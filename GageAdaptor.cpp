/* 

******************************************************************************

Copyright 2008 Universidade Federal do Rio Grande do Sul, Carlos Dietrich

This file is part of Macet.

Macet is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published
by the Free Software Foundation; either version 2 of the License,
or (at your option) any later version.

Macet is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301, USA

******************************************************************************

If you use this work in academic papers, we would really appreciate if
you cited either of these two:

Dietrich et al. Edge Groups: an approach to understanding the mesh
quality of marching methods. IEEE Trans. Vis. Comp. Graph. 2008

Dietrich et al. Edge transformations for improving the quality of
marching methods. IEEE Trans. Vis Comp. Graph. 2009

******************************************************************************

*/

#define MY_LEAN_AND_MEAN_GAGEADAPTOR

#include <cmath>
#include <iostream>

#include "GageAdaptor.h"

/**
*/
CGageAdaptor::CGageAdaptor(void)
{
	Create();
}

/**
*/
CGageAdaptor::CGageAdaptor(const std::string& path)
{
	Create();

	Open(path);
}

/**
*/
CGageAdaptor::~CGageAdaptor(void)
{
	if (IsOpen())
		Close();
}

/**
*/
bool CGageAdaptor::Open(const std::string& path)
{
	if (IsOpen())
		Close();

	if (!OpenNrrd(path))
	{
		std::cerr << "opennrrd failed..." << std::endl;
		return false;
	}

	if (!CreateDefaultContext())
	{
		std::cerr << "create default context failed..." << std::endl;
		return false;
	}

	m_isOpen = true;

	m_isCopy = false;

	return true;
}

/**
*/
bool CGageAdaptor::OpenFromMemory(void *data, VALUE_TYPE type, unsigned int width, unsigned int height, unsigned int depth)
{
	if (IsOpen())
		Close();

	if (!OpenNrrdFromMemory(data, type, width, height, depth))
	{
		return false;
	}

	if (!CreateDefaultContext())
	{
		return false;
	}

	m_isOpen = true;

	m_isCopy = false;

	return true;
}

/**
*/
void CGageAdaptor::Close(void)
{
	if (m_imageInfo && !m_isCopy) 
	{ 
		if (gagePerVolumeDetach(m_measurementContext, m_imageInfo))

		free(gagePerVolumeNix(m_imageInfo));

		m_imageInfo = 0;
	}

	if (m_measurementContext)
	{
		free(gageContextNix(m_measurementContext));

		m_measurementContext = 0;
	}

	if (m_imageHandle)
	{
		nrrdNuke(m_imageHandle);
		
		m_imageHandle = 0;
	}

	m_isOpen = false;

	m_isCopy = false;
}

/**
*/
bool CGageAdaptor::IsOpen(void) const
{
	return m_isOpen;
}

/**
*/
bool CGageAdaptor::EnableQuery(int item)
{
	if (m_isCopy)
	{
		return false;
	}

	if (!m_measurementContext || !m_imageInfo)
	{
		return false;
	}

	if (gageQueryItemOn(m_measurementContext, m_imageInfo, item))
	{
		return false;
	}

	switch (item) {
		case gageSclValue:
			if (!(m_valuePointer = gageAnswerPointer(m_measurementContext, m_imageInfo, gageSclValue)))
			{
				return false;
			}
			break;
		case gageSclNormal:
			if (!(m_normalPointer = gageAnswerPointer(m_measurementContext, m_imageInfo, gageSclNormal)))
			{
				return false;
			}
			break;
		case gageSclGradVec:
			if (!(m_gradientPointer = gageAnswerPointer(m_measurementContext, m_imageInfo, gageSclGradVec)))
			{
				return false;
			}
			break;
		case gageSclGradMag:
			if (!(m_gradientMagnitudePointer = gageAnswerPointer(m_measurementContext, m_imageInfo, gageSclGradMag)))
			{
				return false;
			}
			break;
		case gageSclHessian:
			if (!(m_hessianPointer = gageAnswerPointer(m_measurementContext, m_imageInfo, gageSclHessian)))
			{
				return false;
			}
			break;
		case gageSclLaplacian:
			if (!(m_laplacianPointer = gageAnswerPointer(m_measurementContext, m_imageInfo, gageSclLaplacian)))
			{
				return false;
			}
			break;
		case gageSclHessEval0:
			if (!(m_hessian1stEigenvaluePointer = gageAnswerPointer(m_measurementContext, m_imageInfo, gageSclHessEval0)))
			{
				return false;
			}
			break;
		case gageSclHessEval1:
			if (!(m_hessian2ndEigenvaluePointer = gageAnswerPointer(m_measurementContext, m_imageInfo, gageSclHessEval1)))
			{
				return false;
			}
			break;
		case gageSclHessEval2:
			if (!(m_hessian3rdEigenvaluePointer = gageAnswerPointer(m_measurementContext, m_imageInfo, gageSclHessEval2)))
			{
				return false;
			}
			break;
		case gageSclK1:
			if (!(m_1stPrincipleCurvaturePointer = gageAnswerPointer(m_measurementContext, m_imageInfo, gageSclK1)))
			{
				return false;
			}
			break;
		default:
			return false;
	}

	if (!UpdateKernel())
	{
		return false;
	}

	return true;
}

/**
*/
bool CGageAdaptor::ResetKernel(void)
{
	if (m_isCopy)
	{
		return false;
	}

	if (m_measurementContext)
	{
		gageKernelReset(m_measurementContext);
		
		if (!UpdateKernel())
		{
			return false;
		}
	}

	return true;
}

/**
This method is used to indicate that the kernel is intended to be used with 
clamping.
*/
void CGageAdaptor::SetClamp(bool doClamp)
{
	m_doClamp = doClamp;
}

/**
*/
bool CGageAdaptor::SetValueKernel(const NrrdKernel *type, const double *parameters)
{
	if (!m_measurementContext)
	{
		std::cerr << "m_measurementContext is null: failed." << std::endl;
		return false;
	}

	if (m_isCopy)
	{
		std::cerr << "m_isCopy is true: failed." << std::endl;
		return false;
	}

	if (gageKernelSet(m_measurementContext, gageKernel00, type, parameters))
	{
		std::cerr << "gageKernelSet failed." << std::endl;
		return false;
	}

	// cscheid 20081027 With teem 1.10, UpdateKernel does not work
	// here when context is being initialized, since gage needs the query items which
        // will not have been set

	return true;
}

/**
*/
bool CGageAdaptor::Set1stDerivativeKernel(const NrrdKernel *type, const double *parameters)
{
	if (!m_measurementContext)
	{
		return false;
	}

	if (m_isCopy)
	{
		return false;
	}

	if (gageKernelSet(m_measurementContext, gageKernel11, type, parameters))
	{
		return false;
	}

	return UpdateKernel();
}

/**
*/
bool CGageAdaptor::Set2ndDerivativeKernel(const NrrdKernel *type, const double *parameters)
{
	if (!m_measurementContext)
	{
		return false;
	}

	if (m_isCopy)
	{
		return false;
	}

	if (gageKernelSet(m_measurementContext, gageKernel22, type, parameters))
	{
        return false;
	}

	return UpdateKernel();
}

/**
*/
int CGageAdaptor::GetWidth(void) const
{
	/*if (!m_imageHandle)
	{
		MarkError();
	
		return 0;
	}

	return (int)m_imageHandle->axis[0].size;*/
	if (!m_measurementContext)
	{
		return 0;
	}

	return (int)m_measurementContext->shape->size[0];
}

/**
*/
int CGageAdaptor::GetHeight(void) const
{
	/*if (!m_imageHandle)
	{
		MarkError();
	
		return 0;
	}

	return (int)m_imageHandle->axis[1].size;*/
	if (!m_measurementContext)
	{
		return 0;
	}

	return (int)m_measurementContext->shape->size[1];
}

/**
*/
int CGageAdaptor::GetDepth(void) const
{
	/*if (!m_imageHandle)
	{
		MarkError();
	
		return 0;
	}

	return (int)m_imageHandle->axis[2].size;*/
	if (!m_measurementContext)
	{
		return 0;
	}

	return (int)m_measurementContext->shape->size[2];
}

/**
*/
CGageAdaptor::VALUE_TYPE CGageAdaptor::GetType(void) const
{
	switch (m_imageHandle->type) {
		case nrrdTypeChar:
			return BYTE;
			break;
		case nrrdTypeUChar:
			return UNSIGNED_BYTE;
			break;
		case nrrdTypeShort:
			return SHORT;
		case nrrdTypeUShort:
			return UNSIGNED_SHORT;
			break;
		case nrrdTypeInt:
			return INT;
			break;
		case nrrdTypeUInt:
			return UNSIGNED_INT;
			break;
		case nrrdTypeFloat:
			return FLOAT;
			break;
		case nrrdTypeDouble:
			return DOUBLE;
			break;
	}

	return UNKNOWN_TYPE;
}

/**
*/
const GAGE_TYPE *CGageAdaptor::GetValueArray(void) const
{
#ifndef MY_LEAN_AND_MEAN_GAGEADAPTOR
	if (!m_imageInfo)
	{
		MarkError();
	
		return 0;
	}
#endif // #ifndef MY_LEAN_AND_MEAN_GAGEADAPTOR

	return (GAGE_TYPE*)m_imageInfo->nin->data;
}

/**
*/
GAGE_TYPE CGageAdaptor::GetValue(float x, float y, float z) const
{
#ifndef MY_LEAN_AND_MEAN_GAGEADAPTOR
	if (!m_measurementContext || !m_valuePointer)
	{
		MarkError();
	
		return 0;
	}

	if (m_doClamp)
		Clamp(&x, &y, &z);
#endif // #ifndef MY_LEAN_AND_MEAN_GAGEADAPTOR

    gageProbe(m_measurementContext, x, y, z);

    return *m_valuePointer;
}

/**
*/
const GAGE_TYPE *CGageAdaptor::GetNormal(float x, float y, float z) const
{
#ifndef MY_LEAN_AND_MEAN_GAGEADAPTOR
	if (!m_measurementContext || !m_normalPointer)
	{
		MarkError();
	
		return 0;
	}

	if (m_doClamp)
		Clamp(&x, &y, &z);
#endif // #ifndef MY_LEAN_AND_MEAN_GAGEADAPTOR

    gageProbe(m_measurementContext, x, y, z);

	return m_normalPointer;
}

/**
*/
const GAGE_TYPE *CGageAdaptor::GetNormal(void) const
{
	return m_normalPointer;
}

/**
*/
const GAGE_TYPE *CGageAdaptor::GetGradient(float x, float y, float z) const
{
	if (!m_measurementContext || !m_gradientPointer)
	{
		return 0;
	}

	if (m_doClamp)
		Clamp(&x, &y, &z);

    gageProbe(m_measurementContext, x, y, z);
	return m_gradientPointer;
}

/**
*/
GAGE_TYPE CGageAdaptor::GetGradientMagnitude(float x, float y, float z) const
{
	if (!m_measurementContext || !m_gradientMagnitudePointer)
	{
		return 0;
	}

	if (m_doClamp)
		Clamp(&x, &y, &z);

    gageProbe(m_measurementContext, x, y, z);

	return *m_gradientMagnitudePointer;
}

/**
*/
const GAGE_TYPE *CGageAdaptor::GetHessian(float x, float y, float z) const
{
	if (!m_measurementContext || !m_hessianPointer)
	{
		return 0;
	}

	if (m_doClamp)
		Clamp(&x, &y, &z);

    gageProbe(m_measurementContext, x, y, z);

	return m_hessianPointer;
}

/**
*/
GAGE_TYPE CGageAdaptor::GetLaplacian(float x, float y, float z) const
{
	if (!m_measurementContext || !m_laplacianPointer)
	{
		return 0;
	}

	if (m_doClamp)
		Clamp(&x, &y, &z);

    gageProbe(m_measurementContext, x, y, z);

	return *m_laplacianPointer;
}

/**
*/
GAGE_TYPE CGageAdaptor::GetHessian1stEigenvalue(float x, float y, float z) const
{
	if (!m_measurementContext || !m_hessian1stEigenvaluePointer)
	{
		return 0;
	}

	if (m_doClamp)
		Clamp(&x, &y, &z);

    gageProbe(m_measurementContext, x, y, z);

	return *m_hessian1stEigenvaluePointer;
}

/**
*/
GAGE_TYPE CGageAdaptor::GetHessian2ndEigenvalue(float x, float y, float z) const
{
	if (!m_measurementContext || !m_hessian2ndEigenvaluePointer)
	{
		return 0;
	}

	if (m_doClamp)
		Clamp(&x, &y, &z);

    gageProbe(m_measurementContext, x, y, z);

	return *m_hessian2ndEigenvaluePointer;
}

/**
*/
GAGE_TYPE CGageAdaptor::GetHessian3rdEigenvalue(float x, float y, float z) const
{
	if (!m_measurementContext || !m_hessian3rdEigenvaluePointer)
	{
		return 0;
	}

	if (m_doClamp)
		Clamp(&x, &y, &z);

    gageProbe(m_measurementContext, x, y, z);

	return *m_hessian3rdEigenvaluePointer;
}

/**
*/
GAGE_TYPE CGageAdaptor::Get1stPrincipalCurvature(float x, float y, float z) const
{
	if (!m_measurementContext || !m_1stPrincipleCurvaturePointer)
	{
		return 0;
	}

	if (m_doClamp)
		Clamp(&x, &y, &z);

    gageProbe(m_measurementContext, x, y, z);

	return *m_1stPrincipleCurvaturePointer;
}

/**
*/
bool CGageAdaptor::OpenNrrd(const std::string& path)
{
	if (IsOpen())
	{
		return false;
	}

	m_imageHandle = nrrdNew();

	nrrdStateDisableContent = AIR_TRUE;

	if (nrrdLoad(m_imageHandle, path.c_str(), NULL)) 
	{
		Close();
		
		return false;
	}

	return true;
}

/**
*/
bool CGageAdaptor::OpenNrrdFromMemory(void *data, VALUE_TYPE type, unsigned int width, unsigned int height, unsigned int depth)
{
	if (IsOpen())
	{
		return false;
	}

	m_imageHandle = nrrdNew();

	nrrdStateDisableContent = AIR_TRUE;

	if (nrrdAlloc_va(m_imageHandle, type, 3, width, height, depth))
	{
		Close();
		
		return false;
	}

	switch (type) {
		case nrrdTypeFloat:
			memcpy(m_imageHandle->data, data, width*height*depth*sizeof(float));
			break;
		case nrrdTypeDouble:
			memcpy(m_imageHandle->data, data, width*height*depth*sizeof(double));
			break;
		default:
			return false;
	};

	// Why do I have to do it? Gage cannot do it automatically?
	m_imageHandle->axis[0].spacing = 1.0;
	m_imageHandle->axis[1].spacing = 1.0;
	m_imageHandle->axis[2].spacing = 1.0;

	return true;
}

/**
*/
bool CGageAdaptor::CreateDefaultContext(void)
{
	double parameters[3];
	bool status;
	
	// Scale.
	parameters[0] = 1.0f;
	// Don't care.
	parameters[1] = 0.0f;
	// Don't care.
	parameters[2] = 0.0f;

	

	if (!m_imageHandle)
	{
		std::cerr << "No image handle..." << std::endl;
		return false;
	}

	if (!(m_measurementContext = gageContextNew()))
	{
		std::cerr << "gageContextNew failed..." << std::endl;
		return false;
	}

	if (!(m_imageInfo = gagePerVolumeNew(m_measurementContext, m_imageHandle, gageKindScl)))
	{
        std::cerr << "gagePerVolumeNew failed..." << std::endl;
		return false;
	}

	if (gagePerVolumeAttach(m_measurementContext, m_imageInfo))
	{
		std::cerr << "gagePerVolumeAttach failed..." << std::endl;
		Close();
		return false;
	}

	// The tent function: f(-1)=0, f(0)=1, f(1)=0, with linear ramps in 
	// between, and zero elsewhere. Used for linear (and bilinear and 
	// trilinear) interpolation.

	if (!SetValueKernel(nrrdKernelTent, parameters))
	{
		std::cerr << "SetValueKernel failed..." << std::endl;
		Close();
		return false;
	}
	// Piecewise-linear ramps that implement forward-difference 
	// differentiation.
	//if (status)
	//	status = Set1stDerivativeKernel(nrrdKernelForwDiff, parameters);

	if (!EnableQuery(gageSclValue))
	{
		std::cerr << "EnableQuery failed..." << std::endl;
		Close();
		return false;
	}

	//if (status)
	//	status = EnableQuery(gageSclNormal);

	if (!UpdateKernel())
	{
		std::cerr << "UpdateKernel failed..." << std::endl;
		Close();
		return false;
	}

	if (!(m_valuePointer = gageAnswerPointer(m_measurementContext, m_imageInfo, gageSclValue)))
	{
		std::cerr << "gageAnswerPointer failed..." << std::endl;
		return false;
	}

	//if (!(m_normalPointer = gageAnswerPointer(m_measurementContext, m_imageInfo, gageSclNormal)))
	//{
	//	MarkError();
	//
	//	return false;
	//}

	return true;
}

/**
*/
bool CGageAdaptor::UpdateKernel(void)
{
	if (!m_measurementContext)
	{
		std::cerr << "m_measurementContext is NULL: UpdateKernel failed" << std::endl;
		return false;
	}

	if (m_isCopy)
	{
		std::cerr << "m_isCopy: UpdateKernel failed" << std::endl;
		return false;
	}

	if (gageUpdate(m_measurementContext))
	{
		std::cerr << "gageUpdate failed: UpdateKernel failed" << std::endl;
		std::cerr << biffGetDone(GAGE) << std::endl;

		return false;
	}

	return true;
}

/**
*/
void CGageAdaptor::Clamp(float *x, float *y, float *z) const
{
	if (*x < FLT_EPSILON)
		*x = FLT_EPSILON;
	else if (*x > (m_measurementContext->shape->size[0] - 1.0f - FLT_EPSILON))
		*x = m_measurementContext->shape->size[0] - 1.0f - FLT_EPSILON;

	if (*y < FLT_EPSILON)
		*y = FLT_EPSILON;
	else if (*y > (m_measurementContext->shape->size[1] - 1.0f - FLT_EPSILON))
		*y = m_measurementContext->shape->size[1] - 1.0f - FLT_EPSILON;

	if (*z < FLT_EPSILON)
		*z = FLT_EPSILON;
	else if (*z > (m_measurementContext->shape->size[2] - 1.0f - FLT_EPSILON))
		*z = m_measurementContext->shape->size[2] - 1.0f - FLT_EPSILON;
}

/**
*/
void CGageAdaptor::Create(void)
{
	m_imageHandle = 0;

	m_measurementContext = 0;
	m_imageInfo = 0;
	
	m_valuePointer = 0;
	m_normalPointer = 0;

	m_gradientPointer = 0;
	m_gradientMagnitudePointer = 0;

	m_hessianPointer = 0;

	m_laplacianPointer = 0;

	m_hessian1stEigenvaluePointer = 0;
	m_hessian2ndEigenvaluePointer = 0;
	m_hessian3rdEigenvaluePointer = 0;

	m_1stPrincipleCurvaturePointer = 0;

	m_isOpen = false;

	m_isCopy = false;

	m_doClamp = false;
}

