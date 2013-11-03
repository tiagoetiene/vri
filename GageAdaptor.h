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

#ifndef GAGEADAPTOR_INCLUDED
#define GAGEADAPTOR_INCLUDED

#include <string>

#include <teem/nrrd.h>
#include <teem/gage.h>

#if TEEM_VERSION >= 11000
#define GAGE_TYPE double
#elif TEEM_VERSION >= 10900
#define GAGE_TYPE float
#else
#error "Need teem version >= 1.09, got" TEEM_VERSION_STRING
#endif



#include <boost/utility.hpp>
#include <boost/shared_ptr.hpp>


class CGageAdaptor
	: boost::noncopyable
{
public:
	enum QUERY_ITEM {
		// Data value.
		VALUE = gageSclValue,
		// Gradient vector, normalized.
		NORMAL = gageSclNormal,
		// Gradient vector, un-normalized.
		GRADIENT = gageSclGradVec,
		// Gradient magnitude.
		GRADIENT_MAGNITUDE = gageSclGradMag,
		// Hessian (column-order).
		HESSIAN = gageSclHessian,
		// Laplacian: Dxx + Dyy + Dzz.
		LAPLACIAN = gageSclLaplacian,
		// Hessian's 1st eigenvalue.
		HESSIAN_1ST_EIGENVALUE = gageSclHessEval0,
		// Hessian's 2nd eigenvalue.
		HESSIAN_2ND_EIGENVALUE = gageSclHessEval1,
		// Hessian's 3rd eigenvalue.
		HESSIAN_3RD_EIGENVALUE = gageSclHessEval2,
		// 1st principle curvature.
		PRINCIPAL_CURVATURE = gageSclK1
	};
	enum VALUE_TYPE {
		// Signifies "type is unset/unknown".
		UNKNOWN_TYPE = nrrdTypeUnknown,
		// Signed 1-byte integer.
		BYTE = nrrdTypeChar,
		// Unsigned 1-byte integer.
		UNSIGNED_BYTE = nrrdTypeUChar,
		// Signed 2-byte integer.
		SHORT = nrrdTypeShort,
		// Unsigned 2-byte integer.
		UNSIGNED_SHORT = nrrdTypeUShort,
		// Signed 4-byte integer.
		INT = nrrdTypeInt,
		// Unsigned 4-byte integer.
		UNSIGNED_INT = nrrdTypeUInt,
		// 4-byte floating point.
		FLOAT = nrrdTypeFloat,
		// 8-byte floating point.
		DOUBLE = nrrdTypeDouble
	};
	CGageAdaptor(void);
	CGageAdaptor(const std::string& path);
	virtual ~CGageAdaptor(void);
	virtual bool Open(const std::string& path);
	virtual bool OpenFromMemory(void *data, VALUE_TYPE type, unsigned int width, unsigned int height, unsigned int depth);
	virtual void Close(void);
	virtual bool IsOpen(void) const;
	virtual bool EnableQuery(int item);
	virtual bool ResetKernel(void);
	virtual void SetClamp(bool doClamp);
	virtual bool SetValueKernel(const NrrdKernel *type, const double *parameters);
	virtual bool Set1stDerivativeKernel(const NrrdKernel *type, const double *parameters);
	virtual bool Set2ndDerivativeKernel(const NrrdKernel *type, const double *parameters);
	virtual int GetWidth(void) const;
	virtual int GetHeight(void) const;
	virtual int GetDepth(void) const;
	virtual VALUE_TYPE GetType(void) const;
	virtual const GAGE_TYPE *GetValueArray(void) const;
	virtual GAGE_TYPE GetValue(float x, float y, float z) const;
	virtual const GAGE_TYPE *GetNormal(float x, float y, float z) const;
	const GAGE_TYPE *GetNormal(void) const;
	virtual const GAGE_TYPE *GetGradient(float x, float y, float z) const;
	virtual GAGE_TYPE GetGradientMagnitude(float x, float y, float z) const;
	virtual const GAGE_TYPE *GetHessian(float x, float y, float z) const;
	virtual GAGE_TYPE GetLaplacian(float x, float y, float z) const;
	virtual GAGE_TYPE GetHessian1stEigenvalue(float x, float y, float z) const;
	virtual GAGE_TYPE GetHessian2ndEigenvalue(float x, float y, float z) const;
	virtual GAGE_TYPE GetHessian3rdEigenvalue(float x, float y, float z) const;
	virtual GAGE_TYPE Get1stPrincipalCurvature(float x, float y, float z) const;
private:
	bool OpenNrrd(const std::string& path);
	bool OpenNrrdFromMemory(void *data, VALUE_TYPE type, unsigned int width, unsigned int height, unsigned int depth);
	bool CreateDefaultContext(void);
	bool UpdateKernel(void);
	inline void Clamp(float *x, float *y, float *z) const;
protected:
	void Create(void);
protected:
	Nrrd *m_imageHandle;
	gageContext *m_measurementContext;
	gagePerVolume *m_imageInfo;

	const GAGE_TYPE *m_valuePointer;
	const GAGE_TYPE *m_normalPointer;
	const GAGE_TYPE *m_gradientPointer;
	const GAGE_TYPE *m_gradientMagnitudePointer;
	const GAGE_TYPE *m_hessianPointer;
	const GAGE_TYPE *m_laplacianPointer;
	const GAGE_TYPE *m_hessian1stEigenvaluePointer;
	const GAGE_TYPE *m_hessian2ndEigenvaluePointer;
	const GAGE_TYPE *m_hessian3rdEigenvaluePointer;
	const GAGE_TYPE *m_1stPrincipleCurvaturePointer;
	bool m_isOpen;
	bool m_isCopy;
	bool m_doClamp;
};

#endif // GAGEADAPTOR_INCLUDED

