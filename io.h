#ifndef IO_H
#define IO_H

#include <teem/nrrd.h>
#include <string>

#include "GageAdaptor.h"

/**
*/
boost::shared_ptr<CGageAdaptor> LoadImage(const std::string& imageFileName)
{
    boost::shared_ptr<CGageAdaptor> m_image;
    double kernelParam[3];

    // Scale parameter, in units of samples.
    kernelParam[0] = 1.0;
    // The two-parameter family of first-order continuous, 4-sample support,
    // piece-wise cubic splines, described by "D. P. Mitchell and A. N.
    // Netravali. Reconstruction filters in computer graphics. In Computer
    // Graphics (SIGGRAPH '88 Proceedings), volume 22, pages 221--228, August
    // 1988." Includes the Catmull-Rom spline at (B,C)=(0,0.5) and the uniform
    // cubic B-spline at (B,C)=(1,0).
    kernelParam[1] = 0.0;
    kernelParam[2] = 0.5;

    m_image.reset(new CGageAdaptor);

    if (!m_image.get())
    {
        std::cerr << "get image failed..." << std::endl;
        return false;
    }

    if (!m_image->Open(imageFileName))
    {
        std::cerr << "open image failed..." << std::endl;
        return false;
    }

    // BC family of cubic polynomial splines.
    if (!m_image->SetValueKernel(nrrdKernelBCCubic, kernelParam))
    {
        std::cerr << "set value kernel failed..." << std::endl;
    }

    // 1st deriv. of BC cubic family.
    if (!m_image->Set1stDerivativeKernel(nrrdKernelBCCubicD, kernelParam))
    {
        std::cerr << "set derivative kernel failed..." << std::endl;
    }

    if (!m_image->EnableQuery(CGageAdaptor::NORMAL))
    {
        std::cerr << "enable query failed..." << std::endl;
    }

    return m_image;
}

//Nrrd* open(char *filename)
//{
//  char *err;
//  Nrrd *nin;

//  /* create a nrrd; at this point this is just an empty container */
//  nin = nrrdNew();

//  /* read in the nrrd from file */
//  if (nrrdLoad(nin, filename, NULL))
//  {
//    err = biffGetDone(NRRD);
//    fprintf(stderr, "Trouble reading \"%s\":\n%s", filename, err);
//    free(err);
//    return;
//  }

//  /* say something about the array */
//  printf("\"%s\" is a %d-dimensional nrrd of type %d (%s)\n",
//         filename, nin->dim, nin->type,
//         airEnumStr(nrrdType, nin->type));
//  printf("The array contains %d elements, each %d bytes in size\n",
//         (int)nrrdElementNumber(nin), (int)nrrdElementSize(nin));

//  return nin;
//}

#endif // IO_H
