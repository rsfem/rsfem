/* Copyright (C) 2022 M. A. Yalçın, İ. Temizer*/
/* This file is distributed under the terms of the GNU General Public License.*/
/* See the file COPYING for license details.*/

#include "xc.h"
#include <stdio.h>
#include <mex.h>
#include <matrix.h>

//LIBXC
void callLibxcforLDA(double *rho, double *exc,double *vxc, mwSize m, mwSize n){
xc_func_type func;
  int i, vmajor, vminor, vmicro, func_id = 1;
  /* Get the libxc version */
  xc_version(&vmajor, &vminor, &vmicro);
  /* Initialize the functional */
  if(xc_func_init(&func, func_id, XC_UNPOLARIZED) != 0){
    fprintf(stderr, "Functional '%d' not found\n", func_id);
  }
  /* Evaluate the energy density, depending on the family */
  switch(func.info->family)
    {
      case XC_FAMILY_LDA:
      xc_lda_exc(&func, (mwSize)m, rho, exc);
      xc_lda_vxc(&func, (mwSize)m, rho, vxc);
      break;
    case XC_FAMILY_GGA:
    case XC_FAMILY_HYB_GGA:
      break;
    }

  xc_func_end(&func);



}

void mexFunction(
    int           nlhs,           /* number of expected outputs */
    mxArray       *plhs[],        /* array of pointers to output arguments */
    int           nrhs,           /* number of inputs */
    const mxArray *prhs[]         /* array of pointers to input arguments */
    )
{

    double *rho;
    size_t sz1;
    size_t sz2;
    
    #if MX_HAS_INTERLEAVED_COMPLEX
    rho = mxGetDoubles(prhs[0]);
    #else
    rho = mxGetPr(prhs[0]);
    #endif
    
    sz1 = mxGetM(prhs[0]);
    sz2 = mxGetN(prhs[0]);
    
    double *exc;
    double *vxc;
    
    plhs[0] = mxCreateDoubleMatrix((mwSize)sz1,(mwSize)sz2,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(sz1,sz2,mxREAL);
    
    #if MX_HAS_INTERLEAVED_COMPLEX
    exc = mxGetDoubles(plhs[0]);
    vxc = mxGetDoubles(plhs[1]);
    #else
    exc = mxGetPr(plhs[0]);
    vxc = mxGetPr(plhs[1]);
    #endif
    
    callLibxcforLDA(rho, exc,vxc,(mwSize)sz1, (mwSize)sz2);
    //exc = rho;
    

  
}

