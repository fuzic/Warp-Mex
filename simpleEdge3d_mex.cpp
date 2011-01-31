#include <stdio.h>
#include <math.h>
#include "mex.h"
#include "matrix.h"
#include "simpleEdge3d_l.h"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
	double *orig_graphPatch, *simple;
	int i, dir;
	if (nrhs!=2)
		mexErrMsgTxt("simple=simpleEdge3d_mex(graphPatch,v)");

    dir=(int)mxGetScalar(prhs[1]);
    
	orig_graphPatch=mxGetPr(prhs[0]);
	plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
	simple=mxGetPr(plhs[0]);
    if (dir<=0 || dir>3)
        mexErrMsgTxt("Expecting v equal to 1, 2, or 3 (y,x,z accordingly)");

	bool *graphPatch=(bool *)mxCalloc(81,sizeof(bool));
	for(i=0;i<81;i++)
		graphPatch[i]=(bool)orig_graphPatch[i]!=0;

    *simple=(double)simpleEdge3d(graphPatch,dir,6);
	
    mxFree(graphPatch);
	
	return;
}
