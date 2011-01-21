#include <math.h>
#include "mex.h"
#include "matrix.h"
#include "num_comps_edge.h"

static int patchsize[4] = {3,3,3,3};

int NCCNh1(bool *graphPatch, unsigned int v)
{
    int i;
	for(i=0;i<81;i++)
		graphPatch[i] &= mask10e[v-1][i];
	
	int num_comps = num_comps_edge(graphPatch,patchsize,6,v,true);
	return num_comps;
	/* This is the problem */
}


int NCCNh2(bool *graphPatch, unsigned int v)
{
	int i;
	for(i=0;i<81;i++)
		graphPatch[i] &= mask12e[v-1][i];
	
	int num_comps = num_comps_edge(graphPatch,patchsize,26,v,true);
    return num_comps;
}



bool simpleEdge3d(bool *orig_graphPatch, unsigned int v, unsigned int fg_conn)
{
	int i, dir;
	bool simple;
	bool *graphPatch;
	
	graphPatch=(bool *)mxCalloc(81,sizeof(bool));
	
	if(fg_conn==6)
		for(i=0;i<81;i++)
			graphPatch[i]=orig_graphPatch[i];
	else
		for(i=0;i<81;i++)
			graphPatch[i]=!orig_graphPatch[i];

    int T1=NCCNh1(graphPatch,v);

	if(fg_conn==6)
		for(i=0;i<81;i++)
			graphPatch[i]=!orig_graphPatch[i];
	else
		for(i=0;i<81;i++)
			graphPatch[i]=orig_graphPatch[i];
    
    int T2=NCCNh2(graphPatch,v);
	simple=(T1 == 1 && T2 ==1)?true:false;
	mxFree(graphPatch);
	
	return simple;
}
