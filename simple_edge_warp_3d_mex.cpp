#include <math.h>
#include "mex.h"
#include "matrix.h"
#include "simpleEdge3d_l.h"

struct imorder
{
	unsigned int coord;
	float val;
};

int queue_comp(const void *queue1_, const void *queue2_)
{
	struct imorder *queue1 = (struct imorder *)queue1_;
	struct imorder *queue2 = (struct imorder *)queue2_;
	
	if(queue1->val > queue2->val)
		return -1;
	else if(queue1->val < queue2->val)
		return 1;
	return 0;
}

void construct_patch(unsigned int *source, bool *patch, int i, int xsize, int ysize, int zsize)
{
	int x, y, z, d;
	for(d=0;d<3;d++)
		for(x=0;x<3;x++)
			for(y=0;y<3;y++)
				for(z=0;z<3;z++)
					patch[x + y*3 + z*3*3 + d*3*3*3] = (bool)(source[ i%(xsize*ysize*zsize) + (x-1) + (y-1)*xsize + (z-1)*xsize*ysize + d*xsize*ysize*zsize]!=0);
}

unsigned int unique_neighbor(unsigned int *source, int i, int xsize, int ysize, int zsize, unsigned int fg_conn)
{
	int it, j;
	unsigned int neighbor = 0;
	unsigned int v = (unsigned int)(i/(xsize*ysize*zsize));
	if(fg_conn==6)
	{
		for(j=0;j<10;j++)
		{
			it = i%(xsize*ysize*zsize) + offset10[j][v][0] + offset10[j][v][1]*xsize + offset10[j][v][2]*xsize*ysize + offset10[j][v][3]*xsize*ysize*zsize;
			
			if(source[it]!=0)
			{
				neighbor=(unsigned int)source[it];
				break;
			}
		}
	}
	else
	{
		for(j=0;j<12;j++)
		{
			it = i%(xsize*ysize*zsize) + offset12[j][v][0] + offset12[j][v][1]*xsize + offset12[j][v][2]*xsize*ysize + offset12[j][v][3]*xsize*ysize*zsize;
			if(source[it]!=0)
			{
				neighbor=(unsigned int)source[it];
				break;
			}
		}
	}
	return neighbor;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int i, j, k, it, x, y, z, d, linind, linind_pad, imsize, xsize, ysize, zsize, dsize, imsize_pad, erroro, errorn, delregct, delregbct;
	float binary_threshold;
	unsigned int fg_conn, bg_conn, v;
	bool *mask, *target, *orig_mask, *missclass_points_image, *patch, *queued;
	unsigned int *orig_source, *source, *sourceout;
	float *orig_target, *target_real;
	struct imorder *delreg, *delregb;
	
	imsize = (int)mxGetNumberOfElements(prhs[0]);
	xsize = (int)mxGetDimensions(prhs[0])[0];
	ysize = (int)mxGetDimensions(prhs[0])[1];
	zsize = (int)mxGetDimensions(prhs[0])[2];
	dsize = (int)mxGetDimensions(prhs[0])[3];
	imsize_pad = (int)(xsize+2)*(ysize+2)*(zsize+2)*dsize;
		
	binary_threshold = (float)((nrhs==3)?0.5:mxGetScalar(prhs[3]));
	fg_conn = (unsigned int)((nrhs==4)?6:mxGetScalar(prhs[4]));
	
	if(fg_conn==6)
		bg_conn=26;
	else
		bg_conn=6;
	
	orig_source = (unsigned int *)mxGetData(prhs[0]);
	orig_target = (float *)mxGetData(prhs[1]);
	orig_mask = (bool *)mxGetData(prhs[2]);
		
	source = (unsigned int *)mxCalloc(imsize_pad,sizeof(unsigned int));
	target_real = (float *)mxCalloc(imsize_pad,sizeof(float));
	mask = (bool *)mxCalloc(imsize_pad,sizeof(bool));
	target = (bool *)mxCalloc(imsize_pad,sizeof(bool));
	missclass_points_image = (bool *)mxCalloc(imsize_pad,sizeof(bool));
	delreg = (struct imorder *)mxCalloc(imsize_pad,sizeof(struct imorder));
	delregb = (struct imorder *)mxCalloc(imsize_pad,sizeof(struct imorder));
	patch = (bool *)mxCalloc(81,sizeof(bool));
	queued = (bool *)mxCalloc(imsize_pad,sizeof(bool));
		
	for(d=0;d<dsize;d++)
	{
		for(x=0;x<xsize;x++)
		{
			for(y=0;y<ysize;y++)
			{
				for(z=0;z<zsize;z++)
				{
					linind_pad = (x+1) + (y+1)*(xsize+2) + (z+1)*(xsize+2)*(ysize+2) + d*(xsize+2)*(ysize+2)*(zsize+2);
					linind = x + y*xsize + z*xsize*ysize + d*xsize*ysize*zsize;
					source[linind_pad] = (unsigned int)(orig_source[linind]);
					target_real[linind_pad] = (float)(orig_target[linind]);
					mask[linind_pad] = (bool)(orig_mask[linind]);
				}
			}
		}
	}
		
	erroro = -1;
	errorn = 0;
	
	bool bg,found;
	
	j=0;
	for(i=0;i<imsize_pad;i++)
	{
		target[i] = (bool)(target_real[i]>binary_threshold);
		
		missclass_points_image[i] = (bool)(mask[i] && ((source[i]!=0)!=target[i]));
		
		if(missclass_points_image[i])
		{
			v = (unsigned int)(i/((xsize+2)*(ysize+2)*(zsize+2)));
			errorn++;
			it = i%((xsize+2)*(ysize+2)*(zsize+2)) + offset10e[0][v][0] + offset10e[0][v][1]*(xsize+2) + offset10e[0][v][2]*(xsize+2)*(ysize+2) + offset10e[0][v][3]*(xsize+2)*(ysize+2)*(zsize+2);
			bg = (source[it] == 0);
			found = false;
			if((bg && fg_conn==6) || (!bg && fg_conn==26))
			{
				for(k=1;k<10;k++)
				{
					it = i%((xsize+2)*(ysize+2)*(zsize+2)) + offset10[k][v][0] + offset10[k][v][1]*(xsize+2) + offset10[k][v][2]*(xsize+2)*(ysize+2) + offset10[k][v][3]*(xsize+2)*(ysize+2)*(zsize+2);
					if((bg && source[it]!=0) || (!bg && source[it]==0))
					{
						found=true;
						break;
					}
				}
			}
			else
			{
				for(k=1;k<12;k++)
				{
					it = i%((xsize+2)*(ysize+2)*(zsize+2)) + offset12[k][v][0] + offset12[k][v][1]*(xsize+2) + offset12[k][v][2]*(xsize+2)*(ysize+2) + offset12[k][v][3]*(xsize+2)*(ysize+2)*(zsize+2);
					if((bg && source[it]!=0) || (!bg && source[it]==0))
					{
						found=true;
						break;
					}
				}
			}
			if(found)
			{
				delreg[j].coord = i;
				delreg[j].val = (float)fabs(target_real[i]-binary_threshold);
				j++;
			}
		}
	}
	delregct = j;
	
	
	bool changed=true;
	while( changed )
	{
		changed=false;
		qsort(delreg,delregct,sizeof(struct imorder),queue_comp);
		for(i=0;i<delregct;i++)
			delregb[i] = delreg[i];
		
		for(i=0;i<imsize_pad;i++)
			queued[i] = false;
		
		delregbct = delregct;
		delregct = 0;
		for(k=0;k<delregbct;k++)
		{
			i = delregb[k].coord;
			if(missclass_points_image[i])
			{
				construct_patch(source,patch,i,xsize+2,ysize+2,zsize+2);
				if(simpleEdge3d(patch,(unsigned int)(i/((xsize+2)*(ysize+2)*(zsize+2)) +1),fg_conn))
				{
					source[i] = source[i]==0?unique_neighbor(source,i,xsize+2,ysize+2,zsize+2,fg_conn):0;
					
					missclass_points_image[i] = false;
					errorn--;
					changed = true;
					v = (unsigned int)(i/((xsize+2)*(ysize+2)*(zsize+2)));
					for(j=0;j<34;j++) /* ??  10 aligned neighbors, 12 and 12 antialigned neighbors  ?? */
					{
						it = i%((xsize+2)*(ysize+2)*(zsize+2)) + offset12e[j][v][0] + offset12e[j][v][1]*(xsize+2) + offset12e[j][v][2]*(xsize+2)*(ysize+2) + offset12e[j][v][3]*(xsize+2)*(ysize+2)*(zsize+2);
						if(!queued[it])
						{
							queued[it] = true;
							delreg[delregct].coord = it;
							delreg[delregct].val = (float)fabs(target_real[it]-binary_threshold);
							delregct++;
						}
					}
				}
			}
		}
	}
	
	plhs[0] = mxCreateNumericArray(mxGetNumberOfDimensions(prhs[0]),mxGetDimensions(prhs[0]),mxUINT32_CLASS,mxREAL);
	sourceout = (unsigned int *)mxGetData(plhs[0]);
	
	for(d=0;d<dsize;d++)
		for(x=0;x<xsize;x++)
			for(y=0;y<ysize;y++)
				for(z=0;z<zsize;z++)
					sourceout[x + y*xsize + z*xsize*ysize + d*xsize*ysize*zsize] = source[(x+1) + (y+1)*(xsize+2) + (z+1)*(xsize+2)*(ysize+2) + d*(xsize+2)*(ysize+2)*(zsize+2)];
	
	mxFree(source);
	mxFree(target_real);
	mxFree(mask);
	mxFree(target);
	mxFree(missclass_points_image);
	mxFree(patch);
	mxFree(queued);
	mxFree(delreg);
	mxFree(delregb);
	
	return;
}
