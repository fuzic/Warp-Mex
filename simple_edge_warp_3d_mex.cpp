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
		for(j=0;j<32;j++)
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
	int i, j, k, it, x, y, z, d, linind, linind_pad, imsize, imsize_pad, delregct, delregbct, mask_dim_ct, max_obj_id,;
	int *size, *mask_dims;
	float binary_threshold;
	unsigned int fg_conn, bg_conn, v;
	bool *mask, *target, *missclass_points_image, *patch, *queued;
	unsigned int *orig_source, *source, *sourceout;
	float *orig_target, *target_real;
	struct imorder *delreg, *delregb;
	
	imsize = (int)mxGetNumberOfElements(prhs[0]);
	size = (int *)mxGetDimensions(prhs[0]);
	imsize_pad = (int)(size[0]+2)*(size[1]+2)*(size[2]+2)*size[3];
		
	binary_threshold = (float)((nrhs==3)?0.5:mxGetScalar(prhs[3]));
	fg_conn = (unsigned int)((nrhs==4)?6:mxGetScalar(prhs[4]));
	
	if(fg_conn==6)
		bg_conn=26;
	else
		bg_conn=6;
	
	orig_source = (unsigned int *)mxGetData(prhs[0]);
	orig_target = (float *)mxGetData(prhs[1]);
	mask = (bool *)mxGetData(prhs[2]);
	
	mask_dims = (int *)mxGetDimensions(prhs[2]);
	mask_dim_ct = (int)mxGetNumberOfDimensions(prhs[2]);
		
	source = (unsigned int *)mxCalloc(imsize_pad,sizeof(unsigned int));
	target_real = (float *)mxCalloc(imsize_pad,sizeof(float));
	target = (bool *)mxCalloc(imsize_pad,sizeof(bool));
	missclass_points_image = (bool *)mxCalloc(imsize_pad,sizeof(bool));
	delreg = (struct imorder *)mxCalloc(imsize_pad,sizeof(struct imorder));
	delregb = (struct imorder *)mxCalloc(imsize_pad,sizeof(struct imorder));
	patch = (bool *)mxCalloc(81,sizeof(bool));
	queued = (bool *)mxCalloc(imsize_pad,sizeof(bool));
		
	for(d=0;d<size[3];d++)
	{
		for(x=0;x<size[0];x++)
		{
			for(y=0;y<size[1];y++)
			{
				for(z=0;z<size[2];z++)
				{
					linind_pad = (x+1) + (y+1)*(size[0]+2) + (z+1)*(size[0]+2)*(size[1]+2) + d*(size[0]+2)*(size[1]+2)*(size[2]+2);
					linind = x + y*size[0] + z*size[0]*size[1] + d*size[0]*size[1]*size[2];
					source[linind_pad] = (unsigned int)(orig_source[linind]);
					target_real[linind_pad] = (float)(orig_target[linind]);
					target[linind_pad] = (bool)(target_real[linind_pad]>binary_threshold);
				}
			}
		}
	}
		
	bool bg,found;
	int i_pad;
	bool object_wise_warp=false;
	int obj_id;
	if(mask_dim_ct==5)
	{
		// We should really check the maximum value of the segmented components now to see how many objects we need to allocate
		max_obj_id=0;
		for(i=0;i<linind;i++)
		{
			if(orig_source[i] > max_obj_id)
				max_obj_id = (int)orig_source[i];
		}
		if(mxGetDimensions(prhs[2])[4]<max_obj_id)
			mexErrMsgTxt("You must provide more masks than the highest labeled object in the image.");
		object_wise_warp=true;
	}
	
	int max_obj_ct = object_wise_warp?max_obj_id:1;
	
	for(obj_id=0;obj_id<max_obj_ct;obj_id++)
	{
		j=0;
		for(i=0;i<imsize;i++)// We can redo this so the mask neighborhood doesn't need to be padded
		{
			d=(int)(i/(size[0]*size[1]*size[2]));
			x=i%size[0]+1;
			y=(int)(i/size[0])%size[1]+1;
			z=(int)(i/(size[0]*size[1]))%size[2]+1;
			i_pad = x + y*(size[0]+2) + z*(size[0]+2)*(size[1]+2) + d*(size[0]+2)*(size[1]+2)*(size[2]+2);
			missclass_points_image[i_pad] = (bool)(mask[i+obj_id*imsize] && ((source[i_pad]!=0)!=target[i_pad]));
			
			if(object_wise_warp)
				missclass_points_image[i_pad] &= (source[i_pad]==i+1) || (source[i_pad]==0);
			
			if(missclass_points_image[i_pad])
			{
				v = (unsigned int)(i_pad/((size[0]+2)*(size[1]+2)*(size[2]+2)));
				it = i_pad%((size[0]+2)*(size[1]+2)*(size[2]+2)) + offset10e[0][v][0] + offset10e[0][v][1]*(size[0]+2) + offset10e[0][v][2]*(size[0]+2)*(size[1]+2) + offset10e[0][v][3]*(size[0]+2)*(size[1]+2)*(size[2]+2);
				bg = (source[it] == 0);
				found = false;
				if((bg && fg_conn==6) || (!bg && fg_conn==26))
				{
					for(k=1;k<10;k++)
					{
						it = i_pad%((size[0]+2)*(size[1]+2)*(size[2]+2)) + offset10[k][v][0] + offset10[k][v][1]*(size[0]+2) + offset10[k][v][2]*(size[0]+2)*(size[1]+2) + offset10[k][v][3]*(size[0]+2)*(size[1]+2)*(size[2]+2);
						if(object_wise_warp)
						{
							if((bg && source[it]==(obj_id+1)) || (!bg && source[it]==0))
							{
								found=true;
								break;
							}
						}
						else if((bg && source[it]==(obj_id+1)) || (!bg && source[it]==0))
						{
							found=true;
							break;
						}
					}
				}
				else
				{
					for(k=1;k<32;k++)
					{
						it = i_pad%((size[0]+2)*(size[1]+2)*(size[2]+2)) + offset12[k][v][0] + offset12[k][v][1]*(size[0]+2) + offset12[k][v][2]*(size[0]+2)*(size[1]+2) + offset12[k][v][3]*(size[0]+2)*(size[1]+2)*(size[2]+2);
						if(object_wise_warp)
						{
							if((bg && source[it]!=0) || (!bg && source[it]==0))
							{
								found=true;
								break;
							}
						}
						else if((bg && source[it]!=0) || (!bg && source[it]==0))
						{
							found=true;
							break;
						}
					}
				}
				// We are still only considering i_pad, which we are already sure belongs to the object of interest.
				if(found)
				{
					delreg[j].coord = i_pad;
					delreg[j].val = (float)fabs(target_real[i_pad]-binary_threshold);
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
					construct_patch(source,patch,i,size[0]+2,size[1]+2,size[2]+2);
					if(simpleEdge3d(patch,(unsigned int)(i/((size[0]+2)*(size[1]+2)*(size[2]+2)) +1),fg_conn))
					{
						source[i] = source[i]==0?unique_neighbor(source,i,size[0]+2,size[1]+2,size[2]+2,fg_conn):0;
					
						missclass_points_image[i] = false;
						changed = true;
						v = (unsigned int)(i/((size[0]+2)*(size[1]+2)*(size[2]+2)));
						for(j=0;j<34;j++) /* all potential neighbors */
						{
							it = i%((size[0]+2)*(size[1]+2)*(size[2]+2)) + offset12e[j][v][0] + offset12e[j][v][1]*(size[0]+2) + offset12e[j][v][2]*(size[0]+2)*(size[1]+2) + offset12e[j][v][3]*(size[0]+2)*(size[1]+2)*(size[2]+2);
							if(!queued[it])
							{
								if(object_wise_warp)
								{
									if(source[it]==(obj_id+1) || source[it]==0)
									{
										queued[it] = true;
										delreg[delregct].coord = it;
										delreg[delregct].val = (float)fabs(target_real[it]-binary_threshold);
										delregct++;
									}
								}
								else
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
		}
	}
	
	plhs[0] = mxCreateNumericArray(mxGetNumberOfDimensions(prhs[0]),mxGetDimensions(prhs[0]),mxUINT32_CLASS,mxREAL);
	sourceout = (unsigned int *)mxGetData(plhs[0]);
	
	for(d=0;d<size[3];d++)
		for(x=0;x<size[0];x++)
			for(y=0;y<size[1];y++)
				for(z=0;z<size[2];z++)
					sourceout[x + y*size[0] + z*size[0]*size[1] + d*size[0]*size[1]*size[2]] = source[(x+1) + (y+1)*(size[0]+2) + (z+1)*(size[0]+2)*(size[1]+2) + d*(size[0]+2)*(size[1]+2)*(size[2]+2)];
	
	mxFree(source);
	mxFree(target_real);
	mxFree(target);
	mxFree(missclass_points_image);
	mxFree(patch);
	mxFree(queued);
	mxFree(delreg);
	mxFree(delregb);
	
	return;
}
