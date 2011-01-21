#include <math.h>
#include "mex.h"
#include "matrix.h"
#include "num_comps_edge.h"
#include <pthread.h>
#include <algorithm>

#define NUM_THREADS 8

static int offset6[10][3][4] = {{ {-1,0,0,0},	{0,-1,0,1},	{0,0,-1,2} },
								{ {1,0,0,0},	{0,1,0,1},	{0,0,1,2} },
								{ {0,0,0,1},	{0,0,0,0},	{0,0,0,0} },
								{ {0,0,0,2},	{0,0,0,2},	{0,0,0,1} },
								{ {0,1,0,1},	{1,0,0,0},	{0,1,0,1} },
								{ {0,0,1,2},	{0,0,1,2},	{1,0,0,0} },
								{ {-1,0,0,1},	{0,-1,0,0},	{0,0,-1,0} },
								{ {-1,0,0,2},	{0,-1,0,2},	{0,0,-1,1} },
								{ {-1,1,0,1},	{1,-1,0,0},	{0,1,-1,1} },
								{ {-1,0,1,2},	{0,-1,1,2},	{1,0,-1,0} }};

static int offset26[12][3][4] = {{ {0,0,0,1},	{0,0,0,0},	{0,0,0,0} },
								 { {0,0,0,2},	{0,0,0,2},	{0,0,0,1} },
								 { {0,1,0,1},	{1,0,0,0},	{0,1,0,1} },
								 { {0,0,1,2},	{0,0,1,2},	{1,0,0,0} },
								 { {-1,0,0,1},	{0,-1,0,0},	{0,0,-1,0} },
								 { {-1,0,0,2},	{0,-1,0,2},	{0,0,-1,1} },
								 { {-1,1,0,1},	{1,-1,0,0},	{0,1,-1,1} },
								 { {-1,0,1,2},	{0,-1,1,2},	{1,0,-1,0} },
								 { {0,1,0,0},	{1,0,0,1},	{1,0,0,2} },
								 { {0,-1,0,0},	{-1,0,0,1},	{-1,0,0,2} },
								 { {0,0,1,0},	{0,0,1,1},	{0,1,0,2} },
								 { {0,0,-1,0},	{0,0,-1,1},	{0,-1,0,2} }};

static int offset_all[14][3][4] = {{ {-1,0,0,0},{0,-1,0,1},	{0,0,-1,2} },
								  { {1,0,0,0},	{0,1,0,1},	{0,0,1,2} },
								  { {0,1,0,0},	{1,0,0,1},	{1,0,0,2} },
								  { {0,-1,0,0},	{-1,0,0,1},	{-1,0,0,2} },
								  { {0,0,0,1},	{0,0,0,0},	{0,0,0,0} },
								  { {0,0,0,2},	{0,0,0,2},	{0,0,0,1} },
								  { {0,1,0,1},	{1,0,0,0},	{0,1,0,1} },
								  { {0,0,1,2},	{0,0,1,2},	{1,0,0,0} },
								  { {-1,0,0,1},	{0,-1,0,0},	{0,0,-1,0} },
								  { {-1,0,0,2},	{0,-1,0,2},	{0,0,-1,1} },
								  { {-1,1,0,1},	{1,-1,0,0},	{0,1,-1,1} },
								  { {-1,0,1,2},	{0,-1,1,2},	{1,0,-1,0} },
								  { {0,0,1,0},	{0,0,1,1},	{0,1,0,2} },
								  { {0,0,-1,0},	{0,0,-1,1},	{0,-1,0,2} }};

void construct_patch(unsigned int *source, unsigned int *patch, int i, int xsize, int ysize, int zsize)
{
	int x, y, z, d;
	for(d=0;d<3;d++)
		for(x=0;x<3;x++)
			for(y=0;y<3;y++)
				for(z=0;z<3;z++)
					patch[x + y*3 + z*3*3 + d*3*3*3] = source[ i%(xsize*ysize*zsize) + (x-1) + (y-1)*xsize + (z-1)*xsize*ysize + d*xsize*ysize*zsize];
}

bool neighbor_exists(unsigned int *patch, unsigned int v, bool foreground, int fg_conn)
{
	int i;
	if(fg_conn==6)
	{
		if(foreground)
		{
			for(i=0;i<10;i++)
				if(patch[1+offset6[i][v][0] + (1+offset6[i][v][1])*3 + (1+offset6[i][v][2])*3*3 + offset6[i][v][3]*3*3*3])
					return true;
		}
		else
		{
			for(i=0;i<12;i++)
				if(!patch[1+offset26[i][v][0] + (1+offset26[i][v][1])*3 + (1+offset26[i][v][2])*3*3 + offset26[i][v][3]*3*3*3])
					return true;
		}
	}
	else
	{
		if(foreground)
		{
			for(i=0;i<12;i++)
				if(patch[1+offset26[i][v][0] + (1+offset26[i][v][1])*3 + (1+offset26[i][v][2])*3*3 + offset26[i][v][3]*3*3*3]!=0)
					return true;
		}
		else
		{
			for(i=0;i<10;i++)
				if(patch[1+offset6[i][v][0] + (1+offset6[i][v][1])*3 + (1+offset6[i][v][2])*3*3 + offset6[i][v][3]*3*3*3]==0)
					return true;
		}
	}
	return false;
}

bool multiple_neighbors(unsigned int *patch, unsigned int v, unsigned int fg_conn)
{
	int i;
	unsigned int nbor=0;
	if(fg_conn==6)
	{
		for(i=0;i<10;i++)
		{
			if(patch[(1+offset6[i][v][0]) + (1+offset6[i][v][1])*3 + (1+offset6[i][v][2])*3*3 + offset6[i][v][3]*3*3*3]!=0)
			{
				if(nbor!=0 && nbor!=patch[1+offset6[i][v][0] + (1+offset6[i][v][1])*3 + (1+offset6[i][v][2])*3*3 + offset6[i][v][3]*3*3*3])
					return true;
				else if(nbor==0)
					nbor = patch[1+offset6[i][v][0] + (1+offset6[i][v][1])*3 + (1+offset6[i][v][2])*3*3 + offset6[i][v][3]*3*3*3];
			}
		}
	}
	else
	{
		for(i=0;i<12;i++)
		{
			if(patch[1+offset26[i][v][0] + (1+offset26[i][v][1])*3 + (1+offset26[i][v][2])*3*3 + offset26[i][v][3]*3*3*3]==0)
			{
				if(nbor!=0 && nbor!=patch[1+offset26[i][v][0] + (1+offset26[i][v][1])*3 + (1+offset26[i][v][2])*3*3 + offset26[i][v][3]*3*3*3])
					return true;
				else if(nbor==0)
					nbor = patch[1+offset26[i][v][0] + (1+offset26[i][v][1])*3 + (1+offset26[i][v][2])*3*3 + offset26[i][v][3]*3*3*3];
			}
		}
	}
	return false;
}

pthread_mutex_t iter_mutex;

struct thread_data{
	bool *edges_to_classify;
	unsigned int *segmented_source;
	unsigned int *classify_im;
	int xsize;
	int ysize;
	int zsize;
	int dsize;
	int imsize_pad;
	int radius;
	int fg_conn;
};

int i=0;

void *classify_edges(void *threadarg)
{
	int j, xc, yc, zc, dc, x_begin, x_end, y_begin, y_end, z_begin, z_end, x, y, z, d, num_comps, num_comps_alt;
	struct thread_data *my_data;
	bool *edges_to_classify, *source_alt;
	unsigned int *segmented_source, *classify_im, *patch;
	int radius, fg_conn, xsize, ysize, zsize, dsize, imsize_pad;
	my_data = (struct thread_data *) threadarg;
	edges_to_classify = my_data->edges_to_classify;
	segmented_source = my_data->segmented_source;
	radius = my_data->radius;
	fg_conn = my_data->fg_conn;
	xsize = my_data->xsize;
	ysize = my_data->ysize;
	zsize = my_data->zsize;
	dsize = my_data->dsize;
	imsize_pad = my_data->imsize_pad;
	classify_im = my_data->classify_im;
	
	if(radius==0)
		source_alt = (bool *)calloc((xsize+2)*(ysize+2)*(zsize+2)*dsize,sizeof(bool));
	else
		source_alt = (bool *)calloc((2*radius+1)*(2*radius+1)*(2*radius+1)*(dsize),sizeof(bool));
	
	patch = (unsigned int *)calloc(81,sizeof(unsigned int));
	
	
	unsigned int v;
	int this_iter;
	pthread_mutex_lock(&iter_mutex);
	this_iter=i;
	i++;
	pthread_mutex_unlock(&iter_mutex);
	while(this_iter<imsize_pad)
	{
		if(edges_to_classify[this_iter])
		{

			v=(unsigned int)(this_iter/((xsize+2)*(ysize+2)*(zsize+2)));
			construct_patch(segmented_source,patch,this_iter,xsize+2,ysize+2,zsize+2);
			if(segmented_source[this_iter]==0)
			{
				if(!neighbor_exists(patch,v,true,fg_conn))
					classify_im[this_iter] = 1; /* flipping would create an object */
				else if(!neighbor_exists(patch,v,false,fg_conn))
					classify_im[this_iter] = 2; /* flipping would delete a hole */
				else if(multiple_neighbors(patch,v,fg_conn))
					classify_im[this_iter] = 3; /* flipping would cause a merger */
				else
					classify_im[this_iter] = 4; /* flipping would create a hole */
			}
			else
			{
				if(!neighbor_exists(patch,v,true,fg_conn))
				{
					classify_im[this_iter] = 5; /* flipping would delete an object */
				}
				else if(!neighbor_exists(patch,v,false,fg_conn))
				{
					classify_im[this_iter] = 6; /* flipping would create a hole */
				}
				else
				{
					if(radius==0) 
					{
						for(j=0;j<imsize_pad;j++)
							source_alt[j] = (bool)(segmented_source[j]==segmented_source[i]);
						
						source_alt[this_iter]=!source_alt[this_iter];
						int patchsize[4] = {xsize+2,ysize+2,zsize+2,dsize};
						num_comps_alt=num_comps_edge(source_alt,patchsize,fg_conn,1,false);
					}
					else
					{
						xc = this_iter%(xsize+2);
						yc = ((int)(this_iter/(xsize+2)))%(ysize+2);
						zc = (int)(this_iter/((xsize+2)*(ysize+2)))%(zsize+2);
						dc = (int)(this_iter/((xsize+2)*(ysize+2)*(zsize+2)));
						x_begin = xc-radius<=0?0:xc-radius;
						x_end = xsize+2<xc+radius+1?xsize+2:xc+radius+1;
						
						y_begin = yc-radius<=0?0:yc-radius;
						y_end = ysize+2<yc+radius+1?ysize+2:yc+radius+1;
						
						z_begin = zc-radius<=0?0:zc-radius;
						z_end = zsize+2<zc+radius+1?zsize+2:zc+radius+1;
						
						for(j=0;j<(2*radius+1)*(2*radius+1)*(2*radius+1)*dsize;j++)
							source_alt[j] = false;
						
						for(x=0;x<(x_end-x_begin);x++)
							for(y=0;y<(y_end-y_begin);y++)
								for(z=0;z<(z_end-z_begin);z++)
									for(d=0;d<dsize;d++)
										source_alt[x + (y)*(2*radius+1) + (z)*(2*radius+1)*(2*radius+1) + d*(2*radius+1)*(2*radius+1)*(2*radius+1)] = (bool)(segmented_source[(x_begin+x) + (y_begin+y)*(xsize+2) + (z_begin+z)*(xsize+2)*(ysize+2) + (d)*(xsize+2)*(ysize+2)*(zsize+2)]==segmented_source[this_iter]);
						
						int patchsize[4] = {2*radius+1,2*radius+1,2*radius+1,dsize};
						num_comps = num_comps_edge(source_alt,patchsize,fg_conn,1,false);
						xc = (int)(xc-x_begin);
						yc = (int)(yc-y_begin);
						zc = (int)(zc-z_begin);

						source_alt[xc + yc*(2*radius+1) + zc*(2*radius+1)*(2*radius+1) + dc*(2*radius+1)*(2*radius+1)*(2*radius+1)] = !source_alt[xc + yc*(2*radius+1) + zc*(2*radius+1)*(2*radius+1) + dc*(2*radius+1)*(2*radius+1)*(2*radius+1)];
						num_comps_alt = num_comps_edge(source_alt,patchsize,fg_conn,1,false);
						num_comps_alt = num_comps_alt-num_comps+1;
					}
					if(num_comps_alt>1)
						classify_im[this_iter] = 7; /* flipping would cause a split */
					else if(num_comps_alt==0)
						printf("shouldnt happen\n");
					else
						classify_im[this_iter] = 8; /* flipping would delete a hole */
				}
			}

		}
		pthread_mutex_lock(&iter_mutex);
		this_iter=i;
		i++;
		pthread_mutex_unlock(&iter_mutex);
		
		
	}
	free(patch);
	free(source_alt);
	pthread_exit((void *) 0);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int j, x, y, z, d, xc, yc, zc, dc, x_begin, x_end, y_begin, y_end, z_begin, z_end, xsize, ysize, zsize, dsize, linind, linind_pad, imsize, imsize_pad;
	int radius, fg_conn;
	int num_comps, num_comps_alt;
	bool *edges_to_classify, *orig_edges_to_classify;
	unsigned int *orig_segmented_source, *segmented_source, *classify_im, *patch;
	bool *source_crop_alt, *source_crop;
	unsigned int *classify_im_out;
	double *comps;
	double numcomps;
	
	imsize = (int)mxGetNumberOfElements(prhs[0]);
	xsize = (int)mxGetDimensions(prhs[0])[0];
	ysize = (int)mxGetDimensions(prhs[0])[1];
	zsize = (int)mxGetDimensions(prhs[0])[2];
	dsize = (int)mxGetDimensions(prhs[0])[3];
	imsize_pad = (int)(xsize+2)*(ysize+2)*(zsize+2)*dsize;
	
	radius = (int)((nrhs==2)?0:mxGetScalar(prhs[2]));
	fg_conn = (int)((nrhs==3)?6:mxGetScalar(prhs[3]));
	if(fg_conn!=6 && fg_conn!=26)
		mexErrMsgTxt("fg_conn must be 6 or 26");
	
	orig_segmented_source = (unsigned int *)mxGetData(prhs[0]);
	orig_edges_to_classify = (bool *)mxGetData(prhs[1]);
	
	segmented_source = (unsigned int *)mxCalloc(imsize_pad,sizeof(unsigned int));
	edges_to_classify = (bool *)mxCalloc(imsize_pad,sizeof(bool));
	classify_im = (unsigned int *)mxCalloc(imsize_pad,sizeof(unsigned int));
	
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
					segmented_source[linind_pad] = orig_segmented_source[linind];
					edges_to_classify[linind_pad] = orig_edges_to_classify[linind];
				}
			}
		}
	}
	
	/* this whole inner loop gets sent to the threads which process sequential */
	pthread_t threads[NUM_THREADS];
	int rc;
	int t;
	struct thread_data thread_data_all;
	thread_data_all.edges_to_classify = edges_to_classify;
	thread_data_all.segmented_source = segmented_source;
	thread_data_all.classify_im = classify_im;
	thread_data_all.xsize = xsize;
	thread_data_all.ysize = ysize;
	thread_data_all.zsize = zsize;
	thread_data_all.dsize = dsize;
	thread_data_all.fg_conn = fg_conn;
	thread_data_all.radius = radius;
	thread_data_all.imsize_pad = imsize_pad;
	pthread_attr_t attr;
	void *status;
	
	pthread_mutex_init(&iter_mutex, NULL);
	
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	i=0;
	for(t=0;t<NUM_THREADS;t++)
		pthread_create(&threads[t], &attr, classify_edges, (void *)&thread_data_all);
	
	pthread_attr_destroy(&attr);
	for(t=0; t<NUM_THREADS; t++)
		pthread_join(threads[t], &status);
	
	pthread_mutex_destroy(&iter_mutex);
		
	plhs[0] = mxCreateNumericArray(mxGetNumberOfDimensions(prhs[0]),mxGetDimensions(prhs[0]),mxUINT32_CLASS,mxREAL);
	classify_im_out = (unsigned int *)mxGetData(plhs[0]);
		
	for(x=0;x<xsize;x++)
		for(y=0;y<ysize;y++)
			for(z=0;z<zsize;z++)
				for(d=0;d<dsize;d++)
					classify_im_out[x + y*xsize + z*xsize*ysize + d*xsize*ysize*zsize] = classify_im[(x+1) + (y+1)*(xsize+2) + (z+1)*(xsize+2)*(ysize+2) + d*(xsize+2)*(ysize+2)*(zsize+2)];
		
	mxFree(segmented_source);
	mxFree(edges_to_classify);
	mxFree(classify_im);
	return;
}
