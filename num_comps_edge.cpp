#include "num_comps_edge.h"
#include <stdlib.h>

#define REALLOC_FACTOR 2

ufType *uf_init(int hint)
{
	ufType *uf = NULL;
	uf = (ufType *) calloc(1, sizeof(*uf));
	uf->id = (int *) calloc(hint, sizeof(*(uf->id)));
	uf->sz = (int *) calloc(hint, sizeof(*(uf->sz)));
	uf->allocated_length = hint;
	uf->num_nodes = 0;
	uf->num_sets = 0;
	uf->finalized = false;

	return uf;
}

void uf_new_node(ufType *uf)
{
	int init_idx, i;
	if (uf->num_nodes >= uf->allocated_length)
	{
		init_idx = uf->allocated_length;
		uf->allocated_length *= REALLOC_FACTOR;
		uf->id = (int *) realloc(uf->id, uf->allocated_length * sizeof(int));
		uf->sz = (int *) realloc(uf->sz, uf->allocated_length * sizeof(int));
	}
	uf->num_nodes++;
	uf->id[uf->num_nodes - 1] = uf->num_nodes - 1;
	uf->sz[uf->num_nodes - 1] = 1;
}

int uf_find(ufType *uf, int p)
{
	int i, t;
	int *id = uf->id;
	
	for (i = p; i != id[i]; i = id[i])
	{
		t = i;
		i = id[id[t]];
		id[t] = i;
	}
	
	return i;
}

void uf_union(ufType *uf, int p, int q)
{
	if (uf->sz[p] < uf->sz[q])
	{
		uf->id[p] = q;
		uf->sz[q] += uf->sz[p];
	}
	else
	{
		uf->id[q] = p;
		uf->sz[p] += uf->sz[q];
	}
}

void uf_new_pair(ufType *uf, int p, int q)
{
	int i, j;
	
	i = uf_find(uf, p);
	j = uf_find(uf, q);
	if (i != j)
		uf_union(uf, i, j);
}

int uf_renumber(ufType *uf, int first)
{
	int k;
	int counter = first;

	uf->finalized = true;

	for (k = 0; k < uf->num_nodes; k++)
		if (uf->id[k] == k)
			uf->sz[k] = counter++;
	
	uf->num_sets = counter - first;
	
	return uf->num_sets;
}

int uf_query_set(ufType *uf, int p)
{
	int k;

	k = uf_find(uf, p);
	return uf->sz[k];
}

void uf_destroy(ufType *uf)
{
	free(uf->id);
	free(uf->sz);
}

bool nbor(int *p, int q, int *size, int fg_conn, int j)
{
	bool inbounds=true;
 	int coord[4];
	int cumsz=1;
	
	for(int i=0; i<4; cumsz*=size[i], i++)
		coord[i] = ((int)(q/cumsz))%size[i];
	
	for(int i=0; i<3; i++)
		coord[i] += (fg_conn==6)?offset10[j][coord[3]][i]:offset12[j][coord[3]][i];
	coord[3] = (fg_conn==6)?offset10[j][coord[3]][3]:offset12[j][coord[3]][3];
	
	for(int i=0; i<3; i++)
		inbounds *= coord[i]>=0 && coord[i]<size[i];
	
	if(inbounds)
		*p = coord[0] + coord[1]*size[0] + coord[2]*size[0]*size[1] + coord[3]*size[0]*size[1]*size[2];
	return inbounds;
}

int num_comps_edge(bool *BW, int *size, int fg_conn, int dir, bool reduced_out)
{
	int num_elements=1;
	int next_label=1;
	int comp_ct=0;
	int p, nct, num_sets, i, q, killidx, found_comp;
	int *L;
	bool found;
	ufType *uf;
		
	for(i=0;i<4;i++)
		num_elements *= size[i];
		
	L = (int *)calloc(num_elements,sizeof(int));
	
	nct = (fg_conn==6)?10:32;
	
	uf = uf_init(1000);
		
	for(q=0; q<num_elements; q++)
	{
		if(BW[q])
		{
			found = false;
			for(i=0; i<nct; i++)
			{
				if(nbor(&p,q,size,fg_conn,i))
				{
					if(L[p]!=0)
					{
						L[q] = L[p];
						found = true;
						killidx = i;
						break;
					}
				}
			}
			if(!found)
			{
				L[q] = next_label;
				next_label++;
				uf_new_node(uf);
			}
			else
			{
				for(i = killidx+1; i<nct; i++)
					if(nbor(&p,q,size,fg_conn,i))
						if((L[p]!=0) && (L[p]!=L[q]))
							uf_new_pair(uf, (L[q]-1), (L[p]-1));
			}
		}
	}
	num_sets = uf_renumber(uf, 1);
	
	if(reduced_out && num_elements==81)
	{
		for(q=0; q<num_elements; q++)
			if(BW[q] != 0)
				L[q] = uf_query_set(uf, (L[q]-1))+1;
		
		found_comp=0;
		
		for(q=0;q<num_elements;q++)
		{
			if((fg_conn==26 && mask12[dir-1][q]) || (fg_conn==6 && mask10[dir-1][q]))
			{
				if(found_comp==0 && L[q]!=0)
				{
					comp_ct=1;
					found_comp=L[q];
				}
				else if(comp_ct==1 && found_comp!=L[q] && L[q]!=0)
				{
					comp_ct=2;
					break;
				}
			}
		}
	}
	else if(!reduced_out)
	{
		comp_ct=num_sets;
	}
	free(L);
	uf_destroy(uf);
	free(uf);
	return comp_ct;
}
