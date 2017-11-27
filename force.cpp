#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <stdlib.h>
#include "stdint.h"
#include "variables.h"
#include "force.h"

int initialize1()
{
   for(int i=0;i<=NX;i++)
   {
	for(int j=0;j<=RUNS;j++)
   	{
	   Xi[i][j]=0;
	   Yi[i][j]=0;
	   Zi[i][j]=0;
	   Fx[i][j]=0;
	   Fy[i][j]=0;
	}
   }
   return 0;
}

int initialize2()
{
   for(int i=0;i<=NX;i++)
   {
	for(int j=0;j<MAXBONDS;j++)
        {
		for(int k=0;k<=RUNS;k++)
        	{
			Xj[i][j][k]=0;
			Yj[i][j][k]=0;
			Zj[i][j][k]=0;
		}
	}
    
   }
   return 0;
}

int numneigbhorbonds(int node)
{
 int i,numbonds=0;
 for(i=0;i<Nb;i++)
 {
 	if(bondGroup[i*2]==node || bondGroup[i*2+1]==node)
		numbonds++;
 }
 return numbonds;
}

int neigbhor(int node, int run, int boundarynodenum,int frames)
{
  int i,numbonds=0,nodej;
  for(i=0;i<Nb;i++)
  {
	if(bondGroup[i*2]==node)
		nodej=bondGroup[i*2+1];
	else if(bondGroup[i*2+1]==node)
		nodej=bondGroup[i*2];
	else
		nodej=-1;//no match
	if(nodej != -1)
	{
		Xj[boundarynodenum][numbonds][run]+=position[3*nodej];
		Yj[boundarynodenum][numbonds][run]+=position[3*nodej+1];
		Zj[boundarynodenum][numbonds][run]+=position[3*nodej+2];
		//if(run == 1 && frames==FRAMES/2)
               		//printf("Nodei %d Node j %d #neigbors %d\n",node,nodej,numbonds);
		numbonds++;
		nodej=-1;
	}
  }
  return 0;
}

int normalizeNeighbor(int iblo)
{
  int j,k;
	for(j=0;j<numneigbhorbonds(iblo);j++)
	{
		for(k=0;k<RUNS;k++)
		{
			Xj[iblo][j][k]/=(FRAMES/2);
			Yj[iblo][j][k]/=(FRAMES/2);
			Zj[iblo][j][k]/=(FRAMES/2);
		}
	}
  return 0;
}

int sumNeighbor(int iblo)
{
  int j,k;
        for(j=0;j<numneigbhorbonds(iblo);j++)
        {
                for(k=0;k<RUNS;k++)
                {
                        Xj[iblo][j][RUNS]+=Xj[iblo][j][k];
                        Yj[iblo][j][RUNS]+=Yj[iblo][j][k];
                        Zj[iblo][j][RUNS]+=Zj[iblo][j][k];
                }
        }
  return 0;
}

int jackknifeNeighbor(int iblo)
{
  int j,k;
        for(j=0;j<numneigbhorbonds(iblo);j++)
        {
                for(k=0;k<RUNS;k++)
                {
                        Xj[iblo][j][k]=(Xj[iblo][j][RUNS]-Xj[iblo][j][k])/(RUNS-1.);
                        Yj[iblo][j][k]=(Yj[iblo][j][RUNS]-Yj[iblo][j][k])/(RUNS-1.);
                        Zj[iblo][j][k]=(Zj[iblo][j][RUNS]-Zj[iblo][j][k])/(RUNS-1.);
			//if(k==0)
				//printf("Nodei %d Bondj# %d Run %d %.8f %.8f %.8f\n",iblo,j,k,Xj[iblo][j][k],Yj[iblo][j][k],Zj[iblo][j][k]);
                }
        }
  return 0;
}
