#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <stdlib.h>
#include "stdint.h"
#include "variables.h"
#include "force.h"


int neighbor_node[NXMAX][6];

int initialize1()
{
   for(int i=0;i<=NY;i++)
   {
	for(int j=0;j<=RUNS;j++)
   	{
	   Fx[i][j]=0;
	   Fy[i][j]=0;
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

int neigbhor(int node,int col)
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
          neighbor_node[(node-col)/NX][numbonds]=nodej;
	  numbonds++;
	  nodej=-1;
	}
  }
  return 0;
}

