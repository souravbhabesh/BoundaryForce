#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <math.h>
#include "stdint.h"
#include "gsd.h"
#include "variables.h"
#include "gsd_read.h"
#include "force.h"


int NX,NY,LEN,RUNS,JK_BIN_COUNT;
double KAPPA,EPSILON;
int STEPS,FRAMES;
//Jack Knife blocking variables
double Xi[NXMAX][MAXRUN];
double Yi[NXMAX][MAXRUN];
double Zi[NXMAX][MAXRUN];
double Xj[NXMAX][MAXBONDS][MAXRUN];
double Yj[NXMAX][MAXBONDS][MAXRUN];
double Zj[NXMAX][MAXBONDS][MAXRUN];

// Ribbon Boundary Inplane Stretching Forces
double Fx[NXMAX][MAXRUN];
double Fy[NXMAX][MAXRUN];

int main(int argc, char **argv)
{
  FILE *pf,*vf,*rf;
  char init_strip[256],trajectory_file[256],validrunfile[256],observable_file[1024],numvalidrun_file[256],runnum[256];
  int frame_cnt=0;
  double NORM;

   switch (argc){
     case 5:
       sscanf(argv[1],"%d",&NX);    
       sscanf(argv[2],"%d",&NY);
       sscanf(argv[3],"%lf",&KAPPA);
       sscanf(argv[4],"%d",&STEPS); 
       break;
     default:
       print_and_exit("Usage: %s NX NY KAPPA STEPS\n");
   }
 
  sprintf(validrunfile,"../Sim_dump_ribbon/L%d/W%d/k%.1f/validruns.log",NX,NY,KAPPA);
  sprintf(numvalidrun_file,"../Sim_dump_ribbon/L%d/W%d/k%.1f/numvalidruns.log",NX,NY,KAPPA);

  if(NULL==(vf=fopen(validrunfile,"r")))
	print_and_exit("I could not open file with simulation run numbers %s\n",validrunfile);

  if(NULL==(rf=fopen(numvalidrun_file,"r")))
        print_and_exit("I could not open file with simulation run numbers %s\n",numvalidrun_file);

  fscanf(rf, "%s", runnum);
  RUNS = atoi(runnum);
  JK_BIN_COUNT=RUNS;
  EPSILON = 720.0 * KAPPA;
  //printf("NX = %d\n",NX);
  //printf("RUNS = %d\n",RUNS);

  FRAMES=STEPS/PERIOD;

  initialize1();
  initialize2();

  // Init_strip.gsd filepath
  sprintf(init_strip,"../Sim_dump_ribbon/init_strip_L%d_W%d.gsd",NX,NY);
  //printf("Init_strip.gsd : %s\n",init_strip);

  load_gsd(init_strip,0); 
  //printf("NX = %d\n",NX);

  //Data Blocking
  int iblo,r=0,frames,imes;
  double Rij;

  while (fscanf(vf, "%s", observable_file) == 1)// 1 is returned if fscanf reads a valid filepath
  {
	  // Trajectory.gsd filepath
	  sprintf(trajectory_file,"%s/traj_thermal.gsd",observable_file);
	  //printf("Trajectory file being read : %s\n",trajectory_file);
		
	  for(iblo=0; iblo<NX; iblo++)
  	  {
		  //Looping through frames
		  for(frames=FRAMES/2;frames<FRAMES;frames++)
		  {
			load_gsd(trajectory_file,frames);
			//Data Blocking for node i
			Xi[iblo][r]+=position[3*iblo];
			Yi[iblo][r]+=position[3*iblo+1];
			Zi[iblo][r]+=position[3*iblo+2];
			//Data blocking for node j		
			neigbhor(iblo,r,iblo,frames);
		  }
	  }
	  r++;//counter through the runs
   }
   fclose(vf);

  for(iblo=0; iblo<NX; iblo++)//counter over the ribbon edge nodes
  {
	//printf("iblo %d\n",iblo);
	for(imes=0;imes<RUNS;imes++)
        {
		Xi[iblo][imes]/=(FRAMES/2);
                Yi[iblo][imes]/=(FRAMES/2);
                Zi[iblo][imes]/=(FRAMES/2);
        }
	normalizeNeighbor(iblo);
	for(imes=0;imes<RUNS;imes++)
        {
		Xi[iblo][RUNS]+=Xi[iblo][imes];
		Yi[iblo][RUNS]+=Yi[iblo][imes];
		Zi[iblo][RUNS]+=Zi[iblo][imes];
	}
	sumNeighbor(iblo);
	
	//Jack Knife Blocking
	for(imes=0;imes<RUNS;imes++)
        {
		Xi[iblo][imes]=(Xi[iblo][RUNS]-Xi[iblo][imes])/(RUNS-1.);
		Yi[iblo][imes]=(Yi[iblo][RUNS]-Yi[iblo][imes])/(RUNS-1.);
		Zi[iblo][imes]=(Zi[iblo][RUNS]-Zi[iblo][imes])/(RUNS-1.);
		//if(imes == 0)
			//printf("Nodei %d run %d %.8f %.8f %.8f\n",iblo,imes,Xi[iblo][imes],Yi[iblo][imes],Zi[iblo][imes]);
	}
	jackknifeNeighbor(iblo);

	//Evaluate Force acting on boundary nodes using Jack Knife blocked data
	for(imes=0;imes<RUNS;imes++)
        {
		for(int j=0;j<numneigbhorbonds(iblo);j++)
		{
			Rij = pow((Xi[iblo][imes]-Xj[iblo][j][imes]),2);
			Rij += pow((Yi[iblo][imes]-Yj[iblo][j][imes]),2);
			Rij += pow((Zi[iblo][imes]-Zj[iblo][j][imes]),2);
			Rij = sqrt(Rij);
			//if(imes == 0)
				//printf("Nodei %d Rij %.8g %.8f %.8f %.8f %.8f %.8f\n",iblo,Rij,(Xi[iblo][imes]-Xj[iblo][j][imes]),(Yi[iblo][imes]-Yj[iblo][j][imes]),(Zi[iblo][imes]-Zj[iblo][j][imes]),Xi[iblo][imes],Xj[iblo][j][imes]);
			Fx[iblo][imes]+= -EPSILON * (Rij-a)*(Xi[iblo][imes]-Xj[iblo][j][imes])/Rij;
			Fy[iblo][imes]+= -EPSILON * (Rij-a)*(Yi[iblo][imes]-Yj[iblo][j][imes])/Rij;
		}
  	}

	Fx[iblo][RUNS]=0;
	Fy[iblo][RUNS]=0;
	for(imes=0;imes<RUNS;imes++)
        {
		Fx[iblo][RUNS]+=Fx[iblo][imes];
		Fy[iblo][RUNS]+=Fy[iblo][imes];
	}
	Fx[iblo][RUNS]/=(RUNS*1.0);
	Fy[iblo][RUNS]/=(RUNS*1.0);
   }
	
   //Normalizing using the left clamped width
   //NORM = rwidth[0][RUNS];
   //for(iblo=0; iblo<NX; iblo++)//counter over the ribbon edge nodes
   //{
//	rwidth[iblo][RUNS]/=NORM;
	//printf("%d\t%.8f\n",iblo,rwidth[iblo][RUNS]);
  // }  

   //Jack knife error
       double jk_error[NXMAX],jk_error_term1[NXMAX],jk_error_term2[NXMAX];
    //          Jack Knife Error        
    for(int i=0;i<NX;i++)
    {
        jk_error[i]=0;
        jk_error_term1[i]=0;
        jk_error_term2[i]=0;
    }

    JK_BIN_COUNT=RUNS;
    for(int i=0;i<NX;i++)
    {
        for(int j=0;j<JK_BIN_COUNT;j++)
        {
                jk_error_term1[i] += Fy[i][j] * Fy[i][j];
        }
        jk_error_term1[i] = (1.0/JK_BIN_COUNT) * jk_error_term1[i];
        
	for(int j=0;j<JK_BIN_COUNT;j++)
        {
                jk_error_term2[i] += (1.0/JK_BIN_COUNT) * Fy[i][j];
        }
        jk_error_term2[i] = jk_error_term2[i] * jk_error_term2[i];
        //      JK Error        
        jk_error[i] = sqrt((JK_BIN_COUNT-1)*(jk_error_term1[i] - jk_error_term2[i]));
        printf ("%d\t%.8g\t%.8g\n",i,Fy[i][RUNS],jk_error[i]);
    }


	  
  return 0;
}
