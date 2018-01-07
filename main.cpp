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
/*
double Xi[NXMAX][MAXRUN];
double Yi[NXMAX][MAXRUN];
double Zi[NXMAX][MAXRUN];
double Xj[NXMAX][MAXBONDS][MAXRUN];
double Yj[NXMAX][MAXBONDS][MAXRUN];
double Zj[NXMAX][MAXBONDS][MAXRUN];
*/

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
       print_and_exit("Usage: %s NX NY KAPPA STEPS\n",argv[0]);
   }
 
  sprintf(validrunfile,"../Sim_dump_stretched/L%d/W%d/k%.1f/valid_thermal.log",NX,NY,KAPPA);
  sprintf(numvalidrun_file,"../Sim_dump_stretched/L%d/W%d/k%.1f/numvalidruns.log",NX,NY,KAPPA);

  if(NULL==(vf=fopen(validrunfile,"r")))
	print_and_exit("I could not open file with simulation run numbers %s\n",validrunfile);

  if(NULL==(rf=fopen(numvalidrun_file,"r")))
        print_and_exit("I could not open file with number of valid simulation runs %s\n",numvalidrun_file);

  fscanf(rf, "%s", runnum);
  RUNS = atoi(runnum);
  printf("RUNS %d\n",RUNS);
  JK_BIN_COUNT=RUNS;
  EPSILON = 720.0 * KAPPA;
  //printf("NX = %d\n",NX);
  //printf("RUNS = %d\n",RUNS);

  FRAMES=STEPS/PERIOD;

//  initialize1();
//  initialize2();

  // Init_strip.gsd filepath
  sprintf(init_strip,"../Sim_dump_ribbon/init_strip_L%d_W%d.gsd",NX,NY);
  printf("Init_strip.gsd : %s\n",init_strip);
  sprintf(observable_file,"../Sim_dump_stretched/L%d/W%d/k%.1f/r",NX,NY,KAPPA);
  printf("observable_file: %s\n",observable_file);

  load_gsd(init_strip,0); 
  //printf("NX = %d\n",NX);


  //Data Blocking
  int iblo,r=0,frames,imes,col=0;
  char runnum_str[20];
  double Rij;
  //Building the neghbor node list
  for(iblo=0; iblo<NY; iblo++)
  {
    neigbhor(col+iblo*NX,col);
  }
  for(iblo=0; iblo<NY; iblo++)
  {
    printf("%d    ",col+iblo*NX);
    for(int j=0;j<numneigbhorbonds(col+iblo*NX);j++)
    {
      printf("%d ",neighbor_node[iblo][j]);
    }
    printf("\n");
  }

  initialize1();

  while (fscanf(vf, "%s", runnum_str) == 1)// 1 is returned if fscanf reads a valid filepath
  {
          //printf("%d\n",r);
	  // Trajectory.gsd filepath
          //printf("%s\n",observable_file);
	  sprintf(trajectory_file,"%s%d/traj_stretched.gsd",observable_file,atoi(runnum_str));
	  //printf("Trajectory file being read : %s\n",trajectory_file);
		
	  for(iblo=0; iblo<NY; iblo++)
  	  {
		  //Looping through frames
		  for(frames=FRAMES/2;frames<FRAMES;frames++)
		  {
			load_gsd(trajectory_file,frames);
		        for(int j=0;j<numneigbhorbonds(iblo);j++)
		        {
			  Rij = pow((position[3*(col+iblo*NX)]- position[3*neighbor_node[iblo][j]]),2);
			  Rij += pow((position[3*(col+iblo*NX)+1]-position[(3*neighbor_node[iblo][j])+1]),2);
			  Rij += pow((position[3*(col+iblo*NX)+2]-position[(3*neighbor_node[iblo][j])+2]),2);
			  Rij = sqrt(Rij);
                          //printf("%.8f\n",Rij);
			  //if(imes == 0)
				//printf("Nodei %d Rij %.8g %.8f %.8f %.8f %.8f %.8f\n",iblo,Rij,(Xi[iblo][imes]-Xj[iblo][j][imes]),(Yi[iblo][imes]-Yj[iblo][j][imes]),(Zi[iblo][imes]-Zj[iblo][j][imes]),Xi[iblo][imes],Xj[iblo][j][imes]);
			  Fx[iblo][r]+= -EPSILON * (Rij-a)*(position[3*(col+iblo*NX)]- position[3*neighbor_node[iblo][j]])/Rij;
			  Fy[iblo][r]+= -EPSILON * (Rij-a)*(position[3*(col+iblo*NX)+1]-position[(3*neighbor_node[iblo][j])+1])/Rij;
		        }
		  }
                  Fx[iblo][r]/=FRAMES/2; //Averaging the force by frames at each node
                  Fy[iblo][r]/=FRAMES/2;
                  if(iblo==3)
                    printf("Node %d Run %d Fx %.8f Fy %.8f\n",col+iblo*NX,r,Fx[iblo][r],Fy[iblo][r]);

                  //Sum of all the avg forces at each node for all the runs
                  Fx[iblo][RUNS] += Fx[iblo][r];
                  Fy[iblo][RUNS] += Fy[iblo][r];
                  //printf("Node Num %d  %.8f   %.8f\n",iblo,Fx[iblo][RUNS],Fy[iblo][RUNS]);
	  }
	  r++;//counter through the runs
   }
   fclose(vf);

  for(iblo=0; iblo<NY; iblo++)//counter over the ribbon edge nodes
  {
        //printf("Node Num %d  %.8f   %.8f\n",iblo,Fx[iblo][RUNS],Fy[iblo][RUNS]);
	//Jack Knife Blocking
	for(imes=0;imes<RUNS;imes++)
        {
		Fx[iblo][imes]=(Fx[iblo][RUNS]-Fx[iblo][imes])/(RUNS-1.);
                Fy[iblo][imes]=(Fy[iblo][RUNS]-Fy[iblo][imes])/(RUNS-1.);
	}
   }
	
   //Jack knife error
   double jk_error_Fx[NXMAX],jk_error_term1_Fx[NXMAX],jk_error_term2_Fx[NXMAX];
   double jk_error_Fy[NXMAX],jk_error_term1_Fy[NXMAX],jk_error_term2_Fy[NXMAX];
   //Jack Knife Error 
    for(int i=0;i<NX;i++)
    {
        jk_error_Fx[i]=0;
        jk_error_term1_Fx[i]=0;
        jk_error_term2_Fx[i]=0;
        jk_error_Fy[i]=0;
        jk_error_term1_Fy[i]=0;
        jk_error_term2_Fy[i]=0;
    }

    //JK_BIN_COUNT=RUNS;
    for(int i=0;i<NY;i++)
    {
        for(int j=0;j<JK_BIN_COUNT;j++)
        {
                jk_error_term1_Fy[i] += Fy[i][j] * Fy[i][j];
                jk_error_term1_Fx[i] += Fx[i][j] * Fx[i][j];
        }
        jk_error_term1_Fy[i] = (1.0/JK_BIN_COUNT) * jk_error_term1_Fy[i];
        jk_error_term1_Fx[i] = (1.0/JK_BIN_COUNT) * jk_error_term1_Fx[i];
        
	for(int j=0;j<JK_BIN_COUNT;j++)
        {
                jk_error_term2_Fy[i] += (1.0/JK_BIN_COUNT) * Fy[i][j];
                jk_error_term2_Fx[i] += (1.0/JK_BIN_COUNT) * Fx[i][j];
        }
        jk_error_term2_Fy[i] = jk_error_term2_Fy[i] * jk_error_term2_Fy[i];
        jk_error_term2_Fx[i] = jk_error_term2_Fx[i] * jk_error_term2_Fx[i];

        //      JK Error        
        jk_error_Fy[i] = sqrt((JK_BIN_COUNT-1)*(jk_error_term1_Fy[i] - jk_error_term2_Fy[i]));
        jk_error_Fx[i] = sqrt((JK_BIN_COUNT-1)*(jk_error_term1_Fx[i] - jk_error_term2_Fx[i]));

        printf ("%d\t%.8g\t%.8g\t%.8g\t%.8g\n",i,Fy[i][RUNS]/RUNS,jk_error_Fy[i],Fx[i][RUNS]/RUNS,jk_error_Fx[i]);
  }
  return 0;
}
