#define PERIOD 10000
#define NMAX 20000
#define a 1.0
#define NXMAX 201
#define MAXFRAMES 20001
#define MAXRUN 201
#define MAXPARTICLETYPE 10
#define MAXBONDS 6 //Maximum bonds of a node

extern int N,Nb,Nd,bondGroup[NMAX*2],dihedralGroup[NMAX*4];
//N:#particles, Nb:#bonds, Nd:#dihedrals
extern float position[NMAX*3];
extern uint32_t particleID[NMAX];
extern char particleType[MAXPARTICLETYPE][2];

extern int NX,NY,RUNS,STEPS,LEN,FRAMES,JK_BIN_COUNT;
extern double KAPPA,EPSILON;

//Jack Knife blocking variables
extern double Xi[NXMAX][MAXRUN];
extern double Yi[NXMAX][MAXRUN];
extern double Zi[NXMAX][MAXRUN];
extern double Xj[NXMAX][MAXBONDS][MAXRUN];
extern double Yj[NXMAX][MAXBONDS][MAXRUN];
extern double Zj[NXMAX][MAXBONDS][MAXRUN];

// Ribbon Boundary Inplane Stretching Forces
extern double Fx[NXMAX][MAXRUN];
extern double Fy[NXMAX][MAXRUN];



