/*
 * This file declares all global variables and functions. The existence of these are provided by "all.c". 
 */

#ifndef ALL_H
#define ALL_H

/*Libraries:*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

/*
Definitions:
*/

#define G 1.0 
#define dt 0.003125 /* Has to be a decimal number (write 14.0 instead of 14)*/
#define tmax 0.8 /* Has to be a decimal number*/
#define OPENING_ANGLE 0.5 /* Opening Angle used in the Barnes-Hut algorithm*/
#define SOFTENING 0.004 /* plummer-softening (of the gravitational potential). It has to be a decimal number. May not be 0 */
#define TYPICAL_MASS 0.000679348
#define INPUT_FILE_NAME "gassphere.ascii"
#define PRINT_EVERY_N_TIMESTEP 32


#define Nsph 1472/*Number of SPH-particles*/
#define Ndm 0/*Number of Dark Matter-particles*/
#define Nparticles (Ndm+Nsph) /*total number of particles*/

#define DESIRED_NUMBER_OF_NEIGHBORS	50
#define	NUMBER_OF_NEIGHBORS_TOLLERANCE	1
#define OUTPUT "./outputfiles/file" /* A number and ".txt" will be added to the filename. Make sure that the directory exists */
#define LOGFILE_NAME "logfile.log"
#define MAX_NR_OF_OUTPUT_FILES 100000
#define GAMMA  1.66666666666667 /* =5.0/3.0 Adiabatic index for SPH-particles*/
#define GAMMA_MINUS1  (GAMMA-1)
#define ALPHA 0.8 /* SPH-parameter - ArtBulkViscConst */
#define PI 3.14159265358979323846
#define SQRT3 1.732050807568877
#define SQRT2 1.41421356237310
#define A_HALF_TIMES_SQRT3 0.866025403784440 /* 0.5*sqrt(3) */
#define MAX_NUMBER_OF_DIVISIONS_IN_TREE_BUILD 200 /* Maximum number of divisions in "CreateNodes" with unchanged number of particles in the box in order to avoid errors when different particles have identical positions */

/* Structures */
extern struct Particle
{
	double Pos[3]; /* The position of the particle after the n'th timestep */
	double Vel[3]; /* The velocity of the particle after the (n+1/2)'th timestep */
	double PredictedVel[3]; /* The predicted velocity at the n'th timestep */
	double GravAcc[3]; /* The gravity-acceleration of the particle after the n'th timestep */
	double HydroAcc[3]; /* The hydro-acceleration of the particle after the n'th timestep */
	double f; /* eq. (8) in the Gadget2-paper. h*/
	double drhodh;
	double Mass;
	double GravPotential; 
	double Density;
	double Energy;
	double Entropy;/*Entropy Function defined as Pressure / Density^gamma */
	double PredictedEntropy;
	double DEntropyDt; /*Time-derivative of the entropy function*/
	double Pressure;
	double SmoothingLength; /*SPH smoothing length*/
	int *ListOfNeighbors; /* The ID's of the neighbor-particles width a distance less than the SmoothinLength */
	int NumberOfNeighbors; /* Number of particles in ListOfNeighbors[]: the particles with distance smaller than max(h_i,h_j) */
	int NumNgb2; /* The number of neighbors inside a sphere with a radius given by the SmoothingLength */
} Particles[Nparticles];

extern struct Node
{
	int NumberOfParticlesInBox;
	int *ParticleIDs;
	double Mass;
	double CM[3];
	double xMin[3];
	double LengthOfBox;
	double MaxSmoothingLength;
	struct Node *Children[8];
	unsigned short int EndNode; /* 1 if it is an EndNote. else 0.*/
} *Root;

/*File: */
extern FILE *LogFile;

/* Variables */
extern int FileCounter;
extern double t;

/* Functions */

/* File: density.c */
void DetermineSmoothingLength_Newton(int );
void CreateNeighborList(int , struct Node * );
void CreateNeighborList1(int , struct Node *, double );
int BoxIsTotallyOutsideSphere(int , struct Node * , double );
int BoxIsInsideSphere(int , struct Node *, double );
void InsertParticlesInBoxInNeighborList(int ,struct Node *);
void InsertParticlesInBoxInNeighborList1(int , struct Node *);
void UpdateSmoothing(int , int );
void DensityEvaluate(int , double);
double Constraint(double ,double );
void ResetNeighborList(int );

/* File: gravity.c */
void GravityAcceleration(int );
void ComputeAcc(int , struct Node *);
void ComputeAcc_DirectSummation(void);
void GetGravPotential(void);
void ComputeGravPotential(int, struct Node *);
double DistanceSquared(double *, double *);

/* File: initialize.c */
void Initialize(void);
void ImportData(void);

/* File: print.c */
void PrintToFile(int );
char* itoa (int );

/* File: run.c */
void Run(void);

/* File: sph.c*/
double dWdh(double , double );
double dmin(double , double );
double dmax(double , double );
void HydroAccAndEntropyDerivative(int );
double SPH_Kernel(double ,double );
double NablaKernel(int , int , double );
double SpeedOfSound(int );
void ComputeInitialPressure(void);
void ComputeInitialEntropy(void);
double ViscosityTensor(int , int );

/* File: tree.c*/
void BuildTree(void);
void CreateNodes(struct Node *, int, int );
void DestroyTree(struct Node *);
struct Node *AllocateMemoryForNode(void);
void CreateRoot(struct Node *);
struct Node *CreateChildren(struct Node *);
double ComputeMass(struct Node *);
double ComputeCM(int ,struct Node *);
void GetMaxSmoothingLength(struct Node *);

#endif
