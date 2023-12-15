/*
 * This file defines all global variables and functions.
 */

/*Libraries:*/
#include "all.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

struct Particle Particles[Nparticles	];

struct Node *Root;

/*File: */
FILE *LogFile;

/* Variables */
int FileCounter;
double t;


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
