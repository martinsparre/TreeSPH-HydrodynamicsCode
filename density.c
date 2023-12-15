#include "all.h"


void CreateNeighborList(int aParticleID, struct Node *A)
{
	int i;
	double h = dmax(A->MaxSmoothingLength,Particles[aParticleID].SmoothingLength);

	if ( A->NumberOfParticlesInBox == 0 )
	{
		/* do nothing */
	}
	else if(A->EndNode == 1)
	{
		if(DistanceSquared(Particles[aParticleID].Pos ,  Particles[ A->ParticleIDs[0] ].Pos ) <= h*h)
		InsertParticlesInBoxInNeighborList(aParticleID, A);
	}
	else if( BoxIsTotallyOutsideSphere(aParticleID, A , h ) )
	{
		/* do nothing */
	}
	else if( BoxIsInsideSphere(aParticleID, A , h  ) )
	{
		InsertParticlesInBoxInNeighborList(aParticleID, A);
	}

	else
	{
		for(i=0;i<8;i++)
		{
			CreateNeighborList(aParticleID,A->Children[i]);	
		}
	}
}

void CreateNeighborList1(int aParticleID, struct Node *A, double h)
{
	int i;

	if ( A->NumberOfParticlesInBox == 0 )
	{
		/* do nothing */
	}
	else if(A->EndNode == 1)
	{
		if(DistanceSquared(Particles[aParticleID].Pos ,  Particles[ A->ParticleIDs[0] ].Pos ) <= h*h)
		InsertParticlesInBoxInNeighborList1(aParticleID, A);
	}
	else if( BoxIsTotallyOutsideSphere(aParticleID, A , h ) )
	{
		/* do nothing */
	}
	else if( BoxIsInsideSphere(aParticleID, A , h  ) )
	{
		InsertParticlesInBoxInNeighborList1(aParticleID, A);
	}

	else
	{
		for(i=0;i<8;i++)
		{
			CreateNeighborList1(aParticleID,A->Children[i],h);	
		}
	}
}


int BoxIsTotallyOutsideSphere(int aParticleID, struct Node *A , double h)
{
	if(

	  pow(   Particles[aParticleID].Pos[0] - (A->xMin[0] + 0.5 * A->LengthOfBox)   ,2)
	+ pow(   Particles[aParticleID].Pos[1] - (A->xMin[1] + 0.5 * A->LengthOfBox)   ,2)
	+ pow(   Particles[aParticleID].Pos[2] - (A->xMin[2] + 0.5 * A->LengthOfBox)   ,2)
	
	 > pow(h + A_HALF_TIMES_SQRT3 * A->LengthOfBox,2)

	)
	return 1;
	else
	return 0;
}

int BoxIsInsideSphere(int aParticleID, struct Node *A , double h)
{
	int T=0;
	double temp = h - A_HALF_TIMES_SQRT3*A->LengthOfBox;

	if(	temp > 0.0   )
	{
		if(
		
		  pow(  Particles[aParticleID].Pos[0] - (A->xMin[0] + 0.5 * A->LengthOfBox)   ,2)
		+ pow(  Particles[aParticleID].Pos[1] - (A->xMin[1] + 0.5 * A->LengthOfBox)   ,2)
		+ pow(  Particles[aParticleID].Pos[2] - (A->xMin[2] + 0.5 * A->LengthOfBox)   ,2)
		
		 < pow(temp,2)	
		 
		)
		T=1;
	}
	
	return T;
}


void InsertParticlesInBoxInNeighborList(int aParticleID, struct Node *A)
{
	int i, NewParticle = 0;
	double r2;

	/* This statement allocates memory for the new particles in the Neighbor list */
	
	for(i=0; i <  A->NumberOfParticlesInBox; i++)
	{
		r2 = DistanceSquared( Particles[A->ParticleIDs[i]].Pos , Particles[aParticleID].Pos );

		/* A particle is inserted in the neighbor-list */

		if(r2 < pow(Particles[aParticleID].SmoothingLength,2) && A->ParticleIDs[i] < Nsph)
		{
			Particles[aParticleID].NumNgb2 += 1;
			Particles[aParticleID].NumberOfNeighbors += 1;
			NewParticle = 1;
		}
		else if( r2 < pow(Particles[A->ParticleIDs[i]].SmoothingLength,2)  && A->ParticleIDs[i] < Nsph)
		{
			Particles[aParticleID].NumberOfNeighbors += 1;
			NewParticle = 1;
		}

		if(NewParticle)
		{
			Particles[aParticleID].ListOfNeighbors = realloc(Particles[aParticleID].ListOfNeighbors, (Particles[aParticleID].NumberOfNeighbors ) * sizeof(int));

			Particles[aParticleID].ListOfNeighbors[  Particles[aParticleID].NumberOfNeighbors-1] = A->ParticleIDs[i];

			NewParticle = 0;
		}
	}
}

void InsertParticlesInBoxInNeighborList1(int aParticleID, struct Node *A)
{
	int i;
	double r2;

	for(i=0; i <  A->NumberOfParticlesInBox; i++)
	{
		r2 = DistanceSquared( Particles[A->ParticleIDs[i]].Pos , Particles[aParticleID].Pos );

		/* A particle is inserted in the neighbor-list */

		if(r2 < pow(Particles[aParticleID].SmoothingLength,2)  && A->ParticleIDs[i] < Nsph)
		{
			Particles[aParticleID].NumNgb2 += 1;
			Particles[aParticleID].NumberOfNeighbors += 1;

			Particles[aParticleID].ListOfNeighbors = realloc(Particles[aParticleID].ListOfNeighbors, (Particles[aParticleID].NumberOfNeighbors ) * sizeof(int));

			Particles[aParticleID].ListOfNeighbors[  Particles[aParticleID].NumberOfNeighbors-1] = A->ParticleIDs[i];
		}
	}
}	

/*
void DetermineSmoothingLength_Bisection(int i)
{
	double CM,CL,CR,New;

	if(t!=0.0)
	{		
		Particles[i].Left = Particles[i].SmoothingLength/1.1;
		Particles[i].Right = 1.1 * Particles[i].SmoothingLength;
	}
	else
	{
		Particles[i].Left = Particles[i].SmoothingLength/2.;
		Particles[i].Right = 2.0 * Particles[i].SmoothingLength;			
	}
	
	ResetNeighborList(i);
	CreateNeighborList1(i,Root,Particles[i].Left);
	DensityEvaluate(i, Particles[i].Left);
	
	CL = Constraint(Particles[i].Left,Particles[i].Density);
		
	while(CL > 0)
	{
		Particles[i].Left /= 1.26;
			
		ResetNeighborList(i);
		CreateNeighborList1(i,Root,Particles[i].Left);
		
		DensityEvaluate(i, Particles[i].Left);
		CL = Constraint(Particles[i].Left,Particles[i].Density);
	}

	ResetNeighborList(i);
	CreateNeighborList1(i,Root,Particles[i].Right);
		
	DensityEvaluate(i, Particles[i].Right);
	CR = Constraint(Particles[i].Right,Particles[i].Density);
	while(CR < 0)
	{
		Particles[i].Right *= 1.26;
		ResetNeighborList(i);
		CreateNeighborList1(i,Root,Particles[i].Right);			
			
		DensityEvaluate(i, Particles[i].Right);
		CR = Constraint(Particles[i].Right,Particles[i].Density);
	}
		
	New = Particles[i].Right - CR / (   Particles[i].drhodh*pow(Particles[i].Right,3)+3.0*pow(Particles[i].Right,2)*Particles[i].Density );

	if(New>Particles[i].Left || New < Particles[i].Right)
	Particles[i].SmoothingLength = New;
	else
	Particles[i].SmoothingLength = 0.5 * (Particles[i].Left + Particles[i].Right);
		
	ResetNeighborList(i);
	CreateNeighborList1(i,Root,Particles[i].SmoothingLength);
		
	DensityEvaluate(i, Particles[i].SmoothingLength);
	CM = Constraint(Particles[i].SmoothingLength,Particles[i].Density);

		if(CM>0)
		{
			Particles[i].Right = Particles[i].SmoothingLength;
		}
		else if(CM<0)
		{
			Particles[i].Left = Particles[i].SmoothingLength;
		}
		
	while(fabs( (Particles[i].Right - Particles[i].Left) / (Particles[i].SmoothingLength) ) > 0.01 )
//	while(fabs( CM ) > DESIRED_NUMBER_OF_NEIGHBORS*TYPICAL_MASS*3./(4.*PI)*0.01 )
	{
		Particles[i].SmoothingLength = 0.5 * (Particles[i].Left + Particles[i].Right);

		ResetNeighborList(i);
		CreateNeighborList1(i,Root,Particles[i].SmoothingLength);

		DensityEvaluate(i, Particles[i].SmoothingLength);
		CM = Constraint(Particles[i].SmoothingLength,Particles[i].Density);
		
		
		if(CM>0)
		{
			Particles[i].Right = Particles[i].SmoothingLength;
		}
		else if(CM<0)
		{
			Particles[i].Left = Particles[i].SmoothingLength;
		}
		else
		{
		//	printf("CM=0.0000\n");
		break;
		}	
	}
}
*/

void DetermineSmoothingLength_Newton(int i)
{
		double delta=Particles[i].SmoothingLength;
//		int N=0;
		while(delta > 0.05*Particles[i].SmoothingLength )
//		while(N<3)
		{
			ResetNeighborList(i);
			CreateNeighborList1(i,Root,Particles[i].SmoothingLength);
			DensityEvaluate(i, Particles[i].SmoothingLength);
			
			delta = Constraint(Particles[i].SmoothingLength,Particles[i].Density)
			/ ( Particles[i].drhodh*pow(Particles[i].SmoothingLength,3)  + Particles[i].Density*3*pow(Particles[i].SmoothingLength,2) );
			
			Particles[i].SmoothingLength -= delta ;
//			N++;
		}
}


double Constraint(double h,double rho)
{
	return h*h*h * rho - DESIRED_NUMBER_OF_NEIGHBORS*TYPICAL_MASS*3./(4.*PI);
}

void DensityEvaluate(int i, double h)
{
	int j,k;
	double R;
	double drdh=0.0,dens=0.0;

	for(j=0;j< Particles[i].NumberOfNeighbors; j++)
	{
		k=Particles[i].ListOfNeighbors[j];
		
		R = sqrt(
			  pow( Particles[i].Pos[0] - Particles[k].Pos[0] ,2)
			+ pow( Particles[i].Pos[1] - Particles[k].Pos[1] ,2)
			+ pow( Particles[i].Pos[2] - Particles[k].Pos[2] ,2)
			);

		dens += Particles[ k ].Mass 
		* ( SPH_Kernel(R,h) );
		
		drdh += Particles[ k ].Mass * dWdh(R,h);
	}
	
	Particles[i].f = 1./(1.+ drdh*h / (3.*dens));
	Particles[i].Density = dens;
	Particles[i].drhodh = drdh;
}


void ResetNeighborList(int i)
{
		free(Particles[i].ListOfNeighbors);
		Particles[i].NumberOfNeighbors = 0;
		Particles[i].NumNgb2 = 0;
		Particles[i].ListOfNeighbors = NULL;
}


void UpdateSmoothing(int i, int N)
{
	
	if(Particles[i].NumNgb2 > DESIRED_NUMBER_OF_NEIGHBORS + NUMBER_OF_NEIGHBORS_TOLLERANCE 
	|| Particles[i].NumNgb2 < DESIRED_NUMBER_OF_NEIGHBORS - NUMBER_OF_NEIGHBORS_TOLLERANCE)
	{
		Particles[i].SmoothingLength = Particles[i].SmoothingLength * 0.5 
		*(1.0 +   pow(       ( (double)DESIRED_NUMBER_OF_NEIGHBORS) / ( ((double)Particles[i].NumNgb2 ) + 0.01*DESIRED_NUMBER_OF_NEIGHBORS ), 0.33333333333333 )       );
		
		Particles[i].NumberOfNeighbors = 0;
		Particles[i].NumNgb2 = 0;
		free(Particles[i].ListOfNeighbors);
		Particles[i].ListOfNeighbors = NULL;
		
		CreateNeighborList(i, Root);
		
		if(N<15)
		{
			UpdateSmoothing( i, N+1);
		}	
	}
}

