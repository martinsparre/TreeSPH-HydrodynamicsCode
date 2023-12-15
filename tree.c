#include "all.h"

void BuildTree(void)
{
	Root = AllocateMemoryForNode();

	CreateRoot(Root);

	CreateNodes(Root,0,Nparticles+1);

	ComputeMass(Root);

	#pragma omp parallel
	{
	#pragma omp sections
		{
			#pragma omp section
			ComputeCM(0,Root);
			#pragma omp section
			ComputeCM(1,Root);
			#pragma omp section
			ComputeCM(2,Root);
		}
	}
}


void CreateNodes(struct Node *A, int Divisions, int N) /* Divisions : number of divisions with unchanged number of particles in box. N: number of particles in parent node */
{
	int i;
		if(A->NumberOfParticlesInBox > 1 && Divisions < MAX_NUMBER_OF_DIVISIONS_IN_TREE_BUILD)
		{
			A=CreateChildren(A);
		
			if(A==Root)
			{
				#pragma omp parallel for private(i)
				for(i=0;i<8;i++)
				CreateNodes(A->Children[i], 0,A->NumberOfParticlesInBox);
			}
			else if(A->NumberOfParticlesInBox < N)
			{
				for(i=0;i<8;i++)
				CreateNodes(A->Children[i], 0,A->NumberOfParticlesInBox);
			}
			else if (A->NumberOfParticlesInBox == N)
			{
				for(i=0;i<8;i++)
				CreateNodes(A->Children[i], Divisions+1,A->NumberOfParticlesInBox);
			}

		}
		else
		{
			A->EndNode = 1;
		}
}

void DestroyTree(struct Node *A)
{
	int j;

	if(A != NULL)
	{

		for(j=0;j<8;j++)
				DestroyTree(A->Children[j]);
		
		free(A->ParticleIDs);
		free(A);
	}
}

struct Node *AllocateMemoryForNode(void)
{
	return (struct Node *) malloc(sizeof(struct Node));	
}

void CreateRoot(struct Node *R)
{
	int i;

	double xMin,yMin,zMin,xMax,yMax,zMax, LengthOfBox;
	
	if(R == NULL)
	{
		printf("ERROR: Not enough RAM... Error: 32\n");
	    fprintf(LogFile,"ERROR: Not enough RAM... Error: 32\n");
	}
	
	R->NumberOfParticlesInBox = Nparticles;

	xMin = Particles[0].Pos[0];
	xMax = Particles[0].Pos[0];
	yMin = Particles[0].Pos[1];
	yMax = Particles[0].Pos[1];
	zMin = Particles[0].Pos[2];
	zMax = Particles[0].Pos[2];

	R->ParticleIDs = (int *) malloc(Nparticles*sizeof(int));
	
	for(i=0;i<Nparticles;i++)
	{
		R->ParticleIDs[i] = i;

		if(Particles[i].Pos[0] < xMin)
		{
			xMin = Particles[i].Pos[0];
		}

		if(Particles[i].Pos[1] < yMin)
		{
			yMin = Particles[i].Pos[1];
		}

		if(Particles[i].Pos[2] < zMin)
		{
			zMin = Particles[i].Pos[2];
		}

		if(Particles[i].Pos[0] > xMax)
		{
			xMax = Particles[i].Pos[0];
		}

		if(Particles[i].Pos[1] > yMax)
		{
			yMax = Particles[i].Pos[1];
		}

		if(Particles[i].Pos[2] > zMax)
		{
			zMax = Particles[i].Pos[2];
		}
	}

	LengthOfBox = xMax - xMin;
	
	if( (yMax-yMin) > LengthOfBox )
	{
			LengthOfBox = yMax - yMin;
	}
	
	if( (zMax-zMin) > LengthOfBox )
	{
			LengthOfBox = zMax - zMin;
	}

	R->xMin[0] = xMin;
	R->xMin[1] = yMin;
	R->xMin[2] = zMin;

	R->LengthOfBox = LengthOfBox;
	R->EndNode = 0;


	for(i=0;i<8;i++)
	R->Children[i] = NULL;

}

struct Node *CreateChildren(struct Node *ParentNode)
{
	int j,w,q, i, nX, nY, nZ, BoxNumber;
	
	for(j=0;j<8;j++)
	ParentNode->Children[j] = AllocateMemoryForNode();


	if(ParentNode->Children[7] == NULL)
	{
		printf("ERROR: Not enough RAM...\n");
		fprintf(LogFile,"ERROR: Not enough RAM...\n");
	}
	
	for(w=0;w<8;w++)
	{
		ParentNode->Children[w]->ParticleIDs = (int *) malloc(sizeof(int));

		for(q=0;q<8;q++)
		{
			ParentNode->Children[w]->Children[q] = NULL;
		}
	}

	for(j=0;j<8;j++)
	{
		ParentNode->Children[j]->LengthOfBox = ParentNode->LengthOfBox / 2.0;
		ParentNode->Children[j]->NumberOfParticlesInBox = 0;
		ParentNode->Children[j]->EndNode = 0;		
	}		
	
		ParentNode->Children[0]->xMin[0] = ParentNode->xMin[0];
		ParentNode->Children[0]->xMin[1] = ParentNode->xMin[1];
		ParentNode->Children[0]->xMin[2] = ParentNode->xMin[2];

		ParentNode->Children[1]->xMin[0] = ParentNode->xMin[0] + ParentNode->Children[1]->LengthOfBox;
		ParentNode->Children[1]->xMin[1] = ParentNode->xMin[1];
		ParentNode->Children[1]->xMin[2] = ParentNode->xMin[2];

		ParentNode->Children[2]->xMin[0] = ParentNode->xMin[0];
		ParentNode->Children[2]->xMin[1] = ParentNode->xMin[1] + ParentNode->Children[2]->LengthOfBox;
		ParentNode->Children[2]->xMin[2] = ParentNode->xMin[2];

		ParentNode->Children[3]->xMin[0] = ParentNode->xMin[0];
		ParentNode->Children[3]->xMin[1] = ParentNode->xMin[1];
		ParentNode->Children[3]->xMin[2] = ParentNode->xMin[2] + ParentNode->Children[3]->LengthOfBox;

		ParentNode->Children[4]->xMin[0] = ParentNode->xMin[0] + ParentNode->Children[4]->LengthOfBox;
		ParentNode->Children[4]->xMin[1] = ParentNode->xMin[1] + ParentNode->Children[4]->LengthOfBox;
		ParentNode->Children[4]->xMin[2] = ParentNode->xMin[2] + ParentNode->Children[4]->LengthOfBox;
		
		ParentNode->Children[5]->xMin[0] = ParentNode->xMin[0];
		ParentNode->Children[5]->xMin[1] = ParentNode->xMin[1] + ParentNode->Children[5]->LengthOfBox;
		ParentNode->Children[5]->xMin[2] = ParentNode->xMin[2] + ParentNode->Children[5]->LengthOfBox;

		ParentNode->Children[6]->xMin[0] = ParentNode->xMin[0] + ParentNode->Children[6]->LengthOfBox;
		ParentNode->Children[6]->xMin[1] = ParentNode->xMin[1];
		ParentNode->Children[6]->xMin[2] = ParentNode->xMin[2] + ParentNode->Children[6]->LengthOfBox;
	
		ParentNode->Children[7]->xMin[0] = ParentNode->xMin[0] + ParentNode->Children[7]->LengthOfBox;
		ParentNode->Children[7]->xMin[1] = ParentNode->xMin[1] + ParentNode->Children[7]->LengthOfBox;
		ParentNode->Children[7]->xMin[2] = ParentNode->xMin[2];

	for(i=0;i<ParentNode->NumberOfParticlesInBox;i++)
	{
		nX = 2.0*(Particles[ ParentNode->ParticleIDs[i] ].Pos[0]-ParentNode->xMin[0]) / ParentNode->LengthOfBox; /*the particle is in the left half of the box if nX=0. If nX!=0 the particle is in the right half of the box*/
		nY = 2.0*(Particles[ ParentNode->ParticleIDs[i] ].Pos[1]-ParentNode->xMin[1]) / ParentNode->LengthOfBox;
		nZ = 2.0*(Particles[ ParentNode->ParticleIDs[i] ].Pos[2]-ParentNode->xMin[2]) / ParentNode->LengthOfBox;		

		if(nX==0 && nY==0 && nZ==0)
		{
			BoxNumber = 0;
		}
		else if(nX!=0 && nY==0 && nZ==0)
		{
			BoxNumber = 1;
		}
		else if(nX==0 && nY!=0 && nZ==0)
		{
			BoxNumber = 2;
		}
		else if(nX==0 && nY==0 && nZ!=0)
		{
			BoxNumber = 3;
		}
		else if(nX!=0 && nY!=0 && nZ!=0)
		{
			BoxNumber = 4;
		}
		else if(nX==0 && nY!=0 && nZ!=0)
		{
			BoxNumber = 5;
		}
		else if(nX!=0 && nY==0 && nZ!=0)
		{
			BoxNumber = 6;
		}
		else if(nX!=0 && nY!=0 && nZ==0)
		{
			BoxNumber = 7;
		}
		else
		{
			BoxNumber = 42;
			printf("Error - 42\n");
		}

		/*This statement allocates memory for a new particle in the array, ParentNode->Children[BoxNumber]->ParticleIDs*/
		ParentNode->Children[BoxNumber]->ParticleIDs = realloc(ParentNode->Children[BoxNumber]->ParticleIDs, (ParentNode->Children[BoxNumber]->NumberOfParticlesInBox+1) * sizeof(int));
	
		/*The particle is inserted in child-node..:*/
		ParentNode->Children[BoxNumber]->ParticleIDs[ ParentNode->Children[BoxNumber]->NumberOfParticlesInBox ] = ParentNode->ParticleIDs[i];
		
		/*The node now contains one more particle:*/
		ParentNode->Children[BoxNumber]->NumberOfParticlesInBox++;
		
	}

	return ParentNode;
}

double ComputeMass(struct Node *aNode)
{
	int i;

	if(aNode->EndNode != 1)
	{
		if(aNode==Root) /* parallelization */
		{
			aNode->Mass = 0.0;
			#pragma omp parallel for private(i)
			for(i=0;i<8;i++)
			aNode->Mass += ComputeMass(aNode->Children[i]);
			
		}
		else
		aNode->Mass = ComputeMass(aNode->Children[0]) + ComputeMass(aNode->Children[1])+ComputeMass(aNode->Children[2])+ComputeMass(aNode->Children[3])+ComputeMass(aNode->Children[4])+ComputeMass(aNode->Children[5])+ComputeMass(aNode->Children[6])+ComputeMass(aNode->Children[7]);

		return aNode->Mass;
	}
	else if( aNode->NumberOfParticlesInBox == 1 )
	{
		aNode->Mass = Particles[aNode->ParticleIDs[0]].Mass;
		return aNode->Mass;
	}
	else if(aNode->NumberOfParticlesInBox == 0)
	{
		aNode->Mass = 0.0;
		return 0.0;
	}
	else
	{
		aNode->Mass = 0.0;
		for(i=0;i<aNode->NumberOfParticlesInBox;i++)
		aNode->Mass += Particles[aNode->ParticleIDs[i]].Mass;

		return aNode->Mass;
	}
}

double ComputeCM(int j,struct Node *aNode)/* j=0: x. j=1: y. j=2: z */
{
	int i;

	if(aNode->EndNode != 1)
	{
		aNode->CM[j] = 1.0/aNode->Mass * (
		ComputeCM(j,aNode->Children[0])*aNode->Children[0]->Mass + 
		ComputeCM(j,aNode->Children[1])*aNode->Children[1]->Mass + 
		ComputeCM(j,aNode->Children[2])*aNode->Children[2]->Mass + 
		ComputeCM(j,aNode->Children[3])*aNode->Children[3]->Mass + 
		ComputeCM(j,aNode->Children[4])*aNode->Children[4]->Mass + 
		ComputeCM(j,aNode->Children[5])*aNode->Children[5]->Mass + 
		ComputeCM(j,aNode->Children[6])*aNode->Children[6]->Mass +
		ComputeCM(j,aNode->Children[7])*aNode->Children[7]->Mass );
		return aNode->CM[j];
	}
	else if( aNode->NumberOfParticlesInBox == 1 )
	{
		aNode->CM[j] = Particles[aNode->ParticleIDs[0]].Pos[j];
		return aNode->CM[j];
	}
	else if(aNode->NumberOfParticlesInBox == 0)
	{
		aNode->CM[j] = 0.0;
		return 0.0;
	}
	else
	{
		aNode->CM[j] = 0.0;
		for(i=0;i<aNode->NumberOfParticlesInBox;i++)
		aNode->CM[j] += Particles[aNode->ParticleIDs[i]].Pos[j]  *  Particles[aNode->ParticleIDs[i]].Mass;

		aNode->CM[j] = aNode->CM[j] / aNode->Mass;
		return aNode->CM[j];
	}
}


void GetMaxSmoothingLength(struct Node *aNode)
{
	int i;
	double temp;

	if(aNode->EndNode != 1)
	{
		aNode->MaxSmoothingLength = 0.0;
		for(i=0;i<8;i++)
		{
			GetMaxSmoothingLength(aNode->Children[i]);
			if(aNode->Children[i]->MaxSmoothingLength > aNode->MaxSmoothingLength)
			aNode->MaxSmoothingLength = aNode->Children[i]->MaxSmoothingLength;
		}
	}
	else if( aNode->NumberOfParticlesInBox == 1 )
	{
		if(aNode->ParticleIDs[0]<Nsph)
		aNode->MaxSmoothingLength = Particles[aNode->ParticleIDs[0]].SmoothingLength;
	}
	else if(aNode->NumberOfParticlesInBox == 0)
	{
		aNode->MaxSmoothingLength = 0.0;
	}
	else
	{
		if(aNode->ParticleIDs[0]<Nsph)
		temp = Particles[aNode->ParticleIDs[0]].SmoothingLength;
		else
		temp=0.0;
				
		for(i=0;i<aNode->NumberOfParticlesInBox;i++)
		{
			if(Particles[aNode->ParticleIDs[i]].SmoothingLength > temp && aNode->ParticleIDs[i] < Nsph)
			temp = Particles[aNode->ParticleIDs[i]].SmoothingLength;
		}

		aNode->MaxSmoothingLength = temp;
	}
}
