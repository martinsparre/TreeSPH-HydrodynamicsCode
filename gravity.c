#include "all.h"


void GravityAcceleration(int i)
{
	Particles[i].GravAcc[0] = 0.0;
	Particles[i].GravAcc[1] = 0.0;
	Particles[i].GravAcc[2] = 0.0;
	
	#ifdef DIRECTSUMMATION
		ComputeAcc_DirectSummation();
	#else
		ComputeAcc(i,Root);
	#endif
}


void ComputeAcc(int aParticleID, struct Node *A)
{
	double TEMP;
	double DistSquared = DistanceSquared( Particles[aParticleID].Pos , A->CM );
	
	if( DistSquared * pow(OPENING_ANGLE,2) > pow(A->LengthOfBox,2) || A->NumberOfParticlesInBox < 2)
	{
		TEMP = - G*A->Mass/pow(  DistSquared + pow(SOFTENING,2)  ,1.5);

		Particles[aParticleID].GravAcc[0] = Particles[aParticleID].GravAcc[0] + TEMP * (Particles[aParticleID].Pos[0] - A->CM[0]);
		Particles[aParticleID].GravAcc[1] = Particles[aParticleID].GravAcc[1] + TEMP * (Particles[aParticleID].Pos[1] - A->CM[1]);
		Particles[aParticleID].GravAcc[2] = Particles[aParticleID].GravAcc[2] + TEMP * (Particles[aParticleID].Pos[2] - A->CM[2]);
	}
	else if( A->EndNode == 1)
	{
		TEMP = - G*A->Mass/pow(  DistSquared + pow(SOFTENING,2)  ,1.5);

		Particles[aParticleID].GravAcc[0] = Particles[aParticleID].GravAcc[0] + TEMP * (Particles[aParticleID].Pos[0] - A->CM[0]);
		Particles[aParticleID].GravAcc[1] = Particles[aParticleID].GravAcc[1] + TEMP * (Particles[aParticleID].Pos[1] - A->CM[1]);
		Particles[aParticleID].GravAcc[2] = Particles[aParticleID].GravAcc[2] + TEMP * (Particles[aParticleID].Pos[2] - A->CM[2]);
		
	}
	else
	{
		int i;
		for(i=0;i<8;i++)
		{
		ComputeAcc(aParticleID,A->Children[i]);			
		}
	}
}

#ifdef DIRECTSUMMATION
void ComputeAcc_DirectSummation(void)
{
    int k;
    for(k=0; k<NumberOfParticles;k++)
    {
        Particles[k].GravAcc[0] = 0.0;
        Particles[k].GravAcc[1] = 0.0;
        Particles[k].GravAcc[2] = 0.0;
    }

    double Ax,Ay,Az,m_i,m_j,temp;
    Ax=0;Ay=0;Az=0;
    int i,j;
 
    for(i=0;i<NumberOfParticles;i++)
    {
         m_i = Particles[i].Mass;

        for(j=i+1;j<NumberOfParticles;j++)
        {
            m_j = Particles[j].Mass;
            temp = G*m_j*m_i/(pow(  DistanceSquared(Particles[i].Pos,Particles[j].Pos) + pow(SOFTENING,2.0) ,1.5));

            Ax = temp*(Particles[j].Pos[0]-Particles[i].Pos[0]);
            Ay = temp*(Particles[j].Pos[1]-Particles[i].Pos[1]);
            Az = temp*(Particles[j].Pos[2]-Particles[i].Pos[2]);

			Particles[i].GravAcc[0] = Particles[i].GravAcc[0] + Ax;
			Particles[j].GravAcc[0] = Particles[j].GravAcc[0] - Ax;

			Particles[i].GravAcc[1] = Particles[i].GravAcc[1] + Ay;
			Particles[j].GravAcc[1] = Particles[j].GravAcc[1] - Ay;
			
			Particles[i].GravAcc[2] = Particles[i].GravAcc[2] + Az;
			Particles[j].GravAcc[2] = Particles[j].GravAcc[2] - Az;
        }
    }
    
    for(k=0; k<NumberOfParticles;k++)
    {
        Particles[k].GravAcc[0] = Particles[k].GravAcc[0] /  Particles[k].Mass;
        Particles[k].GravAcc[1] = Particles[k].GravAcc[1] /  Particles[k].Mass;
        Particles[k].GravAcc[2] = Particles[k].GravAcc[2] /  Particles[k].Mass;
    }
}
#endif

void GetGravPotential(void)
{
	int i;
	
	for(i=0;i<Nparticles;i++)
	{
		Particles[i].GravPotential = 0.0;
	}

	#pragma omp parallel for
	for(i=0;i<Nparticles;i++)	
		ComputeGravPotential(i,Root);
}

void ComputeGravPotential(int aParticleID, struct Node *A)
{
	int i;
	double DistSquared = DistanceSquared( Particles[aParticleID].Pos , A->CM );
	
	if( DistSquared * pow(OPENING_ANGLE,2) > pow(A->LengthOfBox,2) || A->NumberOfParticlesInBox < 2)
	{
		Particles[aParticleID].GravPotential = Particles[aParticleID].GravPotential - G*A->Mass/pow( DistSquared + pow(SOFTENING,2)  ,0.5);

	}
	else if( A->EndNode == 1)
	{
		Particles[aParticleID].GravPotential = Particles[aParticleID].GravPotential - G*A->Mass/pow( DistSquared + pow(SOFTENING,2)  ,0.5);
	}
	else
	{
		for(i=0;i<8;i++)
		{
		ComputeGravPotential(aParticleID,A->Children[i]);			
		}
	}
}

double DistanceSquared(double *Pos1, double *Pos2)
{
	return ( pow(Pos2[0]-Pos1[0],2) + pow(Pos2[1]-Pos1[1],2) + pow(Pos2[2]-Pos1[2],2) );
}
