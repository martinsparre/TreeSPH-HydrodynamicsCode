#include "all.h"

double dWdh(double R, double h) /* dW/dh */
{
	double hinv5,hinv4;
	double R2,R3,temp;
	
	if(R <= h/2.0)
	{
		hinv4 = 1.0/(h*h*h*h);
		R2=R*R;
		R3=R*R*R;
	return -24.*hinv4 / PI * ( 1.-6.*R2*h*h*hinv4+6.*R3*h*hinv4 ) + 8./PI *  hinv4 * h * (12.*R2 *hinv4*h-18.*R3*hinv4 ); 
	}
	else if(R>h/2.0 && R<h)
	{
		hinv5 = 1.0/(h*h*h*h*h);
		hinv4 = h*hinv5;
		temp = 1.-R/h;
		return -24.*hinv4 / PI * ( 2. * temp*temp*temp ) + 8./(PI)  * hinv5 * 6*temp*temp*R;
	}
	else
	return 0.;
}

double dmin(double x, double y)
{
	if(x < y)
		return x;
	else
		return y;
}

double dmax(double x, double y)
{
	if(x > y)
		return x;
	else
		return y;
}

void HydroAccAndEntropyDerivative(int i)
{
	int j,k,NeighborID;	
	double acc,R,AccVisc, vdotr, NablaKernel1, NablaKernel2,ViscTensor,DEntropy;
	

	Particles[i].HydroAcc[0] = 0.0;
	Particles[i].HydroAcc[1] = 0.0;
	Particles[i].HydroAcc[2] = 0.0;
		
	DEntropy = 0.0;

	/* the accelerations */
	for(j=0;j<Particles[i].NumberOfNeighbors;j++)
	{
		NeighborID = Particles[i].ListOfNeighbors[j];
		R = sqrt(DistanceSquared( Particles[i].Pos, Particles[NeighborID].Pos)) ;
		
		if(R>0.0)
		{
			ViscTensor = ViscosityTensor(i,NeighborID);
				
			NablaKernel1 = NablaKernel(i,NeighborID, Particles[i].SmoothingLength);
			NablaKernel2 = NablaKernel(i,NeighborID, Particles[NeighborID].SmoothingLength);

			vdotr=0.0;
				
			for(k=0;k<3;k++)
			{
				vdotr += (Particles[i].Pos[k]-Particles[NeighborID].Pos[k])*(Particles[i].PredictedVel[k]-Particles[NeighborID].PredictedVel[k]);	
			}
			
			//TEST:
/*			if(vdotr<0 && (NablaKernel1 + NablaKernel2) < 0)
			{
			  ViscTensor = dmin(ViscTensor, 0.5 * 1 * vdotr /
				      (0.5 * (Particles[i].Mass+Particles[NeighborID].Mass) * (NablaKernel1 + NablaKernel2) * R * 0.02));
			}*/
				
			/* Acceleration due to to viscosity divided by R*/				
			AccVisc = Particles[NeighborID].Mass * ViscTensor * 0.5 * (NablaKernel1+NablaKernel2) / R;
					
			/* Acceleration due to "ordinary" hydrodynamics divided by R */
			acc = -AccVisc - Particles[NeighborID].Mass *
			(
			Particles[i].f*
			Particles[i].Pressure / pow(Particles[i].Density , 2) 
			* NablaKernel1
			+
			Particles[NeighborID].f*
			Particles[NeighborID].Pressure / pow(Particles[NeighborID].Density , 2)
			* NablaKernel2
			)/R;				

			for(k=0;k<3;k++)
			{
				Particles[i].HydroAcc[k] += acc*(Particles[i].Pos[k]-Particles[NeighborID].Pos[k]);
			}
			/* Entropy Derivative */
			DEntropy += 0.5*Particles[NeighborID].Mass * ViscTensor * vdotr * (NablaKernel1 + NablaKernel2)/R;
		}
	}

	Particles[i].DEntropyDt = DEntropy * 0.5 * (GAMMA - 1.0) / pow( Particles[i].Density , GAMMA-1.0);

}

double ViscosityTensor(int i, int j)
{
	int k;
	double vdotr = 0.0;
	
	for(k=0;k<3;k++)
	{
		vdotr += ( Particles[i].Pos[k] - Particles[j].Pos[k] ) * ( Particles[i].PredictedVel[k] - Particles[j].PredictedVel[k] );
	}

	if(vdotr >= 0.0)
	{
		return 0.0;	
	}
	else
	{
		vdotr = vdotr / sqrt(DistanceSquared(Particles[i].Pos,Particles[j].Pos));

		return (- ALPHA* ( SpeedOfSound(i) + SpeedOfSound(j) - 3*vdotr )*vdotr 
		/ ( 
		( Particles[i].Density + Particles[j].Density)
		));
	}
}

double SPH_Kernel(double R, double h)
{
	double W,u;
	
	u=R/h;
	
	if(u <= 0.5)
		W = 1.0-6.0*pow(u,2)+6.0*pow(u,3);
	else if(u > 0.5 &&  u < 1.0)
		W = 2.0*pow(1.0-u,3);
	else
		W=0.0;

	if(W>0)
	W *= 2.5464790894703253722 /( pow(h,3) ); /* 8.0 / Pi = 2.5464790894703253722 */
	
	return W;
}

double NablaKernel(int Particle1, int Particle2, double h) /* h = smoothinglength. */
{
	double R,u;

	R = sqrt(DistanceSquared( Particles[ Particle1 ].Pos , Particles[ Particle2 ].Pos ));
	u = R/h;

	if(u <= 0.5)
		return  u * ((45.836623610466) * u - (30.557749073644)) / pow(h,4);
	else if(u > 0.5 && u < 1.0)
		return -15.278874536822 * (1.0 - u) * (1.0 - u) / pow(h,4);
	else
		return 0.0;
}

double SpeedOfSound(int i) /* The speed of sound of the i'th particle */
{
 return sqrt( GAMMA * Particles[i].Pressure / Particles[i].Density   );
}
