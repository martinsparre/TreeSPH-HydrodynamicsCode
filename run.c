#include "all.h"

void Run(void)
{
	int i,Counter=PRINT_EVERY_N_TIMESTEP;
	
	Root = NULL;
	BuildTree(); /*The tree is used to calculate the density*/

	/*SPH-stuff:*/
	#pragma omp parallel for
	for(i=0;i<Nsph;i++)	
	UpdateSmoothing( i, 0); /* Initial smoothing length */
	
	printf("1\n");
	
	#pragma omp parallel for
	for(i=0;i<Nsph;i++)
	DetermineSmoothingLength_Newton(i);
	
	GetMaxSmoothingLength(Root);

	printf("2\n");

	#pragma omp parallel for
	for(i=0;i<Nsph;i++)	
	{
		ResetNeighborList(i);
		CreateNeighborList(i,Root);
	}
	
	#pragma omp parallel for
	for(i=0;i<Nsph;i++)
	DensityEvaluate(i, Particles[i].SmoothingLength);
	
		printf("3\n");
	#pragma omp parallel for
	for(i=0; i<Nsph;i++)/* Initial Pressure */
	{
		Particles[i].Pressure = (GAMMA-1.0) * Particles[i].Density * Particles[i].Energy;
	}
	
	#pragma omp parallel for
	for(i=0; i<Nsph;i++) /* Initial Entropy */
	{
		Particles[i].Entropy = Particles[i].Pressure *  pow(Particles[i].Density,- GAMMA);
		Particles[i].PredictedEntropy = Particles[i].Pressure *  pow(Particles[i].Density,- GAMMA);
	}

   	/*
	 * The Velocity and the entropy are advanced to t=0.5*dt 
	 * First the accelerion is calculated at t=0 and then the quantities are computed
	 */
	#pragma omp parallel for
	for(i=0;i<Nparticles;i++)
	GravityAcceleration(i);
		
			printf("3\n");
	#pragma omp parallel for
	for(i=0;i<Nsph;i++)	
	HydroAccAndEntropyDerivative(i); /*The Densities must be known when this function is called. */

	#pragma omp parallel for
	for(i=0;i<Nsph;i++)/* computation of velocity and Entropy at t=dt/2.0 */
	{
		Particles[i].Vel[0] += ( Particles[i].GravAcc[0] + Particles[i].HydroAcc[0] ) * dt / 2.0;
		Particles[i].Vel[1] += ( Particles[i].GravAcc[1] + Particles[i].HydroAcc[1] ) * dt / 2.0;
		Particles[i].Vel[2] += ( Particles[i].GravAcc[2] + Particles[i].HydroAcc[2] ) * dt / 2.0;
		
		Particles[i].Entropy += Particles[i].DEntropyDt * dt / 2.0;
	}

	#pragma omp parallel for
	for(i=Nsph;i<Nparticles;i++)/* computation of velocity and Entropy at t=dt/2.0 */
	{
		Particles[i].Vel[0] += Particles[i].GravAcc[0] * dt / 2.0;
		Particles[i].Vel[1] += Particles[i].GravAcc[1] * dt / 2.0;
		Particles[i].Vel[2] += Particles[i].GravAcc[2] * dt / 2.0;
	}

    printf("The first step has been completed\n");
	printf("FILE\t \t \t \t TIME \n");
    fprintf(LogFile,"The first step has been completed\n");
	fprintf(LogFile,"FILE\t \t \t \t TIME \n");

	printf("4\n");
    /* MAIN LOOP */
    while(t <= tmax)
    {
    	/* Print to files */
        if(Counter==PRINT_EVERY_N_TIMESTEP)
        {
            PrintToFile(FileCounter);
            FileCounter++;
            Counter=0;
        }

		/* x(n+1) is calculated: x(n+1) = x(n) + v(n+1/2)*dt */       
		#pragma omp parallel for
		for(i=0;i<Nparticles;i++)
		{
			Particles[i].Pos[0] = Particles[i].Pos[0] + Particles[i].Vel[0]*dt;
			Particles[i].Pos[1] = Particles[i].Pos[1] + Particles[i].Vel[1]*dt;
			Particles[i].Pos[2] = Particles[i].Pos[2] + Particles[i].Vel[2]*dt;
		}   
        
        /* The tree is destroyed and a new tree is build. Then the densities are computed */
        DestroyTree(Root);
        BuildTree();
        
       	#pragma omp parallel for
		for(i=0;i<Nsph;i++)
		DetermineSmoothingLength_Newton(i);
		
		GetMaxSmoothingLength(Root);

		#pragma omp parallel for
		for(i=0;i<Nsph;i++)	
		{
			ResetNeighborList(i);
			CreateNeighborList(i,Root);
		}
		
		#pragma omp parallel for
		for(i=0;i<Nsph;i++)
		DensityEvaluate(i, Particles[i].SmoothingLength);
        
        /* The velocity and entropy are predicted at timestep n+1:  */
		#pragma omp parallel for
		for(i=0;i<Nsph;i++)
		{

			Particles[i].PredictedVel[0] = Particles[i].Vel[0] + (Particles[i].GravAcc[0] + Particles[i].HydroAcc[0]) * dt / 2.0;
			Particles[i].PredictedVel[1] = Particles[i].Vel[1] + (Particles[i].GravAcc[1] + Particles[i].HydroAcc[1]) * dt / 2.0;
			Particles[i].PredictedVel[2] = Particles[i].Vel[2] + (Particles[i].GravAcc[2] + Particles[i].HydroAcc[2]) * dt / 2.0;
			
			Particles[i].PredictedEntropy = Particles[i].Entropy + Particles[i].DEntropyDt * dt / 2.0;		
		}

		#pragma omp parallel for
		for(i=Nsph;i<Nparticles;i++)
		{
			Particles[i].PredictedVel[0] = Particles[i].Vel[0] + (Particles[i].GravAcc[0]) * dt / 2.0;
			Particles[i].PredictedVel[1] = Particles[i].Vel[1] + (Particles[i].GravAcc[1]) * dt / 2.0;
			Particles[i].PredictedVel[2] = Particles[i].Vel[2] + (Particles[i].GravAcc[2]) * dt / 2.0;		
		}

		
		#pragma omp parallel for
		for(i=0;i<Nparticles;i++)
		GravityAcceleration(i);
		
		#pragma omp parallel for
		for(i=0;i<Nsph;i++)
		Particles[i].Pressure = Particles[i].PredictedEntropy * pow( Particles[i].Density , GAMMA );
		
		#pragma omp parallel for
		for(i=0;i<Nsph;i++)	
		HydroAccAndEntropyDerivative(i);
		
		/* Velocity and entropy are calculated at step n+3/2 */
		#pragma omp parallel for
		for(i=0;i<Nsph;i++)
		{
			Particles[i].Vel[0] = Particles[i].Vel[0] + (Particles[i].GravAcc[0] + Particles[i].HydroAcc[0]) * dt;
			Particles[i].Vel[1] = Particles[i].Vel[1] + (Particles[i].GravAcc[1] + Particles[i].HydroAcc[1]) * dt;
			Particles[i].Vel[2] = Particles[i].Vel[2] + (Particles[i].GravAcc[2] + Particles[i].HydroAcc[2]) * dt;

			Particles[i].Entropy += Particles[i].DEntropyDt * dt;
		}
		
		#pragma omp parallel for
		for(i=Nsph;i<Nparticles;i++)
		{
			Particles[i].Vel[0] = Particles[i].Vel[0] + (Particles[i].GravAcc[0]) * dt;
			Particles[i].Vel[1] = Particles[i].Vel[1] + (Particles[i].GravAcc[1]) * dt;
			Particles[i].Vel[2] = Particles[i].Vel[2] + (Particles[i].GravAcc[2]) * dt;
		}
		
		t += dt;
		Counter++;
	}

	PrintToFile(FileCounter);

	printf("Run completed\n");
	fprintf(LogFile,"Run completed\n");

	printf("A logfile has been printed to LOGFILE_NAME\n");
}
