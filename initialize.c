#include "all.h"

void Initialize(void)
{
	int i;

	FileCounter=MAX_NR_OF_OUTPUT_FILES;
	t=0.0;

	printf("Start\n");
	printf("END TIME: %lf\t TIME STEP: %lf\n",tmax,dt);
	fprintf(LogFile,"Start\n");
	fprintf(LogFile,"END TIME: %lf\t TIME STEP: %lf\n",tmax,dt);
	
	if(SOFTENING == 0.0)
	{
	printf("Error - softening is 0.0\n");
	fprintf(LogFile,"Error - softening is 0.0\n");
	}
	
	ImportData();
	
	for(i=0;i<Nsph;i++) /* Initial smoothong length */
	{
		Particles[i].SmoothingLength = 0.5;
	}

	
	/* velocities at t=0 */
	for(i=0; i<Nparticles;i++)
	{
		Particles[i].PredictedVel[0] = Particles[i].Vel[0];
		Particles[i].PredictedVel[1] = Particles[i].Vel[1];
		Particles[i].PredictedVel[2] = Particles[i].Vel[2];
	}

}

void ImportData(void)
{
    int i;
    double TempValue;
    FILE *Filename;

    Filename = fopen(INPUT_FILE_NAME,"r+");
    
    for (i=0;i<Nsph;i++)
    {
			fscanf(Filename, "%lf", &TempValue);
			Particles[i].Pos[0] = TempValue;

			fscanf(Filename, "%lf", &TempValue);
			Particles[i].Pos[1] = TempValue;
			
			fscanf(Filename, "%lf", &TempValue);
			Particles[i].Pos[2] = TempValue;
			
			fscanf(Filename, "%lf", &TempValue);
			Particles[i].Vel[0] = TempValue;

			fscanf(Filename, "%lf", &TempValue);
			Particles[i].Vel[1] = TempValue;

			fscanf(Filename, "%lf", &TempValue);
			Particles[i].Vel[2] = TempValue;

			fscanf(Filename, "%lf", &TempValue);
			Particles[i].Mass = TempValue;
			
			fscanf(Filename, "%lf", &TempValue);
			Particles[i].Energy = TempValue;
    }

    for (i=Nsph;i<Nparticles;i++)
    {
			fscanf(Filename, "%lf", &TempValue);
			Particles[i].Pos[0] = TempValue;

			fscanf(Filename, "%lf", &TempValue);
			Particles[i].Pos[1] = TempValue;
			
			fscanf(Filename, "%lf", &TempValue);
			Particles[i].Pos[2] = TempValue;
			
			fscanf(Filename, "%lf", &TempValue);
			Particles[i].Vel[0] = TempValue;

			fscanf(Filename, "%lf", &TempValue);
			Particles[i].Vel[1] = TempValue;

			fscanf(Filename, "%lf", &TempValue);
			Particles[i].Vel[2] = TempValue;

			fscanf(Filename, "%lf", &TempValue);
			Particles[i].Mass = TempValue;
			
//			fscanf(Filename, "%lf", &TempValue);
//			Particles[i].Energy = TempValue;
    }


    fclose(Filename);

	printf("Data has been imported\n");
	fprintf(LogFile,"Data has been imported\n");
}
