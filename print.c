#include "all.h"

void PrintToFile(int i)
{
	FILE *FileA;
	int j,k;
	char string2[100];
	
	GetGravPotential();
	
	for(j=0;j<Nsph;j++)
	Particles[j].Energy = Particles[j].PredictedEntropy / (GAMMA-1.0) * pow(Particles[j].Density, GAMMA-1);
	
	    
    strcpy(string2,OUTPUT);
    strcat(string2,itoa(i));
    strcat(string2,".ascii");
    

    FileA = fopen(string2,"w+");    
    
    fprintf(FileA,"x\ty\tz\tvx\tvy\tvz\tM\tu\trho\tA\th\tV\n"); 

    for(k=0;k<Nsph;k++)
    {
 fprintf(FileA,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", Particles[k].Pos[0],Particles[k].Pos[1],Particles[k].Pos[2],Particles[k].PredictedVel[0],Particles[k].PredictedVel[1],Particles[k].PredictedVel[2],Particles[k].Mass,Particles[k].Energy,Particles[k].Density,Particles[k].PredictedEntropy,Particles[k].SmoothingLength,Particles[k].GravPotential); 
    }
    
    for(;k<Nparticles;k++)
    {
 fprintf(FileA,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", Particles[k].Pos[0],Particles[k].Pos[1],Particles[k].Pos[2],Particles[k].PredictedVel[0],Particles[k].PredictedVel[1],Particles[k].PredictedVel[2],Particles[k].Mass,Particles[k].GravPotential); 
    }
    
    printf("%s \t %lf \n",string2,t);
    fprintf(LogFile,"%s \t %lf \n",string2,t);
	
    fclose(FileA);
}

/*itoa - integer to char (http://www.freebookzone.com/others/itoa.h) */
char* itoa (int n)
{
    int i=0,j;
    char* s;
    char* u;
    
    s= (char*) malloc(17);
    u= (char*) malloc(17);
    
    do{
    s[i++]=(char)( n%10+48 );
    n-=n%10;
    }
    while((n/=10)>0);
    for (j=0;j<i;j++)
    u[i-1-j]=s[j];
    
    u[j]='\0';
    return u;
}
