#include "all.h"

int a;

int main()
{		
	LogFile=fopen(LOGFILE_NAME,"w+");

	Initialize();
	Run();

	fclose(LogFile);

return 0;
}
