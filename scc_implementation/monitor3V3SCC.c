#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/uio.h>
#include <sys/wait.h>
#include <fcntl.h> 
#include <signal.h>

 struct timeval ttstart,ttemp;

 char monitoring_command[1000]="sccBmc -c status | grep 3V3SCC >> monitor_3V3SCC.txt";
  char temp[1000];



void monitor_voltages(){
   
    gettimeofday(&ttstart,NULL);
    while(1){
	system(monitoring_command);
	gettimeofday(&ttemp,NULL);	
	sprintf(temp,"echo command executed timestep= %lf>>monitor_3V3SCC.txt",(ttemp.tv_sec - ttstart.tv_sec) + (ttemp.tv_usec - ttstart.tv_usec) * 0.000001);
	system(temp);
    }
}


int main (int argc, char *argv[]){
   
   monitor_voltages();

return 0;
}
