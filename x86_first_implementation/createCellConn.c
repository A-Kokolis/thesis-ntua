#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>

int main(int argc, char* argv[]) {

	int IO_NETWORK_DIM1= atoi(argv[1]);
	int IO_NETWORK_DIM2= atoi(argv[2]);
	int IO_NETWORK_SIZE= IO_NETWORK_DIM1*IO_NETWORK_DIM2;
	char conFile[100];
	sprintf(conFile,"/home/apostolis/Desktop/mpi_programs/ckpt_infoli/input/cellConnections.txt");
	FILE *pConFile = fopen(conFile, "w+");
	printf("Open Sesame\n");
	int i, j, k, x, y;
	for (i=0; i<IO_NETWORK_SIZE; i++) {
		x = i / IO_NETWORK_DIM2;
		y = i % IO_NETWORK_DIM2;
		
		for(j=x-1; j<=x+1; j++){
                        for(k=y-1; k<=y+1; k++){
				if ((j==x)&&(k==y))
					;
				else if ((j<0)||(k<0)||(j>=IO_NETWORK_DIM1)||(k>=IO_NETWORK_DIM2))
					;
				else
					fprintf(pConFile,"%d B ",(j*IO_NETWORK_DIM2+k));
			}
		}
		
		fprintf(pConFile,"|\n");		
	}
	fclose(pConFile);
	return 1;
	
}
