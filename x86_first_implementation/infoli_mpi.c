/*
 *
 * Copyright (c) 2012, Neurasmus B.V., The Netherlands,
 * web: www.neurasmus.com email: voltage@neurasmus.com
 *
 * Any use or reproduction in whole or in parts is prohibited
 * without the written consent of the copyright owner.
 *
 * All Rights Reserved.
 *
 *
 * Author: Sebastian Isaza
 * Created: 19-01-2012
 * Modified: 07-08-2012
 *
 * Description: Top source file of the Inferior Olive model, originally written
 * in Matlab by Jornt De Gruijl. It contains the implementation of all functions.
 * The main function allocates the necessary memory, initializes the system
 * state and runs the model calculations.
 *
 */

/*




 * we assume that dim1 refers to the number of ROWS
 * and dim2 refers to COLUMNS (cell network size).
 */


#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <signal.h>
#include <stddef.h>
#include <math.h>
#include <sys/time.h>
#include <sys/file.h>
#include <time.h>
#include <mpi.h>
#include "infoli_mpi.h"

#include <errno.h>

int core_id, cores, cellCount;
int IO_NETWORK_DIM1, IO_NETWORK_DIM2, IO_NETWORK_SIZE,CKPT_INTERVAL;
FILE  *ckptfd;

struct timeval tic, toc, intime;

typedef unsigned long long timestamp_t;

static timestamp_t get_timestamp ()                                             
{                                                                               
    struct timeval now;                                                         
    gettimeofday (&now, NULL);                                                  
    return now.tv_usec + (timestamp_t)now.tv_sec * 1000000;                     
} 

char fpeek(FILE *stream)
{
	char c;

	c = (char)fgetc(stream);
	ungetc(((int)c), stream);

	return c;
}


int main(int argc, char *argv[]){

	int i=0, j, k=0, line_count,l, p, q, x, y, targetCore;
	char c;
	char *inFileName;
	char outFileName[100];
	FILE *pInFile, *coreF;
	char conFile[100],core_file[100];
	FILE *pConFile,*pOutFile;
    FILE *restartTimeFile;
    char restartFileName[100];
    FILE *ckptTimeFile;
    char ckptTimeFileName[100];

    sprintf(ckptTimeFileName,"/home/apostolis//Desktop/mpi_programs/ckpt_infoli/ckpTime/ckpTime.txt");
    sprintf(restartFileName,"/home/apostolis/Desktop/mpi_programs/ckpt_infoli/ckpTime/applicationRestartTime.txt");
	
    mod_prec* iAppArray= NULL;
	int simStep = 0, inputFromFile = 0, initSteps,total_simulation_steps;
	float simTime = 0;
	cellState **cellPtr;
	cellCompParams *cellParamsPtr;
	int seedvar;	
	mod_prec iApp;
	mod_prec *gcal;
	int *init_steps;

	//CHECKPOINT DEFINITIONS
    char tempbuf[OUTFILE_MAXLINE];		
	char ckptFileName[100];
	char hostsFile[100];
    FILE * hosts;
    void * ckptBuffer;
    
    FILE * zero_ckptfd;
    int prev_cores,prev_cellCount;

	long double seconds;
	int simulation_array_ID,target_cell;
	communication_node *communicating_cell,*next_in_list,*communication_list_head;
	sending_node *sending_list_head = NULL;
	receiving_node *receiving_list_head = NULL;
	clock_t start, end;
	struct timeval computingtime, messagingtime, outtime;
    struct timeval tts,ttf;


	MPI_Init(&argc, &argv);


	start = clock();
	MPI_Comm_rank(MPI_COMM_WORLD, &core_id);
	MPI_Comm_size(MPI_COMM_WORLD, &cores);


    if (MAX_CORES % cores != 0) {
        printf("Error: %d cores specified\n", cores);
        exit(EXIT_NUM_UES);
    }


	/* outFileName is file for output
	 * conFile is input of the cell's connections

	 */
	sprintf(core_file,"core%d",core_id);		
	sprintf(outFileName,"/home/apostolis/Desktop/mpi_programs/ckpt_infoli/output/InferiorOlive_Output%d.txt",core_id);		
	sprintf(conFile,"/home/apostolis/Desktop/mpi_programs/ckpt_infoli/input/cellConnections.txt");


	IO_NETWORK_DIM1= atoi(argv[1]);
	IO_NETWORK_DIM2= atoi(argv[2]);
	IO_NETWORK_SIZE= IO_NETWORK_DIM1*IO_NETWORK_DIM2;
	CKPT_INTERVAL= atoi(argv[3]);
	
    /* TEST - timer initialization
     */
     timestamp_t timestamp, timediff, now;
     int counter = 0;

	/* compute how many grid cells 
	 * are assigned to each core
	 */
        cellCount= (IO_NETWORK_SIZE / cores);

    /* ckptFilename is the core's checkpoint file name.
	 * if restarting won't be performed, clear the file and write the header
	 */
	sprintf(ckptFileName, "/home/apostolis/Desktop/mpi_programs/ckpt_infoli/output/ckptFile%d.bin", core_id);
	if (RESTARTING) {
		if (core_id == 0)
			zero_ckptfd = fopen(ckptFileName, "r+");
		else {
			sprintf(ckptFileName, "/home/apostolis/Desktop/mpi_programs/ckpt_infoli/output/ckptFile0.bin");
			zero_ckptfd = fopen(ckptFileName,"r+");
            if (zero_ckptfd == NULL) {
                printf("Error: Couldn't open %s core_id=%d\n", ckptFileName,core_id);
                exit(EXIT_OPEN_FAIL);
            }
		}

        	/* Read the cellCount and cores variables of the checkpointed simulation */
        	flock(fileno(zero_ckptfd), LOCK_SH);
			fread(&prev_cores, sizeof(int), 1, zero_ckptfd);
			fread(&prev_cellCount, sizeof(int), 1, zero_ckptfd);
        	flock(fileno(zero_ckptfd), LOCK_UN);
	
		sprintf(ckptFileName, "/home/apostolis/Desktop/mpi_programs/ckpt_infoli/output/ckptFile%d.bin", core_id);
		if (core_id>=prev_cores){
			ckptfd = fopen(ckptFileName, "w+");
			printf("Error %d \n", errno);
			printf("the ckptfd is %p\n",ckptfd);
		}else{
			ckptfd = fopen(ckptFileName, "r+");
			printf("Error %d \n", errno);
			printf("the ckptfd is %p\n",ckptfd);
		}
	} else if (CHECKPOINTING) {
		ckptfd = fopen(ckptFileName, "w");
		fwrite(&cores, sizeof(int), 1, ckptfd);
		fwrite(&cellCount, sizeof(int), 1, ckptfd);
	} 
    if ((RESTARTING || CHECKPOINTING) && (ckptfd == NULL)) {
        printf("Error: Couldn't open %s\n", ckptFileName);
    	printf("here1\n");
        exit(EXIT_OPEN_FAIL);
    }
    if (CHECKPOINTING) {
        ckptBuffer = malloc(cellCount * (sizeof(cellState) + sizeof(int)) + sizeof(int));
		//chmod(ckptFileName, 0x77777);
        chmod(ckptFileName,S_IRUSR | S_IWUSR);
        //sprintf("chmod 777 %s",ckptFileName);
        //system(tempbuf);
		sprintf(tempbuf, "chown apostolis %s", ckptFileName);
        system(tempbuf);
        //chown(ckptFileName, NULL, (gid_t) (getgrnam("users"))->gr_gid);
    }


	/* Process command line arguments
	 * Case argc = 4 then a one-pulse input is stimulated.
	 * Otherwise we receive input from a specified file in argv[4] and simulation runs accordingly
	 * in the case of inputFromFile, we will also need a buffer to hold the stimulus
	 * for each cell on each step
  	 */
	
	if(argc == 4) {
		inputFromFile = 0;
	}
	else if(argc == 5) {
		inputFromFile = 1;
		inFileName = argv[4];
		pInFile = fopen(inFileName,"r");
		if(pInFile==NULL) {
			printf("Error: Couldn't open %s\n", inFileName);
			exit(EXIT_OPEN_FAIL);
		}
		iAppArray= (mod_prec*) malloc(cellCount*(sizeof(mod_prec)));
	}	
	else {
		printf("Error: Too many arguments.\nUsage: ./InferiorOlive <dim1> <dim2> <CKPT_INTERVAL> <Iapp_input_file> or ./InferiorOlive <dim1> <dim2> <CKPT_INTERVAL>\n");
		exit(EXIT_ARGUMENTS);
	}

	
	/* PRINTING is true in case we want to write the
	 * output to a specified file (outFileName)
	 */
	if (PRINTING && !RESTARTING) {
		pOutFile = fopen(outFileName,"w");
		if(pOutFile==NULL){
			printf("Error: Couldn't create %s\n", outFileName);
			exit(EXIT_OPEN_FAIL);
		}
		//chmod(outFileName, 0x77777);
        chmod(outFileName,S_IRUSR | S_IWUSR);
		sprintf(tempbuf, "chown apostolis %s", ckptFileName);
        system(tempbuf);
		sprintf(tempbuf, "#simSteps Time(ms) Input(Iapp) Output(V_soma)");
		fputs(tempbuf, pOutFile);

	}

	
	/* CellPtr is a 2-D array of size 2*CellCount 
	 * and containd cell states used in every simulation 
	 * step accordingly
	 */	
	cellPtr = malloc(2*sizeof(cellState *));
	if(cellPtr==NULL){
		printf("Error: Couldn't malloc for cellPtr\n");
		exit(EXIT_MALLOC);
	}
	
	cellPtr[0] = malloc(cellCount*sizeof(cellState));
	cellPtr[1] = malloc(cellCount*sizeof(cellState));
	if ((!cellPtr[0])||(!cellPtr[1])) {
		printf("Error: Couldn't malloc the array for cellStates\n");
		exit(EXIT_MALLOC);
	}



	/* cellCompParams struct is used to update a cell state.
	 * It contains the current flowing to this specific cell
	 * voltages from communicating cells , and pointers to the previous
	 * cell state and to the next cell state
	 * We allocate cellCount cellCompParams for all the core's cells
	 * Pointer to previous and next states are pointers to CellPtr[0] and CellPtr[1]
	 * elements.
	 */
	cellParamsPtr = malloc(cellCount*sizeof(cellCompParams));
	if(cellParamsPtr==NULL){
		printf("Error: Couldn't malloc for cellParamsPtr\n");
		exit(EXIT_MALLOC);
	}


	for (i=0;i<cellCount;i++)					//initial amount of neighbours for each cell is 0 (bug fix in case the cell stays isolated)
		cellParamsPtr[i].total_amount_of_neighbours = 0;

	/* Make_Core_Communication_List creates a list that details information that this

	 * core needs to exchange with every other core employed.
	 * !!!  The function also takes care of important buffers, such as
	 * the conductance buffer for every cell  !!!
	 */

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	sending_list_head = Make_Core_Communication_List_newest_format(conFile, cellParamsPtr);

	//communication_list_head = Make_Core_Communication_List_new_format(conFile, core_id, cellParamsPtr);

	/* Reckon_phase is an initial communicational process that exchanges information
	 * between cores which describes which cells each core will send to each other core.
	 * We can do this only once in the beginning since this scheme does not change -for now-
	 * to avoid exchanging unnecessary information during the simulation
	 */ 
	
	MPI_Barrier(MPI_COMM_WORLD);
	receiving_list_head = reckon_phase(sending_list_head, cellParamsPtr);
        MPI_Barrier(MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);


   /* Restart procedure, performs alternative variable initializations, recovers
     * output files and jumps inside the appropriate loop.

	 */
	if (RESTARTING) {
		int first_file, last_file, cur_cell, range, files, \
                    min_max_step, local_max, arranged_step, tempint, k;
        mod_prec templong;
		
        FILE ** ckptFiles;
        FILE ** outFiles;
        FILE * new_outfile;
        int ** prev_simSteps;
        cellState *cellbuff;
               
        size_t snapshot_size = sizeof(cellState) + sizeof(int);

        now = get_timestamp();

        /* DEBUG */
        if (core_id == 0)
            printf("prev_cores: %d, prev_cellCount: %d\n", prev_cores, prev_cellCount);

		if ((prev_cores * prev_cellCount) != (cores * cellCount)) {
			printf("Error: Attempting to restart from a checkpoint with a different network size\n");
            exit(EXIT_ARGUMENTS);
		}

        /* Initializations */
		first_file = (core_id * cellCount) / prev_cellCount;
		last_file = ((core_id + 1) * cellCount - 1) / prev_cellCount;
        files = last_file - first_file + 1;
        ckptFiles = (FILE **) malloc(files * sizeof(FILE *));
        outFiles = ckptFiles;
        prev_simSteps = (int **) malloc(2 * sizeof(int *));
        prev_simSteps[0] = (int *) malloc(files * sizeof(int));
        prev_simSteps[1] = (int *) malloc(files * sizeof(int));

        /* Open all required checkpoint files for this core and retrieve the simSteps */
        MPI_Barrier(MPI_COMM_WORLD);
        for (j = 0; j < files; j++) {
            if ((j + first_file) != 0) {
                sprintf(ckptFileName,"/home/apostolis/Desktop/mpi_programs/ckpt_infoli/output/ckptFile%d.bin", j + first_file);
                if ((j + first_file) == core_id)
                    ckptFiles[j] = ckptfd;
                else {
                    ckptFiles[j] = fopen(ckptFileName, "r+");
                    if (ckptFiles[j] == NULL) {
                        printf("Error: Couldn't open %s\n", ckptFileName);
                        exit(EXIT_OPEN_FAIL);
                    }
                }
            } else {
                ckptFiles[j] = zero_ckptfd;
            }
            flock(fileno(ckptFiles[j]), LOCK_SH);
            fseek(ckptFiles[j], 2 * sizeof(int), SEEK_SET);
            fread(&(prev_simSteps[0][j]), sizeof(int), 1, ckptFiles[j]);

            fseek(ckptFiles[j], (3 * sizeof(int)) + prev_cellCount * snapshot_size, SEEK_SET);
            fread(&(prev_simSteps[1][j]), sizeof(int), 1, ckptFiles[j]);
            flock(fileno(ckptFiles[j]), LOCK_UN);
            /* DEBUG */
            printf("core: %d, file: %d simstep1: %d, simpstep2: %d\n", core_id, (j+first_file), \
                    prev_simSteps[0][j], prev_simSteps[1][j]);
        }

        /* Find the largest recoverable simStep for this core.
         * First, find the max simStep for each open file
         * then, take the minimum of those max steps */
        for (j=0; j < files; j++) {
            if (prev_simSteps[0][j] >= prev_simSteps[1][j])
                local_max = prev_simSteps[0][j];
            else
                local_max = prev_simSteps[1][j];
            if ((j == 0) || (local_max < min_max_step))
                min_max_step = local_max;
        }


        /* Determine the max recoverable simStep for all cores 
         * and communicate it to all other cores.
         * If another core's max simStep is smaller, switch to his step */
        int step_buf;
        for (i=0;i<cores;i++) {
            if (core_id == i)
                memcpy(&step_buf, &min_max_step, sizeof(int));
			MPI_Bcast(&step_buf,1,MPI_INT,i,MPI_COMM_WORLD);
            //RCCE_bcast((char *) &step_buf, sizeof(int), i, RCCE_COMM_WORLD);
            if (step_buf < min_max_step)
                min_max_step = step_buf;
        }

		/* set the arranged simStep */
        simStep = min_max_step;

        /* Recover the cellStates and Params for the arranged simStep */
        cur_cell = core_id * cellCount;
        for (j = 0; j < files; j++) {
            range = min(cellCount - (cur_cell % cellCount), \
                        prev_cellCount - (cur_cell % prev_cellCount));
            flock(fileno(ckptFiles[j]), LOCK_SH);
            if (prev_simSteps[0][j] == simStep) {
                fseek(ckptFiles[j], \
                    (3*sizeof(int)) + (cur_cell%prev_cellCount) * snapshot_size, SEEK_SET); 
            } else {
                fseek(ckptFiles[j], 4 * sizeof(int) + \
                    ((cur_cell % prev_cellCount) + prev_cellCount) * snapshot_size, SEEK_SET);
            }
            for (k=0; k<range; k++) {
                fread(&(cellPtr[(simStep%2)][k + (cur_cell % cellCount)]), sizeof(cellState), 1, \
                        ckptFiles[j]);
                fread(&step_buf, sizeof(int), 1, ckptFiles[j]);
                cellParamsPtr[k + (cur_cell % cellCount)].index_of_neighVdend = step_buf;
            }
            flock(fileno(ckptFiles[j]), LOCK_UN);
            cur_cell += range;
            if (ckptFiles[j] != ckptfd)
                fclose(ckptFiles[j]);
        }
        
        /* Initialize the constants for the next state */
        initState(cellPtr[(simStep%2)^1]);

        /* DEBUG Message: cells recovered per core */
        MPI_Barrier(MPI_COMM_WORLD);
        for (j=0; j<cellCount; j++) {
            printf("core id: %d recovered cell: %d\n", core_id, cellPtr[simStep%2][j].cellID, cellPtr[simStep%2][j].cell_x, cellPtr[simStep%2][j].cell_y);
        }

        if (cellCount != prev_cellCount) {
            /* Create the new output file */
            sprintf(outFileName,"/home/apostolis/Desktop/mpi_programs/ckpt_infoli/output/new_InferiorOlive_Output%d.txt", core_id);
            new_outfile = fopen(outFileName, "w");

            /* Open the required previous output files */
            for (j = 0; j < files; j++) {
                sprintf(outFileName, "/home/apostolis/Desktop/mpi_programs/ckpt_infoli/output/InferiorOlive_Output%d.txt",\
                        j + first_file);
                printf("core: %d, opening output file %d\n",core_id, (j + first_file));
                outFiles[j] = fopen(outFileName, "r");
                if (outFiles[j] == NULL) {
                    printf("Error: Couldn't open %s\n", ckptFileName);
                    exit(EXIT_OPEN_FAIL);
                }

                flock(fileno(outFiles[j]), LOCK_SH);
                fgets(tempbuf, OUTFILE_MAXLINE, outFiles[j]); // skip the header
                flock(fileno(outFiles[j]), LOCK_UN);
            }

            /* Write the new output file's header */
            sprintf(tempbuf, "#simSteps Time(ms) Input(Iapp) Output(V_soma)");
            fputs(tempbuf, new_outfile);
            for (i=0; i<cellCount; i++) {
                sprintf(tempbuf, "[%d][%d] ", \
                        cellPtr[simStep%2][i].cell_x, cellPtr[simStep%2][i].cell_y);
                fputs(tempbuf, new_outfile);
            }
            sprintf(tempbuf, "\n");
            fputs(tempbuf, new_outfile);

            /* Populate the new output file */
            int leftover;
            for (i=0; i<simStep; i++) {
                /* Write simStep, time and input for this simStep */
                iApp = ((i>20000-1)&&(i<20500-1)) ? 6 : 0; 
                if (inputFromFile)
                    sprintf(tempbuf, "%d %.2f ", i + 1, i*0.05); 
                else
                    sprintf(tempbuf, "%d %.2f %.1f ", i + 1, i*0.05, iApp); 
                fputs(tempbuf, new_outfile);

                /* Read the cell outputs from the previous output files */
                cur_cell = core_id * cellCount;
                for (j=0; j<files; j++) {
                    // scroll the file past the simstep, time and input
                    flock(fileno(outFiles[j]), LOCK_SH);
                    fscanf(outFiles[j], "%d", &tempint);
                    fscanf(outFiles[j], "%lf", &templong);
                    if (!inputFromFile)
                        fscanf(outFiles[j], "%lf", &templong);

                    // skip to the first element needed from this file
                    for (k=0; k<(cur_cell % prev_cellCount); k++)
                        fscanf(outFiles[j], "%lf ", &templong);

                    range = min(cellCount - (cur_cell % cellCount), \
                            prev_cellCount - (cur_cell % prev_cellCount));
                    
                    for (k=0; k<range; k++) {
                        fscanf(outFiles[j], "%lf ", &templong);
                        sprintf(tempbuf,"%.8f ", templong);
                        fputs(tempbuf, new_outfile);
                    }

                    // skip to the next line if there are leftover elements
                    leftover = prev_cellCount - (range + cur_cell % prev_cellCount);
                    for (k=0; k < leftover; k++) {
                        fscanf(outFiles[j], "%lf ", &templong);
                    }
                    flock(fileno(outFiles[j]), LOCK_UN);
                    cur_cell += range;
                }
                sprintf(tempbuf, "\n"); 
                fputs(tempbuf, new_outfile);
            }

            /* close the new file */
            MPI_Barrier(MPI_COMM_WORLD);
            fclose(new_outfile);

            /* Replace the old output file with the new one */
            sprintf(tempbuf, \
                "mv /home/apostolis/Desktop/mpi_programs/ckpt_infoli/output/new_InferiorOlive_Output%d.txt /home/apostolis/Desktop/mpi_programs/ckpt_infoli/output/InferiorOlive_Output%d.txt", core_id, core_id);
            system(tempbuf);
           //chmod(outFileName, 0x77777);
            chmod(outFileName,S_IRUSR | S_IWUSR);
            sprintf(tempbuf, "chown apostolis %s", ckptFileName);
            system(tempbuf);
            MPI_Barrier(MPI_COMM_WORLD);
        }

        /* Open the output file to write on */
        sprintf(outFileName,"/home/apostolis/Desktop/mpi_programs/ckpt_infoli/output/InferiorOlive_Output%d.txt", core_id);
        pOutFile = fopen(outFileName, "a");

        /* Truncate this core's checkpoint file 
         * After this point, the simulation is unrecoverable for two CKPT_INTERVALs */
        ckptfd = freopen(NULL, "w", ckptfd);
		fwrite(&cores, sizeof(int), 1, ckptfd);
		fwrite(&cellCount, sizeof(int), 1, ckptfd);

        if (core_id == 0)
            printf("Restarting from simStep %d\n", simStep);

		fsync(fileno(stdout));

        /* Free restart-only data structures */
        free(cellbuff);

		/* jump in the appropriate loop */
        MPI_Barrier(MPI_COMM_WORLD);

        timediff = get_timestamp() - now;
        if (core_id==0){
            restartTimeFile=fopen(restartFileName,"w");
            fprintf(restartTimeFile,"%lf\n",(double)timediff/1000000);
            fclose(restartTimeFile);
            printf("Restart Time: %lf\n", (double)timediff/1000000);
        }
       

	if (inputFromFile) {
            mod_prec tmpnum;
            int counter = 1;
            while (counter < simStep) {
                ReadFileLine(pInFile, iAppArray);
                counter++;
            }
			goto restartFile;
		} else {
			simTime = SIMTIME; 
			total_simulation_steps = ceil(simTime/DELTA);
			goto restartRandom;
		}

	}

	/* Initialise cellPtr[0] with appropriate values.
	 */
	initState(cellPtr[0]);
	
	if (PRINTING) {
		for (i=0;i<cellCount;i++) {
			sprintf(tempbuf, "[%d][%d] ", cellPtr[0][i].cell_x, cellPtr[0][i].cell_y);
			fputs(tempbuf, pOutFile);
		}
		sprintf(tempbuf, "\n");
		fputs(tempbuf, pOutFile);
	}
	
	//Initialize g_CaL
	seedvar = time(NULL)+core_id;		//seedvar will be time- and core-dependent now!
	srand(seedvar);

	for(i=0;i<cellCount;i++){

		cellPtr[1][i].soma.g_CaL = cellPtr[0][i].soma.g_CaL;

		if (RAND_INIT) {
			cellPtr[0][i].soma.g_CaL = 0.6+(0.2*(rand()%100)/100);
			cellPtr[1][i].soma.g_CaL = cellPtr[0][i].soma.g_CaL;
		}
	}

	//random initialization process
	if (RAND_INIT) {
		for(i=0;i<cellCount;i++) {
			initSteps = rand()%(int)ceil(100/DELTA);
			initSteps = initSteps | 0x00000001;//make it odd, so that the final state is in prevCellState

			for(j=0;j<initSteps;j++){
				//Arrange inputs
				cellParamsPtr[i].iAppIn = 0;//No stimulus
				cellParamsPtr[i].prevCellState = &cellPtr[j%2][i];
				cellParamsPtr[i].newCellState = &cellPtr[(j%2)^1][i];

				CompDend(&cellParamsPtr[i], 1);
				CompSoma(&cellParamsPtr[i]);
				CompAxon(&cellParamsPtr[i]);
			}

		}
	}
	

	/* start of the simulation
	 * In case we want to read the stimulus from file inputFromFile = true
	 */

	if(inputFromFile){
		simStep = 0;

		/* Read full lines until end of file. 
		 * Every iteration (line) is one simulation step.
		 */
		while(ReadFileLine(pInFile, iAppArray)) {
restartFile:
			simulation_array_ID = simStep%2;

            if (simStep % SYNC_INTERVAL == 0) {
                fflush(pOutFile);
                if (fsync(fileno(pOutFile)) != 0)
					printf("Error: Outfile file sync failed\n");
                printf("%d\n", simStep);
            }

			if ((CHECKPOINTING) && (simStep % CKPT_INTERVAL == 0)) {
                MPI_Barrier(MPI_COMM_WORLD);
                gettimeofday(&tts,NULL);
                take_checkpoint(simStep, cellPtr, cellParamsPtr, cellCount, ckptfd, ckptBuffer);
                MPI_Barrier(MPI_COMM_WORLD);
                gettimeofday(&ttf,NULL); 
             if (core_id==0){
                ckptTimeFile=fopen(ckptTimeFileName,"w");
                fprintf(ckptTimeFile,"%lf\n",(ttf.tv_sec-tts.tv_sec)+(ttf.tv_usec-tts.tv_usec)*0.000001);
                fclose(ckptTimeFile);
                //printf("ckptTime= %lf ",(ttf.tv_sec-tts.tv_sec)+(ttf.tv_usec-tts.tv_usec)*0.000001);
             }
            
            }

			if (PRINTING) {
				sprintf(tempbuf, "%d %.2f ", (simStep+1), simStep*0.05); // start @ 1 because skipping initial values
				fputs(tempbuf, pOutFile);
			}
			
			/* Perform_Communication() performs the inter core
  			 * core dendrite communication with neighbouring cells
			 * See definition for more details
  			 */

			perform_communication_step(sending_list_head, receiving_list_head, cellParamsPtr, cellPtr[simulation_array_ID]);
			//Perform_Communication_Version_2(communication_list_head, cellParamsPtr, cellPtr, simulation_array_ID);

			for (target_cell=0;target_cell<cellCount;target_cell++) {

				cellParamsPtr[target_cell].iAppIn = iAppArray[target_cell];
				cellParamsPtr[target_cell].prevCellState = &cellPtr[simulation_array_ID][target_cell];
				cellParamsPtr[target_cell].newCellState = &cellPtr[simulation_array_ID^1][target_cell];

				CompDend(&cellParamsPtr[target_cell], 0);

				CompSoma(&cellParamsPtr[target_cell]);
				CompAxon(&cellParamsPtr[target_cell]);

				if (PRINTING) {
					sprintf(tempbuf, "%.16f ", cellPtr[(simulation_array_ID)^1][target_cell].axon.V_axon);
					fputs(tempbuf, pOutFile);
				}

			}
			if (PRINTING) {
				sprintf(tempbuf, "\n");
				fputs(tempbuf,pOutFile);
			}

			simStep++;
		}
		
		if (core_id==0)
			printf("End of simulation reached.\n");
	}

	else {	
		simTime = SIMTIME; 
		total_simulation_steps = (int)(simTime/DELTA);

		now = get_timestamp();		
		for(simStep=0;simStep<total_simulation_steps;simStep++) {
								
restartRandom:
			simulation_array_ID = simStep%2;

            if (simStep % SYNC_INTERVAL == 0) {
                fflush(pOutFile);
                if (fsync(fileno(pOutFile)) != 0)
					printf("Error: Outfile file sync failed\n");
            }

			if ((CHECKPOINTING) && (simStep % CKPT_INTERVAL == 0)) {

                MPI_Barrier(MPI_COMM_WORLD);
                gettimeofday(&tts,NULL);
                take_checkpoint(simStep, cellPtr, cellParamsPtr, cellCount, ckptfd, ckptBuffer);
		        MPI_Barrier(MPI_COMM_WORLD);
                gettimeofday(&ttf,NULL);	
             if (core_id==0){
                ckptTimeFile=fopen(ckptTimeFileName,"w");
                fprintf(ckptTimeFile,"%lf\n",(ttf.tv_sec-tts.tv_sec)+(ttf.tv_usec-tts.tv_usec)*0.000001);
                fclose(ckptTimeFile);
                //printf("ckptTime= %lf ",(ttf.tv_sec-tts.tv_sec)+(ttf.tv_usec-tts.tv_usec)*0.000001);
             }
 
                fflush(stdout);
                counter++;
                now = get_timestamp();
            }

			if ((simStep>=20000)&&(simStep<20500-1))
			//if ((simStep>=20000)&&((simStep%20000)>=0)&&((simStep%20000)<500))		//re-occuring spikes after first sec
				iApp = 6; 
			else
		   		iApp = 0;
			
			if (PRINTING) {
				sprintf(tempbuf, " %d %.2f ", simStep+1, iApp);
				fputs(tempbuf, pOutFile);
			}
			
			/* Perform_Communication() performs the inter core
  			 * core dendrite communication with neighbouring cells
			 * See definition for more details
  			 */
			
			//Perform_Communication(communication_list_head, cellParamsPtr, cellPtr, simulation_array_ID);

			
			perform_communication_step(sending_list_head, receiving_list_head, cellParamsPtr, cellPtr[simulation_array_ID]);
			//printf("OK1\n");
			//Perform_Communication_Version_2(communication_list_head, cellParamsPtr, cellPtr, simulation_array_ID);
			//printf("OK2\n");
			for (target_cell=0;target_cell<cellCount;target_cell++) {
				/* we simulate a hardcoded input pulse here 
				 * that differs from step to step 
				 */

				cellParamsPtr[target_cell].iAppIn = iApp;			
				cellParamsPtr[target_cell].prevCellState = &cellPtr[simulation_array_ID][target_cell];
				cellParamsPtr[target_cell].newCellState = &cellPtr[simulation_array_ID^1][target_cell];
	
				CompDend(&cellParamsPtr[target_cell], 0);
				CompSoma(&cellParamsPtr[target_cell]);
				CompAxon(&cellParamsPtr[target_cell]);
				
				if (PRINTING) {
					sprintf(tempbuf, " %d : %.8f ", (core_id*cellCount+target_cell+1), cellPtr[(simulation_array_ID)^1][target_cell].axon.V_axon);
					fputs(tempbuf, pOutFile);
				}

			}
			
			if (PRINTING) {
				sprintf(tempbuf, "\n");
				fputs(tempbuf, pOutFile);
			}

		}

		simStep = total_simulation_steps;		//so that simStep in the end has the exact value of how many steps we had in this sim, regardless of input method (useful to know which cellPtr has what)
	}
		
	/* Free  memory and close files
	 */

	if (PRINTSTATE) {
		char resultFileName[50];
		FILE *resultFile;
		sprintf(resultFileName, "/home/apostolis/Desktop/mpi_programs/ckpt_infoli/results/lastStateDump%d.txt", core_id);
		resultFile=fopen(resultFileName, "w+");
		sprintf(tempbuf, " %d : %.8f ", (core_id*cellCount+target_cell+1), cellPtr[(simulation_array_ID)^1][target_cell].axon.V_axon);
		fputs(tempbuf, resultFile);
		
		//printState(cellPtr[simStep%2], resultFileName);		//simStep%2 here should refer to the cellPtr which has the last state of the network that we calculated
	}

	free(cellPtr[0]);
	free(cellPtr[1]);
	free(cellPtr);
	free(cellParamsPtr);
	/* reminder, maybe free communication list */
	
	
	if (PRINTING) {
		fclose (pOutFile);
		chmod(outFileName,0x01B6);
	}

	if(inputFromFile)
		fclose (pInFile);


	end = clock() - start;
	if (core_id==0) {
		printf("Total time:\t%.3f secs\n", ((double) end)/((double)CLOCKS_PER_SEC));
//		printf("Input:\t\t%d usecs\n", (intime.tv_sec)*1000000+intime.tv_usec);
//		printf("Output:\t\t%d usecs\n", (outtime.tv_sec)*1000000+outtime.tv_usec);
//		printf("Processing:\t%d usecs\n", (computingtime.tv_sec+messagingtime.tv_sec)*1000000+computingtime.tv_usec+messagingtime.tv_sec);
//		printf("Split in:\n");
//		printf("Computations:\t%d usecs\n", (computingtime.tv_sec)*1000000+computingtime.tv_usec);
//		printf("Messaging:\t%d usecs\n", (messagingtime.tv_sec)*1000000+messagingtime.tv_usec);
	}
	MPI_Finalize();	
	return 0;
}




//DENDRITIC COMPUTATIONAL PART ------------------

void CompDend(cellCompParams *cellParamsPtr, int randomness){

	struct channelParams chPrms;
	struct dendCurrVoltPrms chComps;

	//printf("Dendrite ");

	//Prepare pointers to inputs/outputs
	chPrms.v = &cellParamsPtr->prevCellState->dend.V_dend;
	chPrms.prevComp1 = &cellParamsPtr->prevCellState->dend.Hcurrent_q;
	chPrms.newComp1 = &cellParamsPtr->newCellState->dend.Hcurrent_q;
	//Compute
	DendHCurr(&chPrms);

	//Prepare pointers to inputs/outputs
	chPrms.v = &cellParamsPtr->prevCellState->dend.V_dend;
	chPrms.prevComp1 = &cellParamsPtr->prevCellState->dend.Calcium_r;
	chPrms.newComp1 = &cellParamsPtr->newCellState->dend.Calcium_r;
	//Compute
	DendCaCurr(&chPrms);

	//Prepare pointers to inputs/outputs
	chPrms.prevComp1 = &cellParamsPtr->prevCellState->dend.Potassium_s;
	chPrms.prevComp2 = &cellParamsPtr->prevCellState->dend.Ca2Plus;
	chPrms.newComp1 = &cellParamsPtr->newCellState->dend.Potassium_s;
	//Compute
	DendKCurr(&chPrms);

	//Prepare pointers to inputs/outputs
	chPrms.prevComp1 = &cellParamsPtr->prevCellState->dend.Ca2Plus;
	chPrms.prevComp2 = &cellParamsPtr->prevCellState->dend.I_CaH;
	chPrms.newComp1 = &cellParamsPtr->newCellState->dend.Ca2Plus;
	//Compute
	DendCal(&chPrms);

	/* ANDREAS change here the last parameter of IcNeighbors, we no longer use number_of_neighbours 
	 * instead we have an index in struct cellParameters which tells us where the cell's neighbouring
	 * voltages end in the array
	 */


	//in random initialization mode, cells run in closed circuit mode so neighboring current equals zero
	if (randomness==1)
		chComps.iC = 0;
	else
		chComps.iC = IcNeighbors(cellParamsPtr->neighVdend, cellParamsPtr->neighConductances, cellParamsPtr->prevCellState->dend.V_dend, cellParamsPtr->total_amount_of_neighbours);

	
	chComps.iApp = &cellParamsPtr->iAppIn;
	chComps.vDend = &cellParamsPtr->prevCellState->dend.V_dend;
	chComps.newVDend = &cellParamsPtr->newCellState->dend.V_dend;
	chComps.vSoma = &cellParamsPtr->prevCellState->soma.V_soma;
	chComps.q = &cellParamsPtr->newCellState->dend.Hcurrent_q;
	chComps.r = &cellParamsPtr->newCellState->dend.Calcium_r;
	chComps.s = &cellParamsPtr->newCellState->dend.Potassium_s;
	chComps.newI_CaH = &cellParamsPtr->newCellState->dend.I_CaH;
	DendCurrVolt(&chComps);

	return;
}

void DendHCurr(struct channelParams *chPrms){

	mod_prec q_inf, tau_q, dq_dt, q_local;

	//Get inputs
	mod_prec prevV_dend = *chPrms->v;
	mod_prec prevHcurrent_q = *chPrms->prevComp1;

	// Update dendritic H current component
	q_inf = 1 /(1 + exp((prevV_dend + 80) / 4));
	tau_q = 1 /(exp(-0.086 * prevV_dend - 14.6) + exp(0.070 * prevV_dend - 1.87));
	dq_dt = (q_inf - prevHcurrent_q) / tau_q;
	q_local = DELTA * dq_dt + prevHcurrent_q;
	//Put result
	*chPrms->newComp1 = q_local;

	return;
}

void DendCaCurr(struct channelParams *chPrms){

	mod_prec alpha_r, beta_r, r_inf, tau_r, dr_dt, r_local;

	//Get inputs
	mod_prec prevV_dend = *chPrms->v;
	mod_prec prevCalcium_r = *chPrms->prevComp1;

	// Update dendritic high-threshold Ca current component
	alpha_r = 1.7 / (1 + exp( -(prevV_dend - 5) / 13.9));
	beta_r = 0.02 * (prevV_dend + 8.5) / (exp((prevV_dend + 8.5) / 5) - 1);
	r_inf = alpha_r / (alpha_r + beta_r);
	tau_r = 5 / (alpha_r + beta_r);
	dr_dt = (r_inf - prevCalcium_r) / tau_r;
	r_local = DELTA * dr_dt + prevCalcium_r;
	//Put result
	*chPrms->newComp1 = r_local;

	return;
}

void DendKCurr(struct channelParams *chPrms){

	mod_prec alpha_s, beta_s, s_inf, tau_s, ds_dt, s_local;

	//Get inputs
	mod_prec prevPotassium_s = *chPrms->prevComp1;
	mod_prec prevCa2Plus = *chPrms->prevComp2;

	// Update dendritic Ca-dependent K current component
	alpha_s = min((0.00002*prevCa2Plus), 0.01);
	beta_s = 0.015;
	s_inf = alpha_s / (alpha_s + beta_s);
	tau_s = 1 / (alpha_s + beta_s);
	ds_dt = (s_inf - prevPotassium_s) / tau_s;
	s_local = DELTA * ds_dt + prevPotassium_s;
	//Put result
	*chPrms->newComp1 = s_local;

	return;
}

//Consider merging DendCal into DendKCurr since DendCal's output doesn't go to DendCurrVolt but to DendKCurr
void DendCal(struct channelParams *chPrms){

	mod_prec dCa_dt, Ca2Plus_local;

	//Get inputs
	mod_prec prevCa2Plus = *chPrms->prevComp1;
	mod_prec prevI_CaH = *chPrms->prevComp2;

	// update Calcium concentration
	dCa_dt = -3 * prevI_CaH - 0.075 * prevCa2Plus;
	Ca2Plus_local = DELTA * dCa_dt + prevCa2Plus;
	//Put result
	*chPrms->newComp1 = Ca2Plus_local;//This state value is read in DendKCurr

	return;
}

void DendCurrVolt(struct dendCurrVoltPrms *chComps){

	//Loca variables
	mod_prec I_sd, I_CaH, I_K_Ca, I_ld, I_h, dVd_dt;

	//Get inputs
	mod_prec I_c = chComps->iC;
	mod_prec I_app = *chComps->iApp;
	mod_prec prevV_dend = *chComps->vDend;
	mod_prec prevV_soma = *chComps->vSoma;
	mod_prec q = *chComps->q;
	mod_prec r = *chComps->r;
	mod_prec s = *chComps->s;

	// DENDRITIC CURRENTS

	// Soma-dendrite interaction current I_sd
	I_sd = (G_INT / (1 - P1)) * (prevV_dend - prevV_soma);
	// Inward high-threshold Ca current I_CaH
	I_CaH=G_CAH * r * r * (prevV_dend - V_CA);
	// Outward Ca-dependent K current I_K_Ca
	I_K_Ca =G_K_CA * s * (prevV_dend - V_K);
	// Leakage current I_ld
	I_ld =G_LD * (prevV_dend - V_L);
	// Inward anomalous rectifier I_h
	I_h=G_H * q * (prevV_dend - V_H);

	dVd_dt = (-(I_CaH + I_sd+ I_ld + I_K_Ca + I_c + I_h) + I_app) / C_M;

	//Put result (update V_dend)
	*chComps->newVDend = DELTA * dVd_dt + prevV_dend;
	*chComps->newI_CaH = I_CaH;//This is a state value read in DendCal
	return;
}

mod_prec IcNeighbors(mod_prec *neighVdend, mod_prec *neighConductances, mod_prec prevV_dend, int neighbors){

	int i;
	mod_prec f, V, Cond, I_c=0;






	for(i=0;i<neighbors;i++){
		V = prevV_dend - neighVdend[i];
		f = 0.8 * exp(-1*pow(V, 2)/100) + 0.2;// SCHWEIGHOFER 2004 VERSION
		Cond = neighConductances[i];
		I_c = I_c + (Cond * f * V);
//		V = 1+1;
//		f = 3+3;
	}

	return I_c;
}

//SOMATIC COMPUTATIONAL PART -----------------

void CompSoma(cellCompParams *cellParamsPtr){

	struct channelParams chPrms;
	struct somaCurrVoltPrms chComps;

	// update somatic components
	// SCHWEIGHOFER:

	//Prepare pointers to inputs/outputs
	chPrms.v = &cellParamsPtr->prevCellState->soma.V_soma;
	chPrms.prevComp1 = &cellParamsPtr->prevCellState->soma.Calcium_k;
	chPrms.prevComp2 = &cellParamsPtr->prevCellState->soma.Calcium_l;
	chPrms.newComp1 = &cellParamsPtr->newCellState->soma.Calcium_k;
	chPrms.newComp2 = &cellParamsPtr->newCellState->soma.Calcium_l;
	//Compute
	SomaCalcium(&chPrms);

	//Prepare pointers to inputs/outputs
	chPrms.v = &cellParamsPtr->prevCellState->soma.V_soma;
	chPrms.prevComp1 = &cellParamsPtr->prevCellState->soma.Sodium_m;
	chPrms.prevComp2 = &cellParamsPtr->prevCellState->soma.Sodium_h;
	chPrms.newComp1 = &cellParamsPtr->newCellState->soma.Sodium_m;
	chPrms.newComp2 = &cellParamsPtr->newCellState->soma.Sodium_h;
	//Compute
	SomaSodium(&chPrms);

	//Prepare pointers to inputs/outputs
	chPrms.v = &cellParamsPtr->prevCellState->soma.V_soma;
	chPrms.prevComp1 = &cellParamsPtr->prevCellState->soma.Potassium_n;
	chPrms.prevComp2 = &cellParamsPtr->prevCellState->soma.Potassium_p;
	chPrms.newComp1 = &cellParamsPtr->newCellState->soma.Potassium_n;
	chPrms.newComp2 = &cellParamsPtr->newCellState->soma.Potassium_p;
	//Compute
	SomaPotassium(&chPrms);

	//Prepare pointers to inputs/outputs
	chPrms.v = &cellParamsPtr->prevCellState->soma.V_soma;
	chPrms.prevComp1 = &cellParamsPtr->prevCellState->soma.Potassium_x_s;
	chPrms.newComp1 = &cellParamsPtr->newCellState->soma.Potassium_x_s;
	//Compute
	SomaPotassiumX(&chPrms);

	chComps.g_CaL = &cellParamsPtr->prevCellState->soma.g_CaL;
	chComps.vDend = &cellParamsPtr->prevCellState->dend.V_dend;
	chComps.vSoma = &cellParamsPtr->prevCellState->soma.V_soma;
	chComps.newVSoma = &cellParamsPtr->newCellState->soma.V_soma;
	chComps.vAxon = &cellParamsPtr->prevCellState->axon.V_axon;
	chComps.k = &cellParamsPtr->newCellState->soma.Calcium_k;
	chComps.l = &cellParamsPtr->newCellState->soma.Calcium_l;
	chComps.m = &cellParamsPtr->newCellState->soma.Sodium_m;
	chComps.h = &cellParamsPtr->newCellState->soma.Sodium_h;
	chComps.n = &cellParamsPtr->newCellState->soma.Potassium_n;
	chComps.x_s = &cellParamsPtr->newCellState->soma.Potassium_x_s;
	SomaCurrVolt(&chComps);

	return;
}

void SomaCalcium(struct channelParams *chPrms){

	mod_prec k_inf, l_inf, tau_k, tau_l, dk_dt, dl_dt, k_local, l_local;

	//Get inputs
	mod_prec prevV_soma = *chPrms->v;
	mod_prec prevCalcium_k = *chPrms->prevComp1;
	mod_prec prevCalcium_l = *chPrms->prevComp2;

	k_inf = (1 / (1 + exp(-1 * (prevV_soma + 61) / 4.2)));
	l_inf = (1 / (1 + exp(( prevV_soma + 85.5) / 8.5)));
	tau_k = 1;
	tau_l = ((20 * exp((prevV_soma + 160) / 30) / (1 + exp((prevV_soma + 84) / 7.3))) +35);
	dk_dt = (k_inf - prevCalcium_k) / tau_k;
	dl_dt = (l_inf - prevCalcium_l) / tau_l;
	k_local = DELTA * dk_dt + prevCalcium_k;
	l_local = DELTA * dl_dt + prevCalcium_l;
	//Put result
	*chPrms->newComp1= k_local;
	*chPrms->newComp2= l_local;

	return;
}

void SomaSodium(struct channelParams *chPrms){

	mod_prec m_inf, h_inf, tau_h, dh_dt, m_local, h_local;

	//Get inputs
	mod_prec prevV_soma = *chPrms->v;
	//mod_prec prevSodium_m = *chPrms->prevComp1;
	mod_prec prevSodium_h = *chPrms->prevComp2;

	// RAT THALAMOCORTICAL SODIUM:
	m_inf = 1 / (1 + (exp((-30 - prevV_soma)/ 5.5)));
	h_inf = 1 / (1 + (exp((-70 - prevV_soma)/-5.8)));
	tau_h = 3 * exp((-40 - prevV_soma)/33);
	dh_dt = (h_inf - prevSodium_h)/tau_h;
	m_local = m_inf;
	h_local = prevSodium_h + DELTA * dh_dt;
	//Put result
	*chPrms->newComp1 = m_local;
	*chPrms->newComp2 = h_local;

	return;
}

void SomaPotassium(struct channelParams *chPrms){

	mod_prec n_inf, p_inf, tau_n, tau_p, dn_dt, dp_dt, n_local, p_local;

	//Get inputs
	mod_prec prevV_soma = *chPrms->v;
	mod_prec prevPotassium_n = *chPrms->prevComp1;
	mod_prec prevPotassium_p = *chPrms->prevComp2;

	// NEOCORTICAL
	n_inf = 1 / (1 + exp( ( -3 - prevV_soma) /10));
	p_inf = 1 / (1 + exp( (-51 - prevV_soma) / -12));
	tau_n = 5 + (47 * exp( -(-50 - prevV_soma) /900));
	tau_p = tau_n;
	dn_dt = (n_inf - prevPotassium_n) / tau_n;
	dp_dt = (p_inf - prevPotassium_p) / tau_p;
	n_local = DELTA * dn_dt + prevPotassium_n;
	p_local = DELTA * dp_dt + prevPotassium_p;
	//Put result
	*chPrms->newComp1 = n_local;
	*chPrms->newComp2 = p_local;

	return;
}

void SomaPotassiumX(struct channelParams *chPrms){

	mod_prec alpha_x_s, beta_x_s, x_inf_s, tau_x_s, dx_dt_s, x_s_local;

	//Get inputs
	mod_prec prevV_soma = *chPrms->v;
	mod_prec prevPotassium_x_s = *chPrms->prevComp1;

	// Voltage-dependent (fast) potassium
	alpha_x_s = 0.13 * (prevV_soma + 25) / (1 - exp(-(prevV_soma + 25) / 10));
	beta_x_s= 1.69 * exp(-0.0125 * (prevV_soma + 35));
	x_inf_s = alpha_x_s / (alpha_x_s + beta_x_s);
	tau_x_s = 1 / (alpha_x_s + beta_x_s);
	dx_dt_s = (x_inf_s - prevPotassium_x_s) / tau_x_s;
	x_s_local = 0.05 * dx_dt_s + prevPotassium_x_s;
	//Put result
	*chPrms->newComp1 = x_s_local;

	return;
}

void SomaCurrVolt(struct somaCurrVoltPrms *chComps){

	//Local variables
	mod_prec I_ds, I_CaL, I_Na_s, I_ls, I_Kdr_s, I_K_s, I_as, dVs_dt;

	//Get inputs
	mod_prec g_CaL = *chComps->g_CaL;
	mod_prec prevV_dend = *chComps->vDend;
	mod_prec prevV_soma = *chComps->vSoma;
	mod_prec prevV_axon = *chComps->vAxon;
	mod_prec k = *chComps->k;
	mod_prec l = *chComps->l;
	mod_prec m = *chComps->m;
	mod_prec h = *chComps->h;
	mod_prec n = *chComps->n;
	mod_prec x_s = *chComps->x_s;

	// SOMATIC CURRENTS

	// Dendrite-soma interaction current I_ds
	I_ds= (G_INT / P1) * (prevV_soma - prevV_dend);
	// Inward low-threshold Ca current I_CaL
	I_CaL = g_CaL * k * k * k * l * (prevV_soma - V_CA); //k^3
	// Inward Na current I_Na_s
	I_Na_s= G_NA_S * m * m * m * h * (prevV_soma - V_NA);
	// Leakage current I_ls
	I_ls= G_LS * (prevV_soma - V_L);
	// Outward delayed potassium current I_Kdr
	I_Kdr_s = G_KDR_S * n * n * n * n * (prevV_soma - V_K); // SCHWEIGHOFER
	// I_K_s
	I_K_s = G_K_S * pow(x_s, 4) * (prevV_soma - V_K);
	// Axon-soma interaction current I_as
	I_as= (G_INT / (1 - P2)) * (prevV_soma - prevV_axon);

	dVs_dt = (-(I_CaL + I_ds+ I_as + I_Na_s + I_ls + I_Kdr_s + I_K_s)) / C_M;
	*chComps->newVSoma = DELTA * dVs_dt + prevV_soma;

	return;
}

//AXONAL COMPUTATIONAL PART -------------------

void CompAxon(cellCompParams *cellParamsPtr){

	struct channelParams chPrms;
	struct axonCurrVoltPrms chComps;

	// update somatic components
	// SCHWEIGHOFER:

	//Prepare pointers to inputs/outputs
	chPrms.v = &cellParamsPtr->prevCellState->axon.V_axon;
	chPrms.prevComp1 = &cellParamsPtr->prevCellState->axon.Sodium_h_a;
	chPrms.newComp1 = &cellParamsPtr->newCellState->axon.Sodium_h_a;
	chPrms.newComp2 = &cellParamsPtr->newCellState->axon.Sodium_m_a;
	//Compute
	AxonSodium(&chPrms);

	//Prepare pointers to inputs/outputs
	chPrms.v = &cellParamsPtr->prevCellState->axon.V_axon;
	chPrms.prevComp1 = &cellParamsPtr->prevCellState->axon.Potassium_x_a;
	chPrms.newComp1 = &cellParamsPtr->newCellState->axon.Potassium_x_a;
	//Compute
	AxonPotassium(&chPrms);

	//Get inputs
	chComps.vSoma = &cellParamsPtr->prevCellState->soma.V_soma;
	chComps.vAxon = &cellParamsPtr->prevCellState->axon.V_axon;
	chComps.newVAxon = &cellParamsPtr->newCellState->axon.V_axon;
	chComps.m_a = &cellParamsPtr->newCellState->axon.Sodium_m_a;
	chComps.h_a = &cellParamsPtr->newCellState->axon.Sodium_h_a;
	chComps.x_a = &cellParamsPtr->newCellState->axon.Potassium_x_a;
	AxonCurrVolt(&chComps);

	return;
}

void AxonSodium(struct channelParams *chPrms){

	mod_prec m_inf_a, h_inf_a, tau_h_a, dh_dt_a, m_a_local, h_a_local;

	//Get inputs
	mod_prec prevV_axon = *chPrms->v;
	mod_prec prevSodium_h_a = *chPrms->prevComp1;

	// Update axonal Na components
	// NOTE: current has shortened inactivation to account for high
	// firing frequencies in axon hillock
	m_inf_a = 1 / (1 + (exp((-30 - prevV_axon)/ 5.5)));
	h_inf_a = 1 / (1 + (exp((-60 - prevV_axon)/-5.8)));
	tau_h_a = 1.5 * exp((-40 - prevV_axon)/33);
	dh_dt_a = (h_inf_a - prevSodium_h_a)/tau_h_a;
	m_a_local = m_inf_a;
	h_a_local = prevSodium_h_a + DELTA * dh_dt_a;
	//Put result
	*chPrms->newComp1 = h_a_local;
	*chPrms->newComp2 = m_a_local;

	return;
}

void AxonPotassium(struct channelParams *chPrms){

	mod_prec alpha_x_a, beta_x_a, x_inf_a, tau_x_a, dx_dt_a, x_a_local;

	//Get inputs
	mod_prec prevV_axon = *chPrms->v;
	mod_prec prevPotassium_x_a = *chPrms->prevComp1;

	// D'ANGELO 2001 -- Voltage-dependent potassium
	alpha_x_a = 0.13 * (prevV_axon + 25) / (1 - exp(-(prevV_axon + 25) / 10));
	beta_x_a= 1.69 * exp(-0.0125 * (prevV_axon + 35));
	x_inf_a = alpha_x_a / (alpha_x_a + beta_x_a);
	tau_x_a = 1 / (alpha_x_a + beta_x_a);
	dx_dt_a = (x_inf_a - prevPotassium_x_a) / tau_x_a;
	x_a_local = 0.05 * dx_dt_a + prevPotassium_x_a;
	//Put result
	*chPrms->newComp1 = x_a_local;

	return;
}

void AxonCurrVolt(struct axonCurrVoltPrms *chComps){

	//Local variable
	mod_prec I_Na_a, I_la, I_sa, I_K_a, dVa_dt;

	//Get inputs
	mod_prec prevV_soma = *chComps->vSoma;
	mod_prec prevV_axon = *chComps->vAxon;
	mod_prec m_a = *chComps->m_a;
	mod_prec h_a = *chComps->h_a;
	mod_prec x_a = *chComps->x_a;

	// AXONAL CURRENTS
	// Sodium
	I_Na_a= G_NA_A* m_a * m_a * m_a * h_a * (prevV_axon - V_NA);
	// Leak
	I_la= G_LA* (prevV_axon - V_L);
	// Soma-axon interaction current I_sa
	I_sa= (G_INT / P2) * (prevV_axon - prevV_soma);
	// Potassium (transient)
	I_K_a = G_K_A * pow(x_a, 4) * (prevV_axon - V_K);
	dVa_dt = (-(I_K_a + I_sa + I_la + I_Na_a)) / C_M;
	*chComps->newVAxon = DELTA * dVa_dt + prevV_axon;

	return;
}

//Initialization Function, important ! --------------------------

void initState(cellState *cellPtr){

	int i, j;
	cellState initState;
	//Initial dendritic parameters
	initState.dend.V_dend = -60;
	initState.dend.Calcium_r = 0.0112788;// High-threshold calcium
	initState.dend.Potassium_s = 0.0049291;// Calcium-dependent potassium
	initState.dend.Hcurrent_q = 0.0337836;// H current
	initState.dend.Ca2Plus = 3.7152;// Calcium concentration
	initState.dend.I_CaH = 0.5;// High-threshold calcium current
	//Initial somatic parameters
	initState.soma.g_CaL = 0.68; //default arbitrary value but it should be randomized per cell
	initState.soma.V_soma = -60;
	initState.soma.Sodium_m = 1.0127807;// Sodium (artificial)
	initState.soma.Sodium_h = 0.3596066;
	initState.soma.Potassium_n = 0.2369847;// Potassium (delayed rectifier)
	initState.soma.Potassium_p = 0.2369847;
	initState.soma.Potassium_x_s = 0.1;// Potassium (voltage-dependent)
	initState.soma.Calcium_k = 0.7423159;// Low-threshold calcium
	initState.soma.Calcium_l = 0.0321349;
	// Initial axonal parameters
	initState.axon.V_axon = -60;
	//sisaza: Sodium_m_a doesn't have a state, therefore this assignment doesn'thave any effect
	initState.axon.Sodium_m_a = 0.003596066;// Sodium (thalamocortical)
	initState.axon.Sodium_h_a = 0.9;
	initState.axon.Potassium_x_a = 0.2369847;// Potassium (transient)

	initState.cell_x=0;
	initState.cell_y=0; //Initial irrelevant value of compartment's ID and coords
	initState.cellID= core_id * cellCount;	//Compute the cellID of the first cell in this core

	//Copy init state to all cell states and calculate their coords
	for(j=0;j<cellCount;j++){
		initState.cell_x= initState.cellID / IO_NETWORK_DIM2;		//we assume that dim1 refers to the number of ROWS and dim2 refers to COLUMNS !!!!
		initState.cell_y= initState.cellID % IO_NETWORK_DIM2;	
		memcpy(&cellPtr[j], &initState, sizeof(cellState));
		initState.cellID ++;				//next cell's ID is increased by 1
	}

	if (G_CAL_FROM_FILE)
		read_g_CaL_from_file(cellPtr);
	
	return;
}

int ReadFileLine(FILE *pInFile, mod_prec *iAppArray){



//	gettimeofday(&tic, NULL);
	char c= fpeek(pInFile);
	if (c==EOF)
		return 0;


	char *strNumber;
	int bufSize = cores*cellCount*20;	//more than enough but will do for nao
	char* buffer = (char*) malloc(bufSize*sizeof(char));

	int floats_ignored = core_id*cellCount;
	int floats_needed = cellCount;
	int useful_element_found = 0;
	//int i will stand for elements already processed, hence it starts at 0 even after the first strtok
	int i = 0;
	
	//Get one line
	if(fgets(buffer, bufSize, pInFile)){
		//Convert the ASCII string of one element to a double precision floating point value
		strNumber = strtok(buffer," ");
		i = 0;

		//printf("Line:\n");
		while ((strNumber != NULL) && (i<IO_NETWORK_SIZE)){
			if (i>=floats_ignored)
				useful_element_found = 1;

			if (i>=floats_ignored+floats_needed)
				useful_element_found = 0;

			if (useful_element_found)
				iAppArray[i-floats_ignored] = atof(strNumber);	//atof() should change if using integers or fixed point
			//printf ("(%s) %0.2f ", strNumber, iAppArray[i]);
			strNumber = strtok(NULL, " ");
			i++; 
		}
		//printf("i: %d\n", i);
		if(i<IO_NETWORK_SIZE){
			//BUG: if only one element is missing but the line ends in a space, the error is not detected
			printf("Error: Input line doesn't have enough elements, only %d\n", i);
			exit(EXIT_FAILURE);
		}
		free(buffer);
		return 1;//success
	}else{
		if(!feof(pInFile)){
			printf("Error: Reading from input file didn't finish successfully\n");
			exit(EXIT_FAILURE);
		}
		free(buffer);
		return 0;//end of file
	}
}


void read_g_CaL_from_file(cellState* cellPtr) {

	mod_prec trash=0;
	int i;

	FILE* fd = fopen("gcal_file.txt","r");
	for (i=0;i<(cellCount*core_id);i++)
		fscanf(fd, "%lf ", &trash);
	for (i=0;i<cellCount;i++) 
		fscanf(fd, "%lf ", &cellPtr[i].soma.g_CaL);
	fclose(fd);

	return;

}

void printState(cellState* cellPtr, char* filename) {

	FILE* fd = fopen(filename, "w");
	mod_prec* s = (mod_prec*) malloc(16*sizeof(mod_prec));
	int i, j;

	for (i=0; i<cellCount; i++) {

		s[0] = cellPtr[i].soma.V_soma;
		s[1] = cellPtr[i].soma.Sodium_m;
		s[2] = cellPtr[i].soma.Potassium_n;
		s[3] = cellPtr[i].soma.Potassium_x_s;
		s[4] = cellPtr[i].soma.Calcium_k;
		s[5] = cellPtr[i].soma.Calcium_l;
		s[6] = cellPtr[i].dend.V_dend;
		s[7] = cellPtr[i].dend.Calcium_r;
		s[8] = cellPtr[i].dend.Potassium_s;
		s[9] = cellPtr[i].dend.Hcurrent_q;
		s[10] = cellPtr[i].dend.Ca2Plus;
		s[11] = cellPtr[i].dend.I_CaH;
		s[12] = cellPtr[i].axon.V_axon;
		s[13] = cellPtr[i].axon.Sodium_m_a;
		s[14] = cellPtr[i].axon.Sodium_h_a;
		s[15] = cellPtr[i].axon.Potassium_x_a;

		for (j=0; j<16; j++)
			fprintf(fd, "%.8lf ", s[j]);
		fprintf(fd, "\n");

	}

	fclose(fd);
	return;

}

inline mod_prec min(mod_prec a, mod_prec b){

	return (a < b) ? a : b;
}


/* Creates a double linked list for each core according to the in and out connections of 
 * each core cell. List is sorted by communicating cells involved with core cells.
 * Arguments are the connectivity filename, id of the core that will run and cellCompParams 
 * array in order to store the number of neighbors for each core cell (core_cells[i]->total_amount_of_neighbors)
 * ---LEGACY FUNCTION, NO LONGER IN USE----
 */


communication_node* Make_Core_Communication_List(char *filename, int core_id,cellCompParams *core_cells) {

	FILE *input_file = fopen(filename,"r");
	communication_node *temporary = NULL,*list_head = NULL;
	int line_counter = 0, starting_line = (core_id*cellCount),ending_line = (core_id*cellCount) + cellCount;
	int in_cell_id, out_cell_id,count_neighbors;
	char sample_character;

	/* skip cell conectivity lines 
	 * with cells simulated by other cores
	 */
	for(;line_counter<starting_line;line_counter++)
		while (fgetc(input_file)!='\n');
	
	/* read only cell lines for this specific core 
	 * and make the communication list while parsing the file
	 * we read first the format i cell_id cell_id ... that correspond to incoming connections of this cell
	 * then o cell_id cell_id that corresponds tou outgoing connections of this cell
	 * Decomment if neccessary lines with RCCE api,for SCC.
	 */
	while (line_counter<ending_line) {
		count_neighbors = 0;
		fscanf(input_file,"%c ",&sample_character);
		if (sample_character!='i') {
			printf("Error : Wrong file format of file describing connectivity. Leaving the core..");

			MPI_Finalize();
			exit(1);
		}
		
		sample_character = fpeek(input_file);
		while (sample_character != 'o') {
			fscanf(input_file,"%d ",&in_cell_id);
			temporary = (communication_node *) malloc(sizeof(communication_node));
			temporary->core_cell = line_counter;
			temporary->other_cell = in_cell_id;
			temporary->connection = 'r';
			list_head = Insert_Into_Communication_List(list_head,temporary);
			sample_character = fpeek(input_file);
			count_neighbors++;		
		}
		
		core_cells[line_counter%cellCount].total_amount_of_neighbours = count_neighbors;
		
		fscanf(input_file,"%c ",&sample_character);
		if (sample_character != 'o') {
			printf("Error : Wrong file format of file describing connectivity. Leaving the core..");

			MPI_Finalize();
			exit(1);
		}
		
		
		while (sample_character != '|') {
			fscanf(input_file,"%d ",&out_cell_id);
			temporary = (communication_node *) malloc(sizeof(communication_node));
			temporary->core_cell = line_counter;
			temporary->other_cell = out_cell_id;
			temporary->connection = 's';
			list_head = Insert_Into_Communication_List(list_head,temporary);
			sample_character = fpeek(input_file);
		}
 

		fscanf(input_file,"%c",&sample_character);
		if (sample_character != '|') {
			printf("Error : Wrong file format of file describing connectivity. Leaving the core..");

			MPI_Finalize();
			exit(1);
		}

		fscanf(input_file,"%c",&sample_character);
		if (sample_character==EOF)
			return temporary;

		line_counter++;
	}

	/* return the head of the communication 
	 * list to the caller
	 */
	return list_head;

}


/* Function that inserts a node into the list 
 * sorted by node->other_cell, in case where two nodes have same 
 * other cell then a secondary sort is performed by node->core_cell
 * --LEGACY FUNCTION, NO LONGER IN USE--
 */

communication_node* Insert_Into_Communication_List (communication_node *head,communication_node *node) {

	communication_node *next, *previous;

	/* head insertion special case where 
	 * the list is empty 
	 */
	if (!head) {
		node->next = NULL;
		node->previous = NULL;
		return node;
	}

	/* insert the node to the right place
	 */
	
	next = head;
	while ( (node->other_cell > next->other_cell) || ( (node->other_cell == next->other_cell) && (node->core_cell > next->core_cell)) ) {
	
		previous = next;
		next = next->next;
		/* case where node must be
		 * inserted at the end of the list
		 */

		if (!next) {
			previous->next = node;
			node->previous = previous;
			node->next = NULL;
			return head;
		}
	}

	/* if node already exists change
	 * communicaton_type to SEND/RECEIVE 
	 * because node already added at inocoming
	 * connections
	 */
	if ( (node->core_cell == next->core_cell) && (node->other_cell == next->other_cell) ){
	
		next->connection = 'b';
		free(node);
		return head;
	}
	

	/* case where node must be 
	 * the new head of the list 
	 */
	if (next == head){
		node->next = head;
		head->previous = node;
		node->previous = NULL;
		return node;
	} 

	/* simple insertion case */
	node->next = next;
	node->previous = next->previous;
	(next->previous)->next = node;
	next->previous = node;
	return head;
}




communication_node* Make_Core_Communication_List_new_format(char *filename, int core_id, cellCompParams *core_cells) {

	FILE *input_file = fopen(filename,"r");
	communication_node *temporary = NULL,*list_head = NULL;
	int line_counter = 0, starting_line = (core_id*cellCount),ending_line = (core_id*cellCount) + cellCount, number_of_conductances;
	int communicating_cell_id,count_neighbors;
	char sample_character, sample_character_2;
	mod_prec conductance_read;

	/* skip cell conectivity lines 
	 * with cells simulated by other cores
	 */
	
	for(;line_counter<starting_line;line_counter++)
		while (fgetc(input_file)!='\n');
	
	/* read only cell lines for this specific core 
	 * and make the communication list while parsing the file
	 * we read first the format i cell_id cell_id ... that correspond to incoming connections of this cell
	 * then o cell_id cell_id that corresponds tou outgoing connections of this cell
	 * Decomment if neccessary lines with RCCE api,for SCC.
	 */
	while (line_counter<ending_line) {
		
		count_neighbors = 0;
		number_of_conductances = 0;
		sample_character = fpeek(input_file);
		
		while (sample_character != '|') {
			fscanf(input_file,"%d %c ",&communicating_cell_id/* ,&conductance_read*/, &sample_character);
			//core_cells[line_counter%cellCount].neighConductances[number_of_conductances] = conductance_read;
			//number_of_conductances++;
			
			temporary = (communication_node *) malloc(sizeof(communication_node));
			temporary->core_cell = line_counter;
			temporary->other_cell = communicating_cell_id;

		//	printf("Core %d List : %d -> %d ", core_id, temporary->core_cell, temporary->other_cell);

			switch(sample_character) {
				case ('B'):
					temporary->connection = 'b';
					count_neighbors++;
					break;
				case ('S'):
					temporary->connection = 's';
					break;
				case ('R'):
					temporary->connection = 'r';
					count_neighbors++;
					break;
			}

			list_head = Insert_Into_Communication_List_format_2 (list_head,temporary);	
			sample_character = fpeek(input_file);	
		}
		core_cells[line_counter%cellCount].total_amount_of_neighbours = count_neighbors;
		
		fscanf(input_file,"|%c",&sample_character);

		if (sample_character==EOF)
			return list_head;

		line_counter++;
	}

	/* return the head of the communication 
	 * list to the caller
	 */
	return list_head;

}

sending_node* Make_Core_Communication_List_newest_format(char *filename, cellCompParams *core_cells) {

/*
 * Fair Warning: This function needs heavy testing and is currently being designed as a way to use conductivity matrixes, which seem to be mandatory now, instead of the old
 * communication format. A conductivity matrix basically describes that cell x in rows sends to cell y in columns via a conductivity value z = conductivity_matrix[x][y].
 * The matrix does not have to be symmetrical. Actually the only limit assumed here is that conductivity_matrix[i][i]==0 (no self-feeding, that would make no sense).
 * This function builds a sending list and the necessary receiving buffers of the cellCompParams structure. The new sending list should behave much better since now
 * all cells that need to be sent from one core to the other are BUNDLED TOGETHER, so that we can send all of them in one go (and not sending
 * duplicates of the same cell to the same core either). Read below for functions completing communication based on the structures designed here.
*/
	FILE *input_file = fopen(filename,"r");
	sending_node *temporary = NULL, *list_head = NULL;
	int i, k, my_cell, my_sending_cell, core_to_send, line_counter, my_cell_lower_bound = cellCount*core_id, my_cell_upper_bound = my_cell_lower_bound + cellCount;
	int* my_neighbour_count = (int*) calloc(cellCount, sizeof(int));
	mod_prec cond_value;

	//this list has one node for every core employed, detailing which cells must be sent to each core. The node that corresponds to this core will be kept blank
	list_head = (sending_node*) malloc(1*sizeof(sending_node));
	list_head->target_core = 0;
	list_head->next = NULL;
	list_head->total_cells_to_send = 0;
	list_head->cells_to_send = NULL;
	list_head->voltages_to_send = NULL;
	temporary = list_head;

	for (i=1;i<cores;i++) {
		temporary->next = (sending_node*) malloc(1*sizeof(sending_node));
		temporary = temporary->next;
		temporary->target_core = i;
		temporary->next = NULL;
            	temporary->total_cells_to_send = 0;
               	temporary->cells_to_send = NULL;
		temporary->voltages_to_send = NULL;
	}

	//it is important to keep in mind that the matrix is read as if cell in ROW is sending information to cell in COLUMN!

	for (line_counter=0;line_counter<IO_NETWORK_SIZE;line_counter++) {

		temporary = list_head;				//before examining every (supposedly sending) line, we reset temporary to the beginning of our sending list
		for (i=0; i<IO_NETWORK_SIZE; i++) {

			fscanf(input_file, "%f ", &cond_value);
			if (cond_value==0)						//this connection is considered not existing if conductance = 0
				;
			else {

			//part of the code handling RECEIVING and noting which of my cells needs input from which other cells, from ANY core

				if ((i>=my_cell_lower_bound)&&(i<my_cell_upper_bound)) {				//these are this core's cells' COLUMNS (incoming)
					my_cell = i%cellCount;
					if (my_neighbour_count[my_cell]==0) {                  //if this is the first neighbour, initialize buffers
						core_cells[my_cell].neighVdend = NULL;
						core_cells[my_cell].neighConductances = NULL;
						core_cells[my_cell].neighId = NULL;
					}

					core_cells[my_cell].neighId = allocate_space_int(core_cells[my_cell].neighId, my_neighbour_count[my_cell]);
					core_cells[my_cell].neighId[my_neighbour_count[my_cell]] = line_counter;		//which cell sends this voltage to us (GLOBAL ID)

					core_cells[my_cell].neighConductances = allocate_space_mod(core_cells[my_cell].neighConductances, my_neighbour_count[my_cell]);
					core_cells[my_cell].neighConductances[my_neighbour_count[my_cell]] = cond_value;	//what conductance we use to calculate its impact

					//allocate space for storing this voltage
					core_cells[my_cell].neighVdend = allocate_space_mod(core_cells[my_cell].neighVdend, my_neighbour_count[my_cell]);

					my_neighbour_count[my_cell]++;						//how many neighbours this cell has so far (from ANY core)
					core_cells[my_cell].total_amount_of_neighbours = my_neighbour_count[my_cell];
				}

			//part of the code handling SENDING and bundling together my cells I need to send each core

				if ((line_counter>=my_cell_lower_bound)&&(line_counter<my_cell_upper_bound)) {		//these are this core's cells' ROWS (outgoing)
					core_to_send = i/cellCount;							//which core I need to send to
					while ((temporary!=NULL)&&(temporary->target_core!=core_to_send))		//search for the right list node
						temporary = temporary->next;

					if (core_to_send==core_id)							//Obviously I do not send to my own core
						;
					else if (temporary!=NULL) {
						k = temporary->total_cells_to_send;
						my_sending_cell = line_counter%cellCount;

						if ((k>0)&&(temporary->cells_to_send[k-1]==my_sending_cell))		//I have already marked this cell to send core_to_send
							;
						else {
							//add this cell to the core_to_send's cells_to_send-array ->WITH LOCAL ID<-
							temporary->cells_to_send = allocate_space_int(temporary->cells_to_send, temporary->total_cells_to_send);
							temporary->cells_to_send[k] = my_sending_cell;
							temporary->voltages_to_send = allocate_space_mod(temporary->voltages_to_send, temporary->total_cells_to_send);
							temporary->total_cells_to_send++;
						}
					}
				}

			//end of the code handling proper list and buffer creation concerning communication

			}
		}

	}

//	for (i=0; i<cellCount; i++)
//		printf("%d) [%d] = (%d-%lf,%d-%lf,%d-%lf,%d-%lf)\n", core_id, i, core_cells[i].neighId[0], core_cells[i].neighConductances[0], core_cells[i].neighId[1], core_cells[i].neighConductances[1], core_cells[i].neighId[2], core_cells[i].neighConductances[2], core_cells[i].neighId[3], core_cells[i].neighConductances[3]);

	fclose(input_file);
	return list_head;

}

receiving_node* reckon_phase(sending_node *sending_list_head, cellCompParams *core_cells) {

/* In this function, we will make the opposite of the sending_list: We will create a list detailing what
 * information we will receive from each other core. Since cell mapping and connections do not change,
 * we will create this list once in the beginning so that during simulation, only voltages are to be exchanged
 * This list will be created via some information exchange between cores (although it could be built during
 * the conductivity matrix parsing, this seems simpler)
 */
	receiving_node *r_temporary = NULL, *receiving_list_head = NULL; 
	sending_node *s_temporary = NULL;
        int i, k, targetCore;
	MPI_Request s_request;	//irrelevant
	MPI_Request* request1 = (MPI_Request*) malloc(cores*sizeof(MPI_Request));
	MPI_Request* request2 = (MPI_Request*) malloc(cores*sizeof(MPI_Request));

	receiving_list_head = (receiving_node*) malloc(1*sizeof(receiving_node));
	receiving_list_head->target_core = 0;
	receiving_list_head->next = NULL;
	receiving_list_head->total_cells_to_receive = 0;
	receiving_list_head->cells_to_receive = NULL;
	receiving_list_head->voltages_to_receive = NULL;
	r_temporary = receiving_list_head;

	for (i=1;i<cores;i++) {
		r_temporary->next = (receiving_node*) malloc(1*sizeof(receiving_node));
		r_temporary = r_temporary->next;
		r_temporary->target_core = i;
		r_temporary->next = NULL;
		r_temporary->total_cells_to_receive = 0;
		r_temporary->cells_to_receive = NULL;
		r_temporary->voltages_to_receive = NULL;
	}

	s_temporary = sending_list_head;
	r_temporary = receiving_list_head;

	//phase 1: inform cores how many cells need to be exchanged

	k = 0;
	
	MPI_Status status;
	
	while ((r_temporary!=NULL)&&(s_temporary!=NULL)) {
		targetCore = r_temporary->target_core;
		if (targetCore == core_id)
			request1[k] = MPI_REQUEST_NULL;
		else {
			//MPI_Isend(&(s_temporary->total_cells_to_send), 1, MPI_INT, targetCore, 0, MPI_COMM_WORLD, &s_request);
			MPI_Send(&(s_temporary->total_cells_to_send), 1, MPI_INT, targetCore, 0, MPI_COMM_WORLD);
			//MPI_Irecv(&(r_temporary->total_cells_to_receive), 1, MPI_INT, targetCore, 0, MPI_COMM_WORLD, &request1[k]);
			MPI_Recv(&(r_temporary->total_cells_to_receive), 1, MPI_INT, targetCore, 0, MPI_COMM_WORLD,&status);
		}
		r_temporary = r_temporary->next;
		s_temporary = s_temporary->next;
		k++;
	}

	//syncing(request1);		//phase 1 complete, waiting for sync

	// phase 2: allocate buffers to hold the incoming cells and exchange cell ids to know who sends what
	// WARNING: the receiving buffer here gets filled with ids as the sending core sends them - they are NOT global ids
	// and thus will NEED "translation"

	s_temporary = sending_list_head;
        r_temporary = receiving_list_head;

	k = 0;
	while ((r_temporary!=NULL)&&(s_temporary!=NULL)) {
                targetCore = r_temporary->target_core;
                if (targetCore == core_id)
                        request2[k] = MPI_REQUEST_NULL;
                else {

			if (s_temporary->total_cells_to_send != 0)
				//MPI_Isend(&(s_temporary->cells_to_send[0]), s_temporary->total_cells_to_send, MPI_INT, targetCore, 0, MPI_COMM_WORLD, &s_request);
				MPI_Send(&(s_temporary->cells_to_send[0]), s_temporary->total_cells_to_send, MPI_INT, targetCore, 0, MPI_COMM_WORLD);
			if (r_temporary->total_cells_to_receive != 0) {
				r_temporary->cells_to_receive = (int*) malloc(r_temporary->total_cells_to_receive * sizeof(int));
				r_temporary->voltages_to_receive = (mod_prec*) malloc(r_temporary->total_cells_to_receive * sizeof(mod_prec));
				//MPI_Irecv(&(r_temporary->cells_to_receive[0]), r_temporary->total_cells_to_receive, MPI_INT, targetCore, 0, MPI_COMM_WORLD, &request2[k]);
				MPI_Recv(&(r_temporary->cells_to_receive[0]), r_temporary->total_cells_to_receive, MPI_INT, targetCore, 0, MPI_COMM_WORLD, &status);
			} else {
				request2[k] = MPI_REQUEST_NULL;
			}
		}
                r_temporary = r_temporary->next;
                s_temporary = s_temporary->next;
		k++;
        }

	//syncing(request2);		//phase 2 complete, waiting for sync

	return receiving_list_head;

}

void perform_communication_step(sending_node* sending_list_head, receiving_node* receiving_list_head, cellCompParams* params, cellState* cells) {

/* The function where the magic happens: using the structures created during initialization, exchange necessary dendritic voltage
 * between all cores. This time, we use non-blocking functions and voltages are bundled together per core-target (and this time,
 * no duplicates when a core needs a voltage for more than one of its own cells).
 */

	sending_node* s_temp;
	receiving_node* r_temp;
	int* processing_pointer = (int*) calloc(cellCount, sizeof(int));
	int i, j, k=0, cell_id, translated_cell_id, requested_neighbour, my_requested_cell;
	MPI_Request s_request;		//irrelevant
	MPI_Request* request = (MPI_Request*) malloc(cores*sizeof(MPI_Request));	//array of request handles for MPI_Irecv
		
	MPI_Status status;	
	//phase 1: we send all of our necessary voltages
	
	s_temp = sending_list_head;
	while (s_temp!=NULL) {
		if (s_temp->total_cells_to_send==0)
			;
		else {
			for (i=0; i<s_temp->total_cells_to_send; i++)			//fill the voltage buffer we shall send with the cells we need to send
				s_temp->voltages_to_send[i] = cells[s_temp->cells_to_send[i]].dend.V_dend;
			MPI_Isend(&(s_temp->voltages_to_send[0]), s_temp->total_cells_to_send, MPI_FLOAT, s_temp->target_core, 0, MPI_COMM_WORLD, &s_request);
		}
		s_temp = s_temp->next;
	}
	
	//phase 2: we request all of the voltages we need and we sync the cores together
	
	r_temp = receiving_list_head;
	while (r_temp!=NULL) {
		if (r_temp->total_cells_to_receive==0)
			request[k] = MPI_REQUEST_NULL;
		else
			//MPI_Irecv(&(r_temp->voltages_to_receive[0]), r_temp->total_cells_to_receive, MPI_FLOAT, r_temp->target_core, 0, MPI_COMM_WORLD, &request[k]);
			MPI_Recv(&(r_temp->voltages_to_receive[0]), r_temp->total_cells_to_receive, MPI_FLOAT, r_temp->target_core, 0, MPI_COMM_WORLD, &status);
		k++;
		r_temp = r_temp->next;
	}

	//syncing(request);

	//phase 3: we will now process what we received from all cores in one passing of our structs altogether

	r_temp = receiving_list_head;
	while (r_temp!=NULL) {
		if (r_temp->total_cells_to_receive == 0)
			;
		else
			for (i=0; i<r_temp->total_cells_to_receive; i++) {		//we will run through the list we received
				cell_id = r_temp->cells_to_receive[i];
				translated_cell_id = cell_id + cellCount*(r_temp->target_core);			//global id of the cell under examination

				for (j=0; j<cellCount; j++) {							//we shall examine each one of our cells to see which one needs it
					if (processing_pointer[j]>=params[j].total_amount_of_neighbours)
						;
					else {
						requested_neighbour = params[j].neighId[processing_pointer[j]];		//id of the next neighbour this cell needs
						while (((requested_neighbour/cellCount)==core_id)&&(requested_neighbour!=-1)) {			//if the requested neighbour belongs to this core

							my_requested_cell = requested_neighbour%cellCount;		//which of my cells is needed
							params[j].neighVdend[processing_pointer[j]] = cells[my_requested_cell].dend.V_dend;	//the voltage of the cell requested

							processing_pointer[j]++;
							if (processing_pointer[j]<params[j].total_amount_of_neighbours)
								requested_neighbour = params[j].neighId[processing_pointer[j]];		//request next neighbour
							else
								requested_neighbour = -1;						//if we are done then make sure we do not request any more cells
						}
						if (requested_neighbour==translated_cell_id) {						//we found a cell we need
							params[j].neighVdend[processing_pointer[j]] = r_temp->voltages_to_receive[i];
							processing_pointer[j]++;
						}
					}
				}
			}

		r_temp = r_temp->next;
	}

	//phase 4: check whether all cells have filled up their neighVdend buffers (some of them may still need info from our own cells)
	
	for (j=0; j<cellCount; j++) {
		if (processing_pointer[j]>=params[j].total_amount_of_neighbours)
			;
		else {
			requested_neighbour = params[j].neighId[processing_pointer[j]];
			while (((requested_neighbour/cellCount)==core_id)&&(requested_neighbour!=-1)) {

				my_requested_cell = requested_neighbour%cellCount;
				params[j].neighVdend[processing_pointer[j]] = cells[my_requested_cell].dend.V_dend;

				processing_pointer[j]++;
				if (processing_pointer[j]<params[j].total_amount_of_neighbours)
					requested_neighbour = params[j].neighId[processing_pointer[j]];
				else
					requested_neighbour = -1;
			}
		}
	}
	
	//all tasks have been completed and processed so we can move out of the communication function
	return;

}

int* allocate_space_int(int* pointer, int existing_slots) {

	int new_total_slots = existing_slots + 1;
	int* new_pointer = (int*) realloc(pointer, new_total_slots*sizeof(int));
	return new_pointer;
}

mod_prec* allocate_space_mod(mod_prec* pointer, int existing_slots) {

	int new_total_slots = existing_slots + 1;
        mod_prec* new_pointer = (mod_prec*) realloc(pointer, new_total_slots*sizeof(mod_prec));
	return new_pointer;
}

void syncing(MPI_Request* request) {

/* a function that first makes sure all requests passed on to it are completed
 * and then calls an MPI_Barrier to make sure all processes are in sync. We expect
 * one request from each core employed by the app
 */

	int i, done = 0, flag, no_of_reqs = cores;

	while (!done) {
		done = 1;
		for (i=0;i<no_of_reqs;i++) {
			if (request[i] == MPI_REQUEST_NULL)
				;
			else {
				MPI_Test(&(request[i]), &flag, MPI_STATUS_IGNORE);
				if (!flag)
					done = 0;
			}
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
	return;

}

communication_node* Insert_Into_Communication_List_format_2 (communication_node *head,communication_node *node) {

	communication_node *next, *previous;

	/* head insertion special case where 
	 * the list is empty 
	 */
	if (!head) {
		node->next = NULL;
		node->previous = NULL;
		return node;
	}

	/* insert the node to the right place
	 */
	
	next = head;
	while ( (node->other_cell > next->other_cell) || ( (node->other_cell == next->other_cell) && (node->core_cell > next->core_cell)) ) {
	
		previous = next;
		next = next->next;
		/* case where node must be
		 * inserted at the end of the list
		 */

		if (!next) {
			previous->next = node;
			node->previous = previous;
			node->next = NULL;
			return head;
		}
	}

	/* case where node must be 
	 * the new head of the list 
	 */
	if (next == head){
		node->next = head;
		head->previous = node;
		node->previous = NULL;
		return node;
	} 

	/* simple insertion case */
	node->next = next;
	node->previous = next->previous;
	(next->previous)->next = node;
	next->previous = node;
	return head;
}


/* Check if it works */

void Perform_Communication_Version_2(communication_node *communication_list_head, cellCompParams *cellParamsPtr, cellState **cellPtr, int simulation_array_ID) {

	
	
	int targetCore,index,core_cell_id,other_cell_id,my_storing_cell,other_storing_cell;
	communication_node *communicating_cell;
	mod_prec received_voltage, sent_voltage;
	
	communicating_cell = communication_list_head;


   	while (communicating_cell) {	
   		targetCore = (communicating_cell->other_cell)/cellCount;	

   		/* communication with cell simulated 
   		 * by another core , if core_id != targetCore
   		 * else the cell is simulated by the same cores
   		 */
   		 MPI_Status status;
   		if (core_id != targetCore) {	
			printf("OK3\n");
			core_cell_id = (communicating_cell->core_cell)%cellCount;
			switch(communicating_cell->connection) {
				

				case('s'):
					sent_voltage = cellPtr[simulation_array_ID][core_cell_id].dend.V_dend;
					MPI_Send(&sent_voltage, 1, MPI_DOUBLE, targetCore, 0, MPI_COMM_WORLD);
					break;

				case('r'):
					MPI_Recv(&received_voltage, 1, MPI_DOUBLE, targetCore, 0, MPI_COMM_WORLD, &status);
					index = cellParamsPtr[core_cell_id].index_of_neighVdend;
					cellParamsPtr[core_cell_id].neighVdend[index] = received_voltage;
					cellParamsPtr[core_cell_id].index_of_neighVdend++;
					if (cellParamsPtr[core_cell_id].index_of_neighVdend == cellParamsPtr[core_cell_id].total_amount_of_neighbours)
						cellParamsPtr[core_cell_id].index_of_neighVdend = 0; 
					break ;

				case('b'):
					if (core_id < targetCore) {
						sent_voltage = cellPtr[simulation_array_ID][core_cell_id].dend.V_dend;
						other_storing_cell = (communicating_cell->other_cell)%cellCount;


						MPI_Recv(&received_voltage, 1, MPI_DOUBLE, targetCore, 0, MPI_COMM_WORLD, &status);
						MPI_Send(&sent_voltage, 1, MPI_DOUBLE, targetCore, 0, MPI_COMM_WORLD);

						index = cellParamsPtr[my_storing_cell].index_of_neighVdend;
						cellParamsPtr[my_storing_cell].neighVdend[index] = received_voltage;
						cellParamsPtr[my_storing_cell].index_of_neighVdend++;
						if (cellParamsPtr[my_storing_cell].index_of_neighVdend == cellParamsPtr[my_storing_cell].total_amount_of_neighbours)
							cellParamsPtr[my_storing_cell].index_of_neighVdend = 0;

						break ;
					}
					else {
						sent_voltage = cellPtr[simulation_array_ID][core_cell_id].dend.V_dend;
						other_storing_cell = (communicating_cell->other_cell)%cellCount;
						MPI_Send(&sent_voltage, 1, MPI_DOUBLE, targetCore, 0, MPI_COMM_WORLD);
						MPI_Recv(&received_voltage, 1, MPI_DOUBLE, targetCore, 0, MPI_COMM_WORLD, &status);

						index = cellParamsPtr[my_storing_cell].index_of_neighVdend;
						cellParamsPtr[my_storing_cell].neighVdend[index] = received_voltage;
						cellParamsPtr[my_storing_cell].index_of_neighVdend++;
						if (cellParamsPtr[my_storing_cell].index_of_neighVdend == cellParamsPtr[my_storing_cell].total_amount_of_neighbours)
							cellParamsPtr[my_storing_cell].index_of_neighVdend = 0;
						
						break;

					}
							
			}
		}	
		else {			
			core_cell_id  = (communicating_cell->core_cell)%cellCount;
			other_cell_id = (communicating_cell->other_cell)%cellCount;

			/* no core message passing 
			 * change only values of cellParamsPtr array
			 */
			switch(communicating_cell->connection) {
				case('s'):
					index = cellParamsPtr[other_cell_id].index_of_neighVdend;
					cellParamsPtr[other_cell_id].neighVdend[index] = cellPtr[simulation_array_ID][core_cell_id].dend.V_dend;
					cellParamsPtr[other_cell_id].index_of_neighVdend++;
					if (cellParamsPtr[other_cell_id].index_of_neighVdend == cellParamsPtr[other_cell_id].total_amount_of_neighbours)
						cellParamsPtr[other_cell_id].index_of_neighVdend = 0;
					break;

				case('r'):
					index = cellParamsPtr[core_cell_id].index_of_neighVdend;
					cellParamsPtr[core_cell_id].neighVdend[index] = cellPtr[simulation_array_ID][other_cell_id].dend.V_dend;
					cellParamsPtr[core_cell_id].index_of_neighVdend++;
					if (cellParamsPtr[core_cell_id].index_of_neighVdend == cellParamsPtr[core_cell_id].total_amount_of_neighbours)
						cellParamsPtr[core_cell_id].index_of_neighVdend = 0;
					break;

				case('b'):
					index = cellParamsPtr[other_cell_id].index_of_neighVdend;
					cellParamsPtr[other_cell_id].neighVdend[index] = cellPtr[simulation_array_ID][core_cell_id].dend.V_dend;
					cellParamsPtr[other_cell_id].index_of_neighVdend++;

					index = cellParamsPtr[core_cell_id].index_of_neighVdend;
					cellParamsPtr[core_cell_id].neighVdend[index] = cellPtr[simulation_array_ID][other_cell_id].dend.V_dend;


					cellParamsPtr[core_cell_id].index_of_neighVdend++;

					if (cellParamsPtr[core_cell_id].index_of_neighVdend == cellParamsPtr[core_cell_id].total_amount_of_neighbours)
						cellParamsPtr[core_cell_id].index_of_neighVdend = 0;
					if (cellParamsPtr[other_cell_id].index_of_neighVdend == cellParamsPtr[other_cell_id].total_amount_of_neighbours)
						cellParamsPtr[other_cell_id].index_of_neighVdend = 0;
		
					break;
						
			}
		}
		
		communicating_cell = communicating_cell->next;
	}



	return;


}

void print_cellState(cellState *a) {
    printf("cellID: %d\n", a->cellID);
    printf("cell_x: %d\n", a->cell_x);
    printf("cell_y: %d\n", a->cell_y);
    printf("dend.V_dend %lf\n", a->dend.V_dend);
    printf("dend.Hcurrent_q %lf\n", a->dend.Hcurrent_q);
    printf("dend.Calcium_r: %lf\n", a->dend.Calcium_r);
    printf("dend.Potassium_s: %lf\n", a->dend.Potassium_s);
    printf("dend.I_CaH: %lf\n", a->dend.I_CaH);
    printf("dend.Ca2Plus: %lf\n", a->dend.Ca2Plus);
    printf("soma.g_CaL: %lf\n", a->soma.g_CaL);
    printf("soma.V_soma: %lf\n", a->soma.V_soma);
    printf("soma.Sodium_m: %lf\n", a->soma.Sodium_m);
    printf("soma.Sodium_h: %lf\n", a->soma.Sodium_h);
    printf("soma.Calcium_k: %lf\n", a->soma.Calcium_k);
    printf("soma.Calcium_l: %lf\n", a->soma.Calcium_l);
    printf("soma.Potassium_n: %lf\n", a->soma.Potassium_n);
    printf("soma.Potassium_p: %lf\n", a->soma.Potassium_p);
    printf("soma.Potassium_x_s: %lf\n", a->soma.Potassium_x_s);
    printf("axon.V_axon: %lf\n", a->axon.V_axon);
    printf("axon.Sodium_m_a: %lf\n", a->axon.Sodium_m_a);
    printf("axon.Sodium_h_a: %lf\n", a->axon.Sodium_h_a);
    printf("axon.Potassium_x_a: %lf\n", a->axon.Potassium_x_a);
}

void print_cellCompParams(cellCompParams *a) {
    printf("index_of_neighVdend: %d\n", a->index_of_neighVdend);
    printf("total_amount_of_neighbors: %d\n", a->total_amount_of_neighbours);
    printf("iAppIn: %lf\n", a->iAppIn);
    printf("neighVdend[0]: %lf\n", a->neighVdend[0]);
}


void take_checkpoint(int simStep, cellState **cellPtr, cellCompParams * cellParamsPtr, \
                    int cellCount, FILE *ckptfd, void * ckptBuffer) {
    int i;
    void * location;

    if (core_id==0){
	FILE *fp;
	char ckpTime[100];

	sprintf(ckpTime,"/home/apostolis/Desktop/mpi_programs/ckpt_infoli/ckpTime/%d.txt",simStep);
	fp=fopen(ckpTime,"w");
	fprintf(fp,"%d\n",simStep);
	if (fsync(fileno(fp)) !=0){
		printf("Error syncing ckpTime file...\n");
	}
	fclose(fp);
    }


    memcpy(ckptBuffer, &simStep, sizeof(int));
    for (i=0; i<cellCount; i++) {
        location = ckptBuffer + sizeof(int) + i * (sizeof(cellState) + sizeof(int));
        memcpy(location, &(cellPtr[simStep % 2][i]), sizeof(cellState));
        memcpy(location + sizeof(cellState), &(cellParamsPtr[i].index_of_neighVdend), sizeof(int));
    }

    if ((simStep / CKPT_INTERVAL) % 2 == 0) {
        fseek(ckptfd, 2*sizeof(int), SEEK_SET);
    } else {
        fseek(ckptfd, 3*sizeof(int) + (cellCount * (sizeof(cellState) + sizeof(int))), SEEK_SET);
    }
    if (fwrite(ckptBuffer, sizeof(int) + (cellCount * (sizeof(cellState) + sizeof(int))), 1, ckptfd) != 1)
        fprintf(stderr, "Error: Checkpoint write failed");
    fflush(ckptfd);
    if (fsync(fileno(ckptfd)) != 0)
        printf("Error: Checkpoint file sync failed\n");

}

/*
void take_checkpoint(int simStep, cellState **cellPtr, cellCompParams * cellParamsPtr, \
                    int cellCount, FILE *ckptfd, void * ckptBuffer) {
    int i;
    void * location;
    memcpy(ckptBuffer, &simStep, sizeof(int));
    for (i=0; i<cellCount; i++) {
        location = ckptBuffer + sizeof(int) + i * (sizeof(cellState) + sizeof(int));
        memcpy(location, &(cellPtr[simStep % 2][i]), sizeof(cellState));
        memcpy(location + sizeof(cellState), &(cellParamsPtr[i].index_of_neighVdend), sizeof(int));
    }

    if ((simStep / CKPT_INTERVAL) % 2 == 0) {
        fseek(ckptfd, 2*sizeof(int), SEEK_SET);
    } else {
        fseek(ckptfd, 3*sizeof(int) + (cellCount * (sizeof(cellState) + sizeof(int))), SEEK_SET);
    }
    if (fwrite(ckptBuffer, sizeof(int) + (cellCount * (sizeof(cellState) + sizeof(int))), 1, ckptfd) != 1)
        fprintf(stderr, "Error: Checkpoint write failed");
    fflush(ckptfd);
    if (fsync(fileno(ckptfd)) != 0)
        printf("Error: Checkpoint file sync failed\n");
}
*/



