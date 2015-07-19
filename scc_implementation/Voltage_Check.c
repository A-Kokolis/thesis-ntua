#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/time.h>
#include "Voltage_Check.h"


void
init_voltage_files (void)
{
  int res, loop_count = 0;
  int v[7] = { 0, 1, 2, 3, 4, 5, 7 };
  char temp[200];
  while (loop_count < 7)
    {
      sprintf (temp, "sccBmc -c status | grep VCC%d > voltage_level_%d.txt",
	       v[loop_count], v[loop_count]);
      res = system (temp);
      loop_count++;
    }
}


double
voltage_check (double difference)
{
  int token_count = 0;
  int changed[6]; //keep track of the voltage islands that have already change

  int v[6] = { 0, 1, 3, 4, 5, 7 }; //these are the voltage islands of SCC that i test
  int loop_count = 0;
  char temp[1024];

  double time = 0, voltage[6], current[6], prev_voltage[6], prev_current[6];
  char line[1024];
  char *tok = NULL;
  struct timeval tts, ttf,ttemp;   //inorder to check the overhead of voltage change and the timesteps between each sccBmc status command
  FILE *file;
  gettimeofday (&tts, NULL);
  
  //read prev_values from files
  
  while (loop_count<6){
  	changed[loop_count]=0;
	sprintf (temp, "voltage_level_%d.txt", v[loop_count]);
	file = fopen (temp, "r");
	if (file != NULL)
	{
	  if (fgets (line, sizeof line, file) != NULL)
	    {
	      token_count = 0;
	      tok = strtok (line, " ");
	      while (token_count < 4)
		{
//		  printf ("token= %s\n", tok);
		  if (token_count == 2)
		    {		//it was ==1
		      prev_voltage[loop_count] = atof (tok);
		    }
		  //    if (token_count==3){
		  //           prev_current=atof(tok);
		  //    }
		  token_count++;
		  tok = strtok (NULL, " ");
		}
	    }
	    else
	    {
	      printf ("could not locate file\n");
	      return -1;
	    }
	}
      loop_count++;	
    //  printf ("prev_voltage[%d]=%lf loop_count=%d\n",loop_count-1,prev_voltage[loop_count-1],loop_count);
      fclose (file);
 }
	
int	external_loop_count=0; // set external_loop_count to 0, to check for voltage changes  



  while (external_loop_count < 6)
    {

	  //printf("Entered !change loop_count=%d\n",loop_count);
	  sprintf (temp,"./execute_sccBmc.sh");  //IMPORTANT "execute_sccBmc.sh script needed for the execution of commands
	  if (!system (temp))
	    {
//	    	  gettimeofday (&ttemp, NULL);
//	          sprintf(temp,"echo command executed timestep= %lf>>voltage_level_all.txt",(ttemp.tv_sec - tts.tv_sec) + (ttemp.tv_usec - tts.tv_usec) * 0.000001);
//	          system(temp);
	    	  loop_count=0;
		  while (loop_count<6){
		  	
		     if (!changed[loop_count]){
			sprintf (temp, "voltage_level_%d.txt", v[loop_count]);
			file = fopen (temp, "r");
			if (file != NULL)
			{
			  if (fgets (line, sizeof line, file) != NULL)
			    {
			      token_count = 0;
			      tok = strtok (line, " ");
			      while (token_count < 4)
				{
		//		  printf ("token= %s\n", tok);
				  if (token_count == 2)
				    {		//it was ==1
				      voltage[loop_count] = atof (tok);
				    }
				  //    if (token_count==3){
				  //           current=atof(tok);
				  //    }
				  token_count++;
				  tok = strtok (NULL, " ");
				}
				
				
				if ((fabs (prev_voltage[loop_count] - voltage[loop_count]) > difference))
				{
				  printf ("changed loop_count=%d\n",loop_count);
				  changed[loop_count]=1;
				  external_loop_count++;
				}
			    }
			    else
			    {
			      printf("nothing was written on test_voltage_change.txt\n");
			    }
			}
		        else
			{
			  printf ("could not locate file\n");
			  return -1;
			}
		
		      fclose (file);
		    }
		  //  printf("not entered loop_count=%d\n",loop_count);
		   loop_count++;
		 }

	    }
	  else
	    {
	      printf ("sccBmc command failed....trying again\n");
	    }

    }

  gettimeofday (&ttf, NULL);
  time = (ttf.tv_sec - tts.tv_sec) + (ttf.tv_usec - tts.tv_usec) * 0.000001;
  return time;
}
