/*Last updated 04/10/2020 by Erdong*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

int main (int argc, char * agrv[]){
  FILE * JOB_Info; 
  FILE * JOB_ID;
  char JOB_Info_full[100];
  char JOB_ID_only[100];
  int info_len;
  int i;
  //printf("Getting Job info\n");

  JOB_Info=fopen("JOBInfo","r");
  JOB_ID=fopen("JOBID","w");

  fscanf(JOB_Info,"%s",JOB_Info_full);
  info_len=strlen(JOB_Info_full);

  for(i=0;i<info_len;i++){
    if(JOB_Info_full[i] != '.'){
      JOB_ID_only[i]=JOB_Info_full[i];
    }
    else{
      break;
    }
  }

  //printf("JOB ID is: %s\n",JOB_Info_full);
  fprintf(JOB_ID,"%s",JOB_ID_only);
  //printf("DONE \n");
}
