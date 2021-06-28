/*Last updated 04/10/2020 by Erdong*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

int main (int argc, char * argv[]){
  //declare variables
  int ngfft[3];
  char dummy[50];
  int status=0;
  FILE * info_file;
  FILE * status_file;
  int i=0;
  int divisiable=0;
  memset(dummy,0,50);
  info_file=fopen("ngfft_Info","r");
  status_file=fopen("ngfft_status","w");
  while(fscanf(info_file,"%s",&dummy)==1){
    //printf("dummy is: %s\n",dummy);
    if(strcmp(dummy,"ngfft")==0){
      //printf("ngfft accpted.\n");
      fscanf(info_file,"%d %d %d",&ngfft[0],&ngfft[1],&ngfft[2]);
    }
  }
  for(i=0;i<3;i++){
    if((ngfft[i] % 12) == 0) divisiable++;
  }

  if(divisiable==3){
    //printf("ngfft not accepted.\n");
    status=1;
  }else{
    status=2;
  }

  if(status==1){
    fprintf(status_file,"%d\n",status);
    fprintf(status_file,"ngfft accepted\n");
  }else if(status==2){
    fprintf(status_file,"%d\n",status);
    fprintf(status_file,"ngfft not accepted\n");
  }

}
