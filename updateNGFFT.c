/*Last updated 04/10/2020 by Erdong*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

char * get_in_name (FILE * files_file){
  char in_file_name_read[100];
  char * in_file_name;
  //printf("In function: get_in \n");
  fscanf(files_file,"%s",in_file_name_read);
  //printf("in_file_name_read: %s \n",in_file_name_read);
  in_file_name=in_file_name_read;
  printf("in_file_name copied: %s \n",in_file_name);
  return in_file_name;
}

void update_in_file (char * old_in_name_pt, char * new_in_name_pt, long * ngfft){
  char old_in_name[100], new_in_name[100];
  char dummy1[200], word1[20], word2[20];
  int dummy1_len;
  int ngfft_adjust=0;
  FILE * old_in; 
  FILE * new_in;
  printf("Inside function: update_in_file\n");

  strcpy(old_in_name,old_in_name_pt);
  strcpy(new_in_name,new_in_name_pt);

  old_in=fopen(old_in_name,"r");
  new_in=fopen(new_in_name,"w");

  while(fgets(dummy1,200,old_in) != NULL){
    //printf("dummy1 is: %s \n",dummy1);
    //dummy1_len=strlen(dummy1);

    sscanf(dummy1,"%s %s",word1,word2);

    if(ngfft_adjust==0){
      if((strcmp(word1,"!") == 0) && (strcmp(word2,"ngfft") == 0)){
        //printf("! ngfft line found  \n");
        ngfft_adjust=2;
      }else if(strcmp(word1,"ngfft") == 0){
        //printf("ngfft line found  \n");
        ngfft_adjust=2;
      }
    }

    if(ngfft_adjust == 2){
      printf("Performing ngfft adjust. \n");
      if(ngfft[0]==1 && ngfft[1]==1 && ngfft[2]==1){
        fprintf(new_in,"!    ngfft 72 72 72      ! Updated ngfft value\n",ngfft[0],ngfft[1],ngfft[2]);
      }else{
        fprintf(new_in,"     ngfft %ld %ld %ld      ! Updated ngfft value\n",ngfft[0],ngfft[1],ngfft[2]);
      }
      ngfft_adjust=1;
    }else{
      fprintf(new_in,"%s",dummy1);  
    }
    //printf("Next fget. \n"); 
  }
  printf("Exiting function update.\n");
}

int main (int argc, char * argv[]){
  //declare variables
  long ngfft[3];
  char files_file_name[100];
  char * in_file_name_pt;
  char in_file_name[100];
  char in_file_name_base[100];
  char in_file_name_new[100];
  int in_file_name_len;
  FILE * files_file;
  char * temp;
  char dummy[100];
  int i;
  memset(in_file_name_base,0,100);

  if(argc==5){
    strcpy(files_file_name,argv[1]);
    ngfft[0]=strtol(argv[2],&temp,10);
    ngfft[1]=strtol(argv[3],&temp,10);
    ngfft[2]=strtol(argv[4],&temp,10);
    for(i=0;i<3;i++){
      if(ngfft[i]==0){
        printf("Looks like ngfft[%d] you specify: %s is 0 or is not a number. This is not allowed. \n",i+1,argv[i+3]);
        exit(0);
      }else{
        printf("ngfft[%d]=%ld  ",i+1,ngfft[i]);
      }
    }
    printf("\n");
  }else{
    printf("updateNGFFT usage: <.files file name> <ngfft[1]> <ngfft[2]> <ngfft[3]>\n");
    printf("Enter 1 1 1 for the three ngfft values to get ABINIT default. \n");
    exit(0);
  }
  files_file=fopen(files_file_name,"r");
  in_file_name_pt=get_in_name(files_file);
  strcpy(in_file_name,in_file_name_pt);
  printf("In file name: %s \n",in_file_name);
  fclose(files_file);
  in_file_name_len=strlen(in_file_name);
  for(i=0;i<in_file_name_len;i++){
    if(in_file_name[i] != '.') in_file_name_base[i]=in_file_name[i];
  } 
  printf("in_file_name_base: %s \n",in_file_name_base);
  sprintf(in_file_name_new,"%s_new.in",in_file_name_base);
  //strcpy(in_file_name_new,"test_output.in");
  update_in_file(in_file_name,in_file_name_new,ngfft);
  printf("DONE generating the new input file. \n");
  printf("Now cleaning up the files. \n");
  remove(in_file_name);
  rename(in_file_name_new,in_file_name);
  printf("DONE \n");
}
