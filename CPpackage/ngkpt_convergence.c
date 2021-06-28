/* Last updated 03/26/2020 by Erdong */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#define n_dataset_max 100
#define Ha2meV 27211.3825

int main (int argc, char * argv[]){
  //declare variables
  int n_dataset=0;
  double etotal[n_dataset_max];
  double etotal_diff[n_dataset_max];
  double ediff_per_atom[n_dataset_max];
  int ngkpt[n_dataset_max][3];
  int nkpt[n_dataset_max];
  int natom;
  char etot_file_name[100];
  char abinitInName[100];
  char abinitOutName[100];
  char dummy[100], dummy2[100];
  char ngkpt_name[n_dataset_max][100];
  char nkpt_name[n_dataset_max][100];
  char output_name[100];
  int i, j; //looper
  int stop=0; //used in read abinit output file
  int dummy2_len;
  int found=0;
  FILE * abinit_in_file;
  FILE * etot_file;
  FILE * output_file;
  FILE * abinit_out_file;
  //Parsing input arguments

  if(argc==4){
    strcpy(etot_file_name,argv[1]);
    strcpy(abinitInName,argv[2]);
    strcpy(abinitOutName,argv[3]);
    strcpy(output_name,"ngkpt_ETOT_CONVERGENCE");
    if(access(abinitInName,F_OK) != -1){
      printf("Abinit old input file name: %s \n", abinitInName);
    }else{
      printf("Abinit old input file name: %s does not exist. Exiting this program. \n", abinitInName);
      exit(0);
    }
    if(access(abinitOutName,F_OK) != -1){
      printf("Abinit output file name: %s \n", abinitOutName);
    }else{
      printf("Abinit output file name: %s does not exist. Exiting this program. \n", abinitOutName);
      exit(0);
    }
  }else{
    printf("Usage: <etotal energy file name> <ABINIT input file name> <ABINIT output name>\n");
    exit(0);
  }
  
  //Now read the input file
  printf("Begin reading input\n"); 
  abinit_in_file = fopen(abinitInName,"r");
  while(fscanf(abinit_in_file,"%s",dummy) != EOF ){
    //if(strcmp(dummy,"FILE----!")==0) break;
    //printf("dummy is: %s \n",dummy);
    if(strcmp(dummy,"natom") == 0){
      fscanf(abinit_in_file,"%d",&natom);
      found++;
      //printf("NATOM = %d \n", natom);
    }
    if(strcmp(dummy,"!") != 0){
      fscanf(abinit_in_file,"%s",dummy2);
      dummy2_len=strlen(dummy2);
      if(strcmp(dummy2,"ndtset") == 0){
        fscanf(abinit_in_file,"%d",&n_dataset);
        found++;
        //printf("ndtset = %d \n", n_dataset);
      }else{
        fseek(abinit_in_file,-dummy2_len,SEEK_CUR);
      }
    }
    for(i=0;i<n_dataset;i++){
      sprintf(ngkpt_name[i],"ngkpt%d",i+1);
      if(strcmp(dummy,ngkpt_name[i]) == 0){
        fscanf(abinit_in_file,"%d %d %d",&ngkpt[i][0],&ngkpt[i][1],&ngkpt[i][2]);
      found++;
      }
    }
    if(found==n_dataset+2) break;
  }
  fclose(abinit_in_file);

  //Now read the output file
  printf("Begin reading output\n");
  abinit_out_file = fopen(abinitOutName,"r");
  found=0;
  while(fscanf(abinit_in_file,"%s",dummy) != EOF ){
    //if(strcmp(dummy,"FILE----!")==0) break;
    //printf("dummy is: %s \n",dummy);
    for(i=0;i<n_dataset;i++){
      sprintf(nkpt_name[i],"nkpt%d",i+1); 
      if(strcmp(dummy,nkpt_name[i]) == 0){
        fscanf(abinit_out_file,"%d",&nkpt[i]);
        found++;
      }
    }
    if(found==n_dataset) break;
  }
  fclose(abinit_in_file);



  // Now read the total energy file
  printf("Begin reading ETOT\n");
  etot_file = fopen(etot_file_name,"r");
  for(i=0;i<n_dataset;i++){
    fscanf(etot_file,"%s %lf",&dummy, &etotal[i]);
    //printf("etotal %d: %lf \n",i,etotal[i]);
  }
  fclose(etot_file);

  //Now processing data
  
  for(i=0;i<n_dataset-1;i++){
    //printf("i=%d\n",i);
    etotal_diff[i]=fabs(etotal[i+1]-etotal[i]);
    ediff_per_atom[i]=etotal_diff[i]/natom;
  }
  
  output_file = fopen(output_name,"w");
  fprintf(output_file,"Your original calculation has %d ngkpt (dataset). \n",n_dataset);
  fprintf(output_file,"Therefore, the convergence of the last %d ngkpt will be tested. \n\n",n_dataset-1);
  for(i=0;i<n_dataset-1;i++){
    fprintf(output_file,"For ngkpt grid %d: %d %d %d, there are %d kpoints, etotal = %.8lf Ha, convergence: %.8lf Ha per atom or %.8lf meV per atom\n",i+2,ngkpt[i][0],ngkpt[i][1],ngkpt[i][2],nkpt[i],etotal[i],ediff_per_atom[i],ediff_per_atom[i]*Ha2meV);
  }
  fprintf(output_file,"\nSelect your ngkpt accordingly.\n");
  fprintf(output_file,"\nEND OF FILE");
  printf("DONE!\n");
  fclose(output_file);
}
