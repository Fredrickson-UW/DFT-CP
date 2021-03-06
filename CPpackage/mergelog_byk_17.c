/* Last updated 07/31/2019 by Erdong */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define NGX_MAX 500
#define NGY_MAX 500
#define NGZ_MAX 500
#define NATOMS_MAX 300
#define K_MAX 100
double XSF_sum_ds1[NGX_MAX][NGY_MAX][NGZ_MAX];
double XSF_sum_ds3[NGX_MAX][NGY_MAX][NGZ_MAX];
double LOG_sum_ds1[NATOMS_MAX][5];
double LOG_sum_ds3[NATOMS_MAX][5];


struct XSFfile{
  char dummy[100];
  double cella_x, cella_y, cella_z;
  double cellb_x, cellb_y, cellb_z;
  double cellc_x, cellc_y, cellc_z;
  int atomicno [NATOMS_MAX];
  double Xcart [NATOMS_MAX];
  double Ycart [NATOMS_MAX];
  double Zcart [NATOMS_MAX];
  int NIONS;
  int NGX, NGY, NGZ;
  double grid[NGX_MAX][NGY_MAX][NGZ_MAX];
} inputXSF;


void finish_line (FILE * f3){
    int status=0;
    char dummy;
    while (status == 0){
       dummy=getc(f3);
       if((dummy==EOF) || (dummy==10)){
           status=1;
       }
    }
}

void readXSF(struct XSFfile * inputXSF, FILE * f2, int ndt){
   /*declare variables*/
   int i,j,k;
   /*double volume1, volume2, volume3;*/
   finish_line(f2);
   finish_line(f2);
   finish_line(f2);
   fscanf(f2,"%lf %lf %lf",&inputXSF->cella_x,&inputXSF->cella_y,&inputXSF->cella_z);
   fscanf(f2,"%lf %lf %lf",&inputXSF->cellb_x,&inputXSF->cellb_y,&inputXSF->cellb_z);
   fscanf(f2,"%lf %lf %lf",&inputXSF->cellc_x,&inputXSF->cellc_y,&inputXSF->cellc_z);
   /*printf("%lf %lf %lf \n",inputXSF->cellb_x,inputXSF->cellb_y,inputXSF->cellb_z);*/
   finish_line(f2);
   finish_line(f2);
   fscanf(f2,"%d",&inputXSF->NIONS);
   /*printf("%d \n",inputXSF->NIONS);*/
   finish_line(f2);
   /*missing cell volume calculation... to be finished if needed*/
   for (i=0; i<inputXSF->NIONS; i++){
      fscanf(f2,"%d %lf %lf %lf",&inputXSF->atomicno[i],&inputXSF->Xcart[i],&inputXSF->Ycart[i],&inputXSF->Zcart[i]);
      /*printf("%d %lf %lf %lf \n",inputXSF->atomicno[i],inputXSF->Xcart[i],inputXSF->Ycart[i],inputXSF->Zcart[i]);*/
      finish_line(f2);
   }
   finish_line(f2);
   for (i=0; i<inputXSF->NIONS; i++){
      fscanf(f2,"%d %lf %lf %lf",&inputXSF->atomicno[i],&inputXSF->Xcart[i],&inputXSF->Ycart[i],&inputXSF->Zcart[i]);
      /*printf("%d %lf %lf %lf \n",inputXSF->atomicno[i],inputXSF->Xcart[i],inputXSF->Ycart[i],inputXSF->Zcart[i]);*/
      finish_line(f2);
   }
   finish_line(f2);
   finish_line(f2);
   finish_line(f2);
   fscanf(f2,"%d %d %d", &inputXSF->NGX,&inputXSF->NGY,&inputXSF->NGZ);
   finish_line(f2);
   finish_line(f2);

   finish_line(f2);
   finish_line(f2);
   finish_line(f2);

   for(i=0; i<inputXSF->NGZ; i++) {
      for(j=0; j<inputXSF->NGY; j++) {
         for(k=0; k<inputXSF->NGX; k++) {
              fscanf(f2,"%lf",&inputXSF->grid[i][j][k]);
              /*printf("%.9lf \n",inputXSF->grid[i][j][k]);*/
              if(ndt == 1){
                 XSF_sum_ds1[i][j][k]=XSF_sum_ds1[i][j][k]+inputXSF->grid[i][j][k];
              }
              if(ndt == 3){
                 XSF_sum_ds3[i][j][k]=XSF_sum_ds3[i][j][k]+inputXSF->grid[i][j][k];
              }
         }
      }
   }
   /*printf("%lf \n",inputXSF->grid[0][0][0]);*/
}

void readLOG(FILE * f6, int ndt, int NATOMS){
   /*declare variables*/
   int i;
   double NL_tot, NL_s, NL_p, NL_d, NL_f;
   char dummy1[500];
   char dummy2[500];
   char NL_tot_str[500];
   char NL_s_str[500];
   char NL_p_str[500];
   char NL_d_str[500];
   char NL_f_str[500];
   /*skip the one line header*/
   finish_line(f6);
   /*Reading data*/
   for(i=0; i<NATOMS; i++){
      fscanf(f6,"%s %s %s %s %s %s %s",&dummy1,&dummy2,&NL_tot_str,&NL_s_str,&NL_p_str,&NL_d_str,&NL_f_str);
      /*printf("Read from file: for atom %i, dummy = %s, dummy = %s, NL_tot = %s NL_s = %s, NL_p = %s, NL_d = %s, NL_f = %s \n", (i+1),NL_tot_str,NL_s_str,NL_p_str,NL_d_str,NL_f_str);*/
      sscanf(NL_tot_str,"%lf",&NL_tot);
      sscanf(NL_s_str,"%lf",&NL_s);
      sscanf(NL_p_str,"%lf",&NL_p);
      sscanf(NL_d_str,"%lf",&NL_d);
      sscanf(NL_f_str,"%lf",&NL_f);
      if(ndt==1){
        LOG_sum_ds1[i][0]=LOG_sum_ds1[i][0]+NL_tot;
        LOG_sum_ds1[i][1]=LOG_sum_ds1[i][1]+NL_s;
        LOG_sum_ds1[i][2]=LOG_sum_ds1[i][2]+NL_p;
        LOG_sum_ds1[i][3]=LOG_sum_ds1[i][3]+NL_d;
        LOG_sum_ds1[i][4]=LOG_sum_ds1[i][4]+NL_f;
      }
      if(ndt==3){
        LOG_sum_ds3[i][0]=LOG_sum_ds3[i][0]+NL_tot;
        LOG_sum_ds3[i][1]=LOG_sum_ds3[i][1]+NL_s;
        LOG_sum_ds3[i][2]=LOG_sum_ds3[i][2]+NL_p;
        LOG_sum_ds3[i][3]=LOG_sum_ds3[i][3]+NL_d;
        LOG_sum_ds3[i][4]=LOG_sum_ds3[i][4]+NL_f;
      } 
   }
}
void writeLOG(FILE * f8, int ndt, int NATOMS){
  /*declare variables*/
  int i;
  char word[10];
  /*done*/
  sprintf(word,"  atom ");
  /*Write file header*/
  fprintf(f8,"Nonlocal psp contribution to total energy from DS%i_WFK generated by ***mergeXSF*** \n",ndt);
  if (ndt==1){
    for (i=0; i<NATOMS ;i++){
      fprintf(f8,"%s%i: %20.17lf  %20.17lf  %20.17lf  %20.17lf  %20.17lf\n",word,(i+1),LOG_sum_ds1[i][0],LOG_sum_ds1[i][1],LOG_sum_ds1[i][2],LOG_sum_ds1[i][3],LOG_sum_ds1[i][4]);
    }
  }
  if (ndt==3){
    for (i=0; i<NATOMS ;i++){
      fprintf(f8,"%s%i: %20.17lf  %20.17lf  %20.17lf  %20.17lf  %20.17lf\n",word,(i+1),LOG_sum_ds3[i][0],LOG_sum_ds3[i][1],LOG_sum_ds3[i][2],LOG_sum_ds3[i][3],LOG_sum_ds3[i][4]);
    }
  }
}

void writeXSF(struct XSFfile * inputXSF, FILE * f4, int ndt)
{
   /* declare variables*/
   int i, j, k;
   int line_counter;
   fprintf(f4, " DIM-GROUP\n");
   fprintf(f4, " 3 1\n");
   fprintf(f4, " PRIMVEC\n");
   fprintf(f4, "%19.10lf%19.10lf%19.10lf\n",inputXSF->cella_x,inputXSF->cella_y,inputXSF->cella_z);
   fprintf(f4, "%19.10lf%19.10lf%19.10lf\n",inputXSF->cellb_x,inputXSF->cellb_y,inputXSF->cellb_z);
   fprintf(f4, "%19.10lf%19.10lf%19.10lf\n",inputXSF->cellc_x,inputXSF->cellc_y,inputXSF->cellc_z);
   fprintf(f4, " PRIMCOORD\n");
   fprintf(f4, "%12d%3d\n", inputXSF->NIONS, 1);
   for (i=0; i<inputXSF->NIONS; i++){
      fprintf(f4, "%9d%20.10lf%20.10lf%20.10lf\n", inputXSF->atomicno[i], inputXSF->Xcart[i], inputXSF->Ycart[i], inputXSF->Zcart[i]);
   }
   fprintf(f4, " ATOMS\n");
   for (i=0; i<inputXSF->NIONS; i++){
      fprintf(f4, "%9d%20.10lf%20.10lf%20.10lf\n", inputXSF->atomicno[i], inputXSF->Xcart[i], inputXSF->Ycart[i], inputXSF->Zcart[i]);
   }
   fprintf(f4, " BEGIN_BLOCK_DATAGRID3D\n");
   fprintf(f4, " datagrids\n");
   fprintf(f4, " DATAGRID_3D_DENSITY\n");
   fprintf(f4, "%12d%12d%12d\n", inputXSF->NGX,inputXSF->NGY,inputXSF->NGZ);
   fprintf(f4, " 0.0 0.0 0.0\n");
   fprintf(f4, "%19.10lf%19.10lf%19.10lf\n",inputXSF->cella_x,inputXSF->cella_y,inputXSF->cella_z);
   fprintf(f4, "%19.10lf%19.10lf%19.10lf\n",inputXSF->cellb_x,inputXSF->cellb_y,inputXSF->cellb_z);
   fprintf(f4, "%19.10lf%19.10lf%19.10lf\n",inputXSF->cellc_x,inputXSF->cellc_y,inputXSF->cellc_z);
   for(i=0; i<inputXSF->NGZ; i++) {
      for(j=0; j<inputXSF->NGY; j++) {
         for(k=0; k<inputXSF->NGX; k++) {
            if(ndt==1){
               fprintf(f4,"%20.10lf", XSF_sum_ds1[i][j][k]);
            }
            if(ndt==3){
               fprintf(f4,"%20.10lf", XSF_sum_ds3[i][j][k]);
            }
            line_counter++;
            if(line_counter==6){
               fprintf(f4,"\n");
               line_counter=0;
            }
         }     
      }
   }
   fprintf(f4, " END_DATAGRID_3D\n");
   fprintf(f4, " END_BLOCK_DATAGRID3D\n");
}

   
int main (int argc, char * argv[])
{
    /* declare variables */
    char outfilename[100];
    char nset_str[100];
    char CurrentXSFname[100];
    char OutputXSFname[100];
    char CurrentLOGname[100];
    char OutputLOGname[100];
    char NATOM_str1[5], NATOM_str2[5];
    int Currentset;
    int NATOM;
    char XSF1[25];
    char XSF2[25];
    char XSF3[25];
    char XSF4[25];
    char XSF5[25];
    char LOG1[25];
    char LOG2[25];
    char LOG3[25];
    char LOG4[25];
    char LOG5[25];
    char dummy[100];
    char sampleLOGname[25];
    int nset;
    int i,j,k; /*loopers*/
    int NATOM_len;
    FILE * CurrentXSF;
    FILE * outputXSF1;
    FILE * outputXSF2;
    FILE * CurrentLOG;
    FILE * outputLOG1;
    FILE * outputLOG2;
    FILE * sampleLOG;
    strcpy(XSF1,"abinitout_DS");
    strcpy(XSF2,"_NL_kpt");
    strcpy(XSF3,".xsf");
    /*Reading argument*/
    
    NATOM=0;
    memset(NATOM_str2,0,5);
    if(argc==3){
       strcpy(nset_str,argv[1]);
       sscanf(nset_str,"%d",&nset);
       printf("Number of k points: %i. \n", nset);
       strcpy(outfilename,argv[2]);
    }
    else{
       printf("Usage: <number of k points> <ABINIT output file name> \n");
       exit(0);
    }
    /* Reading XSFs */
    for (i = 0; i<NATOMS_MAX; i++){
      for (j = 0; j<5; j++){
        LOG_sum_ds1[i][j]=0;
        LOG_sum_ds3[i][j]=0;
      }
    }
    sprintf(sampleLOGname,"%s_o_DS%d_NL_kpt%d.log",outfilename,1,1);
    sampleLOG = fopen(sampleLOGname,"r");
    printf("Sample log file: %s opened! \n", sampleLOGname);
    while(fscanf(sampleLOG,"%s",dummy)==1){
      printf("Dummy: %s \n",dummy);
      if(strcmp(dummy,"atom")==0){
        fscanf(sampleLOG,"%s",NATOM_str1);
        printf("NATOM_str found: %s\n", NATOM_str1);
      }
    }
    NATOM_len=strlen(NATOM_str1);
    for(i=0;i<NATOM_len-1;i++){
      NATOM_str2[i]=NATOM_str1[i];
    }
    printf("NATOM_str2: %s\n",NATOM_str2);
    sscanf(NATOM_str2,"%d",&NATOM);
    fclose(sampleLOG);
    printf("Number of atoms found: %i. \n", NATOM);
    /* Now begin to merge logs */
    printf("Now begin to merge log files. \n");
    printf("Total number of atoms = %i \n", NATOM);
    for (i = 1; i<= 3; i = i+2){
       for ( j = 1; j <= nset; j++){
          sprintf(CurrentLOGname,"%s_o_DS%d_NL_kpt%d.log",outfilename,i,j);
          printf("Reading log file name:%s \n", CurrentLOGname);
          CurrentLOG = fopen(CurrentLOGname,"r");
          if(CurrentLOG==NULL){
             printf("Files %s not found. Please check the completeness of your nonlocal calculation. \n", CurrentLOGname);
             exit(0);
          }
          readLOG(CurrentLOG,i,NATOM);
          printf("Done reading %s. \n", CurrentLOGname);
       }
    }
    strcpy(LOG4,".out_o_DS");
    strcpy(LOG5,"_NL.log");
    sprintf(OutputLOGname,"%s_o_DS%d_NL.log",outfilename,1);
    outputLOG1 = fopen(OutputLOGname,"w");
    printf("Writing log file:%s \n",OutputLOGname);
    writeLOG(outputLOG1,1,NATOM);
    fclose(outputLOG1);
    sprintf(OutputLOGname,"%s_o_DS%d_NL.log",outfilename,3);
    outputLOG2 = fopen(OutputLOGname,"w");
    printf("Writing log file:%s \n",OutputLOGname);
    writeLOG(outputLOG2,3,NATOM);
    fclose(outputLOG2);
}
   
