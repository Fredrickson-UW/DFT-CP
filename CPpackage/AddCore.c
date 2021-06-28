#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define PI 3.14159265
#define NPOINTS_MAX 30000
#define NGX_MAX 400
#define NATOMS_MAX 300
#define LONGEST_FILENAME 100

struct XSFfile{
  char systemname[100];
  double cellvolume;
  double scale_xyz;
  double cella_x, cella_y, cella_z;
  double cellb_x, cellb_y, cellb_z;
  double cellc_x, cellc_y, cellc_z;
  int atomicno [NATOMS_MAX];
  double Xcart [NATOMS_MAX];
  double Ycart [NATOMS_MAX];
  double Zcart [NATOMS_MAX];
  int NIONS;
  int NGX, NGY, NGZ;
  double grid[NGX_MAX][NGX_MAX][NGX_MAX];
  double VoxelV;
  double NELECT;
} den;

FILE * inputfile;
FILE * outputfile;

struct radial_profile{
     double r[120][NPOINTS_MAX];
     double rho[120][NPOINTS_MAX];
     int use[120];
     int npoints[120];
} cores;

void finish_line (FILE * f3)
{
     int cont=0;
     char check;
     while (cont==0) {
         check=getc(f3);
	 if ((check==10)||(check==EOF)) cont=1;	    
     }
}

void copy_line (FILE * f3,FILE * f2)
{
     int cont=0;
     char check;
     while (cont==0) {
         check=getc(f3);
	 putc(check,f2);
	 if ((check==10)||(check==EOF)) cont=1;	    
     }
}


void read_line (FILE * f3, char *newstring)
{
     int cont=0;
     int counter=0;
     char check;
     while (cont==0) {
         check=getc(f3);
         if ((check==10)||(check==EOF)) cont=1;
         *(newstring+counter)=check;
         counter++;         
     }
     *(newstring+counter-1)=0;
}



void readXSF(struct XSFfile * XSFIN, FILE *f2)
{
     FILE *f3;
     char jk[100];
     int j,jx, jy, jz;
     int line_counter;
     double CHGvalue;
     char systemname[100];
     finish_line(f2);
     finish_line(f2);
     finish_line(f2);
     fscanf(f2,"%lf %lf %lf",&XSFIN->cella_x, &XSFIN->cella_y, &XSFIN->cella_z);
     fscanf(f2,"%lf %lf %lf",&XSFIN->cellb_x, &XSFIN->cellb_y, &XSFIN->cellb_z);
     fscanf(f2,"%lf %lf %lf",&XSFIN->cellc_x, &XSFIN->cellc_y, &XSFIN->cellc_z);
     finish_line(f2);
     finish_line(f2);
     fscanf(f2,"%d",&XSFIN->NIONS);
     finish_line(f2);
     XSFIN->cellvolume=(XSFIN->cella_x*(XSFIN->cellb_y*XSFIN->cellc_z-XSFIN->cellb_z*XSFIN->cellc_y)-XSFIN->cella_y*(XSFIN->cellb_x*XSFIN->cellc_z-XSFIN->cellb_z*XSFIN->cellc_x)+XSFIN->cella_z*(XSFIN->cellb_x*XSFIN->cellc_y-XSFIN->cellb_y*XSFIN->cellc_x));
     for(j=0; j<XSFIN->NIONS; j++) { 
        fscanf(f2,"%d %lf %lf %lf",&XSFIN->atomicno[j],&XSFIN->Xcart[j],&XSFIN->Ycart[j],&XSFIN->Zcart[j]);
        finish_line(f2);
     }
     finish_line(f2);
     for(j=0; j<XSFIN->NIONS; j++) { 
        fscanf(f2,"%d %lf %lf %lf",&XSFIN->atomicno[j],&XSFIN->Xcart[j],&XSFIN->Ycart[j],&XSFIN->Zcart[j]);
        finish_line(f2);
     }
     finish_line(f2);
     finish_line(f2);
     finish_line(f2);

     fscanf(f2,"%d %d %d", &XSFIN->NGX, &XSFIN->NGY, &XSFIN->NGZ);
     finish_line(f2);
     finish_line(f2);
     fscanf(f2,"%lf %lf %lf",&XSFIN->cella_x, &XSFIN->cella_y, &XSFIN->cella_z);
     fscanf(f2,"%lf %lf %lf",&XSFIN->cellb_x, &XSFIN->cellb_y, &XSFIN->cellb_z);
     fscanf(f2,"%lf %lf %lf",&XSFIN->cellc_x, &XSFIN->cellc_y, &XSFIN->cellc_z);
     for(jz=0;jz<XSFIN->NGZ; jz++) {
         for(jy=0;jy<XSFIN->NGY; jy++) {
            for(jx=0;jx<XSFIN->NGX; jx++) {
                 fscanf(f2,"%lf",&CHGvalue);
                 XSFIN->grid[jx][jy][jz]=CHGvalue;
            }
         }
     }
     XSFIN->VoxelV=XSFIN->cellvolume/((XSFIN->NGX-1)*(XSFIN->NGY-1)*(XSFIN->NGZ-1))/pow(0.52917720859,3.0);
}



void loadProfiles(struct XSFfile * XSFIN, struct radial_profile *CORE)
{
     FILE *f3;
     char firstline[300];
     char profilename[LONGEST_FILENAME];     
     int j,k,atomicno, stop, check;
     int line_counter;
     double integral;
     double CHGvalue;
     double r, rho, drho_dr,d2rho_dr2;
     for(j=0;j<130;j++) CORE->use[j]=0;
     for(j=0;j<XSFIN->NIONS;j++) {
            atomicno=XSFIN->atomicno[j];
	    CORE->use[atomicno]=1;
     }  
     for(j=0;j<130;j++) {
         if(CORE->use[j]==1) {
	     printf("Enter name of file containing radial profile for atomic number %d:  ",j);
	     scanf("%s",profilename);
	     f3=fopen(profilename,"r");
	     if(f3==NULL) {
	         printf("Could not open file %s.\n",profilename);
	         exit(0);
	     }
	     stop=0;
	     line_counter=0;
             integral=0.0;
	     while(stop==0) {
	         check=fscanf(f3,"%lf %lf",&r, &rho);
		 if((check==0)||(check==EOF)) stop=1;
		 else {
		     CORE->r[j][line_counter]=r;
		     CORE->rho[j][line_counter]=rho;
		     line_counter++;
		 }
	     }
             for(k=1;k<line_counter-1;k++) {
                    integral+=CORE->r[j][k]*CORE->r[j][k]*CORE->rho[j][k]*(CORE->r[j][k+1]-CORE->r[j][k-1])/2.0;
             }
             CORE->npoints[j]=line_counter;
	     printf("%d data points read. Approximately %lf electrons.\n",line_counter,4*PI*integral);
	     fclose(f3);
	 }
     }
}

double readProfile(struct radial_profile *CORE, int atomicno, double r)
{
     int j,stop;
     int line_counter;
     double r1, r2,rho, rho1, rho2,deltaR;
     stop =0;
     line_counter=0;
     /* Perform linear interpolation between available data points */
     if(r<=CORE->r[atomicno][0]) rho = CORE->rho[atomicno][0];
     else {
       while (stop==0) {
            line_counter++;
            if(r<=CORE->r[atomicno][line_counter]) {
	        stop = 1;
	    }
            if(line_counter==CORE->npoints[atomicno]) stop=1;
       }
       r1 = CORE->r[atomicno][line_counter-1];
       r2 = CORE->r[atomicno][line_counter];
       rho1 = CORE->rho[atomicno][line_counter-1];
       rho2 = CORE->rho[atomicno][line_counter];
       deltaR = r2-r1;
       rho = ((r2-r)/deltaR)*rho1 + ((r-r1)/deltaR)*rho2;
       if(line_counter==CORE->npoints[atomicno]) {
                /* Distance is outside of range of profile */
                rho = 0.0;
       }
     }
     return(rho);
}



void outputXSF(struct XSFfile * XSFIN, struct XSFfile * XSFOUT, FILE *f2)
{
     int jx,jy,jz;
     int j;
     int line_counter;
     fprintf(f2, " DIM-GROUP\n");
     fprintf(f2, " 3  1\n");
     fprintf(f2, " PRIMVEC\n");
     fprintf(f2, "%17.10lf%19.10lf%19.10lf\n",XSFIN->cella_x,XSFIN->cella_y,XSFIN->cella_z);
     fprintf(f2, "%17.10lf%19.10lf%19.10lf\n",XSFIN->cellb_x,XSFIN->cellb_y,XSFIN->cellb_z);
     fprintf(f2, "%17.10lf%19.10lf%19.10lf\n",XSFIN->cellc_x,XSFIN->cellc_y,XSFIN->cellc_z);
     fprintf(f2, " PRIMCOORD\n");
     fprintf(f2, "%12d%3d\n",XSFIN->NIONS,1);
     for(j=0;j<XSFIN->NIONS;j++) {
         fprintf(f2, "%9d%20.10lf%20.10lf%20.10lf\n",XSFIN->atomicno[j],XSFIN->Xcart[j],XSFIN->Ycart[j],XSFIN->Zcart[j]);
     }
     fprintf(f2, " ATOMS\n");
     for(j=0;j<XSFIN->NIONS;j++) {
         fprintf(f2, "%9d%20.10lf%20.10lf%20.10lf\n",XSFIN->atomicno[j],XSFIN->Xcart[j],XSFIN->Ycart[j],XSFIN->Zcart[j]);
     }
     fprintf(f2," BEGIN_BLOCK_DATAGRID3D\n");
     fprintf(f2," datagrids\n");
     fprintf(f2," DATAGRID_3D_DENSITY\n");
     fprintf(f2,"%12d%12d%12d\n",XSFIN->NGX,XSFIN->NGY,XSFIN->NGZ);
     fprintf(f2," 0.0 0.0 0.0\n");
     fprintf(f2, "%17.10lf%19.10lf%19.10lf\n",XSFIN->cella_x,XSFIN->cella_y,XSFIN->cella_z);
     fprintf(f2, "%17.10lf%19.10lf%19.10lf\n",XSFIN->cellb_x,XSFIN->cellb_y,XSFIN->cellb_z);
     fprintf(f2, "%17.10lf%19.10lf%19.10lf\n",XSFIN->cellc_x,XSFIN->cellc_y,XSFIN->cellc_z);
     line_counter=0;
     for(jz=0;jz<XSFIN->NGZ; jz++) {
        for(jy=0;jy<XSFIN->NGY; jy++) {
           for(jx=0;jx<XSFIN->NGX; jx++) {
              line_counter++;
              fprintf(f2,"%20.5lf",XSFOUT->grid[jx][jy][jz]);
              if(line_counter==6) {
                  fprintf(f2,"\n");
                  line_counter=0;
              }
           }
        }
     }
     fprintf(f2," END_DATAGRID_3D\n");
     fprintf(f2," END_BLOCK_DATAGRID3D\n");

}


void addCores(struct XSFfile * XSFIN, struct radial_profile * CORE)
{
     int jx,jy,jz;
     int n_a, n_b, n_c;
     int j;
     int line_counter;
     int nvoxels,voxel_counter;
     double voxel_x,voxel_y,voxel_z;
     double xf,yf,zf;
     double r_max = 5.0;
     double rho,dist;
     printf("Adding core densities to map...\n");
     nvoxels = (XSFIN->NGZ)*(XSFIN->NGY)*(XSFIN->NGX);
     voxel_counter=0;
     for(jz=0;jz<XSFIN->NGZ; jz++) {
        for(jy=0;jy<XSFIN->NGY; jy++) {
           for(jx=0;jx<XSFIN->NGX; jx++) {
              voxel_counter++;
              for(n_a=-1;n_a<2; n_a++) {
                 xf = (1.0*jx)/(1.0*XSFIN->NGX-1.0)+1.0*n_a;
                 for(n_b=-1;n_b<2; n_b++) {
                    yf = (1.0*jy)/(1.0*XSFIN->NGY-1.0)+1.0*n_b;
                    for(n_c=-1;n_c<2; n_c++) {
                        zf = (1.0*jz)/(1.0*XSFIN->NGZ-1.0)+1.0*n_c;
                        voxel_x = xf*XSFIN->cella_x + yf*XSFIN->cellb_x + zf*XSFIN->cellc_x;
                        voxel_y = xf*XSFIN->cella_y + yf*XSFIN->cellb_y + zf*XSFIN->cellc_y;
                        voxel_z = xf*XSFIN->cella_z + yf*XSFIN->cellb_z + zf*XSFIN->cellc_z;
                        for(j=0;j<XSFIN->NIONS;j++) {
			    dist = pow((voxel_x-XSFIN->Xcart[j])*(voxel_x-XSFIN->Xcart[j])+(voxel_y-XSFIN->Ycart[j])*(voxel_y-XSFIN->Ycart[j])+(voxel_z-XSFIN->Zcart[j])*(voxel_z-XSFIN->Zcart[j]),0.5);
                            
			    if(dist < r_max) {
			          rho = readProfile(CORE,XSFIN->atomicno[j],dist/0.52917721);
				  XSFIN->grid[jx][jy][jz]+=rho;
			    }
			}

                    }
		 }
              }
           }
        }
        printf("\r%d%%",(voxel_counter*100)/nvoxels);
        fflush(stdout);
     }
     printf("\nDone!\n");
}


void outputCHGCAR(struct XSFfile * XSFIN, struct XSFfile * CHGCAROUT, FILE *f2, char name [100],char outfilename [100])
{
     int jx,jy,jz;
     int j;
     int line_counter;
     int stop=0,check;
     FILE * f3;
     char instring [1000];
     double xf,yf,zf;
     f3=fopen(outfilename,"r");
     if(f3==NULL) {
           printf("%s not found. Please supply name of _out file from ABINIT with \n",outfilename);
           printf("fractional coordinates of atoms prefaced by xred keyword. \n");
           exit(0);
     }
     fprintf(f2, "CHGCAR file converted from XSF file %s. \n",name);
     fprintf(f2, "%19.14lf\n",1.000);
     fprintf(f2, "%13.6lf%12.6lf%12.6lf\n",XSFIN->cella_x,XSFIN->cella_y,XSFIN->cella_z);
     fprintf(f2, "%13.6lf%12.6lf%12.6lf\n",XSFIN->cellb_x,XSFIN->cellb_y,XSFIN->cellb_z);
     fprintf(f2, "%13.6lf%12.6lf%12.6lf\n",XSFIN->cellc_x,XSFIN->cellc_y,XSFIN->cellc_z);
     fprintf(f2,"%3d\n",XSFIN->NIONS);
     fprintf(f2,"Direct\n");
     while(stop==0) {
        if(stop==0) check=fscanf(f3,"%s",instring);
        if(strcmp(instring,"xred")==0) stop=1;
        if(check==EOF) { 
             stop=1;
             XSFIN->NIONS=1;
             fprintf(f2,"%10.6lf%10.6lf%10.6lf\n",0.0,0.0,0.0);
        }
     }
     for(j=0;j<XSFIN->NIONS;j++) {
          if(check!=EOF) {
            fscanf(f3,"%lf %lf %lf",&xf,&yf,&zf);
            fprintf(f2,"%10.6lf%10.6lf%10.6lf\n",xf,yf,zf);
          }
     }
     fclose(f3);
     fprintf(f2,"\n");
     fprintf(f2,"%5d %4d %4d\n",XSFIN->NGX-1,XSFIN->NGY-1,XSFIN->NGZ-1);
     line_counter=0;
     for(jz=0;jz<XSFIN->NGZ-1; jz++) {
        for(jy=0;jy<XSFIN->NGY-1; jy++) {
           for(jx=0;jx<XSFIN->NGX-1; jx++) {
              line_counter++;
              fprintf(f2," %17.11lf",XSFIN->cellvolume*CHGCAROUT->grid[jx][jy][jz]/pow(0.52917720859,3.0));
              if(line_counter==5) {
                fprintf(f2,"\n");
                line_counter=0;
              }
           }
        }
     }
}



int main (int argc, char * argv[])
{
    FILE * f2;
    FILE * f4;
    char outfilename[100];
    char filename [LONGEST_FILENAME];
    char chgcarname[LONGEST_FILENAME];
    if (argc > 2) {
          strcpy(filename,argv[1]);
          strcpy(outfilename,argv[2]);
    }
    else {
       printf("Usage:  AddCore <XSF filename without extension> <abinitout file>\n");
       exit(0);
    }        
    /* CHECK THAT NEEDED FILES ARE PRESENT, AND LOAD DATA FOR DATASET 1*/
    strcat(filename,".xsf");
    f2=fopen(filename,"r");
    if(f2==NULL) {
      printf("file %s not found.\n", filename);
      exit(0);
    }
    readXSF(&den, f2);
    fclose(f2);
    f2=fopen("CHGCAR","w");
    outputCHGCAR(&den,&den,f2,filename,outfilename);
    fclose(f2);
    loadProfiles(&den,&cores);
    addCores(&den,&cores);
    f4=fopen("CHGCAR_sum","w");
    outputCHGCAR(&den,&den,f4,filename,outfilename);
    fclose(f4);
}
