/*     bin2xsf, part of the Fredrickson Group Chemical Pressure Package  

                  Copyright (C) 2012, by Daniel C. Fredrickson

                    Last modified:  Mar. 28, 2012

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

/*  NEWS  

    3/28/12  Added more flexibility as to what is included in P_remainder, and how it is calculated.  -DCF

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_sf_legendre.h>


#define CUT3D_COMMAND "cut3d"
#define PI 3.14159265
#define MAX_YAEHMOPFILES 10
#define LONGEST_FILENAME 100
#define DOSPOINTS_MAX 30000
#define NGX_MAX 305
#define NGY_MAX 305
#define NGZ_MAX 335
#define NATOMS_MAX 300

#define TBANDS_MAX 3000
#define BANDS_MAX 80
#define NORBS_MAX 3000
#define MAX_SYM 10000
#define KPOINTS_MAX 100
#define WAVES_MAX 6200
#define NTYPES_MAX 5

struct psp_file{
     double rloc[NTYPES_MAX];
     double rrs[NTYPES_MAX];
     double rrp[NTYPES_MAX];
     double rrd[NTYPES_MAX];
     double cc1[NTYPES_MAX];
     double cc2[NTYPES_MAX];
     double cc3[NTYPES_MAX];
     double cc4[NTYPES_MAX];
     double h[3][3][3][NTYPES_MAX];
     double h11s[NTYPES_MAX];
     double h22s[NTYPES_MAX];
     double h33s[NTYPES_MAX];
     double h11p[NTYPES_MAX];
     double h22p[NTYPES_MAX];
     double h33p[NTYPES_MAX];
     double h11d[NTYPES_MAX];
     double h22d[NTYPES_MAX];
     double h33d[NTYPES_MAX];
     double k11p[NTYPES_MAX];
     double k22p[NTYPES_MAX];
     double k33p[NTYPES_MAX];
     double k11d[NTYPES_MAX];
     double k22d[NTYPES_MAX];
     double k33d[NTYPES_MAX];
} psp;

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
  double grid[NGX_MAX][NGY_MAX][NGZ_MAX];
  double VoxelV;
  double NELECT;
  double epsatm[NATOMS_MAX];
  int epsatm_map[NATOMS_MAX];
  int epsatm_map_type[NATOMS_MAX];
  double nelectrons;
} den0;

struct wfk_file {
   int nbands;
   int nkpts;
   int npw[2][KPOINTS_MAX];
   int nspinor[2][KPOINTS_MAX];
   int nband[2][KPOINTS_MAX];
   double kg[WAVES_MAX][3];
   double eigen[BANDS_MAX];
   double occ[BANDS_MAX];
   double cg[WAVES_MAX][2];
} wfk1;


FILE * inputfile;
FILE * outputfile;
char filename [LONGEST_FILENAME];

void finish_line (FILE * f3)
{
     int cont=0;
     char check;
     while (cont==0) {
         check=getc(f3);
	 if ((check==10)||(check==EOF)) cont=1;	    
     }
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
              fprintf(f2,"%20.10lf  ",XSFOUT->grid[jx][jy][jz]);
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





void read_calc_nonlocal(char filename[200], struct XSFfile * XSFIN)
{
    FILE * f2;
    int k;
    char codvsn [110];
    char title [132];
    char test;
    int headform;
    int fform;
    int bandtot;    
    int date;
    int intxc;               
    int ixc;
    int natom;
    int atomno,pw1,pw2;
    int ngfftx;
    int ngffty;
    int ngfftz;
    int nkpt;
    int npsp;
    int j=0;
    int i,m,l;
    int jx,jy,jz;
    int stop=0;
    int nspden;
    int nsppol;
    int nsym;
    int ntypat;
    int occopt;
    int pertcase;
    int usepaw;
    double ecut;
    double ecutdg;
    double ecutsm;
    double ecut_eff;
    double qptnx;
    double qptny;
    double qptnz;
    double rprimd_ax; 
    double rprimd_ay; 
    double rprimd_az;
    double rprimd_bx; 
    double rprimd_by; 
    double rprimd_bz;
    double rprimd_cx; 
    double rprimd_cy; 
    double rprimd_cz;
    double ax_star; 
    double ay_star; 
    double az_star;
    double bx_star; 
    double by_star; 
    double bz_star;
    double cx_star; 
    double cy_star; 
    double cz_star;
    double ox,oy,oz;
    double theta1,phi1,g1;
    double theta2,phi2,g2;
    double stmbias;
    double tphysel;
    double tsmear;
    double znuclpsp;
    double zionpsp;
    double nonlocalE_byatom[NATOMS_MAX];
    double nonlocalE_byatom_IM[NATOMS_MAX];
    int pspso;
    int pspdat; 
    int pspcod; 
    int pspxc;
    int type;
    int lmn_size;
    int atom_type;
    int usewvl;
    int istwfk [KPOINTS_MAX];
    int istwfkv;
    int nband [BANDS_MAX];
    int nbandv;
    int npwarr [KPOINTS_MAX];
    int npwarrv;
    int so_psp [NATOMS_MAX];
    int symafm [MAX_SYM];
    int symrel [3][3][MAX_SYM];
    int typat [NATOMS_MAX];
    double kpt [3][KPOINTS_MAX];
    double occ;
    double cgt[2][WAVES_MAX][2];
    double tnons [3][MAX_SYM];
    double znucltypat[NATOMS_MAX];
    double wtk[KPOINTS_MAX];
    double wtkv;
    double residm;
    double x,y,z;
    double etotal,fermie;
    double Enonlocal_temp;
    double ga1, gb1, gc1;
    double gx1, gy1, gz1;
    double ga2, gb2, gc2;
    double gx2, gy2, gz2;
    double xred[3][NATOMS_MAX];
    double kptv;
    double php;
    int npw;
    int pw;
    int nspinor;
    int nband_temp;
    int kptno;
    int kx,ky,kz;
    int band;
    double normalization;
    double eigen,occ_temp,cg;
    double cellvolume,p1,p2;
    gsl_complex c1;
    gsl_complex c2;
    gsl_complex c1_star;
    gsl_complex c1c2;
    gsl_complex atom_phase;
    double atom_phase0;
    gsl_complex vg1g2;
    gsl_complex vg1g2_lm;
    gsl_complex vg1g2_temp;
    
    for(j=0;j<NATOMS_MAX;j++) {
            nonlocalE_byatom[j]=0.0;
            nonlocalE_byatom_IM[j]=0.0;
    }
    f2=fopen(filename,"rb+");
    if(f2==NULL) {
        printf("%s not found.\n",filename);
	exit(0);
    }
    /*   READ HEADER OF WFK FILE    */    
    j=0;
    fread(&j, sizeof(int), 1, f2);
    j=fread(codvsn, sizeof(char), 6, f2);
    fread(&headform, sizeof(int), 1, f2);
    fread(&fform, sizeof(int), 1, f2);
//    printf("%s %d %d\n",codvsn,headform,fform);
    fread(&j, sizeof(int), 1, f2);
    fread(&j, sizeof(int), 1, f2);
    fread(&bandtot, sizeof(int), 1, f2);
    /*    int date;*/
    fread(&date, sizeof(int), 1, f2);
//    printf("%d %d\n",bandtot,date);
    /*    int intxc;  */             
    fread(&intxc, sizeof(int), 1, f2);
    /*    int ixc;*/
    fread(&ixc, sizeof(int), 1, f2);
    /*    int natom;*/
    fread(&natom, sizeof(int), 1, f2);
    XSFIN->NIONS=natom;
    /*    int ngfftx;*/
    fread(&ngfftx, sizeof(int), 1, f2);
    /*    int ngffty;*/
    fread(&ngffty, sizeof(int), 1, f2);
    /*    int ngfftz;*/
    fread(&ngfftz, sizeof(int), 1, f2);
//    printf("ngfft = %d x %d x %d\n",ngfftx,ngffty,ngfftz);
    XSFIN->NGX=ngfftx+1;
    XSFIN->NGY=ngffty+1;
    XSFIN->NGZ=ngfftz+1; 
    /*    int nkpt;*/
    fread(&nkpt, sizeof(int), 1, f2);
    /*    int nspden;*/
    fread(&nspden, sizeof(int), 1, f2);
    /*    int nspinor;*/
    fread(&nspinor, sizeof(int), 1, f2);
    /*    int nsppol;*/
    fread(&nsppol, sizeof(int), 1, f2);
    /*    int nsym;*/
    fread(&nsym, sizeof(int), 1, f2);
    /*    int npsp;*/
    fread(&npsp, sizeof(int), 1, f2);
    /*    int ntypat;*/
    fread(&ntypat, sizeof(int), 1, f2);
    /*    int occopt;*/
    fread(&occopt, sizeof(int), 1, f2);
    /*    int pertcase;*/
    fread(&pertcase, sizeof(int), 1, f2);
    /*    int usepaw;*/
    fread(&usepaw, sizeof(int), 1, f2);
    /*    double ecut;*/
    fread(&ecut, sizeof(double), 1, f2);
    /*     double ecutdg;*/
    fread(&ecutdg, sizeof(double), 1, f2);
    /*    double ecutsm;*/
    fread(&ecutsm, sizeof(double), 1, f2);
    /*    double ecut_eff;*/
    fread(&ecut_eff, sizeof(double), 1, f2);
    /*    double qptnx;*/
    fread(&qptnx, sizeof(double), 1, f2);
    /*    double qptny;*/
    fread(&qptny, sizeof(double), 1, f2);
    /*    double qptnz;*/
    fread(&qptnz, sizeof(double), 1, f2);
    /*    double rprimd_ax;*/
    fread(&rprimd_ax, sizeof(double), 1, f2); 
    XSFIN->cella_x=rprimd_ax*0.52917720859;
    /*    double rprimd_ay;*/ 
    fread(&rprimd_ay, sizeof(double), 1, f2); 
    XSFIN->cella_y=rprimd_ay*0.52917720859;
    /*    double rprimd_az;*/
    fread(&rprimd_az, sizeof(double), 1, f2); 
    XSFIN->cella_z=rprimd_az*0.52917720859;
    /*    double rprimd_bx;*/ 
    fread(&rprimd_bx, sizeof(double), 1, f2);
    XSFIN->cellb_x=rprimd_bx*0.52917720859;
    /*    double rprimd_by;*/ 
    fread(&rprimd_by, sizeof(double), 1, f2);
    XSFIN->cellb_y=rprimd_by*0.52917720859;
    /*    double rprimd_bz;*/
    fread(&rprimd_bz, sizeof(double), 1, f2);
    XSFIN->cellb_z=rprimd_bz*0.52917720859;
    /*    double rprimd_cx;*/ 
    fread(&rprimd_cx, sizeof(double), 1, f2);
    XSFIN->cellc_x=rprimd_cx*0.52917720859;
    /*    double rprimd_cy;*/ 
    fread(&rprimd_cy, sizeof(double), 1, f2);
    XSFIN->cellc_y=rprimd_cy*0.52917720859;
    /*    double rprimd_cz;*/
    fread(&rprimd_cz, sizeof(double), 1, f2);
    XSFIN->cellc_z=rprimd_cz*0.52917720859;
    /*    double stmbias;*/
    cellvolume=(rprimd_ax*(rprimd_by*rprimd_cz-rprimd_bz*rprimd_cy)-rprimd_ay*(rprimd_bx*rprimd_cz-rprimd_bz*rprimd_cx)+rprimd_az*(rprimd_bx*rprimd_cy-rprimd_by*rprimd_cx));
//    printf(" cella = %lf %lf %lf \n",rprimd_ax,rprimd_ay,rprimd_az);
//    printf(" cellb = %lf %lf %lf \n",rprimd_bx,rprimd_by,rprimd_bz);
//    printf(" cellc = %lf %lf %lf \n\n",rprimd_cx,rprimd_cy,rprimd_cz);
//    printf(" cell volume = %lf\n\n",cellvolume);
    XSFIN->cellvolume=cellvolume;
    XSFIN->VoxelV=XSFIN->cellvolume/((XSFIN->NGX-1)*(XSFIN->NGY-1)*(XSFIN->NGZ-1))/pow(0.52917720859,3.0);
    ax_star=2*PI*(rprimd_by*rprimd_cz-rprimd_cy*rprimd_bz)/cellvolume;
    ay_star=-2*PI*(rprimd_bx*rprimd_cz-rprimd_cx*rprimd_bz)/cellvolume;
    az_star=2*PI*(rprimd_bx*rprimd_cy-rprimd_cx*rprimd_by)/cellvolume;
    bx_star=2*PI*(rprimd_cy*rprimd_az-rprimd_ay*rprimd_cz)/cellvolume;
    by_star=-2*PI*(rprimd_cx*rprimd_az-rprimd_ax*rprimd_cz)/cellvolume;
    bz_star=2*PI*(rprimd_cx*rprimd_ay-rprimd_ax*rprimd_cy)/cellvolume;
    cx_star=2*PI*(rprimd_ay*rprimd_bz-rprimd_by*rprimd_az)/cellvolume;
    cy_star=-2*PI*(rprimd_ax*rprimd_bz-rprimd_bx*rprimd_az)/cellvolume;
    cz_star=2*PI*(rprimd_ax*rprimd_by-rprimd_bx*rprimd_ay)/cellvolume;
//    printf(" cella* = %lf %lf %lf \n",ax_star,ay_star,az_star);
//    printf(" cellb* = %lf %lf %lf \n",bx_star,by_star,bz_star);
//    printf(" cellc* = %lf %lf %lf \n\n",cx_star,cy_star,cz_star);
//    printf(" a.a*=%lf, a.b*=%lf, a.c*=%lf\n",rprimd_ax*ax_star+rprimd_ay*ay_star+rprimd_az*az_star,rprimd_ax*bx_star+rprimd_ay*by_star+rprimd_az*bz_star,rprimd_ax*cx_star+rprimd_ay*cy_star+rprimd_az*cz_star);
//    printf(" b.a*=%lf, b.b*=%lf, b.c*=%lf\n",rprimd_bx*ax_star+rprimd_by*ay_star+rprimd_bz*az_star,rprimd_bx*bx_star+rprimd_by*by_star+rprimd_bz*bz_star,rprimd_bx*cx_star+rprimd_by*cy_star+rprimd_bz*cz_star);
//    printf(" c.a*=%lf, c.b*=%lf, c.c*=%lf\n\n",rprimd_cx*ax_star+rprimd_cy*ay_star+rprimd_cz*az_star,rprimd_cx*bx_star+rprimd_cy*by_star+rprimd_cz*bz_star,rprimd_cx*cx_star+rprimd_cy*cy_star+rprimd_cz*cz_star);

    fread(&stmbias, sizeof(double), 1, f2);
    /*    double tphysel;*/
    fread(&tphysel, sizeof(double), 1, f2);
    /*    double tsmear;*/
    fread(&tsmear, sizeof(double), 1, f2);
    /*    int usewvl;*/
    fread(&usewvl, sizeof(int), 1, f2);
//    printf("natoms = %d   ecut = %lf  tsmear = %lf  occopt = %d \n",natom,ecut,tsmear, occopt);   
    fread(&j, sizeof(int), 1, f2);
    fread(&j, sizeof(int), 1, f2);
    /*    int istwfk [KPOINTS_MAX]; */
//    printf("nkpt = %d    nsppol = %d   npsp = %d  ntypat = %d \n",nkpt,nsppol,npsp,ntypat);
    for(j=0;j<nkpt;j++) {
         fread(&istwfkv, sizeof(int), 1, f2);
    }
    /*    int nband [BANDS_MAX]; */
    for(j=0;j<(nkpt*nsppol);j++) {
         fread(&nbandv, sizeof(int), 1, f2);
    }
    /*    int npwarr [KPOINTS_MAX]; */
    for(j=0;j<(nkpt);j++) {
         fread(&npwarrv, sizeof(int), 1, f2);
    }
    /*    int so_psp [NATOMS_MAX]; */
    for(j=0;j<(npsp);j++) {
         fread(&so_psp[j], sizeof(int), 1, f2);
    }
    /*    int symafm [MAX_SYM]; */
    for(j=0;j<(nsym);j++) {
         fread(&symafm[j], sizeof(int), 1, f2);
    }
    /*    int symrel [3][3][MAX_SYM]; */
    for(j=0;j<(nsym);j++) {
         fread(&symrel[0][0][j], sizeof(int), 1, f2);
         fread(&symrel[1][0][j], sizeof(int), 1, f2);
         fread(&symrel[2][0][j], sizeof(int), 1, f2);
         fread(&symrel[0][1][j], sizeof(int), 1, f2);
         fread(&symrel[1][1][j], sizeof(int), 1, f2);
         fread(&symrel[2][1][j], sizeof(int), 1, f2);
         fread(&symrel[0][2][j], sizeof(int), 1, f2);
         fread(&symrel[1][2][j], sizeof(int), 1, f2);
         fread(&symrel[2][2][j], sizeof(int), 1, f2);
    }
    /*    int typat [NATOMS_MAX]; */
    for(j=0;j<(natom);j++) {
         fread(&typat[j], sizeof(int), 1, f2);
//	 printf("typeat = %d\n",typat[j]);         
    }    
    /*    double kpt [3][KPOINTS_MAX]; */
    for(j=0;j<(nkpt);j++) {
         fread(&kptv, sizeof(double), 1, f2);         
         fread(&kptv, sizeof(double), 1, f2);         
         fread(&kptv, sizeof(double), 1, f2);  
//         printf("kpoint %d:  %lf %lf %lf \n",j+1,kpt[0][j],kpt[1][j],kpt[2][j]);
    }    
    /*    double occ(TBANDS_MAX); */
    for(j=0;j<(bandtot);j++) {
         fread(&occ, sizeof(double), 1, f2);         
    }    
    /*    double tnons [3][MAX_SYM]; */
    for(j=0;j<(nsym);j++) {
         fread(&tnons[0][j], sizeof(double), 1, f2);         
         fread(&tnons[1][j], sizeof(double), 1, f2);         
         fread(&tnons[2][j], sizeof(double), 1, f2);         
    }    

    /*    double znucltypat[NATOMS_MAX];  */
    for(j=0;j<(ntypat);j++) {
         fread(&znucltypat[j], sizeof(double), 1, f2);
    }    

    for(j=0;j<natom;j++) {
          type=typat[j]-1;
//          printf("atom %d type=%d\n",j, type);
          XSFIN->atomicno[j]=(int)znucltypat[type];
    } 
    /*    double wtk[KPOINTS_MAX]; */
    for(j=0;j<(nkpt);j++) {
         fread(&wtkv, sizeof(double), 1, f2);         
    }    
    fread(&j, sizeof(int), 1, f2);
    for(k=0;k<(npsp);k++) {
          fread(&j, sizeof(int), 1, f2);
          fread(title, sizeof(char), 132, f2);
//          printf("     %s\n",title);
          fread(&znuclpsp, sizeof(double), 1, f2);         
//          printf("     %lf ",znuclpsp);
          fread(&zionpsp, sizeof(double), 1, f2);         
//          printf("%lf ",zionpsp);
          fread(&pspso, sizeof(int), 1, f2);
          fread(&pspdat, sizeof(int), 1, f2);
          fread(&pspcod, sizeof(int), 1, f2);
          fread(&pspxc, sizeof(int), 1, f2);
          fread(&lmn_size, sizeof(int), 1, f2);
//          printf("%d %d %d %d %d \n",pspso,pspdat,pspcod,pspxc,lmn_size);
          fread(&j, sizeof(int), 1, f2);
    }    
    if(usepaw==0) {
          fread(&j, sizeof(int), 1, f2);
          fread(&residm, sizeof(double), 1, f2);
//        printf("     residm = %lf\n",residm);
//          printf("Enter origin x y z: ");
	  for(k=0;k<natom;k++) {
	      fread(&x, sizeof(double), 1, f2);
	      fread(&y, sizeof(double), 1, f2);
	      fread(&z, sizeof(double), 1, f2);
//              printf("     Atom %d:  %lf %lf %lf \n",k+1,x,y,z);
              XSFIN->Xcart[k]=(x)*XSFIN->cella_x+(y)*XSFIN->cellb_x+(z)*XSFIN->cellc_x;
              XSFIN->Ycart[k]=(x)*XSFIN->cella_y+(y)*XSFIN->cellb_y+(z)*XSFIN->cellc_y;
              XSFIN->Zcart[k]=(x)*XSFIN->cella_z+(y)*XSFIN->cellb_z+(z)*XSFIN->cellc_z;
	  }    
	  fread(&etotal, sizeof(double), 1, f2);
	  fread(&fermie, sizeof(double), 1, f2);
//          printf("     Etotal = %lf    FermiE = %lf\n\n",etotal,fermie);
          fread(&j, sizeof(int), 1, f2);
    }
    else {
       printf("Yikes!  usepaw!=0.  We haven't written code for this case yet. \n");
    }
    /*   HEADER FINISHED   */    

    /*   READ DENSITY DATA   */
    fread(&j, sizeof(int), 1, f2);
    for(jz=0;jz<XSFIN->NGZ;jz++) {
       for(jy=0;jy<XSFIN->NGY;jy++) {
         for(jx=0;jx<XSFIN->NGX;jx++) {

           if((jx<XSFIN->NGX-1)&&(jy<XSFIN->NGY-1)&&(jz<XSFIN->NGZ-1)) {
               fread(&eigen,sizeof(double),1,f2);
               XSFIN->grid[jx][jy][jz]=eigen;
           }
         }
       }
    }
    fread(&j, sizeof(int), 1, f2);
    fclose(f2);    
    for(jz=0;jz<XSFIN->NGZ;jz++) {
       for(jy=0;jy<XSFIN->NGY;jy++) {
         for(jx=0;jx<XSFIN->NGX;jx++) {
            if(jx==XSFIN->NGX-1) XSFIN->grid[jx][jy][jz]=XSFIN->grid[0][jy][jz];
            if(jy==XSFIN->NGY-1) XSFIN->grid[jx][jy][jz]=XSFIN->grid[jx][0][jz];
            if(jz==XSFIN->NGZ-1) XSFIN->grid[jx][jy][jz]=XSFIN->grid[jx][jy][0];
            if((jx==XSFIN->NGX-1)&&(jy==XSFIN->NGY-1)) XSFIN->grid[jx][jy][jz]=XSFIN->grid[0][0][jz];
            if((jx==XSFIN->NGX-1)&&(jz==XSFIN->NGZ-1)) XSFIN->grid[jx][jy][jz]=XSFIN->grid[0][jy][0];
            if((jy==XSFIN->NGY-1)&&(jz==XSFIN->NGZ-1)) XSFIN->grid[jx][jy][jz]=XSFIN->grid[jx][0][0];

            if((jx==XSFIN->NGX-1)&&(jy==XSFIN->NGY-1)&&(jz==XSFIN->NGZ-1)) XSFIN->grid[jx][jy][jz]=XSFIN->grid[0][0][0];

         }
       }
    }
}


int main (int argc, char * argv[])
{
    int j;
    int nparams;
    FILE * f2;
    FILE * f4;
    double dV;
    int choice;
    char tim0 [100];
    char tim1 [100];
    char tim2 [100];
    char KDENfile1 [100];
    char KDENfile2 [100];
    char POTfile1 [100];
    char POTfile2 [100];
    char DENfile0 [100];
    char DENfile1 [100];
    char DENfile2 [100];
    char VHXCfile1 [100];
    char VHXCfile2 [100];
    char VHAfile1 [100];
    char VHAfile2 [100];
    char systemcmd [200];
    char WFKfilename [100];
    char outfilename[100];
    double res;
    double V1;
    double V2;
    double Ealpha_map1,Ealpha_map2,Ealpha_nomap1,Ealpha_nomap2;
    double P_entropy;
    double P_Ewald;
    double P_nonlocal;
    double P_Ealpha; 
    if (argc > 1) {
          strcpy(WFKfilename,argv[1]);
    }
    else {
       printf("Usage:  nonlocal_byatom filename \n");
       exit(0);
    }        

    read_calc_nonlocal(WFKfilename,&den0);
    strcat(WFKfilename,".xsf");
    f4=fopen(WFKfilename,"w");
    outputXSF(&den0,&den0,f4);
    fclose(f4);
}
