/*
  CPmap, part of the Fredrickson Group DFT Chemical Pressure Package
  Copyright (C) 2012-2018 by the Fredrickson Group at University of Wisconsin

  Last modified: July 31, 2018

  This program is free software: you can redistribute it or modify it under the
  terms of the GNU General Public License as published by the Free Software
  Foundation, either version 3 of the License, or any later version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANT-ABILITY or FITNESS
  FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

  You should have received a copy of the GNU General Public License along with
  this program. If not, see <http://www.gnu.org/licenses/>.

  To cite DFT-Chemical Pressure or for further reading, see the following papers:
  V.M. Berns, J. Engelkemier, Y. Guo, B.J. Kilduff, D.C. Fredrickson. JCTC. 2014, 10, 3380-3392.
  J. Engelkemier, V.M. Berns, D.C. Fredrickson. J. Chem. Theory Comput. 2013, 9, 3170-3180.
  D.C. Fredrickson. J. Am. Chem. Soc. 2012, 134, 5991-5999.
*/


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_erf.h>

#define AU2GPA 29421.912
#define L_MAX 7  /* higher than intended use */
#define NEQVOX 350
#define NIONMAX 200
#define NPOINTMAX 2000  /* radial density points */
#define ONETHIRD 0.333333333333333333333333333333
#define PI 3.14159265358979323846264338328
#define R_BOHR 0.52917721092  /* angstrom */
#define R_MAX 10.0  /* bohr */
#define STRMAX 100

/* GLOBAL VARIABLES AND STRUCTURES */

char abinitname[STRMAX], abinitout[STRMAX], cpoutname[STRMAX];
int dshi=1, dseq=2, dslo=3, kam=0, kbm=0, kcm=0, ngx=0, ngy=0, ngz=0, prof_nmax[NIONMAX];
int errcount=0, isradii=0, ixc=0, lmax=6, nspin=1, occopt=0, rescp=1, scheme=1, xcpot[2];
int mapcore=1, maphart=1, mapkin=1, maploc=1, mapkinoption=1, mapsym=1, mapxc=1, mapnonloc=1, pspoption=0;
int Ewald_map_done=0, E_alpha_option=0, E_Ewald_option=0, outputewaldalpha=0, vosctype[NIONMAX];
int CV_mode=0,nspinor;
int printbin=0, printen=0, printhmap=0, printvmap=0, standard=1;
double Ewald[NIONMAX], Ewald_vo[NIONMAX], Ewald_sc[NIONMAX], sc_elec[NIONMAX], vo_elec[NIONMAX];
double E_ewald[4], E_core[4], E_alpha[4], en_core[4], p_ealpha[4], p_ewald[4], p_hart[4], p_kin[4], p_loc[4], p_nonloc[4], p_xc[4];
double E_alpha_vo[NIONMAX], E_alpha_sc[NIONMAX], epsatm[NIONMAX];
double zion[NIONMAX], logint[NIONMAX], logslope[NIONMAX], tolerance=0.01, volhi=0.0, voleq=0.0, vollo=0.0;
double rhoprofile[NIONMAX][NPOINTMAX], rprofile[NIONMAX][NPOINTMAX];
double E_Ewald_by_atom[NIONMAX][4], E_alpha_by_atom[NIONMAX][4], E_NL[NIONMAX][5][2], E_NL_equiv[NIONMAX][2];
char profile_filenames[NIONMAX][STRMAX];
int psptypes[NIONMAX];
double sccounts[NIONMAX];
double NLremainder;
double *** xnewup, *** ynewup, *** znewup, *** xnewdn, *** ynewdn, *** znewdn;
FILE * cplog, * htmlfile, * javafile;
char logname[STRMAX];
char htmlname[STRMAX];
char javaname[STRMAX];
char templatefiles[NIONMAX][STRMAX];
int ntemplates;

/* Jonathan adding code */

// variables for cp_html
int factorial;
char bonds[NIONMAX][STRMAX];
char element1[STRMAX];
char element2[STRMAX];
char element3[STRMAX];
const char * colors[] = {"#ff0000", "#005aff", "#00bb00", "#cc00cc"};
int calcdate = 00000000;
int buttoncounter;
int count1;
char bond[STRMAX];
int l;
int m;

// variables for autocalibration
int it = 0; //iteration variable for 'for' loops
int autocali_stopper = 0;
int autocali_min_index = 100;
int autocali_num_site = 0;
double autocali_min_then = 0;
double autocali_min_now = 0;
double autocali_percent = 100;
double autocali_diff1 = 0;
double autocali_diff2 = 0;
double autocali_slope = 0;
double autocali_yint = 0;
double ps1[NIONMAX];
double ps2[NIONMAX];
double es0[NIONMAX];
double es1[NIONMAX];
double es2[NIONMAX];
double ps2_temp[NIONMAX];
double es2_temp[NIONMAX];


/* Jonathan done adding code */

struct CrystData {
  int nion, zatomic[NIONMAX];
  double cella_x, cella_y, cella_z, cellb_x, cellb_y, cellb_z, cellc_x, cellc_y, cellc_z;
  double corerad[NIONMAX], intCP[NIONMAX], intYlm[NIONMAX][L_MAX][2*L_MAX+1], voxcount[NIONMAX];
  double *** grid, entot, volcell, volvox, xcart[NIONMAX], ycart[NIONMAX], zcart[NIONMAX];
  double x[NIONMAX],y[NIONMAX],z[NIONMAX];
} denhi, deneq, denlo, etothi,etothi_temp, etoteq, etotlo, etotlo_temp,kdenhi, kdeneq, kdenlo, pothi, poteq,
  potlo, vhahi, vhaeq, vhalo, vhxchi, vhxceq, vhxclo, ewaldhi, ewaldeq, ewaldlo, alphahi, alphalo,
  core, cp, cp_Y, temp, vxc,
  gdenhi1, gdeneq1, gdenlo1, gdenhi2, gdeneq2, gdenlo2, gdenhi3, gdeneq3, gdenlo3,
  ldenhi, ldenhi2, ldeneq, ldeneq2, ldenlo, ldenlo2,
  denhi2, deneq2, denlo2, kdenhi2, kdeneq2, kdenlo2, pothi2, poteq2, potlo2, vhahi2,
  vhaeq2, vhalo2, vhxchi2, vhxceq2, vhxclo2, core2, temp2, vxc2, gdenhi4, gdeneq4,
  gdenlo4, gdenhi5, gdeneq5, gdenlo5, gdenhi6, gdeneq6, gdenlo6, lochi, loceq, loclo,
  nonlochi, nonloclo, nonlochi2, nonloclo2, denschi, densclo, densceq, bader,bader_temp,delta_pot_map,delta_pot_map2;

struct SymMap {
  int equiv[NIONMAX][NIONMAX], nequiv[NIONMAX], nsymel, symrel[1000][3][3];
  double tnons[1000][3];
} smap;

struct HirshMap {
  int atomid[NEQVOX], neighcount;
  double chg[NIONMAX], wj[NEQVOX], xcart[NEQVOX], ycart[NEQVOX], zcart[NEQVOX], *** hirsh_weight[NIONMAX];
} hmap;

struct ContactVol {
  int *** neighcount, *** neighcount2, *** neighkey[NEQVOX], *** ionmap[NEQVOX];
  double average[7][7][7][NIONMAX][NIONMAX], count[7][7][7][NIONMAX][NIONMAX];
  double total[7][7][7][NIONMAX][NIONMAX], *** swj, *** swjk, *** wj[NEQVOX];
  double *** wmax, *** wmax2;
} vmap;

/* SUPPORT FUNCTIONS */

int ElementName(int Z, char * name) {
  /* called by: PrintAverage, PrintResults, ReadProfile, SetBubbles, CalcEalpha */
  /* calls: none */
  const char * element[119] = {
    "&","H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar",
    "K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",
    "Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe",
    "Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu",
    "Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn",
    "Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr",
    "Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn","Uut","Fl","Uup","Lv","Uus","Uuo"};
  if (Z>=0 && Z<119) strncpy(name, element[Z], 4);
  else strncpy(name, "Un", 4);
  return 0;
}

int FinishLine(FILE * fptr) {
  /* called by: CoreCorrection, ReadProfile, ReadXSF, ReadNonlocAtom */
  /* calls: none */
  char check;
  int cont=0;
  while (cont==0) {
    check = fgetc(fptr);
    if (check==10 || check==EOF) cont = 1;
  }
  if (check==EOF) return 1;
  else return 0;
}

int ReadLine(FILE * fptr, char * newstr) {
  /* called by: CalcKineticTF, PressureContrib, ReadEwald */
  /* calls: none */
  char check;
  int cont=0, counter=0;
  while (cont==0) {
    check = getc(fptr);
    if (check==10) cont = 1;
    else if (check==EOF) return 1;
    *(newstr+counter) = check;
    counter++;
  }
  *(newstr+counter-1) = 0;
  return 0;
}

int Getkm(struct CrystData * gridin) {
  /* called by main */
  /* calls: none */
  double fa=0.0, fb=0.0, fc=0.0;
  fa = pow((gridin->cellb_y*gridin->cellc_z-gridin->cellb_z*gridin->cellc_y), 2)+
    pow((gridin->cellb_z*gridin->cellc_x-gridin->cellb_x*gridin->cellc_z), 2)+
    pow((gridin->cellb_x*gridin->cellc_y-gridin->cellb_y*gridin->cellc_x), 2);
  fb = pow((gridin->cellc_y*gridin->cella_z-gridin->cellc_z*gridin->cella_y), 2)+
    pow((gridin->cellc_z*gridin->cella_x-gridin->cellc_x*gridin->cella_z), 2)+
    pow((gridin->cellc_x*gridin->cella_y-gridin->cellc_y*gridin->cella_x), 2);
  fc = pow((gridin->cella_y*gridin->cellb_z-gridin->cella_z*gridin->cellb_y), 2)+
    pow((gridin->cella_z*gridin->cellb_x-gridin->cella_x*gridin->cellb_z), 2)+
    pow((gridin->cella_x*gridin->cellb_y-gridin->cella_y*gridin->cellb_x), 2);
  kam = (int)ceil(R_MAX/(gridin->volcell/sqrt(fa)));
  kbm = (int)ceil(R_MAX/(gridin->volcell/sqrt(fb)));
  kcm = (int)ceil(R_MAX/(gridin->volcell/sqrt(fc)));
  printf("  Using supercell range %d %d %d\n", kam, kbm, kcm);
  fprintf(cplog, "Supercell range: %d %d %d\n", kam, kbm, kcm);
  if(kam>3 || kbm>3 || kcm>3) {
    printf("\n  BAD NEWS: Unit cell too small or distorted!\n");
    fprintf(cplog, "\nTerminated because an appropriate supercell could not be found\n");
    fprintf(cplog, "Suggestion: check your files or decrease R_MAX and recompile\n");
    errcount++;
    return 1;
  }
  return 0;
}

double Getwj(int atom, double dist) {
  /* called by: CoordSearch, CoreUnwarp */
  /* calls: none */
  int start=0;
  double xmax=0.0, xmin=0.0;
  start = (int)ceil((log(dist+1.0e-6)-logint[atom])/logslope[atom]);
  if (start<0) start = 0;
  while (start<prof_nmax[atom] && rprofile[atom][start]<dist) start++;
  while (rprofile[atom][start]>dist && start>0) start--;
  if (start==0) {
    if (rhoprofile[atom][0]<=0.0) {
      printf("\n  BAD NEWS: Electron density profile for atom #%d is abnormal!\n", atom+1);
      fprintf(cplog, "\nTerminated because density=%f at dist=%f for atom #%d is <= 0\n", rhoprofile[atom][0], dist, atom+1);
      fprintf(cplog, "Suggestion: check your profiles for non-numerical data\n");
      errcount++;
      return -1000.0;
    } else return rhoprofile[atom][0];
  }
  else if (start>=prof_nmax[atom]) {
    printf("\n  CAUTION: The atomic density profiles are too short! Continuing anyway...\n");
    fprintf(cplog, "WARNING: atomic radial density profile for atom #%d should be longer than %.6e bohr\n", atom+1, dist);
    errcount++;
    return 0.0;  /* assume no density at long distance */
  } else {  /* interpolation between closest two points in radial profile */
    xmin = rhoprofile[atom][start];
    xmax = rhoprofile[atom][start+1];
    return (xmin*(rprofile[atom][start+1]-dist)+xmax*(dist-rprofile[atom][start]))/
      (rprofile[atom][start+1]-rprofile[atom][start]);
  }
}

int Cart2Sph(double x, double y, double z, double * r, double * theta, double * phi) {
  /* called by: AverageContact */
  /* calls: none */
  *r = sqrt(x*x+y*y+z*z);
  if (*r>0.0) {
    *theta = acos(z / *r);
    if (fabs(sin(*theta))>0.0) {
      *phi = acos(x/(*r * sin(*theta)));
      if(x/(*r * sin(*theta))>1.0) *phi = acos(1.0);
      if(x/(*r * sin(*theta))<-1.0) *phi = acos(-1.0);
      if(y<0.0) *phi = -*phi;
    } else *phi = 0;
  } else {
    *phi = 0.0;
    *theta = 0.0;
  }
  return 0;
}


/* VOXEL GRID FUNCTIONS */

int FixEdges(struct CrystData * gridin) {
  /* called by: AverageContact, Bin2XSF, CalcVxc, CalcVxc1, CoreUnwarp, OutputXSF, SymmetrizeGrid, OutputWeight, CoreCorrection */
  /* calls: none */
  int jx=0, jy=0, jz=0;
  for (jy=0; jy<=ngy; jy++) {
    for (jx=0; jx<=ngx; jx++) {
      gridin->grid[jx][jy][ngz] = gridin->grid[jx][jy][0];
    }
  }
  for (jz=0; jz<=ngz; jz++) {
    for (jx=0; jx<=ngx; jx++) {
      gridin->grid[jx][ngy][jz] = gridin->grid[jx][0][jz];
    }
  }
  for (jz=0; jz<=ngz; jz++) {
    for (jy=0; jy<=ngy; jy++) {
      gridin->grid[ngx][jy][jz] = gridin->grid[0][jy][jz];
    }
  }
  for (jz=0; jz<=ngz; jz++) {
    gridin->grid[ngx][ngy][jz] = gridin->grid[0][0][jz];
  }
  for (jy=0; jy<=ngy; jy++) {
    gridin->grid[ngx][jy][ngz] = gridin->grid[0][jy][0];
  }
  for (jx=0; jx<=ngx; jx++) {
    gridin->grid[jx][ngy][ngz] = gridin->grid[jx][0][0];
  }
  gridin->grid[ngx][ngy][ngz] = gridin->grid[0][0][0];
  return 0;
}

double IntegrateGrid(struct CrystData * gridin) {
  /* called by: CalcCP, GridStats, MapEntot, MapEwald */
  /* calls: none */
  int jx=0, jy=0, jz=0;
  double sum=0.0;
  for (jz=0; jz<ngz; jz++) {
    for (jy=0; jy<ngy; jy++) {
      for (jx=0; jx<ngx; jx++) {
        sum += gridin->grid[jx][jy][jz];
      }
    }
  }
  return (sum*gridin->volvox);
}

int ShiftGrid(struct CrystData * gridin, double num, struct CrystData * gridout) {
  /* called by: CalcCP, CoreCorrection, GridStats */
  /* calls: none */
  int jx=0, jy=0, jz=0;
  gridout->volvox = gridin->volvox;
  for (jz=0; jz<=ngz; jz++) {
    for (jy=0; jy<=ngy; jy++) {
      for (jx=0; jx<=ngx; jx++) {
        gridout->grid[jx][jy][jz] = gridin->grid[jx][jy][jz]+num;
      }
    }
  }
  return 0;
}

int ScaleGrid(struct CrystData * gridin, double factor, struct CrystData * gridout) {
  /* called by: CalcCP, MapEntot */
  /* calls: none */
  int jx=0, jy=0, jz=0;
  gridout->volvox = gridin->volvox;
  for (jz=0; jz<=ngz; jz++) {
    for (jy=0; jy<=ngy; jy++) {
      for (jx=0; jx<=ngx; jx++) {
        gridout->grid[jx][jy][jz] = gridin->grid[jx][jy][jz]*factor;
      }
    }
  }
  return 0;
}

int PowerGrid(struct CrystData * gridin, int expon, struct CrystData * gridout) {
  /* called by: GridStats */
  /* calls: none */
  int jx=0, jy=0, jz=0;
  gridout->volvox = gridin->volvox;
  for (jz=0; jz<=ngz; jz++) {
    for (jy=0; jy<=ngy; jy++) {
      for (jx=0; jx<=ngx; jx++) {
        gridout->grid[jx][jy][jz] = pow(gridin->grid[jx][jy][jz], expon);
      }
    }
  }
  return 0;
}

int AddGrid(struct CrystData * gridin, struct CrystData * gridin2, struct CrystData * gridout) {
  /* called by: MapEntot */
  /* calls: none */
  int jx=0, jy=0, jz=0;
  gridout->volvox = gridin->volvox;
  for (jz=0; jz<=ngz; jz++) {
    for (jy=0; jy<=ngy; jy++) {
      for (jx=0; jx<=ngx; jx++) {
        gridout->grid[jx][jy][jz] = gridin->grid[jx][jy][jz]+gridin2->grid[jx][jy][jz];
      }
    }
  }
  return 0;
}

int SubtractGrid(struct CrystData * gridin, struct CrystData * gridin2, struct CrystData * gridout) {
  /* called by: CalcCP, MapEntot */
  /* calls: none */
  int jx=0, jy=0, jz=0;
  gridout->volvox = gridin->volvox;
  for (jz=0; jz<=ngz; jz++) {
    for (jy=0; jy<=ngy; jy++) {
      for (jx=0; jx<=ngx; jx++) {
        gridout->grid[jx][jy][jz] = gridin->grid[jx][jy][jz]-gridin2->grid[jx][jy][jz];
      }
    }
  }
  return 0;
}

int MultiplyGrid(struct CrystData * gridin, struct CrystData * gridin2, struct CrystData * gridout) {
  /* called by: MapEntot */
  /* calls: none */
  int jx=0, jy=0, jz=0;
  gridout->volvox = gridin->volvox;
  for (jz=0; jz<=ngz; jz++) {
    for (jy=0; jy<=ngy; jy++) {
      for (jx=0; jx<=ngx; jx++) {
        gridout->grid[jx][jy][jz] = gridin->grid[jx][jy][jz]*gridin2->grid[jx][jy][jz];
      }
    }
  }
  return 0;
}

double GridGrad(int jx, int jy, int jz, struct CrystData * gridin, double *gradx, double *grady, double *gradz) {
  /* called by: MapEntot */
  /* calls: none */
  /* Added by DCF, 9/(28-29)/18 */
  double astepx,bstepx,cstepx;
  double astepy,bstepy,cstepy;
  double astepz,bstepz,cstepz;
  double df_a, df_b, df_c;
  int jx_up,jy_up,jz_up;
  int jx_down,jy_down,jz_down;
  gsl_matrix *step = gsl_matrix_alloc(3,3);
  gsl_vector *df = gsl_vector_alloc(3);
  jx_up = jx+1;
  if(jx_up == ngx) jx_up = 0;
  jy_up = jy+1;
  if(jy_up == ngy) jy_up = 0;
  jz_up = jz+1;
  if(jz_up == ngz) jz_up = 0;
  jx_down = jx-1;
  if(jx_down == -1) jx_down = ngx-1;
  jy_down = jy-1;
  if(jy_down == -1) jy_down = ngy-1;
  jz_down = jz-1;
  if(jz_down == -1) jz_down = ngz-1;
  astepx = 2.0*gridin->cella_x/(1.0*ngx);
  astepy = 2.0*gridin->cella_y/(1.0*ngx);
  astepz = 2.0*gridin->cella_z/(1.0*ngx);
  bstepx = 2.0*gridin->cellb_x/(1.0*ngy);
  bstepy = 2.0*gridin->cellb_y/(1.0*ngy);
  bstepz = 2.0*gridin->cellb_z/(1.0*ngy);
  cstepx = 2.0*gridin->cellc_x/(1.0*ngz);
  cstepy = 2.0*gridin->cellc_y/(1.0*ngz);
  cstepz = 2.0*gridin->cellc_z/(1.0*ngz);
  df_a = (gridin->grid[jx_up][jy][jz]-gridin->grid[jx_down][jy][jz]);
  df_b = (gridin->grid[jx][jy_up][jz]-gridin->grid[jx][jy_down][jz]);
  df_c = (gridin->grid[jx][jy][jz_up]-gridin->grid[jx][jy][jz_down]);
  gsl_vector_set(df, 0, df_a);
  gsl_vector_set(df, 1, df_b);
  gsl_vector_set(df, 2, df_c);
  gsl_matrix_set(step,0,0, astepx);
  gsl_matrix_set(step,0,1, astepy);
  gsl_matrix_set(step,0,2, astepz);
  gsl_matrix_set(step,1,0, bstepx);
  gsl_matrix_set(step,1,1, bstepy);
  gsl_matrix_set(step,1,2, bstepz);
  gsl_matrix_set(step,2,0, cstepx);
  gsl_matrix_set(step,2,1, cstepy);
  gsl_matrix_set(step,2,2, cstepz);
  gsl_linalg_HH_svx(step,df);
  *gradx = gsl_vector_get(df,0);
  *grady = gsl_vector_get(df,1);
  *gradz = gsl_vector_get(df,2);
  gsl_vector_free(df);
  gsl_matrix_free(step);
  return pow(*gradx*(*gradx)+*grady*(*grady)+*gradz*(*gradz),0.5);
}




int CopyStruct(struct CrystData * gridin, struct CrystData * gridout) {
  /* called by: CalcCP, MapEntot, MapEwald */
  /* calls: none */
  int i=0;
  gridout->cella_x = gridin->cella_x;
  gridout->cella_y = gridin->cella_y;
  gridout->cella_z = gridin->cella_z;
  gridout->cellb_x = gridin->cellb_x;
  gridout->cellb_y = gridin->cellb_y;
  gridout->cellb_z = gridin->cellb_z;
  gridout->cellc_x = gridin->cellc_x;
  gridout->cellc_y = gridin->cellc_y;
  gridout->cellc_z = gridin->cellc_z;
  gridout->nion = gridin->nion;
  gridout->volcell = gridin->volcell;
  gridout->volvox = gridin->volvox;
  for (i=0; i<gridin->nion; i++) {
    gridout->corerad[i] = gridin->corerad[i];
    gridout->xcart[i] = gridin->xcart[i];
    gridout->ycart[i] = gridin->ycart[i];
    gridout->zcart[i] = gridin->zcart[i];
    gridout->zatomic[i] = gridin->zatomic[i];
  }
  return 0;
}

int GridStats(struct CrystData * gridin, const char * str) {
  /* called by: main */
  /* calls: IntegrateGrid, PowerGrid, ShiftGrid */
  int jx=0, jy=0, jz=0, nvox=0;
  double min=1.0e7, max=-1.0e7, kurt=0.0, mean=0.0, sdev=0.0, skew=0.0, sum=0.0, var=0.0;
  for (jz=0; jz<ngz; jz++) {
    for (jy=0; jy<ngy; jy++) {
      for (jx=0; jx<ngx; jx++) {
        if (gridin->grid[jx][jy][jz]>max) max=gridin->grid[jx][jy][jz];
        if (gridin->grid[jx][jy][jz]<min) min=gridin->grid[jx][jy][jz];
      }
    }
  }
  nvox = ngx*ngy*ngz;
  sum = IntegrateGrid(gridin)/gridin->volvox;
  mean = sum/(double)nvox;
  ShiftGrid(gridin, -mean, &temp);
  PowerGrid(&temp, 2, &temp);
  sum = IntegrateGrid(&temp)/gridin->volvox;
  var = sum/(double)(nvox-1);
  sdev = sqrt(var);
  ShiftGrid(gridin, -mean, &temp);
  PowerGrid(&temp, 3, &temp);
  sum = IntegrateGrid(&temp)/gridin->volvox;
  skew = sum/(double)(nvox-1)/pow(sdev, 3);
  ShiftGrid(gridin, -mean, &temp);
  PowerGrid(&temp, 4, &temp);
  sum = IntegrateGrid(&temp)/gridin->volvox;
  kurt = sum/(double)(nvox-1)/pow(sdev, 4)-3;
  sum = IntegrateGrid(gridin)/gridin->volvox;
  fprintf(cplog, "%s:\n    voxel volume: %12.6f\n max voxel value: %12.6f\n min voxel value: %12.6f\n       voxel sum: %12.6f\n            mean: %12.6f\n        variance: %20.16f\n        skewness: %12.6f\n        kurtosis: %12.6f\n\n",
    str, gridin->volvox, max, min, sum, mean, var, skew, kurt);
  return 0;
}

/* SYMMETRY FUNCTIONS */

int SymAtoms(struct CrystData * gridin, struct SymMap * map) {
  /* Calls: none */
  /* called by: main */
  int i=0, j=0, k=0;
  double ax=gridin->cella_x, ay=gridin->cella_y, az=gridin->cella_z;
  double bx=gridin->cellb_x, by=gridin->cellb_y, bz=gridin->cellb_z;
  double cx=gridin->cellc_x, cy=gridin->cellc_y, cz=gridin->cellc_z;
  double det=0.0, xf=0.0, yf=0.0, zf=0.0, xf2=0.0, yf2=0.0, zf2=0.0;
  double coord_eqx[NIONMAX][1000], coord_eqy[NIONMAX][1000], coord_eqz[NIONMAX][1000];
  for (i=0; i<NIONMAX; i++) map->nequiv[i] = 0;
  for (i=0; i<gridin->nion; i++) {
    for (j=0; j<gridin->nion; j++) {
      map->equiv[i][j] = 0;
    }
  }
  det = 1.0/(ax*by*cz-ax*cy*bz-bx*ay*cz+bx*cy*az+cx*ay*bz-cx*by*az);
  for (k=0; k<map->nsymel; k++) {
    for (i=0; i<gridin->nion; i++) {
      xf = det*((by*cz-cy*bz)*gridin->xcart[i]+(cx*bz-bx*cz)*gridin->ycart[i]+(bx*cy-cx*by)*gridin->zcart[i]);
      yf = det*((cy*az-ay*cz)*gridin->xcart[i]+(ax*cz-cx*az)*gridin->ycart[i]+(cx*ay-ax*cy)*gridin->zcart[i]);
      zf = det*((ay*bz-by*az)*gridin->xcart[i]+(bx*az-ax*bz)*gridin->ycart[i]+(ax*by-bx*ay)*gridin->zcart[i]);
      xf2 = xf*map->symrel[k][0][0]+yf*map->symrel[k][0][1]+zf*map->symrel[k][0][2]+map->tnons[k][0];
      yf2 = xf*map->symrel[k][1][0]+yf*map->symrel[k][1][1]+zf*map->symrel[k][1][2]+map->tnons[k][1];
      zf2 = xf*map->symrel[k][2][0]+yf*map->symrel[k][2][1]+zf*map->symrel[k][2][2]+map->tnons[k][2];
      while (xf2<0.0)  xf2 += 1.0;
      while (xf2>=1.0) xf2 -= 1.0;
      while (yf2<0.0)  yf2 += 1.0;
      while (yf2>=1.0) yf2 -= 1.0;
      while (zf2<0.0)  zf2 += 1.0;
      while (zf2>=1.0) zf2 -= 1.0;
      /* Fixed treatment of coordinates that need to be shifted by negative integers into home cell.   -DCF 6/22/19 */
      coord_eqx[i][k] = xf2;
      coord_eqy[i][k] = yf2;
      coord_eqz[i][k] = zf2;
    }
  }
  for (i=0; i<gridin->nion; i++) {
    for (k=0; k<map->nsymel; k++) {
      for (j=0; j<gridin->nion; j++) {
        xf = det*((by*cz-cy*bz)*gridin->xcart[j]+(cx*bz-bx*cz)*gridin->ycart[j]+(bx*cy-cx*by)*gridin->zcart[j]);
        yf = det*((cy*az-ay*cz)*gridin->xcart[j]+(ax*cz-cx*az)*gridin->ycart[j]+(cx*ay-ax*cy)*gridin->zcart[j]);
        zf = det*((ay*bz-by*az)*gridin->xcart[j]+(bx*az-ax*bz)*gridin->ycart[j]+(ax*by-bx*ay)*gridin->zcart[j]);
        while (xf<0.0)  xf += 1.0;
        while (xf>=1.0) xf -= 1.0;
        while (yf<0.0)  yf += 1.0;
        while (yf>=1.0) yf -= 1.0;
        while (zf<0.0)  zf += 1.0;
        while (zf>=1.0) zf -= 1.0;
        if (fabs(coord_eqx[i][k]-xf)<0.01 && fabs(coord_eqy[i][k]-yf)<0.01 && fabs(coord_eqz[i][k]-zf)<0.01) {
          map->equiv[i][j] = 1;
        }
      }
    }
  }
  for (i=0; i<gridin->nion; i++) {
    for (j=0; j<gridin->nion; j++) {
      if (map->equiv[i][j]==1) map->nequiv[i]++;
    }
  }
  return 0;
}

int CheckSym(struct SymMap * map) {
  /* called by: SymmetrizeGrid */
  /* calls: none */
  char ngn[4];
  int i=0, j=0, halfflag[3], quartflag[3];
  for (i=0; i<map->nsymel; i++) {
    for (j=0; j<3; j++) {
      if (fabs(map->tnons[i][j]-0.5)<1.0e-10) halfflag[j] = 1;
      else if (fabs(map->tnons[i][j]-0.25)<1.0e-10) quartflag[j] = 1;
    }
  }
  for (j=0; j<3; j++) {
    if (j==0) strncpy(ngn, "ngx", 3);
    else if (j==1) strncpy(ngn, "ngy", 3);
    else if (j==2) strncpy(ngn, "ngz", 3);
    if (halfflag[j]==1 && quartflag[j]==1) fprintf(cplog, "%s FFT index must be divisible by 4\n", ngn);
    else if (quartflag[j]==1) fprintf(cplog, "%s FFT index must be divisible by 4 if it is divisible by 2\n", ngn);
    else if (halfflag[j]==1) fprintf(cplog, "%s FFT index must be divisible by 2\n", ngn);
    else fprintf(cplog, "%s FFT index has no symmetry restrictions\n", ngn);
  }
  return 0;
}

int SymmetrizeGrid(struct CrystData * gridin, struct SymMap * map) {
  /* called by: CalcCP */
  /* calls: CheckSym, FixEdges, MapNonloc */
  int flag=0, i=0, jx=0, jy=0, jz=0, jx2=0, jy2=0, jz2=0;
  double xf=0.0, yf=0.0, zf=0.0, xf2=0.0, yf2=0.0, zf2=0.0;
  double voxspacing_xf=0.0, voxspacing_yf=0.0, voxspacing_zf=0.0;
  voxspacing_xf = 1.0/(double)ngx;
  voxspacing_yf = 1.0/(double)ngy;
  voxspacing_zf = 1.0/(double)ngz;
  for (jz=0; jz<=ngz; jz++) {
    for (jy=0; jy<=ngy; jy++) {
      for (jx=0; jx<=ngx; jx++) {
        core.grid[jx][jy][jz] = 0.0;  /* temp voxel counter */
      }
    }
  }
  for (i=0; i<map->nsymel; i++) {
    for (jz=0; jz<ngz; jz++) {
      zf = (double)jz/(double)ngz;
      for (jy=0; jy<ngy; jy++) {
        yf = (double)jy/(double)ngy;
        for (jx=0; jx<ngx; jx++) {
          xf = (double)jx/(double)ngx;
          xf2 = xf*map->symrel[i][0][0]+yf*map->symrel[i][0][1]+zf*map->symrel[i][0][2]+map->tnons[i][0];
          yf2 = xf*map->symrel[i][1][0]+yf*map->symrel[i][1][1]+zf*map->symrel[i][1][2]+map->tnons[i][1];
          zf2 = xf*map->symrel[i][2][0]+yf*map->symrel[i][2][1]+zf*map->symrel[i][2][2]+map->tnons[i][2];
          jx2 = (int)(floor(xf2/voxspacing_xf+0.5));
          jy2 = (int)(floor(yf2/voxspacing_yf+0.5));
          jz2 = (int)(floor(zf2/voxspacing_zf+0.5));
          while (jx2<0) { jx2 += ngx; }
          while (jx2>=ngx) { jx2 -= ngx; }
          while (jy2<0) { jy2 += ngy; }
          while (jy2>=ngy) { jy2 -= ngy; }
          while (jz2<0) { jz2 += ngz; }
          while (jz2>=ngz) { jz2 -= ngz; }
          core.grid[jx][jy][jz] += gridin->grid[jx2][jy2][jz2];
          if (fabs(xf2/voxspacing_xf-0.5)<0.00001) flag = 1;
          else if (fabs(yf2/voxspacing_yf-0.5)<0.00001) flag = 1;
          else if (fabs(zf2/voxspacing_zf-0.5)<0.00001) flag = 1;
          if (flag==1) {
            printf("\n  BAD NEWS: FFT grid is incommensurate with symmetry operation %d!\n", i+1);
            CheckSym(&smap);
            fprintf(cplog, "\nTerminated because FFT grid is not compatible with symmetry operation %d\n", i+1);
            fprintf(cplog, "Suggestion: re-run Abinit with another ngfft based on the recommendations above\n");
            fprintf(cplog, "Or re-run CP with symmetry restoration turned off (not recommended)\n");
            errcount++;
            return 2;
          }
        }
      }
    }
  }
  for (jz=0; jz<ngz; jz++) {
    for (jy=0; jy<ngy; jy++) {
      for (jx=0; jx<ngx; jx++) {
        gridin->grid[jx][jy][jz] = core.grid[jx][jy][jz]/(double)map->nsymel;
      }
    }
  }
  FixEdges(gridin);
  return 0;
}


/* I/O FUNCTIONS */

int Bin2XSF(char binname[STRMAX], struct CrystData * gridin, struct CrystData * gridin2) {
  /* called by: Den2XSF */
  /* calls: FixEdges */
  char codvsn[10], title[150];
  int bandtot=0, date=0, fform=0, headform=0, lmn_size=0, occopt_this=0;
  int intxc=0, istwfkv=0, ixc_this=0, i=0, jx=0, jy=0, jz=0, k=0;
  int natom=0, nbandv=0, ngfftx=0, ngffty=0, ngfftz=0, nkpt=0, npsp=0;
  int npwarrv=0, nspden=0, nsppol=0, nsym=0, ntypat=0;
  int pertcase=0, pspcod=0, pspdat=0, pspso=0, pspxc=0, type=0, usepaw=0, usewvl=0;
  int so_psp[NIONMAX], symafm[1000], symrel[3][3][1000], typat[NIONMAX];
  double cellvol=0.0, ecut=0.0, ecutdg=0.0, ecut_eff=0.0, ecutsm=0.0, eigen=0.0, etotal=0.0;
  double efermi=0.0, kptv=0.0, occ=0.0, qptnx=0.0, qptny=0.0, qptnz=0.0, residm=0.0;
  double rprimd_ax=0.0, rprimd_ay=0.0, rprimd_az=0.0, rprimd_bx=0.0;
  double rprimd_by=0.0, rprimd_bz=0.0, rprimd_cx=0.0, rprimd_cy=0.0, rprimd_cz=0.0;
  double stmbias=0.0, tphysel=0.0, tsmear=0.0, wtkv=0.0, x=0.0, y=0.0, z=0.0, zionpsp=0.0, znuclpsp=0.0;
  double tnons[3][1000], znucltypat[NIONMAX];
  FILE * fptr;
  fptr = fopen(binname, "rb+");
  if (fptr==NULL) {
    printf("\n  BAD NEWS: File %s not found!\n", binname);
    fprintf(cplog, "\nTerminated because binary file %s not found\n", binname);
    fprintf(cplog, "Suggestion: check your Abinit files or input options\n");
    errcount++;
    return 1;
  }
  fread(&i, sizeof(int), 1, fptr);
  i = fread(codvsn, sizeof(char), 6, fptr);
  fread(&headform, sizeof(int), 1, fptr);
  fread(&fform, sizeof(int), 1, fptr);
  fread(&i, sizeof(int), 1, fptr);
  fread(&i, sizeof(int), 1, fptr);
  fread(&bandtot, sizeof(int), 1, fptr);
  fread(&date, sizeof(int), 1, fptr);
  fprintf(cplog, "Start date of %s: %d\n", binname, date);
  calcdate = date;
  fread(&intxc, sizeof(int), 1, fptr);
  fread(&ixc_this, sizeof(int), 1, fptr);
  ixc = ixc_this;  /* global */
  fread(&natom, sizeof(int), 1, fptr);
  gridin->nion = natom;
  fread(&ngfftx, sizeof(int), 1, fptr);
  fread(&ngffty, sizeof(int), 1, fptr);
  fread(&ngfftz, sizeof(int), 1, fptr);
  if (ngx!=0 || ngy!=0 || ngz!=0) {
    if (ngx!=ngfftx || ngy!=ngffty || ngz!=ngfftz) {
    printf("\n  BAD NEWS: Different numbers of voxels between datasets!\n");
    fprintf(cplog, "\nTerminated because not all FFT grids are equal\n");
    fprintf(cplog, "\nPrevious grid: %d %d %d, New grid for %s: %d %d %d\n",
      binname, ngx, ngy, ngz, ngfftx, ngffty, ngfftz);
    fprintf(cplog, "Suggestion: specify ngfft in Abinit input such that boxcut is ~2\n");
    errcount++;
    return 2;
    }
  }
  ngx = ngfftx;  /* global */
  ngy = ngffty;  /* global */
  ngz = ngfftz;  /* global */
  fread(&nkpt, sizeof(int), 1, fptr);
  fread(&nspden, sizeof(int), 1, fptr);
  fread(&nspinor, sizeof(int), 1, fptr);
  fread(&nsppol, sizeof(int), 1, fptr);
  nspin = nsppol;  /* global */
  if ((nspin==2 && nspden!=2)) {
    printf("\n  BAD NEWS: Spin-orbit coupling, non-scalar magnetism, or non-collinear magnetism detected!\n");
    fprintf(cplog, "\nTerminated because nspinor=1, or nsppol=2 and nspden is not 2\n");
    fprintf(cplog, "Suggestion: Check your files or re-run Abinit with different input options\n");
    errcount++;
    return 3;
  }
  fread(&nsym, sizeof(int), 1, fptr);
  smap.nsymel = nsym;
  fread(&npsp, sizeof(int), 1, fptr);
  fread(&ntypat, sizeof(int), 1, fptr);
  fread(&occopt_this, sizeof(int), 1, fptr);
  occopt = occopt_this;  /* global */
  fread(&pertcase, sizeof(int), 1, fptr);
  fread(&usepaw, sizeof(int), 1, fptr);
  fread(&ecut, sizeof(double), 1, fptr);
  fread(&ecutdg, sizeof(double), 1, fptr);
  fread(&ecutsm, sizeof(double), 1, fptr);
  fread(&ecut_eff, sizeof(double), 1, fptr);
  fread(&qptnx, sizeof(double), 1, fptr);
  fread(&qptny, sizeof(double), 1, fptr);
  fread(&qptnz, sizeof(double), 1, fptr);
  fread(&rprimd_ax, sizeof(double), 1, fptr);
  gridin->cella_x = rprimd_ax;
  fread(&rprimd_ay, sizeof(double), 1, fptr);
  gridin->cella_y = rprimd_ay;
  fread(&rprimd_az, sizeof(double), 1, fptr);
  gridin->cella_z = rprimd_az;
  fread(&rprimd_bx, sizeof(double), 1, fptr);
  gridin->cellb_x = rprimd_bx;
  fread(&rprimd_by, sizeof(double), 1, fptr);
  gridin->cellb_y = rprimd_by;
  fread(&rprimd_bz, sizeof(double), 1, fptr);
  gridin->cellb_z = rprimd_bz;
  fread(&rprimd_cx, sizeof(double), 1, fptr);
  gridin->cellc_x = rprimd_cx;
  fread(&rprimd_cy, sizeof(double), 1, fptr);
  gridin->cellc_y = rprimd_cy;
  fread(&rprimd_cz, sizeof(double), 1, fptr);
  gridin->cellc_z = rprimd_cz;
  cellvol = (rprimd_ax*(rprimd_by*rprimd_cz-rprimd_bz*rprimd_cy)-rprimd_ay*
    (rprimd_bx*rprimd_cz-rprimd_bz*rprimd_cx)+rprimd_az*(rprimd_bx*rprimd_cy-rprimd_by*rprimd_cx));
  gridin->volcell = cellvol;
  gridin->volvox = cellvol/(ngx*ngy*ngz);
  fread(&stmbias, sizeof(double), 1, fptr);
  fread(&tphysel, sizeof(double), 1, fptr);
  fread(&tsmear, sizeof(double), 1, fptr);
  fread(&usewvl, sizeof(int), 1, fptr);
  fread(&i, sizeof(int), 1, fptr);
  fread(&i, sizeof(int), 1, fptr);
  for (i=0; i<nkpt; i++) fread(&istwfkv, sizeof(int), 1, fptr);
  for (i=0; i<(nkpt*nsppol); i++) fread(&nbandv, sizeof(int), 1, fptr);
  for (i=0; i<(nkpt); i++) fread(&npwarrv, sizeof(int), 1, fptr);
  for (i=0; i<(npsp); i++) fread(&so_psp[i], sizeof(int), 1, fptr);
  for (i=0; i<(nsym); i++) fread(&symafm[i], sizeof(int), 1, fptr);
  for (i=0; i<(nsym); i++) {
    fread(&symrel[0][0][i], sizeof(int), 1, fptr);
    smap.symrel[i][0][0] = symrel[0][0][i];
    fread(&symrel[1][0][i], sizeof(int), 1, fptr);
    smap.symrel[i][1][0] = symrel[1][0][i];
    fread(&symrel[2][0][i], sizeof(int), 1, fptr);
    smap.symrel[i][2][0] = symrel[2][0][i];
    fread(&symrel[0][1][i], sizeof(int), 1, fptr);
    smap.symrel[i][0][1] = symrel[0][1][i];
    fread(&symrel[1][1][i], sizeof(int), 1, fptr);
    smap.symrel[i][1][1] = symrel[1][1][i];
    fread(&symrel[2][1][i], sizeof(int), 1, fptr);
    smap.symrel[i][2][1] = symrel[2][1][i];
    fread(&symrel[0][2][i], sizeof(int), 1, fptr);
    smap.symrel[i][0][2] = symrel[0][2][i];
    fread(&symrel[1][2][i], sizeof(int), 1, fptr);
    smap.symrel[i][1][2] = symrel[1][2][i];
    fread(&symrel[2][2][i], sizeof(int), 1, fptr);
    smap.symrel[i][2][2] = symrel[2][2][i];
  }
  for (i=0; i<(natom); i++) fread(&typat[i], sizeof(int), 1, fptr);
  for (i=0; i<(nkpt); i++) {
    fread(&kptv, sizeof(double), 1, fptr);
    fread(&kptv, sizeof(double), 1, fptr);
    fread(&kptv, sizeof(double), 1, fptr);
  }
  for (i=0; i<bandtot; i++) fread(&occ, sizeof(double), 1, fptr);
  for (i=0; i<nsym; i++) {
    fread(&tnons[0][i], sizeof(double), 1, fptr);
    smap.tnons[i][0] = tnons[0][i];
    fread(&tnons[1][i], sizeof(double), 1, fptr);
    smap.tnons[i][1] = tnons[1][i];
    fread(&tnons[2][i], sizeof(double), 1, fptr);
    smap.tnons[i][2] = tnons[2][i];
  }
  for (i=0; i<ntypat; i++) fread(&znucltypat[i], sizeof(double), 1, fptr);
  for(i=0; i<natom; i++) {
    type = typat[i]-1;
    gridin->zatomic[i] = (int)znucltypat[type];
  }
  for (i=0; i<nkpt; i++) fread(&wtkv, sizeof(double), 1, fptr);
  fread(&i, sizeof(int), 1, fptr);
  for (k=0; k<npsp; k++) {
    fread(&i, sizeof(int), 1, fptr);
    fread(title, sizeof(char), 132, fptr);
    fread(&znuclpsp, sizeof(double), 1, fptr);
    fread(&zionpsp, sizeof(double), 1, fptr);
    fread(&pspso, sizeof(int), 1, fptr);
    fread(&pspdat, sizeof(int), 1, fptr);
    fread(&pspcod, sizeof(int), 1, fptr);
    fread(&pspxc, sizeof(int), 1, fptr);
    fread(&lmn_size, sizeof(int), 1, fptr);
    fread(&i, sizeof(int), 1, fptr);
  }
  if (usepaw==0) {
    fread(&i, sizeof(int), 1, fptr);
    fread(&residm, sizeof(double), 1, fptr);
    for (i=0; i<natom; i++) {
      fread(&x, sizeof(double), 1, fptr);
      fread(&y, sizeof(double), 1, fptr);
      fread(&z, sizeof(double), 1, fptr);
      if((x<0)&&(-x < 1E-15)) x =0;
      if((y<0)&&(-y < 1E-15)) y =0;
      if((z<0)&&(-z < 1E-15)) z =0;
      gridin->x[i] = x;
      gridin->y[i] = y;
      gridin->z[i] = z;
      gridin->xcart[i] = x*gridin->cella_x+y*gridin->cellb_x+z*gridin->cellc_x;
      gridin->ycart[i] = x*gridin->cella_y+y*gridin->cellb_y+z*gridin->cellc_y;
      gridin->zcart[i] = x*gridin->cella_z+y*gridin->cellb_z+z*gridin->cellc_z;
    }
    fread(&etotal, sizeof(double), 1, fptr);
    gridin->entot = etotal;
    fread(&efermi, sizeof(double), 1, fptr);
    fread(&i, sizeof(int), 1, fptr);
  } else {
    printf("\n  BAD NEWS: PAW pseudopotential detected!\n");
    fprintf(cplog, "\nTerminated because PAW pseudopotentials are not yet supported\n");
    fprintf(cplog, "Suggestion: re-run Abinit with norm-conserving pseudopotentials\n");
    errcount++;
    return 4;
  }
  /* allocating memory for the voxel grid */
  gridin->grid = (double***)malloc((ngx+1)*sizeof(double**));
  for (jx=0; jx<=ngx; jx++) {
    gridin->grid[jx] = (double**)malloc((ngy+1)*sizeof(double*));
    for (jy=0; jy<=ngy; jy++) {
      gridin->grid[jx][jy] = (double*)malloc((ngz+1)*sizeof(double));
    }
  }
  if (nspin==2) {
    gridin2->grid = (double***)malloc((ngx+1)*sizeof(double**));
    for (jx=0; jx<=ngx; jx++) {
      gridin2->grid[jx] = (double**)malloc((ngy+1)*sizeof(double*));
      for (jy=0; jy<=ngy; jy++) {
        gridin2->grid[jx][jy] = (double*)malloc((ngz+1)*sizeof(double));
      }
    }
  }
  /* reading in the voxel grid values */
  fread(&i, sizeof(int), 1, fptr);
  for (jz=0; jz<ngz; jz++) {
    for (jy=0; jy<ngy; jy++) {
      for (jx=0; jx<ngx; jx++) {
        fread(&eigen, sizeof(double), 1, fptr);
        gridin->grid[jx][jy][jz] = eigen;
      }
    }
  }
  FixEdges(gridin);
  fread(&i, sizeof(int), 1, fptr);
  if (nspin==2) {
    fread(&i, sizeof(int), 1, fptr);
    gridin2->volvox = gridin->volvox;
    for (jz=0; jz<ngz; jz++) {
      for (jy=0; jy<ngy; jy++) {
        for (jx=0; jx<ngx; jx++) {
          fread(&eigen, sizeof(double), 1, fptr);
          gridin2->grid[jx][jy][jz] = eigen;
          /* down density = total density - up density */
          gridin->grid[jx][jy][jz] = gridin->grid[jx][jy][jz]-gridin2->grid[jx][jy][jz];
        }
      }
    }
    FixEdges(gridin);
    FixEdges(gridin2);
    fread(&i, sizeof(int), 1, fptr);
  }
  fclose(fptr);
  return 0;
}

int ReadXSF(char binname[STRMAX], struct CrystData * gridin) {
  /* called by: MapNonloc, ReadAll*/
  /* calls: FinishLine */
  int j, jx, jy, jz,temp1,temp2,temp3;
  double voxelValue, cellvol=0.0;
  FILE * f1;
  f1 = fopen (binname,"r");
  FinishLine(f1);
  FinishLine(f1);
  FinishLine(f1);
  fscanf(f1,"%lf %lf %lf",&gridin->cella_x, &gridin->cella_y, &gridin->cella_z);
  fscanf(f1,"%lf %lf %lf",&gridin->cellb_x, &gridin->cellb_y, &gridin->cellb_z);
  fscanf(f1,"%lf %lf %lf",&gridin->cellc_x, &gridin->cellc_y, &gridin->cellc_z);
  FinishLine(f1);
  FinishLine(f1);
  fscanf(f1,"%d",&gridin->nion);
  FinishLine(f1);
  /* Cross product (Angstrom^3) */
  cellvol = (gridin->cella_x*(gridin->cellb_y*gridin->cellc_z-gridin->cellb_z*gridin->cellc_y)-gridin->cella_y*(gridin->cellb_x*gridin->cellc_z-gridin->cellb_z*gridin->cellc_x)+gridin->cella_z*(gridin->cellb_x*gridin->cellc_y-gridin->cellb_y*gridin->cellc_x));
  for(j=0; j<gridin->nion; j++) {
    fscanf(f1,"%d %lf %lf %lf",&gridin->zatomic[j], &gridin->xcart[j], &gridin->ycart[j], &gridin->zcart[j]);
    FinishLine(f1);
  }
  FinishLine(f1);
  for(j=0; j<gridin->nion; j++) {
    fscanf(f1,"%d %lf %lf %lf",&gridin->zatomic[j], &gridin->xcart[j], &gridin->ycart[j], &gridin->zcart[j]);
    FinishLine(f1);
  }
  FinishLine(f1);
  FinishLine(f1);
  FinishLine(f1);
  fscanf(f1,"%d %d %d",&temp1,&temp2,&temp3);
  ngx = temp1-1;
  ngy = temp2-1;
  ngz = temp3-1;
  gridin->volvox = cellvol/(ngx*ngy*ngz)/pow(R_BOHR, 3);    /* bohr^3 */
  FinishLine(f1);
  FinishLine(f1);
  fscanf(f1,"%lf %lf %lf",&gridin->cella_x, &gridin->cella_y, &gridin->cella_z);
  fscanf(f1,"%lf %lf %lf",&gridin->cellb_x, &gridin->cellb_y, &gridin->cellb_z);
  fscanf(f1,"%lf %lf %lf",&gridin->cellc_x, &gridin->cellc_y, &gridin->cellc_z);
  /* Allocating memory for the voxel grid */
  gridin->grid = (double***)malloc((ngx+1)*sizeof(double**));
  for (jx=0; jx<=ngx; jx++) {
    gridin->grid[jx] = (double**)malloc((ngy+1)*sizeof(double*));
    for (jy=0; jy<=ngy; jy++) {
      gridin->grid[jx][jy] = (double*)malloc((ngz+1)*sizeof(double));
    }
  }
  /* Reading in voxel values */
  for(jz=0; jz<=ngz; jz++) {
    for(jy=0; jy<=ngy; jy++) {
      for(jx=0; jx<=ngx; jx++) {
        fscanf(f1,"%lf",&voxelValue);
        gridin->grid[jx][jy][jz] = voxelValue;
      }
    }
  }
  fclose(f1);
  return 0;
}


int ReadCHGCAR(char binname[STRMAX], struct CrystData * gridin) {
  int j, jx, jy, jz,temp1,temp2,temp3;
  int stop = 0,check,nions_type,nions=0;
  double voxelValue, voxelsum=0.0, cellvol=0.0;
  char str [100];
  FILE * f1;
  f1 = fopen (binname,"r");
  if(f1==NULL) {
        printf("%s file not found.\n",binname);
        exit(0);
  }
  FinishLine(f1);
  FinishLine(f1);
  fscanf(f1,"%lf %lf %lf",&gridin->cella_x, &gridin->cella_y, &gridin->cella_z);
  fscanf(f1,"%lf %lf %lf",&gridin->cellb_x, &gridin->cellb_y, &gridin->cellb_z);
  fscanf(f1,"%lf %lf %lf",&gridin->cellc_x, &gridin->cellc_y, &gridin->cellc_z);
  while(stop==0) {
    check = fscanf(f1,"%d",&nions_type);
    if(check==1) nions+=nions_type;
    else stop=1;
  }
  gridin->nion=nions;
  fscanf(f1,"%s",str);
  /* Cross product (Angstrom^3) */
  cellvol = (gridin->cella_x*(gridin->cellb_y*gridin->cellc_z-gridin->cellb_z*gridin->cellc_y)-gridin->cella_y*(gridin->cellb_x*gridin->cellc_z-gridin->cellb_z*gridin->cellc_x)+gridin->cella_z*(gridin->cellb_x*gridin->cellc_y-gridin->cellb_y*gridin->cellc_x));
  for(j=0; j<gridin->nion; j++) {
    fscanf(f1,"%lf %lf %lf", &gridin->xcart[j], &gridin->ycart[j], &gridin->zcart[j]);
  }
  fscanf(f1,"%d %d %d",&temp1,&temp2,&temp3);
  ngx = temp1;
  ngy = temp2;
  ngz = temp3;
  gridin->volvox = cellvol/(ngx*ngy*ngz)/pow(R_BOHR, 3);    /* bohr^3 */
  /* Reading in voxel values */
  for(jz=0; jz<ngz; jz++) {
    for(jy=0; jy<ngy; jy++) {
      for(jx=0; jx<ngx; jx++) {
        fscanf(f1,"%lf",&voxelValue);
        gridin->grid[jx][jy][jz] = voxelValue;
      }
    }
  }
  fclose(f1);
  return 0;
}


int OutputXSF(FILE * fptr, struct CrystData * gridref, struct CrystData * gridout) {
  /* called by: main, CoreUnwarp, Den2XSF, MapEntot, OutputWeight, MapEwald, CalcCP, PrintResults */
  /* calls: FixEdges */
  int i=0, jx=0, jy=0, jz=0, linecount=0;
  FixEdges(gridout);
  fprintf(fptr, " DIM-GROUP\n");
  fprintf(fptr, " 3  1\n");
  fprintf(fptr, " PRIMVEC\n");
  fprintf(fptr, "%20.14f  %20.14f  %20.14f\n", gridref->cella_x*R_BOHR, gridref->cella_y*R_BOHR, gridref->cella_z*R_BOHR);
  fprintf(fptr, "%20.14f  %20.14f  %20.14f\n", gridref->cellb_x*R_BOHR, gridref->cellb_y*R_BOHR, gridref->cellb_z*R_BOHR);
  fprintf(fptr, "%20.14f  %20.14f  %20.14f\n", gridref->cellc_x*R_BOHR, gridref->cellc_y*R_BOHR, gridref->cellc_z*R_BOHR);
  fprintf(fptr, " PRIMCOORD\n");
  fprintf(fptr, "%12d  1\n", gridref->nion);
  for (i=0; i<gridref->nion; i++) fprintf(fptr, "%9d  %20.14f  %20.14f  %20.14f\n",
    gridref->zatomic[i], gridref->xcart[i]*R_BOHR, gridref->ycart[i]*R_BOHR, gridref->zcart[i]*R_BOHR);
  fprintf(fptr, " ATOMS\n");
  for (i=0; i<gridref->nion; i++) fprintf(fptr, "%9d  %20.14f  %20.14f  %20.14f\n",
    gridref->zatomic[i], gridref->xcart[i]*R_BOHR, gridref->ycart[i]*R_BOHR, gridref->zcart[i]*R_BOHR);
  fprintf(fptr, " BEGIN_BLOCK_DATAGRID3D\n");
  fprintf(fptr, " datagrids\n");
  fprintf(fptr, " DATAGRID_3D_DENSITY\n");
  fprintf(fptr, "%12d  %12d  %12d\n", ngx+1, ngy+1, ngz+1);
  fprintf(fptr, " 0.0  0.0  0.0\n");
  fprintf(fptr, "%20.14f  %20.14f  %20.14f\n", gridref->cella_x*R_BOHR, gridref->cella_y*R_BOHR, gridref->cella_z*R_BOHR);
  fprintf(fptr, "%20.14f  %20.14f  %20.14f\n", gridref->cellb_x*R_BOHR, gridref->cellb_y*R_BOHR, gridref->cellb_z*R_BOHR);
  fprintf(fptr, "%20.14f  %20.14f  %20.14f\n", gridref->cellc_x*R_BOHR, gridref->cellc_y*R_BOHR, gridref->cellc_z*R_BOHR);
  for (jz=0; jz<=ngz; jz++) {
    for (jy=0; jy<=ngy; jy++) {
      for (jx=0; jx<=ngx; jx++) {
        linecount++;
        fprintf(fptr, "  %20.14f", gridout->grid[jx][jy][jz]);
        if (linecount==6) {
          fprintf(fptr, "\n");
          linecount = 0;
        }
      }
    }
  }
  fprintf(fptr, " END_DATAGRID_3D\n");
  fprintf(fptr, " END_BLOCK_DATAGRID3D\n");
  return 0;
}

int Den2XSF(char name[STRMAX], int num, char type[STRMAX], struct CrystData * gridin, struct CrystData * gridin2) {
  /* called by: ReadAll, ReadGradient */
  /* calls: Bin2XSF, OutputXSF */
  char denfile[STRMAX], xsffile[STRMAX];
  int check=0;
  FILE * fptr;
  snprintf(denfile, STRMAX, "%s_DS%d_%s", name, num, type);
  check = Bin2XSF(denfile, gridin, gridin2);  /* always read from binary file */
  if (check!=0) return 1;
  if (printbin==1) {
    strncpy(xsffile, denfile, STRMAX);
    strncat(xsffile, ".xsf", STRMAX);
    fptr = fopen(xsffile, "w");
    OutputXSF(fptr, gridin, gridin);
    fclose(fptr);
    if (nspin==2) {
      strncpy(xsffile, denfile, STRMAX);
      strncat(xsffile, "2.xsf", STRMAX);
      fptr = fopen(xsffile, "w");
      OutputXSF(fptr, gridin, gridin2);
      fclose(fptr);
    }
  }
  return 0;
}

int ReadNonlocAtom(struct CrystData * gridin) {
  /* called by: main */
  /* calls: FinishLine */
  int i=0, j=0;
  char filename[STRMAX], str[STRMAX];
  FILE * fptr;
  printf("  Enter the name of the file containing nonlocal atomic energies: ");
  scanf("%s", filename);
  fptr = fopen(filename, "r");
  FinishLine(fptr);
  /* data for DS1 */
  for (i=0; i<gridin->nion; i++) {
    fscanf(fptr, "%s %lf %lf %lf %lf %lf", str, &E_NL[i][0][0], &E_NL[i][1][0], &E_NL[i][2][0], &E_NL[i][3][0], &E_NL[i][4][0]);
  }
  FinishLine(fptr);
  FinishLine(fptr);
  /* data for DS3 */
  for (i=0; i<gridin->nion; i++) {
    fscanf(fptr, "%s %lf %lf %lf %lf %lf", str, &E_NL[i][0][1], &E_NL[i][1][1], &E_NL[i][2][1], &E_NL[i][3][1], &E_NL[i][4][1]);
  }
  fclose(fptr);
  for (i=0; i<gridin->nion; i++) {
    for (j=i; j<gridin->nion; j++) {
      if (smap.equiv[i][j]==1) {
        E_NL_equiv[i][0] += E_NL[j][0][0];
        E_NL_equiv[i][1] += E_NL[j][0][1];
      }
    }
    E_NL_equiv[i][0] = E_NL_equiv[i][0]/smap.nequiv[i];
    E_NL_equiv[i][1] = E_NL_equiv[i][1]/smap.nequiv[i];
  }

  return 0;
}

double pspElec(int Z) {
  /* called by: ReadProfile */
  /* calls: none */
  double electrons[119] = {
    0.0,1.0,2.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,
    1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,
    1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,
    1.0,2.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0,21.0,22.0,23.0,24.0,25.0,
    12.0,5.0,6.0,7.0,8.0,9.0,10.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,
    0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
    0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  return electrons[Z];
}

double pspElecSC(int Z) {
  /* called by: ReadProfile */
  /* calls: none */
  double electrons[119] = {
    0.0,1.0,2.0,3.0,4.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,3.0,4.0,5.0,6.0,7.0,8.0,
    9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,11.0,12.0,13.0,4.0,5.0,6.0,7.0,8.0,
    9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,11.0,12.0,13.0,4.0,5.0,6.0,7.0,8.0,
    9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0,21.0,22.0,23.0,24.0,25.0,
    12.0,13.0,14.0,15.0,16.0,17.0,18.0,11.0,12.0,13.0,4.0,5.0,6.0,7.0,8.0,
    0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
    0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  return electrons[Z];
}



int ReadProfile(struct CrystData * gridin) {
  /* called by: main */
  /* calls: ElementName, FinishLine, pspElec, pspElecSC */
  char element[STRMAX], profname[STRMAX];
  int i=0, count,j=0, k=0, stop[NIONMAX];
  double elec=0.0, temp_zion[NIONMAX];
  FILE * fptr;
  count=0;
  for (i=0; i<gridin->nion; i++) stop[i] = 0;
  for (i=0; i<gridin->nion; i++) {
    if (stop[i]==1) continue;
    ElementName(gridin->zatomic[i], element);
    if (smap.nequiv[i]>1) {
      fprintf(cplog, "Name of the Hirshfeld profile for atom #%d (%s, %d equivalent sites): ",
        i+1, element, smap.nequiv[i]);
    } else {
      fprintf(cplog, "Name of the Hirshfeld profile for atom #%d (%s, %d site): ",
        i+1, element, smap.nequiv[i]);
    }
    retry_profile:
    if (smap.nequiv[i]>1) {
      printf("  Name of the Hirshfeld profile for atom #%d (%s, %d equivalent sites): ",
        i+1, element, smap.nequiv[i]);
    } else {
      printf("  Name of the Hirshfeld profile for atom #%d (%s, %d site): ",
        i+1, element, smap.nequiv[i]);
    }
    fptr = fopen(profile_filenames[count], "r");
    if (fptr==NULL) {
      printf("  Atomic radial density profile %s not found!\n", profile_filenames[count]);
      exit(1);
    }
    fprintf(cplog, "%s\n", profile_filenames[count]);
    printf("%s\n", profile_filenames[count]);
    if (E_Ewald_option !=0) {
      vosctype[i]=psptypes[count];
      printf("  Pseudopotential type for atom #%d (%s) [0=valence-only] [1=semicore]: %d \n", i+1, element,vosctype[i]);
      fprintf(cplog,"  Pseudopotential type for atom #%d (%s) [0=valence-only] [1=semicore]: %d \n", i+1, element,vosctype[i]);
      elec=0.0;
      if (vosctype[i] == 0) {
        elec = pspElec(deneq.zatomic[i]);
      }
      else {
        elec = pspElecSC(deneq.zatomic[i]);
      }
      temp_zion[i] = elec;
      sc_elec[i]=sccounts[count];
      printf("  Number of the %.1lf electrons on atom %d treated as semicore: %lf\n", temp_zion[i], i,sc_elec[i]);
      fprintf(cplog,"  Number of the %.1lf electrons on atom %d treated as semicore: %lf\n", temp_zion[i], i,sc_elec[i]);
      vo_elec[i] = temp_zion[i] - sc_elec[i];
    }
    for (j=i; j<gridin->nion; j++) {
      if (smap.equiv[i][j]!=1) continue;
      vosctype[j] = vosctype[i];
      sc_elec[j] = sc_elec[i];
      vo_elec[j] = vo_elec[i];
      k = 0;
      while (FinishLine(fptr)==0) {
        fscanf(fptr, "%lf %lf", &rprofile[j][k], &rhoprofile[j][k]);
        k++;
        if (k==NPOINTMAX) {
          printf("\n  BAD NEWS: Density profile %s is bigger than expected!\n", profname);
          fprintf(cplog, "\nTerminated because file %s is too long\n", profname);
          fprintf(cplog, "Suggestion: increase NPOINTMAX and recompile\n");
          errcount++;
          return 1;
        }
      }
      fclose(fptr);
      fptr = fopen(profile_filenames[count], "r");
      logint[j] = log(rprofile[j][0]);
      logslope[j] = (log(rprofile[j][k-1])-log(rprofile[j][0]))/((double)k-1.0);
      prof_nmax[j] = k;
      stop[j] = 1;
    }
    fclose(fptr);
    count++;
  }
  return 0;
}

int OutputWeight(struct CrystData * gridin, struct ContactVol * map, struct CrystData * gridout) {
  /* called by: PrintResults */
  /* calls: FixEdges, OutputXSF */
  char filename[STRMAX];
  int atom=0, i=0, jx=0, jy=0, jz=0;
  double tempcp=0.0;
  FILE * fptr;
  for(i=0; i<map->neighcount[jx][jy][jz]; i++) {
    printf("  Creating voxel weight map for atom %d\n", i+1);
    for(jz=0; jz<ngz; jz++) {
      for(jy=0; jy<ngy; jy++) {
        for(jx=0; jx<ngx; jx++) {
          atom = map->ionmap[i][jx][jy][jz]&127;
          tempcp = 0.5*(map->swj[jx][jy][jz]-map->wj[atom][jx][jy][jz])*map->wj[atom][jx][jy][jz];
          gridout->grid[jx][jy][jz] = tempcp/map->swjk[jx][jy][jz];
        }
      }
    }
    FixEdges(gridout);
    snprintf(filename, STRMAX, "%s-voxelweight-%d.xsf", cpoutname, i+1);
    fptr = fopen(filename, "w");
    OutputXSF(fptr, gridin, gridout);
    fclose(fptr);
  }
  return 0;
}


/* CORE UNWARP FUNCTIONS */

double CubicInterpolation(double x, double y1, double y2, double y3, double y4) {
  /* called by: TricubicInterpolation */
  /* calls: none */
  double a=0.0, b=0.0, c=0.0, d=y2, apc=0.0, qapc=0.0;
  b = 0.5*(y1+y3)-y2;
  apc = 0.5*(y3-y1);
  qapc = 0.5*(y4-4*b-d);
  a = (qapc-apc)/3.0;
  c = apc-a;
  return (((a*x+b)*x+c)*x+d);  /* y = ax^3+bx^2+cx+d */
}

double TricubicInterpolation(struct CrystData * gridin, struct CrystData * gridref, double x, double y, double z) {
  /* called by: CoreUnwarp */
  /* calls: CubicInterpolation */
  int a1=0, a2=0, a3=0, a4=0, b1=0, b2=0, b3=0, b4=0, c1=0, c2=0, c3=0, c4=0;
  double xf=0.0, yf=0.0, zf=0.0, y1=0.0, y2=0.0, y3=0.0, y4=0.0, y11=0.0, y12=0.0, y13=0.0, y14=0.0, y21=0.0;
  double y22=0.0, y23=0.0, y24=0.0, y31=0.0, y32=0.0, y33=0.0, y34=0.0, y41=0.0, y42=0.0, y43=0.0, y44=0.0;
  double delta=0.0, interpvalue=0.0, voxspacing_xf=0.0, voxspacing_yf=0.0, voxspacing_zf=0.0;
  gsl_matrix * xf_to_r = gsl_matrix_alloc(3, 3);
  gsl_vector * rf = gsl_vector_alloc(3);
  gsl_vector * rcart = gsl_vector_alloc(3);
  /* convert cartesian x, y, z to fractional xf, yf, zf */
  gsl_matrix_set(xf_to_r, 0, 0, gridref->cella_x);
  gsl_matrix_set(xf_to_r, 0, 1, gridref->cellb_x);
  gsl_matrix_set(xf_to_r, 0, 2, gridref->cellc_x);
  gsl_matrix_set(xf_to_r, 1, 0, gridref->cella_y);
  gsl_matrix_set(xf_to_r, 1, 1, gridref->cellb_y);
  gsl_matrix_set(xf_to_r, 1, 2, gridref->cellc_y);
  gsl_matrix_set(xf_to_r, 2, 0, gridref->cella_z);
  gsl_matrix_set(xf_to_r, 2, 1, gridref->cellb_z);
  gsl_matrix_set(xf_to_r, 2, 2, gridref->cellc_z);
  gsl_vector_set(rcart, 0, x);
  gsl_vector_set(rcart, 1, y);
  gsl_vector_set(rcart, 2, z);
  gsl_linalg_HH_solve(xf_to_r, rcart, rf);
  xf = gsl_vector_get(rf, 0);
  yf = gsl_vector_get(rf, 1);
  zf = gsl_vector_get(rf, 2);
  gsl_matrix_free(xf_to_r);
  gsl_vector_free(rf);
  gsl_vector_free(rcart);
  /* translate xf, yf, zf into central unit cell */
  xf -= floor(xf);
  yf -= floor(yf);
  zf -= floor(zf);
  if (xf>0.999999) xf = 0.0;
  if (yf>0.999999) yf = 0.0;
  if (zf>0.999999) zf = 0.0;
  /* determine voxels around xf, yf, zf */
  voxspacing_xf = 1.0/(double)ngx;
  voxspacing_yf = 1.0/(double)ngy;
  voxspacing_zf = 1.0/(double)ngz;
  a2 = (int)(xf/voxspacing_xf);
  a1 = (a2-1+ngx)%ngx;
  a3 = (a2+1+ngx)%ngx;
  a4 = (a2+2+ngx)%ngx;
  b2 = (int)(yf/voxspacing_yf);
  b1 = (b2-1+ngy)%ngy;
  b3 = (b2+1+ngy)%ngy;
  b4 = (b2+2+ngy)%ngy;
  c2 = (int)(zf/voxspacing_zf);
  c1 = (c2-1+ngz)%ngz;
  c3 = (c2+1+ngz)%ngz;
  c4 = (c2+2+ngz)%ngz;
  delta = (xf-(double)a2*voxspacing_xf)/voxspacing_xf;
  y11 = CubicInterpolation(delta, gridin->grid[a1][b1][c1], gridin->grid[a2][b1][c1],
    gridin->grid[a3][b1][c1], gridin->grid[a4][b1][c1]);
  y12 = CubicInterpolation(delta, gridin->grid[a1][b1][c2], gridin->grid[a2][b1][c2],
    gridin->grid[a3][b1][c2], gridin->grid[a4][b1][c2]);
  y13 = CubicInterpolation(delta, gridin->grid[a1][b1][c3], gridin->grid[a2][b1][c3],
    gridin->grid[a3][b1][c3], gridin->grid[a4][b1][c3]);
  y14 = CubicInterpolation(delta, gridin->grid[a1][b1][c4], gridin->grid[a2][b1][c4],
    gridin->grid[a3][b1][c4], gridin->grid[a4][b1][c4]);
  y21 = CubicInterpolation(delta, gridin->grid[a1][b2][c1], gridin->grid[a2][b2][c1],
    gridin->grid[a3][b2][c1], gridin->grid[a4][b2][c1]);
  y22 = CubicInterpolation(delta, gridin->grid[a1][b2][c2], gridin->grid[a2][b2][c2],
    gridin->grid[a3][b2][c2], gridin->grid[a4][b2][c2]);
  y23 = CubicInterpolation(delta, gridin->grid[a1][b2][c3], gridin->grid[a2][b2][c3],
    gridin->grid[a3][b2][c3], gridin->grid[a4][b2][c3]);
  y24 = CubicInterpolation(delta, gridin->grid[a1][b2][c4], gridin->grid[a2][b2][c4],
    gridin->grid[a3][b2][c4], gridin->grid[a4][b2][c4]);
  y31 = CubicInterpolation(delta, gridin->grid[a1][b3][c1], gridin->grid[a2][b3][c1],
    gridin->grid[a3][b3][c1], gridin->grid[a4][b3][c1]);
  y32 = CubicInterpolation(delta, gridin->grid[a1][b3][c2], gridin->grid[a2][b3][c2],
    gridin->grid[a3][b3][c2], gridin->grid[a4][b3][c2]);
  y33 = CubicInterpolation(delta, gridin->grid[a1][b3][c3], gridin->grid[a2][b3][c3],
    gridin->grid[a3][b3][c3], gridin->grid[a4][b3][c3]);
  y34 = CubicInterpolation(delta, gridin->grid[a1][b3][c4], gridin->grid[a2][b3][c4],
    gridin->grid[a3][b3][c4], gridin->grid[a4][b3][c4]);
  y41 = CubicInterpolation(delta, gridin->grid[a1][b4][c1], gridin->grid[a2][b4][c1],
    gridin->grid[a3][b4][c1], gridin->grid[a4][b4][c1]);
  y42 = CubicInterpolation(delta, gridin->grid[a1][b4][c2], gridin->grid[a2][b4][c2],
    gridin->grid[a3][b4][c2], gridin->grid[a4][b4][c2]);
  y43 = CubicInterpolation(delta, gridin->grid[a1][b4][c3], gridin->grid[a2][b4][c3],
    gridin->grid[a3][b4][c3], gridin->grid[a4][b4][c3]);
  y44 = CubicInterpolation(delta, gridin->grid[a1][b4][c4], gridin->grid[a2][b4][c4],
    gridin->grid[a3][b4][c4], gridin->grid[a4][b4][c4]);
  delta = (yf-b2*voxspacing_yf)/voxspacing_yf;
  y1 = CubicInterpolation(delta, y11, y21, y31, y41);
  y2 = CubicInterpolation(delta, y12, y22, y32, y42);
  y3 = CubicInterpolation(delta, y13, y23, y33, y43);
  y4 = CubicInterpolation(delta, y14, y24, y34, y44);
  delta = (zf-(double)c2*voxspacing_zf)/voxspacing_zf;
  interpvalue = CubicInterpolation(delta, (double)y1, (double)y2, (double)y3, (double)y4);
  return interpvalue;
}

int CoreUnwarp(struct CrystData * gridinup, struct CrystData * gridin, struct CrystData * gridindn, struct CrystData * gridoutup, struct CrystData * gridoutdn) {
  /* called by: main */
  /* calls: FixEdges, Getwj, OutputXSF, TricubicInterpolation */
  char str[STRMAX];
  int atom=0, count=0, ngcount=0, ngp=0, ngp0=0;
  int i=0, jx=0, jy=0, jz=0, jx1=0, jy1=0, jz1=0, jx2=0, jy2=0, jz2=0, ka=0, kb=0, kc=0;
  int gridx=ngx+1, gridy=ngy+1, gridz=ngz+1;
  double dist=0.0, encoreup=0.0, encoredn=0.0, endelta=0.0, wfrac=0.0, wsum=0.0;
  double upscale=0.0, dnscale=0.0, dxnewup=0.0, dynewup=0.0, dznewup=0.0, dxnewdn=0.0, dynewdn=0.0, dznewdn=0.0;
  double volnew=0.0, voltotup=0.0, voltotdn=0.0, xtran=0.0, ytran=0.0, ztran=0.0;
  double voxcenter_x=0.0, voxcenter_y=0.0, voxcenter_z=0.0, xc=0.0, yc=0.0, zc=0.0, xf=0.0, yf=0.0, zf=0.0;
  double xdist1=0.0, ydist1=0.0, zdist1=0.0, xdist2=0.0, ydist2=0.0, zdist2=0.0, xdist3=0.0, ydist3=0.0, zdist3=0.0;
  FILE * fptr;

  xnewup = (double***)malloc(gridx*sizeof(double**));
  ynewup = (double***)malloc(gridx*sizeof(double**));
  znewup = (double***)malloc(gridx*sizeof(double**));
  xnewdn = (double***)malloc(gridx*sizeof(double**));
  ynewdn = (double***)malloc(gridx*sizeof(double**));
  znewdn = (double***)malloc(gridx*sizeof(double**));
  for (jx=0; jx<gridx; jx++) {
    xnewup[jx] = (double**)malloc(gridy*sizeof(double*));
    ynewup[jx] = (double**)malloc(gridy*sizeof(double*));
    znewup[jx] = (double**)malloc(gridy*sizeof(double*));
    xnewdn[jx] = (double**)malloc(gridy*sizeof(double*));
    ynewdn[jx] = (double**)malloc(gridy*sizeof(double*));
    znewdn[jx] = (double**)malloc(gridy*sizeof(double*));
    for(jy=0; jy<gridy; jy++) {
      xnewup[jx][jy] = (double*)malloc(gridz*sizeof(double));
      ynewup[jx][jy] = (double*)malloc(gridz*sizeof(double));
      znewup[jx][jy] = (double*)malloc(gridz*sizeof(double));
      xnewdn[jx][jy] = (double*)malloc(gridz*sizeof(double));
      ynewdn[jx][jy] = (double*)malloc(gridz*sizeof(double));
      znewdn[jx][jy] = (double*)malloc(gridz*sizeof(double));
    }
  }


  for (jz=0; jz<=ngz; jz++) {
    for (jy=0; jy<=ngy; jy++) {
      for (jx=0; jx<=ngx; jx++) {
        core.grid[jx][jy][jz] = 0.0;   /* placeholder grid */
        pothi.grid[jx][jy][jz] = 0.0;  /* placeholder grid */
        potlo.grid[jx][jy][jz] = 0.0;  /* placeholder grid */
        temp.grid[jx][jy][jz] = 0.0;   /* placeholder grid */
      }
    }
  }
  upscale = pow(gridinup->volvox/gridin->volvox, ONETHIRD)-1.0;
  dnscale = pow(gridindn->volvox/gridin->volvox, ONETHIRD)-1.0;
  fprintf(cplog, "Began core unwarping\n");
  printf("  Interpolating between %d x %d x %d = %d voxels\n", ngx, ngy, ngz, ngx*ngy*ngz);
  fprintf(cplog, "Interpolated %d x %d x %d = %d voxels\n", ngx, ngy, ngz, ngx*ngy*ngz);
  printf("0%%");
  fflush(stdout);
  /* for every voxel in the unit cell */
  for (jz=0; jz<ngz; jz++) {
    /* xf, yf, zf fractional coordinates */
    zf = (double)jz/(double)ngz;
    for (jy=0; jy<ngy; jy++) {
      yf = (double)jy/(double)ngy;
      for (jx=0; jx<ngx; jx++) {
        xf = (double)jx/(double)ngx;
        /* voxel centers in cartesian coordinates */
        voxcenter_x = xf*gridin->cella_x+yf*gridin->cellb_x+zf*gridin->cellc_x;
        voxcenter_y = xf*gridin->cella_y+yf*gridin->cellb_y+zf*gridin->cellc_y;
        voxcenter_z = xf*gridin->cella_z+yf*gridin->cellb_z+zf*gridin->cellc_z;
        /* prints percent completion */
        ngcount++;
        ngp = ngcount*100/(ngx*ngy*ngz);
        if (ngp!=ngp0) {
          printf("\r%d%%", ngp);
          fflush(stdout);
        }
        ngp0 = ngp;
        /* for every unit cell in the supercell */
        hmap.neighcount = 0;
        wsum = 0.0;
        for (ka=-kam; ka<=kam; ka++) {
          for (kb=-kbm; kb<=kbm; kb++) {
            for (kc=-kcm; kc<=kcm; kc++) {
              /* for every atom in cartesian coordinates */
              for (atom=0; atom<gridin->nion; atom++) {
                xc = gridin->xcart[atom]+ka*gridin->cella_x+kb*gridin->cellb_x+kc*gridin->cellc_x;
                yc = gridin->ycart[atom]+ka*gridin->cella_y+kb*gridin->cellb_y+kc*gridin->cellc_y;
                zc = gridin->zcart[atom]+ka*gridin->cella_z+kb*gridin->cellb_z+kc*gridin->cellc_z;
                dist = sqrt((voxcenter_x-xc)*(voxcenter_x-xc)+(voxcenter_y-yc)*
                  (voxcenter_y-yc)+(voxcenter_z-zc)*(voxcenter_z-zc));
                /* if voxel is close enough to translated atom */
                if (dist<R_MAX) {
                  /* determine hirshfeld-like weight for atom */
                  count = hmap.neighcount;
                  hmap.wj[count] = Getwj(atom, dist);
                  if (hmap.wj[count]>hmap.hirsh_weight[atom][jx][jy][jz]) {
                    hmap.hirsh_weight[atom][jx][jy][jz] = hmap.wj[count]*(4*PI);
                  }
                  if (hmap.wj[count]==-1000.0) return 1;
                  wsum += hmap.wj[count];
                  hmap.atomid[count] = atom;
                  hmap.xcart[count] = xc;
                  hmap.ycart[count] = yc;
                  hmap.zcart[count] = zc;
                  hmap.neighcount++;
                  if (hmap.neighcount==NEQVOX) {
                    printf("\n  BAD NEWS: The number of nearby atoms exceeds %d!\n", NEQVOX);
                    fprintf(cplog, "\nTerminated because the number of atoms near voxel %d %d %d is larger than %d\n",
                      jx, jy, jz, NEQVOX);
                    fprintf(cplog, "Suggestion: increase NEQVOX or decrease R_MAX and recompile\n");
                    errcount++;
                    return 2;
                  }
                }
              }
            }
          }
        }
        /* determine shift in voxel position for expanded and contracted datasets */
        dxnewup = 0.0;
        dynewup = 0.0;
        dznewup = 0.0;
        dxnewdn = 0.0;
        dynewdn = 0.0;
        dznewdn = 0.0;
        /* for a voxel near an atomic nuclei, wfrac=1 for that atom and dxnewup is the
           x-vector between the equilibrium voxel center and the upscaled voxel such that
           xnewup = equilibrium voxel position (cf. JCTC 2014, Eq 6) */
        for (i=0; i<hmap.neighcount; i++) {
          wfrac = hmap.wj[i]/wsum;  /* cf. JCTC 2014, Eq 8 */
          if(wsum==0.0) {
            printf("\n  BAD NEWS: A grid voxel was not interpolated!\n");
            fprintf(cplog, "\nTerminated because sum of hirshfeld-like atomic weights for voxel %d %d %d = zero\n", jx, jy, jz);
            fprintf(cplog, "Suggestion: check the atomic density profiles or increase R_MAX and recompile\n");
            errcount++;
            return 4;
          }
          /* vector distance between each voxel center and nearby atoms (cf. delta_r in Eq 7)
             as a hirshfeld-weighted sum of each atom's influence, scaled either up or down
             for the expanded or contracted dataset, respectively */
          dxnewup += wfrac*(hmap.xcart[i]-voxcenter_x)*upscale;
          dynewup += wfrac*(hmap.ycart[i]-voxcenter_y)*upscale;
          dznewup += wfrac*(hmap.zcart[i]-voxcenter_z)*upscale;
          dxnewdn += wfrac*(hmap.xcart[i]-voxcenter_x)*dnscale;
          dynewdn += wfrac*(hmap.ycart[i]-voxcenter_y)*dnscale;
          dznewdn += wfrac*(hmap.zcart[i]-voxcenter_z)*dnscale;
          /* hirshfeld charge analysis */
          hmap.chg[hmap.atomid[i]] += wfrac*gridin->grid[jx][jy][jz]*gridin->volvox;
        }
        /* expanded-volume cartesian voxel position plus the hirshfeld-weighted correction (cf. JCTC 2014, Eq 9a) */
        xnewup[jx][jy][jz] = xf*gridinup->cella_x+yf*gridinup->cellb_x+zf*gridinup->cellc_x+dxnewup;
        ynewup[jx][jy][jz] = xf*gridinup->cella_y+yf*gridinup->cellb_y+zf*gridinup->cellc_y+dynewup;
        znewup[jx][jy][jz] = xf*gridinup->cella_z+yf*gridinup->cellb_z+zf*gridinup->cellc_z+dznewup;
        /* new interpolated energy at shifted voxel center */
        pothi.grid[jx][jy][jz] = TricubicInterpolation(gridoutup, gridinup,
          xnewup[jx][jy][jz], ynewup[jx][jy][jz], znewup[jx][jy][jz]);

        /* cf. JCTC 2014, Eq 9b */
        xnewdn[jx][jy][jz] = xf*gridindn->cella_x+yf*gridindn->cellb_x+zf*gridindn->cellc_x+dxnewdn;
        ynewdn[jx][jy][jz] = xf*gridindn->cella_y+yf*gridindn->cellb_y+zf*gridindn->cellc_y+dynewdn;
        znewdn[jx][jy][jz] = xf*gridindn->cella_z+yf*gridindn->cellb_z+zf*gridindn->cellc_z+dznewdn;
        potlo.grid[jx][jy][jz] = TricubicInterpolation(gridoutdn, gridindn,
          xnewdn[jx][jy][jz], ynewdn[jx][jy][jz], znewdn[jx][jy][jz]);
      }
    }
  }
  printf(" Finished\n");
  /* reassign voxels with interpolated energies and shifted-grid volumes */
  for (jz=0; jz<ngz; jz++) {
    for (jy=0; jy<ngy; jy++) {
      for (jx=0; jx<ngx; jx++) {
        xtran = 0;
        ytran = 0;
        ztran = 0;
        /* jx1 and jx2 are neighboring voxel indices to jx */
        jx1 = (jx==0)? (xtran = 1, ngx-1) : (jx-1);
        jy1 = (jy==0)? (ytran = 1, ngy-1) : (jy-1);
        jz1 = (jz==0)? (ztran = 1, ngz-1) : (jz-1);
        jx2 = (jx==ngx-1)? (xtran = 1, 0) : (jx+1);
        jy2 = (jy==ngy-1)? (ytran = 1, 0) : (jy+1);
        jz2 = (jz==ngz-1)? (ztran = 1, 0) : (jz+1);
        /* assuming insignificant change in angles, the new volume of a voxel is
           based on the shifted distance between neighboring voxel centers */
        xdist1 = xnewup[jx2][jy][jz]-xnewup[jx1][jy][jz]+xtran*gridinup->cella_x;
        ydist1 = ynewup[jx2][jy][jz]-ynewup[jx1][jy][jz]+xtran*gridinup->cella_y;
        zdist1 = znewup[jx2][jy][jz]-znewup[jx1][jy][jz]+xtran*gridinup->cella_z;
        xdist2 = xnewup[jx][jy2][jz]-xnewup[jx][jy1][jz]+ytran*gridinup->cellb_x;
        ydist2 = ynewup[jx][jy2][jz]-ynewup[jx][jy1][jz]+ytran*gridinup->cellb_y;
        zdist2 = znewup[jx][jy2][jz]-znewup[jx][jy1][jz]+ytran*gridinup->cellb_z;
        xdist3 = xnewup[jx][jy][jz2]-xnewup[jx][jy][jz1]+ztran*gridinup->cellc_x;
        ydist3 = ynewup[jx][jy][jz2]-ynewup[jx][jy][jz1]+ztran*gridinup->cellc_y;
        zdist3 = znewup[jx][jy][jz2]-znewup[jx][jy][jz1]+ztran*gridinup->cellc_z;
        volnew = 0.125*(xdist1*(ydist2*zdist3-ydist3*zdist2)+xdist2*
          (ydist3*zdist1-ydist1*zdist3)+xdist3*(ydist1*zdist2-ydist2*zdist1));
        voltotup += volnew;
        core.grid[jx][jy][jz] = volnew;
        /* divide by old volume and subtract old energy to replace their values,
           cf. the numerator of Eq 10 in JCTC 2014 */
        endelta = pothi.grid[jx][jy][jz]*(volnew/gridinup->volvox)-gridoutup->grid[jx][jy][jz];
        gridoutup->grid[jx][jy][jz] += endelta;
        encoreup -= endelta;
        xdist1 = xnewdn[jx2][jy][jz]-xnewdn[jx1][jy][jz]+xtran*gridindn->cella_x;
        ydist1 = ynewdn[jx2][jy][jz]-ynewdn[jx1][jy][jz]+xtran*gridindn->cella_y;
        zdist1 = znewdn[jx2][jy][jz]-znewdn[jx1][jy][jz]+xtran*gridindn->cella_z;
        xdist2 = xnewdn[jx][jy2][jz]-xnewdn[jx][jy1][jz]+ytran*gridindn->cellb_x;
        ydist2 = ynewdn[jx][jy2][jz]-ynewdn[jx][jy1][jz]+ytran*gridindn->cellb_y;
        zdist2 = znewdn[jx][jy2][jz]-znewdn[jx][jy1][jz]+ytran*gridindn->cellb_z;
        xdist3 = xnewdn[jx][jy][jz2]-xnewdn[jx][jy][jz1]+ztran*gridindn->cellc_x;
        ydist3 = ynewdn[jx][jy][jz2]-ynewdn[jx][jy][jz1]+ztran*gridindn->cellc_y;
        zdist3 = znewdn[jx][jy][jz2]-znewdn[jx][jy][jz1]+ztran*gridindn->cellc_z;
        volnew = 0.125*(xdist1*(ydist2*zdist3-ydist3*zdist2)+xdist2*
          (ydist3*zdist1-ydist1*zdist3)+xdist3*(ydist1*zdist2-ydist2*zdist1));
        voltotdn += volnew;
        temp.grid[jx][jy][jz] = volnew;
        endelta = potlo.grid[jx][jy][jz]*(volnew/gridindn->volvox)-gridoutdn->grid[jx][jy][jz];
        gridoutdn->grid[jx][jy][jz] += endelta;
        encoredn -= endelta;
      }
    }
  }

  fprintf(cplog, "Interpolated unit cell vol: %.6e, %.6e\n", voltotup, voltotdn);
  if (fabs(gridinup->volcell-voltotup)>1.0e-5 || fabs(gridindn->volcell-voltotdn)>1.0e-5) {
    printf("\n  CAUTION: Interpolated unit cell volume has changed a lot! Continuing anyway...\n");
    fprintf(cplog, "WARNING: interpolated unit cell volumes: %f and %f\n", voltotup, voltotdn);
    fprintf(cplog, "                      real cell volumes: %f and %f\n", gridinup->volcell, gridindn->volcell);
    errcount++;
  }
  fprintf(cplog, "\n");
  /* record total change in average map value */
  en_core[0] = encoreup;
  en_core[2] = encoredn;
  FixEdges(gridoutup);
  FixEdges(gridoutdn);
  for (i=0; i<gridin->nion; i++) fprintf(cplog, "Hirshfeld charge on atom %d: %.6e\n", i+1, hmap.chg[i]);
  fprintf(cplog, "\n");
  if (printhmap==1) {  /* output voxel map of new volumes */
    FixEdges(&core);
    snprintf(str, STRMAX, "%s-upvox.xsf", cpoutname);
    fptr = fopen(str, "w");
    OutputXSF(fptr, gridin, &core);
    fclose(fptr);
    FixEdges(&temp);
    snprintf(str, STRMAX, "%s-dnvox.xsf", cpoutname);
    fptr = fopen(str, "w");
    OutputXSF(fptr, gridin, &temp);
    fclose(fptr);
  }
  return 0;
}


/* XC FUNCTIONS */

int ReadGradient() {
  /* called by: IdentifyXC */
  /* calls: Den2XSF */
  fprintf(cplog, "\n");
  Den2XSF(abinitname, dshi, "GDEN1", &gdenhi1, &gdenhi4);
  Den2XSF(abinitname, dshi, "GDEN2", &gdenhi2, &gdenhi5);
  Den2XSF(abinitname, dshi, "GDEN3", &gdenhi3, &gdenhi6);
  Den2XSF(abinitname, dseq, "GDEN1", &gdeneq1, &gdeneq4);
  Den2XSF(abinitname, dseq, "GDEN2", &gdeneq2, &gdeneq5);
  Den2XSF(abinitname, dseq, "GDEN3", &gdeneq3, &gdeneq6);
  Den2XSF(abinitname, dslo, "GDEN1", &gdenlo1, &gdenlo4);
  Den2XSF(abinitname, dslo, "GDEN2", &gdenlo2, &gdenlo5);
  Den2XSF(abinitname, dslo, "GDEN3", &gdenlo3, &gdenlo6);
  if (errcount!=0) return 1;
  return 0;
}

int IdentifyXC(int id) {
  /* called by: main */
  /* calls: ReadGradient */
  /* list of one-to-one equivalent Abinit to LibXC functional IDs */
  int abinit2libxc[43] = { 0,20,10,0,0,0,0,0,0,0,0,101130,101,0,0,0,
    161,162,0,0,0,0,0,0,0,0,163,164,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
  //int check=0, i=0, isgga=0;
  //xc_func_type func;
  xcpot[0] = 0;
  xcpot[1] = 0;
  if (id==0) {
    printf("\n  BAD NEWS: No exchange-correlation potential detected!\n", abinitout);
    fprintf(cplog, "\nTerminated because ixc=0\n");
    fprintf(cplog, "Suggestion: re-run Abinit with a different ixc or turn off exch-corr mapping\n");
    errcount++;
    return 1;
  } else if (id>42) {
    printf("\n  BAD NEWS: Invalid exchange-correlation potential detected!\n", abinitout);
    fprintf(cplog, "\nTerminated because ixc=%d > 42 is not defined by Abinit\n", ixc);
    fprintf(cplog, "Suggestion: check Abinit files for discrepencies\n");
    errcount++;
    return 2;
  } else if (id>0 && abinit2libxc[id]==0) {
    printf("\n  BAD NEWS: Exchange-correlation potential not supported!\n", abinitout);
    fprintf(cplog, "\nTerminated because ixc=%d is not supported by CPpackage\n", ixc);
    fprintf(cplog, "Suggestion: re-run Abinit with a different ixc or contact the Fredrickson Group\n");
    errcount++;
    return 3;
  } else {
    if (id<0) {  /* not a built-in Abinit functional */
      xcpot[0] = -1*id/1000;  /* remove lowest three digits */
      xcpot[1] = -1*id-(xcpot[0]*1000);  /* remove highest three digits */
      return 6;  /* libxc not enabled */
    } else {
      xcpot[0] = abinit2libxc[id]/1000;
      xcpot[1] = abinit2libxc[id]-(xcpot[0]*1000);
    }
    /*
    for (i=0; i<2; i++) {
      if (i==0 && xcpot[0]==0) continue;   only three digits
      if (nspin==2) check = xc_func_init(&func, xcpot[i], XC_POLARIZED);
      else check = xc_func_init(&func, xcpot[i], XC_POLARIZED);
      if (check!=0) {
        printf("\n  BAD NEWS: Unrecognized XC functional!\n");
        fprintf(cplog, "\nTerminated because LibXC does not support the exchange-correlation functional: %d\n", xcpot[i]);
        fprintf(cplog, "Sugggestion: re-run Abinit with different ixc input\n");
        errcount++;
        return 4;
      }
      isgga = xc_family_from_id(xcpot[i], NULL, NULL);
      switch (isgga) {
        case XC_FAMILY_LDA:
          xc_func_end(&func);
          break;   do not need electron density gradient
        default:
          check = ReadGradient();
          if (check!=0) return 5;
          xc_func_end(&func);
          goto bottom;   initialize gradient files only once
      }
    }
    bottom:
    */
    return 0;
  }
}

int CoreCorrection(struct CrystData * gridin, struct CrystData * coreout) {
  /* called by: MapEntot */
  /* calls: FinishLine, ShiftGrid, FixEdges */
  char pspfile[STRMAX], str[200];
  int check=0, i=0, jx=0, jy=0, jz=0, jx1=0, jy1=0, jz1=0, k=0, stop=0;
  int lloc=0, lmax_this=0, mmax=0, pspcod=0, pspxc=0, r2well=0;
  int neps_atomvalues=0, zint=0, zused[120];
  double core_elec=0.0, fchrg=0.0, fraction=0.0, rchrg=0.0;
  double core_r[500][120], core_rho[500][120], ddrho=0.0, drho=0.0, rho_core=0.0;
  double d1=0.0, d2=0.0, delta=0.0, distsq=0.0, rmax[120], z=0.0, zion[NIONMAX];
  double xf=0.0, yf=0.0, zf=0.0, voxcenter_x=0.0, voxcenter_y=0.0, voxcenter_z=0.0;
  FILE * fptr1;
  FILE * fptr2;
  for (i=0; i<120; i++) zused[i] = 0;
  fptr1 = fopen(abinitout, "r");
  if (fptr1==NULL) {
    printf("\n  BAD NEWS: File %s not found!\n", abinitout);
    fprintf(cplog, "\nTerminated because file %s not found\n", abinitout);
    fprintf(cplog, "Suggestion: check your Abinit files or input options\n");
    errcount++;
    return 1;
  }
  while (stop==0) {
    check = fscanf(fptr1, "%s", str);
    if (check==EOF) stop = 1;
    if (strncmp(str, "pspini:", 7)==0) {
      fscanf(fptr1, "%s %s %s %s %s %s %s", str, str, str, str, str, str, pspfile);
      FinishLine(fptr1);
      FinishLine(fptr1);
      FinishLine(fptr1);
      fscanf(fptr1, "%s", str);
      fscanf(fptr1, "%lf %lf", &z, &zion[neps_atomvalues]);
      zint = (int)z;
      FinishLine(fptr1);
      fscanf(fptr1, "%d %d %d %d %d %d", &pspcod, &pspxc, &lmax_this, &lloc, &mmax, &r2well);
      FinishLine(fptr1);
      fscanf(fptr1, "%lf %lf", &rchrg, &fchrg);
      rmax[zint] = rchrg;
      FinishLine(fptr1);
      if (fchrg!=0.0) {
        fprintf(cplog, "Core correction detected on Z=%d pseudopotential\n", zint);
        fprintf(cplog, "Reading core density from %s, lmax=%d, mmax=%d\n", pspfile, lmax_this, mmax);
        zused[zint] = 1;
        fptr2 = fopen(pspfile, "r");
        if (fptr2==NULL) {
          printf("\n  BAD NEWS: File %s not found!\n", pspfile);
          fprintf(cplog, "\nTerminated because file %s not found\n", pspfile);
          fprintf(cplog, "Suggestion: check your Abinit files or input options\n");
          errcount++;
          return 2;
        }
        for (i=0; i<18; i++) FinishLine(fptr2);
        for (i=0; i<lmax_this+1; i++) {
          FinishLine(fptr2);
          for (k=0; k<mmax; k++) FinishLine(fptr2);
        }
        for (i=1; i<mmax+1; i++) {
          fscanf(fptr2, "%lf %lf %lf %lf", &core_r[i][zint], &core_rho[i][zint], &drho, &ddrho);
        }
        fclose(fptr2);
        core_r[0][zint] = 0.0;
        core_rho[0][zint] = core_rho[0][zint];
        fprintf(cplog, "r_max = %.6e angstrom\n", rmax[zint]*R_BOHR);
      }
    }
  }
  fclose(fptr1);
  coreout->volvox = gridin->volvox;
  ShiftGrid(gridin, 0.0, coreout);  /* copies gridin grid onto coreout grid */
  for (i=0; i<gridin->nion; i++) {
    zint = gridin->zatomic[i];
    if (zused[zint]==1) {
      for (jz=-(ngz+1)/2; jz<ngz+(ngz+1)/2; jz++) {
        zf = (double)jz/(double)ngz;
        for (jy=-(ngy+1)/2; jy<ngy+(ngy+1)/2; jy++) {
          yf = (double)jy/(double)ngy;
          for (jx=-(ngx+1)/2; jx<ngx+(ngx+1)/2; jx++) {
            xf = (double)jx/(double)ngx;
            voxcenter_x = (xf)*gridin->cella_x+(yf)*gridin->cellb_x+(zf)*gridin->cellc_x;
            voxcenter_y = (xf)*gridin->cella_y+(yf)*gridin->cellb_y+(zf)*gridin->cellc_y;
            voxcenter_z = (xf)*gridin->cella_z+(yf)*gridin->cellb_z+(zf)*gridin->cellc_z;
            distsq = sqrt((voxcenter_x-gridin->xcart[i])*(voxcenter_x-gridin->xcart[i])+
              (voxcenter_y-gridin->ycart[i])*(voxcenter_y-gridin->ycart[i])+
              (voxcenter_z-gridin->zcart[i])*(voxcenter_z-gridin->zcart[i]));
            if (distsq<rmax[zint]) {
              jx1 = (jx+ngx)%ngx;
              jy1 = (jy+ngy)%ngy;
              jz1 = (jz+ngz)%ngz;
              if (distsq==0.0) rho_core = core_rho[0][zint];
              else {
                stop = 0;
                k = 0;
                while (stop==0) {
                  k++;
                  if (core_r[k][zint]>distsq) {
                    /* linear interpolation between two closest mesh points */
                    d2 = core_r[k][zint];
                    d1 = core_r[k-1][zint];
                    delta = d2-d1;
                    fraction = (distsq-d1)/delta;
                    rho_core = fraction*core_rho[k][zint]+(1.0-fraction)*core_rho[k-1][zint];
                    stop = 1;
                  }
                }
              }
              coreout->grid[jx1][jy1][jz1] += rho_core/(4.0*PI);
              core_elec += rho_core*coreout->volvox/(4.0*PI);
            }
          }
        }
      }
    }
  }
  FixEdges(coreout);
  fprintf(cplog, "   Added to core: %20.14f\n", core_elec);
  return 0;
}

int CalcVxc1(struct CrystData * gridin, struct CrystData * gridout) {
  /* called by: CalcVxc */
  /* calls: FixEdges */
  int jx=0, jy=0, jz=0;
  double denom=0.0, numer=0.0, r_s=0.0;
  static const double A0 = 0.4581652932831429;
  static const double A1 = 2.217058676663745;
  static const double A2 = 0.7405551735357053;
  static const double A3 = 0.01968227878617998;
  static const double B1 = 1.0;
  static const double B2 = 4.504130959426697;
  static const double B3 = 1.110667363742916;
  static const double B4 = 0.02359291751427506;
  for(jz=0; jz<ngz; jz++) {
    for(jy=0; jy<ngy; jy++) {
      for(jx=0; jx<ngx; jx++) {
        /* 4/3*pi*r_s^3 = 1/rho; r_s = (3/(4*pi*rho))^1/3 */
        r_s = pow(3.0/(4.0*PI*gridin->grid[jx][jy][jz]), ONETHIRD);
        numer = A0+A1*r_s+A2*pow(r_s, 2)+A3*pow(r_s, 3);
        denom = B1*r_s+B2*pow(r_s, 2)+B3*pow(r_s, 3)+B4*pow(r_s, 4);
        gridout->grid[jx][jy][jz] = -numer/denom;
      }
    }
  }
  FixEdges(gridout);
  return 0;
}

int CalcVxc2(struct CrystData * gridin, struct CrystData * gridin2, struct CrystData * gridout) {
  /* called by: CalcVxc */
  /* calls: FixEdges */
  int jx=0, jy=0, jz=0;
  double denom=0.0, numer=0.0, r_s=0.0;
  double A0,A1,A2,A3,B1,B2,B3,B4,zeta,f_zeta;
  static const double A0o = 0.4581652932831429;
  static const double A1o = 2.217058676663745;
  static const double A2o = 0.7405551735357053;
  static const double A3o = 0.01968227878617998;
  static const double B1o = 1.0;
  static const double B2o = 4.504130959426697;
  static const double B3o = 1.110667363742916;
  static const double B4o = 0.02359291751427506;
  static const double dA0 = 0.119086804055547;
  static const double dA1 = 0.6157402568883345;
  static const double dA2 = 0.1574201515892867;
  static const double dA3 = 0.003532336663397157;
  static const double dB1 = 0.000000000000;
  static const double dB2 = 0.2673612973836267;
  static const double dB3 = 0.2052004607777787;
  static const double dB4 = 0.004200005045691381;
  for(jz=0; jz<ngz; jz++) {
    for(jy=0; jy<ngy; jy++) {
      for(jx=0; jx<ngx; jx++) {
        /* 4/3*pi*r_s^3 = 1/rho; r_s = (3/(4*pi*rho))^1/3 */
        r_s = pow(3.0/(4.0*PI*(gridin->grid[jx][jy][jz]+gridin2->grid[jx][jy][jz])), ONETHIRD);
        zeta = (gridin->grid[jx][jy][jz]-gridin2->grid[jx][jy][jz])/(gridin->grid[jx][jy][jz]+gridin2->grid[jx][jy][jz]);
        f_zeta = (pow(1.0+zeta,4.0/3.0)+pow(1.0-zeta,4.0/3.0)-2.0)/(2.0*(pow(2.0,1.0/3.0)-1));
        A0 = A0o+f_zeta*dA0;
        A1 = A1o+f_zeta*dA1;
        A2 = A2o+f_zeta*dA2;
        A3 = A3o+f_zeta*dA3;
        B1 = B1o+f_zeta*dB1;
        B2 = B2o+f_zeta*dB2;
        B3 = B3o+f_zeta*dB3;
        B4 = B4o+f_zeta*dB4;
        numer = A0+A1*r_s+A2*pow(r_s, 2)+A3*pow(r_s, 3);
        denom = B1*r_s+B2*pow(r_s, 2)+B3*pow(r_s, 3)+B4*pow(r_s, 4);
        gridout->grid[jx][jy][jz] = -numer/denom;
      }
    }
  }
  FixEdges(gridout);
  return 0;
}

int CalcVxc(struct CrystData * gridin, struct CrystData * gridin2, int dsnum, struct CrystData * gridout, struct CrystData *gridout2) {
  /* called by: main, MapEntot */
  /* calls: FixEdges, CalcVxc1 */
  int check=0, i=0, id=0, jx=0, jy=0, jz=0;
  double grad[1], gradx=0.0, grady=0.0, gradz=0.0, energy[1], sden[2];
  //xc_func_type func;
  for (jz=0; jz<ngz; jz++) {
    for (jy=0; jy<ngy; jy++) {
      for (jx=0; jx<ngx; jx++) {
        gridout->grid[jx][jy][jz] = 0.0;
        if (nspin==2) gridout2->grid[jx][jy][jz] = 0.0;
      }
    }
  }
  if (nspin==1 && xcpot[0]==0 && xcpot[1]==20) {  /* default Teter93 LDA functional */
    CalcVxc1(gridin, gridout);
    return 0;
  }
  if (nspin==2 && xcpot[0]==0 && xcpot[1]==20) {  /* default Teter93 LDA functional */
    CalcVxc2(gridin, gridin2, gridout);
    return 0;
  }
  printf("Appropriate XC functional not found.\n");
  /*
  for (i=0; i<2; i++) {
    if (i==0 && xcpot[0]==0) continue;   only three digits
    else if (nspin==2) xc_func_init(&func, xcpot[i], XC_POLARIZED);
    else xc_func_init(&func, xcpot[i], XC_UNPOLARIZED);
    id = xc_family_from_id(xcpot[i], NULL, NULL);
    switch (id) {
      case XC_FAMILY_LDA:
        for (jz=0; jz<ngz; jz++) {
          for (jy=0; jy<ngy; jy++) {
            for (jx=0; jx<ngx; jx++) {
              if (nspin==2) {
                sden[0] = gridin2->grid[jx][jy][jz];   alpha density
                sden[1] = gridin->grid[jx][jy][jz];   beta density
                xc_lda_exc(&func, 1, sden, energy);
                gridout2->grid[jx][jy][jz] += energy[0];
              } else {
                xc_lda_exc(&func, 1, &gridin->grid[jx][jy][jz], energy);
                gridout->grid[jx][jy][jz] += energy[0];
              }
            }
          }
        }
        break;
      case XC_FAMILY_GGA:
      case XC_FAMILY_HYB_GGA:
        for (jz=0; jz<ngz; jz++) {
          for (jy=0; jy<ngy; jy++) {
            for (jx=0; jx<ngx; jx++) {
              if(dsnum==dshi) {
                if (nspin==2) {
                  sden[0] = gridin2->grid[jx][jy][jz];
                  sden[1] = gridin->grid[jx][jy][jz];
                  gradx = gdenhi4.grid[jx][jy][jz];
                  grady = gdenhi5.grid[jx][jy][jz];
                  gradz = gdenhi6.grid[jx][jy][jz];
                  grad[0] = gradx*gradx+grady*grady+gradz*gradz;
                  xc_gga_exc(&func, 1, sden, &grad[0], &energy[0]);
                  gridout2->grid[jx][jy][jz] += energy[0];
                } else {
                  gradx = gdenhi1.grid[jx][jy][jz];
                  grady = gdenhi2.grid[jx][jy][jz];
                  gradz = gdenhi3.grid[jx][jy][jz];
                  grad[0] = gradx*gradx+grady*grady+gradz*gradz;
                  xc_gga_exc(&func, 1, &gridin->grid[jx][jy][jz], &grad[0], &energy[0]);
                  gridout->grid[jx][jy][jz] += energy[0];
                }
              }
              else if(dsnum==dseq) {
                if (nspin==2) {
                  sden[0] = gridin2->grid[jx][jy][jz];
                  sden[1] = gridin->grid[jx][jy][jz];
                  gradx = gdeneq4.grid[jx][jy][jz];
                  grady = gdeneq5.grid[jx][jy][jz];
                  gradz = gdeneq6.grid[jx][jy][jz];
                  grad[0] = gradx*gradx+grady*grady+gradz*gradz;
                  xc_gga_exc(&func, 1, sden, &grad[0], &energy[0]);
                  gridout2->grid[jx][jy][jz] += energy[0];
                } else {
                  gradx = gdeneq1.grid[jx][jy][jz];
                  grady = gdeneq2.grid[jx][jy][jz];
                  gradz = gdeneq3.grid[jx][jy][jz];
                  grad[0] = gradx*gradx+grady*grady+gradz*gradz;
                  xc_gga_exc(&func, 1, &gridin->grid[jx][jy][jz], &grad[0], &energy[0]);
                  gridout->grid[jx][jy][jz] += energy[0];
                }
              }
              else if(dsnum==dslo) {
                if (nspin==2) {
                  sden[0] = gridin2->grid[jx][jy][jz];
                  sden[1] = gridin->grid[jx][jy][jz];
                  gradx = gdenlo4.grid[jx][jy][jz];
                  grady = gdenlo5.grid[jx][jy][jz];
                  gradz = gdenlo6.grid[jx][jy][jz];
                  grad[0] = gradx*gradx+grady*grady+gradz*gradz;
                  xc_gga_exc(&func, 1, sden, &grad[0], &energy[0]);
                  gridout2->grid[jx][jy][jz] += energy[0];
                } else {
                  gradx = gdenlo1.grid[jx][jy][jz];
                  grady = gdenlo2.grid[jx][jy][jz];
                  gradz = gdenlo3.grid[jx][jy][jz];
                  grad[0] = gradx*gradx+grady*grady+gradz*gradz;
                  xc_gga_exc(&func, 1, &gridin->grid[jx][jy][jz], &grad[0], &energy[0]);
                  gridout->grid[jx][jy][jz] += energy[0];
                }
              }   grad[0] = square of vector magnitude of electron density gradient
            }
          }
        }
        break;
      default:
        printf("\n  BAD NEWS: Unrecognized XC functional!\n");
        fprintf(cplog, "\nTerminated because the XC potential was not recognized as an LDA or GGA functional\n");
        fprintf(cplog, "Sugggestion: re-run Abinit with a different XC functional or check that ixc=%d is correct\n", xcpot[i]);
        errcount++;
        return 1;
    }
    xc_func_end(&func);
  }
  FixEdges(gridout);
  */
  return 2;
}


/* MAP CALCULATION FUNCTIONS */

int ReadEwald(int dsnum) {
  /* called by: main */
  /* calls: ReadLine */
  int check = 0, i=0, stop=0;
  char line[500], str1[30], str2[30], str3[30];
  double E_core[dsnum], E_entropy[dsnum], E_hartree[dsnum], E_int[dsnum];
  double E_kinetic[dsnum], E_locPsp[dsnum], E_nlPsp[dsnum], E_total[dsnum], E_xc[dsnum];
  FILE * fptr;
  fptr = fopen(abinitout, "r");
  if (fptr==NULL) {
    printf("\n  BAD NEWS: File %s not found!\n", abinitout);
    fprintf(cplog, "\nTerminated because file %s not found\n", abinitout);
    fprintf(cplog, "Suggestion: check your Abinit files or input options\n");
    errcount++;
    return 1;
  }
  while (stop==0) {
    check = ReadLine(fptr, line);
    if (check==1) {
      printf("\n  BAD NEWS: Energy data in %s not found!\n", abinitout);
      fprintf(cplog, "\nTerminated because total free energy list #%d not found in %s\n", i+1, abinitout);
      fprintf(cplog, "Suggestion: check %s or check if data matches format in CPpackage source code\n", abinitout);
      errcount++;
      return 2;
    }
    if (strncmp(line, " Components of total free energy (in Hartree) :", 47)==0) {
      fscanf(fptr, "%s %s %s %lf", str1, str2, str3, &E_kinetic[i]);
      fscanf(fptr, "%s %s %s %lf", str1, str2, str3, &E_hartree[i]);
      fscanf(fptr, "%s %s %s %lf", str1, str2, str3, &E_xc[i]);
      fscanf(fptr, "%s %s %s %lf", str1, str2, str3, &E_ewald[i]);
      fscanf(fptr, "%s %s %s %lf", str1, str2, str3, &E_core[i]);
      fscanf(fptr, "%s %s %s %lf", str1, str2, str3, &E_locPsp[i]);
      fscanf(fptr, "%s %s %s %lf", str1, str2, str3, &E_nlPsp[i]);
      fscanf(fptr, "%s %s %s %lf", str1, str2, str3, &E_int[i]);
      fscanf(fptr, "%s %s %lf" ,str1, str2, &E_entropy[i]);
      fscanf(fptr, "%s %s %lf" ,str1, str2, &E_total[i]);
      i++;
    }
    if (i==dsnum) stop = 1;
  }
  fclose(fptr);

  return 0;
}


int CalcKineticTF(struct CrystData * gridin, int dsnum, struct CrystData * gridout) {
  /* called by: MapEntot */
  /* calls: ReadLine */
  char line[500], str1[30], str2[30], str3[30];
  int check=0, i=0, jx=0, jy=0, jz=0, stop=0;
  double E_core[dsnum], E_entropy[dsnum], E_ewald[dsnum], E_hartree[dsnum], E_int[dsnum];
  double E_kinetic[dsnum], E_locPsp[dsnum], E_nlPsp[dsnum], E_total[dsnum], E_xc[dsnum];
  double normKE=0;
  FILE * fptr;
  fptr = fopen(abinitout, "r");
  if (fptr==NULL) {
    printf("\n  BAD NEWS: File %s not found!\n", abinitout);
    fprintf(cplog, "\nTerminated because file %s not found\n", abinitout);
    fprintf(cplog, "Suggestion: check your Abinit files or input options\n");
    errcount++;
    return 1;
  }
  while (stop==0) {
    check = ReadLine(fptr, line);
    if (check==1) {
      printf("\n  BAD NEWS: Energy data in %s not found!\n", abinitout);
      fprintf(cplog, "\nTerminated because total free energy list #%d not found in %s\n", i+1, abinitout);
      fprintf(cplog, "Suggestion: check %s or check if data matches format in CPpackage source code\n", abinitout);
      errcount++;
      return 2;
    }
    if (strncmp(line, " Components of total free energy (in Hartree) :", 47)==0) {
      fscanf(fptr, "%s %s %s %lf", str1, str2, str3, &E_kinetic[i]);
      fscanf(fptr, "%s %s %s %lf", str1, str2, str3, &E_hartree[i]);
      fscanf(fptr, "%s %s %s %lf", str1, str2, str3, &E_xc[i]);
      fscanf(fptr, "%s %s %s %lf", str1, str2, str3, &E_ewald[i]);
      fscanf(fptr, "%s %s %s %lf", str1, str2, str3, &E_core[i]);
      fscanf(fptr, "%s %s %s %lf", str1, str2, str3, &E_locPsp[i]);
      fscanf(fptr, "%s %s %s %lf", str1, str2, str3, &E_nlPsp[i]);
      fscanf(fptr, "%s %s %s %lf", str1, str2, str3, &E_int[i]);
      fscanf(fptr, "%s %s %lf" ,str1, str2, &E_entropy[i]);
      fscanf(fptr, "%s %s %lf" ,str1, str2, &E_total[i]);
      i++;
    }
    if (i==dsnum) stop = 1;
  }
  fclose(fptr);
  for (jz=0; jz<ngz; jz++) {
    for (jy=0; jy<ngy; jy++) {
      for (jx=0; jx<ngx; jx++) {
        normKE += pow(gridin->grid[jx][jy][jz], 5.0/3.0)*gridin->volvox;
      }
    }
  }
  for (jz=0; jz<=ngz; jz++) {
    for (jy=0; jy<=ngy; jy++) {
      for (jx=0; jx<=ngx; jx++) {
        gridout->grid[jx][jy][jz] = E_kinetic[dsnum-1]*(pow(gridin->grid[jx][jy][jz], 5.0/3.0)/normKE);
      }
    }
  }
  return 0;
}

int MapEntot(struct CrystData * denin, struct CrystData * denin2, struct CrystData * kdenin, struct CrystData * kdenin2, struct CrystData * ldenin, struct CrystData * ldenin2, struct CrystData * potin, struct CrystData * potin2, struct CrystData * vhxcin, struct CrystData * vhxcin2, struct CrystData * vhain, struct CrystData * vhain2, int dsnum, double vol, struct CrystData * etotout, struct CrystData * locout) {
  /* called by: main */
  /* calls: AddGrid, CalcKineticTF, CalcVxc, CoreCorrection, IntegrateGrid, MultiplyGrid, OutputXSF, ScaleGrid, SubtractGrid, CopyStruct */
  char enfile[STRMAX];
  int check=0;
  FILE * fptr;
  core.volvox = vol;
  core2.volvox = vol;
  temp.volvox = vol;
  temp2.volvox = vol;
  vxc.volvox = vol;
  vxc2.volvox = vol;

  CopyStruct(denin,locout);

  fprintf(cplog, " Cell volume DS%d: %20.14f\n", dsnum, denin->volcell);
  if (nspin==2) fprintf(cplog, "  Electron count: %20.14f\n", IntegrateGrid(denin)+IntegrateGrid(denin2));
  else fprintf(cplog, "  Electron count: %20.14f\n", IntegrateGrid(denin));
  if (mapkin==2) {
    kdenin->volvox = vol;
    if (nspin==2) {
      AddGrid(denin, denin2, &temp);
      check = CalcKineticTF(&temp, dsnum, kdenin);
    } else check = CalcKineticTF(denin, dsnum, kdenin);
    if (check!=0) return 1;
  }
  if (mapkin==1 || mapkin==2) {
    if (mapkinoption==2) {
      ScaleGrid(ldenin,0.25,ldenin);
      SubtractGrid(kdenin,ldenin,kdenin);
    }
    if (nspin==2 && mapkin!=2) p_kin[dsnum-1] = IntegrateGrid(kdenin)+IntegrateGrid(kdenin2);
    else p_kin[dsnum-1] = IntegrateGrid(kdenin);
    fprintf(cplog, "  Kinetic energy: %20.14f\n", p_kin[dsnum-1]);
  }
  SubtractGrid(potin, vhxcin, potin);
  if (nspin==2) {
    SubtractGrid(potin2, vhxcin2, potin2);
    AddGrid(potin, potin2, potin);  /* combine up and down pot from here on */
  }
  if (maploc!=1) {
    SubtractGrid(potin, potin, potin);
    if (nspin==2) SubtractGrid(potin2, potin2, potin2);
  }
  if (nspin==2) {
    AddGrid(denin, denin2, &temp2);
    MultiplyGrid(potin, &temp2, locout);
  }
  else {
    MultiplyGrid(potin, denin, locout);  /* temp is vden here */
  }
  p_loc[dsnum-1] = IntegrateGrid(locout);
  fprintf(cplog, "  V_local energy: %20.14f\n", p_loc[dsnum-1]);

  if (mapxc==1) {  /* Vxc = epsilon(r) where Exc[rho] = integ(epsilon(r)*rho(r))dV */
    check = CoreCorrection(denin, &core);  /* for psuedopotentials with nonlinear core correction */
    printf("here5!\n");
    if (check!=0) return 2;
    if (nspin==2) {
      check = CoreCorrection(denin2, &core2);
      if (check!=0) return 2;
      printf("here6!\n");
    }
    check = CalcVxc(&core, &core2, dsnum, &vxc, &vxc2);
    if (check!=0) return 3;
    printf("here7!\n");
  } else {
    SubtractGrid(&vxc, &vxc, &vxc);
    if (nspin==2) SubtractGrid(&vxc2, &vxc2, &vxc2);
  }
  printf("here4!\n");
  if (maphart==1) {
    ScaleGrid(vhain, 0.5, vhain);
    AddGrid(potin, vhain, potin);
    if (nspin==2) {
      ScaleGrid(vhain2, 0.5, vhain2);
      AddGrid(potin, vhain2, potin);
    }
  } else {
    SubtractGrid(vhain, vhain, vhain);
    if (nspin==2) SubtractGrid(vhain2, vhain2, vhain2);
  }
  if (nspin==2) {
    AddGrid(denin, denin2, &temp2);
    MultiplyGrid(vhain2, &temp2, &temp);  /* vhain1 is empty */
    printf("Hartree: %lf\n",IntegrateGrid(&temp));
  } else MultiplyGrid(vhain, denin, &temp);
  p_hart[dsnum-1] = IntegrateGrid(&temp);
  fprintf(cplog, "  Hartree energy: %20.14f\n", p_hart[dsnum-1]);

  if (nspin==2) {
    AddGrid(&core, &core2, &temp2);
    MultiplyGrid(&temp2, &vxc, &temp);
  } else MultiplyGrid(&core, &vxc, &temp);  /* core = denin for HGH pseudopotentials */
  p_xc[dsnum-1] = IntegrateGrid(&temp);
  fprintf(cplog, "Exch-corr energy: %20.14f\n", p_xc[dsnum-1]);
  if (p_xc[dsnum-1]!=p_xc[dsnum-1]) {  /* undefined number */
    printf("\n  BAD NEWS: The electron density is sometimes negative!\n");
    fprintf(cplog, "\nTerminated because of negative values in the electron density file\n");
    fprintf(cplog, "Suggestion: re-run Abinit with higher ecut or less empty space in the unit cell\n");
    errcount++;
    return 3;
  }

  if (nspin==2) {
    AddGrid(denin, denin2, &temp2);
    MultiplyGrid(potin, &temp2, potin);
  } else MultiplyGrid(potin, denin, potin);  /* Vden(r) = V(r)*rho(r) potential energy density */
  AddGrid(potin, &temp, potin);  /* full potential including exchange-correlation */
  fprintf(cplog, "Potential energy: %20.14f\n", IntegrateGrid(potin));
  if (maploc==1 && IntegrateGrid(potin)>0.0) {
    printf("\n  BAD NEWS: The potential energy is positive!\n");
    fprintf(cplog, "\nTerminated because total potential energy from DS%d is positive\n", dsnum);
    fprintf(cplog, "Suggestion: check your Abinit files for inconsistencies\n");
    errcount++;
    return 4;
  }

  if (mapkin==1 || mapkin==2) AddGrid(kdenin, potin, potin);  /* total energy density */
  if (nspin==2 && mapkin==1) AddGrid(kdenin2, potin, potin);
  ScaleGrid(potin, denin->volvox, etotout);  /* etot = grid of total voxel energies */
  fprintf(cplog, "Total map energy: %20.14f\n\n", IntegrateGrid(etotout));
  if (printen==1) {
    snprintf(enfile, STRMAX, "%s-emap%d.xsf", cpoutname, dsnum);
    fptr = fopen(enfile, "w");
    OutputXSF(fptr, denin, etotout);
    fclose(fptr);
  }
  return 0;
}


int CalcEwald(struct CrystData * gridin) {
  /* called by: main */
  /* calls: none */
  /* see: Martin, R.M. Electronic Structure: Basic Theory and Practical Methods Eqn. F.5 */
  int i=0, j=0, ka=0, kb=0, kc=0, h=0, k=0, l=0;
  double xc1=0.0, xc2=0.0, yc1=0.0, yc2=0.0, zc1=0.0, zc2=0.0;
  double dist=0.0, normG=0.0;
  double ax_star=0.0, ay_star=0.0, az_star=0.0, bx_star=0.0, by_star=0.0, bz_star=0.0,
         cx_star=0.0, cy_star=0.0, cz_star=0.0;
  double normastar=0.0, normbstar=0.0, normcstar=0.0, eta=0.0;
  double Gx=0.0, Gy=0.0, Gz=0.0, delta_rx=0.0, delta_ry=0.0, delta_rz=0.0,shiftx,shifty,shiftz;
  double temp1=0.0, temp2=0.0;

  ax_star = 2*PI*(gridin->cellb_y*gridin->cellc_z-gridin->cellc_y*gridin->cellb_z)/gridin->volcell;
  ay_star = -2*PI*(gridin->cellb_x*gridin->cellc_z-gridin->cellc_x*gridin->cellb_z)/gridin->volcell;
  az_star = 2*PI*(gridin->cellb_x*gridin->cellc_y-gridin->cellc_x*gridin->cellb_y)/gridin->volcell;
  bx_star = 2*PI*(gridin->cellc_y*gridin->cella_z-gridin->cella_y*gridin->cellc_z)/gridin->volcell;
  by_star = -2*PI*(gridin->cellc_x*gridin->cella_z-gridin->cella_x*gridin->cellc_z)/gridin->volcell;
  bz_star = 2*PI*(gridin->cellc_x*gridin->cella_y-gridin->cella_x*gridin->cellc_y)/gridin->volcell;
  cx_star = 2*PI*(gridin->cella_y*gridin->cellb_z-gridin->cellb_y*gridin->cella_z)/gridin->volcell;
  cy_star = -2*PI*(gridin->cella_x*gridin->cellb_z-gridin->cellb_x*gridin->cella_z)/gridin->volcell;
  cz_star = 2*PI*(gridin->cella_x*gridin->cellb_y-gridin->cellb_x*gridin->cella_y)/gridin->volcell;
  normastar = sqrt(ax_star*ax_star+ay_star*ay_star+az_star*az_star);
  normbstar = sqrt(bx_star*bx_star+by_star*by_star+bz_star*bz_star);
  normcstar = sqrt(cx_star*cx_star+cy_star*cy_star+cz_star*cz_star);
  /* eta is a factor intended to speed convergence */
  eta = 10000.0;
  if (normastar<eta) {
    eta = normastar;
  }
  if (normbstar<eta) {
    eta = normbstar;
  }
  if (normcstar<eta) {
    eta = normcstar;
  }

  for (i=0; i<gridin->nion; i++) {
    Ewald[i] = 0.0;
    Ewald_vo[i] = 0.0;
    Ewald_sc[i] = 0.0;
  }
  /* see Equations XX-XX in XXXX, 2018 for semicore (loc) and valence (itin) Ewald energy */
  for (i=0; i<gridin->nion; i++) {
    xc1 = gridin->xcart[i];
    yc1 = gridin->ycart[i];
    zc1 = gridin->zcart[i];
    for (j=0; j<gridin->nion; j++) {
      xc2 = gridin->xcart[j];
      yc2 = gridin->ycart[j];
      zc2 = gridin->zcart[j];
      delta_rx = xc2-xc1;
      delta_ry = yc2-yc1;
      delta_rz = zc2-zc1;
      for (ka=-8; ka<9; ka++) {
        for (kb=-8; kb<9; kb++) {
          for (kc=-8; kc<9; kc++) {
            shiftx = ka*gridin->cella_x+kb*gridin->cellb_x+kc*gridin->cellc_x;
            shifty = ka*gridin->cella_y+kb*gridin->cellb_y+kc*gridin->cellc_y;
            shiftz = ka*gridin->cella_z+kb*gridin->cellb_z+kc*gridin->cellc_z;
            dist = sqrt((delta_rx-shiftx)*(delta_rx-shiftx)+(delta_ry-shifty)*(delta_ry-shifty)+(delta_rz-shiftz)*(delta_rz-shiftz));
            if (dist>0.0) {
              Ewald_vo[i] += (0.5*vo_elec[j]*vo_elec[i]*gsl_sf_erfc(eta*dist))/dist;
              Ewald_sc[i] += (sc_elec[i]*vo_elec[j]*gsl_sf_erfc(eta*dist))/dist;
              Ewald_sc[i] += (0.5*sc_elec[i]*sc_elec[j]*gsl_sf_erfc(eta*dist))/dist;
            }
          }
        }
      }
      /* these limits are excessive, consider cutting down if taking forever */
      for (h=-20; h<21; h++ ){
        for (k=-20; k<21; k++) {
          for (l=-20; l<21; l++) {
            Gx = h*ax_star+k*bx_star+l*cx_star;
            Gy = h*ay_star+k*by_star+l*cy_star;
            Gz = h*az_star+k*bz_star+l*cz_star;
            /* G = h*astar + k*bstar + l*cstar */
            normG = sqrt(Gx*Gx+Gy*Gy+Gz*Gz);
            if (normG>0.0) {
              temp1 = -1.0*(pow(normG,2)/(4*pow(eta,2)));
              temp2 = Gx*delta_rx + Gy*delta_ry + Gz*delta_rz;
              /* separate vo and sc sums of Ewald energy */
              /* vo_elec and sc_elec based on user inputs when radial density profiles entered if semicore pseudopotentials are involved */
              Ewald_vo[i] += (0.5*(4*PI/gridin->volcell)*vo_elec[i]*vo_elec[j]*exp(temp1)*cos(temp2))/(pow(normG,2));
              Ewald_sc[i] += ((4*PI/gridin->volcell)*sc_elec[i]*vo_elec[j]*exp(temp1)*cos(temp2))/(pow(normG,2));
              Ewald_sc[i] += (0.5*(4*PI/gridin->volcell)*sc_elec[i]*sc_elec[j]*exp(temp1)*cos(temp2))/(pow(normG,2));
            }
          }
        }
      }
    }
    Ewald_vo[i] -= (0.5*vo_elec[i]*vo_elec[i]*2*eta)/(pow(PI,0.5));
    Ewald_sc[i] -= (sc_elec[i]*vo_elec[i]*2*eta)/(pow(PI,0.5));
    Ewald_sc[i] -= (0.5*sc_elec[i]*sc_elec[i]*2*eta)/(pow(PI,0.5));
    for (j=0; j<gridin->nion; j++) {
      Ewald_vo[i] -= (0.5*vo_elec[i]*vo_elec[j]*PI)/(pow(eta,2)*gridin->volcell);
      Ewald_sc[i] -= (sc_elec[i]*vo_elec[j]*PI)/(pow(eta,2)*gridin->volcell);
      Ewald_sc[i] -= (0.5*sc_elec[i]*sc_elec[j]*PI)/(pow(eta,2)*gridin->volcell);
    }
  }
  for (i=0; i<gridin->nion; i++) {
    Ewald[i] = Ewald_sc[i]+Ewald_vo[i];
    printf("Ewald_vo[%d] = %lf   ",i,Ewald_vo[i]);
    printf("Ewald_sc[%d] = %lf   ",i,Ewald_sc[i]);
    printf("Ewald_tot[%d] = %lf\n",i,Ewald_sc[i]+Ewald_vo[i]);
  }
  return 0;
}

int CalcEalpha(struct CrystData * gridin) {
  /* called by: main */
  /* calls: ElementName */
  char element[STRMAX], src[100], str1[100];
  int i=0, j=0, check=0, stop=0, ntypat;
  double epsatm_bytype[NIONMAX],znucl_bytype[NIONMAX], temp0=0.0;
  FILE * fptr2;
  /* See equations XX in XXXX, XXXX */
  fptr2 = fopen(abinitout, "r");
  if (fptr2==NULL) {
    printf("\n  BAD NEWS: File %s not found!\n", abinitout);
    fprintf(cplog, "\nTerminated because file %s not found\n", abinitout);
    fprintf(cplog, "Suggestion: check your Abinit files or input options\n");
    errcount++;
    return 1;
  }
  strcpy(src, "znucl");
  stop = 0;
  ntypat=0;
  while (stop==0) {
    check=fscanf(fptr2, "%s",str1);
    if(check == EOF) {
      printf("\n  BAD NEWS: %s value not found in file %s!\n", src,abinitout);
      fprintf(cplog, "\nTerminated because %s value not found in file %s.d\n",src, abinitout);
      fprintf(cplog, "Suggestion: check your Abinit files or input options\n");
      errcount++;
      return 1;
    }
    if (strcmp(str1, src)==0) {
      while (stop==0) {
        check=fscanf(fptr2, "%lf",&znucl_bytype[ntypat]);
        if(check!=0) {
          ntypat++;
        }
        else {
          stop =1;
        }
      }
    }
  }
  strcpy(src, "epsatm=");
  for(j=0;j<ntypat;j++) {
    stop = 0;
    while (stop==0) {
      check=fscanf(fptr2, "%s",str1);
      if(check == EOF) {
        printf("\n  BAD NEWS: %s value not found in file %s!\n", src,abinitout);
        fprintf(cplog, "\nTerminated because %s value not found in file %s.d\n",src, abinitout);
        fprintf(cplog, "Suggestion: check your Abinit files or input options\n");
        errcount++;
        return 1;
      }
      if (strcmp(str1, src)==0) {
        fscanf(fptr2, "%lf",&temp0);
        epsatm_bytype[j] = temp0;
        stop=1;
      }
    }
  }
  fclose(fptr2);
  for (i=0; i<gridin->nion; i++) {
    for(j=0; j<ntypat; j++) {
      if(znucl_bytype[j]==gridin->zatomic[i]) {
        epsatm[i] = epsatm_bytype[j];
        ElementName(gridin->zatomic[i],element);
      }
    }
  }
  for (i=0; i<gridin->nion; i++) {
    E_alpha_vo[i]=0.0;
    E_alpha_sc[i]=0.0;
    E_alpha[i]=0.0;
    for (j=0; j<gridin->nion; j++) {
      E_alpha_vo[i] += vo_elec[i]*epsatm[j]/gridin->volcell;
      E_alpha_sc[i] += sc_elec[i]*epsatm[j]/gridin->volcell;
    }
  }
  for (i=0; i<gridin->nion; i++) {
    E_alpha[i] = E_alpha_vo[i]+E_alpha_sc[i];
    fprintf(cplog, "E_alpha_vo[%d] = %lf E_alpha_sc[%d] = %lf\n",i,E_alpha_vo[i],i,E_alpha_sc[i]);
  }
  return 0;
}

int MapNonloc(struct CrystData * gridin, struct CrystData * nonlocin, struct CrystData * nonlocin2, struct CrystData * gridout, int dsnum) {
  /* called by: main */
  /* calls: SymmetrizeGrid, ReadXSF */
  char filename[STRMAX];
  int check=0, jx=0, jy=0, jz=0;
  /* XSF files containing spatially mapped nonlocal energies from CPnonlocal are read in and added to the appropriate energy grids */
  snprintf(filename, STRMAX, "%s_DS%d_NL.xsf", abinitname, dsnum);
  ReadXSF(filename, nonlocin);
  check = SymmetrizeGrid(nonlocin, &smap);
  for (jx=0; jx<ngx; jx++) {
    for (jy=0; jy<ngy; jy++) {
      for (jz=0; jz<ngz; jz++) {
        gridout->grid[jx][jy][jz] += nonlocin->grid[jx][jy][jz]*gridin->volvox;
        p_nonloc[dsnum-1] += nonlocin->grid[jx][jy][jz]*gridin->volvox;
      }
    }
  }
  return 0;
}

int MapNonloc2(struct CrystData * gridin, struct CrystData * nonlocin, struct CrystData * nonlocin2, struct CrystData * gridout, int dsnum) {
  /* called by: main */
  /* calls: SymmetrizeGrid, ReadXSF */
  char filename[STRMAX];
  int check=0, jx=0, jy=0, jz=0;
  /* XSF files containing spatially mapped nonlocal energies from CPnonlocal are read in and added to the appropriate energy grids */
  snprintf(filename, STRMAX, "%s_DS%d_NL.xsf", abinitname, dsnum);
  ReadXSF(filename, nonlocin);
  check = SymmetrizeGrid(nonlocin, &smap);
  for (jx=0; jx<ngx; jx++) {
    for (jy=0; jy<ngy; jy++) {
      for (jz=0; jz<ngz; jz++) {
        p_nonloc[dsnum-1] += nonlocin->grid[jx][jy][jz]*gridin->volvox;
      }
    }
  }
  return 0;
}


int MapEwald(struct CrystData * gridin, struct CrystData * grideq, struct CrystData * gridout, struct CrystData * gridoutewald, struct CrystData * gridoutalpha, struct CrystData * nonlocin, int dsnum, double E_ewald, double E_ewald2, struct CrystData * gridsc) {
  /* called by: main */
  /* calls: IntegrateGrid, CopyStruct, OutputXSF */
  FILE * fptr;
  char ewaldfile[STRMAX], ealphafile[STRMAX];
  int i=0, jx=0, jy=0, jz=0, atom=0, aoi=0, nmax_atom=0, nvox[NIONMAX];
  int ngcount=0, ngp=0, ngp0=0, voxtot=ngx*ngy*ngz;
  double sum_E_den, mean_E_den, temp_ealpha;
  double Ewald_frac[NIONMAX], Ewald_bkgd_sum[NIONMAX], Alpha_bkgd_sum[NIONMAX], wsum[NIONMAX], Ewald_bkgd, Alpha_bkgd, sc_den_sum[NIONMAX];
  double hirshmax = 0.0, temp=0.0, atom_ewald[NIONMAX], den_sum[NIONMAX];
  double E_Ewald_atom_local[NIONMAX], E_alpha_atom_local[NIONMAX], E_Ewald_homogeneous, E_alpha_homogeneous, nonloc[NIONMAX], Ewald_tot=0.0, Alpha_tot=0.0, E_alpha_tot=0.0, E_alpha_tot2=0.0;
  double alpha_ratio=0.0;
  /* set up homogeneous (energy) and bkgd (density) variables */
  E_Ewald_homogeneous = 0.0;
  Ewald_bkgd = 0.0;
  E_alpha_homogeneous = 0.0;
  Alpha_bkgd = 0.0;
  sum_E_den = IntegrateGrid(gridin)/gridout->volvox;
  mean_E_den = sum_E_den/(double)voxtot;
  /* Calculating E_alpha for the system at the particular mean electron density (expanded or contracted volume) */
  for (i=0; i<gridin->nion; i++) {
    E_alpha_tot2 += epsatm[i]*mean_E_den;
    E_alpha_tot += E_alpha[i];
  }
  alpha_ratio=0.0;
  alpha_ratio = E_alpha_tot2/E_alpha_tot;
  /* based on Zion(itin) and Zion(loc) split E_Ewald and E_alpha into homogeneous and localized components */
  for (i=0; i<gridin->nion; i++) {
    Ewald_frac[i] = Ewald[i]/E_ewald2;
    E_Ewald_by_atom[i][dsnum-1] = Ewald_frac[i]*E_ewald;
    E_Ewald_atom_local[i] = (Ewald_sc[i]/E_ewald2)*E_ewald;
    fprintf(cplog, "Local Ewald energy for atom %d: %lf\n",i,E_Ewald_atom_local[i]);
    E_Ewald_homogeneous += (Ewald_vo[i]/E_ewald2)*E_ewald;
    Ewald_tot += E_Ewald_atom_local[i] + (Ewald_vo[i]/E_ewald2)*E_ewald;
    if (E_alpha_option==1) {
      E_alpha_by_atom[i][dsnum-1] = E_alpha[i]*alpha_ratio;
      E_alpha_atom_local[i] = E_alpha_sc[i]*alpha_ratio;
      E_alpha_homogeneous += (E_alpha_vo[i]*alpha_ratio);
      Alpha_tot += E_alpha_atom_local[i] + (E_alpha_vo[i]*alpha_ratio);
    }
  }
  fprintf(cplog, "Homogeneous Ewald energy: %lf\n\n",E_Ewald_homogeneous);
  Ewald_bkgd = E_Ewald_homogeneous/voxtot;
  if (E_alpha_option==1) {
    Alpha_bkgd = (E_alpha_homogeneous)/voxtot;
  }
  /* add homogeneous portions of Ewald and E_alpha to energy maps */
  for (jx=0; jx<ngx; jx++) {
    for (jy=0; jy<ngy; jy++) {
      for (jz=0; jz<ngz; jz++) {
        gridout->grid[jx][jy][jz]+=Alpha_bkgd;
        gridout->grid[jx][jy][jz]+=Ewald_bkgd;
        gridoutewald->grid[jx][jy][jz] += Ewald_bkgd;
        gridoutalpha->grid[jx][jy][jz] += Alpha_bkgd;
      }
    }
  }
  gridoutewald->volvox = gridout->volvox;
  gridoutalpha->volvox = gridout->volvox;
/* divide into atomic volumes to map localized components */
  if(Ewald_map_done==0) {
    for (jx=0; jx<ngx; jx++) {
      for (jy=0; jy<ngy; jy++) {
        for (jz=0; jz<ngz; jz++) {
          ngcount++;
          ngp = ngcount*100/voxtot;
          if (ngp!=ngp0) {
            printf("\r%d%%", ngp);
            fflush(stdout);
          }
          ngp0=ngp;
          nmax_atom = 0;
          hirshmax = 0.0;
          /* calculate hirshfeld weights for atoms at this voxel */
          /* find hirshmax for voxel jx jy jz */
          for (atom=0; atom<gridin->nion; atom++) {
            if (hmap.hirsh_weight[atom][jx][jy][jz] > hirshmax){
              hirshmax = hmap.hirsh_weight[atom][jx][jy][jz];
            }
            else continue;
          }
          /* add voxel to atom with hirshmax */
          for (atom=0; atom<gridin->nion; atom++) {
            if ((hirshmax - hmap.hirsh_weight[atom][jx][jy][jz])<tolerance) {
              nvox[atom]++;
              nmax_atom++;
            }
          }
          /* assign binary hirshfeld weights */
          for (atom=0; atom<gridin->nion; atom++) {
            if ((hirshmax - hmap.hirsh_weight[atom][jx][jy][jz])<tolerance) {
              hmap.hirsh_weight[atom][jx][jy][jz] = (1/(double)nmax_atom);
            }
            else {
              hmap.hirsh_weight[atom][jx][jy][jz] = 0.0;
            }
          }
        }
      }
    }
    Ewald_map_done=1;
  }
  /* variables to be used in mapping localized Ewald components */
  for (atom = 0; atom<gridin->nion; atom++) {
    wsum[atom] = 0.0;
    atom_ewald[atom] = 0.0;
    den_sum[atom] = 0.0;
    Ewald_bkgd_sum[atom] = 0.0;
    Alpha_bkgd_sum[atom] = 0.0;
    nonloc[atom] = 0.0;
    sc_den_sum[atom] = 0.0;
  }
  for (jx=0; jx<ngx; jx++) {
    for (jy=0; jy<ngy; jy++) {
      for (jz=0; jz<ngz; jz++) {
        for (atom = 0; atom<gridin->nion; atom++) {
          wsum[atom] += hmap.hirsh_weight[atom][jx][jy][jz];
          den_sum[atom] += hmap.hirsh_weight[atom][jx][jy][jz]*grideq->grid[jx][jy][jz];
          Ewald_bkgd_sum[atom] +=hmap.hirsh_weight[atom][jx][jy][jz]*Ewald_bkgd;
          Alpha_bkgd_sum[atom] +=hmap.hirsh_weight[atom][jx][jy][jz]*Alpha_bkgd;
          if (E_Ewald_option==2) {
            sc_den_sum[atom] += hmap.hirsh_weight[atom][jx][jy][jz]*gridsc->grid[jx][jy][jz];
          }
          if (E_Ewald_option==3) {
            nonloc[atom] += hmap.hirsh_weight[atom][jx][jy][jz]*fabs(nonlocin->grid[jx][jy][jz]);
          }
          if(mapnonloc==2) {
            nonloc[atom] += hmap.hirsh_weight[atom][jx][jy][jz]*nonlocin->grid[jx][jy][jz]*nonlocin->volvox;
          }
        }
      }
    }
  }
  printf("  Ewald Background comparisons:\n");
  for (atom = 0; atom<gridin->nion; atom++) {
     printf("    Atom %d.  Background sum:  %13.8lf,  E_Ewald_vo:  %13.8lf\n",atom+1, Ewald_bkgd_sum[atom], (Ewald_vo[atom]/E_ewald2)*E_ewald);
  }

  /* allocate localized Ewald/E_alpha energy in proportion to selected variable */
  for (jx=0; jx<ngx; jx++) {
    for (jy=0; jy<ngy; jy++) {
      for (jz=0; jz<ngz; jz++) {
        for (atom=0; atom<gridin->nion; atom++) {
          if (hmap.hirsh_weight[atom][jx][jy][jz] > 0.0) {
            aoi = atom;
            if (E_Ewald_option == 1) { /* map localized Ewald in proportion to total density */
              temp = E_Ewald_atom_local[aoi]*hmap.hirsh_weight[aoi][jx][jy][jz]*grideq->grid[jx][jy][jz]/den_sum[aoi];
              gridout->grid[jx][jy][jz] +=(temp);
              gridoutewald->grid[jx][jy][jz] += temp;
              if (E_alpha_option == 0) { /* homogeneous E_alpha, no localized component */
                temp_ealpha = 0.0;
              }
              if (E_alpha_option == 1) { /* localized E_alpha mapped */
                temp_ealpha = E_alpha_atom_local[aoi]*hmap.hirsh_weight[aoi][jx][jy][jz]*grideq->grid[jx][jy][jz]/den_sum[aoi];
                gridout->grid[jx][jy][jz] +=temp_ealpha;
                gridoutalpha->grid[jx][jy][jz] += temp_ealpha;
              }
              if(mapnonloc==2) {
                temp = nonloc[aoi]*hmap.hirsh_weight[aoi][jx][jy][jz]*grideq->grid[jx][jy][jz]/den_sum[aoi];
                gridout->grid[jx][jy][jz] +=temp;
              }
            }

            if (E_Ewald_option == 2) { /* map localized Ewald in proportion to semicore electron density */
              temp = E_Ewald_atom_local[aoi]*hmap.hirsh_weight[aoi][jx][jy][jz]*gridsc->grid[jx][jy][jz]/sc_den_sum[aoi];
              gridout->grid[jx][jy][jz] +=(temp);
              gridoutewald->grid[jx][jy][jz] += temp;
              if (E_alpha_option == 0) {
                temp_ealpha = 0.0;
              }
              if (E_alpha_option == 1) {
                temp_ealpha = E_alpha_atom_local[aoi]*hmap.hirsh_weight[aoi][jx][jy][jz]*gridsc->grid[jx][jy][jz]/sc_den_sum[aoi];
                gridout->grid[jx][jy][jz] +=temp_ealpha;
                gridoutalpha->grid[jx][jy][jz] += temp_ealpha;
              }
            }
           if (E_Ewald_option == 3) { /* map localized Ewald in proportion to nonlocal energy */
              temp = E_Ewald_atom_local[aoi]*hmap.hirsh_weight[aoi][jx][jy][jz]*fabs(nonlocin->grid[jx][jy][jz])/nonloc[aoi];
              gridout->grid[jx][jy][jz] +=(temp);
              gridoutewald->grid[jx][jy][jz] += temp;
              if (E_alpha_option == 0) {
                temp_ealpha = 0.0;
              }
              if (E_alpha_option == 1) {
                temp_ealpha = E_alpha_atom_local[aoi]*hmap.hirsh_weight[aoi][jx][jy][jz]*fabs(nonlocin->grid[jx][jy][jz])/nonloc[aoi];
                gridout->grid[jx][jy][jz] +=temp_ealpha;
                gridoutalpha->grid[jx][jy][jz] += temp_ealpha;
              }
            }
            if (E_Ewald_option == 4) { /* map localized Ewald homogeneously within atomic volume */
              temp = E_Ewald_atom_local[aoi]*hmap.hirsh_weight[aoi][jx][jy][jz]/wsum[aoi];
              gridout->grid[jx][jy][jz] +=(temp);
              gridoutewald->grid[jx][jy][jz] += temp;
              if (E_alpha_option == 0) {
                temp_ealpha = 0.0;
              }
              if (E_alpha_option == 1) {
                temp_ealpha = E_alpha_atom_local[aoi]*hmap.hirsh_weight[aoi][jx][jy][jz]/wsum[aoi];
                gridout->grid[jx][jy][jz] +=temp_ealpha;
                gridoutalpha->grid[jx][jy][jz] += temp_ealpha;
              }
            }
            if (E_alpha_option == 0) {
              temp_ealpha = 0.0;
            }
            atom_ewald[aoi] += temp;
          }
        }
      }
    }
  }
  p_ewald[dsnum-1] = IntegrateGrid(gridoutewald)/(gridoutewald->volvox);
  p_ealpha[dsnum-1] = IntegrateGrid(gridoutalpha)/(gridoutewald->volvox);
  fprintf(cplog, "Ewald energy apportioned to atoms:\n");
  for (i=0; i<gridin->nion; i++) {
    fprintf(cplog,"Atom %i: %lf\n",i,atom_ewald[i]);
  }
  if (outputewaldalpha == 1) {
    snprintf(ewaldfile, STRMAX, "%s-ewald%i.xsf", cpoutname,dsnum);
    CopyStruct(gridin,gridoutewald);
    fptr = fopen(ewaldfile, "w");
    OutputXSF(fptr, gridoutewald, gridoutewald);
    fclose(fptr);
    if (E_alpha_option == 1) {
      snprintf(ealphafile, STRMAX, "%s-ealpha%i.xsf", cpoutname, dsnum);
      CopyStruct(gridin, gridoutalpha);
      fptr = fopen(ealphafile, "w");
      OutputXSF(fptr, gridoutalpha, gridoutalpha);
      fclose(fptr);
    }
  }
  return 0;
}


int MapEwald_Bader(struct CrystData * gridin, struct CrystData * grideq, struct CrystData * gridout, struct CrystData * gridoutewald, struct CrystData * gridoutalpha, struct CrystData * nonlocin, int dsnum, double E_ewald, double E_ewald2, struct CrystData * gridsc,  struct ContactVol * map) {
  /* called by: main */
  /* calls: IntegrateGrid, CopyStruct, OutputXSF */
  FILE * fptr;
  char ewaldfile[STRMAX], ealphafile[STRMAX];
  int j,i=0, jx=0, jy=0, jz=0, atom=0, aoi=0, nmax_atom=0, nvox[NIONMAX];
  int ngcount=0, ngp=0, ngp0=0, voxtot=ngx*ngy*ngz;
  double sum_E_den, mean_E_den, temp_ealpha;
  double Ewald_frac[NIONMAX], Ewald_bkgd_sum[NIONMAX], Alpha_bkgd_sum[NIONMAX], wsum[NIONMAX], Ewald_bkgd, Alpha_bkgd, sc_den_sum[NIONMAX];
  double hirshmax = 0.0, temp=0.0, atom_ewald[NIONMAX], den_sum[NIONMAX];
  double E_Ewald_atom_local[NIONMAX], E_alpha_atom_local[NIONMAX], E_Ewald_homogeneous, E_alpha_homogeneous, nonloc[NIONMAX], Ewald_tot=0.0, Alpha_tot=0.0, E_alpha_tot=0.0, E_alpha_tot2=0.0;
  double alpha_ratio=0.0;
  /* set up homogeneous (energy) and bkgd (density) variables */
  E_Ewald_homogeneous = 0.0;
  Ewald_bkgd = 0.0;
  E_alpha_homogeneous = 0.0;
  Alpha_bkgd = 0.0;
  sum_E_den = IntegrateGrid(gridin)/gridout->volvox;
  mean_E_den = sum_E_den/(double)voxtot;
  /* Calculating E_alpha for the system at the particular mean electron density (expanded or contracted volume) */
  for (i=0; i<gridin->nion; i++) {
    E_alpha_tot2 += epsatm[i]*mean_E_den;
    E_alpha_tot += E_alpha[i];
  }
  alpha_ratio=0.0;
  alpha_ratio = E_alpha_tot2/E_alpha_tot;
  /* based on Zion(itin) and Zion(loc) split E_Ewald and E_alpha into homogeneous and localized components */
  for (i=0; i<gridin->nion; i++) {
    Ewald_frac[i] = Ewald[i]/E_ewald2;
    E_Ewald_by_atom[i][dsnum-1] = Ewald_frac[i]*E_ewald;
    E_Ewald_atom_local[i] = (Ewald_sc[i]/E_ewald2)*E_ewald;
    fprintf(cplog, "Local Ewald energy for atom %d: %lf\n",i,E_Ewald_atom_local[i]);
    E_Ewald_homogeneous += (Ewald_vo[i]/E_ewald2)*E_ewald;
    Ewald_tot += E_Ewald_atom_local[i] + (Ewald_vo[i]/E_ewald2)*E_ewald;
    if ((E_alpha_option==1)||(E_Ewald_option==5)) {
      E_alpha_by_atom[i][dsnum-1] = E_alpha[i]*alpha_ratio;
      E_alpha_atom_local[i] = E_alpha_sc[i]*alpha_ratio;
      E_alpha_homogeneous += (E_alpha_vo[i]*alpha_ratio);
      Alpha_tot += E_alpha_atom_local[i] + (E_alpha_vo[i]*alpha_ratio);
    }
    printf("  Alpha components for atom %d:  %lf + %lf \n",i,E_alpha_atom_local[i],E_alpha_vo[i]*alpha_ratio);
  }
  printf("   --------------------------------  %lf\n",Alpha_tot);
  fprintf(cplog, "Homogeneous Ewald energy: %lf\n\n",E_Ewald_homogeneous);
  Ewald_bkgd = E_Ewald_homogeneous/voxtot;
  if (E_alpha_option==1) {
    Alpha_bkgd = (E_alpha_homogeneous)/voxtot;
  }
  if(E_Ewald_option==5) {
         printf("  Entering beta-stage Ewald_option 5 mode...\n");
         Alpha_bkgd = (E_alpha_homogeneous)/voxtot;
         E_alpha_option=1;
  }
  /* add homogeneous portions of Ewald and E_alpha to energy maps */
  for (jx=0; jx<ngx; jx++) {
    for (jy=0; jy<ngy; jy++) {
      for (jz=0; jz<ngz; jz++) {
        gridoutewald->grid[jx][jy][jz] = 0.0;
        gridoutalpha->grid[jx][jy][jz] = 0.0;
      }
    }
  }
  for (jx=0; jx<ngx; jx++) {
    for (jy=0; jy<ngy; jy++) {
      for (jz=0; jz<ngz; jz++) {
        gridout->grid[jx][jy][jz]+=Alpha_bkgd;
        gridout->grid[jx][jy][jz]+=Ewald_bkgd;
        gridoutewald->grid[jx][jy][jz] += Ewald_bkgd;
        gridoutalpha->grid[jx][jy][jz] += Alpha_bkgd;
      }
    }
  }
  gridoutewald->volvox = gridout->volvox;
  gridoutalpha->volvox = gridout->volvox;
/* divide into atomic volumes to map localized components */
  for (atom = 0; atom<gridin->nion; atom++) {
    wsum[atom] = 0.0;
    atom_ewald[atom] = 0.0;
    den_sum[atom] = 0.0;
    Ewald_bkgd_sum[atom] = 0.0;
    Alpha_bkgd_sum[atom] = 0.0;
    nonloc[atom] = 0.0;
    sc_den_sum[atom] = 0.0;
  }
  for (jx=0; jx<ngx; jx++) {
    for (jy=0; jy<ngy; jy++) {
      for (jz=0; jz<ngz; jz++) {
        for(j=0;j<map->neighcount2[jx][jy][jz];j++) {
          atom = map->ionmap[j][jx][jy][jz]&127;
          wsum[atom] += 1.0/map->swj[jx][jy][jz];
          den_sum[atom] += (1.0/map->swj[jx][jy][jz])*grideq->grid[jx][jy][jz];
          Ewald_bkgd_sum[atom] +=Ewald_bkgd/map->swj[jx][jy][jz];
          Alpha_bkgd_sum[atom] +=(1.0/map->swj[jx][jy][jz])*Alpha_bkgd;
          if (E_Ewald_option==2) {
            sc_den_sum[atom] += (1.0/map->swj[jx][jy][jz])*gridsc->grid[jx][jy][jz];
          }
          if (E_Ewald_option==3) {
            nonloc[atom] +=(1.0/map->swj[jx][jy][jz])*fabs(nonlocin->grid[jx][jy][jz]);
          }
          if(mapnonloc==2) {
            nonloc[atom] += (1.0/map->swj[jx][jy][jz])*nonlocin->grid[jx][jy][jz]*nonlocin->volvox;
          }
        }
      }
    }
  }
  printf("  Ewald Background comparisons:\n");
  for (atom = 0; atom<gridin->nion; atom++) {
     printf("    Atom %d.  Ewald background sum:  %13.8lf,  E_Ewald_vo:  %13.8lf   Alpha background sum:  %13.8lf,  alpha_vo:  %13.8lf\n",atom+1, Ewald_bkgd_sum[atom], (Ewald_vo[atom]/E_ewald2)*E_ewald, Alpha_bkgd_sum[atom],E_alpha[atom]*alpha_ratio);
  }
  /* allocate localized Ewald/E_alpha energy in proportion to selected variable */
  for (jx=0; jx<ngx; jx++) {
    for (jy=0; jy<ngy; jy++) {
      for (jz=0; jz<ngz; jz++) {
        for(j=0;j<map->neighcount2[jx][jy][jz];j++) {
            atom = map->ionmap[j][jx][jy][jz]&127;
            aoi = atom;
            if (E_Ewald_option == 1) { /* map localized Ewald in proportion to total density */
              temp = E_Ewald_atom_local[aoi]*(1.0/map->swj[jx][jy][jz])*grideq->grid[jx][jy][jz]/den_sum[aoi];
              gridout->grid[jx][jy][jz] +=(temp);
              gridoutewald->grid[jx][jy][jz] += temp;
              if (E_alpha_option == 0) { /* homogeneous E_alpha, no localized component */
                temp_ealpha = 0.0;
              }
              if (E_alpha_option == 1) { /* localized E_alpha mapped */
                temp_ealpha = E_alpha_atom_local[aoi]*(1.0/map->swj[jx][jy][jz])*grideq->grid[jx][jy][jz]/den_sum[aoi];
                gridout->grid[jx][jy][jz] +=temp_ealpha;
                gridoutalpha->grid[jx][jy][jz] += temp_ealpha;
              }
              if(mapnonloc==2) {
                temp = nonloc[aoi]*(1.0/map->swj[jx][jy][jz])*grideq->grid[jx][jy][jz]/den_sum[aoi];
                gridout->grid[jx][jy][jz] +=temp;
              }
            }
            if (E_Ewald_option == 5) { /* map localized Ewald in proportion to total density */
//              if(Ewald_bkgd_sum[atom] > (Ewald_vo[atom]/E_ewald2)*E_ewald) {
                temp = (((Ewald_vo[atom]/E_ewald2)*E_ewald)-Ewald_bkgd_sum[atom])*(1.0/map->swj[jx][jy][jz])*grideq->grid[jx][jy][jz]/den_sum[aoi];
                gridout->grid[jx][jy][jz] +=(temp);
                gridoutewald->grid[jx][jy][jz] += temp;
                temp_ealpha = ((E_alpha[atom]*alpha_ratio)-Alpha_bkgd_sum[atom])*(1.0/map->swj[jx][jy][jz])*grideq->grid[jx][jy][jz]/den_sum[aoi];
                gridout->grid[jx][jy][jz] +=temp_ealpha;
                gridoutalpha->grid[jx][jy][jz] += temp_ealpha;
//              }
/*              else {
                temp = (((Ewald_vo[atom]/E_ewald2)*E_ewald)-Ewald_bkgd_sum[atom])*(1.0/map->swj[jx][jy][jz])/wsum[aoi];
                gridout->grid[jx][jy][jz] +=(temp);
                gridoutewald->grid[jx][jy][jz] += temp;
                temp_ealpha = ((E_alpha[atom]*alpha_ratio)-Alpha_bkgd_sum[atom])*(1.0/map->swj[jx][jy][jz])/wsum[aoi];
                gridout->grid[jx][jy][jz] +=temp_ealpha;
                gridoutalpha->grid[jx][jy][jz] += temp_ealpha;
              }
*/
	    }
            if (E_Ewald_option == 2) { /* map localized Ewald in proportion to semicore electron density */
              temp = E_Ewald_atom_local[aoi]*(1.0/map->swj[jx][jy][jz])*gridsc->grid[jx][jy][jz]/sc_den_sum[aoi];
              gridout->grid[jx][jy][jz] +=(temp);
              gridoutewald->grid[jx][jy][jz] += temp;
              if (E_alpha_option == 0) {
                temp_ealpha = 0.0;
              }
              if (E_alpha_option == 1) {
                temp_ealpha = E_alpha_atom_local[aoi]*(1.0/map->swj[jx][jy][jz])*gridsc->grid[jx][jy][jz]/sc_den_sum[aoi];
                gridout->grid[jx][jy][jz] +=temp_ealpha;
                gridoutalpha->grid[jx][jy][jz] += temp_ealpha;
              }
            }
           if (E_Ewald_option == 3) { /* map localized Ewald in proportion to nonlocal energy */
              temp = E_Ewald_atom_local[aoi]*(1.0/map->swj[jx][jy][jz])*fabs(nonlocin->grid[jx][jy][jz])/nonloc[aoi];
              gridout->grid[jx][jy][jz] +=(temp);
              gridoutewald->grid[jx][jy][jz] += temp;
              if (E_alpha_option == 0) {
                temp_ealpha = 0.0;
              }
              if (E_alpha_option == 1) {
                temp_ealpha = E_alpha_atom_local[aoi]*(1.0/map->swj[jx][jy][jz])*fabs(nonlocin->grid[jx][jy][jz])/nonloc[aoi];
                gridout->grid[jx][jy][jz] +=temp_ealpha;
                gridoutalpha->grid[jx][jy][jz] += temp_ealpha;
              }
            }
            if (E_Ewald_option == 4) { /* map localized Ewald homogeneously within atomic volume */
              temp = E_Ewald_atom_local[aoi]*(1.0/map->swj[jx][jy][jz])/wsum[aoi];
              gridout->grid[jx][jy][jz] +=(temp);
              gridoutewald->grid[jx][jy][jz] += temp;
              if (E_alpha_option == 0) {
                temp_ealpha = 0.0;
              }
              if (E_alpha_option == 1) {
                temp_ealpha = E_alpha_atom_local[aoi]*(1.0/map->swj[jx][jy][jz])/wsum[aoi];
                gridout->grid[jx][jy][jz] +=temp_ealpha;
                gridoutalpha->grid[jx][jy][jz] += temp_ealpha;
              }
            }
            if (E_alpha_option == 0) {
              temp_ealpha = 0.0;
            }
            atom_ewald[aoi] += temp;
        }
      }
    }
  }
  p_ewald[dsnum-1] = IntegrateGrid(gridoutewald)/(gridoutewald->volvox);
  p_ealpha[dsnum-1] = IntegrateGrid(gridoutalpha)/(gridoutewald->volvox);
  fprintf(cplog, "Ewald energy apportioned to atoms:\n");
  for (i=0; i<gridin->nion; i++) {
    fprintf(cplog,"Atom %i: %lf\n",i,atom_ewald[i]);
  }
  if (outputewaldalpha == 1) {
    snprintf(ewaldfile, STRMAX, "%s-ewald%i.xsf", cpoutname,dsnum);
    CopyStruct(gridin,gridoutewald);
    fptr = fopen(ewaldfile, "w");
    OutputXSF(fptr, gridoutewald, gridoutewald);
    fclose(fptr);
    if (E_alpha_option == 1) {
      snprintf(ealphafile, STRMAX, "%s-ealpha%i.xsf", cpoutname, dsnum);
      CopyStruct(gridin, gridoutalpha);
      fptr = fopen(ealphafile, "w");
      OutputXSF(fptr, gridoutalpha, gridoutalpha);
      fclose(fptr);
    }
  }
  return 0;
}

/* Jonathan adding code */

double autocali_get_slope(float x1, float y1, float x2, float y2) {
  double slope;
  slope = (y2 - y1) / (x2 - x1);
  return(slope);
}

double autocali_stop() {
  float percent;
  float min = 10000;
  float max = -10000;
  for (it=0; it<autocali_num_site; it++) {
    if (min > ps2[it]) {
      min = ps2[it];
    }
    if (max < ps2[it]) {
      max = ps2[it];
    }
  }
  percent = fabs(100 * (max - min) / min);
  printf("Percent difference in pressures is %f %\n", percent);
  return(percent);
}

int autocali_iteration() {
  double percent;
  percent = autocali_stop();
  autocali_min_then = ps1[autocali_min_index];
  autocali_min_now = ps2[autocali_min_index];
  if (percent > 1) {
    for (it=0; it<autocali_num_site; it++) {
      if (it != autocali_min_index) {
        autocali_diff1 = ps1[it] - autocali_min_then;
        autocali_diff2 = ps2[it] - autocali_min_now;
        if (es2[it] != 0) {
          autocali_slope = autocali_get_slope(autocali_diff1, es1[it], autocali_diff2, es2[it]);
          autocali_yint = es1[it] - autocali_slope * autocali_diff1;
          es2_temp[it] = autocali_yint;
        }
        else {
          es2_temp[it] = autocali_diff2 * 0.001;
        }
      }
      else {
        es2_temp[it] = 0.000000;
      }
    }
  }
  else {
    printf("CONVERGENCE REACHED! All pressures are within %lf percent of the smallest pressure\n", percent);
    printf("Localized electron counts are as follows:\n");
    for (it=0; it<autocali_num_site; it++) {
      printf("  Site %d: %lf localized electrons\n", it, es2_temp[it]);
    }
    exit(0);
  }
  return(0);
}

int autocali_start() {
  double min = 10000;
  autocali_num_site = deneq.nion;
  memcpy(ps1, ps2, sizeof(ps1));
  for (it=0; it<autocali_num_site; it++) {
    es1[it] = 0.0;
    es2[it] = 0.0;
  }
  for (it=0; it<autocali_num_site; it++) {
    if (min > ps2[it]) {
      min = ps2[it];
      autocali_min_index = it;
    }
  }
  return(0);
}

/* Jonathan stop adding code */



double PressureContrib(struct CrystData * gridinup, struct CrystData * gridin, struct CrystData * gridindn, double pin_kin, double pin_loc, double pin_hart, double pin_xc, double pin_ew, double pin_ealpha, double pin_nonloc, double * p_entropy, double * p_mapcore, double * pout_nonloc, double * pout_ewald, double * pout_ealpha, double * pout_kin, double * pout_loc, double * pout_hart, double * pout_xc) {
  /* called by: CalcCP */
  /* calls: ReadLine */
  char check=0, i=0, line[500], stop=0, str1[30], str2[30], str3[30];
  double dV=0.0, pin_map=0.0, pin_nonmap=0.0, pin_tot=0.0;
  double p_mapped=0.0, p_nonmapped=0.0, p_total=0.0;
  double E_core[3], E_entropy[3], E_ewald[3], E_hartree[3], E_int[3];
  double E_kinetic[3], E_locpsp[3], E_nonlocpsp[3], E_total[3], E_xc[3];
  FILE * fptr;
  fptr = fopen(abinitout, "r");
  if (fptr==NULL) {
    printf("\n  BAD NEWS: File %s not found!\n", abinitout);
    fprintf(cplog, "\nTerminated because file %s not found\n", abinitout);
    fprintf(cplog, "Suggestion: check your files and input options\n");
    errcount++;
    return 0.0;
  }
  while (stop==0) {
    check = ReadLine(fptr, line);
    if (check==1) {
      printf("\n  BAD NEWS: Energy data in %s not found!\n", abinitout);
      fprintf(cplog, "\nTerminated because total free energy list #%d not found in %s\n", i+1, abinitout);
      fprintf(cplog, "Suggestion: check %s or check if data matches format in CPpackage source code\n", abinitout);
      errcount++;
      return 0.0;
    }
    if (strncmp(line, " Components of total free energy (in Hartree) :", 47)==0) {
      fscanf(fptr, "%s %s %s %lf", str1, str2, str3, &E_kinetic[i]);
      fscanf(fptr, "%s %s %s %lf", str1, str2, str3, &E_hartree[i]);
      fscanf(fptr, "%s %s %s %lf", str1, str2, str3, &E_xc[i]);
      fscanf(fptr, "%s %s %s %lf", str1, str2, str3, &E_ewald[i]);
      fscanf(fptr, "%s %s %s %lf", str1, str2, str3, &E_core[i]);
      fscanf(fptr, "%s %s %s %lf", str1, str2, str3, &E_locpsp[i]);
      fscanf(fptr, "%s %s %s %lf", str1, str2, str3, &E_nonlocpsp[i]);
      if (occopt==3) {
        fscanf(fptr, "%s %s %s %lf", str1, str2, str3, &E_int[i]);
        fscanf(fptr, "%s %s %lf", str1, str2, &E_entropy[i]);
      }
      fscanf(fptr, "%s %s %lf", str1, str2, &E_total[i]);
      i++;
    }
    if (i==3) stop = 1;
  }
  fclose(fptr);
  dV = gridindn->volcell-gridinup->volcell;  /* negative sign */
  *pout_kin = (E_kinetic[0]-E_kinetic[2])/dV;
  *pout_hart = (E_hartree[0]-E_hartree[2])/dV;
  *pout_xc = (E_xc[0]-E_xc[2])/dV;
  *pout_ewald = (E_ewald[0]-E_ewald[2])/dV;
  *pout_ealpha = (E_core[0]-E_core[2])/dV;
  *pout_loc = (E_locpsp[0]-E_locpsp[2])/dV;
  *pout_nonloc = (E_nonlocpsp[0]-E_nonlocpsp[2])/dV;
  *p_entropy = (E_entropy[0]-E_entropy[2])/dV;
  *p_mapcore = (en_core[0]-en_core[2])/dV;
  printf("\n  -----------SUMMARY OF CP PRESSURE CONTRIBUTIONS------------\n");
  fprintf(cplog, "-----------SUMMARY OF CP PRESSURE CONTRIBUTIONS------------\n");
  printf("  Mapped Contributions                 ABINIT         CPmap\n");
  fprintf(cplog, "Mapped Contributions                 ABINIT         CPmap\n");
  if (mapkin==1 || mapkin==2) {
    printf("               kinetic E pressure  %12.8f  %12.8f\n", *pout_kin, pin_kin);
    fprintf(cplog, "             kinetic E pressure  %12.8f  %12.8f\n", *pout_kin, pin_kin);
    p_mapped += *pout_kin;
    pin_map += pin_kin;
  }
  if (maploc==1) {
    printf("             local psp E pressure  %12.8f  %12.8f\n", *pout_loc, pin_loc);
    fprintf(cplog, "           local psp E pressure  %12.8f  %12.8f\n", *pout_loc, pin_loc);
    p_mapped += *pout_loc;
    pin_map += pin_loc;
  }
  if (maphart==1) {
    printf("               Hartree E pressure  %12.8f  %12.8f\n", *pout_hart, pin_hart);
    fprintf(cplog, "             Hartree E pressure  %12.8f  %12.8f\n", *pout_hart, pin_hart);
    p_mapped += *pout_hart;
    pin_map += pin_hart;
  }
  if (mapxc==1) {
    printf("             exch-corr E pressure  %12.8f  %12.8f\n", *pout_xc, pin_xc);
    fprintf(cplog, "           exch-corr E pressure  %12.8f  %12.8f\n", *pout_xc, pin_xc);
    p_mapped += *pout_xc;
    pin_map += pin_xc;
  }
  if (mapcore==1) {
    printf("         -restore core E pressure    same as ->  %12.8f\n", -*p_mapcore);
    fprintf(cplog, "       -restore core E pressure    same as ->  %12.8f\n", -*p_mapcore);
    p_mapped += -*p_mapcore;
    pin_map += -*p_mapcore;
  }
  if (E_Ewald_option != 0) {
    printf("                 Ewald E pressure  %12.8f  %12.8f\n", *pout_ewald, pin_ew);
    fprintf(cplog, "               Ewald E pressure  %12.8f  %12.8f\n", *pout_ewald, pin_ew);
    p_mapped += *pout_ewald;
    pin_map += pin_ew;
  }
  if (E_alpha_option != 0) {
    printf("                 psp core E pressure  %12.8f  %12.8f\n", *pout_ealpha, pin_ealpha);
    fprintf(cplog, "               psp core E pressure  %12.8f  %12.8f\n", *pout_ealpha, pin_ealpha);
    p_mapped += *pout_ealpha;
    pin_map += pin_ealpha;
  }
  if (mapnonloc!=0) {
    printf("               nonlocal E pressure  %12.8f  %12.8f\n", *pout_nonloc, pin_nonloc);
    fprintf(cplog, "             nonlocal E pressure  %12.8f  %12.8f\n", *pout_nonloc, pin_nonloc);
    p_mapped += pin_nonloc;
    pin_map += pin_nonloc;
    if (nspinor==2) {
          NLremainder = *pout_nonloc - pin_nonloc;
    }
  }
 printf("  -----------------------------------------------------------\n");
  fprintf(cplog, "-----------------------------------------------------------\n");
  printf("            Total mapped pressure  %12.8f  %12.8f\n\n", p_mapped, pin_map);
  fprintf(cplog, "          Total mapped pressure  %12.8f  %12.8f\n\n", p_mapped, pin_map);
  printf("  Non-mapped Contributions\n");
  fprintf(cplog, "Non-mapped Contributions\n");
  if (mapcore == 1) {
    printf("          restore core E pressure    same as ->  %12.8f\n", *p_mapcore);
    fprintf(cplog, "        restore core E pressure    same as ->  %12.8f\n", *p_mapcore);
    p_nonmapped += *p_mapcore;
    pin_nonmap += *p_mapcore;
  }
  if (mapkin!=1 && mapkin!=2) {
    printf("               kinetic E pressure  %12.8f    <- same as\n", *pout_kin);
    fprintf(cplog, "             kinetic E pressure  %12.8f    <- same as\n", *pout_kin);
    p_nonmapped += *pout_kin;
    pin_nonmap += *pout_kin;
  }
  if (maploc!=1) {
    printf("             local psp E pressure  %12.8f    <- same as\n", *pout_loc);
    fprintf(cplog, "           local psp E pressure  %12.8f    <- same as\n", *pout_loc);
    p_nonmapped += *pout_loc;
    pin_nonmap += *pout_loc;
  }
  if (maphart!=1) {
    printf("               Hartree E pressure  %12.8f    <- same as\n", *pout_hart);
    fprintf(cplog, "             Hartree E pressure  %12.8f    <- same as\n", *pout_hart);
    p_nonmapped += *pout_hart;
    pin_nonmap += *pout_hart;
  }
  if (mapxc!=1) {
    printf("             exch-corr E pressure  %12.8f    <- same as\n", *pout_xc);
    fprintf(cplog, "           exch-corr E pressure  %12.8f    <- same as\n", *pout_xc);
    p_nonmapped += *pout_xc;
    pin_nonmap += *pout_xc;
  }
  if (E_Ewald_option==0) {
    printf("      non-mapped Ewald E pressure  %12.8f    <- same as\n", *pout_ewald);
    fprintf(cplog, "    non-mapped Ewald E pressure  %12.8f    <- same as\n", *pout_ewald);
    p_nonmapped += *pout_ewald;
    pin_nonmap += *pout_ewald;
  }
  if (E_alpha_option==0) {
    printf("              psp core E pressure  %12.8f    <- same as\n", *pout_ealpha);
    fprintf(cplog, "            psp core E pressure  %12.8f    <- same as\n", *pout_ealpha);
    p_nonmapped += *pout_ealpha;
    pin_nonmap += *pout_ealpha;
  }
  if (mapnonloc==0) {
    printf("          nonlocal psp E pressure  %12.8f    <- same as\n", *pout_nonloc);
    fprintf(cplog, "        nonlocal psp E pressure  %12.8f    <- same as\n", *pout_nonloc);
    p_nonmapped += *pout_nonloc;
    pin_nonmap += *pout_nonloc;
  }
  if (occopt==3) {
    printf("             -kt*entropy pressure  %12.8f    <- same as\n", *p_entropy);
    fprintf(cplog, "           -kt*entropy pressure  %12.8f    <- same as\n", *p_entropy);
    p_nonmapped += *p_entropy;
    pin_nonmap += *p_entropy;
  }
  if(nspinor==2) {
    printf("             nonlocal remainder from SOC %12.8f\n",NLremainder);
    p_nonmapped+=NLremainder;
    pin_nonmap+=NLremainder;
  }
  printf("  -----------------------------------------------------------\n");
  fprintf(cplog, "-----------------------------------------------------------\n");
  printf("        Total non-mapped pressure  %12.8f  %12.8f\n\n", p_nonmapped, pin_nonmap);
  fprintf(cplog, "      Total non-mapped pressure  %12.8f  %12.8f\n\n", p_nonmapped, pin_nonmap);
  p_total = (E_total[0]-E_total[2])/dV;
  pin_tot = pin_map+pin_nonmap;
  printf("  Total Pressure                   %12.8f  %12.8f\n", p_total, pin_tot);
  fprintf(cplog, "Total Pressure                   %12.8f  %12.8f\n", p_total, pin_tot);
  printf("  -----------------------------------------------------------\n");
  fprintf(cplog, "-----------------------------------------------------------\n");
  if (fabs((p_total-pin_tot))>1.0e-5) {
    printf("\n  BAD NEWS: Pressures from Abinit and CP disagree!\n");
    fprintf(cplog, "\nTerminated because pressures from Abinit and CPmap disagree\n");
    fprintf(cplog, "Suggestion: check your files for odd behavior\n");
    errcount++;
    return 0.0;
  }
  return pin_tot;
}

int CalcCP(struct CrystData * denin_up, struct CrystData * denin, struct CrystData * denin_dn, struct CrystData * etotin_up, struct CrystData * etotin_dn, struct CrystData * cpout) {
  /* called by: main */
  /* calls: CopyStruct, IntegrateGrid, PressureContrib, ScaleGrid, ShiftGrid, SubtractGrid, SymmetrizeGrid, OutputXSF */
  char cpfile[STRMAX];
  int check=0, voxtot=ngx*ngy*ngz;
  double dV=0.0, hart=0.0, kin=0.0, loc=0.0, ptot=0.0, xc=0.0, ew=0.0, alpha=0.0, nonloc=0.0;
  double p_check=0.0, p_entropy=0.0, p_mapcore=0, p_remain=0.0;
  FILE * fptr;
  CopyStruct(denin, cpout);  /* copies geometric information */
  dV = etotin_dn->volvox-etotin_up->volvox;  /* smaller minus larger because p = -dE/dV */
  SubtractGrid(etotin_up, etotin_dn, cpout);  /* cpout is energy difference here */
  cpout->volvox = (etotin_up->volvox+etotin_dn->volvox)/2.0;
  ScaleGrid(cpout, 1.0/dV, cpout);
  p_remain = IntegrateGrid(cpout)/(cpout->volvox*voxtot);
  if (mapcore==1 && mapsym==1) {
    printf("  Restoring symmetry after interpolation\n");
    check = SymmetrizeGrid(cpout, &smap);
    if (check!=0) return 1;
    p_check = IntegrateGrid(cpout)/(cpout->volvox*voxtot);
    if (fabs(p_remain-p_check)>1.0e-12) {
      fprintf(cplog, "WARNING: total pressure changed after symmetrization\n");
      fprintf(cplog, " before: %20.14f\n", p_remain);
      fprintf(cplog, "  after: %20.14f\n", p_check);
      errcount++;
      if (fabs(p_remain-p_check)<1.0e-5) printf("\n  CAUTION: Total pressure has changed! Continuing anyway...\n");
      else {
        printf("\n  BAD NEWS: Total pressure has changed!\n");
        fprintf(cplog, "\nTerminated because symmetrization caused a change in pressure\n");
        fprintf(cplog, "Suggestion: check Abinit outfile for correct nsym, symrel, and tnons\n");
        return 2;
      }
    }
  }
  fprintf(cplog, "        dV/volume: %.6e\n", -dV/cpout->volvox);
  fprintf(cplog, "CP before mapping: %.6e\n\n", p_remain);
  p_kin[3] = p_kin[0]-p_kin[2];
  p_loc[3] = p_loc[0]-p_loc[2];
  p_hart[3] = p_hart[0]-p_hart[2];
  p_xc[3] = p_xc[0]-p_xc[2];
  p_ewald[3] = p_ewald[0]-p_ewald[2];
  p_ealpha[3] = p_ealpha[0]-p_ealpha[2];
  p_nonloc[3] = p_nonloc[0]-p_nonloc[2];
  kin = p_kin[3]/(dV*voxtot);
  loc = p_loc[3]/(dV*voxtot);
  hart = p_hart[3]/(dV*voxtot);
  xc = p_xc[3]/(dV*voxtot);
  ew = p_ewald[3]/(dV*voxtot);
  alpha = p_ealpha[3]/(dV*voxtot);
  nonloc = p_nonloc[3]/(dV*voxtot);
  ptot = PressureContrib(denin_up, denin, denin_dn, kin, loc, hart, xc, ew, alpha, nonloc, &p_entropy, &p_mapcore,
    &p_nonloc[3], &p_ewald[3], &p_ealpha[3], &p_kin[3], &p_loc[3], &p_hart[3], &p_xc[3]);
  if (ptot==0.0) return 3;
  if (rescp==3 && mapcore==1) p_remain = p_mapcore;
  else if (rescp==3) p_remain = 0.0;
  else if (rescp==2) {
    p_remain = -1.0*IntegrateGrid(cpout)/(cpout->volvox*voxtot);
    fprintf(cplog, "Homogenous background pressure = %.6e\n", p_remain);
  } else if (rescp==1) {
    p_remain = 0.0;
    if (mapkin!=1 && mapkin!=2) p_remain += p_kin[3];
    if (maploc!=1) p_remain += p_loc[3];
    if (maphart!=1) p_remain += p_hart[3];
    if (mapxc!=1) p_remain += p_xc[3];
    if (E_Ewald_option==0) p_remain += p_ewald[3];
    if (E_alpha_option==0) p_remain += p_ealpha[3];
    if (mapnonloc==0) p_remain += p_nonloc[3];
    if (mapcore==1) p_remain += p_mapcore;
    if (nspinor==2) p_remain += NLremainder;
    p_remain += p_entropy;
  }
  ShiftGrid(cpout, p_remain, cpout);
  p_remain = IntegrateGrid(cpout)/(cpout->volvox*voxtot);
  printf("  Average CP of map = %.6e  au\n", p_remain);
  printf("                    = %.6e GPa\n\n", p_remain*AU2GPA);
  fprintf(cplog, "Average CP = %.6e  au\n", p_remain);
  fprintf(cplog, "           = %.6e GPa\n\n", p_remain*AU2GPA);
  if (fabs(ptot-p_remain)>1.0e-5) {
    printf("  BAD NEWS: Average mapped pressure does not match DFT results!\n");
    printf("ptot is %lf p_remain is %lf great thx \n",ptot,p_remain);
    fprintf(cplog, "\nTerminated because averaged mapped pressure does not match total in Abinit outfile\n");
    fprintf(cplog, "Suggestion: check the Abinit outfile and potentials for odd behavior\n");
    errcount++;
    return 4;
  }
  snprintf(cpfile, STRMAX, "%s-CP.xsf", cpoutname);
  fptr = fopen(cpfile, "w");
  OutputXSF(fptr, denin, cpout);
  fclose(fptr);
  printf("  CP file %s is finished\n\n", cpfile);
  return 0;
}

int CalcCP_prelim(struct CrystData * denin_up, struct CrystData * denin, struct CrystData * denin_dn, struct CrystData * etotin_up, struct CrystData * etotin_dn, struct CrystData * cpout) {
  /* called by: main */
  /* calls: CopyStruct, IntegrateGrid, PressureContrib, ScaleGrid, ShiftGrid, SubtractGrid, SymmetrizeGrid, OutputXSF */
  char cpfile[STRMAX];
  int check=0, voxtot=ngx*ngy*ngz;
  double dV=0.0, hart=0.0, kin=0.0, loc=0.0, ptot=0.0, xc=0.0, ew=0.0, alpha=0.0, nonloc=0.0;
  double p_check=0.0, p_entropy=0.0, p_mapcore=0, p_remain=0.0;
  FILE * fptr;
  CopyStruct(denin, cpout);  /* copies geometric information */
  dV = etotin_dn->volvox-etotin_up->volvox;  /* smaller minus larger because p = -dE/dV */
  SubtractGrid(etotin_up, etotin_dn, cpout);  /* cpout is energy difference here */
  cpout->volvox = (etotin_up->volvox+etotin_dn->volvox)/2.0;
  ScaleGrid(cpout, 1.0/dV, cpout);
  p_remain = IntegrateGrid(cpout)/(cpout->volvox*voxtot);
  if (mapcore==1 && mapsym==1) {
    if(CV_mode!=5) {
         printf("  Restoring symmetry after interpolation\n");
         check = SymmetrizeGrid(cpout, &smap);
    }
    if (check!=0) return 1;
    p_check = IntegrateGrid(cpout)/(cpout->volvox*voxtot);
    if (fabs(p_remain-p_check)>1.0e-12) {
      fprintf(cplog, "WARNING: total pressure changed after symmetrization\n");
      fprintf(cplog, " before: %20.14f\n", p_remain);
      fprintf(cplog, "  after: %20.14f\n", p_check);
      errcount++;
      if (fabs(p_remain-p_check)<1.0e-5) printf("\n  CAUTION: Total pressure has changed! Continuing anyway...\n");
      else {
        printf("\n  BAD NEWS: Total pressure has changed!\n");
        fprintf(cplog, "\nTerminated because symmetrization caused a change in pressure\n");
        fprintf(cplog, "Suggestion: check Abinit outfile for correct nsym, symrel, and tnons\n");
        return 2;
      }
    }
  }
  return 0;
}


/* MAP INTEGRATION FUNCTIONS */

int SetBubbles(struct CrystData * gridin) {
  /* called by: main */
  /* calls: ElementName */
  char element[STRMAX];
  int i=0, j=0, stop[NIONMAX];
  for (i=0; i<gridin->nion; i++) {
    if (stop[i]==1) continue;
    ElementName(gridin->zatomic[i], element);
    if (smap.nequiv[i]>1) {
      fprintf(cplog, "Enter a core radius for atom #%d in Angstrom (%s, %d equivalent sites): ",
        i+1, element, smap.nequiv[i]);
    } else {
      fprintf(cplog, "Enter a core radius for atom #%d in Angstrom (%s, %d site): ",
        i+1, element, smap.nequiv[i]);
    }
    retry_radius:
    if (smap.nequiv[i]>1) {
      printf("  Enter a core radius for atom #%d in Angstrom (%s, %d equivalent sites): ",
        i+1, element, smap.nequiv[i]);
    } else {
      printf("  Enter a core radius for atom #%d in Angstrom (%s, %d site): ",
        i+1, element, smap.nequiv[i]);
    }
    scanf("%lf", gridin->corerad[i]);
    if (gridin->corerad[i]<0.0) {
      printf("  Invalid input. Please try again\n");
      goto retry_radius;
    }
    fprintf(cplog, "%f Angstrom = %f bohr\n", gridin->corerad[i], gridin->corerad[i]/R_BOHR);
    gridin->corerad[i] = gridin->corerad[i]/R_BOHR;
    for (j=i; j<gridin->nion; j++) {
      if (smap.equiv[i][j]!=1) continue;
      gridin->corerad[j] = gridin->corerad[i];
      stop[j] = 1;
    }
  }
  return 0;
}

int CoordSearchBaderHirsh(struct CrystData * gridin, struct ContactVol * map,struct SymMap * smap) {
  /* called by: main */
  /* calls: Getwj */
  int atom=0, atom2=0, atom3[NEQVOX], index=0,index2, i=0, j=0, k,jx=0, jy=0, jz=0;
  int ka1=0, kb1=0, kc1=0, ka2=0, kb2=0, kc2=0, ka3[NEQVOX], kb3[NEQVOX], kc3[NEQVOX];
  int ka=0, kb=0, kc=0, check[NEQVOX], ngcount=0, ngp=0, ngp0=0, voxtot=ngx*ngy*ngz;
  double interatom_dist,dist=0.0,distmin, dmax=0.0, voxcenter_x=0.0, voxcenter_y=0.0, voxcenter_z=0.0;
  double wj_temp[NEQVOX], wmax=0.0, wmax2=0.0, xf=0.0, yf=0.0, zf=0.0, xc=0.0, yc=0.0, zc=0.0, xc2=0.0, yc2=0.0, zc2=0.0;
  double cella=0.0, cellb=0.0, cellc=0.0, minstep=0.0, tolerance=0.0, r_max;
  double stepx,stepy,stepz,testx,testy,testz,gradx,grady,gradz,grad_norm;
  double xf2,yf2,zf2;
  int stop,step,jx2,jy2,jz2,jx2c,jy2c,jz2c;
  double weightmatrix[NEQVOX][NEQVOX];
  double weightmatrix2[NEQVOX][NEQVOX];
  int neighborindex[NEQVOX];
  double weightmatrix_max;
  double hirshsum;
  int contactcounter,counter;
  int stepmax;
  int step_upper,step_lower,step_new;
  int ionmap_temp,ionmap_temp2;
  int ionmap_list[NIONMAX*27];
  double ionmap_dist[NIONMAX*27],diffx,diffy,diffz;
  double iondist_temp;
  char bader_name[100];
  int foundit;
  int shared_voxels,symm_op;
  int ka0,kb0,kc0,atom0;
  cella = pow(gridin->cella_x*gridin->cella_x+gridin->cella_y*gridin->cella_y+gridin->cella_z*gridin->cella_z,0.5)/(ngx);
  cellb = pow(gridin->cellb_x*gridin->cellb_x+gridin->cellb_y*gridin->cellb_y+gridin->cellb_z*gridin->cellb_z,0.5)/(ngy);
  cellc = pow(gridin->cellc_x*gridin->cellc_x+gridin->cellc_y*gridin->cellc_y+gridin->cellc_z*gridin->cellc_z,0.5)/(ngz);
  minstep = cella;
  if (cellb<minstep) minstep = cellb;
  if (cellc<minstep) minstep = cellc;
  tolerance = 0.0004*minstep;
  /* every voxel in the unit cell */
  for (jz=0; jz<ngz; jz++) {
       for (jy=0; jy<ngy; jy++) {
         for (jx=0; jx<ngx; jx++) {
            map->neighcount[jx][jy][jz] = 0;
            map->neighcount2[jx][jy][jz] = 0;
            map->swj[jx][jy][jz]=0.0;
            map->swjk[jx][jy][jz]=0.0;
         }
       }
  }
  if(CV_mode==3) {
    for (jz=0; jz<ngz; jz++) {
       for (jy=0; jy<ngy; jy++) {
         for (jx=0; jx<ngx; jx++) {
              delta_pot_map.grid[jx][jy][jz]=0.0;
         }
       }
    }
  }
  printf("Reading and symmetrizing Bader volumes.\n");
  for(atom=0; atom<gridin->nion; atom++) {
    sprintf(bader_name,"BvAt%04d.dat",atom+1);
    printf("    Processing %s for atom %d (%lf, %lf, %lf)\n",bader_name,atom+1,gridin->x[atom],gridin->y[atom],gridin->z[atom]);
    ReadCHGCAR(bader_name,&bader);
    for(symm_op=0;symm_op<smap->nsymel;symm_op++) {
       for(atom2=0;atom2<gridin->nion;atom2++) {
          xf2 = gridin->x[atom2]*smap->symrel[symm_op][0][0]+gridin->y[atom2]*smap->symrel[symm_op][0][1]+gridin->z[atom2]*smap->symrel[symm_op][0][2]+smap->tnons[symm_op][0];
          yf2 = gridin->x[atom2]*smap->symrel[symm_op][1][0]+gridin->y[atom2]*smap->symrel[symm_op][1][1]+gridin->z[atom2]*smap->symrel[symm_op][1][2]+smap->tnons[symm_op][1];
          zf2 = gridin->x[atom2]*smap->symrel[symm_op][2][0]+gridin->y[atom2]*smap->symrel[symm_op][2][1]+gridin->z[atom2]*smap->symrel[symm_op][2][2]+smap->tnons[symm_op][2];
          xc2 = xf2*gridin->cella_x+yf2*gridin->cellb_x+zf2*gridin->cellc_x;
          yc2 = xf2*gridin->cella_y+yf2*gridin->cellb_y+zf2*gridin->cellc_y;
          zc2 = xf2*gridin->cella_z+yf2*gridin->cellb_z+zf2*gridin->cellc_z;
          for (ka=-kam; ka<=kam; ka++) {
           for (kb=-kbm; kb<=kbm; kb++) {
            for (kc=-kcm; kc<=kcm; kc++) {
                     xc = xc2+ka*gridin->cella_x+kb*gridin->cellb_x+kc*gridin->cellc_x;
                     yc = yc2+ka*gridin->cella_y+kb*gridin->cellb_y+kc*gridin->cellc_y;
                     zc = zc2+ka*gridin->cella_z+kb*gridin->cellb_z+kc*gridin->cellc_z;
                     diffx = xc-gridin->xcart[atom];
                     diffy = yc-gridin->ycart[atom];
                     diffz = zc-gridin->zcart[atom];
                     dist=pow(diffx*diffx+diffy*diffy+diffz*diffz,0.5);
                     if(dist < 0.01) {
                              sprintf(bader_name,"BvAt%04d.dat",atom2+1);
                              ReadCHGCAR(bader_name,&bader_temp);
                              for (jz=0; jz<ngz; jz++) {
                                zf = (double)jz/(double)ngz;
                                for (jy=0; jy<ngy; jy++) {
                                  yf = (double)jy/(double)ngy;
                                  for (jx=0; jx<ngx; jx++) {
                                      xf = (double)jx/(double)ngx;
                                      xf2 = xf*smap->symrel[symm_op][0][0]+yf*smap->symrel[symm_op][0][1]+zf*smap->symrel[symm_op][0][2]+smap->tnons[symm_op][0];
                                      yf2 = xf*smap->symrel[symm_op][1][0]+yf*smap->symrel[symm_op][1][1]+zf*smap->symrel[symm_op][1][2]+smap->tnons[symm_op][1];
                                      zf2 = xf*smap->symrel[symm_op][2][0]+yf*smap->symrel[symm_op][2][1]+zf*smap->symrel[symm_op][2][2]+smap->tnons[symm_op][2];
                                      while(xf2<0.0) xf2=xf2+1.000;
                                      while(yf2<0.0) yf2=yf2+1.000;
                                      while(zf2<0.0) zf2=zf2+1.000;
                                      jx2=xf2*ngx+0.5;
                                      jy2=yf2*ngy+0.5;
                                      jz2=zf2*ngz+0.5;
                                      while(jx2>=ngx) jx2-=ngx;
                                      while(jy2>=ngy) jy2-=ngy;
                                      while(jz2>=ngz) jz2-=ngz;
                                      while(jx2<0) jx2+=ngx;
                                      while(jy2<0) jy2+=ngy;
                                      while(jz2<0) jz2+=ngz;
                                      bader.grid[jx2][jy2][jz2]+=bader_temp.grid[jx][jy][jz];
                                  }
                                }
                              }
                     }
            }
           }
          }
        }
    }
    shared_voxels=0;
    for (jz=0; jz<ngz; jz++) {
       /* fractional coordinates */
       zf = (double)jz/(double)ngz;
       for (jy=0; jy<ngy; jy++) {
         yf = (double)jy/(double)ngy;
         for (jx=0; jx<ngx; jx++) {
            xf = (double)jx/(double)ngx;
            if(bader.grid[jx][jy][jz]>0.000000) {
                voxcenter_x = xf*gridin->cella_x+yf*gridin->cellb_x+zf*gridin->cellc_x;
                voxcenter_y = xf*gridin->cella_y+yf*gridin->cellb_y+zf*gridin->cellc_y;
                voxcenter_z = xf*gridin->cella_z+yf*gridin->cellb_z+zf*gridin->cellc_z;
                distmin=100;
                for (ka=-kam; ka<=kam; ka++) {
                 for (kb=-kbm; kb<=kbm; kb++) {
                  for (kc=-kcm; kc<=kcm; kc++) {
                     xc = gridin->xcart[atom]+ka*gridin->cella_x+kb*gridin->cellb_x+kc*gridin->cellc_x;
                     yc = gridin->ycart[atom]+ka*gridin->cella_y+kb*gridin->cellb_y+kc*gridin->cellc_y;
                     zc = gridin->zcart[atom]+ka*gridin->cella_z+kb*gridin->cellb_z+kc*gridin->cellc_z;
                     dist = sqrt((voxcenter_x-xc)*(voxcenter_x-xc)+(voxcenter_y-yc)*(voxcenter_y-yc)+(voxcenter_z-zc)*(voxcenter_z-zc));
                     if(dist<distmin) distmin=dist;
                  }
                 }
                }
                for (ka=-kam; ka<=kam; ka++) {
                 for (kb=-kbm; kb<=kbm; kb++) {
                  for (kc=-kcm; kc<=kcm; kc++) {
                     xc = gridin->xcart[atom]+ka*gridin->cella_x+kb*gridin->cellb_x+kc*gridin->cellc_x;
                     yc = gridin->ycart[atom]+ka*gridin->cella_y+kb*gridin->cellb_y+kc*gridin->cellc_y;
                     zc = gridin->zcart[atom]+ka*gridin->cella_z+kb*gridin->cellb_z+kc*gridin->cellc_z;
                     dist = sqrt((voxcenter_x-xc)*(voxcenter_x-xc)+(voxcenter_y-yc)*(voxcenter_y-yc)+(voxcenter_z-zc)*(voxcenter_z-zc));
                     if(dist<=distmin+tolerance) {
                           index = map->neighcount2[jx][jy][jz];
                           map->ionmap[index][jx][jy][jz] = ((ka+3)<<13)+((kb+3)<<10)+((kc+3)<<7)+atom;
                           map->swj[jx][jy][jz] += 1.0;
                           map->neighcount2[jx][jy][jz]++;
                           if(index>0) shared_voxels++;
                     }
                  }
                 }
                }
            }
         }
       }
    }
    printf("    %d shared voxels for atom %d.\n",shared_voxels,atom+1);
  }
  ngcount=0;
  if(CV_mode==7) {
  printf("Determining neighbor atoms for each voxel...\n");
  printf("0%%");
  fflush(stdout);
  for (jz=0; jz<ngz; jz++) {
     /* fractional coordinates */
     zf = (double)jz/(double)ngz;
     for (jy=0; jy<ngy; jy++) {
       yf = (double)jy/(double)ngy;
       for (jx=0; jx<ngx; jx++) {
        xf = (double)jx/(double)ngx;
          ngcount++;
          ngp = ngcount*100/voxtot;
          if (ngp!=ngp0) {
            printf("\r%d%%", ngp);
            fflush(stdout);
          }
          ngp0=ngp;
          contactcounter=0;
          voxcenter_x = xf*gridin->cella_x+yf*gridin->cellb_x+zf*gridin->cellc_x;
          voxcenter_y = xf*gridin->cella_y+yf*gridin->cellb_y+zf*gridin->cellc_y;
          voxcenter_z = xf*gridin->cella_z+yf*gridin->cellb_z+zf*gridin->cellc_z;
          distmin=1000.00;
          hirshsum=0.0;
          if(map->neighcount2[jx][jy][jz] < 2) {
                atom0 = map->ionmap[0][jx][jy][jz]&127;
                ka0 = (map->ionmap[0][jx][jy][jz]>>13&7)-3;
                kb0 = (map->ionmap[0][jx][jy][jz]>>10&7)-3;
                kc0 = (map->ionmap[0][jx][jy][jz]>>7&7)-3;
                xc2 = gridin->xcart[atom0]+ka0*gridin->cella_x+kb0*gridin->cellb_x+kc0*gridin->cellc_x;
                yc2 = gridin->ycart[atom0]+ka0*gridin->cella_y+kb0*gridin->cellb_y+kc0*gridin->cellc_y;
                zc2 = gridin->zcart[atom0]+ka0*gridin->cella_z+kb0*gridin->cellb_z+kc0*gridin->cellc_z;
                map->neighcount[jx][jy][jz] = 1;
                distmin=1000.00;
                counter=0;
                for (ka=-kam; ka<=kam; ka++) {
                 for (kb=-kbm; kb<=kbm; kb++) {
                  for (kc=-kcm; kc<=kcm; kc++) {
                   for(atom=0; atom<gridin->nion; atom++) {
                      ionmap_temp = ((ka+3)<<13)+((kb+3)<<10)+((kc+3)<<7)+atom;
                      if(ionmap_temp != map->ionmap[0][jx][jy][jz]) {
                         /* found neighbor to test */
                       xc = gridin->xcart[atom]+ka*gridin->cella_x+kb*gridin->cellb_x+kc*gridin->cellc_x;
                       yc = gridin->ycart[atom]+ka*gridin->cella_y+kb*gridin->cellb_y+kc*gridin->cellc_y;
                       zc = gridin->zcart[atom]+ka*gridin->cella_z+kb*gridin->cellb_z+kc*gridin->cellc_z;
                       interatom_dist = sqrt((xc2-xc)*(xc2-xc)+(yc2-yc)*(yc2-yc)+(zc2-zc)*(zc2-zc));
                       if(interatom_dist < 5.2/0.529) {
                         dist = sqrt((voxcenter_x-xc)*(voxcenter_x-xc)+(voxcenter_y-yc)*(voxcenter_y-yc)+(voxcenter_z-zc)*(voxcenter_z-zc));
                         iondist_temp = Getwj(atom, dist);
                         if(iondist_temp > 0.0001) {
                           ionmap_dist[counter] = iondist_temp;
                           hirshsum+=ionmap_dist[counter];
                           ionmap_list[counter]=ionmap_temp;
                           counter++;
                           if(counter>NEQVOX-1) printf("Warning voxel has more than NEQVOX neighbors!\n");
                         }
                       }
                      }
                   }
                  }
                 }
                }
                for(k=0;k<counter;k++) {
                           index=map->neighcount[jx][jy][jz];
                           map->ionmap[index][jx][jy][jz]=ionmap_list[k];
                           map->neighcount[jx][jy][jz]++;
                }
                contactcounter=0;
                for (i=0; i<map->neighcount[jx][jy][jz]; i++) {
                  for (j=i+1; j<map->neighcount[jx][jy][jz]; j++) {
                    map->neighkey[contactcounter][jx][jy][jz]=0;
                    if (i==0) {
                     map->neighkey[contactcounter][jx][jy][jz]=1;
                     map->swjk[jx][jy][jz] += ionmap_dist[j-1]/hirshsum;
                     map->wj[j][jx][jy][jz]=ionmap_dist[j-1]/hirshsum;
                    }
                    contactcounter++;
                    if(contactcounter>NEQVOX-1) {
                         printf("Warning voxel is involved in more than NEQVOX contacts!\n");
                         exit(0);
                    }
                  }
                }
          }
          if(map->neighcount2[jx][jy][jz] > 1) {
              map->neighcount[jx][jy][jz] = map->neighcount2[jx][jy][jz];
              for (i=0; i<map->neighcount[jx][jy][jz]; i++) {
                  for (j=i+1; j<map->neighcount[jx][jy][jz]; j++) {
                    map->neighkey[contactcounter][jx][jy][jz]=1;
                    map->swjk[jx][jy][jz] += 1.0;
                    contactcounter++;
                    if(contactcounter>NEQVOX-1) printf("Warning voxel is involved in more than NEQVOX contacts!\n");
                  }
               }
          }
       }
     }
  }
  printf(" ...Finished\n");
  }
  return 0;
}


int CoordSearchBader(struct CrystData * gridin, struct ContactVol * map,struct SymMap * smap) {
  /* called by: main */
  /* calls: Getwj */
  int atom=0, atom2=0, atom3[NEQVOX], index=0,index2, i=0, j=0, k,jx=0, jy=0, jz=0;
  int ka1=0, kb1=0, kc1=0, ka2=0, kb2=0, kc2=0, ka3[NEQVOX], kb3[NEQVOX], kc3[NEQVOX];
  int ka=0, kb=0, kc=0, check[NEQVOX], ngcount=0, ngp=0, ngp0=0, voxtot=ngx*ngy*ngz;
  double interatom_dist,dist=0.0,distmin, dmax=0.0, voxcenter_x=0.0, voxcenter_y=0.0, voxcenter_z=0.0;
  double wj_temp[NEQVOX], wmax=0.0, wmax2=0.0, xf=0.0, yf=0.0, zf=0.0, xc=0.0, yc=0.0, zc=0.0, xc2=0.0, yc2=0.0, zc2=0.0;
  double cella=0.0, cellb=0.0, cellc=0.0, minstep=0.0, tolerance=0.0, r_max;
  double stepx,stepy,stepz,testx,testy,testz,gradx,grady,gradz,grad_norm;
  double xf2,yf2,zf2;
  int stop,step,jx2,jy2,jz2,jx2c,jy2c,jz2c;
  double weightmatrix[NEQVOX][NEQVOX];
  double weightmatrix2[NEQVOX][NEQVOX];
  int neighborindex[NEQVOX];
  double weightmatrix_max;
  int contactcounter,counter;
  int stepmax;
  int step_upper,step_lower,step_new;
  int ionmap_temp,ionmap_temp2;
  int ionmap_list[NIONMAX*27];
  double ionmap_dist[NIONMAX*27],diffx,diffy,diffz;
  char bader_name[100];
  int foundit;
  int shared_voxels,symm_op;
  int ka0,kb0,kc0,atom0;
  cella = pow(gridin->cella_x*gridin->cella_x+gridin->cella_y*gridin->cella_y+gridin->cella_z*gridin->cella_z,0.5)/(ngx);
  cellb = pow(gridin->cellb_x*gridin->cellb_x+gridin->cellb_y*gridin->cellb_y+gridin->cellb_z*gridin->cellb_z,0.5)/(ngy);
  cellc = pow(gridin->cellc_x*gridin->cellc_x+gridin->cellc_y*gridin->cellc_y+gridin->cellc_z*gridin->cellc_z,0.5)/(ngz);
  minstep = cella;
  if (cellb<minstep) minstep = cellb;
  if (cellc<minstep) minstep = cellc;
  tolerance = 0.0004*minstep;
  /* every voxel in the unit cell */
  for (jz=0; jz<ngz; jz++) {
       for (jy=0; jy<ngy; jy++) {
         for (jx=0; jx<ngx; jx++) {
            map->neighcount[jx][jy][jz] = 0;
            map->neighcount2[jx][jy][jz] = 0;
            map->swj[jx][jy][jz]=0.0;
            map->swjk[jx][jy][jz]=0.0;
         }
       }
  }
  if(CV_mode==3) {
    for (jz=0; jz<ngz; jz++) {
       for (jy=0; jy<ngy; jy++) {
         for (jx=0; jx<ngx; jx++) {
              delta_pot_map.grid[jx][jy][jz]=0.0;
         }
       }
    }
  }
  printf("Reading and symmetrizing Bader volumes.\n");
  for(atom=0; atom<gridin->nion; atom++) {
    sprintf(bader_name,"BvAt%04d.dat",atom+1);
    printf("    Processing %s for atom %d (%lf, %lf, %lf)\n",bader_name,atom+1,gridin->x[atom],gridin->y[atom],gridin->z[atom]);
    ReadCHGCAR(bader_name,&bader);
    for(symm_op=0;symm_op<smap->nsymel;symm_op++) {
       for(atom2=0;atom2<gridin->nion;atom2++) {
          xf2 = gridin->x[atom2]*smap->symrel[symm_op][0][0]+gridin->y[atom2]*smap->symrel[symm_op][0][1]+gridin->z[atom2]*smap->symrel[symm_op][0][2]+smap->tnons[symm_op][0];
          yf2 = gridin->x[atom2]*smap->symrel[symm_op][1][0]+gridin->y[atom2]*smap->symrel[symm_op][1][1]+gridin->z[atom2]*smap->symrel[symm_op][1][2]+smap->tnons[symm_op][1];
          zf2 = gridin->x[atom2]*smap->symrel[symm_op][2][0]+gridin->y[atom2]*smap->symrel[symm_op][2][1]+gridin->z[atom2]*smap->symrel[symm_op][2][2]+smap->tnons[symm_op][2];
          xc2 = xf2*gridin->cella_x+yf2*gridin->cellb_x+zf2*gridin->cellc_x;
          yc2 = xf2*gridin->cella_y+yf2*gridin->cellb_y+zf2*gridin->cellc_y;
          zc2 = xf2*gridin->cella_z+yf2*gridin->cellb_z+zf2*gridin->cellc_z;
          for (ka=-kam; ka<=kam; ka++) {
           for (kb=-kbm; kb<=kbm; kb++) {
            for (kc=-kcm; kc<=kcm; kc++) {
                     xc = xc2+ka*gridin->cella_x+kb*gridin->cellb_x+kc*gridin->cellc_x;
                     yc = yc2+ka*gridin->cella_y+kb*gridin->cellb_y+kc*gridin->cellc_y;
                     zc = zc2+ka*gridin->cella_z+kb*gridin->cellb_z+kc*gridin->cellc_z;
                     diffx = xc-gridin->xcart[atom];
                     diffy = yc-gridin->ycart[atom];
                     diffz = zc-gridin->zcart[atom];
                     dist=pow(diffx*diffx+diffy*diffy+diffz*diffz,0.5);
                     if(dist < 0.01) {
                              sprintf(bader_name,"BvAt%04d.dat",atom2+1);
                              ReadCHGCAR(bader_name,&bader_temp);
                              for (jz=0; jz<ngz; jz++) {
                                zf = (double)jz/(double)ngz;
                                for (jy=0; jy<ngy; jy++) {
                                  yf = (double)jy/(double)ngy;
                                  for (jx=0; jx<ngx; jx++) {
                                      xf = (double)jx/(double)ngx;
                                      xf2 = xf*smap->symrel[symm_op][0][0]+yf*smap->symrel[symm_op][0][1]+zf*smap->symrel[symm_op][0][2]+smap->tnons[symm_op][0];
                                      yf2 = xf*smap->symrel[symm_op][1][0]+yf*smap->symrel[symm_op][1][1]+zf*smap->symrel[symm_op][1][2]+smap->tnons[symm_op][1];
                                      zf2 = xf*smap->symrel[symm_op][2][0]+yf*smap->symrel[symm_op][2][1]+zf*smap->symrel[symm_op][2][2]+smap->tnons[symm_op][2];
                                      while(xf2<0.0) xf2=xf2+1.000;
                                      while(yf2<0.0) yf2=yf2+1.000;
                                      while(zf2<0.0) zf2=zf2+1.000;
                                      jx2=xf2*ngx+0.5;
                                      jy2=yf2*ngy+0.5;
                                      jz2=zf2*ngz+0.5;
                                      while(jx2>=ngx) jx2-=ngx;
                                      while(jy2>=ngy) jy2-=ngy;
                                      while(jz2>=ngz) jz2-=ngz;
                                      while(jx2<0) jx2+=ngx;
                                      while(jy2<0) jy2+=ngy;
                                      while(jz2<0) jz2+=ngz;
                                      bader.grid[jx2][jy2][jz2]+=bader_temp.grid[jx][jy][jz];
                                  }
                                }
                              }
                     }
            }
           }
          }
        }
    }
    shared_voxels=0;
    for (jz=0; jz<ngz; jz++) {
       /* fractional coordinates */
       zf = (double)jz/(double)ngz;
       for (jy=0; jy<ngy; jy++) {
         yf = (double)jy/(double)ngy;
         for (jx=0; jx<ngx; jx++) {
            xf = (double)jx/(double)ngx;
            if(bader.grid[jx][jy][jz]>0.000000) {
                voxcenter_x = xf*gridin->cella_x+yf*gridin->cellb_x+zf*gridin->cellc_x;
                voxcenter_y = xf*gridin->cella_y+yf*gridin->cellb_y+zf*gridin->cellc_y;
                voxcenter_z = xf*gridin->cella_z+yf*gridin->cellb_z+zf*gridin->cellc_z;
                distmin=100;
                for (ka=-kam; ka<=kam; ka++) {
                 for (kb=-kbm; kb<=kbm; kb++) {
                  for (kc=-kcm; kc<=kcm; kc++) {
                     xc = gridin->xcart[atom]+ka*gridin->cella_x+kb*gridin->cellb_x+kc*gridin->cellc_x;
                     yc = gridin->ycart[atom]+ka*gridin->cella_y+kb*gridin->cellb_y+kc*gridin->cellc_y;
                     zc = gridin->zcart[atom]+ka*gridin->cella_z+kb*gridin->cellb_z+kc*gridin->cellc_z;
                     dist = sqrt((voxcenter_x-xc)*(voxcenter_x-xc)+(voxcenter_y-yc)*(voxcenter_y-yc)+(voxcenter_z-zc)*(voxcenter_z-zc));
                     if(dist<distmin) distmin=dist;
                  }
                 }
                }
                for (ka=-kam; ka<=kam; ka++) {
                 for (kb=-kbm; kb<=kbm; kb++) {
                  for (kc=-kcm; kc<=kcm; kc++) {
                     xc = gridin->xcart[atom]+ka*gridin->cella_x+kb*gridin->cellb_x+kc*gridin->cellc_x;
                     yc = gridin->ycart[atom]+ka*gridin->cella_y+kb*gridin->cellb_y+kc*gridin->cellc_y;
                     zc = gridin->zcart[atom]+ka*gridin->cella_z+kb*gridin->cellb_z+kc*gridin->cellc_z;
                     dist = sqrt((voxcenter_x-xc)*(voxcenter_x-xc)+(voxcenter_y-yc)*(voxcenter_y-yc)+(voxcenter_z-zc)*(voxcenter_z-zc));
                     if(dist<=distmin+tolerance) {
                           index = map->neighcount2[jx][jy][jz];
                           map->ionmap[index][jx][jy][jz] = ((ka+3)<<13)+((kb+3)<<10)+((kc+3)<<7)+atom;
                           map->swj[jx][jy][jz] += 1.0;
                           map->neighcount2[jx][jy][jz]++;
                           if(index>0) shared_voxels++;
                     }
                  }
                 }
                }
            }
         }
       }
    }
    printf("    %d shared voxels for atom %d.\n",shared_voxels,atom+1);
  }
  ngcount=0;
  if((CV_mode<5)||(CV_mode==7)) {
  printf("Determining neighbor atoms for each voxel...\n");
  printf("0%%");
  fflush(stdout);
  for (jz=0; jz<ngz; jz++) {
     /* fractional coordinates */
     zf = (double)jz/(double)ngz;
     for (jy=0; jy<ngy; jy++) {
       yf = (double)jy/(double)ngy;
       for (jx=0; jx<ngx; jx++) {
        xf = (double)jx/(double)ngx;
          ngcount++;
          ngp = ngcount*100/voxtot;
          if (ngp!=ngp0) {
            printf("\r%d%%", ngp);
            fflush(stdout);
          }
          ngp0=ngp;
          contactcounter=0;
          voxcenter_x = xf*gridin->cella_x+yf*gridin->cellb_x+zf*gridin->cellc_x;
          voxcenter_y = xf*gridin->cella_y+yf*gridin->cellb_y+zf*gridin->cellc_y;
          voxcenter_z = xf*gridin->cella_z+yf*gridin->cellb_z+zf*gridin->cellc_z;
          distmin=1000.00;
          if(map->neighcount2[jx][jy][jz] < 2) {
                atom0 = map->ionmap[0][jx][jy][jz]&127;
                ka0 = (map->ionmap[0][jx][jy][jz]>>13&7)-3;
                kb0 = (map->ionmap[0][jx][jy][jz]>>10&7)-3;
                kc0 = (map->ionmap[0][jx][jy][jz]>>7&7)-3;
                xc2 = gridin->xcart[atom0]+ka0*gridin->cella_x+kb0*gridin->cellb_x+kc0*gridin->cellc_x;
                yc2 = gridin->ycart[atom0]+ka0*gridin->cella_y+kb0*gridin->cellb_y+kc0*gridin->cellc_y;
                zc2 = gridin->zcart[atom0]+ka0*gridin->cella_z+kb0*gridin->cellb_z+kc0*gridin->cellc_z;
                map->neighcount[jx][jy][jz] = 1;
                distmin=1000.00;
                counter=0;
                if(CV_mode==4) grad_norm = GridGrad(jx,jy,jz,&delta_pot_map, &gradx,&grady,&gradz);
                for (ka=-kam; ka<=kam; ka++) {
                 for (kb=-kbm; kb<=kbm; kb++) {
                  for (kc=-kcm; kc<=kcm; kc++) {
                   for(atom=0; atom<gridin->nion; atom++) {
                      ionmap_temp = ((ka+3)<<13)+((kb+3)<<10)+((kc+3)<<7)+atom;
                      if(ionmap_temp != map->ionmap[0][jx][jy][jz]) {
                         /* found neighbor to test */
                       xc = gridin->xcart[atom]+ka*gridin->cella_x+kb*gridin->cellb_x+kc*gridin->cellc_x;
                       yc = gridin->ycart[atom]+ka*gridin->cella_y+kb*gridin->cellb_y+kc*gridin->cellc_y;
                       zc = gridin->zcart[atom]+ka*gridin->cella_z+kb*gridin->cellb_z+kc*gridin->cellc_z;
                       interatom_dist = sqrt((xc2-xc)*(xc2-xc)+(yc2-yc)*(yc2-yc)+(zc2-zc)*(zc2-zc));
                       if(interatom_dist < 5.2/0.529) {
                         dist = sqrt((voxcenter_x-xc)*(voxcenter_x-xc)+(voxcenter_y-yc)*(voxcenter_y-yc)+(voxcenter_z-zc)*(voxcenter_z-zc));
                         if(CV_mode==1) {
                           stepmax = dist/tolerance+0.5;
                           foundit=0;
                           step_lower=0;
                           step_upper=stepmax;
                           step = 0;
                           while(foundit==0) {
                               if(foundit==1) printf("loop escaped!\n");
                               step = round(1.0*(step_lower+step_upper)/2.0);
                               xf2 = xf + step*(gridin->x[atom]+1.0*ka-xf)/(double)stepmax;
                               yf2 = yf + step*(gridin->y[atom]+1.0*kb-yf)/(double)stepmax;
                               zf2 = zf + step*(gridin->z[atom]+1.0*kc-zf)/(double)stepmax;
                               jx2 = round(xf2*ngx);
                               jy2 = round(yf2*ngy);
                               jz2 = round(zf2*ngz);
                               jx2c = (jx2 + 4*ngx)%ngx;
                               jy2c = (jy2 + 4*ngy)%ngy;
                               jz2c = (jz2 + 4*ngz)%ngz;
                               ka2 = (jx2-jx2c)/ngx;
                               kb2 = (jy2-jy2c)/ngy;
                               kc2 = (jz2-jz2c)/ngz;
                               ionmap_temp2 = ((ka-ka2+3)<<13)+((kb-kb2+3)<<10)+((kc-kc2+3)<<7)+atom;
                               if((map->neighcount2[jx2c][jy2c][jz2c]==1)&&(map->ionmap[0][jx2c][jy2c][jz2c]==ionmap_temp2)) {
                                      step_upper = step;
                               }
                               else {
                                      step_lower = step;
                               }
                               if(step_upper-step_lower<2) foundit=1;
                               if(foundit==1) {
                                           ionmap_dist[counter]=step_upper*tolerance;
                                           ionmap_list[counter]=ionmap_temp;
                                           if(ionmap_dist[counter]<distmin) distmin=ionmap_dist[counter];
                                           counter++;
                               }
                           }
                         }
                         if(CV_mode==2) {
                                           ionmap_dist[counter]=dist;
                                           ionmap_list[counter]=ionmap_temp;
                                           if(ionmap_dist[counter]<distmin) distmin=ionmap_dist[counter];
                                           counter++;
                         }
                         if(CV_mode==3) {
                                         ionmap_dist[counter] = -Getwj(atom, dist);
                                         delta_pot_map.grid[jx][jy][jz]+=ionmap_dist[counter];
                                         ionmap_list[counter]=ionmap_temp;
                                         if(ionmap_dist[counter]<distmin) distmin=ionmap_dist[counter];
                                         counter++;
                         }
                         if((CV_mode==4)) {
                                         if(fabs(grad_norm) > 0.0) {
                                            ionmap_dist[counter] = -(gradx*(xc2-xc) + grady*(yc2-yc) + gradz*(zc2-zc))/grad_norm/interatom_dist;
                                         }
                                         else ionmap_dist[counter]=dist;
                                         ionmap_list[counter]=ionmap_temp;
                                         if(ionmap_dist[counter]<distmin) distmin=ionmap_dist[counter];
                                         counter++;
                         }
                       }
                      }
                   }
                  }
                 }
                }
                for(k=0;k<counter;k++) {
                       if(ionmap_dist[k]<distmin+tolerance) {
                           index=map->neighcount[jx][jy][jz];
                           map->ionmap[index][jx][jy][jz]=ionmap_list[k];
                           map->neighcount[jx][jy][jz]++;
                       }
                }
                contactcounter=0;
                for (i=0; i<map->neighcount[jx][jy][jz]; i++) {
                  for (j=i+1; j<map->neighcount[jx][jy][jz]; j++) {
                    map->neighkey[contactcounter][jx][jy][jz]=0;
                    if (i==0) {
                     map->neighkey[contactcounter][jx][jy][jz]=1;
                     map->swjk[jx][jy][jz] += 1.0;
                    }
                    contactcounter++;
                  }
                }
          }
          if(map->neighcount2[jx][jy][jz] > 1) {
              map->neighcount[jx][jy][jz] = map->neighcount2[jx][jy][jz];
              for (i=0; i<map->neighcount[jx][jy][jz]; i++) {
                  for (j=i+1; j<map->neighcount[jx][jy][jz]; j++) {
                    map->neighkey[contactcounter][jx][jy][jz]=1;
                    map->swjk[jx][jy][jz] += 1.0;
                    contactcounter++;
                  }
               }
          }
       }
     }
  }
  printf(" ...Finished\n");
  }
  return 0;
}


int CoreIntegrate(struct CrystData * gridin, struct CrystData * gridin2,struct ContactVol * map) {
  /* called by: main */
  /* calls: Getwj */
  int atom=0, atom2=0, atom3[NEQVOX], index=0,index2, i=0, j=0, k,jx=0, jy=0, jz=0;
  int ka1=0, kb1=0, kc1=0, ka2=0, kb2=0, kc2=0, ka3[NEQVOX], kb3[NEQVOX], kc3[NEQVOX];
  int ka=0, kb=0, kc=0, check[NEQVOX], ngcount=0, ngp=0, ngp0=0, voxtot=ngx*ngy*ngz;
  double dist=0.0,dist2,distmin, dmax=0.0, voxcenter_x=0.0, voxcenter_y=0.0, voxcenter_z=0.0;
  double voxcenter_x2=0.0, voxcenter_y2=0.0, voxcenter_z2=0.0;
  double wj_temp[NEQVOX], wmax=0.0, wmax2=0.0, xf=0.0, yf=0.0, zf=0.0, xc=0.0, yc=0.0, zc=0.0, xc2=0.0, yc2=0.0, zc2=0.0;
  double cella=0.0, cellb=0.0, cellc=0.0, minstep=0.0, tolerance=0.0, r_max;
  double stepx,stepy,stepz,testx,testy,testz;
  double xf2,yf2,zf2;
  double diffx,diffy,diffz;
  int stop,step,jx2,jy2,jz2,jx2c,jy2c,jz2c;
  double weightmatrix[NEQVOX][NEQVOX];
  double weightmatrix2[NEQVOX][NEQVOX];
  int neighborindex[NEQVOX];
  double weightmatrix_max;
  double sqrt2;
  double den_i,den_i_dV;
  double nelectrons_total=0.0;
  int contactcounter,counter;
  int stepmax;
  int ionmap_temp,ionmap_temp2,type;
  int ionmap_list[NIONMAX*27];
  double rvox,cc1,cc2,cc3,cc4,ionmap_dist[NIONMAX*27];
  double rloc,r_rloc,r_rloc2,r_rloc4,r_rloc6,volvox;
  double charge_fraction[NIONMAX],charge_total;
  double Madelung_fraction[NIONMAX], MappedE[NIONMAX];
  double DFTMadelungE;
  double Hartree_check;
  double lPot_check;
  double Ewald_remainder;
  double CoreCharges[NIONMAX][10];
  double cutoffs[10];
  double positivecutoff;
  double den_up,den_down;
  double temp_zion[NIONMAX];
  double cutoff_profile[NIONMAX][NPOINTMAX];;
  int voxcount,voxelstep,voxelsteps3D;
  FILE * f2;
  char profile_name[100],element[5];
  cella = pow(gridin->cella_x*gridin->cella_x+gridin->cella_y*gridin->cella_y+gridin->cella_z*gridin->cella_z,0.5)/(ngx);
  cellb = pow(gridin->cellb_x*gridin->cellb_x+gridin->cellb_y*gridin->cellb_y+gridin->cellb_z*gridin->cellb_z,0.5)/(ngy);
  cellc = pow(gridin->cellc_x*gridin->cellc_x+gridin->cellc_y*gridin->cellc_y+gridin->cellc_z*gridin->cellc_z,0.5)/(ngz);
  minstep = cella;
  sqrt2=pow(2.0,0.5);
  cutoffs[0] = 0.016;
  cutoffs[1] = 0.018;
  cutoffs[2] = 0.020;
  cutoffs[3] = 0.022;
  cutoffs[4] = 0.024;
  cutoffs[5] = 0.026;
  cutoffs[6] = 0.028;
  cutoffs[7] = 0.030;
  /* Volume of sphere = 4/3 pi r^3  */
  rvox = pow(0.75*gridin->volvox/PI,0.3333333333333333);
  volvox=gridin->volvox;
  printf("\nCore charge tables\n");
  for (atom=0; atom<gridin->nion; atom++) {
    for(j=0;j<10;j++) CoreCharges[atom][j]=0.0;
    for(j=0;j<200;j++) cutoff_profile[atom][j]=0.0;
  }
  for (jz=0; jz<ngz; jz++) {
      zf = (double)jz/(double)ngz;
      for (jy=0; jy<ngy; jy++) {
       yf = (double)jy/(double)ngy;
       for (jx=0; jx<ngx; jx++) {
        xf = (double)jx/(double)ngx;
         gridin2->grid[jx][jy][jz]=0.0;
         voxcenter_x = xf*gridin->cella_x+yf*gridin->cellb_x+zf*gridin->cellc_x;
         voxcenter_y = xf*gridin->cella_y+yf*gridin->cellb_y+zf*gridin->cellc_y;
         voxcenter_z = xf*gridin->cella_z+yf*gridin->cellb_z+zf*gridin->cellc_z;
         if(map->neighcount2[jx][jy][jz]==1) {
            atom = map->ionmap[0][jx][jy][jz]&127;
            ka = (map->ionmap[0][jx][jy][jz]>>13&7)-3;
            kb = (map->ionmap[0][jx][jy][jz]>>10&7)-3;
            kc = (map->ionmap[0][jx][jy][jz]>>7&7)-3;
            xc = gridin->xcart[atom]+ka*gridin->cella_x+kb*gridin->cellb_x+kc*gridin->cellc_x;
            yc = gridin->ycart[atom]+ka*gridin->cella_y+kb*gridin->cellb_y+kc*gridin->cellc_y;
            zc = gridin->zcart[atom]+ka*gridin->cella_z+kb*gridin->cellb_z+kc*gridin->cellc_z;
            dist = sqrt((voxcenter_x-xc)*(voxcenter_x-xc)+(voxcenter_y-yc)*(voxcenter_y-yc)+(voxcenter_z-zc)*(voxcenter_z-zc));
                den_up=TricubicInterpolation(&denhi,&denhi,voxcenter_x-xc+denhi.xcart[atom],voxcenter_y-yc+denhi.ycart[atom],voxcenter_z-zc+denhi.zcart[atom]);
                den_down=TricubicInterpolation(&denlo,&denlo,voxcenter_x-xc+denlo.xcart[atom],voxcenter_y-yc+denlo.ycart[atom],voxcenter_z-zc+denlo.zcart[atom]);
                gridin2->grid[jx][jy][jz]=(den_up-den_down)*gridin->volvox/(gridin->grid[jx][jy][jz])/(denhi.volvox-denlo.volvox);
            for(j=0;j<8;j++) {
               if(fabs(gridin2->grid[jx][jy][jz])<=cutoffs[j]) CoreCharges[atom][j]+=gridin->grid[jx][jy][jz]*gridin->volvox;
            }
            for(j=0;j<200;j++) {
              if(dist<=0.01+j*0.01) cutoff_profile[atom][j]+=gridin->grid[jx][jy][jz]*gridin->volvox;
            }
            if(sc_elec[atom]<0.0) {
                  if(fabs(gridin2->grid[jx][jy][jz])<=-sc_elec[atom]) {
                         CoreCharges[atom][9]+=gridin->grid[jx][jy][jz]*gridin->volvox/map->swj[jx][jy][jz];
                  }
            }
         }
         else gridin2->grid[jx][jy][jz]=0;
       }
      }
   }
   for(j=0;j<8;j++) {
      printf("   CUTOFF = %lf TABLE\n",cutoffs[j]);
      for(atom=0;atom<gridin->nion;atom++) {
         printf("      Atom%d:  %lf\n",atom+1,CoreCharges[atom][j]);
      }
      printf("\n");
   }
   printf("   CUSTOM CUTOFF TABLE\n");
   for(atom=0;atom<gridin->nion;atom++) {
         ElementName(gridin->zatomic[atom], element);
         sprintf(profile_name,"cutoff_profile-%s%d",element,atom+1);
         f2=fopen(profile_name,"w");
         for(j=0;j<200;j++) {
               fprintf(f2,"%lf  %lf\n",0.01+j*0.01,cutoff_profile[atom][j]);
         }
         fclose(f2);
         temp_zion[atom]=vo_elec[atom] + sc_elec[atom];
        if(sc_elec[atom]<0.0)  sc_elec[atom] = CoreCharges[atom][9];
         vo_elec[atom] = temp_zion[atom] - sc_elec[atom];
        printf("      Atom%d localized electrons:  %lf\n",atom+1,sc_elec[atom]);
   }
   printf("\n");
   return 0;
}

int CoordSearchAtom(struct CrystData * gridin, struct ContactVol * map) {
  /* called by: main */
  /* calls: Getwj */
  int atom=0, atom2=0, atom3[NEQVOX], index=0,index2, i=0, j=0, jx=0, jy=0, jz=0;
  int ka1=0, kb1=0, kc1=0, ka2=0, kb2=0, kc2=0, ka3[NEQVOX], kb3[NEQVOX], kc3[NEQVOX];
  int ka=0, kb=0, kc=0, check[NEQVOX], ngcount=0, ngp=0, ngp0=0, voxtot=ngx*ngy*ngz;
  double dist=0.0, dmax=0.0, voxcenter_x=0.0, voxcenter_y=0.0, voxcenter_z=0.0;
  double wj_temp[NEQVOX], wmax=0.0, wmax2=0.0, xf=0.0, yf=0.0, zf=0.0, xc=0.0, yc=0.0, zc=0.0;
  double weightmatrix[NEQVOX][NEQVOX];
  int neighborindex[NEQVOX];
  double weightmatrix_max,rho_max;
  int contactcounter;
  printf("0%%");
  fflush(stdout);
  /* every voxel in the unit cell */
  for (jz=0; jz<ngz; jz++) {
    /* fractional coordinates */
    zf = (double)jz/(double)ngz;
    for (jy=0; jy<ngy; jy++) {
      yf = (double)jy/(double)ngy;
      for (jx=0; jx<ngx; jx++) {
        xf = (double)jx/(double)ngx;
        /* cartesian coordinates */
        voxcenter_x = xf*gridin->cella_x+yf*gridin->cellb_x+zf*gridin->cellc_x;
        voxcenter_y = xf*gridin->cella_y+yf*gridin->cellb_y+zf*gridin->cellc_y;
        voxcenter_z = xf*gridin->cella_z+yf*gridin->cellb_z+zf*gridin->cellc_z;
        for (i=0; i<NIONMAX; i++) check[i] = 0;
        map->neighcount[jx][jy][jz] = 0;
        map->neighcount2[jx][jy][jz] = 0;
        map->swj[jx][jy][jz] = 0.0;
        map->swjk[jx][jy][jz] = 0.0;
        wmax=0.0;
        wmax2=0.0;
        ngcount++;
        ngp = ngcount*100/voxtot;
        if (ngp!=ngp0) {
          printf("\r%d%%", ngp);
          printf("\r%d%%", ngp);
          fflush(stdout);
        }
        ngp0=ngp;
        /* every unit cell in the supercell */
        rho_max=0.0;
        for (ka=-kam; ka<=kam; ka++) {
          for (kb=-kbm; kb<=kbm; kb++) {
            for (kc=-kcm; kc<=kcm; kc++) {
              /* every atom in cartesian coordinates */
              for (atom=0; atom<gridin->nion; atom++) {
                xc = gridin->xcart[atom]+ka*gridin->cella_x+kb*gridin->cellb_x+kc*gridin->cellc_x;
                yc = gridin->ycart[atom]+ka*gridin->cella_y+kb*gridin->cellb_y+kc*gridin->cellc_y;
                zc = gridin->zcart[atom]+ka*gridin->cella_z+kb*gridin->cellb_z+kc*gridin->cellc_z;
                dist = sqrt((voxcenter_x-xc)*(voxcenter_x-xc)+(voxcenter_y-yc)*(voxcenter_y-yc)+
                  (voxcenter_z-zc)*(voxcenter_z-zc));
                if (dist<R_MAX) {  /* if the voxel is close enough to the translated atom */

                  index = map->neighcount2[jx][jy][jz];
                  map->ionmap[index][jx][jy][jz] = ((ka+3)<<13)+((kb+3)<<10)+((kc+3)<<7)+atom;
                  if (scheme==1) {  /* Hirshfeld weight */
                    map->wj[index][jx][jy][jz] = Getwj(atom, dist);
                    if(map->wj[index][jx][jy][jz] > rho_max) rho_max = map->wj[index][jx][jy][jz];
                    if (map->wj[index][jx][jy][jz]==-1000.0) return 1;
                  } else if (scheme==2) map->wj[index][jx][jy][jz] = R_MAX-dist+gridin->corerad[atom];  /* distance weight */
                  map->neighcount2[jx][jy][jz]++;
                  if (map->neighcount2[jx][jy][jz]==NEQVOX) {
                    printf("\n  BAD NEWS: The number of nearby atoms exceeds %d!\n", NEQVOX);
                    fprintf(cplog, "\nTerminated because the number of atoms near voxel %d %d %d is larger than %d\n",
                      jx, jy, jz, NEQVOX);
                    fprintf(cplog, "Suggestion: increase NEQVOX or decrease R_MAX and recompile\n");
                    errcount++;
                    return 2;
                  }
                }
              }
            }
          }
        }
        for(i=0;i<map->neighcount2[jx][jy][jz];i++) {
                 if(map->wj[i][jx][jy][jz] > rho_max*(1.0-tolerance)) {
                     index=map->neighcount[jx][jy][jz];
                     map->ionmap[index][jx][jy][jz]=map->ionmap[i][jx][jy][jz];
                     map->swj[jx][jy][jz]+=1.0;
                     map->neighcount[jx][jy][jz]++;
                 }
        }
        map->neighcount2[jx][jy][jz]=map->neighcount[jx][jy][jz];
       }
      }
     }
}

int CoordSearchDist(struct CrystData * gridin, struct ContactVol * map) {
  /* called by: main */
  /* calls: Getwj */
  int atom=0, atom2=0, atom3[NEQVOX], index=0,index2, i=0, j=0, jx=0, jy=0, jz=0;
  int ka1=0, kb1=0, kc1=0, ka2=0, kb2=0, kc2=0, ka3[NEQVOX], kb3[NEQVOX], kc3[NEQVOX];
  int ka=0, kb=0, kc=0, check[NEQVOX], ngcount=0, ngp=0, ngp0=0, voxtot=ngx*ngy*ngz;
  double dist=0.0, dmax=0.0, voxcenter_x=0.0, voxcenter_y=0.0, voxcenter_z=0.0;
  double wj_temp[NEQVOX], wmax=0.0, wmax2=0.0, xf=0.0, yf=0.0, zf=0.0, xc=0.0, yc=0.0, zc=0.0;
  double weightmatrix[NEQVOX][NEQVOX];
  int neighborindex[NEQVOX];
  double weightmatrix_max;
  int contactcounter;
  printf("0%%");
  fflush(stdout);
  /* every voxel in the unit cell */
  for (jz=0; jz<ngz; jz++) {
    /* fractional coordinates */
    zf = (double)jz/(double)ngz;
    for (jy=0; jy<ngy; jy++) {
      yf = (double)jy/(double)ngy;
      for (jx=0; jx<ngx; jx++) {
        xf = (double)jx/(double)ngx;
        /* cartesian coordinates */
        voxcenter_x = xf*gridin->cella_x+yf*gridin->cellb_x+zf*gridin->cellc_x;
        voxcenter_y = xf*gridin->cella_y+yf*gridin->cellb_y+zf*gridin->cellc_y;
        voxcenter_z = xf*gridin->cella_z+yf*gridin->cellb_z+zf*gridin->cellc_z;
        for (i=0; i<NIONMAX; i++) check[i] = 0;
        map->neighcount[jx][jy][jz] = 0;
        map->neighcount2[jx][jy][jz] = 0;
        map->swj[jx][jy][jz] = 0.0;
        map->swjk[jx][jy][jz] = 0.0;
        wmax=100000.0;
        wmax2=100000.0;
        ngcount++;
        ngp = ngcount*100/voxtot;
        if (ngp!=ngp0) {
          printf("\r%d%%", ngp);
          fflush(stdout);
        }
        ngp0=ngp;
        /* every unit cell in the supercell */
        for (ka=-kam; ka<=kam; ka++) {
          for (kb=-kbm; kb<=kbm; kb++) {
            for (kc=-kcm; kc<=kcm; kc++) {
              /* every atom in cartesian coordinates */
              for (atom=0; atom<gridin->nion; atom++) {
                xc = gridin->xcart[atom]+ka*gridin->cella_x+kb*gridin->cellb_x+kc*gridin->cellc_x;
                yc = gridin->ycart[atom]+ka*gridin->cella_y+kb*gridin->cellb_y+kc*gridin->cellc_y;
                zc = gridin->zcart[atom]+ka*gridin->cella_z+kb*gridin->cellb_z+kc*gridin->cellc_z;
                dist = sqrt((voxcenter_x-xc)*(voxcenter_x-xc)+(voxcenter_y-yc)*(voxcenter_y-yc)+
                  (voxcenter_z-zc)*(voxcenter_z-zc));
                if (dist<R_MAX) {  /* if the voxel is close enough to the translated atom */

                  index = map->neighcount2[jx][jy][jz];
                  map->ionmap[index][jx][jy][jz] = ((ka+3)<<13)+((kb+3)<<10)+((kc+3)<<7)+atom;
                  map->wj[index][jx][jy][jz] =  dist;
                  if (map->wj[index][jx][jy][jz]==-1000.0) return 1;
                  map->neighcount2[jx][jy][jz]++;
                  if (map->neighcount2[jx][jy][jz]==NEQVOX) {
                    printf("\n  BAD NEWS: The number of nearby atoms exceeds %d!\n", NEQVOX);
                    fprintf(cplog, "\nTerminated because the number of atoms near voxel %d %d %d is larger than %d\n",
                      jx, jy, jz, NEQVOX);
                    fprintf(cplog, "Suggestion: increase NEQVOX or decrease R_MAX and recompile\n");
                    errcount++;
                    return 2;
                  }
                }
              }
            }
          }
        }
        for (i=0; i<map->neighcount2[jx][jy][jz]; i++) {
          for (j=0; j<map->neighcount2[jx][jy][jz]; j++) {
              weightmatrix[i][j]=0.0;
          }
        }
        for (i=0; i<map->neighcount2[jx][jy][jz]; i++) {
          if(wmax>map->wj[i][jx][jy][jz]) {
            wmax2 = wmax;  /* second highest weight */
            wmax = map->wj[i][jx][jy][jz];  /* highest weight */
          } else if (wmax2>map->wj[i][jx][jy][jz]) wmax2 = map->wj[i][jx][jy][jz];
          /* temporary value holders to prevent overwriting */
          ka3[i] = (map->ionmap[i][jx][jy][jz]>>13&7)-3;
          kb3[i] = (map->ionmap[i][jx][jy][jz]>>10&7)-3;
          kc3[i] = (map->ionmap[i][jx][jy][jz]>>7&7)-3;
          atom3[i] = map->ionmap[i][jx][jy][jz]&127;
          wj_temp[i] = map->wj[i][jx][jy][jz];
        }
        weightmatrix_max = 100000.0;
        for (i=0; i<map->neighcount2[jx][jy][jz]; i++) {
           for (j=0; j<map->neighcount2[jx][jy][jz]; j++) {
            if(j!=i) {
                  weightmatrix[i][j] = wj_temp[i]*wj_temp[j];
                  if(weightmatrix[i][j] < weightmatrix_max) weightmatrix_max =weightmatrix[i][j];
            }
           }
         }
         wmax=100000.0;
         wmax2=100000.0;
        for (i=0; i<map->neighcount2[jx][jy][jz]; i++) {
                 check[i]=0;
        }
        for (i=0; i<map->neighcount2[jx][jy][jz]; i++) {
          for (j=i+1; j<map->neighcount2[jx][jy][jz]; j++) {
            /* if the weight product of the two atoms is within tolerance */
            if (weightmatrix[i][j]<=weightmatrix_max*(1.0+tolerance)) {
              ka1=ka3[i];
              kb1=kb3[i];
              kc1=kc3[i];
              atom=atom3[i];
              ka2=ka3[j];
              kb2=kb3[j];
              kc2=kc3[j];
              atom2=atom3[j];
              if (check[i]==0) {
                check[i] = 1;
                index = map->neighcount[jx][jy][jz];
                neighborindex[index]=i;
                map->ionmap[index][jx][jy][jz] = ((ka1+3)<<13)+((kb1+3)<<10)+((kc1+3)<<7)+atom;
                map->swj[jx][jy][jz] += 1.0;
                map->neighcount[jx][jy][jz]++;
              }
              if (check[j]==0) {
                check[j] = 1;
                index = map->neighcount[jx][jy][jz];
                neighborindex[index]=j;
                map->ionmap[index][jx][jy][jz] = ((ka2+3)<<13)+((kb2+3)<<10)+((kc2+3)<<7)+atom2;
                map->swj[jx][jy][jz] += 1.0;
                map->neighcount[jx][jy][jz]++;
              }
              map->swjk[jx][jy][jz] += 1.0;
              xc = gridin->xcart[atom]+ka1*gridin->cella_x+kb1*gridin->cellb_x+kc1*gridin->cellc_x;
              yc = gridin->ycart[atom]+ka1*gridin->cella_y+kb1*gridin->cellb_y+kc1*gridin->cellc_y;
              zc = gridin->zcart[atom]+ka1*gridin->cella_z+kb1*gridin->cellb_z+kc1*gridin->cellc_z;
              dist = sqrt((voxcenter_x-xc)*(voxcenter_x-xc)+(voxcenter_y-yc)*(voxcenter_y-yc)+
                (voxcenter_z-zc)*(voxcenter_z-zc));
              if (dist>dmax) dmax = dist;
              xc = gridin->xcart[atom2]+ka2*gridin->cella_x+kb2*gridin->cellb_x+kc2*gridin->cellc_x;
              yc = gridin->ycart[atom2]+ka2*gridin->cella_y+kb2*gridin->cellb_y+kc2*gridin->cellc_y;
              zc = gridin->zcart[atom2]+ka2*gridin->cella_z+kb2*gridin->cellb_z+kc2*gridin->cellc_z;
              dist = sqrt((voxcenter_x-xc)*(voxcenter_x-xc)+(voxcenter_y-yc)*(voxcenter_y-yc)+
                (voxcenter_z-zc)*(voxcenter_z-zc));
              if (dist>dmax) dmax = dist;
            }
          }
        }
        contactcounter=0;
        for (i=0; i<map->neighcount[jx][jy][jz]; i++) {
          for (j=i+1; j<map->neighcount[jx][jy][jz]; j++) {
            /* if the weight product of the two atoms is within tolerance */
            map->neighkey[contactcounter][jx][jy][jz]=0;
            index = neighborindex[i];
            index2 = neighborindex[j];
            if (weightmatrix[index][index2]<=weightmatrix_max*(1.0+tolerance)) {
              map->neighkey[contactcounter][jx][jy][jz]=1;
            }
            contactcounter++;
          }
        }
        if(map->neighcount[jx][jy][jz] > 100) {
            printf("Voxel %d %d %d reputedly has %d neighbors!  Please fix this.\n",jx,jy,jz);
            for(j=0; j< map->neighcount2[jx][jy][jz];j++) {
                     printf("%d %lf\n",j+1,wj_temp[j]);
            }
            exit(0);
        }
      }
    }
  }
  fprintf(cplog, "Longest contact within %.2lf bohr was %.2lf bohr away\n", R_MAX, dmax);
  fprintf(cplog, "Adjust R_MAX and recompile to improve accuracy/speed as necessary\n", R_MAX, dmax);
  printf(" Finished\n");
  return 0;
}

int CoordSearch(struct CrystData * gridin, struct ContactVol * map) {
  /* called by: main */
  /* calls: Getwj */
  int atom=0, atom2=0, atom3[NEQVOX], index=0,index2, i=0, j=0, jx=0, jy=0, jz=0;
  int ka1=0, kb1=0, kc1=0, ka2=0, kb2=0, kc2=0, ka3[NEQVOX], kb3[NEQVOX], kc3[NEQVOX];
  int ka=0, kb=0, kc=0, check[NEQVOX], ngcount=0, ngp=0, ngp0=0, voxtot=ngx*ngy*ngz;
  double dist=0.0, dmax=0.0, voxcenter_x=0.0, voxcenter_y=0.0, voxcenter_z=0.0;
  double wj_temp[NEQVOX], wmax=0.0, wmax2=0.0, xf=0.0, yf=0.0, zf=0.0, xc=0.0, yc=0.0, zc=0.0;
  double weightmatrix[NEQVOX][NEQVOX];
  int neighborindex[NEQVOX];
  double weightmatrix_max;
  int contactcounter;
  printf("0%%");
  fflush(stdout);
  /* every voxel in the unit cell */
  for (jz=0; jz<ngz; jz++) {
    /* fractional coordinates */
    zf = (double)jz/(double)ngz;
    for (jy=0; jy<ngy; jy++) {
      yf = (double)jy/(double)ngy;
      for (jx=0; jx<ngx; jx++) {
        xf = (double)jx/(double)ngx;
        /* cartesian coordinates */
        voxcenter_x = xf*gridin->cella_x+yf*gridin->cellb_x+zf*gridin->cellc_x;
        voxcenter_y = xf*gridin->cella_y+yf*gridin->cellb_y+zf*gridin->cellc_y;
        voxcenter_z = xf*gridin->cella_z+yf*gridin->cellb_z+zf*gridin->cellc_z;
        for (i=0; i<NIONMAX; i++) check[i] = 0;
        map->neighcount[jx][jy][jz] = 0;
        map->neighcount2[jx][jy][jz] = 0;
        map->swj[jx][jy][jz] = 0.0;
        map->swjk[jx][jy][jz] = 0.0;
        wmax=0.0;
        wmax2=0.0;
        ngcount++;
        ngp = ngcount*100/voxtot;
        if (ngp!=ngp0) {
          printf("\r%d%%", ngp);
          fflush(stdout);
        }
        ngp0=ngp;
        /* every unit cell in the supercell */
        for (ka=-kam; ka<=kam; ka++) {
          for (kb=-kbm; kb<=kbm; kb++) {
            for (kc=-kcm; kc<=kcm; kc++) {
              /* every atom in cartesian coordinates */
              for (atom=0; atom<gridin->nion; atom++) {
                xc = gridin->xcart[atom]+ka*gridin->cella_x+kb*gridin->cellb_x+kc*gridin->cellc_x;
                yc = gridin->ycart[atom]+ka*gridin->cella_y+kb*gridin->cellb_y+kc*gridin->cellc_y;
                zc = gridin->zcart[atom]+ka*gridin->cella_z+kb*gridin->cellb_z+kc*gridin->cellc_z;
                dist = sqrt((voxcenter_x-xc)*(voxcenter_x-xc)+(voxcenter_y-yc)*(voxcenter_y-yc)+
                  (voxcenter_z-zc)*(voxcenter_z-zc));
                if (dist<R_MAX) {  /* if the voxel is close enough to the translated atom */

                  index = map->neighcount2[jx][jy][jz];
                  map->ionmap[index][jx][jy][jz] = ((ka+3)<<13)+((kb+3)<<10)+((kc+3)<<7)+atom;
                  if (scheme==1) {  /* Hirshfeld weight */
                    map->wj[index][jx][jy][jz] = Getwj(atom, dist);
                    if (map->wj[index][jx][jy][jz]==-1000.0) return 1;
                  } else if (scheme==2) map->wj[index][jx][jy][jz] = R_MAX-dist+gridin->corerad[atom];  /* distance weight */
                  map->neighcount2[jx][jy][jz]++;
                  if (map->neighcount2[jx][jy][jz]==NEQVOX) {
                    printf("\n  BAD NEWS: The number of nearby atoms exceeds %d!\n", NEQVOX);
                    fprintf(cplog, "\nTerminated because the number of atoms near voxel %d %d %d is larger than %d\n",
                      jx, jy, jz, NEQVOX);
                    fprintf(cplog, "Suggestion: increase NEQVOX or decrease R_MAX and recompile\n");
                    errcount++;
                    return 2;
                  }
                }
              }
            }
          }
        }
        for (i=0; i<map->neighcount2[jx][jy][jz]; i++) {
          for (j=0; j<map->neighcount2[jx][jy][jz]; j++) {
              weightmatrix[i][j]=0.0;
          }
        }
        for (i=0; i<map->neighcount2[jx][jy][jz]; i++) {
          if(wmax<map->wj[i][jx][jy][jz]) {
            wmax2 = wmax;  /* second highest weight */
            wmax = map->wj[i][jx][jy][jz];  /* highest weight */
          } else if (wmax2<map->wj[i][jx][jy][jz]) wmax2 = map->wj[i][jx][jy][jz];
          /* temporary value holders to prevent overwriting */
          ka3[i] = (map->ionmap[i][jx][jy][jz]>>13&7)-3;
          kb3[i] = (map->ionmap[i][jx][jy][jz]>>10&7)-3;
          kc3[i] = (map->ionmap[i][jx][jy][jz]>>7&7)-3;
          atom3[i] = map->ionmap[i][jx][jy][jz]&127;
          wj_temp[i] = map->wj[i][jx][jy][jz];
        }
        weightmatrix_max = -100.0;
        for (i=0; i<map->neighcount2[jx][jy][jz]; i++) {
           for (j=0; j<map->neighcount2[jx][jy][jz]; j++) {
            if(j!=i) {
                  weightmatrix[i][j] = wj_temp[i]*wj_temp[j];
                  if(weightmatrix[i][j] > weightmatrix_max) weightmatrix_max =weightmatrix[i][j];
            }
           }
         }
         wmax=0.0;
         wmax2=0.0;
        for (i=0; i<map->neighcount2[jx][jy][jz]; i++) {
                 check[i]=0;
        }
        for (i=0; i<map->neighcount2[jx][jy][jz]; i++) {
          for (j=i+1; j<map->neighcount2[jx][jy][jz]; j++) {
            /* if the weight product of the two atoms is within tolerance */
            if (weightmatrix[i][j]>=weightmatrix_max*(1.0-tolerance)) {
              ka1=ka3[i];
              kb1=kb3[i];
              kc1=kc3[i];
              atom=atom3[i];
              ka2=ka3[j];
              kb2=kb3[j];
              kc2=kc3[j];
              atom2=atom3[j];
              if (check[i]==0) {
                check[i] = 1;
                index = map->neighcount[jx][jy][jz];
                neighborindex[index]=i;
                map->ionmap[index][jx][jy][jz] = ((ka1+3)<<13)+((kb1+3)<<10)+((kc1+3)<<7)+atom;
                map->swj[jx][jy][jz] += 1.0;
                map->neighcount[jx][jy][jz]++;
              }
              if (check[j]==0) {
                check[j] = 1;
                index = map->neighcount[jx][jy][jz];
                neighborindex[index]=j;
                map->ionmap[index][jx][jy][jz] = ((ka2+3)<<13)+((kb2+3)<<10)+((kc2+3)<<7)+atom2;
                map->swj[jx][jy][jz] += 1.0;
                map->neighcount[jx][jy][jz]++;
              }
              map->swjk[jx][jy][jz] += 1.0;
              xc = gridin->xcart[atom]+ka1*gridin->cella_x+kb1*gridin->cellb_x+kc1*gridin->cellc_x;
              yc = gridin->ycart[atom]+ka1*gridin->cella_y+kb1*gridin->cellb_y+kc1*gridin->cellc_y;
              zc = gridin->zcart[atom]+ka1*gridin->cella_z+kb1*gridin->cellb_z+kc1*gridin->cellc_z;
              dist = sqrt((voxcenter_x-xc)*(voxcenter_x-xc)+(voxcenter_y-yc)*(voxcenter_y-yc)+
                (voxcenter_z-zc)*(voxcenter_z-zc));
              if (dist>dmax) dmax = dist;
              xc = gridin->xcart[atom2]+ka2*gridin->cella_x+kb2*gridin->cellb_x+kc2*gridin->cellc_x;
              yc = gridin->ycart[atom2]+ka2*gridin->cella_y+kb2*gridin->cellb_y+kc2*gridin->cellc_y;
              zc = gridin->zcart[atom2]+ka2*gridin->cella_z+kb2*gridin->cellb_z+kc2*gridin->cellc_z;
              dist = sqrt((voxcenter_x-xc)*(voxcenter_x-xc)+(voxcenter_y-yc)*(voxcenter_y-yc)+
                (voxcenter_z-zc)*(voxcenter_z-zc));
              if (dist>dmax) dmax = dist;
            }
          }
        }
        contactcounter=0;
        for (i=0; i<map->neighcount[jx][jy][jz]; i++) {
          for (j=i+1; j<map->neighcount[jx][jy][jz]; j++) {
            /* if the weight product of the two atoms is within tolerance */
            map->neighkey[contactcounter][jx][jy][jz]=0;
            index = neighborindex[i];
            index2 = neighborindex[j];
            if (weightmatrix[index][index2]>=weightmatrix_max*(1.0-tolerance)) {
              map->neighkey[contactcounter][jx][jy][jz]=1;
            }
            contactcounter++;
          }
        }
        if(map->neighcount[jx][jy][jz] > 100) {
            printf("Voxel %d %d %d reputedly has %d neighbors!  Please fix this.\n",jx,jy,jz);
            for(j=0; j< map->neighcount2[jx][jy][jz];j++) {
                     printf("%d %lf\n",j+1,wj_temp[j]);
            }
            exit(0);
        }
      }
    }
  }
  fprintf(cplog, "Longest contact within %.2lf bohr was %.2lf bohr away\n", R_MAX, dmax);
  fprintf(cplog, "Adjust R_MAX and recompile to improve accuracy/speed as necessary\n", R_MAX, dmax);
  printf(" Finished\n");
  return 0;
}

int AssignContactHirsh(struct ContactVol * map, struct CrystData * gridin) {
  /* called by: main */
  /* calls: none */
  int atom=0, atom2=0, index=0, ngcount=0, ngp=0, ngp0=0, voxtot=ngx*ngy*ngz;
  int i=0, j=0, k=0, jx=0, jy=0, jz=0, ka=0, kb=0, kc=0, ka1=0, kb1=0, kc1=0, ka2=0, kb2=0, kc2=0;
  double coeff=0.0, denom=0.0, numer=0.0;
  int contactcounter;
  for (i=0; i<7; i++) {
    for (j=0; j<7; j++) {
      for (k=0; k<7; k++) {
        for (atom=0; atom<gridin->nion; atom++) {
          for (atom2=0; atom2<gridin->nion; atom2++) {
            map->count[i][j][k][atom][atom2] = 0.0;
            map->total[i][j][k][atom][atom2] = 0.0;
          }
        }
      }
    }
  }
  for(atom=0;atom<gridin->nion; atom++) {
      gridin->voxcount[atom]=0.0;
  }
  printf("0%%");
  fflush(stdout);
  for (jz=0; jz<ngz; jz++) {
    for (jy=0; jy<ngy; jy++) {
      for (jx=0; jx<ngx; jx++) {
        ngcount++;
        ngp = ngcount*100/voxtot;
        if (ngp!=ngp0) {
          printf("\r%d%%", ngp);
          fflush(stdout);
        }
        ngp0 = ngp;
        contactcounter=0;
        for (i=0; i<map->neighcount[jx][jy][jz]; i++) {
          ka1 = (map->ionmap[i][jx][jy][jz]>>13&7)-3;
          kb1 = (map->ionmap[i][jx][jy][jz]>>10&7)-3;
          kc1 = (map->ionmap[i][jx][jy][jz]>>7&7)-3;
          atom = map->ionmap[i][jx][jy][jz]&127;
          for (j=i+1; j<map->neighcount[jx][jy][jz]; j++) {
            index=map->neighkey[contactcounter][jx][jy][jz];
            if(index==1) {
             ka2 = (map->ionmap[j][jx][jy][jz]>>13&7)-3;
             kb2 = (map->ionmap[j][jx][jy][jz]>>10&7)-3;
             kc2 = (map->ionmap[j][jx][jy][jz]>>7&7)-3;
             atom2 = map->ionmap[j][jx][jy][jz]&127;
             /* weight of ions i*j compared to sum of all weight products at voxel jx, jy, jz */
             coeff = map->wj[j][jx][jy][jz];
             if(map->neighcount2[jx][jy][jz]>1) coeff=1.0/map->swjk[jx][jy][jz];
             map->count[ka2-ka1+3][kb2-kb1+3][kc2-kc1+3][atom][atom2] += coeff;
             gridin->voxcount[atom]+=0.5*coeff;
             gridin->voxcount[atom2]+=0.5*coeff;
             map->total[ka2-ka1+3][kb2-kb1+3][kc2-kc1+3][atom][atom2] += coeff*gridin->grid[jx][jy][jz];
            }
            contactcounter++;
          }
        }
      }
    }
  }
  printf(" Finished\n");
  /* every pair of atoms */
  for (atom=0; atom<gridin->nion; atom++) {
    for (atom2=atom; atom2<gridin->nion; atom2++) {
      /* every cell in the supercell */
      for (ka=-kam; ka<=kam; ka++) {
        for (kb=-kbm; kb<=kbm; kb++) {
          for (kc=-kcm; kc<=kcm; kc++) {
            if (atom==atom2 && ka==0 && kb==0 && kc==0) continue;  /* skip self interaction */
            numer = map->total[ka+3][kb+3][kc+3][atom][atom2]+map->total[-ka+3][-kb+3][-kc+3][atom2][atom];
            denom = map->count[ka+3][kb+3][kc+3][atom][atom2]+map->count[-ka+3][-kb+3][-kc+3][atom2][atom];
            /* total CP contribution of every ion pair / total number of voxels within their influence */
            map->average[ka+3][kb+3][kc+3][atom][atom2] = (denom==0.0)? 0.0 : numer/denom;
            map->average[-ka+3][-kb+3][-kc+3][atom2][atom] = (denom==0.0)? 0.0 : numer/denom;
          }
        }
      }
    }
  }
  return 0;
}

int AssignContact(struct ContactVol * map, struct CrystData * gridin) {
  /* called by: main */
  /* calls: none */
  int atom=0, atom2=0, index=0, ngcount=0, ngp=0, ngp0=0, voxtot=ngx*ngy*ngz;
  int i=0, j=0, k=0, jx=0, jy=0, jz=0, ka=0, kb=0, kc=0, ka1=0, kb1=0, kc1=0, ka2=0, kb2=0, kc2=0;
  double coeff=0.0, denom=0.0, numer=0.0;
  int contactcounter;
  for (i=0; i<7; i++) {
    for (j=0; j<7; j++) {
      for (k=0; k<7; k++) {
        for (atom=0; atom<gridin->nion; atom++) {
          for (atom2=0; atom2<gridin->nion; atom2++) {
            map->count[i][j][k][atom][atom2] = 0.0;
            map->total[i][j][k][atom][atom2] = 0.0;
          }
        }
      }
    }
  }
  for(atom=0;atom<gridin->nion; atom++) {
      gridin->voxcount[atom]=0.0;
  }
  printf("0%%");
  fflush(stdout);
  for (jz=0; jz<ngz; jz++) {
    for (jy=0; jy<ngy; jy++) {
      for (jx=0; jx<ngx; jx++) {
        ngcount++;
        ngp = ngcount*100/voxtot;
        if (ngp!=ngp0) {
          printf("\r%d%%", ngp);
          fflush(stdout);
        }
        ngp0 = ngp;
        contactcounter=0;
        for (i=0; i<map->neighcount[jx][jy][jz]; i++) {
          ka1 = (map->ionmap[i][jx][jy][jz]>>13&7)-3;
          kb1 = (map->ionmap[i][jx][jy][jz]>>10&7)-3;
          kc1 = (map->ionmap[i][jx][jy][jz]>>7&7)-3;
          atom = map->ionmap[i][jx][jy][jz]&127;
          for (j=i+1; j<map->neighcount[jx][jy][jz]; j++) {
            index=map->neighkey[contactcounter][jx][jy][jz];
            if(index==1) {
             ka2 = (map->ionmap[j][jx][jy][jz]>>13&7)-3;
             kb2 = (map->ionmap[j][jx][jy][jz]>>10&7)-3;
             kc2 = (map->ionmap[j][jx][jy][jz]>>7&7)-3;
             atom2 = map->ionmap[j][jx][jy][jz]&127;
             /* weight of ions i*j compared to sum of all weight products at voxel jx, jy, jz */
             coeff = 1.0/map->swjk[jx][jy][jz];
             map->count[ka2-ka1+3][kb2-kb1+3][kc2-kc1+3][atom][atom2] += coeff;
             gridin->voxcount[atom]+=0.5*coeff;
             gridin->voxcount[atom2]+=0.5*coeff;
             map->total[ka2-ka1+3][kb2-kb1+3][kc2-kc1+3][atom][atom2] += coeff*gridin->grid[jx][jy][jz];
            }
            contactcounter++;
          }
        }
      }
    }
  }
  printf(" Finished\n");
  /* every pair of atoms */
  for (atom=0; atom<gridin->nion; atom++) {
    for (atom2=atom; atom2<gridin->nion; atom2++) {
      /* every cell in the supercell */
      for (ka=-kam; ka<=kam; ka++) {
        for (kb=-kbm; kb<=kbm; kb++) {
          for (kc=-kcm; kc<=kcm; kc++) {
            if (atom==atom2 && ka==0 && kb==0 && kc==0) continue;  /* skip self interaction */
            numer = map->total[ka+3][kb+3][kc+3][atom][atom2]+map->total[-ka+3][-kb+3][-kc+3][atom2][atom];
            denom = map->count[ka+3][kb+3][kc+3][atom][atom2]+map->count[-ka+3][-kb+3][-kc+3][atom2][atom];
            /* total CP contribution of every ion pair / total number of voxels within their influence */
            map->average[ka+3][kb+3][kc+3][atom][atom2] = (denom==0.0)? 0.0 : numer/denom;
            map->average[-ka+3][-kb+3][-kc+3][atom2][atom] = (denom==0.0)? 0.0 : numer/denom;
          }
        }
      }
    }
  }
  return 0;
}


int AverageAtom(struct ContactVol * map, struct CrystData * gridin, struct CrystData * gridout) {
  /* called by: main */
  /* calls: Cart2Sph, FixEdges */
  int i=0, j=0, jx=0, jy=0, jz=0, ka=0, kb=0, kc=0, ka2=0, kb2=0, kc2=0, l=0, m=0;
  int atom=0, atom2=0, index=0, ngcount=0, ngp=0, ngp0=0, voxtot=ngx*ngy*ngz;
  double coeff=0.0, cosmphi[L_MAX], sinmphi[L_MAX], tempcp=0.0;
  double costheta=0.0, voxel_r=0.0, voxel_theta=0.0, voxel_phi=0.0, Y00=0.0;
  double voxcenter_x=0.0, voxcenter_y=0.0, voxcenter_z=0.0, xf=0.0, yf=0.0, zf=0.0;
  int contactcounter;
  int contactcounts[NEQVOX];
  for (i=0; i<gridin->nion; i++) {
    gridout->voxcount[i] = 0;
    gridout->intCP[i] = 0;
    for (l=0; l<10; l++) {
      for (m=0; m<l+1; m++) {
        gridout->intYlm[i][l][2*m] = 0.0;
        gridout->intYlm[i][l][2*m+1] = 0.0;
      }
    }
  }
  Y00 = gsl_sf_legendre_sphPlm(0,0,0);
  printf("0%%");
  fflush(stdout);
  for (jz=0; jz<ngz; jz++) {
    zf = (double)jz/(double)ngz;
    for (jy=0; jy<ngy; jy++) {
      yf = (double)jy/(double)ngy;
      for (jx=0; jx<ngx; jx++) {
        xf = (double)jx/(double)ngx;
        for(i=0;i<map->neighcount2[jx][jy][jz];i++) contactcounts[i]=0;
        ngcount++;
        ngp = ngcount*100/voxtot;
        if (ngp!=ngp0) {
          printf("\r%d%%", ngp);
          fflush(stdout);
        }
        ngp0 = ngp;
        tempcp = 0.0;
        contactcounter=0;
        gridout->grid[jx][jy][jz] = tempcp/map->swj[jx][jy][jz];
        for (i=0; i<map->neighcount2[jx][jy][jz]; i++) {
          atom = map->ionmap[i][jx][jy][jz]&127;
          ka = (map->ionmap[i][jx][jy][jz]>>13&7)-3;
          kb = (map->ionmap[i][jx][jy][jz]>>10&7)-3;
          kc = (map->ionmap[i][jx][jy][jz]>>7&7)-3;
          voxcenter_x = (xf-ka)*gridin->cella_x+(yf-kb)*gridin->cellb_x+(zf-kc)*gridin->cellc_x;
          voxcenter_y = (xf-ka)*gridin->cella_y+(yf-kb)*gridin->cellb_y+(zf-kc)*gridin->cellc_y;
          voxcenter_z = (xf-ka)*gridin->cella_z+(yf-kb)*gridin->cellb_z+(zf-kc)*gridin->cellc_z;
          Cart2Sph(voxcenter_x-gridin->xcart[atom], voxcenter_y-gridin->ycart[atom],
            voxcenter_z-gridin->zcart[atom], &voxel_r, &voxel_theta, &voxel_phi);
          costheta = cos(voxel_theta);
          for (m=1; m<lmax+1; m++) {
            sinmphi[m] = sin(m*voxel_phi);
            cosmphi[m] = cos(m*voxel_phi);
          }
          /* coeff = 1/2*[S(wj)-wj]*wj/S(wjk); sum of weight products involving ion j / total weight product */
          coeff = 1.0/map->swj[jx][jy][jz];
          gridout->voxcount[atom] += coeff;  /* total number of voxels within atom's influence (divided by 2) */
          gridout->intCP[atom] += coeff*gridin->grid[jx][jy][jz];  /* integrated CP of atom */
          /* l=0, m=0 coefficient of spherical harmonics */
          gridout->intYlm[atom][0][0] += Y00*coeff*gridin->grid[jx][jy][jz];
          /* higher order spherical harmonic coefficients */
          if (voxel_r>0.0) {
            for (l=1; l<lmax+1; l++) {
              m = 0;
              gridout->intYlm[atom][l][m] += gsl_sf_legendre_sphPlm(l, m, costheta)*gridin->grid[jx][jy][jz]*coeff;
              for (m=1; m<l+1; m++) {
                tempcp = 1.4142135*gsl_sf_legendre_sphPlm(l, m, costheta)*gridin->grid[jx][jy][jz]*coeff;
                gridout->intYlm[atom][l][2*m-1] += tempcp*cosmphi[m];
                gridout->intYlm[atom][l][2*m] += tempcp*sinmphi[m];
              }
            }
          }

        }
      }
    }
  }
  gridout->volvox = gridin->volvox;
  FixEdges(gridout);
  printf(" Finished\n\n");
  return 0;
}


int AverageAtomFast(struct ContactVol * map, struct CrystData * gridin, struct CrystData * gridout) {
  /* called by: main */
  /* calls: Cart2Sph, FixEdges */
  int i=0, j=0, jx=0, jy=0, jz=0, ka=0, kb=0, kc=0, ka2=0, kb2=0, kc2=0, l=0, m=0;
  int atom=0, atom2=0, index=0, ngcount=0, ngp=0, ngp0=0, voxtot=ngx*ngy*ngz;
  double coeff=0.0, cosmphi[L_MAX], sinmphi[L_MAX], tempcp=0.0;
  double costheta=0.0, voxel_r=0.0, voxel_theta=0.0, voxel_phi=0.0, Y00=0.0;
  double voxcenter_x=0.0, voxcenter_y=0.0, voxcenter_z=0.0, xf=0.0, yf=0.0, zf=0.0;
  int contactcounter;
  int contactcounts[NEQVOX];
  for (i=0; i<gridin->nion; i++) {
    gridout->voxcount[i] = 0;
    gridout->intCP[i] = 0;
    for (l=0; l<10; l++) {
      for (m=0; m<l+1; m++) {
        gridout->intYlm[i][l][2*m] = 0.0;
        gridout->intYlm[i][l][2*m+1] = 0.0;
      }
    }
  }
  Y00 = gsl_sf_legendre_sphPlm(0,0,0);
  printf("0%%");
  fflush(stdout);
  for (jz=0; jz<ngz; jz++) {
    zf = (double)jz/(double)ngz;
    for (jy=0; jy<ngy; jy++) {
      yf = (double)jy/(double)ngy;
      for (jx=0; jx<ngx; jx++) {
        xf = (double)jx/(double)ngx;
        for(i=0;i<map->neighcount2[jx][jy][jz];i++) contactcounts[i]=0;
        ngcount++;
        ngp = ngcount*100/voxtot;
        if (ngp!=ngp0) {
          printf("\r%d%%", ngp);
          fflush(stdout);
        }
        ngp0 = ngp;
        tempcp = 0.0;
        contactcounter=0;
        gridout->grid[jx][jy][jz] = tempcp/map->swj[jx][jy][jz];
        for (i=0; i<map->neighcount2[jx][jy][jz]; i++) {
          atom = map->ionmap[i][jx][jy][jz]&127;
          coeff = 1.0/map->swj[jx][jy][jz];
          gridout->voxcount[atom] += coeff;  /* total number of voxels within atom's influence (divided by 2) */
          gridout->intCP[atom] += coeff*gridin->grid[jx][jy][jz];  /* integrated CP of atom */
          /* l=0, m=0 coefficient of spherical harmonics */
          gridout->intYlm[atom][0][0] += Y00*coeff*gridin->grid[jx][jy][jz];
          /* higher order spherical harmonic coefficients */
        }
      }
    }
  }
  gridout->volvox = gridin->volvox;
  FixEdges(gridout);
  printf(" Finished\n\n");
  return 0;
}


int AverageContactHirsh(struct ContactVol * map, struct CrystData * gridin, struct CrystData * gridout) {
  /* called by: main */
  /* calls: Cart2Sph, FixEdges */
  int i=0, j=0, jx=0, jy=0, jz=0, ka=0, kb=0, kc=0, ka2=0, kb2=0, kc2=0, l=0, m=0;
  int atom=0, atom2=0, index=0, ngcount=0, ngp=0, ngp0=0, voxtot=ngx*ngy*ngz;
  double coeff=0.0, cosmphi[L_MAX], sinmphi[L_MAX], tempcp=0,neighcp[NEQVOX];
  double costheta=0.0, voxel_r=0.0, voxel_theta=0.0, voxel_phi=0.0, Y00=0.0;
  double voxcenter_x=0.0, voxcenter_y=0.0, voxcenter_z=0.0, xf=0.0, yf=0.0, zf=0.0;
  int contactcounter;
  int contactcounts[NEQVOX];
  for (i=0; i<gridin->nion; i++) {
    gridout->voxcount[i] = 0;
    gridout->intCP[i] = 0;
    for (l=0; l<10; l++) {
      for (m=0; m<l+1; m++) {
        gridout->intYlm[i][l][2*m] = 0.0;
        gridout->intYlm[i][l][2*m+1] = 0.0;
      }
    }
  }
  Y00 = gsl_sf_legendre_sphPlm(0,0,0);
  printf("0%%");
  fflush(stdout);
  for (jz=0; jz<ngz; jz++) {
    zf = (double)jz/(double)ngz;
    for (jy=0; jy<ngy; jy++) {
      yf = (double)jy/(double)ngy;
      for (jx=0; jx<ngx; jx++) {
        xf = (double)jx/(double)ngx;
        for(i=0;i<map->neighcount[jx][jy][jz];i++) contactcounts[i]=0;
        ngcount++;
        ngp = ngcount*100/voxtot;
        if (ngp!=ngp0) {
          printf("\r%d%%", ngp);
          fflush(stdout);
        }
        ngp0 = ngp;
        tempcp = 0.0;
        contactcounter=0;
        neighcp[0]=0.0;
        for (i=0; i<map->neighcount[jx][jy][jz]; i++) {
          ka = (map->ionmap[i][jx][jy][jz]>>13&7)-3;
          kb = (map->ionmap[i][jx][jy][jz]>>10&7)-3;
          kc = (map->ionmap[i][jx][jy][jz]>>7&7)-3;
          atom = map->ionmap[i][jx][jy][jz]&127;
          for (j=i+1; j<map->neighcount[jx][jy][jz]; j++) {
            ka2 = (map->ionmap[j][jx][jy][jz]>>13&7)-3;
            kb2 = (map->ionmap[j][jx][jy][jz]>>10&7)-3;
            kc2 = (map->ionmap[j][jx][jy][jz]>>7&7)-3;
            atom2 = map->ionmap[j][jx][jy][jz]&127;
            index = map->neighkey[contactcounter][jx][jy][jz];
            /* weighted sum of CP from all ion pairs */
            if (index==1) {
                      neighcp[0]= 1.0*map->average[ka2-ka+3][kb2-kb+3][kc2-kc+3][atom][atom2];
                      contactcounts[i]++;
                      contactcounts[j]++;
                      voxcenter_x = (xf-ka2)*gridin->cella_x+(yf-kb2)*gridin->cellb_x+(zf-kc2)*gridin->cellc_x;
                      voxcenter_y = (xf-ka2)*gridin->cella_y+(yf-kb2)*gridin->cellb_y+(zf-kc2)*gridin->cellc_y;
                      voxcenter_z = (xf-ka2)*gridin->cella_z+(yf-kb2)*gridin->cellb_z+(zf-kc2)*gridin->cellc_z;
                      Cart2Sph(voxcenter_x-gridin->xcart[atom2], voxcenter_y-gridin->ycart[atom2],
                            voxcenter_z-gridin->zcart[atom2], &voxel_r, &voxel_theta, &voxel_phi);
                      costheta = cos(voxel_theta);
                      for (m=1; m<lmax+1; m++) {
                         sinmphi[m] = sin(m*voxel_phi);
                         cosmphi[m] = cos(m*voxel_phi);
                      }
                      coeff = 0.5*map->wj[j][jx][jy][jz];
                      if(map->neighcount2[jx][jy][jz]>1) coeff = 0.5/map->swjk[jx][jy][jz];
                      gridout->voxcount[atom2] += coeff;  /* total number of voxels within atom's influence (divided by 2) */
                      gridout->intCP[atom2] += coeff*neighcp[0];  /* integrated CP of atom */
                      /* l=0, m=0 coefficient of spherical harmonics */
                      gridout->intYlm[atom2][0][0] += Y00*coeff*neighcp[0];
                      /* higher order spherical harmonic coefficients */
                      if (voxel_r>0.0) {
                        for (l=1; l<lmax+1; l++) {
                           m = 0;
                           gridout->intYlm[atom2][l][m] += gsl_sf_legendre_sphPlm(l, m, costheta)*neighcp[0]*coeff;
                           for (m=1; m<l+1; m++) {
                                tempcp = 1.4142135*gsl_sf_legendre_sphPlm(l, m, costheta)*neighcp[0]*coeff;
                                gridout->intYlm[atom2][l][2*m-1] += tempcp*cosmphi[m];
                                gridout->intYlm[atom2][l][2*m] += tempcp*sinmphi[m];
                           }
                        }
                      }
                      voxcenter_x = (xf-ka)*gridin->cella_x+(yf-kb)*gridin->cellb_x+(zf-kc)*gridin->cellc_x;
                      voxcenter_y = (xf-ka)*gridin->cella_y+(yf-kb)*gridin->cellb_y+(zf-kc)*gridin->cellc_y;
                      voxcenter_z = (xf-ka)*gridin->cella_z+(yf-kb)*gridin->cellb_z+(zf-kc)*gridin->cellc_z;
                      Cart2Sph(voxcenter_x-gridin->xcart[atom], voxcenter_y-gridin->ycart[atom],
                            voxcenter_z-gridin->zcart[atom], &voxel_r, &voxel_theta, &voxel_phi);
                      costheta = cos(voxel_theta);
                      for (m=1; m<lmax+1; m++) {
                         sinmphi[m] = sin(m*voxel_phi);
                         cosmphi[m] = cos(m*voxel_phi);
                      }
                      coeff = 0.5*map->wj[j][jx][jy][jz];
                      if(map->neighcount2[jx][jy][jz]>1) coeff = 0.5/map->swjk[jx][jy][jz];
                      gridout->voxcount[atom] += coeff;  /* total number of voxels within atom's influence (divided by 2) */
                      gridout->intCP[atom] += coeff*neighcp[0];  /* integrated CP of atom */
                      /* l=0, m=0 coefficient of spherical harmonics */
                      gridout->intYlm[atom][0][0] += Y00*coeff*neighcp[0];
                      /* higher order spherical harmonic coefficients */
                      if (voxel_r>0.0) {
                        for (l=1; l<lmax+1; l++) {
                           m = 0;
                           gridout->intYlm[atom][l][m] += gsl_sf_legendre_sphPlm(l, m, costheta)*neighcp[0]*coeff;
                           for (m=1; m<l+1; m++) {
                                tempcp = 1.4142135*gsl_sf_legendre_sphPlm(l, m, costheta)*neighcp[0]*coeff;
                                gridout->intYlm[atom][l][2*m-1] += tempcp*cosmphi[m];
                                gridout->intYlm[atom][l][2*m] += tempcp*sinmphi[m];
                           }
                        }
                      }
            }
            contactcounter++;
          }
        }
        /* total CP from all ions pairs at voxel jx, jy, jz (normalized by total sum of weights) */
      }
    }
  }
  gridout->volvox = gridin->volvox;
  FixEdges(gridout);
  printf(" Finished\n\n");
  return 0;
}

int AverageContact(struct ContactVol * map, struct CrystData * gridin, struct CrystData * gridout) {
  /* called by: main */
  /* calls: Cart2Sph, FixEdges */
  int i=0, j=0, jx=0, jy=0, jz=0, ka=0, kb=0, kc=0, ka2=0, kb2=0, kc2=0, l=0, m=0;
  int atom=0, atom2=0, index=0, ngcount=0, ngp=0, ngp0=0, voxtot=ngx*ngy*ngz;
  double coeff=0.0, cosmphi[L_MAX], sinmphi[L_MAX], tempcp=0.0;
  double costheta=0.0, voxel_r=0.0, voxel_theta=0.0, voxel_phi=0.0, Y00=0.0;
  double voxcenter_x=0.0, voxcenter_y=0.0, voxcenter_z=0.0, xf=0.0, yf=0.0, zf=0.0;
  int contactcounter;
  int contactcounts[NEQVOX];
  for (i=0; i<gridin->nion; i++) {
    gridout->voxcount[i] = 0;
    gridout->intCP[i] = 0;
    for (l=0; l<10; l++) {
      for (m=0; m<l+1; m++) {
        gridout->intYlm[i][l][2*m] = 0.0;
        gridout->intYlm[i][l][2*m+1] = 0.0;
      }
    }
  }
  Y00 = gsl_sf_legendre_sphPlm(0,0,0);
  printf("0%%");
  fflush(stdout);
  for (jz=0; jz<ngz; jz++) {
    zf = (double)jz/(double)ngz;
    for (jy=0; jy<ngy; jy++) {
      yf = (double)jy/(double)ngy;
      for (jx=0; jx<ngx; jx++) {
        xf = (double)jx/(double)ngx;
        for(i=0;i<map->neighcount[jx][jy][jz];i++) contactcounts[i]=0;
        ngcount++;
        ngp = ngcount*100/voxtot;
        if (ngp!=ngp0) {
          printf("\r%d%%", ngp);
          fflush(stdout);
        }
        ngp0 = ngp;
        tempcp = 0.0;
        contactcounter=0;
        for (i=0; i<map->neighcount[jx][jy][jz]; i++) {
          ka = (map->ionmap[i][jx][jy][jz]>>13&7)-3;
          kb = (map->ionmap[i][jx][jy][jz]>>10&7)-3;
          kc = (map->ionmap[i][jx][jy][jz]>>7&7)-3;
          atom = map->ionmap[i][jx][jy][jz]&127;
          for (j=i+1; j<map->neighcount[jx][jy][jz]; j++) {
            ka2 = (map->ionmap[j][jx][jy][jz]>>13&7)-3;
            kb2 = (map->ionmap[j][jx][jy][jz]>>10&7)-3;
            kc2 = (map->ionmap[j][jx][jy][jz]>>7&7)-3;
            atom2 = map->ionmap[j][jx][jy][jz]&127;
            index = map->neighkey[contactcounter][jx][jy][jz];
            /* weighted sum of CP from all ion pairs */
            if (index==1) {
                      tempcp += 1.0*map->average[ka2-ka+3][kb2-kb+3][kc2-kc+3][atom][atom2];
                      contactcounts[i]++;
                      contactcounts[j]++;
            }
            contactcounter++;
          }
        }
        /* total CP from all ions pairs at voxel jx, jy, jz (normalized by total sum of weights) */
        gridout->grid[jx][jy][jz] = tempcp/map->swjk[jx][jy][jz];
        for (i=0; i<map->neighcount[jx][jy][jz]; i++) {
          atom = map->ionmap[i][jx][jy][jz]&127;
          ka = (map->ionmap[i][jx][jy][jz]>>13&7)-3;
          kb = (map->ionmap[i][jx][jy][jz]>>10&7)-3;
          kc = (map->ionmap[i][jx][jy][jz]>>7&7)-3;
          index = map->neighkey[i][jx][jy][jz];
          voxcenter_x = (xf-ka)*gridin->cella_x+(yf-kb)*gridin->cellb_x+(zf-kc)*gridin->cellc_x;
          voxcenter_y = (xf-ka)*gridin->cella_y+(yf-kb)*gridin->cellb_y+(zf-kc)*gridin->cellc_y;
          voxcenter_z = (xf-ka)*gridin->cella_z+(yf-kb)*gridin->cellb_z+(zf-kc)*gridin->cellc_z;
          Cart2Sph(voxcenter_x-gridin->xcart[atom], voxcenter_y-gridin->ycart[atom],
            voxcenter_z-gridin->zcart[atom], &voxel_r, &voxel_theta, &voxel_phi);
          costheta = cos(voxel_theta);
          for (m=1; m<lmax+1; m++) {
            sinmphi[m] = sin(m*voxel_phi);
            cosmphi[m] = cos(m*voxel_phi);
          }
          /* coeff = 1/2*[S(wj)-wj]*wj/S(wjk); sum of weight products involving ion j / total weight product */
          coeff = 0.5*contactcounts[i]/map->swjk[jx][jy][jz];
          gridout->voxcount[atom] += coeff;  /* total number of voxels within atom's influence (divided by 2) */
          gridout->intCP[atom] += coeff*gridout->grid[jx][jy][jz];  /* integrated CP of atom */
          /* l=0, m=0 coefficient of spherical harmonics */
          gridout->intYlm[atom][0][0] += Y00*coeff*gridout->grid[jx][jy][jz];
          /* higher order spherical harmonic coefficients */
          if (voxel_r>0.0) {
            for (l=1; l<lmax+1; l++) {
              m = 0;
              gridout->intYlm[atom][l][m] += gsl_sf_legendre_sphPlm(l, m, costheta)*gridout->grid[jx][jy][jz]*coeff;
              for (m=1; m<l+1; m++) {
                tempcp = 1.4142135*gsl_sf_legendre_sphPlm(l, m, costheta)*gridout->grid[jx][jy][jz]*coeff;
                gridout->intYlm[atom][l][2*m-1] += tempcp*cosmphi[m];
                gridout->intYlm[atom][l][2*m] += tempcp*sinmphi[m];
              }
            }
          }
        }
      }
    }
  }
  gridout->volvox = gridin->volvox;
  FixEdges(gridout);
  printf(" Finished\n\n");
  return 0;
}



int PrintAverage(struct CrystData * gridin, struct CrystData * gridout) {
  /* called by: main */
  /* calls: ElementName */
  char element[10];
  int atom=0;
  double average=0.0;
  for(atom=0; atom<gridin->nion; atom++) {
    average = gridout->intCP[atom]/gridout->voxcount[atom];
    ElementName(gridin->zatomic[atom], element);
    printf("  Atom %d (%s) Integration Results\n", atom+1, element);
    printf("  around point (%.2f, %.2f, %.2f) Angstrom\n", gridin->xcart[atom]*R_BOHR,
      gridin->ycart[atom]*R_BOHR, gridin->zcart[atom]*R_BOHR);
    printf("  average over: %9.0f voxels\n", gridout->voxcount[atom]);
    printf("  l=0 coefficient: %+.2e a.u.\n", gridout->intYlm[atom][0][0]/gridout->voxcount[atom]);
    printf("  net pressure:    %+.2e a.u.\n", average);
    printf("  %+22.2f     GPa\n\n", average*AU2GPA);
    fprintf(cplog, "Atom %d (%s) Integration Results\n", atom+1, element);
    fprintf(cplog, "around point (%.2f, %.2f, %.2f) Angstrom\n", gridin->xcart[atom]*R_BOHR,
      gridin->ycart[atom]*R_BOHR, gridin->zcart[atom]*R_BOHR);
    fprintf(cplog, "average over: %9.0f voxels\n", gridout->voxcount[atom]);
    fprintf(cplog, "l=0 coefficient: %+.2e a.u.\n", gridout->intYlm[atom][0][0]/gridout->voxcount[atom]);
    fprintf(cplog, "net pressure:    %+.2e a.u.\n", average);
    fprintf(cplog, "%+22.2f     GPa\n\n", average*AU2GPA);
  }
  return 0;
}


int PrintAverage2(struct CrystData * gridin, struct CrystData * gridout) {
  /* called by: main */
  /* calls: ElementName */
  char element[10];
  int atom=0;
  double average=0.0;
  for(atom=0; atom<gridin->nion; atom++) {
    average = gridout->intCP[atom]/gridout->voxcount[atom];
    ElementName(gridin->zatomic[atom], element);
    printf("  Atom %d (%s) ", atom+1, element);
    printf("  net pressure:   ");
    printf("  %+22.2f     GPa\n", average*AU2GPA);
  }
  return 0;
}


// Jonathan adding code //

int PrintAverage3(struct CrystData * gridin, struct CrystData * gridout) {
  /* called by: main */
  /* calls: ElementName */
  char element[10];
  int atom=0;
  double average=0.0;
  double avg=0.0;
  for(atom=0; atom<gridin->nion; atom++) {
    average = gridout->intCP[atom]/gridout->voxcount[atom];
    avg = average*AU2GPA;
    ps2[atom] = avg;
    ElementName(gridin->zatomic[atom], element);
    printf("  Atom %d (%s) ", atom+1, element);
    printf("  net pressure:   ");
    printf("  %+22.2f     GPa\n", average*AU2GPA);
  }
  return(0);
}

// Jonathan stop adding code //

int PrintCoeff(struct CrystData * gridin, struct CrystData * gridout) {
  /* called by: PrintResults */
  /* calls: none */
  char filename[STRMAX];
  int atom=0, l=0, m=0;
  FILE * fptr;
  strncpy(filename, cpoutname, STRMAX);
  strncat(filename, "-coeff", STRMAX);
  fptr = fopen(filename, "w");
  for(atom=0; atom<gridin->nion; atom++) {
    fprintf(fptr, "l_%dm_%d=  %20.14f\n", 0, 0, gridout->intYlm[atom][0][0]/gridout->voxcount[atom]);
    for(l=1; l<lmax+1; l++) {
      fprintf(fptr, "l_%dm_%dp=  %20.14f\n", l, 0, gridout->intYlm[atom][l][0]/gridout->voxcount[atom]);
      for(m=1; m<l+1; m++) {
        fprintf(fptr, "l_%dm_%dp=  %20.14f\n", l, m, gridout->intYlm[atom][l][2*m-1]/gridout->voxcount[atom]);
        fprintf(fptr, "l_%dm_%dm=  %20.14f\n", l, m, gridout->intYlm[atom][l][2*m]/gridout->voxcount[atom]);
      }
    }
  }
  fclose(fptr);
  return 0;
}


/* MAIN HELPER FUNCTIONS */

int SetOptions() {
  /* called by: main */
  /* calls: none */
  int ds_option=0;
  fprintf(cplog, "Selected custom options:\n");
  fprintf(cplog, "Mapped following energy terms: ");
  printf("  Map kinetic energy? [1=Yes] [2=Yes, Thomas-Fermi] ");
  scanf("%d", &mapkin);
  if (mapkin==1) fprintf(cplog, "kinetic ");
  else if (mapkin==2) fprintf(cplog, "Thomas-Fermi kinetic ");
  printf("  Kinetic energy mapping method? [1=Gradient-based] [2=Laplacian-based] ");
  scanf("%d", &mapkinoption);
  if (mapkinoption==1) fprintf(cplog, "Gradient-based ");
  else if (mapkinoption==2) fprintf(cplog, "Laplacian-based ");
  printf("  Map local energy?                    [1=Yes] ");
  scanf("%d", &maploc);
  if (maploc==1) fprintf(cplog, "local ");
  printf("  Map Hartree energy?                  [1=Yes] ");
  scanf("%d", &maphart);
  if (maphart==1) fprintf(cplog, "hartree ");
  printf("  Map exchange-correlation energy?     [1=Yes] ");
  scanf("%d", &mapxc);
  if (mapxc==1) fprintf(cplog, "exchange-correlation ");
  /*NEW IN CP_PACKAGE2*/
  printf("  Map localized Ewald energy?    [0=No] [1=in proportion to total density] [2=in proportion to semicore electron density] [3=in proportion to nonlocal energy] [4=homogeneously in atomic volumes]\n");
  scanf("%d", &E_Ewald_option);
  if (E_Ewald_option==1) fprintf(cplog, "localized Ewald (in proportion to total density) ");
  if (E_Ewald_option==2) fprintf(cplog, "localized Ewald (in proportion to semicore electron density) ");
  if (E_Ewald_option==3) fprintf(cplog, "localized Ewald (in proportion to nonlocal energy) ");
  if (E_Ewald_option==4) fprintf(cplog, "localized Ewald (homogeneously in atomic volumes) ");
  printf("  Map E_alpha? [0=No] [1=with Ewald energy]\n");
  scanf("%d", &E_alpha_option);
  if (E_alpha_option==1) fprintf(cplog, "E_alpha (with localized Ewald energy) ");
  printf("  Output Ewald and E_alpha maps (if generated?) [1=Yes] ");
  scanf("%d", &outputewaldalpha);
  printf("  Map nonlocal pseudopotential energy?  [1=Yes] [2=Yes, read in files by kpt]");
  scanf("%d", &mapnonloc);
  if (mapnonloc==1 || mapnonloc==2) fprintf(cplog, "nonlocal pseudopotential ");
  /**/
  printf("  Use core unwarping (recommended)?    [1=Yes] ");
  scanf("%d", &mapcore);
  if (mapcore==1) {
    printf("  Restore symmetry (recommended)?      [1=Yes] ");
    scanf("%d", &mapsym);
  }
  fprintf(cplog, "\n");
  if (mapsym==1 && mapcore==1) fprintf(cplog, "With hirshfeld-inspired core unwarping and symmetry restoration\n");
  else if (mapsym==1) fprintf(cplog, "With hirshfeld-inspired core unwarping\n");
  retrycp:
  printf("  Add unmapped pressure homogeneously? [1=Yes] ");
  scanf("%d", &rescp);
  if (rescp!=1) {
    printf("  Substract it homogeneously instead?  [2=Yes] ");
    scanf("%d", &rescp);
  }
  if (rescp!=1 && rescp!=2) {
    printf("  Ignore it instead (not recommended)? [3=Yes] ");
    scanf("%d", &rescp);
  }
  if (rescp!=1 && rescp!=2 && rescp!=3) {
    printf("  Invalid option. Please choose again\n");
    goto retrycp;
  }
  fprintf(cplog, "Choice to handle remaining pressure = %d\n", rescp);
  retryscheme:
  printf("  Partition scheme? [1=Hirshfeld] [2=Distance] ");
  scanf("%d", &scheme);
  if (scheme!=1 && scheme!=2) {
    printf("  Invalid choice. Please try again\n");
    goto retryscheme;
  }
  retrytolerance:
  printf("  Percent tolerance?               [1=Default] ");
  scanf("%lf", &tolerance);
  if (tolerance<0.0 || tolerance>100.0) {
    printf("  Tolerance out of bounds. Please try again\n");
    goto retrytolerance;
  } else tolerance = tolerance/100.0;
  fprintf(cplog, "Scheme=%d with l_max=%d and tolerance=%.2f%% was used\n", scheme, lmax, tolerance*100.0);
  retrylmax:
  printf("  Max l for spherical harmonics?   [6=Default] ");
  scanf("%d", &lmax);
  if (lmax<=0) {
    printf("  Invalid input. Please choose again\n");
    goto retrylmax;
  } else if (lmax>L_MAX) {
    printf("\n  BAD NEWS: Max l=%d is higher than expected!\n", lmax);
    fprintf(cplog, "\nTerminated because Max l=%d is >= than L_MAX=%d\n", lmax, L_MAX);
    fprintf(cplog, "Suggestion: increase the value of L_MAX and recompile\n");
    errcount++;
    return 1;
  }
  if (scheme==2) {
    printf("  Enter core 'bubble' radii?           [1=Yes] ");
    scanf("%d", &isradii);
  }
  printf("  Use datasets other than 1 through 3? [1=Yes] ");
  scanf("%d", &ds_option);
  if (ds_option==1) {
    printf("  Enter the number of the expanded dataset:    ");
    scanf("%d", &dshi);
    printf("  Enter the number of the equilibrium dataset: ");
    scanf("%d", &dseq);
    printf("  Enter the number of the contracted dataset:  ");
    scanf("%d", &dslo);
    fprintf(cplog, "Datasets DS%d, DS%d, DS%d were used\n", dshi, dseq, dslo);
  }
  printf("  Output potential maps as XSF files?  [1=Yes] ");
  scanf("%d", &printbin);
  printf("  Output total energy maps?            [1=Yes] ");
  scanf("%d", &printen);
  if (mapcore==1) {
    printf("  Output map of new voxel volumes?     [1=Yes] ");
    scanf("%d", &printhmap);
  }
  printf("  Output voxel weight maps?            [1=Yes] ");
  scanf("%d", &printvmap);
  if (printbin==1) fprintf(cplog, "XSF files from binary potential files were output\n");
  if (printen==1) fprintf(cplog, "Intermediate energy maps were output\n");
  if (printhmap==1) fprintf(cplog, "Interpolated voxel maps were output\n");
  if (printvmap==1) fprintf(cplog, "Voxel weight maps were output\n");
  fprintf(cplog, "\n");
  printf("\n");
  return 0;
}

int ReadAll() {
  /* called by: main */
  /* calls: Den2XSF, ReadXSF */
  char denfile[STRMAX];
  Den2XSF(abinitname, dshi, "DEN", &denhi, &denhi2);
  Den2XSF(abinitname, dseq, "DEN", &deneq, &deneq2);
  Den2XSF(abinitname, dslo, "DEN", &denlo, &denlo2);
  if (mapkin==1) {
    Den2XSF(abinitname, dshi, "KDEN", &kdenhi, &kdenhi2);
    Den2XSF(abinitname, dseq, "KDEN", &kdeneq, &kdeneq2);
    Den2XSF(abinitname, dslo, "KDEN", &kdenlo, &kdenlo2);
  }
  if (mapkinoption==2) {
    Den2XSF(abinitname, dshi, "LDEN", &ldenhi, &ldenhi2);
    Den2XSF(abinitname, dseq, "LDEN", &ldeneq, &ldeneq2);
    Den2XSF(abinitname, dslo, "LDEN", &ldenlo, &ldenlo2);
  }
  Den2XSF(abinitname, dshi, "POT", &pothi, &pothi2);
  Den2XSF(abinitname, dseq, "POT", &poteq, &poteq2);
  Den2XSF(abinitname, dslo, "POT", &potlo, &potlo2);
  Den2XSF(abinitname, dshi, "VHA", &vhahi, &vhahi2);
  Den2XSF(abinitname, dseq, "VHA", &vhaeq, &vhaeq2);
  Den2XSF(abinitname, dslo, "VHA", &vhalo, &vhalo2);
  Den2XSF(abinitname, dshi, "VHXC", &vhxchi, &vhxchi2);
  Den2XSF(abinitname, dseq, "VHXC", &vhxceq, &vhxceq2);
  Den2XSF(abinitname, dslo, "VHXC", &vhxclo, &vhxclo2);
  if (E_Ewald_option == 2) {
    snprintf(denfile, STRMAX, "%s_DS%d_%s.xsf", abinitname, dseq, "pDEN");
    ReadXSF(denfile, &densceq);
    printf("  Semicore density files read in\n");
  }

  fprintf(cplog, "\n");
  if (errcount!=0) return 1;
  volhi = denhi.volvox;  /* global */
  voleq = deneq.volvox;  /* global */
  vollo = denlo.volvox;  /* global */
  return 0;
}

int ErrorCheck(double vol1, double vol2, double vol3) {
  /* called by: main */
  /* calls: none */
  if (NIONMAX < deneq.nion) {
    printf("\n  BAD NEWS: Too many atoms!\n");
    fprintf(cplog, "\nTerminated because nion (%d) > NIONMAX (%d)\n", deneq.nion, NIONMAX);
    fprintf(cplog, "Suggestion: increase NIONMAX and recompile CPpackage\n");
    errcount++;
    return 1;
  } else if (vol1<=vol2 || vol2<=vol3) {
    printf("\n  BAD NEWS: Negative or no change in volume!\n");
    fprintf(cplog, "\nTerminated because change in volume is not positive\n");
    fprintf(cplog, "DS1 vol=%20.14f \nDS2 vol=%20.14f \nDS3 vol=%20.14f\n", vol1, vol2, vol3);
    fprintf(cplog, "Suggestion: check Abinit input file for discrepencies\n");
    errcount++;
    return 2;
  } else if (denhi.nion!=deneq.nion || deneq.nion!=denlo.nion) {
    printf("\n  BAD NEWS: Unequal number of atoms between Abinit datasets!\n");
    fprintf(cplog, "\nTerminated because number of atoms is not equal\n");
    fprintf(cplog, "DS1 nion=%d, DS2 nion=%d, DS3 nion=%d\n", denhi.nion, deneq.nion, denlo.nion);
    fprintf(cplog, "Suggestion: check your files for discrepencies\n");
    errcount++;
    return 3;
  }
  return 0;
}



int AllocDbl(struct CrystData * gridin) {
  /* called by: main */
  /* calls: none */
  int gridx=ngx+1, gridy=ngy+1, gridz=ngz+1, jx=0, jy=0;
  gridin->grid = (double***)malloc(gridx*sizeof(double**));
  for (jx=0; jx<gridx; jx++) {
    gridin->grid[jx] = (double**)malloc(gridy*sizeof(double*));
    for (jy=0; jy<gridy; jy++) {
      gridin->grid[jx][jy] = (double*)malloc(gridz*sizeof(double));
    }
  }
  return 0;
}

int AllocHirsh(struct HirshMap *map) {
  /* called by: main */
  /* calls: none */
  int gridx=ngx+1, gridy=ngy+1, gridz=ngz+1, i=0, jx=0, jy=0;
  for (i=0; i<NIONMAX; i++) {
    hmap.hirsh_weight[i] = (double***)malloc(gridx*sizeof(double**));
    for (jx=0; jx<gridx; jx++) {
      hmap.hirsh_weight[i][jx] = (double**)malloc(gridy*sizeof(double*));
      for (jy=0; jy<gridy; jy++) {
        hmap.hirsh_weight[i][jx][jy] = (double*)malloc(gridz*sizeof(double));
      }
    }
  }
  return 0;
}

int AllocInt(struct ContactVol * map) {
  /* called by: main */
  /* calls: none */
  int gridx=ngx+1, gridy=ngy+1, gridz=ngz+1, i=0, jx=0, jy=0;
  map->neighcount = (int***)malloc(gridx*sizeof(int**));
  map->neighcount2 = (int***)malloc(gridx*sizeof(int**));
  map->swj = (double***)malloc(gridx*sizeof(double**));
  map->swjk = (double***)malloc(gridx*sizeof(double**));;
  map->wmax = (double***)malloc(gridx*sizeof(double**));
  map->wmax2 = (double***)malloc(gridx*sizeof(double**));
  for (jx=0; jx<gridx; jx++) {
    map->neighcount[jx] = (int**)malloc(gridy*sizeof(int*));
    map->neighcount2[jx] = (int**)malloc(gridy*sizeof(int*));
    map->swj[jx] = (double**)malloc(gridy*sizeof(double*));
    map->swjk[jx] = (double**)malloc(gridy*sizeof(double*));
    map->wmax[jx] = (double**)malloc(gridy*sizeof(double*));
    map->wmax2[jx] = (double**)malloc(gridy*sizeof(double*));

    for(jy=0; jy<gridy; jy++) {
      map->neighcount[jx][jy] = (int*)malloc(gridz*sizeof(int));
      map->neighcount2[jx][jy] = (int*)malloc(gridz*sizeof(int));
      map->swj[jx][jy] = (double*)malloc(gridz*sizeof(double));
      map->swjk[jx][jy] = (double*)malloc(gridz*sizeof(double));
      map->wmax[jx][jy] = (double*)malloc(gridz*sizeof(double));
      map->wmax2[jx][jy] = (double*)malloc(gridz*sizeof(double));

    }
  }
  for(i=0; i<NEQVOX; i++) {
    map->neighkey[i] = (int***)malloc(gridx*sizeof(int**));
    map->ionmap[i] = (int***)malloc(gridx*sizeof(int**));
    map->wj[i] = (double***)malloc(gridx*sizeof(double**));
    for(jx=0; jx<gridx; jx++) {
      map->neighkey[i][jx] = (int**)malloc(gridy*sizeof(int*));
      map->ionmap[i][jx] = (int**)malloc(gridy*sizeof(int*));
      map->wj[i][jx] = (double**)malloc(gridy*sizeof(double*));
      for(jy=0; jy<gridy; jy++) {
        map->neighkey[i][jx][jy] = (int*)malloc(gridz*sizeof(int));
        map->ionmap[i][jx][jy] = (int*)malloc(gridz*sizeof(int));
        map->wj[i][jx][jy] = (double*)malloc(gridz*sizeof(double));
      }
    }
  }
  fprintf(cplog, "Allocated memory for integration maps\n");
  return 0;
}

int PrintResults(struct CrystData * gridin, struct CrystData * gridout) {
  /* called by: main */
  /* calls: ElementName, OutputWeight, PrintCoeff, OutputXSF */
  char element[10], element2[10],filename[STRMAX],str[STRMAX];
  int i=0,atom,l,m,j,k,count=0,check;
  int elementlist[120];
  int elementcount=0,foundit;
  int natoms;
  FILE * fptr;
  int stop,stop2[NIONMAX];
  double xf,yf,zf,xc,yc,zc;
  stop=0;

  strncpy(filename, cpoutname, STRMAX);
  strncat(filename, "-averaged.xsf", STRMAX);
  fptr = fopen(filename, "w");
  OutputXSF(fptr, gridin, gridout);
  fclose(fptr);
  PrintCoeff(gridin, gridout);


/* Jonathan adding/editing code */


  for (i=0; i<gridin->nion; i++) {
    ElementName(gridin->zatomic[i], element);
    foundit=0.0;
    for(j=0;j<elementcount;j++) {
       if(gridin->zatomic[i]==elementlist[j]) {
                foundit=1;
       }
    }
    if(foundit==0) {
         elementlist[elementcount]=gridin->zatomic[i];
         elementcount++;
    }
  }

  factorial = 0;
  count = 0;
  for(i=0; i<elementcount; i++) {
    ElementName(elementlist[i], element1);
    for(j=count; j<elementcount; j++) {
      ElementName(elementlist[j], element2);
      strcpy(element3, element1);
      strcat(element3, element2);
      strcpy(bonds[factorial],element3);
      factorial += 1;
    }
    count++;
  }
  fprintf(htmlfile,"<!DOCTYPE html>\n<html lang='en'>\n<head>\n\t<meta charset='utf-8'>\n");
  fprintf(htmlfile,"\t<title> %s </title>\n\t<link rel='stylesheet' href='style.css'>\n", cpoutname);
  fprintf(htmlfile,"</head>\n<body>\n\t<div class='header'>\n\t\t<h1> CP Analysis: %s </h1>\n", cpoutname);
  fprintf(htmlfile,"\t</div>\n\t<div class='topnav'></div>\n");
  fprintf(htmlfile,"\t<div class='columns'>\n\t\t<div class='column right'>\n");
  fprintf(htmlfile,"\t\t\t<div class='tab'>\n");
  fprintf(htmlfile,"\t\t\t\t<button class='tablinks' onclick=\"opentab(event, 'settings')\" id='defaultOpen'>Settings</button>\n");
  fprintf(htmlfile,"\t\t\t\t<button class='tablinks' onclick=\"opentab(event, 'details')\">Details</button>\n");
  fprintf(htmlfile,"\t\t\t</div>\n\t\t\t<div id='settings' class='tabcontent'>\n");
  fprintf(htmlfile,"\t\t\t\t<h4> CP Scale: </h4>\n\t\t\t\t<div class='leftcolumn'>\n");
  fprintf(htmlfile,"\t\t\t\t\t<p> New CP Scale: </p>\n\t\t\t\t\t<p> Current CP Scale: </p>\n\t\t\t\t</div>\n");
  fprintf(htmlfile,"\t\t\t\t<div class='rightcolumn'>\n\t\t\t\t\t<input type='number' min='0' max='5000.0' value ='2500.0' id='cpscale'><input type='button' value='Update' id='myRange'><br>\n");
  fprintf(htmlfile,"\t\t\t\t\t<b><span id='demo'></span></b><br>\n\t\t\t\t</div>\n\t\t\t\t<div class='centercolumn'><hr></div>\n\t\t\t\t<h4> Atom Settings: </h4>\n");
  for(i=0; i<elementcount; i++){
    ElementName(elementlist[i], element);
    fprintf(htmlfile,"\t\t\t\t<p><input class='atom_color' type='color' value='%s' id='color%d'>&nbsp %s Color: <b>(<span id='output%dr'></span>,<span id='output%dg'></span>,<span id='output%db'></span>)</b></p><br>\n", colors[i], i+1, element, i+1, i+1, i+1);
    fprintf(htmlfile,"\t\t\t\t<p><input class='checkbox' type='checkbox' id='atom%d' checked><input class='atom_slider' type='range' min='0' max='1' value='0.2' id='atom%dsize' step'0.01>&nbsp %s Size: <b><span id='size%datom'></span> &#8491;</b></p><br>\n", i+1, i+1, element, i+1);
  }
  fprintf(htmlfile,"\t\t\t\t<hr>\n\t\t\t\t<h4> Connection Settings: </h4>\n\t\t\t\t<div class='smalltab'>\n");
  fprintf(htmlfile,"\t\t\t\t\t<button class=\"smalltablinks_2\" onclick=\"opensmalltab_2(event, 'bond_lengths')\" id='defaultopen_2'>Cutoff Lengths</button>\n");
  fprintf(htmlfile,"\t\t\t\t\t<button class=\"smalltablinks_2\" onclick=\"opensmalltab_2(event, 'bond_radii')\">Cylinder Radii</button>\n");
  fprintf(htmlfile,"\t\t\t\t</div>\n\t\t\t\t<div id='bond_lengths' class='smalltabcontent_2'\n");
  for(i=0; i<factorial; i++){
    fprintf(htmlfile,"\t\t\t\t\t<p> %s Connection Maximum Length: <b><span id='bond%dlen'></span> &#8491;</b></p>\n", bonds[i], i+1);
    fprintf(htmlfile,"\t\t\t\t\t<input type='range' min='0.1' max='5' value='3.2' step='0.1' class='slider' id='bond%dlength; oninput='update_bond_%d()'>\n", i+1, i+1);
    fprintf(htmlfile,"\t\t\t\t\t<input class='checkbox' type='checkbox' id='bond%d' checked>\n", i+1);
  }
  fprintf(htmlfile,"\t\t\t\t</div>\n\t\t\t\t<div id='bond_radii' class='smalltabcontent_2'>\n");
  for(i=0; i<factorial; i++){
    fprintf(htmlfile,"\t\t\t\t\t<p> %s Connection Radius: <b><span id='bond%drad'></span> &#8491;</b></p>\n", bonds[i], i+1);
    fprintf(htmlfile,"\t\t\t\t\t<input type='range' min='0.01' max='0.10' value='0.05' step='0.01' class='slider' id='bond%dradius' oninput='update_bond_1()'>\n", i+1, i+1);
  }
  fprintf(htmlfile,"\t\t\t\t\t<script>\n\t\t\t\t\t\tfunction opensmalltab_2(evt, tabname) {\n\t\t\t\t\t\t\tvar i, tabcontent_2, tablinks_2;\n");
  fprintf(htmlfile,"\t\t\t\t\t\t\ttabcontent_2 = document.getElementsByClassName('smalltabcontent_2');\n\t\t\t\t\t\t\tfor (i = 0; i < tabcontent_2.length; i++) {\n");
  fprintf(htmlfile,"\t\t\t\t\t\t\t\ttabcontent_2[i].style.display = 'none';\n\t\t\t\t\t\t\t}\n\t\t\t\t\t\t\ttablinks_2 = document.getElementsByClassName('smalltablinks_2');\n");
  fprintf(htmlfile,"\t\t\t\t\t\t\tfor (i = 0; i < tablinks_2.length; i++) {\n\t\t\t\t\t\t\t\ttablinks_2[i].className = tablinks_2[i].className.replace(\" active\", \"\");\n");
  fprintf(htmlfile,"\t\t\t\t\t\t\t}\n\t\t\t\t\t\t\tdocument.getElementById(tabname).style.display = 'block';\n\t\t\t\t\t\t\tevt.currentTarget.className += ' active';\n");
  fprintf(htmlfile,"\t\t\t\t\t\t}\n\t\t\t\t\t\tdocument.getElementById('defaultopen_2').click();\n\t\t\t\t\t</script>\n\t\t\t\t\t</div>\n\t\t\t\t\t<br>\n\t\t\t\t\t<hr>\n\t\t\t\t\t<h4> Add Polyhedra: </h4>\n");
  fprintf(htmlfile,"\t\t\t\t<div class='smalltab'>\n");
  fprintf(htmlfile,"\t\t\t\t\t<button class='smalltablinks' onclick='opensmalltab(event, \"atomcentered\")' id='defaultopen'>Atom Centered</button>\n");
  fprintf(htmlfile,"\t\t\t\t\t<button class='smalltablinks' onclick='opensmalltab(event, \"noncentered\")'>Non Centered</button>\n");
  fprintf(htmlfile,"\t\t\t\t</div>\n\t\t\t\t<div id='atomcentered' class='smalltabcontent'>\n\t\t\t\t\t<p> Atom Centered </p>\n");
  fprintf(htmlfile,"\t\t\t\t\t<select id='centeratom'>\n\t\t\t\t\t\t<option value='0'> Center Atom </option>\n");
  for(i=0; i<gridin->nion; i++) {
    fprintf(htmlfile,"\t\t\t\t\t\t<option value='%d'> Atom %d </option>\n", i+1, i+1);
  }
  fprintf(htmlfile,"\t\t\t\t\t\t<option value='216'> All </option>\n");
  fprintf(htmlfile,"\t\t\t\t\t</select>\n\t\t\t\t\t<p> Polyhedron Radius: <b><span id='polyrad'></span> &#8491;</b></p>\n");
  fprintf(htmlfile,"\t\t\t\t\t<input type='range' min='0.1' max='5' value='3' step='0.1' class='slider' id='range'>\n");
  fprintf(htmlfile,"\t\t\t\t\t<p> Material Transparency: <b><span id='trans1out'></span></b></p>\n");
  fprintf(htmlfile,"\t\t\t\t\t<input type='range' min='0.1' max='1' value='0.75' step='0.05' class='slider' id='trans1'>\n");
  fprintf(htmlfile,"\t\t\t\t\t<p> Material Color: <b>(<span id='polyout1r'></span>,<span id='polyout1g'></span>,<span id='polyout1b'></span>)</b><input type='color' value='#00c800' id='polycolor1' style='width:90%;'></p>\n");
  fprintf(htmlfile,"\t\t\t\t\t<p><input value='Add Polyhedron' type='button' id='addpoly1'></p>\n");
  fprintf(htmlfile,"\t\t\t\t\t<p><input value='Remove Polyhedron' type='button' id='removepoly1'></p>\n\t\t\t\t</div>\n");
  fprintf(htmlfile,"\t\t\t\t\t<div id='noncentered' class='smalltabcontent'>\n\t\t\t\t\t<p> Non Centered </p>\n");
  if(gridin->nion > 8) {
    for(i=0; i<gridin->nion; i++) {
      fprintf(htmlfile,"\t\t\t\t\t<select id='polyatom%d'>\n", i+1);
      fprintf(htmlfile,"\t\t\t\t\t\t<option value='0'> Vertex %d </option>\n", i+1);
      for(j=0; j<gridin->nion; j++) {
        fprintf(htmlfile,"\t\t\t\t\t\t<option value='%d'> Atom %d </option>\n", j+1, j+1);
      }
      fprintf(htmlfile,"\t\t\t\t\t</select>\n");
    }
  }
  else {
    for(i=0; i<8; i++) {
      fprintf(htmlfile,"\t\t\t\t\t<select id='polyatom%d'>\n", i+1);
      fprintf(htmlfile,"\t\t\t\t\t\t<option value='0'> Vertex %d </option>\n", i+1);
      for(j=0; j<gridin->nion; j++) {
        fprintf(htmlfile,"\t\t\t\t\t\t<option value='%d'> Atom %d </option>\n", j+1, j+1);
      }
      fprintf(htmlfile,"\t\t\t\t\t</select>\n");
    }
  }
  fprintf(htmlfile,"\t\t\t\t\t<p> Material Transparency: <b><span id='trans2out'></span></b></p>\n");
  fprintf(htmlfile,"\t\t\t\t\t<input type='range' min='0.1' max='1' value='0.75' step='0.05' class='slider' id='trans2'>\n");
  fprintf(htmlfile,"\t\t\t\t\t<p> Material Color: <b>(<span id='polyout2r'></span>,<span id='polyout2g'></span>,<span id='polyout2b'></span>)</b><input type='color' value='#00c800' id='polycolor2' style='width:90%;'></p>\n");
  fprintf(htmlfile,"\t\t\t\t\t<p><input value='Add Polyhedron' type='button' id='addpoly2'></p>\n");
  fprintf(htmlfile,"\t\t\t\t\t<p><input value='Remove Polyhedron' type='button' id='removepoly2'></p>\n");
  fprintf(htmlfile,"\t\t\t\t\t<script>\n\t\t\t\t\t\tfunction opensmalltab(evt, tabname) {\n\t\t\t\t\t\t\tvar i, tabcontent, tablinks;\n");
  fprintf(htmlfile,"\t\t\t\t\t\t\ttabcontent = document.getElementsByClassName(\"smalltabcontent\");\n\t\t\t\t\t\t\tfor (i = 0; i < tabcontent.length; i++) {\n");
  fprintf(htmlfile,"\t\t\t\t\t\t\t\ttabcontent[i].style.display = 'none';\n\t\t\t\t\t\t\t}\n\t\t\t\t\t\t\ttablinks = document.getElementsByClassName(\"smalltablinks\");\n");
  fprintf(htmlfile,"\t\t\t\t\t\t\tfor (i = 0; i < tablinks.length; i++) {\n\t\t\t\t\t\t\t\ttablinks[i].className = tablinks[i].className.replace(\" active\", \"\");\n");
  fprintf(htmlfile,"\t\t\t\t\t\t\t}\n\t\t\t\t\t\t\tdocument.getElementById(tabname).style.display = 'block';\n\t\t\t\t\t\t\tevt.currentTarget.className += \" active\";\n");
  fprintf(htmlfile,"\t\t\t\t\t\t}\n\t\t\t\t\t\tdocument.getElementById(\"defaultopen\").click();\n\t\t\t\t\t</script>\n\t\t\t\t\t</div>\n");
  fprintf(htmlfile,"\t\t\t\t<br>\n\t\t\t\t<hr>\n\t\t\t\t<h4> ScaleBar Settings: </h4>\n\t\t\t\t<p> ScaleBar Postition X-Shift: <b><span id='scalebar_y_shift'></span></b></p>\n");
  fprintf(htmlfile,"\t\t\t\t<input type='range' class='slider' id='scalebar_y' min='-20' max='20' step='0.1' value='8'>\n");
  fprintf(htmlfile,"\t\t\t\t<p> ScaleBar Postition Y-Shift: <b><span id='scalebar_x_shift'></span></b></p>\n");
  fprintf(htmlfile,"\t\t\t\t<input type='range' class='slider' id='scalebar_x' min='-20' max='20' step='0.1' value='0'>\n");
  fprintf(htmlfile,"\t\t\t\t<p> Show ScaleBar: <input class='checkbox' type='checkbox' id='scalebar_checked'></p>\n");
  fprintf(htmlfile,"\t\t\t\t<p> ScaleBar Size: <input type='number' min='0' max='5000.0' value ='100' id='scalebar_scale'><input type='button' value='Update' id='scalebar_update_scale'><b><span id='scalebar_scale_out'></span> GPa </b></p>\n");
  fprintf(htmlfile,"\t\t\t\t<br>\n\t\t\t\t<hr>\n\t\t\t\t<h4> Template Settings: </h4>\n");
  for(i=0; i<ntemplates; i++) {
    fprintf(htmlfile,"\t\t\t\t<button class='button' type='button' id='temp%d'> Toggle Template %d </button>\n", i+1, i+1);
  }
  fprintf(htmlfile,"\t\t\t\t<button class='button' type='button' id='save_image'> Save Image </button>\n");
  fprintf(htmlfile,"\t\t\t\t<br><br>\n\t\t\t</div>\n");
  fprintf(javafile,"\nvar geo = [];\n");
  for (i=0; i<gridin->nion; i++) {
    ElementName(gridin->zatomic[i], element);
    fprintf(javafile,"geo[%d]=['%s', %20.14lf, %20.14lf, %20.14lf ];\n",i, element,gridin->xcart[i]*R_BOHR, gridin->ycart[i]*R_BOHR, gridin->zcart[i]*R_BOHR);
  }
  fprintf(javafile, "\n");
  if(ntemplates<1) {
    fprintf(javafile,"var template1 = [];\n");
    for (i=0; i<gridin->nion; i++) {
        ElementName(gridin->zatomic[i], element);
        fprintf(javafile,"template1[%d]=['%s', %20.14lf, %20.14lf, %20.14lf ];\n",i, element,gridin->xcart[i]*R_BOHR, gridin->ycart[i]*R_BOHR, gridin->zcart[i]*R_BOHR);
    }
    ntemplates=1;
    fprintf(javafile,"\n");
  }
  else {
    for(j=0;j<ntemplates;j++) {
      fprintf(javafile,"var template%d = [];\n",j+1);
      fptr = fopen(templatefiles[j], "r");
      stop=0;
      while(stop==0) {
        check=fscanf(fptr,"%s",str);
        if(check==EOF) {
          printf("Structure data not found in %s.\n",templatefiles[j]);
          exit(1);
        }
        natoms =0;
        if(strcmp(str,"NAME")==0) {
          FinishLine(fptr);
          while(stop==0) {
            check=fscanf(fptr,"%s %lf %lf %lf",str,&xf,&yf,&zf);
            xc = xf*gridin->cella_x*R_BOHR + yf*gridin->cellb_x*R_BOHR + zf*gridin->cellc_x*R_BOHR;
            yc = xf*gridin->cella_y*R_BOHR + yf*gridin->cellb_y*R_BOHR + zf*gridin->cellc_y*R_BOHR;
            zc = xf*gridin->cella_z*R_BOHR + yf*gridin->cellb_z*R_BOHR + zf*gridin->cellc_z*R_BOHR;
            if(check > 1) {
              fprintf(javafile,"template%d[%d]=['%s', %20.14lf, %20.14lf, %20.14lf ];\n",j+1,natoms, str,xc, yc, zc);
              natoms++;
            }
            else stop = 1;
          }
        }
      }
      fprintf(javafile,"\n");
    }
  }
  fprintf(htmlfile,"\t\t\t<div id='details' class='tabcontent'>\n");
  fprintf(htmlfile,"\t\t\t\t<p>\n");
  fprintf(htmlfile,"\t\t\t\tSystem Name: %s <br>\n", cpoutname);
  fprintf(htmlfile,"\t\t\t\tCalculation was completed on: %d <br>\n", calcdate);
  fprintf(htmlfile,"\t\t\t\t</p>\n");
  fprintf(htmlfile,"\t\t\t</div>\n");
  fprintf(htmlfile,"\t\t\t<script>\n");
  fprintf(htmlfile,"\t\t\t\tfunction opentab(evt, tabname) {\n");
  fprintf(htmlfile,"\t\t\t\t\tvar i, tabcontent, tablinks;\n");
  fprintf(htmlfile,"\t\t\t\t\ttabcontent = document.getElementsByClassName('tabcontent');\n");
  fprintf(htmlfile,"\t\t\t\t\tfor (i = 0; i < tabcontent.length; i++) {\n");
  fprintf(htmlfile,"\t\t\t\t\t\ttabcontent[i].style.display = 'none';\n");
  fprintf(htmlfile,"\t\t\t\t\t}\n");
  fprintf(htmlfile,"\t\t\t\t\ttablinks = document.getElementsByClassName('tablinks');\n");
  fprintf(htmlfile,"\t\t\t\t\tfor (i = 0; i < tablinks.length; i++) {\n");
  fprintf(htmlfile,"\t\t\t\t\t\ttablinks[i].className = tablinks[i].className.replace(' active', '');\n");
  fprintf(htmlfile,"\t\t\t\t\t}\n");
  fprintf(htmlfile,"\t\t\t\t\tdocument.getElementById(tabname).style.display = 'block';\n");
  fprintf(htmlfile,"\t\t\t\t\tevt.currentTarget.className += ' active';\n");
  fprintf(htmlfile,"\t\t\t\t}\n");
  fprintf(htmlfile,"\t\t\t\tdocument.getElementById('defaultOpen').click();\n");
  fprintf(htmlfile,"\t\t\t</script>\n\t\t</div>\n");
  fprintf(htmlfile,"\t\t<div class='column middle'>\n");
  fprintf(htmlfile,"\t\t\t<div id='canvasbox' class='canvasbox'></div>\n");
  fprintf(htmlfile,"\t\t\t<script src='three.js'></script>\n");
  fprintf(htmlfile,"\t\t\t<script src='TrackballControls.js'></script>\n");
  fprintf(htmlfile,"\t\t\t<script src='figuretoolweb.js'></script>\n");
  fprintf(htmlfile,"\t\t\t<script src='%s.js'></script>\n", cpoutname);
  fprintf(htmlfile,"\t\t\t<script src='create2.js'></script>\n");
  fprintf(htmlfile,"\t\t</div>\n\t</div>\n\t<div class='footer'>\n");
  fprintf(htmlfile,"\t\t<p> <i><a href='https://www2.chem.wisc.edu/~danny/group/'>FigureToolWeb</a>, Copyright &#169 2021 The Fredrickson Group, Department of Chemistry, University of Wisconsin - Madison</i> </p>\n");
  fprintf(htmlfile,"\t</div>\n</body>\n</html>\n");
  fclose(htmlfile);

  fprintf(javafile,"var cell = [];\n");
  fprintf(javafile, "cell[0]=[%20.14f,  %20.14f,  %20.14f];\n", gridin->cella_x*R_BOHR, gridin->cella_y*R_BOHR, gridin->cella_z*R_BOHR);
  fprintf(javafile, "cell[1]=[%20.14f,  %20.14f,  %20.14f];\n", gridin->cellb_x*R_BOHR, gridin->cellb_y*R_BOHR, gridin->cellb_z*R_BOHR);
  fprintf(javafile, "cell[2]=[%20.14f,  %20.14f,  %20.14f];\n", gridin->cellc_x*R_BOHR, gridin->cellc_y*R_BOHR, gridin->cellc_z*R_BOHR);
  fprintf(javafile,"\n");
  fprintf(javafile,"var CPcoeff = [];\n");
  for(atom=0; atom<gridin->nion; atom++) {
    fprintf(javafile, "CPcoeff[%d]= [%20.14f",atom,gridout->intYlm[atom][0][0]/gridout->voxcount[atom]);
    for(l=1; l<lmax+1; l++) {
      fprintf(javafile, ", %20.14f", gridout->intYlm[atom][l][0]/gridout->voxcount[atom]);
      for(m=1; m<l+1; m++) {
        fprintf(javafile, ", %20.14f", gridout->intYlm[atom][l][2*m-1]/gridout->voxcount[atom]);
        fprintf(javafile, ", %20.14f", gridout->intYlm[atom][l][2*m]/gridout->voxcount[atom]);
      }
    }
    fprintf(javafile,"];\n");
  }
  fprintf(javafile, "\n");
  for(i=0;i<elementcount;i++) {
    fprintf(javafile, "var color%d = document.getElementById('color%d');\n", i+1, i+1);
    fprintf(javafile, "var output%dr = document.getElementById('output%dr');\n", i+1, i+1);
    fprintf(javafile, "var output%dg = document.getElementById('output%dg');\n", i+1, i+1);
    fprintf(javafile, "var output%db = document.getElementById('output%db');\n", i+1, i+1);
    fprintf(javafile, "output%dr.innerHTML = hextorgb(color%d.value).r;\n", i+1, i+1);
    fprintf(javafile, "output%dg.innerHTML = hextorgb(color%d.value).g;\n", i+1, i+1);
    fprintf(javafile, "output%db.innerHTML = hextorgb(color%d.value).b;\n", i+1, i+1);
    fprintf(javafile, "var color%dout = { color: color%d.value };\n", i+1, i+1);
  }
  fprintf(javafile,"\n");
  for(i=0;i<elementcount;i++) {
    ElementName(elementlist[i], element);
    fprintf(javafile,"var %s_material = define_material(color%dout,1);\n", element, i+1);
  }
  fprintf(javafile,"var plus = define_material({ color: 'rgb(230,230,230)'},1);\n");
  fprintf(javafile,"var minus = define_material({ color: 'rgb(0,0,0)'},1);\n");
  fprintf(javafile,"\n");
  for(i=0;i<ntemplates;i++) {
    for(j=0;j<elementcount;j++) {
      ElementName(elementlist[j], element);
      fprintf(javafile,"var %s_atoms_template%d = drawmol_spheres(template%d, '%s',%s_material, 0.2);\n",element, i+1, i+1, element, element);
    }
    count = 0;
    for(j=0;j<elementcount;j++) {
      ElementName(elementlist[j], element1);
      for(k=count;k<elementcount;k++) {
        ElementName(elementlist[k], element2);
        strcpy(element3, element1);
        strcat(element3, element2);
        for(l=0;l<factorial;l++) {
          if(strstr(element3, bonds[l]) != NULL) {
            fprintf(javafile, "var %s_bonds_template%d = drawbonds_cylinder(template%d,'%s',%s_material,'%s',%s_material, 0.1, 3.2, 0.05);\n", element3, i+1, i+1, element1, element1, element2, element2);
          }
        }
      }
      count++;
    }
  }
  fprintf(javafile,"\n// var poly_geo_ref = ['num1', 'num2']; // Replace num1 and num2 with the geo[] indices you wish to add simultaneously, add as many sites as desired\n\n");
  for(i=0;i<elementcount;i++){
    fprintf(javafile, "var atom%dstat = 0;\n", i+1);
  }
  for(i=0;i<factorial;i++) {
    fprintf(javafile, "var bond%dstat = 0;\n", i+1);
  }
  for(i=0;i<ntemplates;i++) {
    fprintf(javafile, "var togg%dstat = 0;\n", i+1);
  }
  fprintf(javafile,"\n");
  for(i=0;i<elementcount;i++) {
    ElementName(elementlist[i], element);
    fprintf(javafile, "var atom%dsize = document.getElementById('atom%dsize');\n", i+1, i+1);
    fprintf(javafile, "var size%datom = document.getElementById('size%datom');\n", i+1, i+1);
    fprintf(javafile, "size%datom.innerHTML = atom%dsize.value;\n", i+1, i+1);
    fprintf(javafile, "var size%d = atom%dsize.value;\n", i+1, i+1);
    fprintf(javafile, "atom%dsize.oninput = function() {\n", i+1);
    fprintf(javafile, "\tsize%datom.innerHTML = this.value;\n", i+1);
    fprintf(javafile, "\tsize%d = this.value;\n", i+1);
    for(j=0;j<ntemplates;j++) {
      fprintf(javafile, "\tif (atom%dstat===0 && togg%dstat===0) {\n", i+1, j+1);
    	fprintf(javafile, "\t\tremove_objects(%s_atoms_template%d);\n", element, j+1);
    	fprintf(javafile, "\t\t%s_atoms_template%d = drawmol_spheres(template%d,'%s',%s_material, size%d);\n", element, j+1, j+1, element, element, i+1);
    	fprintf(javafile, "\t\tshow_objects(%s_atoms_template%d);\n\t}\n", element, j+1);
    	fprintf(javafile, "\telse {\n\t\tremove_objects(%s_atoms_template%d);\n", element, j+1);
    	fprintf(javafile, "\t\t%s_atoms_template%d = drawmol_spheres(template%d,'%s',%s_material, size%d);\n", element, j+1, j+1, element, element, i+1);
    	fprintf(javafile, "\t\tdontshow_objects(%s_atoms_template%d);\n\t}\n", element, j+1);
    }
    fprintf(javafile,"}\n\n");
  }
  count=1;
  count1=0;
  for(i=0;i<elementcount;i++) {
    ElementName(elementlist[i], element1);
    for(j=count1;j<elementcount;j++) {
      ElementName(elementlist[j], element2);
      fprintf(javafile, "var bond%dlength = document.getElementById('bond%dlength');\n", count, count);
      fprintf(javafile, "var bond%dlen = document.getElementById('bond%dlen');\n", count, count);
      fprintf(javafile, "bond%dlen.innerHTML = bond%dlength.value;\n", count, count);
      fprintf(javafile, "var len%d = bond%dlength.value;\n", count, count);
      fprintf(javafile, "var bond%dradius = document.getElementById('bond%dradius');\n", count, count);
      fprintf(javafile, "var bond%drad = document.getElementById('bond%drad');\n", count, count);
      fprintf(javafile, "bond%drad.innerHTML = bond%dradius.value;\n", count, count);
      fprintf(javafile, "var rad%d = bond%dradius.value;\n", count, count);
      fprintf(javafile, "function update_bond_%d() {\n", count);
      fprintf(javafile, "\tbond%dlen.innerHTML = document.getElementById('bond%dlength').value;\n", count, count);
      fprintf(javafile, "\tlen%d = document.getElementById('bond%dlength').value;\n", count, count);
      fprintf(javafile, "\tbond%drad.innerHTML = document.getElementById('bond%dradius').value;\n", count, count);
      fprintf(javafile, "\trad%d = document.getElementById('bond%dradius').value;\n", count, count);
      for(k=0;k<ntemplates;k++) {
        fprintf(javafile, "\tif (bond%dstat===0 && togg%dstat===0) {\n", count, k+1);
    	  fprintf(javafile, "\t\tremove_objects(%s%s_bonds_template%d);\n", element1, element2, k+1);
    	  fprintf(javafile, "\t\t%s%s_bonds_template%d = drawbonds_cylinder(template%d,'%s',%s_material,'%s',%s_material,0.1,len%d,rad%d);\n", element1, element2, k+1, k+1, element1, element1, element2, element2, count, count);
    	  fprintf(javafile, "\t\tshow_objects(%s%s_bonds_template%d);\n\t}\n", element1, element2, k+1);
    	  fprintf(javafile, "\telse {\n\t\tremove_objects(%s%s_bonds_template%d);\n", element1, element2, k+1);
    	  fprintf(javafile, "\t\t%s%s_bonds_template%d = drawbonds_cylinder(template%d,'%s',%s_material,'%s',%s_material,0.1,len%d,rad%d);\n", element1, element2, k+1, k+1, element1, element1, element2, element2, count, count);
    	  fprintf(javafile, "\t\tdontshow_objects(%s%s_bonds_template%d);\n\t}\n", element1, element2, k+1);
      }
      fprintf(javafile, "}\n\n");
      count++;
    }
    count1++;
  }

  for(i=0;i<elementcount;i++) {
    ElementName(elementlist[i], element);
    fprintf(javafile, "color%d.oninput = function() {\n", i+1);
    fprintf(javafile, "\toutput%dr.innerHTML = hextorgb(this.value).r;\n", i+1);
    fprintf(javafile, "\toutput%dg.innerHTML = hextorgb(this.value).g;\n", i+1);
    fprintf(javafile, "\toutput%db.innerHTML = hextorgb(this.value).b;\n", i+1);
  	fprintf(javafile, "\tcolor%dout = { color: this.value };\n", i+1);
  	fprintf(javafile, "\t%s_material = define_material(color%dout,1);\n", element, i+1);
    for(j=0;j<ntemplates;j++) {
      fprintf(javafile, "\tif (atom%dstat===0 && togg%dstat===0) {\n", i+1, j+1);
  	  fprintf(javafile, "\t\tremove_objects(%s_atoms_template%d);\n", element, j+1);
  	  fprintf(javafile, "\t\t%s_atoms_template%d = drawmol_spheres(template%d,'%s',%s_material, size%d);\n", element, j+1, j+1, element, element, i+1);
  	  fprintf(javafile, "\t\tshow_objects(%s_atoms_template%d);\n\t}\n", element, j+1);
  	  fprintf(javafile, "\telse {\n\t\tremove_objects(%s_atoms_template%d);\n", element, j+1);
  	  fprintf(javafile, "\t\t%s_atoms_template%d = drawmol_spheres(template%d,'%s',%s_material, size%d);\n", element, j+1, j+1, element, element, i+1);
  	  fprintf(javafile, "\t\tdontshow_objects(%s_atoms_template%d);\n\t}\n", element, j+1);
      for(k=0;k<factorial;k++) {
        if(strstr(bonds[k], element) != NULL) {
          for(l=0;l<elementcount;l++) {
            ElementName(elementlist[l], element1);
            for(m=0;m<elementcount;m++) {
              ElementName(elementlist[m], element2);
              strcpy(element3, element1);
              strcat(element3, element2);
              if(strcmp(bonds[k], element3) == 0) {
                fprintf(javafile, "\tif (bond%dstat===0 && togg%dstat===0) {\n", k+1, j+1);
                fprintf(javafile, "\t\tremove_objects(%s_bonds_template%d);\n", bonds[k], j+1);
                fprintf(javafile, "\t\t%s_bonds_template%d = drawbonds_cylinder(template%d,'%s',%s_material,'%s',%s_material,0.1,len%d,rad%d);\n", bonds[k], j+1, j+1, element1, element1, element2, element2, i+1, i+1);
                fprintf(javafile, "\t\tshow_objects(%s_bonds_template%d);\n\t}\n", bonds[k], j+1);
                fprintf(javafile, "\telse {\n\t\tremove_objects(%s_bonds_template%d);\n", bonds[k], j+1);
                fprintf(javafile, "\t\t%s_bonds_template%d = drawbonds_cylinder(template%d,'%s',%s_material,'%s',%s_material,0.1,len%d,rad%d);\n", bonds[k], j+1, j+1, element1, element1, element2, element2, i+1, i+1);
                fprintf(javafile, "\t\tdontshow_objects(%s_bonds_template%d);\n\t}\n", bonds[k], j+1);
              }
            }
          }
        }
      }
    }
    fprintf(javafile, "}\n\n");
	}
  for(i=0;i<ntemplates;i++) {
    fprintf(javafile, "var CPs_template%d = drawCPs(geo,cell,template%d,CPcoeff,400.0);\n", i+1, i+1);
  }
  fprintf(javafile, "var slider = document.getElementById('myRange');\n");
  fprintf(javafile, "var value = document.getElementById('cpscale');\n");
  fprintf(javafile, "var scale = document.getElementById('demo');\n");
  fprintf(javafile, "scale.innerHTML = value.value;\n");
  fprintf(javafile, "CPscale = value.value;\n");
  fprintf(javafile, "slider.onclick = function() {\n");
  fprintf(javafile, "\tscale.innerHTML = value.value;\n");
  fprintf(javafile, "\tCPscale = value.value;\n");
  for(i=0;i<ntemplates;i++) {
    fprintf(javafile, "\tif (togg%dstat===0) {\n", i+1);
  	fprintf(javafile, "\t\tremove_objects(CPs_template%d);\n", i+1);
  	fprintf(javafile, "\t\tCPs_template%d = drawCPs(geo,cell,template%d,CPcoeff,CPscale);\n", i+1, i+1);
  	fprintf(javafile, "\t\tshow_objects(CPs_template%d);\n\t}\n", i+1);
  	fprintf(javafile, "\telse {\n\t\tremove_objects(CPs_template%d);\n", i+1);
  	fprintf(javafile, "\t\tCPs_template%d = drawCPs(geo,cell,template%d,CPcoeff,CPscale);\n", i+1, i+1);
  	fprintf(javafile, "\t\tdontshow_objects(CPs_template%d);\n\t}\n", i+1);
  }
  fprintf(javafile, "}\n\n");
  buttoncounter = 1;
  for(i=0;i<elementcount;i++) {
    ElementName(elementlist[i], element);
    fprintf(javafile, "const atom%d = document.getElementById('atom%d');\n", i+1, i+1);
    fprintf(javafile, "atom%d.addEventListener('click', buttonfunction%d);\n", i+1, buttoncounter);
    fprintf(javafile, "function buttonfunction%d() {\n", i+1);
    fprintf(javafile, "\tif (atom%dstat===0) {\n", i+1);
    for(j=0;j<ntemplates;j++) {
      fprintf(javafile, "\t\tdontshow_objects(%s_atoms_template%d);\n", element, j+1);
    }
    fprintf(javafile, "\t\tatom%dstat = 1;\n\t}\n", i+1);
    fprintf(javafile, "\telse {\n");
    for(j=0;j<ntemplates;j++) {
      fprintf(javafile, "\t\tshow_objects(%s_atoms_template%d);\n", element, j+1);
    }
    fprintf(javafile, "\t\tatom%dstat = 0;\n\t}\n}\n\n", i+1);
    buttoncounter++;
  }
  for(i=0;i<factorial;i++) {
    fprintf(javafile, "const bond%d = document.getElementById('bond%d');\n", i+1, i+1);
    fprintf(javafile, "bond%d.addEventListener('click', buttonfunction%d);\n", i+1, buttoncounter);
    fprintf(javafile, "function buttonfunction%d() {\n", buttoncounter);
    fprintf(javafile, "\tif (bond%dstat===0) {\n", i+1);
    for(j=0;j<ntemplates;j++) {
      fprintf(javafile, "\t\tdontshow_objects(%s_bonds_template%d);\n", bonds[i], j+1);
    }
    fprintf(javafile, "\t\tbond%dstat = 1;\n\t}\n", i+1);
    fprintf(javafile, "\telse {\n");
    for(j=0;j<ntemplates;j++) {
      fprintf(javafile, "\t\tshow_objects(%s_bonds_template%d);\n", bonds[i], j+1);
    }
    fprintf(javafile, "\t\tbond%dstat = 0;\n\t}\n}\n\n", i+1);
    buttoncounter++;
  }
  for(i=0;i<ntemplates;i++) {
    fprintf(javafile, "const ttemplate%d = document.getElementById('temp%d');\n", i+1, i+1);
    fprintf(javafile, "ttemplate%d.addEventListener('click', buttonfunction%d);\n", i+1, buttoncounter);
    fprintf(javafile, "function buttonfunction%d() {\n", buttoncounter);
    fprintf(javafile, "\tif (togg%dstat===0) {\n", i+1);
    for(j=0;j<elementcount;j++) {
      ElementName(elementlist[j], element);
      fprintf(javafile, "\t\tdontshow_objects(%s_atoms_template%d);\n", element, i+1);
    }
    for(j=0;j<factorial;j++) {
      fprintf(javafile, "\t\tdontshow_objects(%s_bonds_template%d);\n", bonds[j], i+1);
    }
    fprintf(javafile, "\t\tdontshow_objects(CPs_template%d);\n", i+1);
    fprintf(javafile, "\t\ttogg%dstat = 1;\n\t}\n", i+1);
    fprintf(javafile, "\telse {\n");
    for(j=0;j<elementcount;j++) {
      ElementName(elementlist[j], element);
      fprintf(javafile, "\t\tshow_objects(%s_atoms_template%d);\n", element, i+1);
    }
    for(j=0;j<factorial;j++) {
      fprintf(javafile, "\t\tshow_objects(%s_bonds_template%d);\n", bonds[j], i+1);
    }
    fprintf(javafile, "\t\tshow_objects(CPs_template%d);\n", i+1);
    fprintf(javafile, "\t\ttogg%dstat = 0;\n\t}\n}\n", i+1);
    buttoncounter++;
  }
  if (printvmap==1) OutputWeight(gridin, &vmap, gridout);
  strncpy(filename, cpoutname, STRMAX);
  strncat(filename, "-geo", STRMAX);
  fptr = fopen(filename, "w");


  for (i=0; i<gridin->nion; i++) {
    ElementName(gridin->zatomic[i], element);
    fprintf(fptr, "%s  %20.14f  %20.14f  %20.14f\n", element,gridin->xcart[i]*R_BOHR, gridin->ycart[i]*R_BOHR, gridin->zcart[i]*R_BOHR);
    }
  fclose(fptr);
  strncpy(filename, cpoutname, STRMAX);
  strncat(filename, "-cell", STRMAX);
  fptr = fopen(filename, "w");
  fprintf(fptr, "%20.14f  %20.14f  %20.14f\n", gridin->cella_x*R_BOHR, gridin->cella_y*R_BOHR, gridin->cella_z*R_BOHR);
  fprintf(fptr, "%20.14f  %20.14f  %20.14f\n", gridin->cellb_x*R_BOHR, gridin->cellb_y*R_BOHR, gridin->cellb_z*R_BOHR);
  fprintf(fptr, "%20.14f  %20.14f  %20.14f\n", gridin->cellc_x*R_BOHR, gridin->cellc_y*R_BOHR, gridin->cellc_z*R_BOHR);
  fclose(fptr);
  fclose(javafile);

  /* End Jonathan's interference */


  printf("  Projections finished\n");
  return 0;
}

void write_param(FILE * f2)
{
      int j;
      int count,stop;
      int check;
      char atomname [3];
      char str [STRMAX];
      FILE * f3;
      char element[STRMAX], profname[STRMAX];
      int i=0, k=0, stop2[NIONMAX];
      double elec=0.0, temp_zion[NIONMAX];
      FILE * fptr;
      fprintf(f2, "; CPpackage INPUT FILE FOR %s",abinitout);
      fprintf(f2, "\n\n\n");
      fprintf(f2, "ABINIT_OUT     %s \t;  Abinit out file \n",abinitout);
      f3=fopen(abinitout,"r");
      if(f3==NULL) {
          printf("  %s file not found.\n",abinitout);
          exit(1);
      }
      stop=0;
      while(stop==0) {
         check=fscanf(f3,"%s",str);
         if(check==EOF) {
            printf("  %s file incomplete.\n",abinitout);
            exit(1);
         }
         if(strcmp(str,"root")==0) {
             stop=1;
             FinishLine(f3);
             fscanf(f3,"%s %s %s %s %s %s %s",str,str,str,str,str,str,abinitname);
         }
      }
      fclose(f3);
      snprintf(logname, STRMAX, "%s-cplog", cpoutname);
      cplog = fopen(logname, "w");
      fprintf(f2, "ABINIT_ROOT    %s \t;  Abinit root for _o_  files \n",abinitname);
      fprintf(f2, "DEFAULT_MODE   %d \t;  Use defaults?  1=yes, otherwise no.  \n",1);
      fprintf(f2, "EWALD_MODE     %d \t;  Map E_Ewald+E_alpha?  1=yes, otherwise no. \n",1);
      fprintf(f2, "NONLOCAL_MODE  %d \t;  Map exact NL energy?  1=yes, otherwise no. \n",1);
      fprintf(f2, "\n");
      fprintf(f2, ";   Contact volume construction options (CV_MODE): \n");
      fprintf(f2, ";    -2 = Distance-based contact volumes. \n");
      fprintf(f2, ";    -1 = Hirshfeld-inpired atoms. \n");
      fprintf(f2, ";     0 = Hirshfeld-inspired contact volumes. \n");
      fprintf(f2, ";     1 = Closest intercept with Bader volume surface. \n");
      fprintf(f2, ";     2 = Bader volume + Closest neighbor.   \n");
      fprintf(f2, ";     3 = Bader volume + highest Hirshfeld weight. \n");
      fprintf(f2, ";     4 = Bader volume + gradient of _lLOC.xsf file from DFTMad. \n");
      fprintf(f2, ";     5 = Bader atoms (interactive). \n");
      fprintf(f2, ";     6 = Bader atoms (automatic).  \n");
      fprintf(f2, ";     7 = Bader volume + full Hirshfeld neighbors. \n");
      fprintf(f2, "CV_MODE        %d\n\n",CV_mode);
      fprintf(f2, "TEMPLATES      %d \t; number of template .xtl files provided as list below.\n\n",0);

      Den2XSF(abinitname, dseq, "DEN", &deneq, &deneq2);
      SymAtoms(&deneq, &smap);
      fprintf(f2, "ATOM_MODELS\n");
      for (i=0; i<deneq.nion; i++) stop2[i] = 0;
      for (i=0; i<deneq.nion; i++) {
         if (stop2[i]==1) continue;
         ElementName(deneq.zatomic[i], element);
         fprintf(f2,"  %s-0       \t;  ",element);
         if (smap.nequiv[i]>1) {
           fprintf(f2, "density profile file for atom #%d (%s, %d equivalent sites). \n",
             i+1, element, smap.nequiv[i]);
          } else {
           fprintf(f2, "density profile file for atom #%d (%s, %d site).\n",
             i+1, element, smap.nequiv[i]);
          }
          fprintf(f2,"  %d                \t;  0= valence only, 1=semicore.\n",0);
          fprintf(f2,"  %lf             \t;  number of semicore electrons.\n",0.000);
          for(j=i;j<deneq.nion;j++) {
                 if(smap.equiv[i][j]==1) stop2[j]=1;
          }
    }
    fprintf(f2,"END\n");
    fclose(cplog);
}

void read_param(FILE * f2)
{
      int i,j;
      int n_s;
      int n_p;
      int n_d;
      int stop;
      char keyword [200];
      int check;
      FILE * testfile;
      int tag1;
      int tag2;
      double Hii_s;
      double Hii_p;
      double Hii_d;
      double zeta_s;
      double zeta_p;
      double zeta1_d;
      double zeta2_d;
      double zeta1_d_wght;
      double zeta2_d_wght;
      int count;
      double test1;
      int test2;
      double weight;
      double dos_param;
      int stop2[NIONMAX];
      stop=0;
      snprintf(logname, STRMAX, "%s-cplog", cpoutname);
      cplog = fopen(logname, "w");
      snprintf(htmlname, STRMAX, "%s.html", cpoutname);
      htmlfile = fopen(htmlname,"w");
      snprintf(javaname, STRMAX, "%s.js", cpoutname);
      javafile = fopen(javaname,"w");
      while (stop==0) {
         check=fscanf(f2,"%s",keyword);
         if(check==EOF) stop=1;
         if (keyword[0]==59) FinishLine(f2);
         if (strcmp(keyword,"TEMPLATES")==0) {
             check=fscanf(f2,"%d",&ntemplates);
             FinishLine(f2);
             if(ntemplates>0) {
                for(j=0;j<ntemplates;j++) {
                 check=fscanf(f2,"%s",templatefiles[j]);
                 testfile=fopen(templatefiles[j],"r");
                 if(testfile==NULL) {
                      printf("%s not found.\n",templatefiles[j]);
                      exit(1);
                 }
                 fclose(testfile);
                }
             }
         }
         if (strcmp(keyword,"ABINIT_OUT")==0) {
            check=fscanf(f2,"%s", abinitout);
            FinishLine(f2);
         }
         if (strcmp(keyword,"END")==0) stop=1;
         if (strcmp(keyword,"ABINIT_ROOT")==0) {
            check=fscanf(f2,"%s", abinitname);
            printf("  ABINIT_ROOT = %s\n",abinitname);
            FinishLine(f2);
            Den2XSF(abinitname, dseq, "DEN", &deneq, &deneq2);
            SymAtoms(&deneq, &smap);
         }
         if (strcmp(keyword,"DEFAULT_MODE")==0) {
             check=fscanf(f2,"%d",&standard);
             if(standard!=1) standard=0;
             FinishLine(f2);
             printf("  DEFAULT_MODE = %d\n",standard);
         }
         if (strcmp(keyword,"EWALD_MODE")==0) {
             check=fscanf(f2,"%d",&E_Ewald_option);
             FinishLine(f2);
             printf("  EWALD_MODE = %d\n",E_Ewald_option);
         }
         if (strcmp(keyword,"NONLOCAL_MODE")==0) {
             check=fscanf(f2,"%d",&mapnonloc);
             FinishLine(f2);
             printf("  NONLOCAL = %d\n",mapnonloc);
         }
         if (strcmp(keyword,"CV_MODE")==0) {
             check=fscanf(f2,"%d",&CV_mode);
             FinishLine(f2);
             printf("  CV_MODE = %d\n",CV_mode);
         }

/*
char profile_filenames[NIONMAX][STRMAX];
int psptypes[NIONMAX];
double sccounts[NIONMAX];
*/
         if (strcmp(keyword,"ATOM_MODELS")==0) {
             count=0;
             for (i=0; i<deneq.nion; i++) stop2[i] = 0;
             for (i=0; i<deneq.nion; i++) {
               if (stop2[i]==1) continue;
               check=fscanf(f2,"%s",profile_filenames[count]);
               FinishLine(f2);
               check=fscanf(f2,"%d",&psptypes[count]);
               FinishLine(f2);
               check=fscanf(f2,"%lf",&sccounts[count]);
               FinishLine(f2);
               for(j=i;j<deneq.nion;j++) {
                 if(smap.equiv[i][j]==1) stop2[j]=1;
               }
               count++;
             }
         }
      }
}


/* MAIN FUNCTION */

int main(int argc, char * argv[]) {
  /* calls: AllocAll, AllocInt, AllocDbl, AllocHirsh, Getkm, AssignContact, AverageContact,
 *          CalcCP, CoordSearch, CoreUnwarp, ErrorCheck, SymAtoms, SetBubbles,
 *          ReadNonlocAtom, MapNonloc, ReadEwald, CalcEwald, CalcEalpha, MapEwald,
 *          GridStats, IdentifyXC, MapEntot, OutputXSF, PrintAverage, PrintResults,
 *          ReadAll, ReadProfile, SetOptions */
  char option[STRMAX], pspoption_in[STRMAX];
  int errchk=0,i,j,count,check;
  FILE * corefile;
  FILE * inifile;
  double sccounts[NIONMAX];
  double stop2[NIONMAX];
  double z_ion_temp;
  char elementname[STRMAX];
  char inifilename[STRMAX],delta_pot_file[STRMAX];
  printf("\nFREDRICKSON GROUP DFT-CHEMICAL PRESSURE ANALYSIS\n");
  if(argc<3) {
       printf("  Usage:  CPpackage2_BaderCV4 <abinit out file> <cp out file root>\n");
       exit(0);
  }
  strcpy(abinitout,argv[1]);
  strcpy(cpoutname,argv[2]);
  sprintf(inifilename,"%s.ini",cpoutname);
  inifile=fopen(inifilename,"r");
  if(inifile==NULL) {
        printf("  This looks like a new CP calculation.  Preparing %s parameter file...\n",inifilename);
        inifile=fopen(inifilename,"w");
        write_param(inifile);
        fclose(inifile);
        printf("  Done!\n\n");
        printf("  Please edit %s and rerun CPpackage.\n",inifilename);
        exit(0);
  }
  read_param(inifile);
  fclose(inifile);
  fprintf(cplog, "FREDRICKSON GROUP DFT-CHEMICAL PRESSURE ANALYSIS\n");
  fprintf(cplog, "Last modified: XX July 2018\n");
  fprintf(cplog, "Echo of input: %s %s %s %d\n", abinitout, abinitname, cpoutname, standard);
  fprintf(cplog, "All values are in a.u. unless otherwise noted\n\n");
  if (standard==1) {
    fprintf(cplog, "Selected default options:\n");
    fprintf(cplog, "Map kinetic, local, hartree, exchange-correlation energy, and nonlocal pseudopotential terms\n");
    fprintf(cplog, "Map localized Ewald and E_alpha terms if semicore electrons present\n");
    fprintf(cplog, "With hirshfeld-inspired core unwarping and symmetry restoration\n");
    fprintf(cplog, "And hirshfeld-inspired integration: l_max=%d and tolerance=%.2f%%\n\n", lmax, tolerance*100);
  } else {
    errchk = SetOptions();  /* user input of custom settings */
    if (errchk!=0) goto end;
  }

  if (mapcore==1 || scheme==1) printf("  Reading files...\n");
  else printf("  Reading files\n");
  errchk = ReadAll();  /* assigns volhi, voleq, vollo */
  printf("HERE!\n");
  if (errchk!=0) goto end;
  if (mapcore==1 || scheme==1 || isradii==1) SymAtoms(&deneq, &smap);
  if (mapcore==1 || scheme==1) {
    errchk = ReadProfile(&deneq);  /* radial electron density profiles */
    if (errchk!=0) goto end;
  }
//  if (mapnonloc==2) {
//    errchk = ReadNonlocAtom(&deneq);
//  }
  if (isradii==1) SetBubbles(&deneq);  /* bubble radii for distance-based contact volumes */
  if (mapxc==1) {
    errchk = IdentifyXC(ixc);  /* determine which XC functional was used in Abinit */
    if (errchk!=0) goto end;
  }
  errchk = ErrorCheck(volhi, voleq, vollo);  /* check that volumes and ions make sense */
  printf("\n  Creating energy density maps\n");

  /* allocate memory dynamically for voxel grids */
  errchk = AllocDbl(&core);
  errchk = AllocDbl(&bader);
  errchk = AllocDbl(&bader_temp);
  errchk = AllocDbl(&cp);
  errchk = AllocDbl(&cp_Y);
  errchk = AllocDbl(&etothi);
  etothi.volvox = volhi;
  errchk = AllocDbl(&etoteq);
  etoteq.volvox = voleq;
  errchk = AllocDbl(&etotlo);
  etotlo.volvox = vollo;
  errchk = AllocDbl(&temp);
  errchk = AllocDbl(&vxc);
  errchk = AllocDbl(&lochi);
  errchk = AllocDbl(&loceq);
  errchk = AllocDbl(&loclo);
  errchk = AllocDbl(&etotlo_temp);
  errchk = AllocDbl(&etothi_temp);

  if (mapkin==2) {
    errchk = AllocDbl(&kdenhi);
    errchk = AllocDbl(&kdeneq);
    errchk = AllocDbl(&kdenlo);
  }
  if (E_Ewald_option !=0) {
    errchk = AllocDbl(&ewaldhi);
    errchk = AllocDbl(&ewaldeq);
    errchk = AllocDbl(&ewaldlo);
    errchk = AllocDbl(&alphahi);
    errchk = AllocDbl(&alphalo);
  }
  if (mapnonloc!=0) {
    errchk = AllocDbl(&nonlochi);
    errchk = AllocDbl(&nonloclo);
    errchk = AllocDbl(&nonlochi2);
    errchk = AllocDbl(&nonloclo2);
  }
  if (nspin==2) {
//    errchk = AllocDbl(&vhahi2);
//    errchk = AllocDbl(&vhalo2);
//    errchk = AllocDbl(&pothi2);
//    errchk = AllocDbl(&potlo2);
    errchk = AllocDbl(&core2);
    errchk = AllocDbl(&temp2);
    errchk = AllocDbl(&vxc2);
  }
  printf("Here2!\n");
  errchk = MapEntot(&denhi, &denhi2, &kdenhi, &kdenhi2, &ldenhi, &ldenhi2, &pothi, &pothi2,
    &vhxchi, &vhxchi2, &vhahi, &vhahi2, dshi, volhi, &etothi, &lochi);
  if (errchk!=0) goto end;
  errchk = MapEntot(&deneq, &deneq2, &kdeneq, &kdeneq2, &ldeneq, &ldeneq2, &poteq, &poteq2,
    &vhxceq, &vhxceq2, &vhaeq, &vhaeq2, dseq, voleq, &etoteq, &loceq);
  if (errchk!=0) goto end;
  errchk = MapEntot(&denlo, &denlo2, &kdenlo, &kdenlo2, &ldenlo, &ldenlo2, &potlo, &potlo2,
    &vhxclo, &vhxclo2, &vhalo, &vhalo2, dslo, vollo, &etotlo, &loclo);
  if (errchk!=0) goto end;
  printf("here3!\n");
  errchk = Getkm(&deneq);  /* determines appropriate supercell */
  if (errchk!=0) goto end;
  if (mapcore==1) {
    errchk = AllocHirsh(&hmap);
    printf("  Performing core unwarping\n");
    errchk = CoreUnwarp(&denhi, &deneq, &denlo, &etothi, &etotlo);
  }
  errchk = AllocInt(&vmap);  /* allocate memory dynamically for integration grids */
  if(CV_mode == 7) {
        Den2XSF(abinitname, dseq, "DEN", &delta_pot_map, &delta_pot_map2);
        errchk = CoordSearchBaderHirsh(&deneq, &vmap,&smap);
        if (mapnonloc==1) {
          printf("  Mapping nonlocal energy\n");
          errchk = MapNonloc(&deneq, &nonlochi, &nonlochi2, &etothi, 1);
          errchk = MapNonloc(&deneq, &nonloclo, &nonloclo2, &etotlo, 3);
        }
        errchk = ReadEwald(3);
        if (E_Ewald_option != 0) {
            printf("  Mapping Ewald energy\n");
            errchk = CalcCP_prelim(&denhi, &deneq, &denlo, &etothi, &etotlo, &cp);
            errchk = CalcEwald(&deneq);
            errchk = CalcEalpha(&deneq);
            errchk = MapEwald_Bader(&denhi, &cp, &etothi, &ewaldhi, &alphahi, &nonlochi, 1, E_ewald[0], E_ewald[1], &densceq,&vmap);
            errchk = MapEwald_Bader(&denlo, &cp, &etotlo, &ewaldlo, &alphalo, &nonloclo, 3, E_ewald[2], E_ewald[1], &densceq,&vmap);
        }
        sprintf(delta_pot_file,"%s-HirshMap.xsf",cpoutname);
        corefile = fopen(delta_pot_file,"w");
        OutputXSF(corefile,&delta_pot_map,&delta_pot_map);
        fclose(corefile);
        printf("\n  Creating chemical pressure map\n");
        errchk = CalcCP(&denhi, &deneq, &denlo, &etothi, &etotlo, &cp);
        if (errchk!=0) goto end;
        GridStats(&cp, "CP XSFfile stats");  /* calculates statistical values of the voxel grid */
        fclose(cplog);  /* close to ensure writing of CPmap log */
        cplog = fopen(logname, "a");  /* reopen for CPintegrate */
        if (errchk!=0) goto end;
        printf("  Assigning %d x %d x %d = %d voxels to contact volumes: full Hirshfeld neighbors\n", ngx, ngy, ngz, ngx*ngy*ngz);
        AssignContactHirsh(&vmap, &cp);
        printf("  Averaging within contact volumes\n");
        AverageContactHirsh(&vmap, &cp, &cp_Y);
  }
  // Jonathan adding code

  // ADD CV_MODE 6: THE BADER ATOMS AUTOCALIBRATION METHOD

  if (CV_mode == 6) {
        errchk = CoordSearchBader(&deneq, &vmap,&smap);
        if (mapnonloc==1) {
          printf("  Mapping nonlocal energy\n");
          errchk = MapNonloc(&deneq, &nonlochi, &nonlochi2, &etothi, 1);
          errchk = MapNonloc(&deneq, &nonloclo, &nonloclo2, &etotlo, 3);
        }
        errchk = ReadEwald(3);
        ScaleGrid(&etothi,1.0,&etothi_temp);
        ScaleGrid(&etotlo,1.0,&etotlo_temp);


        repeatcp1:
        if (E_Ewald_option != 0) {
            printf("  Mapping Ewald energy\n");
            errchk = CalcCP_prelim(&denhi, &deneq, &denlo, &etothi, &etotlo, &cp);
            errchk = CalcEwald(&deneq);
            errchk = CalcEalpha(&deneq);
            errchk = MapEwald_Bader(&denhi, &cp, &etothi, &ewaldhi, &alphahi, &nonlochi, 1, E_ewald[0], E_ewald[1], &densceq,&vmap);
            errchk = MapEwald_Bader(&denlo, &cp, &etotlo, &ewaldlo, &alphalo, &nonloclo, 3, E_ewald[2], E_ewald[1], &densceq,&vmap);
        }
        printf("\n  Creating chemical pressure map\n");
        errchk = CalcCP(&denhi, &deneq, &denlo, &etothi, &etotlo, &cp);
        if (errchk!=0) goto end;
//        GridStats(&cp, "CP XSFfile stats");  /* calculates statistical values of the voxel grid */
        fclose(cplog);  /* close to ensure writing of CPmap log */
        cplog = fopen(logname, "a");  /* reopen for CPintegrate */
        if (errchk!=0) goto end;
        printf("  Averaging within Bader volumes\n");
        AverageAtomFast(&vmap, &cp, &cp_Y);
        PrintAverage3(&cp, &cp_Y);  /* prints integrated CP to the screen and log file */
        ScaleGrid(&etothi_temp,1.0,&etothi);
        ScaleGrid(&etotlo_temp,1.0,&etotlo);
        count=0;
        if (autocali_stopper == 0) {
          autocali_start();
          autocali_stopper = 1;
        }
        autocali_iteration();
        for (i=0; i<deneq.nion; i++) stop2[i] = 0;
        for (i=0; i<deneq.nion; i++) {
               if(stop2[i]==1) continue;
               ElementName(deneq.zatomic[i], elementname);
               // DONT ASK FOR INPUT, DECIDE YOUR OWN DESTINY!!!
               // printf("Enter new localized electron count for atom %d(%s) [%lf]:",i+1,elementname,sc_elec[i]);
               // check=scanf("%lf",&sccounts[count]); //previously used sccounts[count]
               sccounts[i] = fabs(es2_temp[i]);
               printf("New Elec Val for site %d is %lf, es2_temp of %d is %lf\n", i, sccounts[i], i, es2_temp[i]);
               z_ion_temp = vo_elec[i]+sc_elec[i];
               vo_elec[i] = z_ion_temp - sccounts[i];
               sc_elec[i] = sccounts[i];
               for(j=i;j<deneq.nion;j++) {
                 if(smap.equiv[i][j]==1) {
                    stop2[j]=1;
                    vo_elec[j] = vo_elec[i];
                    sc_elec[j] = sc_elec[i];
                 }
               }
               count++;
        }
        memcpy(ps1, ps2, sizeof(ps1));
        memcpy(es1, es2, sizeof(es1));
        memcpy(es2, es2_temp, sizeof(es2));
        goto repeatcp1;
  }


  //Jonathan stop adding code


  if(CV_mode == 5) {
        errchk = CoordSearchBader(&deneq, &vmap,&smap);
        if (mapnonloc==1) {
          printf("  Mapping nonlocal energy\n");
          errchk = MapNonloc(&deneq, &nonlochi, &nonlochi2, &etothi, 1);
          errchk = MapNonloc(&deneq, &nonloclo, &nonloclo2, &etotlo, 3);
        }
        errchk = ReadEwald(3);
        ScaleGrid(&etothi,1.0,&etothi_temp);
        ScaleGrid(&etotlo,1.0,&etotlo_temp);


        repeatcp2:
        if (E_Ewald_option != 0) {
            printf("  Mapping Ewald energy\n");
            errchk = CalcCP_prelim(&denhi, &deneq, &denlo, &etothi, &etotlo, &cp);
            errchk = CalcEwald(&deneq);
            errchk = CalcEalpha(&deneq);
            errchk = MapEwald_Bader(&denhi, &cp, &etothi, &ewaldhi, &alphahi, &nonlochi, 1, E_ewald[0], E_ewald[1], &densceq,&vmap);
            errchk = MapEwald_Bader(&denlo, &cp, &etotlo, &ewaldlo, &alphalo, &nonloclo, 3, E_ewald[2], E_ewald[1], &densceq,&vmap);
        }
        printf("\n  Creating chemical pressure map\n");
        errchk = CalcCP(&denhi, &deneq, &denlo, &etothi, &etotlo, &cp);
        if (errchk!=0) goto end;
//        GridStats(&cp, "CP XSFfile stats");  /* calculates statistical values of the voxel grid */
        fclose(cplog);  /* close to ensure writing of CPmap log */
        cplog = fopen(logname, "a");  /* reopen for CPintegrate */
        if (errchk!=0) goto end;
        printf("  Averaging within Bader volumes\n");
        AverageAtomFast(&vmap, &cp, &cp_Y);
        PrintAverage2(&cp, &cp_Y);  /* prints integrated CP to the screen and log file */
        ScaleGrid(&etothi_temp,1.0,&etothi);
        ScaleGrid(&etotlo_temp,1.0,&etotlo);
        count=0;
        for (i=0; i<deneq.nion; i++) stop2[i] = 0;
        for (i=0; i<deneq.nion; i++) {
               if(stop2[i]==1) continue;
               ElementName(deneq.zatomic[i], elementname);
               printf("Enter new localized electron count for atom %d(%s) [%lf]:",i+1,elementname,sc_elec[i]);
               check=scanf("%lf",&sccounts[count]);
               z_ion_temp = vo_elec[i]+sc_elec[i];
               vo_elec[i] = z_ion_temp - sccounts[count];
               sc_elec[i] = sccounts[count];
               for(j=i;j<deneq.nion;j++) {
                 if(smap.equiv[i][j]==1) {
                    stop2[j]=1;
                    vo_elec[j] = vo_elec[i];
                    sc_elec[j] = sc_elec[i];
                 }
               }
               count++;
        }
        goto repeatcp2;
  }
  if(CV_mode == 4) {
        sprintf(delta_pot_file,"%s_lLOC.xsf",abinitout);
        ReadXSF(delta_pot_file,&delta_pot_map);
        errchk = CoordSearchBader(&deneq, &vmap,&smap);
        if (mapnonloc==1) {
          printf("  Mapping nonlocal energy\n");
          errchk = MapNonloc(&deneq, &nonlochi, &nonlochi2, &etothi, 1);
          errchk = MapNonloc(&deneq, &nonloclo, &nonloclo2, &etotlo, 3);
        }
        errchk = ReadEwald(3);
        if (E_Ewald_option != 0) {
            printf("  Mapping Ewald energy\n");
            errchk = CalcEwald(&deneq);
            errchk = CalcEalpha(&deneq);
            errchk = MapEwald_Bader(&denhi, &cp, &etothi, &ewaldhi, &alphahi, &nonlochi, 1, E_ewald[0], E_ewald[1], &densceq,&vmap);
            errchk = MapEwald_Bader(&denlo, &cp, &etotlo, &ewaldlo, &alphalo, &nonloclo, 3, E_ewald[2], E_ewald[1], &densceq,&vmap);

        }
        printf("\n  Creating chemical pressure map\n");
        errchk = CalcCP(&denhi, &deneq, &denlo, &etothi, &etotlo, &cp);
        if (errchk!=0) goto end;
        GridStats(&cp, "CP XSFfile stats");  /* calculates statistical values of the voxel grid */
        fclose(cplog);  /* close to ensure writing of CPmap log */
        cplog = fopen(logname, "a");  /* reopen for CPintegrate */
        if (errchk!=0) goto end;
        printf("  Assigning %d x %d x %d = %d voxels to contact volumes:  lLOC-based\n", ngx, ngy, ngz, ngx*ngy*ngz);
        AssignContact(&vmap, &cp);
        printf("  Averaging within contact volumes\n");
        AverageContact(&vmap, &cp, &cp_Y);
  }
  if(CV_mode == 3) {
        Den2XSF(abinitname, dseq, "DEN", &delta_pot_map, &delta_pot_map2);
        errchk = CoordSearchBader(&deneq, &vmap,&smap);
        if (mapnonloc==1) {
          printf("  Mapping nonlocal energy\n");
          errchk = MapNonloc(&deneq, &nonlochi, &nonlochi2, &etothi, 1);
          errchk = MapNonloc(&deneq, &nonloclo, &nonloclo2, &etotlo, 3);
        }
        if (mapnonloc==2) {
          printf("  Mapping nonlocal energy by valence electron density within atomic cells\n");
          errchk = MapNonloc2(&deneq, &nonlochi, &nonlochi2, &etothi, 1);
          errchk = MapNonloc2(&deneq, &nonloclo, &nonloclo2, &etotlo, 3);
        }
        errchk = ReadEwald(3);
        if (E_Ewald_option != 0) {
            printf("  Mapping Ewald energy\n");
            errchk = CalcCP_prelim(&denhi, &deneq, &denlo, &etothi, &etotlo, &cp);
            errchk = CalcEwald(&deneq);
            errchk = CalcEalpha(&deneq);
            errchk = MapEwald_Bader(&denhi, &cp, &etothi, &ewaldhi, &alphahi, &nonlochi, 1, E_ewald[0], E_ewald[1], &densceq,&vmap);
            errchk = MapEwald_Bader(&denlo, &cp, &etotlo, &ewaldlo, &alphalo, &nonloclo, 3, E_ewald[2], E_ewald[1], &densceq,&vmap);
        }
        sprintf(delta_pot_file,"%s-HirshMap.xsf",cpoutname);
        corefile = fopen(delta_pot_file,"w");
        OutputXSF(corefile,&delta_pot_map,&delta_pot_map);
        fclose(corefile);
        printf("\n  Creating chemical pressure map\n");
        errchk = CalcCP(&denhi, &deneq, &denlo, &etothi, &etotlo, &cp);
        if (errchk!=0) goto end;
        GridStats(&cp, "CP XSFfile stats");  /* calculates statistical values of the voxel grid */
        fclose(cplog);  /* close to ensure writing of CPmap log */
        cplog = fopen(logname, "a");  /* reopen for CPintegrate */
        if (errchk!=0) goto end;
        printf("  Assigning %d x %d x %d = %d voxels to contact volumes: Hirshfeld-inspired 2nd neighbors\n", ngx, ngy, ngz, ngx*ngy*ngz);
        AssignContact(&vmap, &cp);
        printf("  Averaging within contact volumes\n");
        AverageContact(&vmap, &cp, &cp_Y);
  }
  if(CV_mode == 2) {
        errchk = CoordSearchBader(&deneq, &vmap,&smap);
        if (mapnonloc==1) {
          printf("  Mapping nonlocal energy\n");
          errchk = MapNonloc(&deneq, &nonlochi, &nonlochi2, &etothi, 1);
          errchk = MapNonloc(&deneq, &nonloclo, &nonloclo2, &etotlo, 3);
        }
        if (mapnonloc==2) {
          printf("  Mapping nonlocal energy by valence electron density within atomic cells\n");
          errchk = MapNonloc2(&deneq, &nonlochi, &nonlochi2, &etothi, 1);
          errchk = MapNonloc2(&deneq, &nonloclo, &nonloclo2, &etotlo, 3);
        }
        errchk = ReadEwald(3);
        if (E_Ewald_option != 0) {
            printf("  Mapping Ewald energy\n");
            errchk = CalcCP_prelim(&denhi, &deneq, &denlo, &etothi, &etotlo, &cp);
            errchk = CalcEwald(&deneq);
            errchk = CalcEalpha(&deneq);
            errchk = MapEwald_Bader(&denhi, &cp, &etothi, &ewaldhi, &alphahi, &nonlochi, 1, E_ewald[0], E_ewald[1], &densceq,&vmap);
            errchk = MapEwald_Bader(&denlo, &cp, &etotlo, &ewaldlo, &alphalo, &nonloclo, 3, E_ewald[2], E_ewald[1], &densceq,&vmap);

        }
        printf("\n  Creating chemical pressure map\n");
        errchk = CalcCP(&denhi, &deneq, &denlo, &etothi, &etotlo, &cp);
        if (errchk!=0) goto end;
        GridStats(&cp, "CP XSFfile stats");  /* calculates statistical values of the voxel grid */
        fclose(cplog);  /* close to ensure writing of CPmap log */
        cplog = fopen(logname, "a");  /* reopen for CPintegrate */
        if (errchk!=0) goto end;
        printf("  Assigning %d x %d x %d = %d voxels to contact volumes:  distance-based 2nd neighbors\n", ngx, ngy, ngz, ngx*ngy*ngz);
        AssignContact(&vmap, &cp);
        printf("  Averaging within contact volumes\n");
        AverageContact(&vmap, &cp, &cp_Y);
  }
  if(CV_mode == 1) {
        errchk = CoordSearchBader(&deneq, &vmap,&smap);
        if (mapnonloc==1) {
          printf("  Mapping nonlocal energy\n");
          errchk = MapNonloc(&deneq, &nonlochi, &nonlochi2, &etothi, 1);
          errchk = MapNonloc(&deneq, &nonloclo, &nonloclo2, &etotlo, 3);
        }
        if (mapnonloc==2) {
          printf("  Mapping nonlocal energy by valence electron density within atomic cells\n");
          errchk = MapNonloc2(&deneq, &nonlochi, &nonlochi2, &etothi, 1);
          errchk = MapNonloc2(&deneq, &nonloclo, &nonloclo2, &etotlo, 3);
        }
        errchk = ReadEwald(3);
        if (E_Ewald_option != 0) {
            printf("  Mapping Ewald energy\n");
            errchk = CalcCP_prelim(&denhi, &deneq, &denlo, &etothi, &etotlo, &cp);
            errchk = CalcEwald(&deneq);
            errchk = CalcEalpha(&deneq);
            errchk = MapEwald_Bader(&denhi, &cp, &etothi, &ewaldhi, &alphahi, &nonlochi, 1, E_ewald[0], E_ewald[1], &densceq,&vmap);
            errchk = MapEwald_Bader(&denlo, &cp, &etotlo, &ewaldlo, &alphalo, &nonloclo, 3, E_ewald[2], E_ewald[1], &densceq,&vmap);
        }
        printf("\n  Creating chemical pressure map\n");
        errchk = CalcCP(&denhi, &deneq, &denlo, &etothi, &etotlo, &cp);
        if (errchk!=0) goto end;
        GridStats(&cp, "CP XSFfile stats");  /* calculates statistical values of the voxel grid */
        fclose(cplog);  /* close to ensure writing of CPmap log */
        cplog = fopen(logname, "a");  /* reopen for CPintegrate */
        if (errchk!=0) goto end;
        printf("  Assigning %d x %d x %d = %d voxels to contact volumes: Bader surface-defined 2nd neighbors\n", ngx, ngy, ngz, ngx*ngy*ngz);
        AssignContact(&vmap, &cp);
        printf("  Averaging within voxels\n");
        AverageContact(&vmap, &cp, &cp_Y);
  }
  if(CV_mode == 0) {
        errchk = CoordSearch(&deneq, &vmap);
        if (mapnonloc==1) {
          printf("  Mapping nonlocal energy\n");
          errchk = MapNonloc(&deneq, &nonlochi, &nonlochi2, &etothi, 1);
          errchk = MapNonloc(&deneq, &nonloclo, &nonloclo2, &etotlo, 3);
        }
       if (mapnonloc==2) {
          printf("  Mapping nonlocal energy by valence electron density within atomic cells\n");
          errchk = MapNonloc2(&deneq, &nonlochi, &nonlochi2, &etothi, 1);
          errchk = MapNonloc2(&deneq, &nonloclo, &nonloclo2, &etotlo, 3);
        }

        errchk = ReadEwald(3);
        if (E_Ewald_option != 0) {
            printf("  Mapping Ewald energy\n");
            errchk = CalcCP_prelim(&denhi, &deneq, &denlo, &etothi, &etotlo, &cp);
            errchk = CalcEwald(&deneq);
            errchk = CalcEalpha(&deneq);
            errchk = MapEwald(&denhi, &cp, &etothi, &ewaldhi, &alphahi, &nonlochi, 1, E_ewald[0], E_ewald[1], &densceq);
            errchk = MapEwald(&denlo, &cp, &etotlo, &ewaldlo, &alphalo, &nonloclo, 3, E_ewald[2], E_ewald[1], &densceq);
        }
        printf("\n  Creating chemical pressure map\n");
        errchk = CalcCP(&denhi, &deneq, &denlo, &etothi, &etotlo, &cp);
        if (errchk!=0) goto end;
        GridStats(&cp, "CP XSFfile stats");  /* calculates statistical values of the voxel grid */
        fclose(cplog);  /* close to ensure writing of CPmap log */
        cplog = fopen(logname, "a");  /* reopen for CPintegrate */
        if (errchk!=0) goto end;
        printf("  Assigning %d x %d x %d = %d voxels to contact volumes:  Hirshfeld-inspired\n", ngx, ngy, ngz, ngx*ngy*ngz);
        AssignContact(&vmap, &cp);
        printf("  Averaging within contact volumes\n");
        AverageContact(&vmap, &cp, &cp_Y);
  }
  if(CV_mode == -1) {
        errchk = CoordSearchAtom(&deneq, &vmap);
        if (mapnonloc==1) {
          printf("  Mapping nonlocal energy\n");
          errchk = MapNonloc(&deneq, &nonlochi, &nonlochi2, &etothi, 1);
          errchk = MapNonloc(&deneq, &nonloclo, &nonloclo2, &etotlo, 3);
        }
        errchk = ReadEwald(3);
        if (E_Ewald_option != 0) {
            printf("  Mapping Ewald energy\n");
            errchk = CalcCP_prelim(&denhi, &deneq, &denlo, &etothi, &etotlo, &cp);
            errchk = CalcEwald(&deneq);
            errchk = CalcEalpha(&deneq);
            errchk = MapEwald(&denhi, &cp, &etothi, &ewaldhi, &alphahi, &nonlochi, 1, E_ewald[0], E_ewald[1], &densceq);
            errchk = MapEwald(&denlo, &cp, &etotlo, &ewaldlo, &alphalo, &nonloclo, 3, E_ewald[2], E_ewald[1], &densceq);
        }
        printf("\n  Creating chemical pressure map\n");
        errchk = CalcCP(&denhi, &deneq, &denlo, &etothi, &etotlo, &cp);
        if (errchk!=0) goto end;
        GridStats(&cp, "CP XSFfile stats");  /* calculates statistical values of the voxel grid */
        fclose(cplog);  /* close to ensure writing of CPmap log */
        cplog = fopen(logname, "a");  /* reopen for CPintegrate */
        if (errchk!=0) goto end;
        printf("  Averaging within Hirshfeld-inspired atomic volumes\n");
        AverageAtom(&vmap, &cp, &cp_Y);
  }
  if(CV_mode == -2) {
        errchk = CoordSearchDist(&deneq, &vmap);
        if (mapnonloc==1) {
          printf("  Mapping nonlocal energy\n");
          errchk = MapNonloc(&deneq, &nonlochi, &nonlochi2, &etothi, 1);
          errchk = MapNonloc(&deneq, &nonloclo, &nonloclo2, &etotlo, 3);
        }
        if (mapnonloc==2) {
          printf("  Mapping nonlocal energy by valence electron density within atomic cells\n");
          errchk = MapNonloc2(&deneq, &nonlochi, &nonlochi2, &etothi, 1);
          errchk = MapNonloc2(&deneq, &nonloclo, &nonloclo2, &etotlo, 3);
        }

        errchk = ReadEwald(3);
        if (E_Ewald_option != 0) {
            printf("  Mapping Ewald energy\n");
            errchk = CalcCP_prelim(&denhi, &deneq, &denlo, &etothi, &etotlo, &cp);
            errchk = CalcEwald(&deneq);
            errchk = CalcEalpha(&deneq);
            errchk = MapEwald(&denhi, &cp, &etothi, &ewaldhi, &alphahi, &nonlochi, 1, E_ewald[0], E_ewald[1], &densceq);
            errchk = MapEwald(&denlo, &cp, &etotlo, &ewaldlo, &alphalo, &nonloclo, 3, E_ewald[2], E_ewald[1], &densceq);

        }
        printf("\n  Creating chemical pressure map\n");
        errchk = CalcCP(&denhi, &deneq, &denlo, &etothi, &etotlo, &cp);
        if (errchk!=0) goto end;
        GridStats(&cp, "CP XSFfile stats");  /* calculates statistical values of the voxel grid */
        fclose(cplog);  /* close to ensure writing of CPmap log */
        cplog = fopen(logname, "a");  /* reopen for CPintegrate */
        if (errchk!=0) goto end;
        printf("  Assigning %d x %d x %d = %d voxels to contact volumes: distance-based contact volumes\n", ngx, ngy, ngz, ngx*ngy*ngz);
        AssignContact(&vmap, &cp);
        printf("  Averaging within contact volumes\n");
        AverageContact(&vmap, &cp, &cp_Y);
  }
  if (errchk!=0) goto end;
  PrintAverage(&cp, &cp_Y);  /* prints integrated CP to the screen and log file */
  PrintResults(&cp, &cp_Y);  /* makes all necessary files for Figuretool */
  GridStats(&cp_Y, "Averaged CP XSFfile stats");
  end:
  if (errcount==0) {
    printf("  Logfile %s is ready to view. Goodbye!\n\n", logname);
    fprintf(cplog, "Normal exit without incident\n");
  } else {
    printf("  Logfile %s is ready to view\n", logname);
    if (errcount==1) printf("  There is %d warning\n\n", errcount);
    else printf("  There are %d warnings\n\n", errcount);
    fprintf(cplog, "\n***********************************************************\n");
    fprintf(cplog, "Random errors can occur when nproc is not a factor of nkpt\n");
    fprintf(cplog, "As a last resort re-run Abinit with 1 processor\n");
    fprintf(cplog, "***********************************************************\n\n");
    fprintf(cplog, "Last errchk value was %d\n", errchk);
    if (errcount==1) fprintf(cplog, "Exited with %d warning\n", errcount);
    else fprintf(cplog, "Exited with %d warnings\n", errcount);
  }
  fclose(cplog);
  return 0;  /* memory deallocation is handled automatically by the OS */
}
