/* Last updated 04/20/2021 by Erdong */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <ctype.h>

#define NTYPES_MAX 10
#define NATOM_MAX 1000
#define SITE_MAX 100
#define PI 3.1415926

int ntype;
int scvo_flag=1; //vo=1; sc=2
int spgorig_read=0; //0=not read; 1=read
int convert_primitive_flag=1; //1=Convert; 2=Do not convert
int ecut_unknown[NTYPES_MAX]; //0=known; 1=unknown
int rprim_found=0; //0=not found; 1=found
int kmesh_flag=1; //1=ngkpt; 2=kptrlatt
int CIF_source_flag=1; //1=ICSD; 2=JANA No longer used
int convert_FAIL=0; //0=success or no conversion; 1=fail
int spgaxor_relevant=0; //0=spgaxor irrelevant; 1 = spgaxor relevant; 2=trigoanl space group
int SOC_flag=1; //1=No spin orbit coupling; 2=Spin orbit coupling; 3=Spin polarization

struct CIF{
  double cell_a, cell_b, cell_c;
  char cell_a_str[20], cell_b_str[20], cell_c_str[20];
  double ang_alpha, ang_beta, ang_gamma;
  char ang_alpha_str[20], ang_beta_str[20], ang_gamma_str[20];
  char chemical_formula[NTYPES_MAX][20];
  char HM_FullspgName[10][50];
  char HM_FullspgName_oneword[50];
  char HM_STDName[50];
  char element_identity[NTYPES_MAX][3];
  int element_No[NTYPES_MAX];
  char site_identity[SITE_MAX][6];
  char site_identity_processed[SITE_MAX][6];
  int site_element_No[SITE_MAX];
  int site_eqv[SITE_MAX];
  int total_site;
  int Z_unit;
  int brvltt, nband, old_brvltt, old_chkprim;
  char spgroup_center[3];
  char spg_center_letter;
  int spgroup, spgorig, spg_type_No;
  int spgaxor;
  char xred_str[SITE_MAX][20], yred_str[SITE_MAX][20], zred_str[SITE_MAX][20];
  char xred_upd[SITE_MAX][20], yred_upd[SITE_MAX][20], zred_upd[SITE_MAX][20];
  double xred[SITE_MAX], yred[SITE_MAX], zred[SITE_MAX];
  int natom, natrd;
  int nelectrons_psp[NTYPES_MAX];
  char psp_name[NTYPES_MAX][50];
  int natom_type[NTYPES_MAX];
  int No_site_by_type[NTYPES_MAX];
  int site_typat[SITE_MAX];
  int ecut[NTYPES_MAX];
  int ecut_final;
  double ecut_lambda_bohr, ecut_lambda_A;
  double occupancy[SITE_MAX];
  char occupancy_str[SITE_MAX][20];
  char rprim_str[3][3][100];
  char kptrlatt_str[9][10];
  int kptrlatt[9];
  double shiftk[3];
  int nkpt;
  int kpt_mesh_accepted;
  double cell_volume;
  int ngkpt[5][3]; //ngkpt1: 0.3 A-1; ngkpt2: 0.2 A-1; ngkpt3: 0.1 A-1; ngkpt: 0.07 A-1; ngkpt: 0.05 A-1.
  int ngkpt_from_in[3];
} inputCIF;

void DOS_calculation_adjustment(struct CIF * inputCIF, int fining_factor){
  int i;
  printf("\nSince you are running DOS or pDOS calculation, your k-point mesh is increased by a factor of %d. \n",fining_factor);
  if(kmesh_flag==1){
    for(i=0;i<3;i++){
      inputCIF->ngkpt_from_in[i]=inputCIF->ngkpt_from_in[i]*fining_factor;
    }
  }else if(kmesh_flag==2){
    for(i=0;i<9;i++){
      inputCIF->kptrlatt[i]=inputCIF->kptrlatt[i]*fining_factor;
    }
  }
}

void TriHex_Adjustment(struct CIF * inputCIF, double tolerance){
  int i;
  double diff;
  double one_third=0.333333333333333;
  double two_third=0.666666666666667;
  if((inputCIF->spgroup>=143)&&(inputCIF->spgroup<=194)){
    printf("\nSince the space group %d is hexagonal or trigonal, will make sure 1/3 and 2/3 have enough sig figs. \n",inputCIF->spgroup);
    for(i=0;i<inputCIF->total_site;i++){

      diff = fabs(inputCIF->xred[i]-one_third);
      //printf("diff for x and 0.333333: %lf\n",diff);
      if(diff<tolerance){
        //printf("Adjusting!\n");
        inputCIF->xred[i]=one_third;
        //printf("Adjusted value: %lf. \n",inputCIF->xred[i]);
      }

      diff = fabs(inputCIF->yred[i]-one_third);
      //printf("diff for y and 0.333333: %lf\n",diff);
      if(diff<tolerance){
        //printf("Adjusting!\n");
        inputCIF->yred[i]=one_third;
        //printf("Adjusted value: %lf. \n",inputCIF->yred[i]);
      }
      
      diff = fabs(inputCIF->zred[i]-one_third);
      //printf("diff for z and 0.333333: %lf\n",diff);
      if(diff<tolerance){
        //printf("Adjusting!\n");
        inputCIF->zred[i]=one_third;
        //printf("Adjusted value: %lf. \n",inputCIF->zred[i]);
      }
      
      diff = fabs(inputCIF->xred[i]-two_third);
      //printf("diff for x and 0.666667: %lf\n",diff);
      if(diff<tolerance){
        //printf("Adjusting!\n");
        inputCIF->xred[i]=two_third;
        //printf("Adjusted value: %lf. \n",inputCIF->xred[i]);
      }

      diff = fabs(inputCIF->yred[i]-two_third);
      //printf("diff for y and 0.666667: %lf\n",diff);
      if(diff<tolerance){
        //printf("Adjusting!\n");
        inputCIF->yred[i]=two_third;
        //printf("Adjusted value: %lf. \n",inputCIF->yred[i]);
      }

      diff = fabs(inputCIF->zred[i]-two_third);
      //printf("diff for z and 0.666667: %lf\n",diff);
      if(diff<tolerance){
        //printf("Adjusting!\n");
        inputCIF->zred[i]=two_third;
        //printf("Adjusted value: %lf. \n",inputCIF->zred[i]);
      }
    }
  }
}

int find_spg_Brav_Latt(char center, int spgroup){
  //LOGIC: find spgroup and corresponding Bravis Lattice. Return value of Brav Latt if found; if not found, return -1 at end.
  //1: C-centered Monoclinic; 2: C-centered Orthorhombic; 3: I-centered Orthorhombic; 4: F-centered Orthorhombic;
  //5: I-centered tetragonal; 6: F-centered cubic; 7: I-centered cubic; 
  if((spgroup>=3)&&(spgroup<=15)){
    printf("Monoclinic: ");
    if((center=='C')||(center=='c')){
      printf("C-centered\n");
      return 1;
    }
    else if((center=='A')||(center=='a')){
      printf("A-centered\n");
      return 11;
    }
    else if((center=='B')||(center=='b')){
      printf("B-centered\n");
      return 12;
    }
  }else if((spgroup>=16)&&(spgroup<=74)){
    printf("Orthorhombic: ");
    if((center=='C')||(center=='c')){
      printf("C-centered\n");
      return 2;
    }else if((center=='A')||(center=='a')){
      printf("A-centered\n");
      return 21;
    }else if((center=='B')||(center=='b')){
      printf("B-centered\n");
      return 22;
    }else if((center=='I')||(center=='i')){
      printf("I-centered\n");
      return 3;
    }else if((center=='F')||(center=='f')){
      printf("F-centered\n");
      return 4;
    }
  }else if((spgroup>=75)&&(spgroup<=142)){
    printf("Tetragonal: ");
    if((center=='I')||(center=='i')){
      printf("I-centered\n");
      return 5;
    }
  }else if((spgroup>=143)&&(spgroup<=167)){
    printf("Trigonal: ");
    if((center=='R')||(center=='r')){
      printf("R-centered\n");
      return 31;
    }  
  }else if((spgroup>=195)&&(spgroup<=230)){
    printf("Cubic: ");
    if((center=='F')||(center=='f')){
      printf("F-centered\n");
      return 6;
    }else if((center=='I')||(center=='i')){
      printf("I-centered\n");
      return 7;
    }
  }
  printf("Crystal structure is NOT of the 14 Bravis lattices!\n\n");
  return -1;
}

void generate_ngkpt(struct CIF * inputCIF){
  //declare variables
  double volume;
  double a_rad, b_rad, r_rad; //Radian of alpha, beta, gamma
  double cos_a, cos_b, cos_r; //COS of alpha, beta, gamma
  double sin_a, sin_b, sin_r; //SIN of alpha, beta, gamma
  double len_a, len_b, len_c; //length in angstrom of alpha, beta, gamma
  double prim_a, prim_b, prim_c; //length of primitive vectors
  double factor1, factor2; //factor1=cosa^2+cosb^2+cosr^2; factor2=2*cosa*cosb*cosr
  double a_star, b_star, c_star;
  double ngkpt_factor[5]; //Used to generate ngkpt
  double nkpt_a, nkpt_b, nkpt_c;
  int Brav_Latt_type=0;
  double Vec[3][3]; //three vectors
  double V1,V2,V3;
  int i; //Looper
  //printf("Inside function: generate_ngkpt. \n");
  //printf("Alpha = %lf, beta = %lf, gamma = %lf. \n", inputCIF->ang_alpha,inputCIF->ang_beta,inputCIF->ang_gamma);
  inputCIF->kpt_mesh_accepted=0;
  ngkpt_factor[0]=0.3;
  ngkpt_factor[1]=0.2;
  ngkpt_factor[2]=0.1;
  ngkpt_factor[3]=0.07;
  ngkpt_factor[4]=0.05;
  len_a=inputCIF->cell_a;
  len_b=inputCIF->cell_b;
  len_c=inputCIF->cell_c;
  a_rad=inputCIF->ang_alpha/180*PI;
  b_rad=inputCIF->ang_beta/180*PI;
  r_rad=inputCIF->ang_gamma/180*PI;
  cos_a=cos(a_rad);
  cos_b=cos(b_rad);
  cos_r=cos(r_rad);
  sin_a=sin(a_rad);
  sin_b=sin(b_rad);
  sin_r=sin(r_rad);
  //printf("cos_a = %lf, cos_b = %lf, cos_r = %lf. \n", cos_a,cos_b,cos_r);
  //printf("len_a = %lf, len_b = %lf, len_c = %lf. \n", len_a,len_b,len_c);
  //printf("space group center letter is: %c \n",inputCIF->spg_center_letter);
  //printf("Convert flag is : %d \n", convert_primitive_flag);
  //if(convert_primitive_flag==1) printf("primitive flag 1\n");
  //if(inputCIF->spg_center_letter!='P') printf("Not P\n");
  //if(inputCIF->spg_center_letter!='p') printf("Not p\n");
  if((convert_primitive_flag==1)&&(inputCIF->spg_center_letter!='P')&&(inputCIF->spg_center_letter!='p')&&(inputCIF->spg_center_letter!='R')&&(inputCIF->spg_center_letter!='r')){
    //printf("This is not a primitive cell, will convert\n");
    Brav_Latt_type=find_spg_Brav_Latt(inputCIF->spg_center_letter,inputCIF->spgroup);
    printf("Converting to primitive cell. \n");
    switch (Brav_Latt_type){
      case 1:  //Monoclinic C
        Vec[0][0]=0.5*len_a*sin_b;
        Vec[0][1]=-0.5*len_a*cos_b;
        Vec[0][2]=0.5*-len_b;
        Vec[1][0]=0.5*len_a*sin_b;
        Vec[1][1]=-0.5*len_a*cos_b;
        Vec[1][2]=0.5*len_b;
        Vec[2][0]=0;
        Vec[2][1]=0;
        Vec[2][2]=len_c;
        V1=(Vec[0][0]*Vec[0][0])+(Vec[0][1]*Vec[0][1])+(Vec[0][2]*Vec[0][2]);
        V2=(Vec[1][0]*Vec[1][0])+(Vec[1][1]*Vec[1][1])+(Vec[1][2]*Vec[1][2]);
        V3=(Vec[2][0]*Vec[2][0])+(Vec[2][1]*Vec[2][1])+(Vec[2][2]*Vec[2][2]);
        break;
      case 2: //Orthorhombic C
        Vec[0][0]=0.5*len_a;
        Vec[0][1]=-0.5*len_b;
        Vec[0][2]=0;
        Vec[1][0]=0.5*len_a;
        Vec[1][1]=0.5*len_b;
        Vec[1][2]=0;
        Vec[2][0]=0;
        Vec[2][1]=0;
        Vec[2][2]=len_c;
        V1=(Vec[0][0]*Vec[0][0])+(Vec[0][1]*Vec[0][1])+(Vec[0][2]*Vec[0][2]);
        V2=(Vec[1][0]*Vec[1][0])+(Vec[1][1]*Vec[1][1])+(Vec[1][2]*Vec[1][2]);
        V3=(Vec[2][0]*Vec[2][0])+(Vec[2][1]*Vec[2][1])+(Vec[2][2]*Vec[2][2]); 
        break;
      case 21: //Orthorhombic A
        Vec[0][0]=len_a;
        Vec[0][1]=0;
        Vec[0][2]=0;
        Vec[1][0]=0;
        Vec[1][1]=0.5*len_b;
        Vec[1][2]=-0.5*len_c;
        Vec[2][0]=0;
        Vec[2][1]=0.5*len_b;
        Vec[2][2]=0.5*len_c;
        V1=(Vec[0][0]*Vec[0][0])+(Vec[0][1]*Vec[0][1])+(Vec[0][2]*Vec[0][2]);
        V2=(Vec[1][0]*Vec[1][0])+(Vec[1][1]*Vec[1][1])+(Vec[1][2]*Vec[1][2]);
        V3=(Vec[2][0]*Vec[2][0])+(Vec[2][1]*Vec[2][1])+(Vec[2][2]*Vec[2][2]);
        break;
      case 22: //Orthorhombic B
        Vec[0][0]=-0.5*len_a;
        Vec[0][1]=0;
        Vec[0][2]=0.5*len_c;
        Vec[1][0]=0;
        Vec[1][1]=len_b;
        Vec[1][2]=0;
        Vec[2][0]=0.5*len_a;
        Vec[2][1]=0;
        Vec[2][2]=0.5*len_c;
        V1=(Vec[0][0]*Vec[0][0])+(Vec[0][1]*Vec[0][1])+(Vec[0][2]*Vec[0][2]);
        V2=(Vec[1][0]*Vec[1][0])+(Vec[1][1]*Vec[1][1])+(Vec[1][2]*Vec[1][2]);
        V3=(Vec[2][0]*Vec[2][0])+(Vec[2][1]*Vec[2][1])+(Vec[2][2]*Vec[2][2]);
        break;
      case 3: //Orthorhombic I
        Vec[0][0]=0.5*len_a;
        Vec[0][1]=0.5*len_b;
        Vec[0][2]=0.5*len_c;
        Vec[1][0]=-0.5*len_a;
        Vec[1][1]=-0.5*len_b;
        Vec[1][2]=0.5*len_c;
        Vec[2][0]=0.5*len_a;
        Vec[2][1]=-0.5*len_b;
        Vec[2][2]=-0.5*len_c;
        V1=(Vec[0][0]*Vec[0][0])+(Vec[0][1]*Vec[0][1])+(Vec[0][2]*Vec[0][2]);
        V2=(Vec[1][0]*Vec[1][0])+(Vec[1][1]*Vec[1][1])+(Vec[1][2]*Vec[1][2]);
        V3=(Vec[2][0]*Vec[2][0])+(Vec[2][1]*Vec[2][1])+(Vec[2][2]*Vec[2][2]); 
        break;
      case 4: //Orthorhombic F
        Vec[0][0]=0.5*len_a;
        Vec[0][1]=0;
        Vec[0][2]=0.5*len_c;
        Vec[1][0]=0;
        Vec[1][1]=-0.5*len_b;
        Vec[1][2]=0.5*len_c;
        Vec[2][0]=0.5*len_a;
        Vec[2][1]=-0.5*len_b;
        Vec[2][2]=0;
        V1=(Vec[0][0]*Vec[0][0])+(Vec[0][1]*Vec[0][1])+(Vec[0][2]*Vec[0][2]);
        V2=(Vec[1][0]*Vec[1][0])+(Vec[1][1]*Vec[1][1])+(Vec[1][2]*Vec[1][2]);
        V3=(Vec[2][0]*Vec[2][0])+(Vec[2][1]*Vec[2][1])+(Vec[2][2]*Vec[2][2]);
        break;
      case 5: //Tetragonal I
        Vec[0][0]=-0.5*len_a;
        Vec[0][1]=0.5*len_a;
        Vec[0][2]=0.5*len_c;
        Vec[1][0]=0.5*len_a;
        Vec[1][1]=-0.5*len_a;
        Vec[1][2]=0.5*len_c;
        Vec[2][0]=0.5*len_a;
        Vec[2][1]=0.5*len_a;
        Vec[2][2]=-0.5*len_c;
        V1=(Vec[0][0]*Vec[0][0])+(Vec[0][1]*Vec[0][1])+(Vec[0][2]*Vec[0][2]);
        V2=(Vec[1][0]*Vec[1][0])+(Vec[1][1]*Vec[1][1])+(Vec[1][2]*Vec[1][2]);
        V3=(Vec[2][0]*Vec[2][0])+(Vec[2][1]*Vec[2][1])+(Vec[2][2]*Vec[2][2]);
        break;
      case 6: //Cubic F
        Vec[0][0]=0;
        Vec[0][1]=0.5*len_a;
        Vec[0][2]=0.5*len_a;
        Vec[1][0]=0.5*len_a;
        Vec[1][1]=0;
        Vec[1][2]=0.5*len_a;
        Vec[2][0]=0.5*len_a;
        Vec[2][1]=0.5*len_a;
        Vec[2][2]=0;
        V1=(Vec[0][0]*Vec[0][0])+(Vec[0][1]*Vec[0][1])+(Vec[0][2]*Vec[0][2]);
        V2=(Vec[1][0]*Vec[1][0])+(Vec[1][1]*Vec[1][1])+(Vec[1][2]*Vec[1][2]);
        V3=(Vec[2][0]*Vec[2][0])+(Vec[2][1]*Vec[2][1])+(Vec[2][2]*Vec[2][2]);
        break;
      case 7: //Cubic I
        Vec[0][0]=-0.5*len_a;
        Vec[0][1]=0.5*len_a;
        Vec[0][2]=0.5*len_a;
        Vec[1][0]=0.5*len_a;
        Vec[1][1]=-0.5*len_a;
        Vec[1][2]=0.5*len_a;
        Vec[2][0]=0.5*len_a;
        Vec[2][1]=0.5*len_a;
        Vec[2][2]=-0.5*len_a;
        V1=(Vec[0][0]*Vec[0][0])+(Vec[0][1]*Vec[0][1])+(Vec[0][2]*Vec[0][2]);
        V2=(Vec[1][0]*Vec[1][0])+(Vec[1][1]*Vec[1][1])+(Vec[1][2]*Vec[1][2]);
        V3=(Vec[2][0]*Vec[2][0])+(Vec[2][1]*Vec[2][1])+(Vec[2][2]*Vec[2][2]);
        break;
      case -1:
        printf("The centering in your CIF file is not one of the 14 Bravis Lattice. \n\n");
        convert_FAIL=1;
        break;
      default:
        printf("Unexpected value of Brav_Latt_Type: %d. \n\n",Brav_Latt_type);
        break;
    }
    prim_a=sqrt(V1);
    prim_b=sqrt(V2);
    prim_c=sqrt(V3);
    printf("Values of primitive cell vectors a, b and c: %lf %lf %lf Angstroms\n",prim_a,prim_b,prim_c);
    a_star=2*PI/prim_a;
    b_star=2*PI/prim_b;
    c_star=2*PI/prim_c;
  }else{
    printf("You select not to convert to primitive cell OR the unit cell is primitive.\n");
    factor1=cos_a*cos_a+cos_b*cos_b+cos_r*cos_r;
    factor2=2*cos_a*cos_b*cos_r;
    volume=len_a*len_b*len_c*sqrt(1-factor1+factor2);
    //printf("The cell volume is: %lf \n",volume);
    inputCIF->cell_volume=volume;
    a_star=len_b*len_c*sin_a/volume*2*PI;
    b_star=len_c*len_a*sin_b/volume*2*PI;
    c_star=len_a*len_b*sin_r/volume*2*PI;
  }
  for(i=0;i<5;i++){
    nkpt_a=ceil(a_star/ngkpt_factor[i]);
    nkpt_b=ceil(b_star/ngkpt_factor[i]);
    nkpt_c=ceil(c_star/ngkpt_factor[i]);
    inputCIF->ngkpt[i][0]=nkpt_a;
    inputCIF->ngkpt[i][1]=nkpt_b;
    inputCIF->ngkpt[i][2]=nkpt_c;
    if(nkpt_a*nkpt_b*nkpt_c<=1730){
      inputCIF->kpt_mesh_accepted++;
    }
    //printf("Value of ngkpt for ngkpt%i: %i, %i, %i \n",(i+1),inputCIF->ngkpt[i][0],inputCIF->ngkpt[i][1],inputCIF->ngkpt[i][2]);
  }
  printf("Values of a*, b* and c*: %lf, %lf, %lf 1/Angstrom. \n\n",a_star,b_star,c_star);
  //printf("Number of accepted k-point mesh: %d. \n", inputCIF->kpt_mesh_accepted);
}

void readCIF(struct CIF * inputCIF, FILE * CIFfile){
    //declare variables//
    char dummy[100];
    char HMtemp[100];
    int stop=0;
    int i=0;
    int all_site_finished=0;
    char dummy2[100]; //used in reading through site
    char temp[50][100];
    int site_count=0;
    int ntype_count=0;
    int HM_count=0;
    int str_length, dummy3_len, dummy4_len, HMtemp_len;
    int xred_keyword_position=0, yred_keyword_position=0, zred_keyword_position=0, occ_keyword_position=0, symmul_keyword_position, label_keyword_position=0;
    int correct_loop=0; //0=incorrect, 1=correct
    int keyword_count;
    int dummy2_length;
    char dummy3[100];
    char dummy4[100];
    /*Reading CIF
      logic: when reach EOF, set stop = 1
      Otherwise, read and find keywords*/
    while(stop == 0){
      if(fscanf(CIFfile,"%s",dummy) != 1){
        stop = 1;
        continue;
      }
      //printf("The dummy word read is: %s \n", dummy);
      if(strcmp(dummy,"_cell_length_a")==0){
        fscanf(CIFfile,"%s",&inputCIF->cell_a_str);
        //printf("Cell parameter a is: %s\n",inputCIF->cell_a_str);
      }
      if(strcmp(dummy,"_cell_length_b")==0){
        fscanf(CIFfile,"%s",&inputCIF->cell_b_str);
        //printf("Cell parameter b is: %s\n",inputCIF->cell_b_str);
      }
      if(strcmp(dummy,"_cell_length_c")==0){
        fscanf(CIFfile,"%s",&inputCIF->cell_c_str);
        //printf("Cell parameter c is: %s\n",inputCIF->cell_c_str);
      }
      if(strcmp(dummy,"_cell_angle_alpha")==0){
        fscanf(CIFfile,"%s",&inputCIF->ang_alpha_str);
        //printf("Alpha: %s\n",inputCIF->ang_alpha_str);
      }
      if(strcmp(dummy,"_cell_angle_beta")==0){
        fscanf(CIFfile,"%s",&inputCIF->ang_beta_str);
        //printf("Beta: %s\n",inputCIF->ang_beta_str);
      }
      if(strcmp(dummy,"_cell_angle_gamma")==0){
        fscanf(CIFfile,"%s",&inputCIF->ang_gamma_str);
        //printf("Gamma: %s\n",inputCIF->ang_gamma_str);
      }
      if(strcmp(dummy,"_cell_formula_units_Z")==0){
        fscanf(CIFfile,"%d",&inputCIF->Z_unit);
        //printf("Z value: %d\n",inputCIF->Z_unit);
      }
      /*if(strcmp(dummy,"_symmetry_space_group_name_H-M")==0){
        fscanf(CIFfile,"%s",&inputCIF->spgroup_center);
        //printf("Space group centering: %s\n",inputCIF->spgroup_center);
      }*/
      if((strcmp(dummy,"_symmetry_Int_Tables_number")==0)||(strcmp(dummy,"_space_group_IT_number")==0)||(strcmp(dummy,"_symmetry_int_tables_number")==0)){
        fscanf(CIFfile,"%d",&inputCIF->spgroup);
        printf("Space group: %d\n",inputCIF->spgroup);
      }
      if(strcmp(dummy,"_chemical_formula_sum")==0){
        while(1==1){
          fscanf(CIFfile,"%s",&inputCIF->chemical_formula[ntype_count]);
          //printf("Element %d: %s\n",(ntype_count+1),inputCIF->chemical_formula[ntype_count]);
          str_length=strlen(inputCIF->chemical_formula[ntype_count]);
          ntype_count++;
          //printf("Current ntype: %d \n",ntype_count);
          if(inputCIF->chemical_formula[ntype_count-1][str_length-1]==39){
            break;
          }
          fscanf(CIFfile,"%s",dummy3);
          dummy3_len=strlen(dummy3);
          if(dummy3[0]==95){
            fseek(CIFfile,-dummy3_len,SEEK_CUR);
            break;
          }
          fseek(CIFfile,-dummy3_len,SEEK_CUR);
        }
      ntype=ntype_count; 
      }
      if((strcmp(dummy,"_symmetry_space_group_name_H-M")==0)||(strcmp(dummy,"_space_group_name_H-M_alt")==0)){
        while(1==1){
          fscanf(CIFfile,"%s",&inputCIF->HM_FullspgName[HM_count]);
          //printf("H-M Full name word %d: %s \n",(HM_count+1),inputCIF->HM_FullspgName[HM_count]);
          str_length=strlen(inputCIF->HM_FullspgName[HM_count]);
          HM_count++;
          //printf("Current ntype: %d \n",ntype_count);
          if(inputCIF->HM_FullspgName[HM_count-1][str_length-1]==39){
            break;
          }
          fscanf(CIFfile,"%s",dummy4);
          dummy4_len=strlen(dummy4);
          if(dummy4[0]==95){
            fseek(CIFfile,-dummy4_len,SEEK_CUR);
            break;
          }
          fseek(CIFfile,-dummy4_len,SEEK_CUR);
        }
        //printf("HM_count: %d \n",HM_count);
        strcpy(HMtemp,inputCIF->HM_FullspgName[0]);
        if(HM_count>4) HM_count=4;
        for(i=1;i<HM_count;i++){
          strcat(HMtemp,inputCIF->HM_FullspgName[i]);
          //printf("HMtemp is: %s \n", HMtemp);
        }
        HMtemp_len=strlen(HMtemp);
        if(HMtemp[0]==39 && HMtemp[HMtemp_len-1]==39){
          for(i=1;i<HMtemp_len-1;i++){
            inputCIF->HM_FullspgName_oneword[i-1]=HMtemp[i];
          }
        }else if(HMtemp[0]==39){
          for(i=1;i<HMtemp_len;i++){
            inputCIF->HM_FullspgName_oneword[i-1]=HMtemp[i];
          }
        }else if(HMtemp[HMtemp_len-1]==39){
          for(i=0;i<HMtemp_len-1;i++){
            inputCIF->HM_FullspgName_oneword[i]=HMtemp[i];
          }
        }else{
          strcpy(inputCIF->HM_FullspgName_oneword,HMtemp);
        }
        strcpy(inputCIF->spgroup_center,inputCIF->HM_FullspgName[0]);
        printf("Full space group name: %s \n",inputCIF->HM_FullspgName_oneword);
      }
      if((strcmp(dummy,"loop_")==0)&&(correct_loop==0)){
        keyword_count=0;
        //printf("Entering a loop in the CIF file.\n");
        while(1==1){
          fscanf(CIFfile,"%s",dummy2);
          if(dummy2[0]==95){
            keyword_count++;
            //printf("The %d keyword is: %s. \n",keyword_count,dummy2);
            if(strcmp(dummy2,"_atom_site_label")==0){
              label_keyword_position=keyword_count;
              //printf("The %d keyword is label. \n",label_keyword_position);
            }else if(strcmp(dummy2,"_atom_site_fract_x")==0){
              xred_keyword_position=keyword_count;
              //printf("The %d keyword is xred. \n",xred_keyword_position);
            }else if(strcmp(dummy2,"_atom_site_fract_y")==0){
              yred_keyword_position=keyword_count;
              //printf("The %d keyword is yred. \n",yred_keyword_position);
            }else if(strcmp(dummy2,"_atom_site_fract_z")==0){
              zred_keyword_position=keyword_count;
              //printf("The %d keyword is zred. \n",zred_keyword_position);
            }else if(strcmp(dummy2,"_atom_site_occupancy")==0){
              occ_keyword_position=keyword_count;
              //printf("The %d keyword is occ. \n",occ_keyword_position);
            }else if(strcmp(dummy2,"_atom_site_symmetry_multiplicity")==0){
              symmul_keyword_position=keyword_count;
              //printf("The %d keyword is symmul. \n",symmul_keyword_position);
            }else if(strcmp(dummy2,"_atom_site_site_symmetry_multiplicity")==0){
              symmul_keyword_position=keyword_count;
              //printf("The %d keyword is symmul. \n",symmul_keyword_position);
            }
          }else{
            //printf("There are %d keywords in this loop in total. \n", keyword_count);
            dummy2_length=strlen(dummy2);
            fseek(CIFfile,-dummy2_length,SEEK_CUR);
            //printf("DONE reading the keyword fields for the loop! \n");
            break;
          }
        }
        if((xred_keyword_position!=0)&&(yred_keyword_position!=0)&&(zred_keyword_position!=0)&&(occ_keyword_position!=0)&&(symmul_keyword_position!=0)&&(label_keyword_position!=0)){
          correct_loop=1;
          printf("Wyckoff position information found! \n");
        }
        if(correct_loop==1){
          //printf("Reading data from the correct loop. \n");
          while(1==1){
          //while(fscanf(CIFfile,"%s",dummy2)!=EOF){
            //fscanf(CIFfile,"%s",dummy2);
            if((fscanf(CIFfile,"%s",dummy2)!=EOF)&&(dummy2[0]!=95)&&(strcmp(dummy2,"loop_")!=0)&&(strcmp(dummy2,"#End")!=0)){
              //printf("dummy 2 is : %s. \n",dummy2);
              printf("Processing data for site %d. \n",site_count+1);
              dummy2_length=strlen(dummy2);
              fseek(CIFfile,-dummy2_length,SEEK_CUR);
              for(i=0;i<keyword_count;i++){
                fscanf(CIFfile,"%s",temp[i]); 
              }
              sscanf(temp[label_keyword_position-1],"%s",&inputCIF->site_identity[site_count]);
              //printf("site identity for site %d is :%s \n", site_count+1, inputCIF->site_identity[site_count]);
              sscanf(temp[symmul_keyword_position-1],"%d",&inputCIF->site_eqv[site_count]);
              //printf("site eqv for site %d is :%d \n", site_count, inputCIF->site_eqv[site_count]);
              sscanf(temp[xred_keyword_position-1],"%s",&inputCIF->xred_str[site_count]);
              //printf("xred for site %d is :%s \n", site_count+1, inputCIF->xred_str[site_count]);
              sscanf(temp[yred_keyword_position-1],"%s",&inputCIF->yred_str[site_count]);
              //printf("yred for site %d is :%s \n", site_count+1, inputCIF->yred_str[site_count]);
              sscanf(temp[zred_keyword_position-1],"%s",&inputCIF->zred_str[site_count]);
              //printf("zred for site %d is :%s \n", site_count+1, inputCIF->zred_str[site_count]);   
              sscanf(temp[occ_keyword_position-1],"%s",&inputCIF->occupancy_str[site_count]);
              //printf("occ for site %d is :%s \n", site_count+1, inputCIF->occupancy_str[site_count]);
              site_count++;
            }else{
              //printf("The current dummy2 is: %s. \n",dummy2);
              printf("There are %d sites in total. \n",site_count);
              inputCIF->total_site=site_count;
              dummy2_length=strlen(dummy2);
              fseek(CIFfile,-dummy2_length,SEEK_CUR);
              printf("DONE processing each position \n");
              break;
            }
          }
        }
      }
    }
    inputCIF->old_brvltt=0; //Set old_brvltt to 0: 
}

void read_abinitout(struct CIF * inputCIF, FILE * abinitout){
    //declare variables
    char dummy[100];
    int i=0;
    int stop=0;
    //Logic: need to update: acell, xred and possibly rprim.
    printf("Begin to read abinit output file from previous calculation.\n");
    while(stop == 0){
      if(fscanf(abinitout,"%s",dummy) != 1){
        stop = 1;
        continue;
      }
      if(strcmp(dummy,"acell") == 0){
        fscanf(abinitout,"%lf %lf %lf",&inputCIF->cell_a,&inputCIF->cell_b,&inputCIF->cell_c);
        //printf("Updated acell: %lf %lf %lf \n", inputCIF->cell_a,inputCIF->cell_b,inputCIF->cell_c);
      }
      if(strcmp(dummy,"xred") == 0){
        for(i=0;i<inputCIF->total_site;i++){
          fscanf(abinitout,"%s %s %s",&inputCIF->xred_upd[i],&inputCIF->yred_upd[i],&inputCIF->zred_upd[i]);
          //printf("Updated xred, y red, zred for site %d: %s %s %s \n",i, inputCIF->xred_upd[i],inputCIF->yred_upd[i],inputCIF->zred_upd[i]);
        }
      }
      if(strcmp(dummy,"typat") == 0){
        for(i=0;i<inputCIF->total_site;i++){
          fscanf(abinitout,"%d",&inputCIF->site_typat[i]);
        }
      }
      if(strcmp(dummy,"rprim") == 0){
        rprim_found=1;
        printf("Found dummy: %s \n",dummy);
        fscanf(abinitout,"%s %s %s",&inputCIF->rprim_str[0][0],&inputCIF->rprim_str[0][1],&inputCIF->rprim_str[0][2]);
        fscanf(abinitout,"%s %s %s",&inputCIF->rprim_str[1][0],&inputCIF->rprim_str[1][1],&inputCIF->rprim_str[1][2]);
        fscanf(abinitout,"%s %s %s",&inputCIF->rprim_str[2][0],&inputCIF->rprim_str[2][1],&inputCIF->rprim_str[2][2]);
        //printf("rprim1: %s %s %s \n",&inputCIF->rprim_str[0][0],&inputCIF->rprim_str[0][1],&inputCIF->rprim_str[0][2]);
        //printf("rprim2: %s %s %s \n",&inputCIF->rprim_str[1][0],&inputCIF->rprim_str[1][1],&inputCIF->rprim_str[1][2]);
        //printf("rprim3: %s %s %s \n",&inputCIF->rprim_str[2][0],&inputCIF->rprim_str[2][1],&inputCIF->rprim_str[2][2]);
      }
      if(strcmp(dummy,"kptrlatt") == 0){
        //printf("Found dummy: %s \n",dummy);
        for(i=0;i<9;i++){
          fscanf(abinitout,"%d",&inputCIF->kptrlatt[i]);
          //printf("kptrlatt[%d]: %d \n",(i+1),inputCIF->kptrlatt[i]);
        }
      }
      if(strcmp(dummy,"nkpt") == 0){
        //printf("Found dummy: %s \n",dummy);
        fscanf(abinitout,"%d",&inputCIF->nkpt);
        //printf("kptrlatt[%d]: %s \n",(i+1),inputCIF->kptrlatt_str[i]);
      } 
    }
    
    printf("Done reading abinit output file. \n\n");
}

void read_old_in(struct CIF * inputCIF, FILE * old_in_file){
    //Logic: if an old calculation exist, then the shiftk and spgorig will be read to use in the subsequent calculation.
    int stop=0;
    char dummy[100];
    int natrd_found=0;
    int i=0;

    inputCIF->old_brvltt=0;
    printf("Entering read old abinit in file function \n");
    //Start by reading the natrd from the input file and then reset the pointer.
    printf("Now trying to retrive the natrd, which is inputCIF->total_site from the old in file.\n");
    while(natrd_found == 0){
      if(fscanf(old_in_file,"%s",dummy) != 1){
        natrd_found = 1;
        continue;
      }
      if(strcmp(dummy,"natrd") == 0){
        fscanf(old_in_file,"%d",&inputCIF->total_site);
        printf("natrd: %d \n", inputCIF->total_site);
      }
    }
    printf("DONE.\n");
    //reset pointer
    rewind(old_in_file);
    //for(i=0;i<10;i++){
      //fscanf(old_in_file,"%s",dummy);
      //printf("Word read %d: %s \n",i,dummy);
    //}
    //start to process the rest of the file
    while(stop == 0){
      if(fscanf(old_in_file,"%s",dummy) != 1){
        stop = 1;
        continue;
      }
      if(strcmp(dummy,"angdeg") == 0){
        fscanf(old_in_file,"%lf %lf %lf",&inputCIF->ang_alpha,&inputCIF->ang_beta,&inputCIF->ang_gamma);
        //printf("Shiftk: %lf %lf %lf \n", inputCIF->shiftk[0],inputCIF->shiftk[1],inputCIF->shiftk[2]);
      }
      if(strcmp(dummy,"brvltt") == 0){
        fscanf(old_in_file,"%d",&inputCIF->old_brvltt);
        //printf("Shiftk: %lf %lf %lf \n", inputCIF->shiftk[0],inputCIF->shiftk[1],inputCIF->shiftk[2]);
      }
      if(strcmp(dummy,"chkprim") == 0){
        inputCIF->old_chkprim=1; //Found chkprim in the old input file
        //printf("Shiftk: %lf %lf %lf \n", inputCIF->shiftk[0],inputCIF->shiftk[1],inputCIF->shiftk[2]);
      }
      if(strcmp(dummy,"ecut") == 0){
        fscanf(old_in_file,"%d",&inputCIF->ecut_final);
        inputCIF->ecut_lambda_bohr=2*PI/sqrt(2*inputCIF->ecut_final);
        //printf("Shiftk: %lf %lf %lf \n", inputCIF->shiftk[0],inputCIF->shiftk[1],inputCIF->shiftk[2]);
      }
      if(strcmp(dummy,"znucl") == 0){
        for(i=0;i<inputCIF->total_site;i++){
          fscanf(old_in_file,"%d",&inputCIF->element_No[i]);
          //printf("Shiftk: %lf %lf %lf \n", inputCIF->shiftk[0],inputCIF->shiftk[1],inputCIF->shiftk[2]);
        }
      }
      if(strcmp(dummy,"ntypat") == 0){
        fscanf(old_in_file,"%d",&ntype);
        //printf("Shiftk: %lf %lf %lf \n", inputCIF->shiftk[0],inputCIF->shiftk[1],inputCIF->shiftk[2]);
      }
      if(strcmp(dummy,"nband") == 0){
        fscanf(old_in_file,"%d",&inputCIF->nband);
        //printf("Shiftk: %lf %lf %lf \n", inputCIF->shiftk[0],inputCIF->shiftk[1],inputCIF->shiftk[2]);
      }
      if(strcmp(dummy,"natom") == 0){
        fscanf(old_in_file,"%d",&inputCIF->natom);
        //printf("Shiftk: %lf %lf %lf \n", inputCIF->shiftk[0],inputCIF->shiftk[1],inputCIF->shiftk[2]);
      }
      if(strcmp(dummy,"shiftk") == 0){
        fscanf(old_in_file,"%lf %lf %lf",&inputCIF->shiftk[0],&inputCIF->shiftk[1],&inputCIF->shiftk[2]);
        //printf("Shiftk: %lf %lf %lf \n", inputCIF->shiftk[0],inputCIF->shiftk[1],inputCIF->shiftk[2]);
      }
      if(strcmp(dummy,"spgroup") == 0){
        fscanf(old_in_file,"%d",&inputCIF->spgroup);
        //printf("Shiftk: %lf %lf %lf \n", inputCIF->shiftk[0],inputCIF->shiftk[1],inputCIF->shiftk[2]);
      }
      if(strcmp(dummy,"spgorig") == 0){
        fscanf(old_in_file,"%d",&inputCIF->spgorig);
        spgorig_read=1;
        //printf("Shiftk: %lf %lf %lf \n", inputCIF->shiftk[0],inputCIF->shiftk[1],inputCIF->shiftk[2]);
      }
      if(strcmp(dummy,"spgaxor") == 0){
        fscanf(old_in_file,"%d",&inputCIF->spgaxor);
        spgaxor_relevant=5;
        //printf("Shiftk: %lf %lf %lf \n", inputCIF->shiftk[0],inputCIF->shiftk[1],inputCIF->shiftk[2]);
      }
      if(kmesh_flag==1){
        if(strcmp(dummy,"ngkpt") == 0){
        fscanf(old_in_file,"%d %d %d",&inputCIF->ngkpt_from_in[0],&inputCIF->ngkpt_from_in[1],&inputCIF->ngkpt_from_in[2]);
        //printf("Shiftk: %lf %lf %lf \n", inputCIF->shiftk[0],inputCIF->shiftk[1],inputCIF->shiftk[2]);
        }
      }
    }
    if(inputCIF->old_brvltt==-1){
      printf("The old abinit calculation has converted your unit cell to primitive cell.\n");
      printf("The new input file will be adjusted accordingly.\n");
      inputCIF->total_site=inputCIF->natom;
      inputCIF->spgroup=0;
    }
    for(i=0;i<ntype;i++){
      inputCIF->nelectrons_psp[i]=PSP_Assigner(inputCIF->element_No[i],inputCIF,i);
    }
    printf("Done reading old input file. \n\n");
}

char * ReducePT(char * input_Str){
    //declare variables
    int input_len;
    char * output_Str;
    char temp_Str[100];
    int i = 0;
    int count=0;
    for(i=0;i<=99;i++){
      temp_Str[i]='\0';
    }
    /*Logic: check for () in input string. If a character is not ()
    Then copy that character to output string. Use count to count number
    of characters already in the output string*/
    input_len=strlen(input_Str);
    //printf("Input string is: %s with length %d \n",input_Str, input_len);
    for(i=0;i<=input_len-1;i++){
      if( (input_Str[i] != '(') && (input_Str[i] != ')')){
        //printf("Input[%d]: %c \n", i, input_Str[i]);
        temp_Str[count] = input_Str[i];
        count++;
      }
    }
    //printf("temp_Str: %s \n", temp_Str);
    output_Str=temp_Str;
    //printf("Output string is: %s \n", output_Str); 
    return output_Str;
}

char * DOT(char * input_Str){
    //declare variables
    int input_len;
    char * output_Str;
    char temp_Str[100];
    int i = 0;
    int count=0;
    //printf("DOT received input_Str is: %s \n", input_Str);
    strcpy(temp_Str,input_Str);
    input_len=strlen(input_Str);
    if(temp_Str[input_len-1]=='.'){
        temp_Str[input_len]='0';
    }
    //printf("temp_Str: %s \n", temp_Str);
    output_Str=temp_Str;
    //printf("Output string is: %s \n", output_Str);
    return output_Str;
}

char * AdjustNegative(char * input_Str){
    //declare variables
    int input_len;
    char * output_Str;
    char temp_Str[100];
    int i = 0;
    int count=0;
    for(i=0;i<=99;i++){
      temp_Str[i]='\0';
    }
    input_len=strlen(input_Str);
    if( (input_Str[0]=='-') && (input_Str[1]=='.') ){
        temp_Str[0]='-';
        temp_Str[1]='0';
        for(i=2;i<input_len+1;i++){
          temp_Str[i]=input_Str[i-1];
        }
      output_Str=temp_Str;
    }else{
      output_Str=input_Str;
    }
    //printf("Output string is: %s \n", output_Str);
    return output_Str;
}


double NumberString(char * input_Str){
    char * str_after_ReducePT;
    char * str_after_DOT;
    char * str_after_AdjustNegative;
    char final_Str[100];
    char final_Str2[100];
    double final_number;
    int i;
    for(i=0;i<=99;i++){
      final_Str[i]='\0';
    }
    str_after_ReducePT=ReducePT(input_Str);
    //printf("str_after_ReducePT: %s \n", str_after_ReducePT);
    str_after_DOT=DOT(str_after_ReducePT);
    //printf("str_after_DOT: %s \n", str_after_DOT);
    strcpy(final_Str2,str_after_DOT);
    str_after_AdjustNegative=AdjustNegative(final_Str2);
    strcpy(final_Str,str_after_AdjustNegative);
    //printf("final str: %s \n", final_Str);
    sscanf(final_Str,"%lf",&final_number);
    //printf("Final number is: %lf \n",final_number);
    return final_number;
}

char * ExtractLetter(char * input_Str){
    //declare variables
    int input_len;
    char * output_Str;
    char temp_Str[100];
    int i = 0;
    int count=0;
    //printf("In ExtractLetter. \n");
    memset(temp_Str,0,100); 
    /*for(i=0;i<100;i++){
      temp_Str[i]='\0';
      output_Str[i]='\0';
    }*/
    
    //printf("Received input: %s \n",input_Str);
    input_len=strlen(input_Str);    
    for(i=0;i<input_len;i++){
        if(isalpha(input_Str[i])){
           temp_Str[count]=input_Str[i];
           count++;
        }
    }
    count=0;
    //printf("Temp is: %s \n", temp_Str);
    output_Str=temp_Str;
    //printf("Extracted words: %s \n",output_Str);
    return output_Str;
}

int AtomNoFinder(char * Ele_abbv){
    if(strcmp(Ele_abbv,"H")==0) return 1;
    if(strcmp(Ele_abbv,"He")==0) return 2;
    if(strcmp(Ele_abbv,"Li")==0) return 3;
    if(strcmp(Ele_abbv,"Be")==0) return 4;
    if(strcmp(Ele_abbv,"B")==0) return 5;
    if(strcmp(Ele_abbv,"C")==0) return 6;
    if(strcmp(Ele_abbv,"N")==0) return 7;
    if(strcmp(Ele_abbv,"O")==0) return 8;
    if(strcmp(Ele_abbv,"F")==0) return 9;
    if(strcmp(Ele_abbv,"Ne")==0) return 10;
    if(strcmp(Ele_abbv,"Na")==0) return 11;
    if(strcmp(Ele_abbv,"Mg")==0) return 12;
    if(strcmp(Ele_abbv,"Al")==0) return 13;
    if(strcmp(Ele_abbv,"Si")==0) return 14;
    if(strcmp(Ele_abbv,"P")==0) return 15;
    if(strcmp(Ele_abbv,"S")==0) return 16;
    if(strcmp(Ele_abbv,"Cl")==0) return 17;
    if(strcmp(Ele_abbv,"Ar")==0) return 18;
    if(strcmp(Ele_abbv,"K")==0) return 19;
    if(strcmp(Ele_abbv,"Ca")==0) return 20;
    if(strcmp(Ele_abbv,"Sc")==0) return 21;
    if(strcmp(Ele_abbv,"Ti")==0) return 22;
    if(strcmp(Ele_abbv,"V")==0) return 23;
    if(strcmp(Ele_abbv,"Cr")==0) return 24;
    if(strcmp(Ele_abbv,"Mn")==0) return 25;
    if(strcmp(Ele_abbv,"Fe")==0) return 26;
    if(strcmp(Ele_abbv,"Co")==0) return 27;
    if(strcmp(Ele_abbv,"Ni")==0) return 28;
    if(strcmp(Ele_abbv,"Cu")==0) return 29;
    if(strcmp(Ele_abbv,"Zn")==0) return 30;
    if(strcmp(Ele_abbv,"Ga")==0) return 31;
    if(strcmp(Ele_abbv,"Ge")==0) return 32;
    if(strcmp(Ele_abbv,"As")==0) return 33;
    if(strcmp(Ele_abbv,"Se")==0) return 34;
    if(strcmp(Ele_abbv,"Br")==0) return 35;
    if(strcmp(Ele_abbv,"Kr")==0) return 36;
    if(strcmp(Ele_abbv,"Rb")==0) return 37;
    if(strcmp(Ele_abbv,"Sr")==0) return 38;
    if(strcmp(Ele_abbv,"Y")==0) return 39;
    if(strcmp(Ele_abbv,"Zr")==0) return 40;
    if(strcmp(Ele_abbv,"Nb")==0) return 41;
    if(strcmp(Ele_abbv,"Mo")==0) return 42;
    if(strcmp(Ele_abbv,"Tc")==0) return 43;
    if(strcmp(Ele_abbv,"Ru")==0) return 44;
    if(strcmp(Ele_abbv,"Rh")==0) return 45;
    if(strcmp(Ele_abbv,"Pd")==0) return 46;
    if(strcmp(Ele_abbv,"Ag")==0) return 47;
    if(strcmp(Ele_abbv,"Cd")==0) return 48;
    if(strcmp(Ele_abbv,"In")==0) return 49;
    if(strcmp(Ele_abbv,"Sn")==0) return 50;
    if(strcmp(Ele_abbv,"Sb")==0) return 51;
    if(strcmp(Ele_abbv,"Te")==0) return 52;
    if(strcmp(Ele_abbv,"I")==0) return 53;
    if(strcmp(Ele_abbv,"Xe")==0) return 54;
    if(strcmp(Ele_abbv,"Cs")==0) return 55;
    if(strcmp(Ele_abbv,"Ba")==0) return 56;
    if(strcmp(Ele_abbv,"La")==0) return 57;
    if(strcmp(Ele_abbv,"Ce")==0) return 58;
    if(strcmp(Ele_abbv,"Pr")==0) return 59;
    if(strcmp(Ele_abbv,"Nd")==0) return 60;
    if(strcmp(Ele_abbv,"Pm")==0) return 61;
    if(strcmp(Ele_abbv,"Sm")==0) return 62;
    if(strcmp(Ele_abbv,"Eu")==0) return 63;
    if(strcmp(Ele_abbv,"Gd")==0) return 64;
    if(strcmp(Ele_abbv,"Tb")==0) return 65;
    if(strcmp(Ele_abbv,"Dy")==0) return 66;
    if(strcmp(Ele_abbv,"Ho")==0) return 67;
    if(strcmp(Ele_abbv,"Er")==0) return 68;
    if(strcmp(Ele_abbv,"Tm")==0) return 69;
    if(strcmp(Ele_abbv,"Yb")==0) return 70;
    if(strcmp(Ele_abbv,"Lu")==0) return 71;
    if(strcmp(Ele_abbv,"Hf")==0) return 72;
    if(strcmp(Ele_abbv,"Ta")==0) return 73;
    if(strcmp(Ele_abbv,"W")==0) return 74;
    if(strcmp(Ele_abbv,"Re")==0) return 75;
    if(strcmp(Ele_abbv,"Os")==0) return 76;
    if(strcmp(Ele_abbv,"Ir")==0) return 77;
    if(strcmp(Ele_abbv,"Pt")==0) return 78;
    if(strcmp(Ele_abbv,"Au")==0) return 79;
    if(strcmp(Ele_abbv,"Hg")==0) return 80;
    if(strcmp(Ele_abbv,"Tl")==0) return 81;
    if(strcmp(Ele_abbv,"Pb")==0) return 82;
    if(strcmp(Ele_abbv,"Bi")==0) return 83;
    if(strcmp(Ele_abbv,"Po")==0) return 84;
    if(strcmp(Ele_abbv,"At")==0) return 85;
    if(strcmp(Ele_abbv,"Rn")==0) return 86;
    if(strcmp(Ele_abbv,"Fr")==0) return 87;
    if(strcmp(Ele_abbv,"Ra")==0) return 88;
    if(strcmp(Ele_abbv,"Ac")==0) return 89;
    if(strcmp(Ele_abbv,"Th")==0) return 90;
    if(strcmp(Ele_abbv,"Pa")==0) return 91;
    if(strcmp(Ele_abbv,"U")==0) return 92;
    if(strcmp(Ele_abbv,"Np")==0) return 93;
    if(strcmp(Ele_abbv,"Pu")==0) return 94;
    if(strcmp(Ele_abbv,"Am")==0) return 95;
    if(strcmp(Ele_abbv,"Cm")==0) return 96;
    if(strcmp(Ele_abbv,"Bk")==0) return 97;
    if(strcmp(Ele_abbv,"Cf")==0) return 98;
    if(strcmp(Ele_abbv,"Es")==0) return 99;
    if(strcmp(Ele_abbv,"Fm")==0) return 100;
    if(strcmp(Ele_abbv,"Md")==0) return 101;
    if(strcmp(Ele_abbv,"No")==0) return 102;
    if(strcmp(Ele_abbv,"Lr")==0) return 103;
    if(strcmp(Ele_abbv,"Rf")==0) return 104;
    if(strcmp(Ele_abbv,"Db")==0) return 105;
    if(strcmp(Ele_abbv,"Sg")==0) return 106;
    if(strcmp(Ele_abbv,"Bh")==0) return 107;
    if(strcmp(Ele_abbv,"Hs")==0) return 108;
    if(strcmp(Ele_abbv,"Mt")==0) return 109;
    if(strcmp(Ele_abbv,"Ds")==0) return 110;
    if(strcmp(Ele_abbv,"Rg")==0) return 111;
    if(strcmp(Ele_abbv,"Cn")==0) return 112;
    if(strcmp(Ele_abbv,"Nh")==0) return 113;
    if(strcmp(Ele_abbv,"Fl")==0) return 114;
    if(strcmp(Ele_abbv,"Mc")==0) return 115;
    if(strcmp(Ele_abbv,"Lv")==0) return 116;
    if(strcmp(Ele_abbv,"Ts")==0) return 117;
    if(strcmp(Ele_abbv,"Og")==0) return 118;

    return 666;
}
int PSP_Assigner(int Ele_No, struct CIF * inputCIF, int typeNo){
    //printf("Element number is: %d \n", Ele_No);
    inputCIF->ecut[typeNo]=1;
    if(scvo_flag==1){
      switch(Ele_No){ 
        case 1:
          strcpy(inputCIF->psp_name[typeNo],"1h.1.hgh");
          inputCIF->ecut[typeNo]=35;
          return 1;
        case 2:
          strcpy(inputCIF->psp_name[typeNo],"2he.2.hgh");
          inputCIF->ecut[typeNo]=110;
          return 2;
        case 3:
          strcpy(inputCIF->psp_name[typeNo],"3li.1.hgh");
          inputCIF->ecut[typeNo]=15;
          return 1;
        case 4:
          strcpy(inputCIF->psp_name[typeNo],"4be.2.hgh");
          inputCIF->ecut[typeNo]=20;
          return 2;
        case 5:
          strcpy(inputCIF->psp_name[typeNo],"5b.3.hgh");
          inputCIF->ecut[typeNo]=40;
          return 3;
        case 6:
          strcpy(inputCIF->psp_name[typeNo],"6c.4.hgh");
          inputCIF->ecut[typeNo]=65;
          return 4;
        case 7:
          strcpy(inputCIF->psp_name[typeNo],"7n.5.hgh");
          inputCIF->ecut[typeNo]=90;
          return 5;
        case 8:
          strcpy(inputCIF->psp_name[typeNo],"8o.6.hgh");
          inputCIF->ecut[typeNo]=130;
          return 6;
        case 9:
          strcpy(inputCIF->psp_name[typeNo],"9f.7.hgh");
          inputCIF->ecut[typeNo]=290;
          return 7;
        case 10:
          strcpy(inputCIF->psp_name[typeNo],"10ne.8.hgh");
          inputCIF->ecut[typeNo]=260;
          return 8;
        case 11:
          strcpy(inputCIF->psp_name[typeNo],"11na.1.hgh");
          inputCIF->ecut[typeNo]=5;
          return 1;
        case 12:
          strcpy(inputCIF->psp_name[typeNo],"12mg.2.hgh");
          inputCIF->ecut[typeNo]=10;
          return 2;
        case 13:
          strcpy(inputCIF->psp_name[typeNo],"13al.3.hgh");
          inputCIF->ecut[typeNo]=15;
          return 3;
        case 14:
          strcpy(inputCIF->psp_name[typeNo],"14si.4.hgh");
          inputCIF->ecut[typeNo]=30;
          return 4;
        case 15:
          strcpy(inputCIF->psp_name[typeNo],"15p.5.hgh");
          inputCIF->ecut[typeNo]=40;
          return 5;
        case 16:
          strcpy(inputCIF->psp_name[typeNo],"16s.6.hgh");
          inputCIF->ecut[typeNo]=40;
          return 6;
        case 17:
          strcpy(inputCIF->psp_name[typeNo],"17cl.7.hgh");
          inputCIF->ecut[typeNo]=55;
          return 7;
        case 18:
          strcpy(inputCIF->psp_name[typeNo],"18ar.8.hgh");
          inputCIF->ecut[typeNo]=60;
          return 8;
        case 19:
          strcpy(inputCIF->psp_name[typeNo],"19k.1.hgh");
          inputCIF->ecut[typeNo]=10;
          return 1;
        case 20:
          strcpy(inputCIF->psp_name[typeNo],"20ca.2.hgh");
          inputCIF->ecut[typeNo]=10;
          return 2;
        case 21:
          strcpy(inputCIF->psp_name[typeNo],"21sc.3.hgh");
          inputCIF->ecut[typeNo]=10;
          return 3;
        case 22:
          strcpy(inputCIF->psp_name[typeNo],"22ti.4.hgh");
          inputCIF->ecut[typeNo]=50;
          return 4;
        case 23:
          strcpy(inputCIF->psp_name[typeNo],"23v.5.hgh");
          inputCIF->ecut[typeNo]=50;
          return 5;
        case 24:
          strcpy(inputCIF->psp_name[typeNo],"24cr.6.hgh");
          inputCIF->ecut[typeNo]=55;
          return 6;
        case 25:
          strcpy(inputCIF->psp_name[typeNo],"25mn.7.hgh");
          inputCIF->ecut[typeNo]=70;
          return 7;
        case 26:
          strcpy(inputCIF->psp_name[typeNo],"26fe.8.hgh");
          inputCIF->ecut[typeNo]=80;
          return 8;
        case 27:
          strcpy(inputCIF->psp_name[typeNo],"27co.9.hgh");
          inputCIF->ecut[typeNo]=90;
          return 9;
        case 28:
          strcpy(inputCIF->psp_name[typeNo],"28ni.10.hgh");
          inputCIF->ecut[typeNo]=95;
          return 10;
        case 29:
          strcpy(inputCIF->psp_name[typeNo],"29cu.11.hgh");
          inputCIF->ecut[typeNo]=105;
          return 11;
        case 30:
          strcpy(inputCIF->psp_name[typeNo],"30zn.12.hgh");
          inputCIF->ecut[typeNo]=110;
          return 12;
        case 31:
          strcpy(inputCIF->psp_name[typeNo],"31ga.3.hgh");
          inputCIF->ecut[typeNo]=15;
          return 3;
        case 32:
          strcpy(inputCIF->psp_name[typeNo],"32ge.4.hgh");
          inputCIF->ecut[typeNo]=25;
          return 4;
        case 33:
          strcpy(inputCIF->psp_name[typeNo],"33as.5.hgh");
          inputCIF->ecut[typeNo]=30;
          return 5;
        case 34:
          strcpy(inputCIF->psp_name[typeNo],"34se.6.hgh");
          inputCIF->ecut[typeNo]=35;
          return 6;
        case 35:
          strcpy(inputCIF->psp_name[typeNo],"35br.7.hgh");
          inputCIF->ecut[typeNo]=40;
          return 7;
        case 36:
          strcpy(inputCIF->psp_name[typeNo],"36kr.8.hgh");
          inputCIF->ecut[typeNo]=40;
          return 8;
        case 37:
          strcpy(inputCIF->psp_name[typeNo],"37rb.1.hgh");
          inputCIF->ecut[typeNo]=5;
          return 1;
        case 38:
          strcpy(inputCIF->psp_name[typeNo],"38sr.2.hgh");
          inputCIF->ecut[typeNo]=10;
          return 2;
        case 39:
          strcpy(inputCIF->psp_name[typeNo],"39y.3.hgh");
          inputCIF->ecut[typeNo]=15;
          return 3;
        case 40:
          strcpy(inputCIF->psp_name[typeNo],"40zr.4.hgh");
          inputCIF->ecut[typeNo]=20;
          return 4;
        case 41:
          strcpy(inputCIF->psp_name[typeNo],"41nb.5.hgh");
          inputCIF->ecut[typeNo]=25;
          return 5;
        case 42:
          strcpy(inputCIF->psp_name[typeNo],"42mo.6.hgh");
          inputCIF->ecut[typeNo]=40;
          return 6;
        case 43:
          strcpy(inputCIF->psp_name[typeNo],"43tc.7.hgh");
          inputCIF->ecut[typeNo]=120;
          return 7;
        case 44:
          strcpy(inputCIF->psp_name[typeNo],"44ru.8.hgh");
          inputCIF->ecut[typeNo]=35;
          return 8;
        case 45:
          strcpy(inputCIF->psp_name[typeNo],"45rh.9.hgh");
          inputCIF->ecut[typeNo]=65;
          return 9;
        case 46:
          strcpy(inputCIF->psp_name[typeNo],"46pd.10.hgh");
          inputCIF->ecut[typeNo]=35;
          return 10;
        case 47:
          strcpy(inputCIF->psp_name[typeNo],"47ag.11.hgh");
          inputCIF->ecut[typeNo]=45;
          return 11;
        case 48:
          strcpy(inputCIF->psp_name[typeNo],"48cd.12.hgh");
          inputCIF->ecut[typeNo]=50;
          return 12;
        case 49:
          strcpy(inputCIF->psp_name[typeNo],"49in.3.hgh");
          inputCIF->ecut[typeNo]=10;
          return 3;
        case 50:
          strcpy(inputCIF->psp_name[typeNo],"50sn.4.hgh");
          inputCIF->ecut[typeNo]=20;
          return 4;
        case 51:
          strcpy(inputCIF->psp_name[typeNo],"51sb.5.hgh");
          inputCIF->ecut[typeNo]=25;
          return 5;
        case 52:
          strcpy(inputCIF->psp_name[typeNo],"52te.6.hgh");
          inputCIF->ecut[typeNo]=25;
          return 6;
        case 53:
          strcpy(inputCIF->psp_name[typeNo],"53i.7.hgh");
          inputCIF->ecut[typeNo]=50;
          return 7;
        case 54:
          strcpy(inputCIF->psp_name[typeNo],"54xe.8.hgh");
          inputCIF->ecut[typeNo]=30;
          return 8;
        case 55:
          strcpy(inputCIF->psp_name[typeNo],"55cs.1.hgh");
          inputCIF->ecut[typeNo]=5;
          return 1;
        case 56:
          strcpy(inputCIF->psp_name[typeNo],"56ba.2.hgh");
          inputCIF->ecut[typeNo]=10;
          return 2;
        case 57:
          strcpy(inputCIF->psp_name[typeNo],"57la.11.hgh");
          inputCIF->ecut[typeNo]=70;
          return 11;
        case 58:
          strcpy(inputCIF->psp_name[typeNo],"58ce.12.hgh");
          return 12;
        case 59:
          strcpy(inputCIF->psp_name[typeNo],"59pr.13.hgh");
          return 13;
        case 60:
          strcpy(inputCIF->psp_name[typeNo],"60nd.14.hgh");
          return 14;
        case 61:
          strcpy(inputCIF->psp_name[typeNo],"61pm.15.hgh");
          return 15;
        case 62:
          strcpy(inputCIF->psp_name[typeNo],"62sm.16.hgh");
          return 16;
        case 63:
          strcpy(inputCIF->psp_name[typeNo],"63eu.17.hgh");
          return 17;
        case 64:
          strcpy(inputCIF->psp_name[typeNo],"64gd.18.hgh");
          return 18;
        case 65:
          strcpy(inputCIF->psp_name[typeNo],"65tb.19.hgh");
          return 19;
        case 66:
          strcpy(inputCIF->psp_name[typeNo],"66dy.20.hgh");
          return 20;
        case 67:
          strcpy(inputCIF->psp_name[typeNo],"67ho.14.hgh");
          return 21;
        case 68:
          strcpy(inputCIF->psp_name[typeNo],"68er.22.hgh");
          return 22;
        case 69:
          strcpy(inputCIF->psp_name[typeNo],"69tm.23.hgh");
          return 23;
        case 70:
          strcpy(inputCIF->psp_name[typeNo],"70yb.24.hgh");
          return 24;
        case 71:
          strcpy(inputCIF->psp_name[typeNo],"71lu.25.hgh");
          return 25;
        case 72:
          strcpy(inputCIF->psp_name[typeNo],"72hf.12.hgh");
          inputCIF->ecut[typeNo]=35;
          return 12;
        case 73:
          strcpy(inputCIF->psp_name[typeNo],"73ta.5.hgh");
          inputCIF->ecut[typeNo]=20;
          return 5;
        case 74:
          strcpy(inputCIF->psp_name[typeNo],"74w.6.hgh");
          inputCIF->ecut[typeNo]=20;
          return 6;
        case 75:
          strcpy(inputCIF->psp_name[typeNo],"75re.7.hgh");
          inputCIF->ecut[typeNo]=25;
          return 7;
        case 76:
          strcpy(inputCIF->psp_name[typeNo],"76os.8.hgh");
          return 8;
        case 77:
          strcpy(inputCIF->psp_name[typeNo],"77ir.9.hgh");
          inputCIF->ecut[typeNo]=35;
          return 9;
        case 78:
          strcpy(inputCIF->psp_name[typeNo],"78pt.10.hgh");
          inputCIF->ecut[typeNo]=30;
          return 10;
        case 79:
          strcpy(inputCIF->psp_name[typeNo],"79au.11.hgh");
          inputCIF->ecut[typeNo]=35;
          return 11;
        case 80:
          strcpy(inputCIF->psp_name[typeNo],"80hg.12.hgh");
          inputCIF->ecut[typeNo]=45;
          return 12;
        case 81:
          strcpy(inputCIF->psp_name[typeNo],"81tl.3.hgh");
          inputCIF->ecut[typeNo]=15;
          return 3;
        case 82:
          strcpy(inputCIF->psp_name[typeNo],"82pb.4.hgh");
          inputCIF->ecut[typeNo]=15;
          return 4;
        case 83:
          strcpy(inputCIF->psp_name[typeNo],"83bi.5.hgh");
          inputCIF->ecut[typeNo]=20;
          return 5;
        case 84:
          strcpy(inputCIF->psp_name[typeNo],"84po.6.hgh");
          return 6;
        case 85:
          strcpy(inputCIF->psp_name[typeNo],"85at.7.hgh");
          return 7;
      }
    }
    if(scvo_flag==2){
      switch(Ele_No){ 
        case 1:
          strcpy(inputCIF->psp_name[typeNo],"1h.1.hgh");
          inputCIF->ecut[typeNo]=35;
          return 1;
        case 2:
          strcpy(inputCIF->psp_name[typeNo],"2he.2.hgh");
          inputCIF->ecut[typeNo]=110;
          return 2;
        case 3:
          strcpy(inputCIF->psp_name[typeNo],"3li.3.hgh");
          inputCIF->ecut[typeNo]=100;
          return 3;
        case 4:
          strcpy(inputCIF->psp_name[typeNo],"4be.4.hgh");
          inputCIF->ecut[typeNo]=170;
          return 4;
        case 5:
          strcpy(inputCIF->psp_name[typeNo],"5b.3.hgh");
          inputCIF->ecut[typeNo]=40;
          return 3;
        case 6:
          strcpy(inputCIF->psp_name[typeNo],"6c.4.hgh");
          inputCIF->ecut[typeNo]=65;
          return 4;
        case 7:
          strcpy(inputCIF->psp_name[typeNo],"7n.5.hgh");
          inputCIF->ecut[typeNo]=90;
          return 5;
        case 8:
          strcpy(inputCIF->psp_name[typeNo],"8o.6.hgh");
          inputCIF->ecut[typeNo]=130;
          return 6;
        case 9:
          strcpy(inputCIF->psp_name[typeNo],"9f.7.hgh");
          inputCIF->ecut[typeNo]=290;
          return 7;
        case 10:
          strcpy(inputCIF->psp_name[typeNo],"10ne.8.hgh");
          inputCIF->ecut[typeNo]=260;
          return 8;
        case 11:
          strcpy(inputCIF->psp_name[typeNo],"11na.9.hgh");
          inputCIF->ecut[typeNo]=335;
          return 9;
        case 12:
          strcpy(inputCIF->psp_name[typeNo],"12mg.10.hgh");
          inputCIF->ecut[typeNo]=430;
          return 10;
        case 13:
          strcpy(inputCIF->psp_name[typeNo],"13al.3.hgh");
          inputCIF->ecut[typeNo]=15;
          return 3;
        case 14:
          strcpy(inputCIF->psp_name[typeNo],"14si.4.hgh");
          inputCIF->ecut[typeNo]=30;
          return 4;
        case 15:
          strcpy(inputCIF->psp_name[typeNo],"15p.5.hgh");
          inputCIF->ecut[typeNo]=40;
          return 5;
        case 16:
          strcpy(inputCIF->psp_name[typeNo],"16s.6.hgh");
          inputCIF->ecut[typeNo]=40;
          return 6;
        case 17:
          strcpy(inputCIF->psp_name[typeNo],"17cl.7.hgh");
          inputCIF->ecut[typeNo]=55;
          return 7;
        case 18:
          strcpy(inputCIF->psp_name[typeNo],"18ar.8.hgh");
          inputCIF->ecut[typeNo]=60;
          return 8;
        case 19:
          strcpy(inputCIF->psp_name[typeNo],"19k.9.hgh");
          inputCIF->ecut[typeNo]=70;
          return 9;
        case 20:
          strcpy(inputCIF->psp_name[typeNo],"20ca.10.hgh");
          inputCIF->ecut[typeNo]=75;
          return 10;
        case 21:
          strcpy(inputCIF->psp_name[typeNo],"21sc.11.hgh");
          inputCIF->ecut[typeNo]=110;
          return 11;
        case 22:
          strcpy(inputCIF->psp_name[typeNo],"22ti.12.hgh");
          inputCIF->ecut[typeNo]=135;
          return 12;
        case 23:
          strcpy(inputCIF->psp_name[typeNo],"23v.13.hgh");
          inputCIF->ecut[typeNo]=125;
          return 13;
        case 24:
          strcpy(inputCIF->psp_name[typeNo],"24cr.14.hgh");
          inputCIF->ecut[typeNo]=145;
          return 14;
        case 25:
          strcpy(inputCIF->psp_name[typeNo],"25mn.15.hgh");
          inputCIF->ecut[typeNo]=150;
          return 15;
        case 26:
          strcpy(inputCIF->psp_name[typeNo],"26fe.16.hgh");
          inputCIF->ecut[typeNo]=140;
          return 16;
        case 27:
          strcpy(inputCIF->psp_name[typeNo],"27co.17.hgh");
          inputCIF->ecut[typeNo]=110;
          return 17;
        case 28:
          strcpy(inputCIF->psp_name[typeNo],"28ni.18.hgh");
          inputCIF->ecut[typeNo]=165;
          return 18;
        case 29:
          strcpy(inputCIF->psp_name[typeNo],"29cu.11.hgh");
          inputCIF->ecut[typeNo]=105;
          return 11;
        case 30:
          strcpy(inputCIF->psp_name[typeNo],"30zn.12.hgh");
          inputCIF->ecut[typeNo]=110;
          return 12;
        case 31:
          strcpy(inputCIF->psp_name[typeNo],"31ga.13.hgh");
          inputCIF->ecut[typeNo]=155;
          return 13;
        case 32:
          strcpy(inputCIF->psp_name[typeNo],"32ge.4.hgh");
          inputCIF->ecut[typeNo]=25;
          return 4;
        case 33:
          strcpy(inputCIF->psp_name[typeNo],"33as.5.hgh");
          inputCIF->ecut[typeNo]=30;
          return 5;
        case 34:
          strcpy(inputCIF->psp_name[typeNo],"34se.6.hgh");
          inputCIF->ecut[typeNo]=35;
          return 6;
        case 35:
          strcpy(inputCIF->psp_name[typeNo],"35br.7.hgh");
          inputCIF->ecut[typeNo]=40;
          return 7;
        case 36:
          strcpy(inputCIF->psp_name[typeNo],"36kr.8.hgh");
          inputCIF->ecut[typeNo]=40;
          return 8;
        case 37:
          strcpy(inputCIF->psp_name[typeNo],"37rb.9.hgh");
          inputCIF->ecut[typeNo]=70;
          return 9;
        case 38:
          strcpy(inputCIF->psp_name[typeNo],"38sr.10.hgh");
          inputCIF->ecut[typeNo]=75;
          return 10;
        case 39:
          strcpy(inputCIF->psp_name[typeNo],"39y.11.hgh");
          inputCIF->ecut[typeNo]=40;
          return 11;
        case 40:
          strcpy(inputCIF->psp_name[typeNo],"40zr.12.hgh");
          inputCIF->ecut[typeNo]=45;
          return 12;
        case 41:
          strcpy(inputCIF->psp_name[typeNo],"41nb.13.hgh");
          inputCIF->ecut[typeNo]=45;
          return 13;
        case 42:
          strcpy(inputCIF->psp_name[typeNo],"42mo.14.hgh");
          inputCIF->ecut[typeNo]=55;
          return 14;
        case 43:
          strcpy(inputCIF->psp_name[typeNo],"43tc.15.hgh");
          inputCIF->ecut[typeNo]=120;
          return 7;
        case 44:
          strcpy(inputCIF->psp_name[typeNo],"44ru.16.hgh");
          inputCIF->ecut[typeNo]=55;
          return 16;
        case 45:
          strcpy(inputCIF->psp_name[typeNo],"45rh.17.hgh");
          inputCIF->ecut[typeNo]=60;
          return 17;
        case 46:
          strcpy(inputCIF->psp_name[typeNo],"46pd.18.hgh");
          inputCIF->ecut[typeNo]=60;
          return 18;
        case 47:
          strcpy(inputCIF->psp_name[typeNo],"47ag.11.hgh");
          inputCIF->ecut[typeNo]=45;
          return 11;
        case 48:
          strcpy(inputCIF->psp_name[typeNo],"48cd.12.hgh");
          inputCIF->ecut[typeNo]=50;
          return 12;
        case 49:
          strcpy(inputCIF->psp_name[typeNo],"49in.13.hgh");
          inputCIF->ecut[typeNo]=55;
          return 13;
        case 50:
          strcpy(inputCIF->psp_name[typeNo],"50sn.4.hgh");
          inputCIF->ecut[typeNo]=20;
          return 4;
        case 51:
          strcpy(inputCIF->psp_name[typeNo],"51sb.5.hgh");
          inputCIF->ecut[typeNo]=25;
          return 5;
        case 52:
          strcpy(inputCIF->psp_name[typeNo],"52te.6.hgh");
          inputCIF->ecut[typeNo]=25;
          return 6;
        case 53:
          strcpy(inputCIF->psp_name[typeNo],"53i.7.hgh");
          inputCIF->ecut[typeNo]=50;
          return 7;
        case 54:
          strcpy(inputCIF->psp_name[typeNo],"54xe.8.hgh");
          inputCIF->ecut[typeNo]=30;
          return 8;
        case 55:
          strcpy(inputCIF->psp_name[typeNo],"55cs.9.hgh");
          inputCIF->ecut[typeNo]=30;
          return 9;
        case 56:
          strcpy(inputCIF->psp_name[typeNo],"56ba.10.hgh");
          inputCIF->ecut[typeNo]=20;
          return 10;
        case 57:
          strcpy(inputCIF->psp_name[typeNo],"57la.11.hgh");
          return 11;
        case 58:
          strcpy(inputCIF->psp_name[typeNo],"58ce.12.hgh");
          return 12;
        case 59:
          strcpy(inputCIF->psp_name[typeNo],"59pr.13.hgh");
          return 13;
        case 60:
          strcpy(inputCIF->psp_name[typeNo],"60nd.14.hgh");
          return 14;
        case 61:
          strcpy(inputCIF->psp_name[typeNo],"61pm.15.hgh");
          return 15;
        case 62:
          strcpy(inputCIF->psp_name[typeNo],"62sm.16.hgh");
          return 16;
        case 63:
          strcpy(inputCIF->psp_name[typeNo],"63eu.17.hgh");
          return 17;
        case 64:
          strcpy(inputCIF->psp_name[typeNo],"64gd.18.hgh");
          return 18;
        case 65:
          strcpy(inputCIF->psp_name[typeNo],"65tb.19.hgh");
          return 19;
        case 66:
          strcpy(inputCIF->psp_name[typeNo],"66dy.20.hgh");
          return 20;
        case 67:
          strcpy(inputCIF->psp_name[typeNo],"67ho.14.hgh");
          return 21;
        case 68:
          strcpy(inputCIF->psp_name[typeNo],"68er.22.hgh");
          return 22;
        case 69:
          strcpy(inputCIF->psp_name[typeNo],"69tm.23.hgh");
          return 23;
        case 70:
          strcpy(inputCIF->psp_name[typeNo],"70yb.24.hgh");
          return 24;
        case 71:
          strcpy(inputCIF->psp_name[typeNo],"71lu.25.hgh");
          return 25;
        case 72:
          strcpy(inputCIF->psp_name[typeNo],"72hf.12.hgh");
          inputCIF->ecut[typeNo]=35;
          return 12;
        case 73:
          strcpy(inputCIF->psp_name[typeNo],"73ta.13.hgh");
          inputCIF->ecut[typeNo]=45;
          return 13;
        case 74:
          strcpy(inputCIF->psp_name[typeNo],"74w.14.hgh");
          inputCIF->ecut[typeNo]=50;
          return 14;
        case 75:
          strcpy(inputCIF->psp_name[typeNo],"75re.15.hgh");
          inputCIF->ecut[typeNo]=55;
          return 15;
        case 76:
          strcpy(inputCIF->psp_name[typeNo],"76os.16.hgh");
          return 16;
        case 77:
          strcpy(inputCIF->psp_name[typeNo],"77ir.17.hgh");
          inputCIF->ecut[typeNo]=65;
          return 17;
        case 78:
          strcpy(inputCIF->psp_name[typeNo],"78pt.18.hgh");
          inputCIF->ecut[typeNo]=65;
          return 18;
        case 79:
          strcpy(inputCIF->psp_name[typeNo],"79au.11.hgh");
          inputCIF->ecut[typeNo]=35;
          return 11;
        case 80:
          strcpy(inputCIF->psp_name[typeNo],"80hg.12.hgh");
          inputCIF->ecut[typeNo]=45;
          return 12;
        case 81:
          strcpy(inputCIF->psp_name[typeNo],"81tl.13.hgh");
          inputCIF->ecut[typeNo]=50;
          return 3;
        case 82:
          strcpy(inputCIF->psp_name[typeNo],"82pb.4.hgh");
          inputCIF->ecut[typeNo]=15;
          return 4;
        case 83:
          strcpy(inputCIF->psp_name[typeNo],"83bi.5.hgh");
          inputCIF->ecut[typeNo]=20;
          return 5;
        case 84:
          strcpy(inputCIF->psp_name[typeNo],"84po.6.hgh");
          return 6;
        case 85:
          strcpy(inputCIF->psp_name[typeNo],"85at.7.hgh");
          return 7;
      }
    }
  
}

int FindSpgorig(int spgroup){
  //Find if Spgorig for the space group will be relevant, if yes, return 3.
  switch (spgroup){
    case 48:
      return 3;
    case 50:
      return 3;
    case 59:
      return 3;
    case 70:
      return 3;
    case 85:
      return 3;
    case 86:
      return 3;
    case 88:
      return 3;
    case 125:
      return 3;
    case 126:
      return 3;
    case 129:
      return 3;
    case 130:
      return 3;
    case 133:
      return 3;
    case 134:
      return 3;
    case 137:
      return 3;
    case 141:
      return 3;
    case 142:
      return 3;
    case 201:
      return 3;
    case 203:
      return 3;
    case 222:
      return 3;
    case 224:
      return 3;
    case 227:
      return 3;
    case 228:
      return 3;
    default:
      return 1;
  }
}

int Findspgaxor(struct CIF * inputCIF){
  int spg_No=0;
  char HM_original[50];
  char HM_simplified[50];
  int ori_len=0;
  int i=0;
  int count=0;
  memset(HM_original,50,0);
  //printf("In new Findspgaxor!!!!!!!!!!!!!!!!!!!\n");
  spg_No=inputCIF->spgroup;
  //printf("spgNo is: %d \n", spg_No);
  if(spg_No==1 || spg_No==2 || spg_No>74){
    printf("Spgaxor irrelevant. \n");
    //printf("Space group > 3. Spgaxor irrelevant \n");
    return 0;
  }
  if(spg_No==16 || spg_No==19 || spg_No==22 || spg_No==23 || spg_No==24 || spg_No==47 || spg_No==48 || spg_No==69 || spg_No==70 || spg_No==71){
    printf("Spgaxor irrelevant \n");
    return 0;
  }
  printf("Spgaxor relevant! \n");
  spgaxor_relevant=1;
  strcpy(HM_original,inputCIF->HM_FullspgName_oneword);
  ori_len=strlen(HM_original);
  //printf("HM original is: %s, length is: %d \n",HM_original, ori_len);

  //Now get HM_simplified
  /*for(i=0;i<ori_len;i++){
    if(isalpha(HM_original[i])){
      HM_simplified[count]=HM_original[i];
      count++;
    }
    if(isalpha(HM_original[i])){
      HM_simplified[count]=HM_original[i];
      count++;
    }
  }*/
  
  //printf("Now begin to compare. \n");

  memset(inputCIF->HM_STDName,50,0);

  switch (spg_No){
     
    case 3:
    strcpy(inputCIF->HM_STDName,"P121");
    if(strcmp(HM_original,"P121")==0){
      printf("Match!\n");
      return 1;
    }else{
      printf("Mismatch!\n");
      return -1;
    }
    
    case 4:
    strcpy(inputCIF->HM_STDName,"P1211");
    if(strcmp(HM_original,"P1211")==0){
      printf("Match!\n");
      return 1;
    }else{
      printf("Mismatch!\n");
      return -1;
    }    

    case 5:
    strcpy(inputCIF->HM_STDName,"C121");
    if(strcmp(HM_original,"C121")==0){
      printf("Match!\n");
      return 1;
    }else{
      printf("Mismatch!\n");
      return -1;
    }

    case 6:
    strcpy(inputCIF->HM_STDName,"P1m1");
    if(strcmp(HM_original,"P1m1")==0){
      printf("Match!\n");
      return 1;
    }else{
      printf("Mismatch!\n");
      return -1;
    }

    case 7:
    strcpy(inputCIF->HM_STDName,"P1c1");
    if(strcmp(HM_original,"P1c1")==0){
      printf("Match!\n");
      return 1;
    }else{
      printf("Mismatch!\n");
      return -1;
    }

    case 8:
    strcpy(inputCIF->HM_STDName,"C1m1");
    if(strcmp(HM_original,"C1m1")==0){
      printf("Match!\n");
      return 1;
    }else{
      printf("Mismatch!\n");
      return -1;
    }

    case 9:
    strcpy(inputCIF->HM_STDName,"C1c1");
    if(strcmp(HM_original,"C1c1")==0){
      printf("Match!\n");
      return 1;
    }else{
      printf("Mismatch!\n");
      return -1;
    }

    case 10:
    strcpy(inputCIF->HM_STDName,"P12/m1");
    if(strcmp(HM_original,"P12/m1")==0){
      printf("Match!\n");
      return 1;
    }else{
      printf("Mismatch!\n");
      return -1;
    }

    case 11:
    strcpy(inputCIF->HM_STDName,"P121/m1");
    if(strcmp(HM_original,"P121/m1")==0){
      printf("Match!\n");
      return 1;
    }else{
      printf("Mismatch!\n");
      return -1;
    }

    case 12: 
    strcpy(inputCIF->HM_STDName,"C12/m1");
    if(strcmp(HM_original,"C12/m1")==0){
      printf("Match!\n");
      return 1;
    }else{
      printf("Mismatch!\n");
      return -1;
    }
    
    case 13:
    strcpy(inputCIF->HM_STDName,"P12/c1");
    if(strcmp(HM_original,"P12/c1")==0){
      printf("Match!\n");
      return 1;
    }else{
      printf("Mismatch!\n");
      return -1;
    }    

    case 14:
    strcpy(inputCIF->HM_STDName,"P121/c1");
    if(strcmp(HM_original,"P121/c1")==0){
      printf("Match!\n");
      return 1;
    }else{
      printf("Mismatch!\n");
      return -1;
    }

    case 15:
    strcpy(inputCIF->HM_STDName,"C12/c1");
    if(strcmp(HM_original,"C12/c1")==0){
      printf("Match!\n");
      return 1;
    }else{
      printf("Mismatch!\n");
      return -1;
    }    

    case 17:
    strcpy(inputCIF->HM_STDName,"P2221");
    if(strcmp(HM_original,"P2221")==0){
      printf("Match!\n");
      return 1;
    }else{
      printf("Mismatch!\n");
      return -1;
    }

    case 18:
    strcpy(inputCIF->HM_STDName,"P21212");
    if(strcmp(HM_original,"P21212")==0){
      printf("Match!\n");
      return 1;
    }else{
      printf("Mismatch!\n");
      return -1;
    }

    case 20:
    strcpy(inputCIF->HM_STDName,"C2221");
    if(strcmp(HM_original,"C2221")==0){
      printf("Match!\n");
      return 1;
    }else{
      printf("Mismatch!\n");
      return -1;
    }

    case 21:
    strcpy(inputCIF->HM_STDName,"C222");
    if(strcmp(HM_original,"C222")==0){
      printf("Match!\n");
      return 1;
    }else{
      printf("Mismatch!\n");
      return -1;
    }

    case 25:
    strcpy(inputCIF->HM_STDName,"Pmm2");
    if(strcmp(HM_original,"Pmm2")==0){
      printf("Match!\n");
      return 1;
    }else{
      printf("Mismatch!\n");
      return -1;
    }

    case 26:
    strcpy(inputCIF->HM_STDName,"Pmc21");
    if(strcmp(HM_original,"Pmc21")==0){
      printf("Match!\n");
      return 1;
    }else{
      printf("Mismatch!\n");
      return -1;
    }

    case 27:
    strcpy(inputCIF->HM_STDName,"Pcc2");
    if(strcmp(HM_original,"Pcc2")==0){
      printf("Match!\n");
      return 1;
    }else{
      printf("Mismatch!\n");
      return -1;
    }

    case 28:
    strcpy(inputCIF->HM_STDName,"Pma2");
    if(strcmp(HM_original,"Pma2")==0){
      printf("Match!\n");
      return 1;
    }else{
      printf("Mismatch!\n");
      return -1;
    }

    case 29:
    strcpy(inputCIF->HM_STDName,"Pca21");
    if(strcmp(HM_original,"Pca21")==0){
      printf("Match!\n");
      return 1;
    }else{
      printf("Mismatch!\n");
      return -1;
    }

    case 30:
    strcpy(inputCIF->HM_STDName,"Pnc2");
    if(strcmp(HM_original,"Pnc2")==0){
      printf("Match!\n");
      return 1;
    }else{
      printf("Mismatch!\n");
      return -1;
    }

    case 31:
    strcpy(inputCIF->HM_STDName,"Pmn21");
    if(strcmp(HM_original,"Pmn21")==0){
      printf("Match!\n");
      return 1;
    }else{
      printf("Mismatch!\n");
      return -1;
    }

    case 32:
    strcpy(inputCIF->HM_STDName,"Pba2");
    if(strcmp(HM_original,"Pba2")==0){
      printf("Match!\n");
      return 1;
    }else{
      printf("Mismatch!\n");
      return -1;
    }

    case 33:
    strcpy(inputCIF->HM_STDName,"Pna21");
    if(strcmp(HM_original,"Pna21")==0){
      printf("Match!\n");
      return 1;
    }else{
      printf("Mismatch!\n");
      return -1;
    }

    case 34:
    strcpy(inputCIF->HM_STDName,"Pnn2");
    if(strcmp(HM_original,"Pnn2")==0){
      printf("Match!\n");
      return 1;
    }else{
      printf("Mismatch!\n");
      return -1;
    }

    case 35:
    strcpy(inputCIF->HM_STDName,"Cmm2");
    if(strcmp(HM_original,"Cmm2")==0){
      printf("Match!\n");
      return 1;
    }else{
      printf("Mismatch!\n");
      return -1;
    }

    case 36:
    strcpy(inputCIF->HM_STDName,"Cmc21");
    if(strcmp(HM_original,"Cmc21")==0){
      printf("Match!\n");
      return 1;
    }else{
      printf("Mismatch!\n");
      return -1;
    }

    case 37:
    strcpy(inputCIF->HM_STDName,"Ccc2");
    if(strcmp(HM_original,"Ccc2")==0){
      printf("Match!\n");
      return 1;
    }else{
      printf("Mismatch!\n");
      return -1;
    }

    case 38:
    strcpy(inputCIF->HM_STDName,"Amm2");
    if(strcmp(HM_original,"Amm2")==0){
      printf("Match!\n");
      return 1;
    }else{
      printf("Mismatch!\n");
      return -1;
    }

    case 39:
    strcpy(inputCIF->HM_STDName,"Abm2");
    if(strcmp(HM_original,"Abm2")==0){
      printf("Match!\n");
      return 1;
    }else{
      printf("Mismatch!\n");
      return -1;
    }

    case 40:
    strcpy(inputCIF->HM_STDName,"Ama2");
    if(strcmp(HM_original,"Ama2")==0){
      printf("Match!\n");
      return 1;
    }else{
      printf("Mismatch!\n");
      return -1;
    }

    case 41:
    strcpy(inputCIF->HM_STDName,"Aba2");
    if(strcmp(HM_original,"Aba2")==0){
      printf("Match!\n");
      return 1;
    }else{
      printf("Mismatch!\n");
      return -1;
    }

    case 42:
    strcpy(inputCIF->HM_STDName,"Fmm2");
    if(strcmp(HM_original,"Fmm2")==0){
      printf("Match!\n");
      return 1;
    }else{
      printf("Mismatch!\n");
      return -1;
    }

    case 43:
    strcpy(inputCIF->HM_STDName,"Fdd2");
    if(strcmp(HM_original,"Fdd2")==0){
      printf("Match!\n");
      return 1;
    }else{
      printf("Mismatch!\n");
      return -1;
    }

    case 44:
    strcpy(inputCIF->HM_STDName,"Imm2");
    if(strcmp(HM_original,"Imm2")==0){
      printf("Match!\n");
      return 1;
    }else{
      printf("Mismatch!\n");
      return -1;
    }

    case 45:
    strcpy(inputCIF->HM_STDName,"Iba2");
    if(strcmp(HM_original,"Iba2")==0){
      printf("Match!\n");
      return 1;
    }else{
      printf("Mismatch!\n");
      return -1;
    }

    case 46:
    strcpy(inputCIF->HM_STDName,"Ima2");
    if(strcmp(HM_original,"Ima2")==0){
      printf("Match!\n");
      return 1;
    }else{
      printf("Mismatch!\n");
      return -1;
    }

    case 49:
    strcpy(inputCIF->HM_STDName,"Pccm");
    if(strcmp(HM_original,"Pccm")==0){
      printf("Match!\n");
      return 1;
    }else{
      printf("Mismatch!\n");
      return -1;
    }

    case 50:
    strcpy(inputCIF->HM_STDName,"Pban");
    if(strcmp(HM_original,"Pban")==0){
      printf("Match!\n");
      return 1;
    }else{
      printf("Mismatch!\n");
      return -1;
    }

    case 51:
    strcpy(inputCIF->HM_STDName,"Pmma");
    if(strcmp(HM_original,"Pmma")==0){
      printf("Match!\n");
      return 1;
    }else{
      printf("Mismatch!\n");
      return -1;
    }

    case 52:
    strcpy(inputCIF->HM_STDName,"Pnna");
    if(strcmp(HM_original,"Pnna")==0){
      printf("Match!\n");
      return 1;
    }else{
      printf("Mismatch!\n");
      return -1;
    }

    case 53:
    strcpy(inputCIF->HM_STDName,"Pmna");
    if(strcmp(HM_original,"Pmna")==0){
      printf("Match!\n");
      return 1;
    }else{
      printf("Mismatch!\n");
      return -1;
    }

    case 54:
    strcpy(inputCIF->HM_STDName,"Pcca");
    if(strcmp(HM_original,"Pcca")==0){
      printf("Match!\n");
      return 1;
    }else{
      printf("Mismatch!\n");
      return -1;
    }

    case 55:
    strcpy(inputCIF->HM_STDName,"Pbam");
    if(strcmp(HM_original,"Pbam")==0){
      printf("Match!\n");
      return 1;
    }else{
      printf("Mismatch!\n");
      return -1;
    }

    case 56:
    strcpy(inputCIF->HM_STDName,"Pccn");
    if(strcmp(HM_original,"Pccn")==0){
      printf("Match!\n");
      return 1;
    }else{
      printf("Mismatch!\n");
      return -1;
    }

    case 57:
    strcpy(inputCIF->HM_STDName,"Pbcm");
    if(strcmp(HM_original,"Pbcm")==0){
      printf("Match!\n");
      return 1;
    }else{
      printf("Mismatch!\n");
      return -1;
    }

    case 58:
    strcpy(inputCIF->HM_STDName,"Pnnm");
    if(strcmp(HM_original,"Pnnm")==0){
      printf("Match!\n");
      return 1;
    }else{
      printf("Mismatch!\n");
      return -1;
    }

    case 59:
    strcpy(inputCIF->HM_STDName,"Pmmn");
    if(strcmp(HM_original,"Pmmn")==0){
      printf("Match!\n");
      return 1;
    }else{
      printf("Mismatch!\n");
      return -1;
    }

    case 60:
    strcpy(inputCIF->HM_STDName,"Pbcn");
    if(strcmp(HM_original,"Pbcn")==0){
      printf("Match!\n");
      return 1;
    }else{
      printf("Mismatch!\n");
      return -1;
    }

    case 61:
    strcpy(inputCIF->HM_STDName,"Pbca");
    if(strcmp(HM_original,"Pbca")==0){
      printf("Match!\n");
      return 1;
    }else{
      printf("Mismatch!\n");
      return -1;
    }
    
    case 62:
    strcpy(inputCIF->HM_STDName,"Pnma");
    if(strcmp(HM_original,"Pnma")==0){
      printf("Match!\n");
      return 1;
    }else{
      printf("Mismatch!\n");
      return -1;
    }

    case 63:
    strcpy(inputCIF->HM_STDName,"Cmcm");
    if(strcmp(HM_original,"Cmcm")==0){
      printf("Match!\n");
      return 1;
    }else{
      printf("Mismatch!\n");
      return -1;
    }

    case 64:
    strcpy(inputCIF->HM_STDName,"Cmca");
    if(strcmp(HM_original,"Cmca")==0){
      printf("Match!\n");
      return 1;
    }else{
      printf("Mismatch!\n");
      return -1;
    }

    case 65:
    strcpy(inputCIF->HM_STDName,"Cmmm");
    if(strcmp(HM_original,"Cmmm")==0){
      printf("Match!\n");
      return 1;
    }else{
      printf("Mismatch!\n");
      return -1;
    }

    case 66:
    strcpy(inputCIF->HM_STDName,"Cccm");
    if(strcmp(HM_original,"Cccm")==0){
      printf("Match!\n");
      return 1;
    }else{
      printf("Mismatch!\n");
      return -1;
    }

    case 67:
    strcpy(inputCIF->HM_STDName,"Cmma");
    if(strcmp(HM_original,"Cmma")==0){
      printf("Match!\n");
      return 1;
    }else{
      printf("Mismatch!\n");
      return -1;
    }

    case 68:
    strcpy(inputCIF->HM_STDName,"Ccca");
    if(strcmp(HM_original,"Ccca")==0){
      printf("Match!\n");
      return 1;
    }else{
      printf("Mismatch!\n");
      return -1;
    }

    case 72:
    strcpy(inputCIF->HM_STDName,"Ibam");
    if(strcmp(HM_original,"Ibam")==0){
      printf("Match!\n");
      return 1;
    }else{
      printf("Mismatch!\n");
      return -1;
    }

    case 73:
    strcpy(inputCIF->HM_STDName,"Ibca");
    if(strcmp(HM_original,"Ibca")==0){
      printf("Match!\n");
      return 1;
    }else{
      printf("Mismatch!\n");
      return -1;
    }

    case 74:
    strcpy(inputCIF->HM_STDName,"Imma");
    if(strcmp(HM_original,"Imma")==0){
      printf("Match!\n");
      return 1;
    }else{
      printf("Mismatch!\n");
      return -1;
    }

    default:
    return 9;
  }
}

int Adjust_nband(int input_nband){
  double input, temp;
  int output;
  input = (double) input_nband;
  temp = input/12.0;
  output = ceil(temp)*12;
  return output;
}

int find_spg_type(int spgNO){
  if(spgNO<3){
    printf("Triclinic space group. \n");
    return 1;
  }else if(spgNO<16){
    printf("Monoclinic space group. \n");
    return 2;
  }else if(spgNO<75){
    printf("Orthorhombic space group. \n");
    return 3;
  }else if(spgNO<143){
    printf("Tetragonal space group. \n");
    return 4;
  }else if(spgNO<168){
    printf("Trigonal space group. \n");
    return 5;
  }else if(spgNO<195){
    printf("Hexagonal space group. \n");
    return 6;
  }else{
    printf("Cubic space group. \n");
    return 7;
  }
}

void ProcessCIFdata(struct CIF * inputCIF, int mode){
    int i=0, j=0;
    char * temp_receiver;
    int nband=0, current_e_nband=0, atom_nband=0;
    double temp_number;
    int ecut_max=3;
    //Before: total site number in CIF and spg group number read.
    //Before: abinit input: spgroup, natrd, ntypat.

    //Initialize the natom_type matrix to all 0s.
    inputCIF->spg_type_No=find_spg_type(inputCIF->spgroup);
    //printf("The space group type number is: %d. \n",inputCIF->spg_type_No);
    for(i=0;i<10;i++){
      inputCIF->natom_type[i]=0;
      inputCIF->No_site_by_type[i]=0;
    }
 
    //This first part deals with cell a, b, c and ang alpha, beta, gamma
    //Generate abinit input: acell, angdeg
    inputCIF->cell_a=NumberString(inputCIF->cell_a_str);
    //printf("Cell parameter a after processing: %lf \n",inputCIF->cell_a);
    inputCIF->cell_b=NumberString(inputCIF->cell_b_str);
    //printf("Cell parameter b after processing: %lf \n",inputCIF->cell_b);
    inputCIF->cell_c=NumberString(inputCIF->cell_c_str);
    //printf("Cell parameter c after processing: %lf \n",inputCIF->cell_c);
    inputCIF->ang_alpha=NumberString(inputCIF->ang_alpha_str);
    //printf("Angle alpha after : %lf \n",inputCIF->ang_alpha);
    inputCIF->ang_beta=NumberString(inputCIF->ang_beta_str);
    //printf("Angle beta after : %lf \n",inputCIF->ang_beta);
    inputCIF->ang_gamma=NumberString(inputCIF->ang_gamma_str);
    //printf("Angle gamma after : %lf \n",inputCIF->ang_gamma);



    if(inputCIF->spg_type_No==4){
      inputCIF->cell_b=inputCIF->cell_a;
    }else if(inputCIF->spg_type_No==5){
      if(inputCIF->ang_gamma>=119.9 && inputCIF->ang_gamma<=120.1){
        printf("Trigonal: gamma = 120 degree. \n");
        inputCIF->cell_b=inputCIF->cell_a;
        spgaxor_relevant=2;
        inputCIF->spgaxor=1;
      }else{
        printf("Trigonal: alpha = beta = gamma. \n");
        inputCIF->cell_b=inputCIF->cell_a;
        inputCIF->cell_c=inputCIF->cell_a;
        inputCIF->spgaxor=2;
        spgaxor_relevant=2;
      }
    }else if(inputCIF->spg_type_No==6){
      inputCIF->cell_b=inputCIF->cell_a;
    }else if(inputCIF->spg_type_No==7){
      inputCIF->cell_b=inputCIF->cell_a;
      inputCIF->cell_c=inputCIF->cell_a;
    }

    //Now generating spgaxor
    if(spgaxor_relevant==0){
      //printf("spg_type_No: %d, spgroup: %d, HM_name: %s \n", inputCIF->spg_type_No,inputCIF->spgroup,inputCIF->HM_FullspgName_oneword);
      inputCIF->spgaxor=Findspgaxor(inputCIF);
      //printf("Out of find spgaxor, value is: %d \n", inputCIF->spgaxor);
    }

    if(spgaxor_relevant==1){
      printf("The standard H-M designation is: %s \n", inputCIF->HM_STDName);
    }

    //This second part deals with element identity and pseudopotential
    //Generate abinit input: znucl, files file
    i=1;
    while(i<=ntype){
      //printf("Before extract letter: input %d: %s. \n",(i+1),inputCIF->chemical_formula[i-1]);
      temp_receiver=ExtractLetter(inputCIF->chemical_formula[i-1]);
      //printf("After extract letter: input %d: %s. \n",(i+1),temp_receiver);
      strcpy(inputCIF->element_identity[i-1],temp_receiver);
      i++;
    }
    i=1;
    while(i<=ntype){
      inputCIF->element_No[i-1]=AtomNoFinder(inputCIF->element_identity[i-1]);
      if(inputCIF->element_No[i-1]==666){
        printf("One of the element in the CIF file does not match with any elements in the periodic table!\n");
        printf("Double check the CIF file!!\n");
      }
      //printf("Element %d: %s, Number %d. \n",i,inputCIF->element_identity[i-1],inputCIF->element_No[i-1]);
      i++;
    }
    i=1;
    while(i<=ntype){
      inputCIF->nelectrons_psp[i-1]=PSP_Assigner(inputCIF->element_No[i-1],inputCIF,i-1); 
      //printf("Element: %s, %d electrons, %s, ecut: %d. \n",inputCIF->element_identity[i-1],inputCIF->nelectrons_psp[i-1],inputCIF->psp_name[i-1],inputCIF->ecut[i-1]);
      i++;
    }
    
    i=0;
    inputCIF->natom=0;
    
    //This third part deals with different sites
    //Generate abinit input: xred, typat, natom
    for(i=0;i<inputCIF->total_site;i++){
      inputCIF->xred[i]=NumberString(inputCIF->xred_str[i]);
      sprintf(inputCIF->xred_upd[i],"%lf",inputCIF->xred[i]);
      inputCIF->yred[i]=NumberString(inputCIF->yred_str[i]);
      sprintf(inputCIF->yred_upd[i],"%lf",inputCIF->yred[i]);
      inputCIF->zred[i]=NumberString(inputCIF->zred_str[i]);
      sprintf(inputCIF->zred_upd[i],"%lf",inputCIF->zred[i]);
      inputCIF->occupancy[i]=NumberString(inputCIF->occupancy_str[i]);
      //printf("The reduced coordinates for site %d is: %lf %lf %lf \n", (i+1), inputCIF->xred[i], inputCIF->yred[i], inputCIF->zred[i]);
      inputCIF->natom=inputCIF->natom+inputCIF->site_eqv[i];
      temp_receiver=ExtractLetter(inputCIF->site_identity[i]);
      strcpy(inputCIF->site_identity_processed[i],temp_receiver);
      inputCIF->site_element_No[i]=AtomNoFinder(inputCIF->site_identity_processed[i]);
      ////printf("The identity of site %d is: %d %s \n", (i+1), inputCIF->site_element_No[i], inputCIF->site_identity_processed[i]);
      for(j=0;j<ntype;j++){
        if(inputCIF->site_element_No[i]==inputCIF->element_No[j]){
          inputCIF->natom_type[j]=inputCIF->natom_type[j]+inputCIF->site_eqv[i];
          inputCIF->No_site_by_type[j]++;
          //printf("Site %d, identity: %s, duplicity: %d belongs to element %d. \n",(i+1),inputCIF->site_identity_processed[i],inputCIF->site_eqv[i],inputCIF->element_No[j]);
        }
      }
    }
    for(i=0;i<inputCIF->total_site;i++){
      for(j=0;j<ntype;j++){
        if(inputCIF->site_element_No[i]==inputCIF->element_No[j]){
          inputCIF->site_typat[i]=j+1;
        }
      }
    }
    //printf("Typat: ");
    for(i=0;i<inputCIF->total_site;i++){
      //printf("%d ",inputCIF->site_typat[i]);
    }
    //printf("\n");
    //printf("Natom is: %d \n",inputCIF->natom); 
    //Now find site element number and number of atom for each type 
    for(j=0;j<ntype;j++){ 
      //printf("There are %d %s atoms in %d differnt sites in one unit cell. \n", inputCIF->natom_type[j],inputCIF->element_identity[j],inputCIF->No_site_by_type[j]);
    }
    i=0;
    j=0;

    //This part deals with other needed input variables
    //Generate abinit input: spgorig, brvltt, nband
    inputCIF->spg_center_letter=inputCIF->spgroup_center[1];
    //printf("Space group centering: %c \n", inputCIF->spg_center_letter);
    inputCIF->brvltt=0; 
    if( (inputCIF->spg_center_letter=='P') || (inputCIF->spg_center_letter=='p')){
      inputCIF->brvltt=0;
    }
    else if( (inputCIF->spg_center_letter=='F') || (inputCIF->spg_center_letter=='f')) { 
      inputCIF->brvltt=-2;
    }
    else if( (inputCIF->spg_center_letter=='R') || (inputCIF->spg_center_letter=='r')) {
      inputCIF->brvltt=3;  //Can be changed later if figured out how does the Rhmohedral centering owrk
    }
    else{
      inputCIF->brvltt=-1;
    }
    //printf("brvltt value: %d \n", inputCIF->brvltt);
    //Now find spgorig
    inputCIF->spgorig=FindSpgorig(inputCIF->spgroup);
    if(inputCIF->spgorig==3){
      printf("spgorig relevant. \n");
    }
    if(inputCIF->spgorig==1){
      printf("spgorig irrelevant. \n");
    }

    //Now find spgaxor

    //Now find nband
    for(i=0;i<ntype;i++){
      current_e_nband=inputCIF->natom_type[i]*inputCIF->nelectrons_psp[i];
      nband=nband+current_e_nband; 
    }
    if(SOC_flag!=1){
      nband=nband*2;
    }
    nband=ceil(nband/2);
    atom_nband=inputCIF->natom*2;
    nband=nband+atom_nband;
    temp_number=(double) nband;
    if(convert_primitive_flag==1){
      if(inputCIF->brvltt==-1){
        temp_number=(double) nband/2;
        inputCIF->natom=inputCIF->natom/2; 
      }
      else if(inputCIF->brvltt==-2){
        temp_number=(double) nband/4;
        inputCIF->natom=inputCIF->natom/4;
      }else if(inputCIF->brvltt==-3){ //NOT useful at this moment
        temp_number=(double) nband/3;
        inputCIF->natom=inputCIF->natom/3;
      }
    }
    nband=Adjust_nband(temp_number);
    inputCIF->nband=nband;
    //printf("Value of nband: %d \n",inputCIF->nband);
    
    //Now find ecut
    for(i=0;i<ntype;i++){
      ecut_unknown[i]=0;
    }
    for(i=0;i<ntype;i++){
      if(ecut_max<inputCIF->ecut[i]) ecut_max=inputCIF->ecut[i];
      if(inputCIF->ecut[i]==1) ecut_unknown[i]=1;
    }
    inputCIF->ecut_final=ecut_max;
    //printf("Ecut value is: %d. \n",inputCIF->ecut_final);    
    inputCIF->ecut_lambda_bohr=2*PI/sqrt(2*ecut_max); 
    inputCIF->ecut_lambda_A=inputCIF->ecut_lambda_bohr*0.529177;

    //Set old checkprim to 0
    inputCIF->old_chkprim=0;

    if(inputCIF->spgaxor==-1 && convert_primitive_flag==1 && kmesh_flag==1){
      if((inputCIF->spg_center_letter != 'R') && (inputCIF->spg_center_letter !='P')){
        printf("Your cell is centered and you chose to convert to primitive. \n");
        printf("You chose to use ngkpt to present your K-point mesh");
        printf("Your space group %d has standard setting %s according to the International Tables for Crystallography. \n",inputCIF->spgroup,inputCIF->HM_STDName);
        printf("However, the setting your CIF has is %s. \n", inputCIF->HM_FullspgName_oneword);
        printf("This means for ABINIT to carry out your calculation, first a permutation of axes is required. Then, reduced unit cell must be calculated to generate ngkpt. \n");
        printf("Such changes cannot be calculated by prepareABINIT because ABINIT used different permutations for different space groups. \n");
        printf("As a result, no .in and .files will be generated. \n\n");
        printf("To do: There are many choices you have here: \n");
        printf("1. Use kptrlatt to specify K-mesh and let ABINIT generate K-mesh for you.\n");
        printf("2. Load your CIF file into a software like Diamond and out put it in the %s standard setting.\n",inputCIF->HM_STDName);
        printf("3. Do not convert to primitive cell.\n\n");
        printf("NOW EXIT prepareABINIT. NO FILES ARE GENERATED. CHANGE OPTIONS OR CIF FILE AND THEN RUN prepareABINIT AGAIN!!! \n");
        exit(0);
      }
    }
}

void Write_files_file(struct CIF * inputCIF, FILE * files_file, char * psp_address_header, char * outputbase){
    int i=0;
    printf("Now writing the ABINIT .files file: %s.files ...", outputbase);
    fprintf(files_file,"%s.in\n",outputbase);
    fprintf(files_file,"%s.out\n",outputbase);
    fprintf(files_file,"%s_i\n",outputbase);
    fprintf(files_file,"%s_o\n",outputbase);
    fprintf(files_file,"\n");
    for(i=0;i<ntype;i++){
      fprintf(files_file,"%s%s",psp_address_header,inputCIF->psp_name[i]);
      if (i != ntype-1) fprintf(files_file,"\n");
    }
    fprintf(files_file,"\n",outputbase);
    printf("Done. \n");
}

void Write_in_file(struct CIF * inputCIF, FILE * in_file, char * outputbase, char * CIFname, int mode){
    //Declare variables
    int i=0;
    int Warning=0;
    int Current_Warning=0;
    FILE * LOGfile;
    char LOGfilename[100];
    //Print header
    printf("Now writing the ABINIT .in file: %s.in ... ", outputbase);
    fprintf(in_file,"!----Abinit input file generated by prepareAbinit----!\n");
    fprintf(in_file,"!----Calculation based on %s----!\n\n",CIFname);

    //Print default options section
    fprintf(in_file,"!----DEFAULT OPTIONS----!\n");
    fprintf(in_file,"     nstep  100         ! max SCF steps\n");
    fprintf(in_file,"     toldfe  1.0d-5 eV  ! energy difference for stop SCF steps\n");
    fprintf(in_file,"     occopt  3          ! metallic smearing of bands\n");
    fprintf(in_file,"     tsmear  0.005      ! smearing factor (temperature)\n");
    fprintf(in_file,"     prtvol  3          ! output file printing option\n");
    if(mode==4 || mode==5 || mode==6 ){
      fprintf(in_file,"     prtwf   1          ! 0 = DO NOT write WFK; 1=write WFK\n");
    }else{
      fprintf(in_file,"     prtwf   0          ! 0 = DO NOT write WFK; 1=write WFK\n");
    }
    if(mode==4 || mode==5 || mode==6 ){
      fprintf(in_file,"     istwfk  %d*1       ! Used for generate WFK for nonlocal\n\n",inputCIF->nkpt);
    }else{
      fprintf(in_file,"!    istwfk  (NKPT)*1   ! Used for generate WFK for nonlocal\n\n");
    }

    //Print chemical system section
    fprintf(in_file,"!----CHEMICAL SYSTEM----!\n");
    fprintf(in_file,"     ntypat  %d         ! number of different element types\n", ntype);
    fprintf(in_file,"     znucl ");
    for(i=0;i<ntype;i++){
      fprintf(in_file," %d",inputCIF->element_No[i]);
    }
    fprintf(in_file,"      ! element number of each element.\n");
    fprintf(in_file,"     natom %d         ! number of atom in unit cell\n", inputCIF->natom);
    fprintf(in_file,"     nband %d         ! number of bands\n\n", inputCIF->nband); 
 
    //Print starting geometry section
    fprintf(in_file,"!----STARTING GEOMETRY----!\n");
    if(mode==1 || mode==2){
      fprintf(in_file,"     acell  %lf %lf %lf angstrom     ! conventional cell parameters\n",inputCIF->cell_a,inputCIF->cell_b,inputCIF->cell_c);
    }else{
      fprintf(in_file,"     acell  %lf %lf %lf bohr         ! conventional cell parameters\n",inputCIF->cell_a,inputCIF->cell_b,inputCIF->cell_c);
    }
    if(rprim_found==0){
      fprintf(in_file,"     angdeg  %lf %lf %lf          ! conventional cell angles\n",inputCIF->ang_alpha,inputCIF->ang_beta,inputCIF->ang_gamma);
    }else{
      fprintf(in_file,"     rprim   %s %s %s          ! rprim from abinit out file\n",inputCIF->rprim_str[0][0],inputCIF->rprim_str[0][1],inputCIF->rprim_str[0][2]);
      fprintf(in_file,"             %s %s %s          \n",inputCIF->rprim_str[1][0],inputCIF->rprim_str[1][1],inputCIF->rprim_str[1][2]);
      fprintf(in_file,"             %s %s %s          \n",inputCIF->rprim_str[2][0],inputCIF->rprim_str[2][1],inputCIF->rprim_str[2][2]);
    }

    if(inputCIF->brvltt == 3){
      fprintf(in_file,"     chkprim  0            ! do not check for primitive cell\n");
    }else if(inputCIF->brvltt != 0){
      if(convert_primitive_flag==1){
        fprintf(in_file,"     brvltt  -1            ! convert to primitive cell\n");
      }
      else{
        fprintf(in_file,"     chkprim  0            ! do not check for primitive cell\n");
      }
    }      
    if(inputCIF->spgroup != 0){
      fprintf(in_file,"     spgroup  %d          ! IUC space group number\n",inputCIF->spgroup);
    }
    else{
      fprintf(in_file,"!    spgroup  %d          ! IUC space group number\n",inputCIF->spgroup);
    }
    //printf("\nvalue of old_brvltt: %d\n\n", inputCIF->old_brvltt);
    
    if(inputCIF->old_chkprim==1){
      fprintf(in_file,"     chkprim  0            ! do not check for primitive cell\n");
    }
    
    if(inputCIF->old_brvltt==0){
      //printf("\nold brvltt is 0\n\n");
      if(inputCIF->spgorig==1 || inputCIF->spgorig==2){
        fprintf(in_file,"     spgorig  %d            ! Space group setting\n",inputCIF->spgorig);
      }
      else if(inputCIF->spgorig==3){
        fprintf(in_file,"     spgorig  1 OR 2        ! Space group setting\n");
      }
    }


    if(inputCIF->old_brvltt==0){
      if(spgaxor_relevant==2){
        if(inputCIF->spgaxor==2){
          fprintf(in_file,"   spgaxor  2             ! Trigonal space group - Rhomohedral axes \n",inputCIF->spgaxor);
        }
      }else if(spgaxor_relevant==1){
        if(inputCIF->spgaxor==-1){
          fprintf(in_file,"   spgaxor  ?             ! Orthorhombic / Monoclinic non-standard setting. Change manually\n",inputCIF->spgaxor);
        }
      }else if(spgaxor_relevant==5){
        fprintf(in_file,"     spgaxor  %d            ! Spgaxor from old input file\n",inputCIF->spgaxor);
      }
    }

    fprintf(in_file,"     natrd  %d           ! Number of atoms listed below - total sites\n",inputCIF->total_site);
    fprintf(in_file,"     typat  ");
    for(i=0;i<inputCIF->total_site;i++){
      fprintf(in_file,"%d ",inputCIF->site_typat[i]);
    }
    fprintf(in_file,"    ! identity of each atom listed\n",inputCIF->site_typat[i]);
    fprintf(in_file,"     xred               ! reduced coordinates of each atom listed\n");
    if(mode==1||mode==2||mode==9){
      for(i=0;i<inputCIF->total_site;i++){
        fprintf(in_file,"     %.10lf %.10lf %.10lf \n",inputCIF->xred[i],inputCIF->yred[i],inputCIF->zred[i]);
      }
    }else if(mode>2&&mode!=9){
      for(i=0;i<inputCIF->total_site;i++){
        fprintf(in_file,"     %s %s %s \n",inputCIF->xred_upd[i],inputCIF->yred_upd[i],inputCIF->zred_upd[i]);
      }
    }
    fprintf(in_file,"\n");
    
    //Print mesh of kpoints section
    if(mode==1&&kmesh_flag==1){
      fprintf(in_file,"!----MESH OF KPOINTS----!\n");
      fprintf(in_file,"     ndtset   %d            ! number of ngkpt to be included \n",inputCIF->kpt_mesh_accepted);
      fprintf(in_file,"     ngkpt1   %d %d %d        ! distance betweeen Kpoints<0.3 A-1\n",inputCIF->ngkpt[0][0],inputCIF->ngkpt[0][1],inputCIF->ngkpt[0][2]);
      if(inputCIF->kpt_mesh_accepted>1){
        fprintf(in_file,"     ngkpt2   %d %d %d        ! distance betweeen Kpoints<0.2 A-1\n",inputCIF->ngkpt[1][0],inputCIF->ngkpt[1][1],inputCIF->ngkpt[1][2]);
      }else{
        fprintf(in_file,"!    ngkpt2   %d %d %d        ! distance betweeen Kpoints<0.2 A-1, not included by default due to too many kpts\n",inputCIF->ngkpt[1][0],inputCIF->ngkpt[1][1],inputCIF->ngkpt[1][2]);
      }
      if(inputCIF->kpt_mesh_accepted>2){
        fprintf(in_file,"     ngkpt3   %d %d %d        ! distance betweeen Kpoints<0.1 A-1\n",inputCIF->ngkpt[2][0],inputCIF->ngkpt[2][1],inputCIF->ngkpt[2][2]);
      }else{
        fprintf(in_file,"!    ngkpt3   %d %d %d        ! distance betweeen Kpoints<0.1 A-1, not included by default due to too many kpts\n",inputCIF->ngkpt[2][0],inputCIF->ngkpt[2][1],inputCIF->ngkpt[2][2]);
      }       
      if(inputCIF->kpt_mesh_accepted>3){
        fprintf(in_file,"     ngkpt4   %d %d %d        ! distance betweeen Kpoints<0.07 A-1\n",inputCIF->ngkpt[3][0],inputCIF->ngkpt[3][1],inputCIF->ngkpt[3][2]);
      }else{
        fprintf(in_file,"!    ngkpt4   %d %d %d        ! distance betweeen Kpoints<0.07 A-1, not included by default due to too many kpts\n",inputCIF->ngkpt[3][0],inputCIF->ngkpt[3][1],inputCIF->ngkpt[3][2]);
      }
      if(inputCIF->kpt_mesh_accepted>4){
        fprintf(in_file,"     ngkpt5   %d %d %d        ! distance betweeen Kpoints<0.05 A-1\n",inputCIF->ngkpt[4][0],inputCIF->ngkpt[4][1],inputCIF->ngkpt[4][2]);
      }else{
        fprintf(in_file,"!    ngkpt5   %d %d %d        ! distance betweeen Kpoints<0.05 A-1, not included by default due to too many kpts\n",inputCIF->ngkpt[4][0],inputCIF->ngkpt[4][1],inputCIF->ngkpt[4][2]);
      }
      fprintf(in_file,"     shiftk  0.0 0.0 0.0   ! Shift of kpoints\n");
      fprintf(in_file,"!     ngkpt   %d %d %d         ! Recommended K-point mesh, K-point distance<0.1 A-1 \n\n",inputCIF->ngkpt[2][0],inputCIF->ngkpt[2][1],inputCIF->ngkpt[2][2]);
    }else if(mode==2&&kmesh_flag==1){
      fprintf(in_file,"!----MESH OF KPOINTS----!\n");
      fprintf(in_file,"!     ndtset   5            ! number of datasets \n");
      fprintf(in_file,"!     ngkpt1   %d %d %d        ! distance betweeen Kpoints<0.3 A-1\n",inputCIF->ngkpt[0][0],inputCIF->ngkpt[0][1],inputCIF->ngkpt[0][2]);
      fprintf(in_file,"!     ngkpt2   %d %d %d        ! distance betweeen Kpoints<0.2 A-1\n",inputCIF->ngkpt[1][0],inputCIF->ngkpt[1][1],inputCIF->ngkpt[1][2]);
      fprintf(in_file,"!     ngkpt3   %d %d %d        ! distance betweeen Kpoints<0.1 A-1\n",inputCIF->ngkpt[2][0],inputCIF->ngkpt[2][1],inputCIF->ngkpt[2][2]);
      fprintf(in_file,"!     ngkpt4   %d %d %d        ! distance betweeen Kpoints<0.07 A-1\n",inputCIF->ngkpt[3][0],inputCIF->ngkpt[3][1],inputCIF->ngkpt[3][2]);
      fprintf(in_file,"!     ngkpt5   %d %d %d        ! distance betweeen Kpoints<0.05 A-1\n",inputCIF->ngkpt[4][0],inputCIF->ngkpt[4][1],inputCIF->ngkpt[4][2]);
      fprintf(in_file,"     ngkpt   %d %d %d         ! Recommended K-point mesh, K-point distance<0.1 A-1 \n",inputCIF->ngkpt[2][0],inputCIF->ngkpt[2][1],inputCIF->ngkpt[2][2]);
      fprintf(in_file,"     shiftk  0.0 0.0 0.0   ! Shift of kpoints\n\n");
    }else if(mode==9&&kmesh_flag==1){
      fprintf(in_file,"!----MESH OF KPOINTS----!\n");
      fprintf(in_file,"!     ndtset   5            ! number of datasets \n");
      fprintf(in_file,"!     ngkpt1   %d %d %d        ! distance betweeen Kpoints<0.3 A-1\n",inputCIF->ngkpt[0][0],inputCIF->ngkpt[0][1],inputCIF->ngkpt[0][2]);
      fprintf(in_file,"!     ngkpt2   %d %d %d        ! distance betweeen Kpoints<0.2 A-1\n",inputCIF->ngkpt[1][0],inputCIF->ngkpt[1][1],inputCIF->ngkpt[1][2]);
      fprintf(in_file,"!     ngkpt3   %d %d %d        ! distance betweeen Kpoints<0.1 A-1\n",inputCIF->ngkpt[2][0],inputCIF->ngkpt[2][1],inputCIF->ngkpt[2][2]);
      fprintf(in_file,"!     ngkpt4   %d %d %d        ! distance betweeen Kpoints<0.07 A-1\n",inputCIF->ngkpt[3][0],inputCIF->ngkpt[3][1],inputCIF->ngkpt[3][2]);
      fprintf(in_file,"!     ngkpt5   %d %d %d        ! distance betweeen Kpoints<0.05 A-1\n",inputCIF->ngkpt[4][0],inputCIF->ngkpt[4][1],inputCIF->ngkpt[4][2]);
      fprintf(in_file,"     ngkpt   %d %d %d         ! Recommended K-point mesh, K-point distance<0.1 A-1 \n",inputCIF->ngkpt[2][0],inputCIF->ngkpt[2][1],inputCIF->ngkpt[2][2]);
      fprintf(in_file,"     shiftk  0.0 0.0 0.0   ! Shift of kpoints\n\n");
    }else if(mode==1&&kmesh_flag==2){
      fprintf(in_file,"!----MESH OF KPOINTS----!\n");
      fprintf(in_file,"     prtkpt  1        ! generate recommended kpt mesh\n");
      fprintf(in_file,"     kptrlen  120     ! kpt mesh quality - larger value = better\n");
      fprintf(in_file,"!     kptrlatt 0 0 0   ! an alternative to ngkpt, generated by ABINIT\n");
      fprintf(in_file,"!              0 0 0  \n");
      fprintf(in_file,"!              0 0 0  \n");
      fprintf(in_file,"!     shiftk  0.0 0.0 0.0   ! Shift of kpoints, generated by ABINIT\n\n");
    }else if(mode==2&&kmesh_flag==2){
      fprintf(in_file,"!----MESH OF KPOINTS----!\n");
      fprintf(in_file,"!     prtkpt  1        ! generate recommended kpt mesh\n");
      fprintf(in_file,"!     kptrlen  120     ! kpt mesh quality - larger value = better\n");
      fprintf(in_file,"     kptrlatt  _ _ _   ! an alternative to ngkpt, generated by ABINIT\n");
      fprintf(in_file,"               _ _ _  \n");
      fprintf(in_file,"               _ _ _  \n");
      fprintf(in_file,"     shiftk    _ _ _   ! Shift of kpoints, generated by ABINIT\n\n");
    }else if(mode==9&&kmesh_flag==2){
      fprintf(in_file,"!----MESH OF KPOINTS----!\n");
      fprintf(in_file,"!     prtkpt  1        ! generate recommended kpt mesh\n");
      fprintf(in_file,"!     kptrlen  120     ! kpt mesh quality - larger value = better\n");
      fprintf(in_file,"     kptrlatt  _ _ _   ! an alternative to ngkpt, generated by ABINIT\n");
      fprintf(in_file,"               _ _ _  \n");
      fprintf(in_file,"               _ _ _  \n");
      fprintf(in_file,"     shiftk    _ _ _   ! Shift of kpoints, generated by ABINIT\n\n");
    }else if(mode>2&&kmesh_flag==1){
      fprintf(in_file,"!----MESH OF KPOINTS----!\n");
      fprintf(in_file,"     ngkpt %d %d %d     ! ngkpt based on old input file\n",inputCIF->ngkpt_from_in[0],inputCIF->ngkpt_from_in[1],inputCIF->ngkpt_from_in[2]);
      fprintf(in_file,"     shiftk  %f %f %f   ! Shift of kpoints based on old input file\n\n",inputCIF->shiftk[0],inputCIF->shiftk[1],inputCIF->shiftk[2]);
    }else if(mode>2&&kmesh_flag==2){
      fprintf(in_file,"!----MESH OF KPOINTS----!\n");
      fprintf(in_file,"!     prtkpt  1        ! generate recommended kpt mesh\n");
      fprintf(in_file,"!     kptrlen  120     ! kpt mesh quality - larger value = better\n");
      fprintf(in_file,"     kptrlatt %d %d %d   ! kpoint mesh based on output file\n",inputCIF->kptrlatt[0],inputCIF->kptrlatt[1],inputCIF->kptrlatt[2]);
      fprintf(in_file,"              %d %d %d  \n",inputCIF->kptrlatt[3],inputCIF->kptrlatt[4],inputCIF->kptrlatt[5]);
      fprintf(in_file,"              %d %d %d  \n",inputCIF->kptrlatt[6],inputCIF->kptrlatt[7],inputCIF->kptrlatt[8]);
      fprintf(in_file,"     shiftk  %f %f %f   ! Shift of kpoints based on old input file\n\n",inputCIF->shiftk[0],inputCIF->shiftk[1],inputCIF->shiftk[2]);
    }
 
    //Print planwave ecut section
    if(mode==9){
      fprintf(in_file,"!----PLANEWAVE AND ECUT----!\n");
      fprintf(in_file,"     ndtset 15        ! number of ecut value to test\n");
      fprintf(in_file,"     ecut: 30         ! starting ecut value\n");
      fprintf(in_file,"     ecut+ 10         ! ecut increment per test\n");
      fprintf(in_file,"     getden -1        ! re-use electron density\n\n");
    }else{
      fprintf(in_file,"!----PLANEWAVE AND ECUT----!\n");
      fprintf(in_file,"!    ndtset 0        ! number of ecut value to test\n");
      fprintf(in_file,"!    ecut: 0         ! starting ecut value\n");
      fprintf(in_file,"!    ecut+ 0         ! ecut increment per test\n");
      fprintf(in_file,"!    getden -1       ! re-use electron density\n");
      fprintf(in_file,"     ecut %d         ! ecut value, corresponding to lambda = %lf bohr\n\n",inputCIF->ecut_final,inputCIF->ecut_lambda_bohr);
    }

    //Print geometry optimization section
    if(mode==2){
      fprintf(in_file,"!----GEOMETRY OPTIMIZATION----!\n");
      fprintf(in_file,"     ionmov 3        ! optimization algorithm, change to 2 if 3 does not work\n");
      fprintf(in_file,"     ntime 50        ! max move steps, normally okay\n");
      fprintf(in_file,"     prtden 0        ! do not output _DEN file\n");
      fprintf(in_file,"     optcell 0       ! 0 = only ionic position; 2 = total optimization\n");
      fprintf(in_file,"     dilatmx 1.225   ! max acell increase. Can set to 1 if optcell = 0\n");
      fprintf(in_file,"     ecutsm 13.6 eV  ! basis set smearing, normally okay\n\n");
    }else if(mode==3){
      fprintf(in_file,"!----GEOMETRY OPTIMIZATION----!\n");
      fprintf(in_file,"     ionmov 3        ! optimization algorithm, change to 2 if 3 does not work\n");
      fprintf(in_file,"     ntime 50        ! max move steps, normally okay\n");
      fprintf(in_file,"     prtden 0        ! do not output _DEN file\n");
      fprintf(in_file,"     optcell 2       ! 0 = only ionic position; 2 = total optimization\n");
      fprintf(in_file,"     dilatmx 1.225   ! max acell increase. Can set to 1 if optcell = 0\n");
      fprintf(in_file,"     ecutsm 13.6 eV  ! basis set smearing, normally okay\n\n");
    }else{
      fprintf(in_file,"!----GEOMETRY OPTIMIZATION----!\n");
      fprintf(in_file,"!    ionmov 3        ! optimization algorithm, change to 2 if 3 does not work\n");
      fprintf(in_file,"!    ntime 50        ! max move steps, normally okay\n");
      fprintf(in_file,"!    prtden 0        ! do not output _DEN file\n");
      fprintf(in_file,"!    optcell 0       ! 0 = only ionic position; 2 = total optimization\n");
      fprintf(in_file,"!    dilatmx 1.225   ! max acell increase. Can set to 1 if optcell = 0\n");
      fprintf(in_file,"!    ecutsm 13.6 eV  ! basis set smearing, normally okay\n\n");
    }
    
    //Print spin polarization section

    if(SOC_flag==2){
      fprintf(in_file,"!----SPIN POLARIZATION----!\n");
      fprintf(in_file,"!    nsppol 2        ! 1 = non spin polarized 2 = spin polarized\n");
      fprintf(in_file,"     nspinor 2       ! 1 = no spin-orbit coupling 2 = spin orbit coupling\n");
      fprintf(in_file,"!    nspden 2        ! 1 = no spin 2 = scalar 4 = vector magnetization\n");
      fprintf(in_file,"!    spinat %d*1.0   ! initial magnetization, unit: hbar/2\n\n",inputCIF->natom);
    }else{
      fprintf(in_file,"!----SPIN POLARIZATION----!\n");
      fprintf(in_file,"!    nsppol 2        ! 1 = non spin polarized 2 = spin polarized\n");
      fprintf(in_file,"!    nspinor 1       ! 1 = no spin-orbit coupling 2 = spin orbit coupling\n");
      fprintf(in_file,"!    nspden 2        ! 1 = no spin 2 = scalar 4 = vector magnetization\n");
      fprintf(in_file,"!    spinat %d*1.0   ! initial magnetization, unit: hbar/2\n\n",inputCIF->natom);
    }

    
    //Print CP calculations section
    if(mode==4){
      fprintf(in_file,"!----CP CALCULATIONS----!\n");
      fprintf(in_file,"!   boxcutmin 3         ! fining ngfft grid; normally unnecessary\n");
      fprintf(in_file,"    ndtset 3            ! number of datasets\n");
      fprintf(in_file,"    scalecart1 3*1.005  ! volume expansion\n");
      fprintf(in_file,"    scalecart2 3*1.000  ! constant volume\n");
      fprintf(in_file,"    scalecart3 3*0.995  ! volume contraction\n");
      fprintf(in_file,"    usekden 1           ! calculate kinetic density KDEN file\n");
      fprintf(in_file,"    prtkden 1           ! print kinetic density KDEN file\n");
      fprintf(in_file,"    prtpot 1            ! print KS potential POT file\n");
      fprintf(in_file,"    prtvha 1            ! print Hartree potential VHA file\n");
      fprintf(in_file,"    prtvxc 1            ! print exchange-correlation potential VXC file\n");
      fprintf(in_file,"    prtvhxc 1           ! print sum of Hatree and XC potential VHXC file\n");
      fprintf(in_file,"!    ngfft 72 72 72      ! fast Fourier transform grid, all numbers must be multiply of 12\n\n");
    }else if(mode==5){
      fprintf(in_file,"!----CP CALCULATIONS----!\n");
      fprintf(in_file,"!   boxcutmin 3         ! fining ngfft grid; normally unnecessary\n");
      fprintf(in_file,"    ndtset 3            ! number of datasets\n");
      fprintf(in_file,"    scalecart1 3*0.805  ! volume expansion\n");
      fprintf(in_file,"    scalecart2 3*0.800  ! constant volume\n");
      fprintf(in_file,"    scalecart3 3*0.795  ! volume contraction\n");
      fprintf(in_file,"    usekden 1           ! calculate kinetic density KDEN file\n");
      fprintf(in_file,"    prtkden 1           ! print kinetic density KDEN file\n");
      fprintf(in_file,"    prtpot 1            ! print KS potential POT file\n");
      fprintf(in_file,"    prtvha 1            ! print Hartree potential VHA file\n");
      fprintf(in_file,"    prtvxc 1            ! print exchange-correlation potential VXC file\n");
      fprintf(in_file,"    prtvhxc 1           ! print sum of Hatree and XC potential VHXC file\n");
      fprintf(in_file,"!    ngfft 72 72 72      ! fast Fourier transform grid, all numbers must be multiply of 12\n\n");
    }else if(mode==6){
      fprintf(in_file,"!----CP CALCULATIONS----!\n");
      fprintf(in_file,"!   boxcutmin 3         ! fining ngfft grid; normally unnecessary\n");
      fprintf(in_file,"    ndtset 3            ! number of datasets\n");
      fprintf(in_file,"    scalecart1 3*1.205  ! volume expansion\n");
      fprintf(in_file,"    scalecart2 3*1.200  ! constant volume\n");
      fprintf(in_file,"    scalecart3 3*1.195  ! volume contraction\n");
      fprintf(in_file,"    usekden 1           ! calculate kinetic density KDEN file\n");
      fprintf(in_file,"    prtkden 1           ! print kinetic density KDEN file\n");
      fprintf(in_file,"    prtpot 1            ! print KS potential POT file\n");
      fprintf(in_file,"    prtvha 1            ! print Hartree potential VHA file\n");
      fprintf(in_file,"    prtvxc 1            ! print exchange-correlation potential VXC file\n");
      fprintf(in_file,"    prtvhxc 1           ! print sum of Hatree and XC potential VHXC file\n");
      fprintf(in_file,"!    ngfft 72 72 72      ! fast Fourier transform grid, all numbers must be multiply of 12\n\n");
    }else{
      fprintf(in_file,"!----CP CALCULATIONS----!\n");
      fprintf(in_file,"!    boxcutmin 3         ! fining ngfft grid; normally unnecessary\n");
      fprintf(in_file,"!    ndtset 3            ! number of datasets\n");
      fprintf(in_file,"!    scalecart1 3*1.005  ! volume expansion\n");
      fprintf(in_file,"!    scalecart2 3*1.000  ! constant volume\n");
      fprintf(in_file,"!    scalecart3 3*0.995  ! volume contraction\n");
      fprintf(in_file,"!    usekden 1           ! calculate kinetic density KDEN file\n");
      fprintf(in_file,"!    prtkden 1           ! print kinetic density KDEN file\n");
      fprintf(in_file,"!    prtpot 1            ! print KS potential POT file\n");
      fprintf(in_file,"!    prtvha 1            ! print Hartree potential VHA file\n");
      fprintf(in_file,"!    prtvxc 1            ! print exchange-correlation potential VXC file\n");
      fprintf(in_file,"!    prtvhxc 1           ! print sum of Hatree and XC potential VHXC file\n");
      fprintf(in_file,"!    ngfft 72 72 72      ! fast Fourier transform grid, all numbers must be multiply of 12\n\n");
    }

    //Print DOS calculation section
    if(mode==7){
      fprintf(in_file,"!----DOS CALCULATION----!\n");
      fprintf(in_file,"!    ndtset 2                 ! number of datasets\n");
      fprintf(in_file,"     prtdos 2                 ! output total DOS, tetrahedron method\n");
      fprintf(in_file,"!    prtdos2 3                ! output dos by atom (projected DOS)\n");
      fprintf(in_file,"!    ratsph (Calc. Radius)    ! calculated atomic radius\n");
      fprintf(in_file,"!----END OF IN FILE----!");
    }else if(mode==8){
      fprintf(in_file,"!----DOS CALCULATION----!\n");
      fprintf(in_file,"     ndtset 2                 ! number of datasets\n");
      fprintf(in_file,"     prtdos1 2                ! output total DOS, tetrahedron method\n");
      fprintf(in_file,"     prtdos2 3                ! output dos by atom, (projected DOS)\n");
      fprintf(in_file,"     ratsph (Calc. Radius)    ! calculated atomic radius\n");
      fprintf(in_file,"!----END OF IN FILE----!");
    }else{
      fprintf(in_file,"!----DOS CALCULATION----!\n");
      fprintf(in_file,"!    ndtset 2                 ! number of datasets\n");
      fprintf(in_file,"!    prtdos1 2                ! output total DOS, tetrahedron method\n");
      fprintf(in_file,"!    prtdos2 3                ! output dos by atom, (projected DOS)\n");
      fprintf(in_file,"!    ratsph (Calc. Radius)    ! calculated atomic radius\n");
      fprintf(in_file,"!----END OF IN FILE----!");
    }
    printf("Done.\n");

    //Finished writing in file, proceed to write log file
    sprintf(LOGfilename,"prepareABINIT-%s.log",outputbase);    
    LOGfile=fopen(LOGfilename,"w");
    fprintf(LOGfile,"The abinit .in and .files file have been successfully created.\n\n");
    //printf("The abinit .in and .files file have been successfully created.\n\n");
    fprintf(LOGfile,"Options specified:\n");
    if(convert_primitive_flag==1){
      fprintf(LOGfile,"Convert your unit cell to a primitive cell if it is centered.\n");
    }else if(convert_primitive_flag==2){
      fprintf(LOGfile,"Retain your unit cell if it is centered.\n");
    }
    if(scvo_flag==1){
      fprintf(LOGfile,"Use semicore pseudopotential only for 11 or 12 electron d-metals.\n");
    }else if(scvo_flag==2){
      fprintf(LOGfile,"Use semicore pseudopotential whenever available.\n");
    }
    if(kmesh_flag==1){
      fprintf(LOGfile,"Use ngkpt for k-point mesh.\n");
    }else if(scvo_flag==2){
      fprintf(LOGfile,"Use kptrlatt for k-point mesh.\n");
    }
    if(SOC_flag==1){
      fprintf(LOGfile,"Do not enable spin-orbit coupling.\n");
    }else if(SOC_flag==2){
      fprintf(LOGfile,"Enable spin-orbit coupling.\n");
    }
    
 
    //fprintf(LOGfile,"Below are the steps you need to follow before you start the calculation. To ensure a smooth calculation, please read carefully.\n\n\n");
    //printf("Below are the steps you need to follow before you start the calculation. To ensure a smooth calculation, please read carefully.\n\n\n");
    if(mode==1||mode==2||mode==9){
      if(mode==1){
        fprintf(LOGfile,"\n****************** NEXT STEPS ******************\n\n");
        printf("\n****************** NEXT STEPS ******************\n\n");
        fprintf(LOGfile,"You chose automatic generation of the k-point mesh for the system.\n\n");
        printf("You chose automatic generation of the k-point mesh for the system.\n\n");
        //fprintf(LOGfile,"You chose to generate K-point mesh for target system.\n\n");
        //printf("You chose to generate K-point mesh for target system.\n\nThis means you are staring a new calculation (No previous calculation should have been done on this system.)\n\n");
        if(kmesh_flag==1){
          fprintf(LOGfile,"You chose to specify the k-point mesh with the variable ngkpt (default). The .in file has been populated with 5 choices.");
          printf("You chose to specify the k-point mesh with the variable ngkpt (default). The .in file has been populated with 5 choices.");
          fprintf(LOGfile,"Some of these might be commented out for practical reasons.");
          printf("Some of these might be commented out for practical reasons.");
          fprintf(LOGfile,"You may uncomment the ones you wish to include, adjusting ndtset accordingly.\n\n");
          printf("You may uncomment the ones you wish to include, adjusting ndtset accordingly.\n\n");
          fprintf(LOGfile,"After you deal with any warnings below, you may submit your ABINIT job.\n\n");
          printf("After you deal with any warnings below, you may submit your ABINIT job.\n\n");
          fprintf(LOGfile,"Once the calculation is finished, you can find energy of different ngkpt grids by using [grep etotal \"base name.log\"].");
          printf("Once the calculation is finished, you can find energy of different ngkpt grids by using [grep etotal \"base name.log\"].");
          fprintf(LOGfile,"You can then select a ngkpt grid based on this data, and update the ABINIT .in file with your choice for future calculations.");
          printf("You can then select a ngkpt grid based on this data, and update the ABINIT .in file with your choice for future calculations.");
          fprintf(LOGfile,"Convergence of the energy differences for different k-point meshes to lower than 5 meV/atom (0.00018Ha/atom) is recommended.\n\n");
          printf("Convergence of the energy differences for different k-point meshes to lower than 5 meV/atom (0.00018Ha/atom) is recommended.\n\n");
        }else if(kmesh_flag==2){
          fprintf(LOGfile,"You chose to specify k-point mesh with the variable kptrlatt.\n\n");
          printf("You chose to specify k-point mesh with the variable kptrlatt.\n\n");
          fprintf(LOGfile,"After you deal with any warnings below, you may submit your abinit job.\n\n");
          printf("After you deal with any warnings below, you may submit your abinit job.\n\n");
          //fprintf(LOGfile,"To avoid a weird bug, it is recommended to use one processor with -np 1.\n\n");
          //printf("To avoid a weird bug, it is recommended to use one processor with -np 1.\n\n");
          fprintf(LOGfile,"Once the calculation is finished, you can find a list of ABINIT recommended grid with merit facotrs in the ABINIT .log file.\n\n");
          printf("Once the calculation is finished, you can find a list of ABINIT recommended grid with merit facotrs in the ABINIT .log file.\n\n");
          fprintf(LOGfile,"You can then select a kptrlatt grid with its shiftk based on this list, and then update ABINIT .in file with your choice for future calculations.\n\n");
          printf("You can then select a kptrlatt grid with its shiftk based on this list, and then update ABINIT .in file with your choice for future calculations.\n\n");
        }
      }else if(mode==2){
        fprintf(LOGfile,"****************** NEXT STEPS ******************\n\n");
        printf("****************** NEXT STEPS ******************\n\n");
        fprintf(LOGfile,"You chose to run the optimizition of the geometry of the cell by only changing the positions of the atoms.\n\n");
        printf("You chose to run the optimizition of the geometry of the cell by only changing the positions of the atoms.\n\n");
        if(kmesh_flag==1){
          fprintf(LOGfile,"You chose to specify K-point mesh with ngkpt (default), Please enter your K-point mesh. Or you can use one of the meshes generated by the program\n\n");
          printf("You chose to specify K-point mesh with ngkpt (default), Please enter your K-point mesh. Or you can use one of the meshes generated by the program\n\n");
          fprintf(LOGfile,"After you deal with all the Warnings (if there is any), you may submit your abinit job.\n\n");
          printf("After you deal with all the Warnings (if there is any), you may submit your abinit job.\n\n");
        }else if(kmesh_flag==2){
          fprintf(LOGfile,"You chose to specify K-point mesh with kptrlatt. Please enter the value of kptrlatt and shiftk you chose based on the ABINIT k-point mesh generation calculation\n\n");
          printf("You chose to specify K-point mesh with kptrlatt. Please enter the value of kptrlatt and shiftk you chose based on the ABINIT k-point mesh generation calculation\n\n");
          fprintf(LOGfile,"After you deal with all the Warnings (if there is any), you may submit your abinit job.\n\n");
          printf("After you deal with all the Warnings (if there is any), you may submit your abinit job.\n\n");
        }
      }else if(mode==9){
        fprintf(LOGfile,"This calculation generate the files for ecut determination.\n\n");
        fprintf(LOGfile,"By default, the calculation will start with ecut = 30 Ha and calculate total energy at increment of 10 Ha for 15 data sets.\n\n");
        fprintf(LOGfile,"Therefore, it covers the range of 30 Ha - 180 Ha, which should be more than sufficient for any pseudopotential\n\n");
        fprintf(LOGfile,"You may reduce the number of data sets to speed up this calculation.\n\n");
        fprintf(LOGfile,"You may now submit your calculation\n");
        fprintf(LOGfile,"\n------END OF LOG FILE------");
        printf("This calculation generate the files for ecut determination.\n\n");
        printf("By default, the calculation will start with ecut = 30 Ha and calculate total energy at increment of 10 Ha for 15 data sets.\n\n");
        printf("Therefore, it covers the range of 30 Ha - 180 Ha, which should be more than sufficient for any pseudopotential\n\n");
        printf("You may reduce the number of data sets to speed up this calculation.\n\n");
        printf("You may now submit your calculation\n");
        if(kmesh_flag==1){
          fprintf(LOGfile,"You chose to specify K-point mesh with ngkpt (default), Please enter your K-point mesh. Or you can use one of the meshes generated by the program\n\n");
          printf("You chose to specify K-point mesh with ngkpt (default), Please enter your K-point mesh. Or you can use one of the meshes generated by the program\n\n");
          fprintf(LOGfile,"After you deal with all the Warnings (if there is any), you may submit your abinit job.\n\n");
          printf("After you deal with all the Warnings (if there is any), you may submit your abinit job.\n\n");
        }else if(kmesh_flag==2){
          fprintf(LOGfile,"You chose to specify K-point mesh with kptrlatt. Please enter the value of kptrlatt and shiftk you chose based on the ABINIT k-point mesh generation calculation\n\n");
          printf("You chose to specify K-point mesh with kptrlatt. Please enter the value of kptrlatt and shiftk you chose based on the ABINIT k-point mesh generation calculation\n\n");
          fprintf(LOGfile,"After you deal with all the Warnings (if there is any), you may submit your abinit job.\n\n");
          printf("After you deal with all the Warnings (if there is any), you may submit your abinit job.\n\n");
        }
      }
      if(convert_FAIL==1) Warning++;
      if(inputCIF->spgaxor==-1) Warning++;
      if(spgorig_read==0){
        if(inputCIF->spgorig==3) Warning++;
      }
      for(i=0;i<ntype;i++){
        if(ecut_unknown[i]==1) Warning++;
      }

      for(i=0;i<inputCIF->total_site;i++){
        if(inputCIF->occupancy[i]<1) Warning++;
      }

      fprintf(LOGfile,"\n****************** WARNINGS ******************\n\n");
      printf("\n\n\n****************** WARNINGS ******************\n\n");
      if(Warning==0){
        fprintf(LOGfile,"There are no warnings!\n\n");
        printf("There are no warnings!\n\n");
        fprintf(LOGfile,"However, since this is a new calculation, it is still a good idea to check the input file to be sure the CCIF was properly read before you proceed.\n\n");
        printf("However, since this is a new calculation, it is still a good idea to check the input file to be sure the CCIF was properly read before you proceed.\n\n");
        //fprintf(LOGfile,"Some of the system info might not be properly acquired due to CIF files can be written in different formats.\n\n");
        //printf("Some of the system info might not be properly acquired due to CIF files can be written in different formats.\n\n");
        fprintf(LOGfile,"You may now submit your calculation.\n\n");
        printf("You may now submit your calculation.\n\n");
        fprintf(LOGfile,"\n------END OF LOG FILE------");
      }else{
        fprintf(LOGfile,"There are %d warnings you need to deal with before start ABINIT calculation. \n\n",Warning);
        printf("There are %d warnings you need to deal with before start ABINIT calculation. \n\n",Warning);
      }
      if(spgorig_read==0){
        if(inputCIF->spgorig==3){
          Current_Warning++;
          fprintf(LOGfile,"Warning %d: Your space group (%d) have two possible origin settings.\n\n",Current_Warning, inputCIF->spgroup);
          printf("Warning %d: Your space group (%d) have two possible origin settings.\n\n",Current_Warning, inputCIF->spgroup);
          fprintf(LOGfile,"To do: change your spgorig value to 1 or 2 before running ABINIT.\n\n");
          printf("To do: change your spgorig value to 1 or 2 before running ABINIT.\n\n");
        }
      }
      if(inputCIF->spgaxor==-1){
        Current_Warning++;
        fprintf(LOGfile, "Warning %d: For your space group %d, the standard setting according to International Tables for Crystallography is %s \n\n", Current_Warning,inputCIF->spgroup,inputCIF->HM_STDName);
        fprintf(LOGfile, "However, you have setting %s. This means the spgaxor cannot be default value (1) in the .in file.\n\n",inputCIF->HM_FullspgName_oneword);
        fprintf(LOGfile, "Please choose the correct spgaxor according to your unit cell setting or generate an output file in the setting %s.\n\n",inputCIF->HM_STDName);
        fprintf(LOGfile, "To do: You have two options:\n");
        fprintf(LOGfile, "1. Choose correct spgaxor and change it manually in the .in file. You probably need information on https://docs.abinit.org/guide/spacegroup/ \n");
        fprintf(LOGfile, "2. Load your CIF into a software like Diamond and then output a new CIF file in the standard setting %s \n\n",inputCIF->HM_STDName);
        printf("Warning %d: For your space group %d, the standard setting according to International Tables for Crystallography is %s \n\n", Current_Warning,inputCIF->spgroup,inputCIF->HM_STDName);
        printf("However, you have setting %s. This means the spgaxor cannot be default value (1) in the .in file.\n\n",inputCIF->HM_FullspgName_oneword);
        printf("Please choose the correct spgaxor according to your unit cell setting or generate an output file in the setting %s.\n\n",inputCIF->HM_STDName);
        printf("To do: You have two options:\n");
        printf("1. Choose correct spgaxor and change it manually in the .in file. You probably need information on https://docs.abinit.org/guide/spacegroup/ \n");
        printf("2. Load your CIF into a software like Diamond and then output a new CIF file in the standard setting %s \n\n",inputCIF->HM_STDName);
      }
      if(convert_FAIL==1){
        Current_Warning++;
        fprintf(LOGfile,"Your CIF file contains a lattice that is not one of the 14 Bravis lattices. This will mess up the ngkpt generation step.\n\n");
        fprintf(LOGfile,"The 14 Bravis lattices are: P-triclinic, P-monoclinic, C-monoclinic, P-orthorhombic, C-orthorhombic, I-orthorhombic, F-orthorhombic, P-tetragonal, I-tetragonal, P-trigonal, P-hexagonal, P-cubic, F-cubic, I-cubic. \n\n");
        fprintf(LOGfile,"The C- means based centered. It includes A-centered, B-centered and C-centered.\n\n");
        fprintf(LOGfile,"Change your CIF file so it contains an unit cell that is one of the 14 Bravis lattices and rerun prepareABINIT.\n\n");
        fprintf(LOGfile,"To do: Check your CIF file. If it is correct, manually calculate the primitive cell or rerun this program and select not convert to primitive");
        printf("Your CIF file contains a lattice that is not one of the 14 Bravis lattices. This will mess up the ngkpt generation step.\n\n");
        printf("The 14 Bravis lattices are: P-triclinic, P-monoclinic, C-monoclinic, P-orthorhombic, C-orthorhombic, I-orthorhombic, F-orthorhombic, P-tetragonal, I-tetragonal, P-trigonal, P-hexagonal, P-cubic, F-cubic, I-cubic. \n\n");
        printf("The C- means based centered. It includes A-centered, B-centered and C-centered.\n\n");
        printf("Change your CIF file so it contains an unit cell that is one of the 14 Bravis lattices and rerun prepareABINIT.\n\n");
        printf("To do: Check your CIF file. If it is correct, manually calculate the primitive cell or rerun this program and select not convert to primitive");
      }
      for(i=0;i<ntype;i++){
        if(ecut_unknown[i]==1){
          Current_Warning++;
          fprintf(LOGfile,"Warning %d: the specified psp %s for element %s does not have a pre-calculated ecut value.\n\n",Current_Warning, inputCIF->psp_name[i], inputCIF->element_identity[i]);
          printf("Warning %d: the specified psp %s for element %s does not have a pre-calculated ecut value.\n\n",Current_Warning, inputCIF->psp_name[i], inputCIF->element_identity[i]);
          fprintf(LOGfile,"To do: Run the ecut determination calculation or make sure your ecut value is large enough for an accurate calculation.\n\n");
          printf("To do: Run the ecut determination calculation or make sure your ecut value is large enough for an accurate calculation.\n\n");
        }
      }

      for(i=0;i<inputCIF->total_site;i++){
        if(inputCIF->occupancy[i]<1){
          Current_Warning++;
          fprintf(LOGfile,"Warning %d: Site %d has non zero occupancy %lf. It is currently written into the .in file. \n\n",Current_Warning, (i+1), inputCIF->occupancy[i]);
          printf("Warning %d: Site %d has non zero occupancy %lf. It is currently written into the .in file.\n\n",Current_Warning, (i+1), inputCIF->occupancy[i]);
          fprintf(LOGfile,"To do: It is recommended to clean up your CIF files and only keep atoms you want and change their occupancies to 1 to create a ordered model for your system.\nAlternativly, you can also make sure you only include atoms you want by manually changing the .in file. However, only do this if you are absolutely confident as you need to change more than one entries.\n\n");
          printf("To do: It is recommended to clean up your CIF files and only keep atoms you want and change their occupancies to 1 to create a ordered model for your system.\nAlternativly, you can also make sure you only include atoms you want by manually changing the .in file. However, only do this if you are absolutely confident as you need to change more than one entries.\n\n");
        }
      }
   
      if(Warning!=0){
        fprintf(LOGfile,"Above are all the warnings. \n\n");
        fprintf(LOGfile,"However, since this is a new calculation (not based on previous input and output files), it is still recommended to carefully check the input file before you proceed.\n\n");
        fprintf(LOGfile,"Some of the system info might not be properly acuqired due to CIF files can be written in different formats\n\n");
        
        fprintf(LOGfile,"You may submit your calculation after resolving all the warnings.\n\n");
        fprintf(LOGfile,"\n------END OF LOG FILE------");
        printf("Above are all the warnings. \n\n");
        printf("However, since this is a new calculation (not based on previous input and output files), it is still recommended to carefully check the input file before you proceed.\n\n");
        printf("Some of the system info might not be properly acuqired due to CIF files can be written in different formats\n\n");
        printf("You may submit your calculation after resolving all the warnings.\n\n");
      }
    }else if(mode>2&&mode!=9){
      fprintf(LOGfile,"Your calculation is based on a previous calculation.\n\n");
      fprintf(LOGfile,"All the info for this calculation are acuqired from the ABINIT output file and previous input file.\n\n");
      fprintf(LOGfile,"Therefore, there will be no warnings because all info must be correct in order for the previous calculation to finish.\n\n");
      fprintf(LOGfile,"******************INSCTRUTION******************\n\n");
      printf("Your calculation is based on a previous calculation.\n\n");
      printf("All the info for this calculation are acuqired from the ABINIT output file and previous input file.\n\n");
      printf("Therefore, there will be no warnings because all info must be correct in order for the previous calculation to finish.\n\n");
      printf("******************INSCTRUTION******************\n\n");
      if(mode==3){
        fprintf(LOGfile,"This calculation fully optimizes the geometry of the unit cell\n\n");
        fprintf(LOGfile,"The xred has been updated based on the output of the previous run, which only optimizes the atomic positions.\n\n");
        fprintf(LOGfile,"You may now submit your calculation");
        fprintf(LOGfile,"\n------END OF LOG FILE------");
        printf("This calculation fully optimizes the geometry of the unit cell\n\n");
        printf("The xred has been updated based on the output of the previous run, which only optimizes the atomic positions.\n\n");
        printf("You may now submit your calculation");
      }else if(mode==4){
        fprintf(LOGfile,"This calculation generate the files for CP calulation at equilibrium volume\n\n");
        fprintf(LOGfile,"The acell and xred have been updated based on the output of the the geometry optimization.\n\n");
        fprintf(LOGfile,"You may now submit your calculation\n");
        fprintf(LOGfile,"\n------END OF LOG FILE------");
        printf("This calculation generate the files for CP calulation at equilibrium volume\n\n");
        printf("The acell and xred have been updated based on the output of the the geometry optimization.\n\n");
        printf("You may now submit your calculation\n");
      }else if(mode==5){
        fprintf(LOGfile,"This calculation generate the files for CP calibration at contrated volume (V_con = 0.8 V_eq)\n\n");
        fprintf(LOGfile,"The acell and xred have been updated based on the output of the the geometry optimization.\n\n");
        fprintf(LOGfile,"You may now submit your calculation\n");
        fprintf(LOGfile,"\n------END OF LOG FILE------");
        printf("This calculation generate the files for CP calibration at contrated volume (V_con = 0.8 V_eq)\n\n");
        printf("The acell and xred have been updated based on the output of the the geometry optimization.\n\n");
        printf("You may now submit your calculation\n");
      }else if(mode==6){
        fprintf(LOGfile,"This calculation generate the files for CP calibration at expanded volume (V_exp = 1.2 V_eq)\n\n");
        fprintf(LOGfile,"The acell and xred have been updated based on the output of the geometry optimization.\n\n");
        fprintf(LOGfile,"You may now submit your calculation\n");
        fprintf(LOGfile,"\n------END OF LOG FILE------");
        printf("This calculation generate the files for CP calibration at expanded volume (V_exp = 1.2 V_eq)\n\n");
        printf("The acell and xred have been updated based on the output of the geometry optimization.\n\n");
        printf("You may now submit your calculation\n");
      }else if(mode==7){
        fprintf(LOGfile,"This calculation generate the files for total DOS calculation, no projected DOS will be calculated.\n\n");
        fprintf(LOGfile,"The acell and xred have been updated based on the output of the geometry optimization.\n\n");
        fprintf(LOGfile,"The K-point mesh has been multiplied by a factor of 3 in order to fully capture the DOS.\n\n");
        fprintf(LOGfile,"You may now submit your calculation\n");
        printf("This calculation generate the files for total DOS calculation, no projected DOS will be calculated.\n\n");
        printf("The acell and xred have been updated based on the output of the geometry optimization.\n\n");
        printf("The K-point mesh has been multiplied by a factor of 3 in order to fully capture the DOS.\n\n");
        printf("You may now submit your calculation\n");
        fprintf(LOGfile,"\n------END OF LOG FILE------");
      }else if(mode==8){
        fprintf(LOGfile,"This calculation generate the files for projected DOS calculation.\n\n");
        fprintf(LOGfile,"You need to manually calculate the ratsph before submitting you calculation.\n\n");
        fprintf(LOGfile,"The K-point mesh has been multiplied by a factor of 3 in order to fully capture the DOS.\n");
        printf("This calculation generate the files for projected DOS calculation.\n\n");
        printf("You need to manually calculate the ratsph before submitting you calculation.\n\n");
        printf("The K-point mesh has been multiplied by a factor of 3 in order to fully capture the DOS.\n");
        fprintf(LOGfile,"\n------END OF LOG FILE------");
      }else if(mode==9){
        fprintf(LOGfile,"This calculation generate the files for ecut determination.\n\n");
        fprintf(LOGfile,"By default, the calculation will start with ecut = 30 Ha and calculate total energy at increment of 10 Ha for 15 data sets.\n\n");
        fprintf(LOGfile,"Therefore, it covers the range of 30 Ha - 180 Ha, which should be more than sufficient for any pseudopotential\n\n");
        fprintf(LOGfile,"You may reduce the number of data sets to speed up this calculation.\n\n");
        fprintf(LOGfile,"You may now submit your calculation\n");
        fprintf(LOGfile,"\n------END OF LOG FILE------");
        printf("This calculation generate the files for ecut determination.\n\n");
        printf("By default, the calculation will start with ecut = 30 Ha and calculate total energy at increment of 10 Ha for 15 data sets.\n\n");
        printf("Therefore, it covers the range of 30 Ha - 180 Ha, which should be more than sufficient for any pseudopotential\n\n");
        printf("You may reduce the number of data sets to speed up this calculation.\n\n");
        printf("You may now submit your calculation\n");
      }
    }
    fclose(LOGfile);
}

int main (int argc, char * argv[]){
    //declare variables
    char CIFname[100];
    char outputbase[100];
    char outputINname[100];
    char outputFILESname[100];
    char PhoenixUserName[100];
    char psp_address_header[200];
    char abinitOutName[100];
    char abinitOldInName[100];
    char optionsName[100];
    FILE * CIFfile;
    FILE * files_file;
    FILE * in_file;
    FILE * old_in_file;
    FILE * abinitout_file;
    FILE * options_file;
    int calculation_mode;  //1~12.
    int default_flag;
    int fining_factor;
    double TriHex_tolerance;
    //Check program input
    fining_factor=3;
    TriHex_tolerance=0.0004;
 
    if(argc>=4){
      strcpy(CIFname,argv[1]);
      strcpy(outputbase,argv[2]);
      strcpy(PhoenixUserName,argv[3]);
      sprintf(outputINname,"%s.in",outputbase);
      sprintf(outputFILESname,"%s.files",outputbase);
      sprintf(psp_address_header,"/home/%s/psp/",PhoenixUserName);
      if(access(CIFname,F_OK) != -1){
        printf("Name of the CIF file is: %s \n", CIFname);
        printf("Name of the output file: %s and %s \n", outputFILESname,outputINname);
      }
      else{
        printf("The CIF file: %s does not exist. Please double check!!\n", CIFname);
        exit(0);
      }
    }
    else{
      printf("Usage: <CIF file name> <output base name> <Phoenix user name> if necessary <abinit output file> <old abinit input file>\n");
      exit(0);
    }
    

    printf("What type of calculation are you going to run?\n");
    printf("[1=K-point mesh determination]\n[2=Geometry optimization - atomic position]\n");
    printf("[3=Geometry optimization - full optimization]\n[4=CP calculation preparation - equilibrium volume]\n");
    printf("[5=CP calculation preparation - contracted]\n[6=CP calculation preparation - expanded]\n");
    printf("[7=Total DOS]\n[8=Projected DOS]\n[9=ecut determination]\n\n");
    printf("Options: ");
    while(1==1){
      scanf("%d",&calculation_mode);
      if(calculation_mode==1){
        printf("Mode: K-point mesh determination\n\n");
        break;
      }
      else if(calculation_mode==2){
        printf("Mode: Geometry optimization - only optimizing atomic positions.\n\n");
        break;
      }
      else if(calculation_mode==3){
        printf("Mode: Geometry optimization - optimizing all cell parameters.\n\n");
        break;
      }
      else if(calculation_mode==4){
        printf("Mode: CP calculation preparation - at equilibrium volume.\n\n");
        break;
      }
      else if(calculation_mode==5){
        printf("Mode: CP calculation preparation - at contracted volume.\n\n");
        sprintf(outputFILESname,"%s_con.files",outputbase);
        break;
      }
      else if(calculation_mode==6){
        printf("Mode: CP calculation preparation - at expanded volume.\n\n");
        sprintf(outputFILESname,"%s_exp.files",outputbase);
        break;
      }
      else if(calculation_mode==7){
        printf("Mode: Total DOS calculation.\n\n");
        sprintf(outputFILESname,"%s_DOS.files",outputbase);
        break;
      }
      else if(calculation_mode==8){
        printf("Mode: Projected DOS calculation.\n\n");
        sprintf(outputFILESname,"%s_DOS.files",outputbase);
        break;
      }
      else if(calculation_mode==9){
        printf("Mode: ecut determination.\n\n");
        break;
      }
      else{
        printf("Input invalid. Please re-enter. Mode must be 1 - 10.\n");
      }
    }

    printf("Do you wish to use the default settings for this calculation? [1=Yes] [2=No]   ");
    while(1==1){
      scanf("%d",&default_flag);
      if(default_flag==1){ //Use default
        printf("\nUsing default for this calculation. Here are the default:\n\n");
        printf("1. Convert to primitive cell if the given geometry is centered.\n");
        printf("2. Use semicore pseudopotential only for 11 or 12 electron d-elements.\n");
        printf("3. Use ngkpt to specify the K-points mesh.\n");
        printf("4. Do not enable spin orbit coupling.\n\n");
        break;
      }
      else if(default_flag==2){
        printf("Please specify your settings for this calculation below.\n");
        break;
      }
      else{
        printf("Input invalid. Please re-enter. [1=Use default] [2=Do not use default] ");
      }
    }
      
    if(default_flag==1){
      convert_primitive_flag=1;
      scvo_flag=1;
      kmesh_flag=1; //ngkpt
      SOC_flag=1;
      //printf("Setting up default, convert primitive flag is: %d. \n",convert_primitive_flag);
    }
    else if(default_flag==2){
      //Read in settings
      
      /*printf("Is your CIF file from ICSD or Jana? [1=ICSD] [2=Jana]   ");
      while(1==1){
        scanf("%d",&CIF_source_flag);
        if(CIF_source_flag==1){
          printf("Your CIF file is from ICSD.\n");
          break;
        }
        else if(CIF_source_flag==2){
          printf("Your CIF file is from Jana.\n");
          break;
        }
        else{
          printf("Input invalid. Please re-enter. [1=ICSD] [2=Jana]    ");
        }
      }
      */

      printf("Do you want ABINIT to convert to primitive cell if current one is centered? [1=Yes] [2=No]   ");
      while(1==1){
        scanf("%d",&convert_primitive_flag);
        if(convert_primitive_flag==1){
          printf("If your unit cell in centered, it will be converted.\n");
          break;
        }
        else if(convert_primitive_flag==2){
          printf("Unit cell will be retained.\n");
          break;
        }
        else{
          printf("Input invalid. Please re-enter. [1=Yes] [2=No]\n");
        }
      }  

      printf("What kind of pseudopotential do you wish to use? [1=semicore only for 11 or 12 electron d-elements] [2=semicore when available]   ");
      while(1==1){
        scanf("%d",&scvo_flag);
        if(scvo_flag==1){
          printf("The calculation will only use semicore psp when necessary. \n");
          break;
        }
        else if(scvo_flag==2){
          printf("The calculation will only use semicore psp when available. \n");
          break;
        }
        else{
          printf("Input invalid. Please re-enter. [1=semicore when necessary] [2=semicore when possible]    ");
        }
      }
     
      printf("How do you want to specify your k-point mesh? [1=ngkpt] [2=kptrlatt]  ");
      while(1==1){
        scanf("%d",&kmesh_flag);
        if(kmesh_flag==1){
          printf("Use ngkpt. \n");
          break;
        }
        else if(kmesh_flag==2){
          printf("Use kptrlatt. \n");
          break;
        }
        else{
          printf("Input invalid. Please re-enter. [1=ngkpt] [2=kptrlatt]    ");
        }
      }
      
      printf("Do you want to enable spin orbit coupling? [1=No] [2=Yes] ");
      while(1==1){
        scanf("%d",&SOC_flag);
        if(SOC_flag==1){
          printf("Do not enable spin orbit coupling. \n\n");
          break;
        }
        else if(SOC_flag==2){
          printf("Enable spin orbit coupling. \n\n");
          break;
        }
        /*else if(SOC_flag==3){
          printf("Enable spin orbit coupling with spin polarization. \n\n");
          break;
        }*/
        else{
          printf("Input invalid. Please re-enter. [1=No] [2=Yes]  \n\n");
        }
      }

    }

    //printf("Convert primitive flag is: %d \n", convert_primitive_flag);

    //Write a basename-options file if this is mode 1 (K-point mesh generation)
    if(calculation_mode==1 || calculation_mode==2 || calculation_mode==9){
      sprintf(optionsName,"prepareABINIT-options-%s",outputbase);
      options_file = fopen(optionsName,"w");
      fprintf(options_file,"2\n%d\n%d\n%d\n%d",convert_primitive_flag,scvo_flag,kmesh_flag,SOC_flag);
    }

   //Begin to read abinit .out file if necessary
   if(calculation_mode>2 && calculation_mode!=9){
      if(argc==6){
        strcpy(abinitOutName,argv[4]);
        printf("Abinit out file name: %s.\n", abinitOutName);
        if(access(abinitOutName, F_OK) == -1){
          printf("The output file you specified does not exit. \n");
          exit(0);
        }
        
        strcpy(abinitOldInName,argv[5]);
        printf("Abinit old input file name: %s. \n", abinitOldInName);
        if(access(abinitOldInName, F_OK) == -1){
          printf("The old input file you specified does not exit. \n");
          exit(0);
        }
      }else{
        printf("Your mode of calculation requires abinit output file and old input file.\n");
        printf("Usage: <CIF file name> <output base name> <Phoenix user name> <abinit output file> <abinit old input file>\n");
        printf("Decision taken to exit. \nFILES ARE NOT BEING GENERATED!!!!\n");
        exit(0);
      }
    }

    //Calculation mode = 1 or = 2 or = 9 means this is a new calculation
    if((calculation_mode==1) || (calculation_mode==2) || (calculation_mode==9)){
      //Begin reading CIF
      CIFfile = fopen(CIFname,"r");
      readCIF(&inputCIF,CIFfile);
      fclose(CIFfile);
      printf("Done reading CIF file. \n\n");
    
      //Begin processing raw CIF format data
      ProcessCIFdata(&inputCIF,calculation_mode);
      //printf("Done processing CIF data. \n\n");
    
      //Invoke the ngkpt generation step
      if(kmesh_flag==1){
        //printf("Entering ngkpt calculation. \n");
        generate_ngkpt(&inputCIF);
      }
      
      TriHex_Adjustment(&inputCIF,TriHex_tolerance);
    }

    //Other calculaiton mode means this is based on an old calculation
    //Read and process data in .in and .out file
    if(calculation_mode>2 && calculation_mode!=9){
      old_in_file = fopen(abinitOldInName,"r");
      read_old_in(&inputCIF,old_in_file);
      fclose(old_in_file);
      abinitout_file = fopen(abinitOutName,"r");
      read_abinitout(&inputCIF,abinitout_file);
      fclose(abinitout_file);
    }
    
    //If the calculaiton is a DOS or pDOS calculation, increase the K-point mesh to make it finer
    if(calculation_mode==7 || calculation_mode==8){
      DOS_calculation_adjustment(&inputCIF,fining_factor);
    }
    //Write .files file
    files_file = fopen(outputFILESname,"w");
    Write_files_file(&inputCIF,files_file,psp_address_header,outputbase);
    fclose(files_file);
    //printf("Finished writing the .files file \n");

    //Write .in file
    in_file = fopen(outputINname,"w");
    Write_in_file(&inputCIF,in_file,outputbase,CIFname,calculation_mode);
    fclose(in_file);
    //printf("prepareABINIT successful. \n");
    printf("Make sure to carefully read the NEXT STEPS and WARNINGS before running ABINIT. \n");
    printf("These information is also available in the prepareABINIT-basename.log file. \n");
}
