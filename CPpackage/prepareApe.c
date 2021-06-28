#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#define PI 3.14159265358979323846264338328
#define ntype_max 20
#define natom_max 100
#define nsym_max 500
#define nsite_max 100

double sym_tol=0.01; //symmetry tolrance in reduced coordinates.

struct Crystal{
  int natom, ntype, nsym, nsite;
  double symrel[nsym_max][3][3];
  double trans_vec[nsym_max][3];
  int znucl[ntype_max];
  double xred[natom_max], yred[natom_max], zred[natom_max];
  double atom_coord[natom_max][3]; //Adjusted to be in-between 0 and 1
  int typat[natom_max];
  char element_name[ntype_max][10];
  char psp_name[ntype_max][50];
  int psp_e[ntype_max];
  double eqv_coord[natom_max][nsym_max][3];  //Logic: all symmetry related position of atom [i] is stored here. 
  int eqv_atom_No[natom_max][natom_max];
  int uneqv_atom[natom_max];  //Atom must be treated individually bc it is not equivlent to any other atom except itself
  int eqv_multi[natom_max];
  double bader_e_atom[natom_max];
  double bader_e_site[nsite_max];
  double bader_charge[nsite_max];
  double psp_e_atom[nsite_max];
  int site_element_No[nsite_max];
  double bader_100[nsite_max];
  double bader_75[nsite_max];
  double bader_50[nsite_max];
  double bader_25[nsite_max];
  int noble_gas_id[nsite_max];
  double e_left[nsite_max];
  double e1s,e2s,e2p,e3s,e3p,e4s,e3d,e4p,e5s,e4d,e5p,e6s,e5d,e4f,e6p,e7s,e6d,e7p;
} Cryst;

double max(double n1, double n2){
  if(n1>=n2){
    return n1;
  }else{
    return n2;
  }
}

double min(double n1, double n2){
  if(n1<=n2){
    return n1;
  }else{
    return n2;
  }
}

void find_sym_atom(struct Crystal * Cryst){
  int eqv_found[natom_max];
  int i=0, j=0, k=0, l=0;
  double dist_total=0.0, dist_final=0.0;
  double dist[3];
  int eqv_count=0;
  int uneqv_count=0;
  double temp_x=0.0, temp_y=0.0, temp_z=0.0;
  for(i=0;i<natom_max;i++){
    eqv_found[i]=0;
  }
  for(i=0;i<3;i++){
    dist[i]=0.0;
  }
  for(i=0;i<Cryst->natom;i++){
    for(j=0;j<Cryst->nsym;j++){
      for(k=0;k<3;k++){
        temp_x=temp_x+Cryst->symrel[j][0][k]*Cryst->atom_coord[i][k]+Cryst->trans_vec[j][0];
        temp_y=temp_y+Cryst->symrel[j][1][k]*Cryst->atom_coord[i][k]+Cryst->trans_vec[j][1];
        temp_z=temp_z+Cryst->symrel[j][2][k]*Cryst->atom_coord[i][k]+Cryst->trans_vec[j][2];
      }
      while(temp_x<0){
        temp_x=temp_x+1;
      }
      while(temp_x>=1){
        temp_x=temp_x-1;
      }
      while(temp_y<0){
        temp_y=temp_y+1;
      }
      while(temp_y>=1){
        temp_y=temp_y-1;
      }
      while(temp_z<0){
        temp_z=temp_z+1;
      }
      while(temp_z>=1){
        temp_z=temp_z-1;
      }
      if(i==0){
      //printf("For atom %d, its position is: %lf %lf %lf, symrel %d is [%lf] [%lf] [%lf]; [%lf] [%lf] [%lf]; [%lf] [%lf] [%lf]\n",(i+1),Cryst->atom_coord[i][0],Cryst->atom_coord[i][1],Cryst->atom_coord[i][2],(j+1),Cryst->symrel[j][0][0],Cryst->symrel[j][0][1],Cryst->symrel[j][0][2],Cryst->symrel[j][1][0],Cryst->symrel[j][1][1],Cryst->symrel[j][1][2],Cryst->symrel[j][2][0],Cryst->symrel[j][2][1],Cryst->symrel[j][2][2]);
      //printf("The corresponding symmetry pos is: %lf %lf %lf \n",temp_x,temp_y,temp_z);
      }
      Cryst->eqv_coord[i][j][0]=temp_x;
      Cryst->eqv_coord[i][j][1]=temp_y;
      Cryst->eqv_coord[i][j][2]=temp_z;
      temp_x=0.0;
      temp_y=0.0;
      temp_z=0.0;
      //}
    }
  }
  for(i=0;i<Cryst->natom;i++){
    if(eqv_found[i]==0){
      //printf("Atom %d has not been found by symmetry in the past.\n",(i+1));
      Cryst->uneqv_atom[uneqv_count]=i;
    }else{
      //printf("Atom %d has already been accounted for by symmetry.\n",(i+1));
      //printf("Proceed to next atom\n");
      continue;
    }
    //printf("Starting to find atom symmetry for atom %d.\n",(i+1));
    for(j=0;j<Cryst->nsym;j++){
      for(k=0;k<Cryst->natom;k++){
        if(eqv_found[k]!=0){
          continue;
        }
        for(l=0;l<3;l++){
          dist[l]=Cryst->eqv_coord[i][j][l]-Cryst->atom_coord[k][l];
          /*if(i==5){
            printf("Cryst->eqv_coord[i][j][l]: %lf\n",Cryst->eqv_coord[i][j][l]);
            printf("Cryst->atom_coord[k][l]: %lf\n",Cryst->atom_coord[k][l]);
            printf("dist[%d] is: %lf\n",(l+1),dist[l]);
          }*/
          dist_total=dist_total+dist[l]*dist[l];
          dist_final=sqrt(dist_total);
        }
        //if(i==5) printf("Distance between atom %d's %d symmetry eqv position and atom %d: %lf.\n",(i+1),(j+1),(k+1),dist_final);
        if(dist_final<0.01){
          Cryst->eqv_atom_No[uneqv_count][eqv_count]=k;
          eqv_count++;
          eqv_found[k]=1;
          //printf("Atom %d and %d are equivlent by symmetry. \n",(i+1),(k+1));
          break;
        }
        dist_total=0.0;
        dist_final=0.0;
      }
    }
    Cryst->eqv_multi[uneqv_count]=eqv_count;
    //printf("There are %d equivlent sites for atom %d in total.\n",Cryst->eqv_multi[uneqv_count],(i+1)); 
    uneqv_count++;
    eqv_count=0;
  }
  Cryst->nsite=uneqv_count;
  printf("There are %d sites in total.\n",Cryst->nsite);
  for(i=0;i<Cryst->nsite;i++){
    //printf("Site %d, multiplicity: %d, eqv atoms:",(i+1),Cryst->eqv_multi[i]);
    for(j=0;j<Cryst->eqv_multi[i];j++){
      //printf(" %d, ",(Cryst->eqv_atom_No[i][j]+1));
    }
    printf("\n");
  }
}

void read_files_file(struct Crystal * Cryst, FILE * files_file){
  char dummy[300];
  int line_count=0;
  int psp_count=0;
  int letter_count=0;
  int i=0,j=0,k=0;
  int psp_len;
  int e_No_found=0; //0=not found; 1=found
  char temp[ntype_max][5];
  int slash_pos=0;
  memset(dummy,0,300);
  while( fgets (dummy,300,files_file) !=NULL ){
    line_count++;
    //printf("dummy: %s",dummy);
    if(line_count>=6){
      //printf("dummy: %s",dummy);
      strcpy(Cryst->psp_name[psp_count],dummy);
      psp_count++;
      if(psp_count==Cryst->ntype) break;
    }  
  }
  
  for(i=0;i<Cryst->ntype;i++){
    printf("element type: %d - psp name: %s \n",(i+1),Cryst->psp_name[i]);
    psp_len=strlen(Cryst->psp_name[i]);
    e_No_found=0;
    letter_count=0;
    for(j=0;j<psp_len;j++){
      if(Cryst->psp_name[i][j]=='.'){
        for(k=j+1;k<psp_len;k++){
          if(isdigit(Cryst->psp_name[i][k])){
            //printf("Current letter: %c \n",Cryst->psp_name[i][k]);
            temp[i][letter_count]=Cryst->psp_name[i][k];
            letter_count++;
          }else{
            e_No_found=1;
            break;
          }
        }
      }
      if(e_No_found==1){
        break;
      }
    }
  }

  letter_count=0;

  for(i=0;i<Cryst->ntype;i++){
    psp_len=strlen(Cryst->psp_name[i]);
    for(j=0;j<psp_len;j++){
      if(Cryst->psp_name[i][j]=='/'){
        slash_pos=j;
      }
    }
    //printf("Slash pos for psp %d is %d",(i+1),slash_pos);
    for(j=slash_pos;j<psp_len;j++){
      if(Cryst->psp_name[i][j]=='.'){
        break;
      }
      if(isalpha(Cryst->psp_name[i][j])){
        Cryst->element_name[i][letter_count]=Cryst->psp_name[i][j];
        letter_count++;
      }
    }
    Cryst->element_name[i][0]=Cryst->element_name[i][0]-32;
    letter_count=0;
    slash_pos=0;
    //printf("Name of element %d: %s\n",(i+1),Cryst->element_name[i]);
  }

  letter_count=0;

  //printf("Element %d name: %s \n",(i+1),Cryst->psp_name[i]);
 
  for(i=0;i<Cryst->ntype;i++){
    //printf("Temp[%d]: %s \n",(i+1),temp[i]);
    sscanf(temp[i],"%d",&Cryst->psp_e[i]);
  }

  //printf("DONE.\n");  

  for(i=0;i<Cryst->ntype;i++){
    //printf("Element type: %d - Electron in psp %s for element %d: %d \n",(i+1),Cryst->psp_name[i],Cryst->znucl[i],Cryst->psp_e[i]);
  }
}

void read_output_file(struct Crystal * Cryst, FILE * abinit_out){
  char dummy[300],dummy2[300];
  double temp=0.0;
  int i=0, j=0;
  int stop=0;
  double row_sum=0.0;
  double norm_factor=0.0;
  memset(dummy,0,300);
  memset(dummy2,0,300);
  while(stop == 0){
    if(fscanf(abinit_out,"%s",dummy) != 1 ){
      stop = 1;
      continue;
    }
    //printf("dummy is: %s \n",dummy);
    if(strcmp(dummy,"nsym")==0){
      //printf("nsym match found. \n");
      fscanf(abinit_out,"%d",&Cryst->nsym);
      //printf("Number of symmetry elements: %d\n",Cryst->nsym);
    }
    if(strcmp(dummy,"natom")==0){
      //printf("natom match found. \n");
      fscanf(abinit_out,"%d",&Cryst->natom);
      //printf("Number of atoms: %d\n",Cryst->natom);
    }
    if(strcmp(dummy,"ntypat")==0){
      //printf("ntypat match found. \n");
      fscanf(abinit_out,"%d",&Cryst->ntype);
      //printf("Number of atom types: %d\n",Cryst->ntype);
    }
  }
  
  printf("Number of symmetry elements: %d\n",Cryst->nsym);
  printf("Number of atoms: %d\n",Cryst->natom);
  printf("Number of atom types: %d\n",Cryst->ntype);

  rewind(abinit_out);
  stop = 0; 
  while(stop == 0){
    if(fscanf(abinit_out,"%s",dummy) != 1){
      stop = 1;
      continue;
    }

    if(strcmp(dummy,"symrel")==0){
      for(i=0;i<Cryst->nsym;i++){
        for(j=0;j<3;j++){
          fscanf(abinit_out,"%lf %lf %lf",&Cryst->symrel[i][0][j], &Cryst->symrel[i][1][j], &Cryst->symrel[i][2][j]);
          //row_sum=Cryst->symrel[i][j][0]*Cryst->symrel[i][j][0]+Cryst->symrel[i][j][1]*Cryst->symrel[i][j][1]+Cryst->symrel[i][j][2]*Cryst->symrel[i][j][2];
          //norm_factor=sqrt(row_sum);
          //Cryst->symrel[i][j][0]=Cryst->symrel[i][j][0]/norm_factor;
          //Cryst->symrel[i][j][1]=Cryst->symrel[i][j][1]/norm_factor;
          //Cryst->symrel[i][j][2]=Cryst->symrel[i][j][2]/norm_factor;
          //printf("Symrel #%d: [%lf] [%lf] [%lf]\n",i+1,Cryst->symrel[i][0][j],Cryst->symrel[i][1][j],Cryst->symrel[i][2][j]);
        }
      }
    }

    if(strcmp(dummy,"tnons")==0){
      for(i=0;i<Cryst->nsym;i++){
        fscanf(abinit_out,"%lf %lf %lf",&Cryst->trans_vec[i][0], &Cryst->trans_vec[i][1], &Cryst->trans_vec[i][2]);
        //printf("Translation vector accompanying symrel #%d: [%lf] [%lf] [%lf]\n",i+1,Cryst->trans_vec[i][0], Cryst->trans_vec[i][1], Cryst->trans_vec[i][2]);
      }
    }

    if(strcmp(dummy,"znucl")==0){
      for(i=0;i<Cryst->ntype;i++){
        fscanf(abinit_out,"%lf",&temp);
        //printf("Temp is: %lf \n",temp);
        Cryst->znucl[i]=(int) temp;
        //printf("Reading %d znucl. \n",i+1);
        //printf("Element %d number: %d\n",i+1, Cryst->znucl[i]);
      }
    }

    if(strcmp(dummy,"xred")==0){
      for(i=0;i<Cryst->natom;i++){
        fscanf(abinit_out,"%lf %lf %lf",&Cryst->xred[i],&Cryst->yred[i],&Cryst->zred[i]);
        Cryst->atom_coord[i][0]=Cryst->xred[i];
        Cryst->atom_coord[i][1]=Cryst->yred[i];
        Cryst->atom_coord[i][2]=Cryst->zred[i];
        while(Cryst->atom_coord[i][0]<0){
          Cryst->atom_coord[i][0]++;
        }
        while(Cryst->atom_coord[i][0]>=1){
          Cryst->atom_coord[i][0]--;
        }
        while(Cryst->atom_coord[i][1]<0){
          Cryst->atom_coord[i][1]++;
        }
        while(Cryst->atom_coord[i][1]>=1){
          Cryst->atom_coord[i][1]--;
        }
        while(Cryst->atom_coord[i][2]<0){
          Cryst->atom_coord[i][2]++;
        }
        while(Cryst->atom_coord[i][2]>=1){
          Cryst->atom_coord[i][2]--;
        }
        //printf("Atom %d coordinates: %lf %lf %lf\n",i+1, Cryst->xred[i],Cryst->yred[i],Cryst->zred[i]);
      }
    }

    if(strcmp(dummy,"typat")==0){
      for(i=0;i<Cryst->natom;i++){
        fscanf(abinit_out,"%d",&Cryst->typat[i]);
        //printf("Atom %d type: %d\n",i+1, Cryst->typat[i]);
      }
    }
  }
  //printf("End of reading abinit out. \n");
}

void read_ACF(struct Crystal * Cryst, FILE * ACF){
  char dummy[200];
  int line_count=0;
  int atom_count=0;
  char temp[20][20];
  int i;
  memset(dummy,0,200);
  //printf("In ACF file. \n");
  while(fgets(dummy,200,ACF) != NULL){
    //printf("Line %d is: %s \n",(line_count+1),dummy);
    line_count++;
    if(line_count>=Cryst->natom+3){
      //printf("Terminating line reached.\n");
      break;
    }
    if(line_count>=3){
      sscanf(dummy,"%s %s %s %s %s %s %s",temp[0],temp[1],temp[2],temp[3],temp[4],temp[5],temp[6]);
      //printf("Line scanning done. \n");
      sscanf(temp[4],"%lf",&Cryst->bader_e_atom[atom_count]);
      //printf("The charge on atom %d is: %lf \n",(atom_count+1),Cryst->bader_e_atom[atom_count]);
      atom_count++;       
    } 
  }
}

void calc_Bader_charge(struct Crystal * Cryst){
  int i=0,j=0;
  int atom_No=0;
  int site_count=0;
  double charge_sum=0.0;
  for(i=0;i<Cryst->nsite;i++){
    atom_No=Cryst->uneqv_atom[i];
    Cryst->site_element_No[i]=Cryst->typat[atom_No]-1;
    //printf("Site %d with multiplicity %d has element No: %d, has %d electrons in pseudopotential.\n",(i+1),Cryst->eqv_multi[i],Cryst->site_element_No[i],Cryst->psp_e[Cryst->site_element_No[i]]);
  }
  
  for(i=0;i<Cryst->natom;i++){
    Cryst->psp_e_atom[i]=(double) Cryst->psp_e[Cryst->typat[i]-1]; 
    //printf("Atom %d has %lf electrons in psp. \n",(i+1),Cryst->psp_e_atom[i]);
    Cryst->bader_charge[i]=Cryst->psp_e_atom[i]-Cryst->bader_e_atom[i];
    //printf("Charge on atom %d: %lf\n",(i+1),Cryst->bader_charge[i]);
  }
 
  for(i=0;i<Cryst->nsite;i++){
    for(j=0;j<Cryst->eqv_multi[i];j++){
      charge_sum=charge_sum+Cryst->bader_charge[Cryst->eqv_atom_No[i][j]];
    }
    Cryst->bader_e_site[i]=charge_sum/Cryst->eqv_multi[i];
    charge_sum=0.0;
    //printf("Average charge on site %d: %lf \n",(i+1),Cryst->bader_e_site[i]); 
  }

  for(i=0;i<Cryst->nsite;i++){
    Cryst->bader_100[i]=Cryst->bader_e_site[i]*1.00;
    Cryst->bader_75[i]=Cryst->bader_e_site[i]*0.75;
    Cryst->bader_50[i]=Cryst->bader_e_site[i]*0.50;
    Cryst->bader_25[i]=Cryst->bader_e_site[i]*0.25;
    printf("Charge on site %d - 100%: %lf; 75%: %lf; 50%: %lf; 25%: %lf\n", (i+1),Cryst->bader_100[i],Cryst->bader_75[i],Cryst->bader_50[i],Cryst->bader_25[i]);
  }
  
}

void write_ape(FILE * ape_file, struct Crystal * Cryst, int element_site_id, int site_No, int ionicity){
  int Z=0;
  char noble_gas_symbol[10][10];
  int i=0, j=0;
  double charge=0.0;

  if(ionicity==0){
    charge=Cryst->bader_100[site_No];
  }else if(ionicity==1){
    charge=Cryst->bader_75[site_No];
  }else if(ionicity==2){
    charge=Cryst->bader_50[site_No];
  }else if(ionicity==3){
    charge=Cryst->bader_25[site_No];
  }

  strcpy(noble_gas_symbol[1],"He");
  strcpy(noble_gas_symbol[2],"Ne");
  strcpy(noble_gas_symbol[3],"Ar");
  strcpy(noble_gas_symbol[4],"Kr");
  strcpy(noble_gas_symbol[5],"Xe");
  strcpy(noble_gas_symbol[6],"Rn");

  fprintf(ape_file,"CalculationMode = ae\n");
  fprintf(ape_file,"Units = 1\n");
  fprintf(ape_file,"Verbose = 30\n");
  fprintf(ape_file,"WaveEquation = scalar_rel\n");
  fprintf(ape_file,"SpinMode = unpolarized\n");
  fprintf(ape_file,"XCFunctional = lda_xc_teter93\n");
  fprintf(ape_file,"XCCorrections = none\n");
  fprintf(ape_file,"TheoryLevel = dft\n");
  fprintf(ape_file,"NuclearCharge = %d\n",Cryst->znucl[element_site_id]);
  fprintf(ape_file,"%%Orbitals\n");
  Z=Cryst->znucl[element_site_id];
  Cryst->noble_gas_id[site_No]=-1;

  if(Z<=2){
    Cryst->noble_gas_id[site_No]=0;
    Cryst->e_left[site_No]=Z-charge;
  }else if(Z<=10){
    Cryst->noble_gas_id[site_No]=1;
    Cryst->e_left[site_No]=Z-2-charge;
  }else if(Z<=18){
    Cryst->noble_gas_id[site_No]=2;
    Cryst->e_left[site_No]=Z-10-charge;
  }else if(Z<=36){
    Cryst->noble_gas_id[site_No]=3;
    Cryst->e_left[site_No]=Z-18-charge;
  }else if(Z<=54){
    Cryst->noble_gas_id[site_No]=4;
    Cryst->e_left[site_No]=Z-36-charge;
  }else if(Z<=86){
    Cryst->noble_gas_id[site_No]=5;
    Cryst->e_left[site_No]=Z-54-charge;
  }else if(Z<=118){
    Cryst->noble_gas_id[site_No]=6;
    Cryst->e_left[site_No]=Z-86-charge;
  }

  if(Cryst->noble_gas_id[site_No]!=0){
    fprintf(ape_file," \"%s\"\n",noble_gas_symbol[Cryst->noble_gas_id[site_No]]);
  }

  //print("Value of i is: %d",i);

  if(Z<=2){
    if(Cryst->e_left[site_No]>2.0){
      fprintf(ape_file,"%%There are %lf electrons to be filled, but only 1s orbital should be used. Put 2 elctrons in 1s and ignore the rest",Cryst->e_left[site_No]);
      Cryst->e1s=2.0;
    }else if(Cryst->e_left[site_No]<0){
      fprintf(ape_file,"%%There are less than 0 electrons. Setting 1s electron to 0.");
      Cryst->e1s=0.0;
    }
    Cryst->e1s=Cryst->e_left[site_No];
    fprintf(ape_file," 1 | 0 | %lf\n",Cryst->e1s);
  }else if(Z<=10){
    Cryst->e2s=min(2.0,Cryst->e_left[site_No]);
    Cryst->e2p=max(0.0,Cryst->e_left[site_No]-2);
    fprintf(ape_file," 2 | 0 | %lf\n",Cryst->e2s);
    if (Cryst->e2p != 0.0) fprintf(ape_file," 2 | 1 | %lf\n",Cryst->e2p);
  }else if(Z<=18){
    Cryst->e3s=min(2.0,Cryst->e_left[site_No]);
    Cryst->e3p=max(0.0,Cryst->e_left[site_No]-2);
    fprintf(ape_file," 3 | 0 | %lf\n",Cryst->e3s);
    if (Cryst->e3p != 0.0) fprintf(ape_file," 3 | 1 | %lf\n",Cryst->e3p);
  }else if(Z<=20){
    Cryst->e4s=min(2.0,Cryst->e_left[site_No]);
    Cryst->e3d=max(0.0,Cryst->e_left[site_No]-2);
    Cryst->e4p=max(0.0,Cryst->e_left[site_No]-12);
    fprintf(ape_file," 4 | 0 | %lf\n",Cryst->e4s);
    if (Cryst->e3d != 0.0) fprintf(ape_file," 3 | 2 | %lf\n",Cryst->e3d);
    if (Cryst->e4p != 0.0) fprintf(ape_file," 4 | 1 | %lf\n",Cryst->e4p);
  }else if(Z<=30 && Z!=24 && Z!=29){
    if(charge>=0){
      Cryst->e4s=max(0.0,2.0-charge);
      Cryst->e3d=Z-20.0-max(0.0,charge-2.0);
    }else if(charge<0){
      Cryst->e4s=2.0;
      Cryst->e3d=min(10.0,Cryst->e_left[site_No]-2.0);
      Cryst->e4p=max(0.0,Cryst->e_left[site_No]-12.0);
    }
    fprintf(ape_file," 4 | 0 | %lf\n",Cryst->e4s);
    fprintf(ape_file," 3 | 2 | %lf\n",Cryst->e3d);
    if (Cryst->e4p!=0) fprintf(ape_file," 4 | 1 | %lf\n",Cryst->e4p);
  }else if(Z==24 || Z==29){
    if(charge>=0){
      Cryst->e4s=max(0.0,1.0-charge);
      Cryst->e3d=Z-19.0-max(0.0,charge-1.0);
    }else if(charge<0){
      Cryst->e4s=min(2.0,1+max(0.0,Cryst->e_left[site_No]-11.0));
      Cryst->e3d=min(10.0,Cryst->e_left[site_No]-1.0);
      Cryst->e4p=max(0.0,Cryst->e_left[site_No]-12.0);
    }
    fprintf(ape_file," 4 | 0 | %lf\n",Cryst->e4s);
    fprintf(ape_file," 3 | 2 | %lf\n",Cryst->e3d);
    if (Cryst->e4p!=0) fprintf(ape_file," 4 | 1 | %lf\n",Cryst->e4p);
  }else if(Z<=36){
    Cryst->e4s=2.0-max(0.0,12.0-Cryst->e_left[site_No]);
    Cryst->e3d=10.0;
    Cryst->e4p=max(0.0,Cryst->e_left[site_No]-12.0);
    fprintf(ape_file," 4 | 0 | %lf\n",Cryst->e4s);
    fprintf(ape_file," 3 | 2 | %lf\n",Cryst->e3d);
    fprintf(ape_file," 4 | 1 | %lf\n",Cryst->e4p);    
  }else if(Z<=38){
    Cryst->e5s=min(2.0,Cryst->e_left[site_No]);
    Cryst->e4d=max(0.0,Cryst->e_left[site_No]-2);
    Cryst->e5p=max(0.0,Cryst->e_left[site_No]-12);
    fprintf(ape_file," 5 | 0 | %lf\n",Cryst->e5s);
    if (Cryst->e4d != 0.0) fprintf(ape_file," 4 | 2 | %lf\n",Cryst->e4d);
    if (Cryst->e5p != 0.0) fprintf(ape_file," 5 | 1 | %lf\n",Cryst->e5p);
  }else if(Z==39 || Z==40 || Z==43 || Z==48){
    if(charge>=0){
      Cryst->e5s=max(0.0,2.0-charge);
      Cryst->e4d=Z-38-max(0.0,charge-2.0);
    }else if(charge<0){
      Cryst->e5s=2.0;
      Cryst->e4d=min(10.0,Cryst->e_left[site_No]-2.0);
      Cryst->e5p=max(0.0,Cryst->e_left[site_No]-12.0);
    }
    fprintf(ape_file," 5 | 0 | %lf\n",Cryst->e5s);
    fprintf(ape_file," 4 | 2 | %lf\n",Cryst->e4d);
    if (Cryst->e5p != 0.0) fprintf(ape_file," 5 | 1 | %lf\n",Cryst->e5p);    
  }else if(Z==41 || Z==42 || Z==44 || Z==45 || Z==47){
    if(charge>=0){
      Cryst->e5s=max(0.0,1.0-charge);
      Cryst->e4d=Z-36-max(0.0,charge-1.0);
    }else if(charge<0){
      Cryst->e5s=min(2.0,1+max(0.0,Cryst->e_left[site_No]-11.0));;
      Cryst->e4d=min(10.0,Cryst->e_left[site_No]-1.0);
      Cryst->e5p=max(0.0,Cryst->e_left[site_No]-12.0);
    }
    fprintf(ape_file," 5 | 0 | %lf\n",Cryst->e5s);
    fprintf(ape_file," 4 | 2 | %lf\n",Cryst->e4d);
    if (Cryst->e5p != 0.0) fprintf(ape_file," 5 | 1 | %lf\n",Cryst->e5p);
  }else if(Z==46){
    if(charge>=0){
      Cryst->e5s=0;
      Cryst->e4d=10-charge;
      Cryst->e5p=0;
    }else if(charge<0){
      Cryst->e5s=min(2,-charge);
      Cryst->e4d=10;
      Cryst->e5p=max(0.0,-charge-2);
    }
    fprintf(ape_file," 5 | 0 | %lf\n",Cryst->e5s);
    fprintf(ape_file," 4 | 2 | %lf\n",Cryst->e4d);
    fprintf(ape_file," 5 | 1 | %lf\n",Cryst->e5p);
  }else if(Z<=54){
    Cryst->e5s=2.0-max(0.0,12.0-Cryst->e_left[site_No]);
    Cryst->e4d=10.0;
    Cryst->e5p=max(0.0,Cryst->e_left[site_No]-12.0);
    fprintf(ape_file," 5 | 0 | %lf\n",Cryst->e5s);
    fprintf(ape_file," 4 | 2 | %lf\n",Cryst->e4d);
    fprintf(ape_file," 5 | 1 | %lf\n",Cryst->e5p);
  }else if(Z<=56){
    Cryst->e6s=min(2.0,Cryst->e_left[site_No]);
    Cryst->e5d=max(0.0,min(1.0,Cryst->e_left[site_No]-2));
    Cryst->e4f=max(0.0,Cryst->e_left[site_No]-3);
    fprintf(ape_file," 6 | 0 | %lf\n",Cryst->e6s);
    if(Cryst->e5d!=0) fprintf(ape_file," 5 | 2 | %lf\n",Cryst->e5d);
    if(Cryst->e4f!=0) fprintf(ape_file," 4 | 3 | %lf\n",Cryst->e4f);
  }else if(Z<=70 && Z!=57 && Z!=58 && Z!=64){
    if(charge>0){
      Cryst->e6s=min(0.0,2.0-charge);
      Cryst->e4f=Z-56-max(0.0,charge-2.0);
      Cryst->e5d=0;
    }else if(charge<0){
      Cryst->e6s=2.0;
      Cryst->e4f=min(14.0,Cryst->e_left[site_No]-2.0);
      Cryst->e5d=max(0.0,Cryst->e_left[site_No]-16.0);
    }
    fprintf(ape_file," 6 | 0 | %lf\n",Cryst->e6s);
    if(Cryst->e5d!=0) fprintf(ape_file," 5 | 2 | %lf\n",Cryst->e5d);
    if(Cryst->e4f!=0) fprintf(ape_file," 4 | 3 | %lf\n",Cryst->e4f);
  }else if(Z==57 || Z==58 || Z==64){
    if(charge>0){
      Cryst->e6s=max(0.0,2.0-charge);
      Cryst->e5d=max(0.0,1-max(0.0,charge-2.0));
      Cryst->e4f=Z-57-max(0.0,charge-3.0);
    }else if(charge<0){
      Cryst->e6s=2.0;
      Cryst->e4f=min(14.0,Cryst->e_left[site_No]-3.0);
      Cryst->e5d=1.0+max(0.0,Cryst->e_left[site_No]-16.0);
    }
    fprintf(ape_file," 6 | 0 | %lf\n",Cryst->e6s);
    if(Cryst->e5d!=0) fprintf(ape_file," 5 | 2 | %lf\n",Cryst->e5d);
    if(Cryst->e4f!=0) fprintf(ape_file," 4 | 3 | %lf\n",Cryst->e4f);
  }else if(Z<=80 && Z!=78 && Z!=79){
    if(charge>0){
      Cryst->e6s=min(0.0,2.0-charge);
      Cryst->e5d=min(0.0,Z-70.0-max(charge-2.0,0.0));
      Cryst->e4f=14.0-(charge-2.0-(Z-70));
    }else if(charge<0){
      Cryst->e6s=2.0;
      Cryst->e5d=min(10.0,Z-70.0-charge);
      Cryst->e4f=14.0;
      Cryst->e6p=max(0.0,Z-70.0-charge-10.0);
    }
    fprintf(ape_file," 6 | 0 | %lf\n",Cryst->e6s);
    if(Cryst->e5d!=0) fprintf(ape_file," 5 | 2 | %lf\n",Cryst->e5d);
    if(Cryst->e4f!=0) fprintf(ape_file," 4 | 3 | %lf\n",Cryst->e4f);
    if(Cryst->e6p!=0) fprintf(ape_file," 6 | 1 | %lf\n",Cryst->e6p);
  }else if(Z==78 || Z==79){
    if(charge>0){
      Cryst->e6s=min(0.0,1.0-charge);
      Cryst->e5d=min(0.0,Z-69.0-max(charge-1.0,0.0));
      Cryst->e4f=14.0;
    }else if(charge<0){
      Cryst->e6s=min(2.0,1.0+max(0.0,Z-69.0-charge-10.0));
      Cryst->e5d=min(10.0,Z-69.0-charge);
      Cryst->e4f=14.0;
      Cryst->e6p=max(0.0,Z-70-charge-10);
    }
    fprintf(ape_file," 6 | 0 | %lf\n",Cryst->e6s);
    if(Cryst->e5d!=0) fprintf(ape_file," 5 | 2 | %lf\n",Cryst->e5d);
    if(Cryst->e4f!=0) fprintf(ape_file," 4 | 3 | %lf\n",Cryst->e4f);
    if(Cryst->e6p!=0) fprintf(ape_file," 6 | 1 | %lf\n",Cryst->e6p);
  }else if(Z<=86){
    Cryst->e6p=max(0.0,Cryst->e_left[site_No]-26.0);
    Cryst->e6s=2.0-max(0.0,26.0-Cryst->e_left[site_No]);
    Cryst->e5d=10.0-max(0.0,24.0-Cryst->e_left[site_No]);
    Cryst->e4f=14;
    fprintf(ape_file," 6 | 0 | %lf\n",Cryst->e6s);
    if(Cryst->e5d!=0) fprintf(ape_file," 5 | 2 | %lf\n",Cryst->e5d);
    if(Cryst->e4f!=0) fprintf(ape_file," 4 | 3 | %lf\n",Cryst->e4f);
    if(Cryst->e6p!=0) fprintf(ape_file," 6 | 1 | %lf\n",Cryst->e6p);
  }else{
    fprintf(ape_file,"Your atom has Z > 86, it is currently not supported by prepareApe.");
  }
  
  fprintf(ape_file,"%%");
}

void generate_inp(struct Crystal * Cryst){
  int i=0, j=0;
  char directory_path[100];
  char inp_path[100];
  char Ionicity[4][10];
  char Ionicity_path[100];
  int element_site_id=0;
  int site_count[nsite_max];
  FILE * current_ape;
  FILE * Bader_Info;
  mkdir("Bader",0755);
  for(i=0;i<nsite_max;i++){
    site_count[i]=1;
  }

  Bader_Info = fopen("Bader_Info","w");
  fprintf(Bader_Info,"%d\n",Cryst->nsite);
  strcpy(Ionicity[0],"I100");
  strcpy(Ionicity[1],"I75");
  strcpy(Ionicity[2],"I50");
  strcpy(Ionicity[3],"I25");
  for(i=0;i<4;i++){
    sprintf(Ionicity_path,"Bader/%s",Ionicity[i]);
    mkdir(Ionicity_path,0755);
    for(j=0;j<Cryst->nsite;j++){
      element_site_id=Cryst->typat[Cryst->uneqv_atom[j]]-1;
      sprintf(directory_path,"Bader/%s/%s%d",Ionicity[i],Cryst->element_name[element_site_id],site_count[element_site_id]);
      fprintf(Bader_Info,"%s%d\n",Cryst->element_name[element_site_id],site_count[element_site_id]);
      site_count[element_site_id]++;
      mkdir(directory_path,0755);
      sprintf(inp_path,"%s/inp.ape",directory_path);
      //printf("Dirctory path is: %s \n",directory_path);
      //printf("inp.ape file path: %s \n",inp_path);
      current_ape=fopen(inp_path,"w");
      write_ape(current_ape,Cryst,element_site_id,j,i);
      fclose(current_ape);
    }
    for(j=0;j<nsite_max;j++){
      site_count[j]=1;
    }
  }
}

int main (int argc, char * argv[]){
  char base_name[100];
  char out_file_name[100];
  char files_file_name[100];
  char ACF_name[100];
  FILE * ACF;
  FILE * abinit_out;
  FILE * files_file;

  if(argc==3 || argc==4){
    strcpy(files_file_name,argv[1]);
    strcpy(out_file_name,argv[2]);
    if(argc==4) strcpy(ACF_name,argv[3]);
    if(argc==3) strcpy(ACF_name,"ACF.dat");
  }
  else{
    printf("Usage: <abinit files file name> <ABINIT output file name> optionally: <Bader output file name>\n");
    exit(0);
  }

  if(access(files_file_name,F_OK) != -1){
    printf("Name of the files file is: %s \n", files_file_name);
  }
  else{
    printf("The files file: %s does not exist. Please double check!!\n", files_file_name);
    exit(0);
  }

  if(access(out_file_name,F_OK) != -1){
    printf("Name of the ABINIT output file is: %s \n", out_file_name);
  }
  else{
    printf("The ABINIT output file: %s does not exist. Please double check!!\n", out_file_name);
    exit(0);
  }

  if(access(ACF_name,F_OK) != -1){
    printf("Name of the Bader output file is: %s \n", ACF_name);
  }
  else{
    printf("The Bader output file: %s does not exist. Please make sure you have ran the Bader calculaiton.\n",ACF_name);
    exit(0);
  }

  printf("Begin reading ABINIT output file.\n");
  abinit_out=fopen(out_file_name,"r");
  read_output_file(&Cryst,abinit_out);
  fclose(abinit_out);
  printf("Done reading ABINIT output file.\n");

  printf("Begin reading ABINIT files file.\n");
  files_file=fopen(files_file_name,"r");
  read_files_file(&Cryst,files_file);
  fclose(files_file);
  printf("Done reading ABINIT files file.\n");
  
  //printf("Proceed to find symmetry equivlent atoms. \n");  
  find_sym_atom(&Cryst);
  //printf("Done finding symmetry in crystal structure. \n");

  //printf("Begin reading Bader profile. \n");
  ACF=fopen(ACF_name,"r");
  read_ACF(&Cryst,ACF);
  fclose(ACF);
  //printf("Done reading Bader profile. \n");

  //printf("Start to calculate Bader charge. \n");
  calc_Bader_charge(&Cryst);
  //printf("Bader charges have been calculated. \n");

  printf("Begin to write inp.ape files. \n");
  generate_inp(&Cryst);
  //printf("Done writing inp.ape files .\n");

  printf("Done!\n");

}
