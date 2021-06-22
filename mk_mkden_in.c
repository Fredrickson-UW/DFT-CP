#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#define PI 3.14159265358979323846264338328
#define ntype_max 20

struct Crystal{
  int ntype;
  int znucl[ntype_max];
  char element_name[ntype_max][10];
  char psp_name_full[ntype_max][50];
  char psp_name[ntype_max][50];
  int psp_e[ntype_max];
  int core_e[ntype_max];
  int ele_No_order[ntype_max];
} Cryst;

void write_mkin_i0(struct Crystal * Cryst, FILE * mkin_file, int current_ntype){
  fprintf(mkin_file,"%d\n",Cryst->znucl[current_ntype]);
  fprintf(mkin_file,"0\n");
  switch(Cryst->znucl[current_ntype]){
    case 1:
      fprintf(mkin_file,"%d\n",1);
      break;
    case 2:
      fprintf(mkin_file,"%d\n",2);
      break;
    case 3:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",1);
      break;
    case 4:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      break;
    case 5:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",1);
      break;
    case 6:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      break;
    case 7:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",3);
      break;
    case 8:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",4);
      break;
    case 9:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",5);
      break;
    case 10:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      break;
    case 11:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",1);
      break;
    case 12:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      break;
    case 13:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",1);
      break;
    case 14:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      break;
    case 15:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",3);
      break;
    case 16:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",4);
      break;
    case 17:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",5);
      break;
    case 18:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      break;
    case 19:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",1);
      break;
    case 20:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      break;
    case 21:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",1);
      break;
    case 22:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      break;
    case 23:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",3);
      break;
    case 24:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",1);
      fprintf(mkin_file,"%d\n",5);
      break;
    case 25:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",5);
      break;
    case 26:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      break;
    case 27:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",7);
      break;
    case 28:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",8);
      break;
    case 29:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",1);
      fprintf(mkin_file,"%d\n",10);
      break;
    case 30:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      break;
    case 31:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",1);
      break;
    case 32:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",2);
      break;
    case 33:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",3);
      break;
    case 34:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",4);
      break;
    case 35:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",5);
      break;
    case 36:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      break;
    case 37:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",1);
      break;
    case 38:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      break;
    case 39:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",1);
      break;
    case 40:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      break;
    case 41:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",1);
      fprintf(mkin_file,"%d\n",4);
      break;
    case 42:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",1);
      fprintf(mkin_file,"%d\n",5);
      break;
    case 43:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",5);
      break;
    case 44:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",1);
      fprintf(mkin_file,"%d\n",7);
      break;
    case 45:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",1);
      fprintf(mkin_file,"%d\n",8);
      break;
    case 46:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",0);
      fprintf(mkin_file,"%d\n",10);
      break;
    case 47:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",1);
      fprintf(mkin_file,"%d\n",10);
      break;
    case 48:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      break;
    case 49:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",1);
      break;
    case 50:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",2);
      break;
    case 51:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",3);
      break;
    case 52:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",4);
      break;
    case 53:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",5);
      break;
    case 54:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      break;
    case 55:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",1);
      break;
    case 56:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      break;
    case 57:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",0);
      fprintf(mkin_file,"%d\n",1);
      break;
    case 58:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",1);
      fprintf(mkin_file,"%d\n",1);
      break;
    case 59:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",3);
      break;
    case 60:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",4);
      break;
    case 61:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",5);
      break;
    case 62:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      break;
    case 63:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",7);
      break;
    case 64:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",7);
      fprintf(mkin_file,"%d\n",1);
      break;
    case 65:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",9);
      break;
    case 66:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      break;
    case 67:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",11);
      break;
    case 68:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",12);
      break;
    case 69:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",13);
      break;
    case 70:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",14);
      break;
    case 71:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",14);
      fprintf(mkin_file,"%d\n",1);
      break;
    case 72:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",14);
      fprintf(mkin_file,"%d\n",2);
      break;
    case 73:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",14);
      fprintf(mkin_file,"%d\n",3);
      break;
    case 74:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",14);
      fprintf(mkin_file,"%d\n",4);
      break;
    case 75:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",14);
      fprintf(mkin_file,"%d\n",5);
      break;
    case 76:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",14);
      fprintf(mkin_file,"%d\n",6);
      break;
    case 77:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",14);
      fprintf(mkin_file,"%d\n",7);
      break;
    case 78:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",1);
      fprintf(mkin_file,"%d\n",14);
      fprintf(mkin_file,"%d\n",9);
      break;
    case 79:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",1);
      fprintf(mkin_file,"%d\n",14);
      fprintf(mkin_file,"%d\n",10);
      break;
    case 80:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",14);
      fprintf(mkin_file,"%d\n",10);
      break;
    case 81:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",14);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",1);
      break;
    case 82:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",14);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",2);
      break;
    case 83:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",14);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",3);
      break;
    case 84:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",14);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",4);
      break;
    case 85:
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",6);
      fprintf(mkin_file,"%d\n",2);
      fprintf(mkin_file,"%d\n",14);
      fprintf(mkin_file,"%d\n",10);
      fprintf(mkin_file,"%d\n",5);
      break;
    default:
      printf("One of the atom type in your crystal structure is not currently supported by this program.\n");
      break;
  }
}

void write_mkin_core(struct Crystal * Cryst, FILE * mkin_file, int current_ntype){
  fprintf(mkin_file,"%d\n",Cryst->znucl[current_ntype]);
  fprintf(mkin_file,"%d\n",Cryst->psp_e[current_ntype]);
  if(Cryst->core_e[current_ntype]==0){
    fprintf(mkin_file,"%d\n",0);
  }else if(Cryst->core_e[current_ntype]==2){
    fprintf(mkin_file,"%d\n",2);
  }else if(Cryst->core_e[current_ntype]==10){
    fprintf(mkin_file,"%d\n",2);
    fprintf(mkin_file,"%d\n",2);
    fprintf(mkin_file,"%d\n",6);
  }else if(Cryst->core_e[current_ntype]==18){
    fprintf(mkin_file,"%d\n",2);
    fprintf(mkin_file,"%d\n",2);
    fprintf(mkin_file,"%d\n",6);
    fprintf(mkin_file,"%d\n",2);
    fprintf(mkin_file,"%d\n",6);
  }else if(Cryst->core_e[current_ntype]==28){
    fprintf(mkin_file,"%d\n",2);
    fprintf(mkin_file,"%d\n",2);
    fprintf(mkin_file,"%d\n",6);
    fprintf(mkin_file,"%d\n",2);
    fprintf(mkin_file,"%d\n",6);
    fprintf(mkin_file,"%d\n",0);
    fprintf(mkin_file,"%d\n",10);
  }else if(Cryst->core_e[current_ntype]==36){
    fprintf(mkin_file,"%d\n",2);
    fprintf(mkin_file,"%d\n",2);
    fprintf(mkin_file,"%d\n",6);
    fprintf(mkin_file,"%d\n",2);
    fprintf(mkin_file,"%d\n",6);
    fprintf(mkin_file,"%d\n",2);
    fprintf(mkin_file,"%d\n",10);
    fprintf(mkin_file,"%d\n",6);
  }else if(Cryst->core_e[current_ntype]==46){
    fprintf(mkin_file,"%d\n",2);
    fprintf(mkin_file,"%d\n",2);
    fprintf(mkin_file,"%d\n",6);
    fprintf(mkin_file,"%d\n",2);
    fprintf(mkin_file,"%d\n",6);
    fprintf(mkin_file,"%d\n",2);
    fprintf(mkin_file,"%d\n",10);
    fprintf(mkin_file,"%d\n",6);
    fprintf(mkin_file,"%d\n",0);
    fprintf(mkin_file,"%d\n",10);
  }else if(Cryst->core_e[current_ntype]==54){
    fprintf(mkin_file,"%d\n",2);
    fprintf(mkin_file,"%d\n",2);
    fprintf(mkin_file,"%d\n",6);
    fprintf(mkin_file,"%d\n",2);
    fprintf(mkin_file,"%d\n",6);
    fprintf(mkin_file,"%d\n",2);
    fprintf(mkin_file,"%d\n",10);
    fprintf(mkin_file,"%d\n",6);
    fprintf(mkin_file,"%d\n",2);
    fprintf(mkin_file,"%d\n",10);
    fprintf(mkin_file,"%d\n",6);
  }else if(Cryst->core_e[current_ntype]==60){
    fprintf(mkin_file,"%d\n",2);
    fprintf(mkin_file,"%d\n",2);
    fprintf(mkin_file,"%d\n",6);
    fprintf(mkin_file,"%d\n",2);
    fprintf(mkin_file,"%d\n",6);
    fprintf(mkin_file,"%d\n",2);
    fprintf(mkin_file,"%d\n",10);
    fprintf(mkin_file,"%d\n",6);
    fprintf(mkin_file,"%d\n",0);
    fprintf(mkin_file,"%d\n",10);
    fprintf(mkin_file,"%d\n",0);
    fprintf(mkin_file,"%d\n",0);
    fprintf(mkin_file,"%d\n",14);
  }else if(Cryst->core_e[current_ntype]==68){
    fprintf(mkin_file,"%d\n",2);
    fprintf(mkin_file,"%d\n",2);
    fprintf(mkin_file,"%d\n",6);
    fprintf(mkin_file,"%d\n",2);
    fprintf(mkin_file,"%d\n",6);
    fprintf(mkin_file,"%d\n",2);
    fprintf(mkin_file,"%d\n",10);
    fprintf(mkin_file,"%d\n",6);
    fprintf(mkin_file,"%d\n",2);
    fprintf(mkin_file,"%d\n",10);
    fprintf(mkin_file,"%d\n",6);
    fprintf(mkin_file,"%d\n",0);
    fprintf(mkin_file,"%d\n",14);
  }else if(Cryst->core_e[current_ntype]==78){
    fprintf(mkin_file,"%d\n",2);
    fprintf(mkin_file,"%d\n",2);
    fprintf(mkin_file,"%d\n",6);
    fprintf(mkin_file,"%d\n",2);
    fprintf(mkin_file,"%d\n",6);
    fprintf(mkin_file,"%d\n",2);
    fprintf(mkin_file,"%d\n",10);
    fprintf(mkin_file,"%d\n",6);
    fprintf(mkin_file,"%d\n",2);
    fprintf(mkin_file,"%d\n",10);
    fprintf(mkin_file,"%d\n",6);
    fprintf(mkin_file,"%d\n",0);
    fprintf(mkin_file,"%d\n",14);
    fprintf(mkin_file,"%d\n",10);
  }
}

void write_addcore_in_file(struct Crystal * Cryst, FILE * addcore_in_file){
  int i=0, j=0;
  int max=0;
  int counted[ntype_max];
  int index=0;
  FILE * AddCore_Info_file;
  for(i=0;i<ntype_max;i++){
    counted[i]=0;
  }
  for(i=0;i<Cryst->ntype;i++){
    //printf("Finding the %d element. \n",(i+1));
    for(j=0;j<Cryst->ntype;j++){
      if(Cryst->znucl[j]>max && counted[j]==0){
        max=Cryst->znucl[j];
        index=j;
        //printf("max is: %d \n",max);
      }
    }
    Cryst->ele_No_order[Cryst->ntype-1-i]=index;
    counted[index]=1;
    max=0;
    index=0; 
  }

  //for(i=0;i<Cryst->ntype;i++){
    printf("Element with %d small No: %d \n",(i+1),Cryst->znucl[Cryst->ele_No_order[i]]);
  //}

  for(i=0;i<Cryst->ntype;i++){
    fprintf(addcore_in_file,"%s-%d.000000\n",Cryst->element_name[Cryst->ele_No_order[i]],Cryst->psp_e[Cryst->ele_No_order[i]]); 
  }

  AddCore_Info_file=fopen("AddCore_Info","w");
  fprintf(AddCore_Info_file,"%d\n",Cryst->ntype);
  for(i=0;i<Cryst->ntype;i++){
    fprintf(AddCore_Info_file,"%s\n",Cryst->element_name[i]);
  }
}

void read_files_file(struct Crystal * Cryst, FILE * files_file){
  char dummy[300];
  int line_count=0;
  int psp_count=0;
  int i=0,j=0,k=0;
  int psp_len;
  char temp1[20], temp2[20];
  int slash_pos=0;
  int dot1_pos=0;
  int temp1_len=0;
  char Z_str[20];
  int Z_no, e_no;
  int letter_count=0;
  char ele_name[20];
  memset(ele_name,0,20);
  memset(temp1,0,20);
  memset(temp2,0,20);
  memset(Z_str,0,20);
  memset(dummy,0,300);
  while( fgets (dummy,300,files_file) !=NULL ){
    line_count++;
    //printf("dummy: %s",dummy);
    if(line_count>=6){
      //printf("dummy: %s",dummy);
      strcpy(Cryst->psp_name_full[psp_count],dummy);
      psp_count++;
    }  
  }
  Cryst->ntype=line_count-5;
  printf("There are %d types of elements.\n",Cryst->ntype); 

  for(i=0;i<Cryst->ntype;i++){
    slash_pos=0;
    printf("element type: %d - psp name full: %s \n",(i+1),Cryst->psp_name_full[i]);
    psp_len=strlen(Cryst->psp_name_full[i]);
    for(j=0;j<psp_len;j++){
      if(Cryst->psp_name_full[i][j]=='/'){
        slash_pos=j;
      } 
    }
    printf("Slash pos is: %d \n",slash_pos);
    for(j=slash_pos+1;j<psp_len;j++){
      Cryst->psp_name[i][j-slash_pos-1]=Cryst->psp_name_full[i][j];
    }
    printf("element type: %d - psp name: %s \n",(i+1),Cryst->psp_name[i]);
  }

  for(i=0;i<Cryst->ntype;i++){
    letter_count=0;
    memset(temp1,0,20);
    memset(temp2,0,20);
    memset(ele_name,0,20);
    psp_len=strlen(Cryst->psp_name[i]);
    for(j=0;j<psp_len;j++){
      if(Cryst->psp_name[i][j]=='.'){
        dot1_pos=j;
        break;
      }else{
        temp1[j]=Cryst->psp_name[i][j];
      }
    }
    printf("Temp1 is: %s \n",temp1);
    for(j=dot1_pos+1;j<psp_len;j++){
      if(Cryst->psp_name[i][j]=='.'){
        break;
      }else{
        temp2[j-dot1_pos-1]=Cryst->psp_name[i][j];
      }
    }
    printf("Temp2 is: %s \n",temp2);
    temp1_len=strlen(temp1);
    for(j=0;j<temp1_len;j++){
      if(isdigit(temp1[j])){
        Z_str[j]=temp1[j];
      }
    }
    for(j=0;j<temp1_len;j++){
      if(isalpha(temp1[j])){
        ele_name[letter_count]=temp1[j];
        letter_count++;
      }
    }
    ele_name[0]=ele_name[0]-32;
    printf("Element name is: %s \n",ele_name);
    printf("Z is: %s \n",Z_str);
    sscanf(Z_str,"%d",&Z_no);
    sscanf(temp2,"%d",&e_no);
    printf("Z=%d, e=%d. \n",Z_no,e_no);
    Cryst->znucl[i]=Z_no;
    Cryst->psp_e[i]=e_no;
    Cryst->core_e[i]=Z_no-e_no;
    strcpy(Cryst->element_name[i],ele_name);
    printf("psp e: %d, core e: %d\n",Cryst->psp_e[i],Cryst->core_e[i]);
  }  
}

int main (int argc, char * argv[]){
  char base_name[100];
  char files_file_name[100];
  char mkin_name[100];
  char addcore_in_name[100];
  FILE * addcore_in_file;
  FILE * files_file;
  FILE * current_mkin_file;
  int i=0;
 
  if(argc==2){
    strcpy(files_file_name,argv[1]);
  }
  else{
    printf("Usage: <abinit files file name>\n");
    exit(0);
  }

  if(access(files_file_name,F_OK) != -1){
    printf("Name of the files file is: %s \n", files_file_name);
  }
  else{
    printf("The files file: %s does not exist. Please double check!!\n", files_file_name);
    exit(0);
  }

  printf("Begin reading ABINIT files file.\n");
  files_file=fopen(files_file_name,"r");
  read_files_file(&Cryst,files_file);
  fclose(files_file);
  printf("Done reading ABINIT files file.\n");
  
  printf("Begin writing mkden input file.\n");
  for(i=0;i<Cryst.ntype;i++){
    sprintf(mkin_name,"mkden_input%d",(2*i+1));
    current_mkin_file=fopen(mkin_name,"w");
    write_mkin_core(&Cryst,current_mkin_file,i);
    fclose(current_mkin_file);
    sprintf(mkin_name,"mkden_input%d",(2*i+2));
    current_mkin_file=fopen(mkin_name,"w");
    write_mkin_i0(&Cryst,current_mkin_file,i);
    fclose(current_mkin_file);
  }
  
  sprintf(addcore_in_name,"AddCore_input");
  addcore_in_file=fopen(addcore_in_name,"w");
  write_addcore_in_file(&Cryst,addcore_in_file);
  printf("Done!\n");
}
