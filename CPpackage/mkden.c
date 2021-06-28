/*
mkden.c: A program that makes an APE profile from <element-orbital> files.

Last modified on June 28, 2019 by Amber Lim.
*/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#define PI 3.14159265
#define NTYPES_MAX 20
#define NPOINTMAX 12000
#define STRMAX 100

int Z = 0;
double charge = 0;

struct Denfile{
	double orbocc[20];
	double wf[NPOINTMAX][NTYPES_MAX],rad[NPOINTMAX],den[NPOINTMAX];
} denfile;

int FinishLine(FILE * fptr) {
	char check;
	int cont = 0;;
	while (cont ==0) {
		check = fgetc(fptr);
		if (check==10 || check==EOF) {
			cont = 1;
		}
	}
	if (check==EOF) {
		return 1;
	}
	else return 0;
}

int elemName(int Z, char * name) {
	char * element[119] = {
	"&","H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar",
    "K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",
    "Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe",
    "Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu",
    "Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn",
    "Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr",
    "Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og"};
	if (Z>=0 && Z<119) {
		strcpy(name, element[Z]);
	}
	else {
		strcpy(name,"N/A");
	}
	return 0;
}

int orbOcc(int Z, double charge) {
	int cumulativeOcc[19] = {2, 4, 10, 12, 18, 20, 30, 36, 38, 48, 54, 56, 70, 80, 86, 88, 102, 112, 118};
	int fullOcc[19] = {2, 2,  6,  2,  6,  2, 10,  6, 2, 10, 6, 2, 14, 10, 6, 2, 14, 10, 6};
	double nElec = (double)Z - charge;
	double nElecRemaining = nElec;
	int i;
	for (i=0; i<19; i++) {
		if (nElecRemaining > fullOcc[i]) {
			denfile.orbocc[i] = (double)fullOcc[i];
			nElecRemaining = nElecRemaining-(double)fullOcc[i];
		}
		else if (nElecRemaining == 0) {
			break;
		}
		else {
			denfile.orbocc[i] = nElecRemaining;
			nElecRemaining = 0;
			break;
		}
	}
	return 0;
}

int Orbital(int ind, char * name) {
	const char * orbital[19] = {	
    	"1s","2s","2p","3s","3p","4s","3d","4p","5s","4d","5p","6s","4f","5d","6p","7s","5f","6d","7p"};
	if (ind>=0 && ind<19) {
		strcpy(name,orbital[ind]);
	}
	return 0;
}

int makedenfile(char * element) {
	char orbital[STRMAX], orb[STRMAX], name[STRMAX], name2[STRMAX];
	int i = 0, k = 0, max_k = 0;
	double temp = 0;
	double counter = 0;
	double nElect = (double)Z - charge;
	FILE * fptr;
	snprintf(name,STRMAX,"%s-",element);
	for (i=0; i<19; i++)
	{
		Orbital(i,orb);
		printf("Enter occupancy of orbital %s (%lf): ",orb,denfile.orbocc[i]);
		scanf("%lf",&denfile.orbocc[i]);
		counter = counter + denfile.orbocc[i];
		if (denfile.orbocc[i]>0) {
			snprintf(name2,STRMAX,"%s%s",name,orb);
			fptr = fopen(name2,"r");
			k = 0;
			while(FinishLine(fptr) ==0) {
				fscanf(fptr, "%lf %lf %lf", &denfile.rad[k],&denfile.wf[k][i],&temp);
				denfile.den[k] +=denfile.wf[k][i]*denfile.wf[k][i]*denfile.orbocc[i]/(4*PI);
				k++;
			}
			fclose(fptr);
			max_k = k;
		}
        if (nElect-counter <= 0.0) {break;}
	}
	return max_k = k;
}

int writedenfile(char * element, double charge, int max_k) {
	char outfilename[STRMAX];
	int i;
	FILE * fptr;
	snprintf(outfilename,STRMAX,"%s-%lf",element,charge);
	fptr = fopen(outfilename,"w");
	for (i=0; i<max_k-1; i++) {
		if (!(denfile.rad[i]==0 && denfile.den[i]==0)) {
		fprintf(fptr,"%lf %lf\n",denfile.rad[i],denfile.den[i]);
		}
	}
	fclose(fptr);
	return 0;
}

int main() {
	char elem[4];
	printf("What is the element number you would like to make a profile for? ");
	scanf("%d",&Z);
	elemName(Z,elem);
	printf("What is the charge on the element (%s)? ",elem);
	scanf("%lf",&charge);
	orbOcc(Z,charge);
	int max_k = 0;
	max_k = makedenfile(elem);
	writedenfile(elem,charge,max_k);
}
