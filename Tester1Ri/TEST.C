#include <conio.h>
#include <stdio.h>
#include <stdlib.h>
#include <process.h>

#define depart 10
#define taille_max 200
#define pas 20

#define iterations 50
#define range 0.6

float temps, temps_moy, temps_max,temps_min;
float tempsP, temps_moyP, temps_maxP,temps_minP;
int nbmin,nbmax;
float dev, dev_min, dev_max,dev_moy,nbmoy;
int UBref,LBref;
float tmpf;
int UB, LB, OPT, nb;

void generer(int nbtrav)
{
 int i;
 int ri,pi[taille_max],som=0;
 FILE *fichier;
 
 fichier=fopen("donnees.dat","wt");
 for (i=0;i<nbtrav;i++)
 {
  pi[i]=((float)rand()/(float)RAND_MAX)*99+1;
  som+=pi[i];
 }
 for (i=0;i<nbtrav;i++)
 {
  ri=((float)rand()/(float)RAND_MAX)*(float)50.5*nbtrav*(float)range;
  fprintf(fichier,"%d %d %d\n",i+1,ri,pi[i]);
 }
 fclose(fichier);
}

void main(void)
{
 int i,j;
 FILE *fichier;
 
 srand(time(NULL));
 fichier=fopen("test.res","wt");
 fprintf(fichier,"Comparison of the Beam search to the Optimal lower bound\n"); 
 fprintf(fichier,"-----------------------------------------------------\n\n"); 
 fprintf(fichier,"These tests are done on a PC Pentium II 366 Mhz\n"); 
 fprintf(fichier,"Range factor : %f\n",range); 
 fprintf(fichier,"Number of iterations for a fixed value of n : %d\n\n",iterations); 
 fprintf(fichier," n   dev_min(%) dev_moy(%) dev_max(%) tps_minUB(s) tps_moyUB(s) tps_maxUB(s) nodes_min nodes_moy nodes_max tps_minUB(s) tps_moyUB(s) tps_maxUB(s)\n");
 fprintf(fichier,"-------------------------------------------------------------------------------------------------------------------------------------------------\n"); 
 fclose(fichier);
 for (i=depart;i<taille_max;i+=pas)
 {
  printf("Nombre de travaux : %d\n",i);
  temps_moy=0;
  temps_min=9999999.9;
  temps_max=0;
  temps_moyP=0;
  temps_minP=9999999.9;
  temps_maxP=0;
  dev_moy=0;
  dev_min=99999999.9;
  dev_max=0;
  nbmin=99999999;
  nbmoy=nbmax=0;
  for (j=0;j<iterations;j++)
  {
   printf("Jeux n°%d\n",j+1);
   generer(i);
   spawnl(P_WAIT,"hds3.exe","hds3.exe",NULL);
   fichier=fopen("hds3.txt","rt");
   fscanf(fichier,"%f\n",&temps);
   fscanf(fichier,"%d\n",&UB);
   fscanf(fichier,"%d\n",&LB);
   fclose(fichier);
   if (LB>UB)
   {
    printf("Erreur : LB=%d et UB=%d\n",LB,UB);
    getch();
    exit(-1);
   }

   // On ecrit la borne supérieure
   fichier=fopen("initpse.txt","wt");
   fprintf(fichier,"%d\n",UB);
   fclose(fichier);

   spawnl(P_WAIT,"pse1.exe","pse1.exe",NULL);
   fichier=fopen("pse1.txt","rt");
   fscanf(fichier,"%f\n",&tempsP);
   fscanf(fichier,"%d\n",&OPT);
   fscanf(fichier,"%d\n",&LB);
   fscanf(fichier,"%d\n",&nb);
   fclose(fichier);
   if (OPT>UB)
   {
    printf("Erreur : OPT=%d et UB=%d\n",OPT,UB);
    getch();
    exit(-1);
   }

   if (temps<temps_min) temps_min=temps;
   if (temps>temps_max) temps_max=temps;
   temps_moy+=temps;
   if (tempsP<temps_minP) temps_minP=tempsP;
   if (tempsP>temps_maxP) temps_maxP=tempsP;
   temps_moyP+=tempsP;

   if (nb<nbmin) nbmin=nb;
   if (nb>nbmax) nbmax=nb;
   nbmoy+=nb;

   dev=(float)(UB-OPT)/(float)OPT;
   if (dev<dev_min) dev_min=dev;
   if (dev>dev_max) dev_max=dev;
   dev_moy+=dev;
  }
  temps_moy/=(float)iterations;
  temps_moyP/=(float)iterations;
  nbmoy/=(float)iterations;
  dev_moy=100.0*dev_moy/(float)iterations;
  dev_min=100.0*dev_min;
  dev_max=100.0*dev_max;
  fichier=fopen("test.res","at");
  fprintf(fichier,"%3d   %3.6f  %3.6f  %3.6f %3.6f %3.6f %3.6f %6d %6.2f %6d %3.6f %3.6f %3.6f\n",i,dev_min,dev_moy,dev_max,
        temps_min,temps_moy,temps_max,nbmin,nbmoy,nbmax,temps_minP,temps_moyP,temps_maxP);
  fclose(fichier);
 }
}
