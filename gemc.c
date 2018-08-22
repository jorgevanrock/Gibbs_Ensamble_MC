#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "GEN_RAN2/ran2.h"
#include "GEN_RAN2/Random.h"

#define SI 1
#define NO 0

//**************************************************** STRUCTS *
struct box{double x;double y;double z;};
struct sim{int nat;double dens;double upot;struct box box;};
struct sys{int nat;double dens;double mcStep;double rcut;double temp;double dr;int dim;char potential[5];struct sim sim1,sim2;struct box box;};
struct lj{double sig;double eps;};
struct vector{double x;double y;double z;};
typedef struct {int esp;struct vector pos;} atomo;

//**************************************************** PROTOTIPOS *
void compBoxAux(int dim,double dens,int nat,double box[]);
void compBox(struct sys *sys);
void readData(char inFile[],struct sys *sys,struct lj *lj);
void minima(int dim,struct box box,double *dx,double *dy,double *dz);
int traslape(int dim,struct sim sim,atomo atom[],int atoms,double xo,double yo,double zo);
void makeAtoms(int dim,struct sim sim,atomo atom[]);
void printXYZ(char outxyz[],int dim,struct sim sim,atomo atom[]);
double potencial(struct sys sys,struct sim sim,atomo atom[],int i,int j,struct lj lj);
double energia_total(struct sys sys,struct sim sim,atomo atom[],struct lj lj);
double ener_atom(int o,struct sys sys,struct sim sim,atomo atom[],struct lj lj);
void mcmove(struct sys sys,atomo atom[],struct sim *sim,struct lj lj);
void gibbs(struct sys sys,atomo atom1[],atomo atom2[],struct lj lj);

//######### MAIN GEMC ###############
void main(int argc, char *argv[]){
  char inFile[15]="";
  struct sys sys;
  struct lj lj;
  
  if(argc < 2 || argc >2){printf("Syntax: ./exe  <run.file>\n");exit(0);}
  strcpy(inFile,argv[1]);

  readData(inFile, &sys, &lj);
  atomo atom1[sys.sim1.nat], atom2[sys.sim2.nat];
  makeAtoms(sys.dim,sys.sim1,atom1);
  makeAtoms(sys.dim,sys.sim2,atom2);
  printXYZ("conf1.xyz",sys.dim,sys.sim1,atom1);
  printXYZ("conf2.xyz",sys.dim,sys.sim2,atom2);
  
  printf("energía total Sistema 1: %lf\n",energia_total(sys,sys.sim1,atom1,lj));
  printf("energía total Sistema 2: %lf\n",energia_total(sys,sys.sim2,atom2,lj));

  //call gibbs
  gibbs(sys,atom1,atom2,lj);
  printXYZ("conf1.2.xyz",sys.dim,sys.sim1,atom1);


}
//###################################
void compBoxAux(int dim,double dens,int nat,double box[]){

  if(dim == 3){
    box[0] = pow((double)nat/dens,1.0f/3.0f);
    box[1] = box[0];
    box[2] = box[0];
  }
  else{
    box[0] = pow((double)nat/dens,1.0f/2.0f);
    box[1] = box[0];
  }
}
//******************************************************
void compBox(struct sys *sys){
  double box[(*sys).dim];

  compBoxAux((*sys).dim, (*sys).sim1.dens, (*sys).sim1.nat, box);
  (*sys).sim1.box.x = box[0];
  (*sys).sim1.box.y = box[1];
  if((*sys).dim == 3)   (*sys).sim1.box.z = box[2];
 
  compBoxAux((*sys).dim, (*sys).sim2.dens, (*sys).sim2.nat, box);
  (*sys).sim2.box.x = box[0];
  (*sys).sim2.box.y = box[1];
  if((*sys).dim == 3)   (*sys).sim2.box.z = box[2];

}
//****************************************************** READ DATA *
void readData(char inFile[],struct sys *sys,struct lj *lj){
  FILE *f;
  char name[20];
  
  printf("%s\n",inFile);

  f = fopen(inFile,"r");
    fscanf(f,"%i\t%s\n",&(sys->dim),name);
    fscanf(f,"%i\t%s\n",&(sys->sim1.nat),name);
    fscanf(f,"%lf\t%s\n",&(sys->sim1.dens),name);
    fscanf(f,"%i\t%s\n",&(sys->sim2.nat),name);
    fscanf(f,"%lf\t%s\n",&(sys->sim2.dens),name);
    fscanf(f,"%lf\t%s\n",&(sys->rcut),name);
    fscanf(f,"%s\t%s\n",sys->potential,name);
    fscanf(f,"%lf\t%lf\t%s\n",&(lj->sig),&(lj->eps),name);
    fscanf(f,"%lf\t%s\n",&(sys->dr),name);
    fscanf(f,"%lf\t%s\n",&(sys->temp),name);
    fscanf(f,"%lf\t%s\n",&(sys->mcStep),name);
  fclose(f);

  //Num atomos total
  (*sys).nat = (*sys).sim1.nat + (*sys).sim2.nat;

  compBox(sys);

  printf("nat1: %i  nat2: %i\n",(*sys).sim1.nat,(*sys).sim2.nat);
  printf("rho1: %lf rho2: %lf\n",(*sys).sim1.dens,(*sys).sim2.dens);
  printf("lx1: %lf ly1: %lf  lz1: %lf\n",(*sys).sim1.box.x,(*sys).sim1.box.y,(*sys).sim1.box.z);
  printf("lx2: %lf ly2: %lf  lz2: %lf\n",(*sys).sim2.box.x,(*sys).sim2.box.y,(*sys).sim2.box.z);
  printf("potencial: %s\n",(*sys).potential);
  printf("epsilon: %lf  sigma: %lf\n",(*lj).eps,(*lj).sig);
}
//*************************************************** MINIMA IMAGEN *
void minima(int dim,struct box box,double *dx,double *dy,double *dz){
  if (*dx> 0.5*box.x )        *dx -= box.x;
  else if (*dx<-0.5*box.x)    *dx += box.x;
  if (*dy>0.5*box.y)          *dy -= box.y;
  else if (*dy<-0.5*box.y)    *dy += box.y;
  if(dim == 3){
    if (*dz>0.5*box.z)        *dz -= box.z;
    else if (*dz<-0.5*box.z)  *dz += box.z;
  }
}

//******************************************************** TRASLAPE *
int traslape(int dim,struct sim sim,atomo atom[],int atoms,double xo,double yo,double zo){
  double dx,dy,dz,rio;
  int i;

  for(i=0; i<atoms; i++){
    dx = atom[i].pos.x - xo;
    dy = atom[i].pos.y - yo;
    if(dim == 3)   dz = atom[i].pos.z - zo;
    else               dz = 0.0f;
    minima(dim,sim.box,&dx,&dy,&dz);
    rio = sqrt(dx*dx + dy*dy + dz*dz);
    if(rio < 1.00f){
      return SI;
    }
  }
  return NO;
}

//**************************************************** MAKE ATOMS *
void makeAtoms(int dim,struct sim sim,atomo atom[]){
  int i, nparti=0;
  double xo,yo,zo;

  while(nparti < sim.nat){
    xo = Random() * sim.box.x;
    yo = Random() * sim.box.y;
    if(dim == 3)   zo = Random() * sim.box.z;
    else               zo = 0.0f;
    if(traslape(dim,sim,atom,nparti,xo,yo,zo) == NO){
      atom[nparti].pos.x = xo;
      atom[nparti].pos.y = yo;
      if(dim == 3)   atom[nparti].pos.z = zo;
      atom[nparti].esp = 1;
      nparti++;
    }
  }
}

//****************************************************** PRINT OUT *
void printXYZ(char outxyz[],int dim,struct sim sim,atomo atom[]){
  int i;
  FILE *f;

  f = fopen(outxyz,"w");
  fprintf(f,"%d\n\n",sim.nat);
  for(i=0;i<sim.nat;i++){
    fprintf(f,"%i %lf %lf ",atom[i].esp, atom[i].pos.x, atom[i].pos.y);
    if(dim == 3)   fprintf(f,"%lf\n",atom[i].pos.z);
    else           fprintf(f,"\n");
  }
}
//********************************************************* POTENCIAL *
double potencial(struct sys sys,struct sim sim,atomo atom[],int i,int j,struct lj lj){
  double dx,dy,dz,sig_rij6,u=0.0,rij;

  dx = atom[i].pos.x - atom[j].pos.x;
  dy = atom[i].pos.y - atom[j].pos.y;
  if(sys.dim == 3)   dz = atom[i].pos.z - atom[j].pos.z;
  else               dz = 0.0f;
  minima(sys.dim,sim.box,&dx,&dy,&dz);
  rij = sqrt(dx*dx + dy*dy + dz*dz);

  if(rij <= sys.rcut){
    if(strcmp(sys.potential,"lj") == 0){
      sig_rij6 = pow((lj.sig/rij),6.0f);
      u = 4.0f * lj.eps * sig_rij6 * (sig_rij6 - 1.0f);
    }
  }
  return u;
}
//****************************************************** ENERGÍA TOTAL *
double energia_total(struct sys sys,struct sim sim,atomo atom[],struct lj lj){
  int i,j;
  double uij,upot=0.0;

  for(i=0;i<sim.nat-1;i++){
    for(j=i+1;j<sim.nat;j++){
      uij = potencial(sys,sim,atom,i,j,lj);
      upot += uij;
    }
  }
  return upot/(double)sim.nat;
}
//****************************************************** ENERGÍA POR ÁTOMO *
double ener_atom(int o,struct sys sys,struct sim sim,atomo atom[],struct lj lj){
  double uoj, uo=0.0;
  int j;

  for(j=0;j<sim.nat;j++){
    if(j != o){
      uoj = potencial(sys,sim,atom,o,j,lj);
      uo += uoj;
    }
  }
  return uo;
}
//**************************************************
int volado(){
  double aux;

  aux = Random();
  if(aux < 0.5f)   return 1;
  else             return 2;
}
//**************************************************
void mcmove(struct sys sys,atomo atom[],struct sim *sim,struct lj lj){
  int o;
  double xold,yold,zold,ener_old,ener_new,dE,beta;

  o = rand()%(*sim).nat; 

  xold = atom[o].pos.x;
  yold = atom[o].pos.y;
  if(sys.dim == 3)   zold = atom[o].pos.z;

  ener_old = ener_atom(o,sys,*sim,atom,lj);

  atom[o].pos.x += (2.0f*Random() - 1.0f)*sys.dr; 
  atom[o].pos.y += (2.0f*Random() - 1.0f)*sys.dr; 
  if(sys.dim == 3)   atom[o].pos.z += (2.0f*Random() - 1.0f)*sys.dr; 

  ener_new = ener_atom(o,sys,*sim,atom,lj);
  dE = ener_new - ener_old;
  
  //metropolis
  beta = sys.temp;
  if( dE < 0.0f  ||  Random() < exp(-beta*dE) ){
    (*sim).upot += dE;  
  }
  else{
    atom[o].pos.x = xold; 
    atom[o].pos.y = yold; 
    if(sys.dim == 3)   atom[o].pos.z = zold; 
  }
}
//**************************************************
void gibbs(struct sys sys,atomo atom1[],atomo atom2[],struct lj lj){
  int step;

  step = 0;
  while(step < sys.mcStep){
    if(volado() == 1)   mcmove(sys,atom1,&(sys.sim1),lj);   
    else                mcmove(sys,atom2,&(sys.sim2),lj);
    step++;
  }
}
