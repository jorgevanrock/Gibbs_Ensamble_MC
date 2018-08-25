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
struct sim{int nat;double dens;double upot;double chemPot;struct box box;double dr;};
struct sys{int nat;double dens;double mcStep;int nAjusta;int nPrint;double accep;double rcut;double temp;double dv;int dim;char potential[5];struct sim sim1,sim2;struct box box;};
struct lj{double sig;double eps;};
struct vector{double x;double y;double z;};
typedef struct {struct vector pos;} atomo;
struct move{int accep; int intentos;};

//**************************************************** PROTOTIPOS *
void clear();
void compBoxAux(int dim,double dens,int nat,double box[]);
void compBox(struct sys *sys);
void readData(char inFile[],struct sys *sys,struct lj *lj);
void minima(int dim,struct box box,double *dx,double *dy,double *dz);
int traslape(int dim,struct sim sim,atomo atom[],int atoms,double xo,double yo,double zo);
void makeAtoms(int dim,struct sim sim,atomo atom[]);
void printXYZ(char outxyz[],int dim,struct sim sim,atomo atom[]);
void superPrint(int dim,int nat1,int nat2,atomo atom1[],atomo atom2[],double lx);
double potencial(struct sys sys,struct sim sim,atomo atom[],int i,int j,struct lj lj);
double energia_total(struct sys sys,struct sim sim,atomo atom[],struct lj lj);
double ener_atom(int o,struct sys sys,struct sim sim,atomo atom[],struct lj lj);
void pbc(double *x,double *y,double *z,struct box box);
void mcmove(struct sys sys,atomo atom[],struct sim *sim,struct lj lj,struct move *move);
void escala(double factor,struct sim *sim,atomo atom[],int dim);
void mcvol(struct sys *sys,atomo atom1[],atomo atom2[],struct lj lj,struct move *volu);
void creaParti(struct sys sys,struct sim *simA,struct sim *simB,atomo atomA[],atomo atomB[],struct lj lj);
void screen(int step,struct sys sys,int flag);
void ajustaDr(struct move *move,struct sys sys,struct sim *sim);
void ajustaDv(struct move *volu,struct sys *sys);
void gibbs(struct sys *sys,atomo atom1[],atomo atom2[],struct lj lj,struct move *move1,struct move *move2,struct move *volu);

//######### MAIN GEMC ###############
void main(int argc, char *argv[]){
  char inFile[15]="";
  struct sys sys;
  struct lj lj;
  struct move move1,move2,volu;
 
  clear(); 
  if(argc < 2 || argc >2){printf("Syntax: ./exe  <run.file>\n");exit(0);}
  strcpy(inFile,argv[1]);
  move1.accep = 0;   move1.intentos = 0;
  move2.accep = 0;   move2.intentos = 0;
  volu.accep  = 0;   volu.intentos  = 0;

  readData(inFile, &sys, &lj);
  atomo atom1[sys.sim1.nat + sys.sim2.nat], atom2[sys.sim1.nat + sys.sim2.nat];
  makeAtoms(sys.dim,sys.sim1,atom1);
  makeAtoms(sys.dim,sys.sim2,atom2);
  sys.sim1.upot = energia_total(sys,sys.sim1,atom1,lj);
  sys.sim2.upot = energia_total(sys,sys.sim2,atom2,lj);
  gibbs(&sys,atom1,atom2,lj,&move1,&move2,&volu);
}
//###################################
void clear(){
  system("rm -f *.xyz *txt");
}
//******************************************************
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
    fscanf(f,"%lf\t%s\n",&(sys->sim1.dr),name);
    fscanf(f,"%lf\t%s\n",&(sys->dv),name);
    fscanf(f,"%lf\t%s\n",&(sys->temp),name);
    fscanf(f,"%lf\t%s\n",&(sys->mcStep),name);
    fscanf(f,"%lf\t%s\n",&(sys->accep),name);
    fscanf(f,"%i\t%s\n",&(sys->nPrint),name);
    fscanf(f,"%i\t%s\n",&(sys->nAjusta),name);
  fclose(f);

  //Num atomos total
  (*sys).nat = (*sys).sim1.nat + (*sys).sim2.nat;
  (*sys).sim2.dr = (*sys).sim1.dr;
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
      atom[nparti].pos.z = zo;
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
    fprintf(f,"1 %lf %lf %lf\n",atom[i].pos.x, atom[i].pos.y, atom[i].pos.z);
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
  return upot;
}
//****************************************************** ENERGÍA POR ÁTOMO *
double ener_atom(int o,struct sys sys,struct sim sim,atomo atom[],struct lj lj){
  double uoj=0.0f, uo=0.0f;
  int j;

  for(j=0;j<sim.nat;j++){
    if(j != o){
      uoj = potencial(sys,sim,atom,o,j,lj);
      uo += uoj;
    }
  }
  return uo;
}
//************************************** SUPER PRINT *
void superPrint(int dim,int nat1,int nat2,atomo atom1[],atomo atom2[],double lx){
  int i,nat;
  double sepp = 8.0;
  FILE *f;

  nat = nat1 + nat2;

  f = fopen("super.xyz","a");
    fprintf(f,"%d\n\n",nat);
    for(i=0;i<nat1;i++){
      fprintf(f,"1 %lf %lf %lf\n",atom1[i].pos.x, atom1[i].pos.y, atom1[i].pos.z);
    }
    for(i=0;i<nat2;i++){
      fprintf(f,"2 %lf %lf %lf\n",atom2[i].pos.x+lx+sepp, atom2[i].pos.y, atom2[i].pos.z);
    }
  fclose(f);
}
//******************************************** PBC *
void pbc(double *x,double *y,double *z,struct box box){
  if(*x > box.x)         *x -= box.x;
  else if(*x <0.0f)      *x += box.x;
  if(*y > box.y)         *y -= box.y;
  else if(*y <0.0f)      *y += box.y;
  if(*z > box.z)         *z -= box.z;
  else if(*z <0.0f)      *z += box.z;
}
//******************************************* MC MOVE *
void mcmove(struct sys sys,atomo atom[],struct sim *sim,struct lj lj,struct move *move){
  int o;
  double xold,yold,zold,ener_old,ener_new,dE,beta;

int N = sim->nat,mueve=0;
while(mueve < N){
  o = rand()%(*sim).nat; 

  xold = atom[o].pos.x;
  yold = atom[o].pos.y;
  if(sys.dim == 3)   zold = atom[o].pos.z;

  ener_old = ener_atom(o,sys,*sim,atom,lj);

  atom[o].pos.x += (2.0f*Random() - 1.0f)*((*sim).dr); 
  atom[o].pos.y += (2.0f*Random() - 1.0f)*((*sim).dr); 
  if(sys.dim == 3)   atom[o].pos.z += (2.0f*Random() - 1.0f)*((*sim).dr); 

  pbc(&(atom[o].pos.x), &(atom[o].pos.y), &(atom[o].pos.z),(*sim).box );

  ener_new = ener_atom(o,sys,*sim,atom,lj);
  dE = ener_new - ener_old;
  //metropolis
  beta = 1.0f/sys.temp;
  if( dE < 0.0f  ||  Random() < exp(-beta*dE) ){
    (*sim).upot += dE;
    (*move).accep++;  
  }
  else{
    atom[o].pos.x = xold; 
    atom[o].pos.y = yold; 
    if(sys.dim == 3)   atom[o].pos.z = zold; 
  }
  (*move).intentos++;
mueve++;
}
}
//********************************************* AJUSTA DR *
void ajustaDr(struct move *move,struct sys sys,struct sim *sim){
  double minL,fracc,dro,fac;

  if((*move).intentos != 0){
    minL = (*sim).box.x;
  
    fracc = (double)((*move).accep) / (double)((*move).intentos);
    dro = (*sim).dr;
    (*sim).dr *= fabs(fracc*100.0f / (sys.accep));
    fac = (*sim).dr/dro;
    if(fac > 1.5f)                   (*sim).dr = dro*1.5f;
    else if(fac < 0.5f)              (*sim).dr = dro*0.5f;
    else if((*sim).dr > 0.5f*minL)   (*sim).dr = 0.5f*minL;
  
    (*move).accep = 0;
    (*move).intentos = 0;
  }
}
//********************************************* AJUSTA DV *
void ajustaDv(struct move *volu,struct sys *sys){
  double minL,fracc,dvo,fac;

  if((*volu).intentos != 0){
    if((*sys).sim1.box.x < (*sys).sim2.box.x)   minL = (*sys).sim1.box.x;
    else                                        minL = (*sys).sim2.box.x;
  
    fracc = (double)((*volu).accep) / (double)((*volu).intentos);
    dvo = (*sys).dv;
    (*sys).dv *= fabs(fracc*100.0f / (sys->accep));
    fac = (*sys).dv/dvo;
    if(fac > 1.5f)                   (*sys).dv = dvo*1.5f;
    else if(fac < 0.5f)              (*sys).dv = dvo*0.5f;
    else if((*sys).dv > 0.5f*minL)   (*sys).dv = 0.5f*minL;
  
    (*volu).accep = 0;
    (*volu).intentos = 0;
  }
}
//*************************************************** ESCALA *
void escala(double factor,struct sim *sim,atomo atom[],int dim){
  int i;

  (*sim).box.y = (*sim).box.x;
  if(dim == 3)   (*sim).box.z = (*sim).box.x;
  else           (*sim).box.z = 0.0f;

  for(i=0;i < (*sim).nat; i++){
     atom[i].pos.x *= factor;
     atom[i].pos.y *= factor;
     if(dim == 3)   atom[i].pos.z *= factor;
     else           atom[i].pos.z = 0.0f;
  }
}
//*************************************************** MC VOLUME *
void mcvol(struct sys *sys,atomo atom1[],atomo atom2[],struct lj lj,struct move *volu){
  double vol1,vol2,volT,lnVol,vol1new,vol2new,fac1,fac2;
  double utn1,utn2,dE1,dE2,arg1,arg2,arg,dv1,dv2;
  int dim = (*sys).dim;
  double beta = (*sys).temp;

  vol1 = (*sys).sim1.box.x * (*sys).sim1.box.y;
  if(dim == 3)   vol1 *= (*sys).sim1.box.z;

  vol2 = (*sys).sim2.box.x * (*sys).sim2.box.y;
  if(dim == 3)   vol2 *= (*sys).sim2.box.z;
  
  volT = vol1 + vol2;
  lnVol = log(vol1/vol2) + (2.0f*Random() - 1.0f)*(*sys).dv;

  vol1new = volT * exp(lnVol) / (exp(lnVol) + 1.0f);
  if(dim == 3)   (*sys).sim1.box.x = pow(vol1new,1.0f/3.0f);
  else                  (*sys).sim1.box.x = pow(vol1new,1.0f/2.0f);
  if(dim == 3)   fac1 = (*sys).sim1.box.x / pow(vol1,1.0f/3.0f);
  else                  fac1 = (*sys).sim1.box.x / pow(vol1,1.0f/2.0f);
  escala(fac1,&(*sys).sim1,atom1,dim);
  utn1 = energia_total(*sys,(*sys).sim1,atom1,lj);

  vol2new = volT - vol1new;
  if(dim == 3)   (*sys).sim2.box.x = pow(vol2new,1.0f/3.0f);
  else           (*sys).sim2.box.x = pow(vol2new,1.0f/2.0f);
  if(dim == 3)   fac2 = (*sys).sim2.box.x / pow(vol2,1.0f/3.0f);
  else           fac2 = (*sys).sim2.box.x / pow(vol2,1.0f/2.0f);
  escala(fac2,&(*sys).sim2,atom2,dim);
  utn2 = energia_total(*sys,(*sys).sim2,atom2,lj);
  
  //Criterio de aceptación
  dE1 = utn1 - (*sys).sim1.upot;
  dv1 = ((double)((*sys).sim1.nat + 1) * log(vol1new/vol1)) / beta;
  arg1 = dE1 - dv1;

  dE2 = utn2 - (*sys).sim2.upot;
  dv2 = ((double)((*sys).sim2.nat + 1) * log(vol2new/vol2)) / beta;
  arg2 = dE2 - dv2;

  arg = arg1 + arg2;
  if(Random() < exp(-arg*beta)){
    vol1 = vol1new;
    vol2 = vol2new;
    (*sys).sim1.upot = utn1;
    (*sys).sim2.upot = utn2;
    (*volu).accep++;
  }
  else{
    fac1 = 1.0f/fac1;
    fac2 = 1.0f/fac2;
    if(dim == 3)   (*sys).sim1.box.x = pow(vol1,1.0f/3.0f);
    else           (*sys).sim1.box.x = pow(vol1,1.0f/2.0f);
    if(dim == 3)   (*sys).sim2.box.x = pow(vol2,1.0f/3.0f);
    else           (*sys).sim2.box.x = pow(vol2,1.0f/2.0f);
    escala(fac1,&(*sys).sim1,atom1,dim);
    escala(fac2,&(*sys).sim2,atom2,dim);
  }
  (*volu).intentos++;  
}
//*************************************************** ELIGE *
int elige(int N){
  double rnd,prob; 
  int resp;

  rnd = Random();
  prob = 1.0f/(double)N;
  resp = (int)(rnd/prob) + 1;

  return resp;
}
//******************************************** IMPRIME PANTALLA *
void screen(int step,struct sys sys,int flag){
  double vol1,vol2,upot1,upot2,dens1,dens2;
  FILE *f;
  
  f = fopen("out.txt","a");
  switch(flag){
    case 1:   printf("STEP\tUPOT1\tUPOT2\tDENS1\tDENS2\tDR1\tDR2\tDV\n"); break;
    case 2:   {
      upot1 = sys.sim1.upot/(double)(sys.sim1.nat);
      upot2 = sys.sim2.upot/(double)(sys.sim2.nat);
      vol1  = sys.sim1.box.x*sys.sim1.box.y;
      if(sys.dim == 3)   vol1 *= sys.sim1.box.z;
      vol2  = sys.sim2.box.x*sys.sim2.box.y;
      if(sys.dim == 3)   vol2 *= sys.sim1.box.z;
      dens1 = (double)(sys.sim1.nat)/vol1;
      dens2 = (double)(sys.sim2.nat)/vol2;
      printf("%i  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf\n",step,upot1,upot2,dens1,dens2,sys.sim1.dr,sys.sim2.dr,sys.dv,vol1+vol2);
      fprintf(f,"%i  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf\n",step,upot1,upot2,dens1,dens2,sys.sim1.dr,sys.sim2.dr,sys.dv,vol1+vol2);
    }; break;
  }
  fclose(f);
}
//********************************************** CREAR PARTÍCULA *
void creaParti(struct sys sys,struct sim *simA,struct sim *simB,atomo atomA[],atomo atomB[],struct lj lj){
  double volA,volB,upnew,upDest,arg,beta = sys.temp,ax,ay,az;
  int nat,oCrea,oDest,dim = sys.dim,nA,nB,jren;

int swap=0,N=50;
  if((*simB).nat != 0){
while(swap < N  &&  (*simB).nat != 0){
    nat = (*simA).nat + (*simB).nat;
    nA = (*simA).nat;
    volA = (*simA).box.x * (*simA).box.y;
    if(dim == 3)   volA *= (*simA).box.z;

    oCrea = (*simA).nat;
    atomA[oCrea].pos.x = Random() * (*simA).box.x;
    atomA[oCrea].pos.y = Random() * (*simA).box.y;
    if(dim == 3)   atomA[oCrea].pos.z = Random() * (*simA).box.z;
    else           atomA[oCrea].pos.z = 0.0f;

    (*simA).nat += 1;

    upnew = ener_atom(oCrea,sys,*simA,atomA,lj);
   
    nB = (*simB).nat;
    volB = (*simB).box.x * (*simB).box.y;
    if(dim == 3)   volB *= (*simB).box.z;
  
    oDest = (int)(rand() % nB);
    upDest = ener_atom(oDest,sys,*simB,atomB,lj);
  
    //Criterio de aceptación
    arg = upnew - upDest + log( (volB*(double)(nA+1)) / (volA*(double)(nB)) )/beta;
    if(Random() < exp(-arg*beta)){
      (*simA).upot += upnew;
      (*simB).upot -= upDest; 
      (*simB).nat -= 1;
      nat = (*simA).nat + (*simB).nat;
      jren = oDest + 1;
      while(jren < nB){
        ax = atomB[jren].pos.x;
        ay = atomB[jren].pos.y;
        az = atomB[jren].pos.z;
        atomB[jren - 1].pos.x = ax;
        atomB[jren - 1].pos.y = ay;
        atomB[jren - 1].pos.z = az;
        jren++;
      }
    }
    else{
      atomA[(*simA).nat].pos.x = 0.0f;
      atomA[(*simA).nat].pos.y = 0.0f;
      atomA[(*simA).nat].pos.z = 0.0f;
      (*simA).nat -= 1;
    }
swap++;}
  }
}
//*************************************************** GIBBS *
void gibbs(struct sys *sys,atomo atom1[],atomo atom2[],struct lj lj,struct move *move1,struct move *move2,struct move *volu){
  int step,caso;

  step = 0;
  screen(step,*sys,1);

  while(step <= (*sys).mcStep){
    caso = elige(5);
    switch(caso){
      case 1:   mcmove(*sys,atom1,&((*sys).sim1),lj,move1); break;   
      case 2:   mcmove(*sys,atom2,&((*sys).sim2),lj,move2); break;
      case 3:   mcvol(sys,atom1,atom2,lj,volu); break;
      case 4:   creaParti(*sys,&((*sys).sim1),&((*sys).sim2),atom1,atom2,lj); break;
      case 5:   creaParti(*sys,&((*sys).sim2),&((*sys).sim1),atom2,atom1,lj); break;
    }
    if(step%(*sys).nPrint == 0){
      superPrint((*sys).dim,(*sys).sim1.nat,(*sys).sim2.nat,atom1,atom2,(*sys).sim1.box.x);
      screen(step,*sys,2);
    }
    if(step%(*sys).nAjusta == 0){
      ajustaDr(move1,*sys,&(*sys).sim1);
      ajustaDr(move2,*sys,&(*sys).sim2);
      ajustaDv(volu,&(*sys));
    }
    step++;
  }
}
//*************************************************************
