/* Encabezado donde se declaran las estructuras de datos 
   y se declaran los prototipos de las subrutinas */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "ran2.h"
#include "Random.h"

#define SI 1
#define NO 0

//**************************************************** STRUCTS *
struct box{double x;double y;double z;};
struct sim{int nat;double dens;double upot;double chemPot;double virial;struct box box;double dr;};
struct sys{int nat;double dens;double mcStep;int nAjusta;int nPrint;double accep;double rcut;double temp;double dv;int dim;char potential[5];struct sim sim1,sim2;struct box box;};
struct lj{double sig;double eps;};
struct vector{double x;double y;double z;};
typedef struct {struct vector pos;} atomo;
struct move{int accep; int intentos;};

#include "rutinas.h"
//**************************************************** PROTOTIPOS *
void clear();
void compBoxAux(int dim,double dens,int nat,double box[]);
void compBox(struct sys *sys);
int elige(int N);
void screen(int step,int samp,struct sys sys,int flag);
void superPrint(int dim,int nat1,int nat2,atomo atom1[],atomo atom2[],double lx);
void printXYZ(char outxyz[],int dim,struct sim sim,atomo atom[]);
void readData(char inFile[],struct sys *sys,struct lj *lj);
void minima(int dim,struct box box,double *dx,double *dy,double *dz);
int traslape(int dim,struct sim sim,atomo atom[],int atoms,double xo,double yo,double zo);
void makeAtoms(int dim,struct sim sim,atomo atom[]);
double chemPot(int nat,double beta,double vol,double U);
void potencial(double *u,double *vi,struct sys sys,struct sim sim,atomo atom[],int i,int j,struct lj lj);
void energia_total(double *u,double *vi,struct sys sys,struct sim sim,atomo atom[],struct lj lj);
void ener_atom(int o,double *uo,double *vio,struct sys sys,struct sim sim,atomo atom[],struct lj lj);
void pbc(double *x,double *y,double *z,struct box box);
void init(double *u,double *vi,atomo atom1[],atomo atom2[],struct move *move1,struct move *move2,struct move *volu,struct sys *sys,struct lj lj);
void mcmove(struct sys sys,atomo atom[],struct sim *sim,struct lj lj,struct move *move);
void ajustaDr(struct move *move,struct sys sys,struct sim *sim);
void ajustaDv(struct move *volu,struct sys *sys);
void escala(double factor,struct sim *sim,atomo atom[],int dim);
void mcvol(struct sys *sys,atomo atom1[],atomo atom2[],struct lj lj,struct move *volu);
void creaParti(struct sys sys,struct sim *simA,struct sim *simB,atomo atomA[],atomo atomB[],struct lj lj,int *samp);
void gibbs(struct sys *sys,atomo atom1[],atomo atom2[],struct lj lj,struct move *move1,struct move *move2,struct move *volu,int *samp);
