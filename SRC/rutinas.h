
void clear(){
  system("rm -f *.xyz *txt");
}
//******************************************************
void cafe(){

  printf("╔═════════════════════════════════════════════════════╗\n");
  printf("║     ( (                                             ║\n");
  printf("║      ) )                                            ║\n");
  printf("║  ...........          GIBBS ENSAMBLE                ║\n");
  printf("║  |         |]           MONTE CARLO                 ║\n");
  printf("║   \\       /                                         ║\n");
  printf("║    `-----´                                          ║\n");
  printf("╚═════════════════════════════════════════════════════╝\n\n");

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
    box[2] = 1.0f;
  }
}
//******************************************************
void compBox(struct sys *sys){
  double box[3];

  compBoxAux((*sys).dim, (*sys).sim1.dens, (*sys).sim1.nat, box);
  (*sys).sim1.box.x = box[0];
  (*sys).sim1.box.y = box[1];
  (*sys).sim1.box.z = box[2];
 
  compBoxAux((*sys).dim, (*sys).sim2.dens, (*sys).sim2.nat, box);
  (*sys).sim2.box.x = box[0];
  (*sys).sim2.box.y = box[1];
  (*sys).sim2.box.z = box[2];
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
void screen(int step,int samp,struct sys sys,int flag){
  double mu1,mu2,vol1,vol2,upot1,upot2,press1,press2,dens1,dens2;
  int Neq = (int)(sys.mcStep/2.0f);
  FILE *f;
  
  f = fopen("out.txt","a");
  switch(flag){
    case 1:   printf("STEP   \tUPOT1   \tUPOT2   \tDENS1   \tDENS2   \tMU1     \tMU2     \tPRESS1   \tPRESS2\n"); break;
    case 2:   {
      upot1 = sys.sim1.upot/(double)(sys.sim1.nat);
      upot2 = sys.sim2.upot/(double)(sys.sim2.nat);
      vol1  = sys.sim1.box.x*sys.sim1.box.y;
      if(sys.dim == 3)   vol1 *= sys.sim1.box.z;
      vol2  = sys.sim2.box.x*sys.sim2.box.y;
      if(sys.dim == 3)   vol2 *= sys.sim1.box.z;
      dens1 = (double)(sys.sim1.nat)/vol1;
      dens2 = (double)(sys.sim2.nat)/vol2;
      press1 = dens1 * sys.temp + sys.sim1.virial/((double)(sys.dim)*vol1);
      press2 = dens2 * sys.temp + sys.sim2.virial/((double)(sys.dim)*vol2);
      mu1 = -log(sys.sim1.chemPot / (double)(samp)) * sys.temp;
      mu2 = -log(sys.sim2.chemPot / (double)(samp)) * sys.temp;

      //printf("%lf\t%lf\t%lf\t%lf\n",dens1,dens2,mu1,mu2);

      printf("%i\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",step,upot1,upot2,dens1,dens2,mu1,mu2,press1,press2);
//printf("lz1: %lf\tlz2: %lf\n",sys.sim1.box.z,sys.sim2.box.z);
      if(step >= Neq){
        fprintf(f,"%i\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",step,upot1,upot2,dens1,dens2,mu1,mu2,press1,press2);
      }
    }; break;
  }
  fclose(f);
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
//**********************************************************
double Sech(double x){
  double secantHip;
  secantHip = 1.0f / cosh(x);
  return secantHip;
}
