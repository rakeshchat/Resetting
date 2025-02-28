//two-point correlation measurement between the different pairs
//using Monte-Carlo Simulation with scaling timestep dt
//tasep with translation(f) and pause(k)
//L(2)trelax(6)tmax(3)ens(2) = 
//L(3)trelax(6)tmax(3)ens(0) = 


# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <time.h>
# include "mt19937ar.c"


# define L          1000
# define trelax     1000000
# define tmax       1000
# define ens        100


int main()
{ clock_t tic = clock();

int    en,z,i,t,r,p,y[L]={0},yy[L]={0},yyy[L]={0},L1,M1,dist,mdist=50;
double rho0,rho1,rho2,rhov,rhoa,rhop,d,j,rho[L]={0},rate[L]={0},dt,max,e,f,tau,ft,A1[L]={0},A2[L]={0},A3[L]={0},A4[L]={0},A5[L]={0},B1[L]={0},B2[L]={0},B3[L]={0},B4[L]={0},B5[L]={0},Ciaa[L]={0},Cipp[L]={0},N_aa,N_pp,N_vv,N_ap,N_av,N_pv,N_pa,N_va,N_vp;
init_genrand(4589);
FILE   *fp;
//===========================================
fp=fopen("pair_dcorln_tau7","w");

//tau loop starts
//for(tau=0.01;tau<100.0;tau*=1.4){
//for(tau=2.0;tau<3.0;tau+=1.0){
N_aa=0;N_pp=0;N_vv=0;N_ap=0;N_av=0;N_pv=0;N_pa=0;N_va=0;N_vp=0;

j=0; rhoa=0; rhop=0; dt=0; max=0; tau=7.0;

//parameters and output file
d=0.5; e=1.0; ft=0.1; f=ft/tau; 

i=1; while (i<4){
rate[i] = e;       i++;
rate[i] = f;       i++;
rate[i] = 1.0/tau; i++;
                }
max = rate[1];
for(i=1;i<4;i++){if(rate[i] > max)max=rate[i]; }
dt = 1.0/max;

L1 = (int)L*max; M1 = (int)L*d;

for(z=0;z<L;z++){y[z]=0; yy[z]=0; yyy[z]=0; rho[z]=0; }
i=0; while (i<M1){z=genrand_real1()*L ; if(!y[z]){y[z]=1; i++ ; }  }

for(t=0;t<trelax;t++){
//---------------------------------------------
for(z=0;z<L1;z++){
    p=genrand_real1()*L;
    if(y[p]==1 && y[(p+1)%L]==0 && genrand_real1()<e*dt){y[p]=0; y[(p+1)%L]=1; }
    r=genrand_real1()*L;
    if(y[r]==1 && genrand_real1()<f*dt){y[r]=2; }
    else if (y[r]==2 && genrand_real1()<1.0/tau*dt){y[r]=1;}
                 }
                     }//trelax ends


//for(z=0;z<L;z++){fprintf(fp,"%d",y[z]); } fprintf(fp,"\n");

//the stationary state is saved
for(z=0;z<L;z++){yy[z] = y[z];}


for(en=0;en<ens;en++){

for(z=0;z<L;z++){y[z] = yy[z];}
for(z=0;z<L;z++){rho[z]=0;}
for(z=0;z<L;z++){A1[z]=0; A2[z]=0; A3[z]=0; A4[z]=0; A5[z]=0; B1[z]=0; B2[z]=0; B3[z]=0; B4[z]=0; B5[z]=0; }

for(t=0;t<tmax;t++){
for(z=0;z<L1;z++){
    p=genrand_real1()*L;
    if(y[p]==1 && y[(p+1)%L]==0 && genrand_real1()<e*dt){y[p]=0; y[(p+1)%L]=1; }
    r=genrand_real1()*L;
    if(y[r]==1 && genrand_real1()<f*dt){y[r]=2; }
    else if (y[r]==2 && genrand_real1()<1.0/tau*dt){y[r]=1;}
                 }//1 MC step ends


//measure number of pair particles
for(z=0;z<L;z++){
     if(y[z]==1 && y[(z+1)%L]==1){N_aa++;}
else if(y[z]==2 && y[(z+1)%L]==2){N_pp++;}
else if(y[z]==0 && y[(z+1)%L]==0){N_vv++;}
else if(y[z]==1 && y[(z+1)%L]==2){N_ap++;}
else if(y[z]==1 && y[(z+1)%L]==0){N_av++;}
else if(y[z]==2 && y[(z+1)%L]==0){N_pv++;}
else if(y[z]==2 && y[(z+1)%L]==1){N_pa++;}
else if(y[z]==0 && y[(z+1)%L]==1){N_va++;}
else if(y[z]==0 && y[(z+1)%L]==2){N_vp++;}
                }

rho0=0; for(z=0;z<L;z++){if(y[z]==0){rho0++ ;} }
rho1=0; for(z=0;z<L;z++){if(y[z]==1){rho1++ ;} }
rho2=0; for(z=0;z<L;z++){if(y[z]==2){rho2++ ;} }
rhov += rho0; rhoa += rho1; rhop += rho2;

                   }//tmax ends

                     }//ens ends

fprintf(fp,"  void= %lf\n  active= %lf\n  paused= %lf\n  N_aa= %lf\n  N_pp= %lf\n  N_vv= %lf\n  N_ap= %lf\n  N_av= %lf\n  N_pv= %lf\n  N_pa= %lf\n  N_va= %lf\n  N_vp= %lf\n",rhov/tmax/ens,rhoa/tmax/ens,rhop/tmax/ens,N_aa/tmax/ens,N_pp/tmax/ens,N_vv/tmax/ens,N_ap/tmax/ens,N_av/tmax/ens,N_pv/tmax/ens,N_pa/tmax/ens,N_va/tmax/ens,N_vp/tmax/ens);

//                }//tau loop ends

clock_t toc = clock(); fprintf(fp,"#CPU RUNTIME: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
fclose(fp);

}//main ends
