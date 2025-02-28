//two-point correlation measurement between two distant sites, correlation as a function of dist between two sites
//using Monte-Carlo Simulation with scaling timestep dt
//tasep with translation(f) and pause(k)
//L(2)trelax(6)tmax(3)ens(2) = 89"
//L(3)trelax(6)tmax(3)ens(1) = (for mdist=50) 2.7'


# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <time.h>
# include "mt19937ar.c"


# define L          1000
# define trelax     1000000
# define tmax       1000
# define ens        10


int main()
{ clock_t tic = clock();

int    en,z,i,t,r,p,y[L]={0},yy[L]={0},L1,M1,dist,mdist=50;
double d,j,rho[L]={0},rate[L]={0},dt,max,e,f,tau,ft,A1[L]={0},A2[L]={0},A3[L]={0},A4[L]={0},A5[L]={0},Ci[L]={0};
init_genrand(4589);
FILE   *fp;
//===========================================
fp=fopen("dist_corln_tau7","w");

//tau loop starts
//for(tau=0.01;tau<100.0;tau*=1.4){
//for(tau=2.0;tau<3.0;tau+=1.0){

j=0; dt=0; max=0; tau=7.0;

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

for(z=0;z<L;z++){y[z]=0; yy[z]=0; rho[z]=0; }
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

//the stationary state is saved
for(z=0;z<L;z++){yy[z] = y[z];}


for(en=0;en<ens;en++){

for(z=0;z<L;z++){y[z] = yy[z];}
for(z=0;z<L;z++){rho[z]=0;}
for(z=0;z<L;z++){A1[z]=0; A2[z]=0; A3[z]=0; A4[z]=0; A5[z]=0; }

for(t=0;t<tmax;t++){
for(z=0;z<L1;z++){
    p=genrand_real1()*L;
    if(y[p]==1 && y[(p+1)%L]==0 && genrand_real1()<e*dt){y[p]=0; y[(p+1)%L]=1; }
    r=genrand_real1()*L;
    if(y[r]==1 && genrand_real1()<f*dt){y[r]=2; }
    else if (y[r]==2 && genrand_real1()<1.0/tau*dt){y[r]=1;}
                 }//1 MC step ends


//measure correlation
for(dist=1;dist<mdist;dist++){
for(z=0;z<L;z++){
A1[dist] += y[z]*y[(z+dist)%L]; A3[dist] += y[(z+dist)%L]; A5[dist] += y[(z+dist)%L]*y[(z+dist)%L];
A2[dist] += y[z];               A4[dist] += y[z]*y[z];
                               
                }
                             }//dist loop ends
                   }//tmax ends

for(dist=1;dist<mdist;dist++){
Ci[dist] += (A1[dist]/tmax/L - A2[dist]/tmax/L*A3[dist]/tmax/L)/
(sqrt((A4[dist]/tmax/L - (A2[dist]/tmax/L)*(A2[dist]/tmax/L) )*(A5[dist]/tmax/L - (A3[(dist)%L]/tmax/L)*(A3[(dist)%L]/tmax/L))));

                             }

                     }//ens ends

for(dist=1;dist<mdist;dist++){fprintf(fp,"%d    %lf\n",dist,Ci[dist]/ens);}

//                }//tau loop ends


clock_t toc = clock(); fprintf(fp,"#CPU RUNTIME: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
fclose(fp);

}//main ends
