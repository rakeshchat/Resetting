//two-point correlation measurement between two neighbouring sites
//using Monte-Carlo Simulation with scaling timestep dt
//tasep with translation(f) and pause(k)
//L(2)trelax(6)tmax(3)ens(4) = 04'
//L(3)trelax(6)tmax(3)ens(4) = 35'


# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <time.h>
# include "mt19937ar.c"


# define L          100
# define trelax     1000000
# define tmax       1000
# define ens        1000


int main()
{ clock_t tic = clock();

int    en,z,i,t,r,p,y[L]={0},yy[L]={0},yyy[L]={0},L1,M1;
double rho1,rho2,rhoa,rhop,d,j,rho[L]={0},rate[L]={0},dt,max,e,f,tau,ft,A1[L]={0},A2[L]={0},A3[L]={0},A4[L]={0},A5[L]={0},Ci1[L]={0},Ci2[L]={0};;
init_genrand(4589);
FILE   *fp;
//===========================================
fp=fopen("site_corln","w");

//tau loop starts
//for(tau=0.01;tau<100.0;tau*=1.4){
//for(tau=2.0;tau<3.0;tau+=1.0){

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
for(z=0;z<L;z++){A1[z]=0; A2[z]=0; A3[z]=0; A4[z]=0; A5[z]=0; }

for(t=0;t<tmax;t++){
for(z=0;z<L1;z++){
    p=genrand_real1()*L;
    if(y[p]==1 && y[(p+1)%L]==0 && genrand_real1()<e*dt){y[p]=0; y[(p+1)%L]=1; }
    r=genrand_real1()*L;
    if(y[r]==1 && genrand_real1()<f*dt){y[r]=2; }
    else if (y[r]==2 && genrand_real1()<1.0/tau*dt){y[r]=1;}
                 }//1 MC step ends

//measuring average occupancy of sites
for(z=0;z<L;z++){yyy[z] = y[z];}
for(z=0;z<L;z++){if(y[z]==2)y[z]=1;}
for(z=0;z<L;z++){rho[z] += y[z];}
for(z=0;z<L;z++){y[z] = yyy[z];}

//measure correlation
for(z=0;z<L;z++){
A1[z] += y[z]*y[(z+1)%L];  A2[z] += y[z];  A3[z] += y[(z+1)%L]; 
A4[z] += y[z]*y[z]; A5[z] += y[(z+1)%L]*y[(z+1)%L];                           
                }
                   }//tmax ends

for(z=0;z<L;z++){
Ci1[z] += (A1[z]/tmax - A2[z]/tmax*A3[z]/tmax)/(sqrt((A4[z]/tmax - (rho[z]/tmax)*(rho[z]/tmax) )*(A5[z]/tmax - (rho[(z+1)%L]/tmax)*(rho[(z+1)%L]/tmax))));
                }
                     }//ens ends

for(z=0;z<L;z++){fprintf(fp,"%lf   %lf\n",z*1.0/L,Ci1[z]/ens);}

//                }//tau loop ends


//rho1=0; for(z=0;z<L;z++){if(y[z]==1){rho1++ ;} }
//rho2=0; for(z=0;z<L;z++){if(y[z]==2){rho2++ ;} }
//rhoa += rho1; rhop += rho2;


clock_t toc = clock(); fprintf(fp,"#CPU RUNTIME: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
fclose(fp);

}//main ends
