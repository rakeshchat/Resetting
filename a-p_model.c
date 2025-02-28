//actual Wang dynamics
//using Monte-Carlo Simulation with scaling timestep dt
//tasep with translation(f) and pause(k)
//L(2)trelax(5)tmax(5)ens(0) = 


# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <time.h>
# include "mt19937ar.c"


# define L          100
# define trelax     1000
# define tmax       1000
# define ens        10


int main()
{ clock_t tic = clock();

int    en,z,i,t,r,p,y[L]={0},yy[L]={0},L1,M1;
double rho1,rho2,rhoa,rhop,d,j,rho[L]={0},rate[L]={0},dt,max,e,f,tau,ft;
init_genrand(4589);
FILE   *fp;
//===========================================
fp=fopen("asd","w");

for(tau=0.01;tau<100;tau*=1.2){j=0; rhoa=0; rhop=0; dt=0; max=0;
//for(d=0.02;d<1.0;d+=0.02){j=0; rhoa=0; rhop=0;

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

for(z=0;z<L;z++){y[z]=0; yy[z]=0;}
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

for(z=0;z<L;z++){yy[z] = y[z];}

for(en=0;en<ens;en++){
for(z=0;z<L;z++){y[z] = yy[z];}

for(t=0;t<tmax;t++){
//---------------------------------------------
for(z=0;z<L1;z++){
    p=genrand_real1()*L;
    if(y[p]==1 && y[(p+1)%L]==0 && genrand_real1()<e*dt){y[p]=0; y[(p+1)%L]=1; j++; }
    r=genrand_real1()*L;
    if(y[r]==1 && genrand_real1()<f*dt){y[r]=2; }
    else if (y[r]==2 && genrand_real1()<1.0/tau*dt){y[r]=1;}
                    }
                  }//tmax ends
rho1=0; for(z=0;z<L;z++){if(y[z]==1){rho1++ ;} }
rho2=0; for(z=0;z<L;z++){if(y[z]==2){rho2++ ;} }
rhoa += rho1; rhop += rho2;
              }//ens ends

//fprintf(fp,"%lf   %lf\n",tau,j/L/tmax/ens);   
fprintf(fp,"%lf   %lf   %lf\n",tau,rhoa/L/ens,rhop/L/ens);
      }//tau loop ends
//=========================================================


clock_t toc = clock(); fprintf(fp,"#CPU RUNTIME: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
fclose(fp);
}//main ends
