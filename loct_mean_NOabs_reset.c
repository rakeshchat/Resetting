//system with NO ABSORBING BOUNDARY: FREE CASE W/WO reset,potential
//local time distribution, mean with local time
//RR(6)t(10)dt(0.0001) = 1hr 20'
//for local time dist. RR=10^6, for mean calc. RR=10^5


# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <time.h>
# include "mt19937ar.c"

# define RR     1000000


int main()
{
clock_t tic = clock();

int    j,i,realz,count;
double time,t,mm,a1,a2,a3,PI,R1,R2,y1,eta,D,dt,r,r1,s,x_0,x_f,lambda,loc_time,total_loc,z_0,mfp,A[RR]={},bin,mean,mean_exact_large=0,mean_exact_small=0,mean_exact_large_r=0;
init_genrand(4089);
FILE   *fp;


fp=fopen("localt_dist_NOabs-s1-r.6-t10","w");
//fp=fopen("t_mean_NOabs-s1-r0","w");
//fp=fopen("asd-pdft10","w");


D=0.5;         //diffusion rate
t=10.0;        //total observed time
dt=0.0001;     //timestep 0.0001
r=0.6;         //reset probability
s=1.0;         //reset position
lambda=0.0;    //lambda positive:bounded potential

bin=0.001;     //bin size


//for(t=5;t<101;t+=5){

for(i=1;i<RR+1;i++){A[i]=0;}
      
//realization loop starts
for(realz=1;realz<RR+1;realz++){

count=0; 
total_loc=0;   //survival time at each realization
x_0=s;
z_0=s;


//timestep loop starts
for(time=0;time<t+dt;time+=dt){ 

//Box-Muller to generate Gaussian Random Variable
a1=genrand_real1(); a2=genrand_real1();
PI=4.0*atan(1.0);
R1=sqrt(-2.0*log(a1));
R2=2.0*PI*a2;
y1=R1*cos(R2);

//Noise in Brownian motion
eta=sqrt(2.0*D*dt)*y1;

//reset rate
a3=genrand_real1();
r1=r*dt;

if(a3<r1){x_f = s; }                      //reset to initial position at 1
    else {x_f = x_0 - lambda*dt + eta; }    //diffusion with a potential


//define mfp
//local time spend within a region should be less than mfp
mfp = sqrt(2.0*D*dt);

//final position update
x_0 = x_f;

//time spend at a position (z_0 +- mfp/2) within t to (t+dt) is local time
if(x_f < (z_0+mfp/2.0) && x_f > (z_0-mfp/2.0) ){count++; }

                 }//timestep ends

total_loc = count*dt;

j    = total_loc/bin;
A[j] = A[j] + 1;

              }//realization ends


//local time distribution
for(i=1;i<RR+1;i++)fprintf(fp,"%lf    %lf\n",i*bin/mfp,A[i]/bin/RR*mfp);
 
//calculation of numerical mean
//mean=0; for(i=1;i<RR+1;i++){mean += i*1.0*A[i]/RR*bin/mfp;}

//FROM EXACT CALCULATION without reset
//mean_exact_small=sqrt(t/(4*3.14*D));
//mean_exact_large=0;

//FROM EXACT CALCULATION with reset
//mean_exact_large_r=sqrt(t/(4*3.14*D))*exp(-r*t) + (1+2*r*t)/(4*sqrt(r*D))*erf(sqrt(r*t));


//printf("%lf    %lf    %lf\n",t,mean,mean_exact_large);
//fprintf(fp,"%lf      %lf\n",t,mean);

//    }//t loop ends


clock_t toc = clock(); fprintf(fp,"#CPU RUNTIME: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
	
fclose(fp);
}//main ends
