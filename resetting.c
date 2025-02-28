//system with ABSORBING BOUNDARY: SEMI-INFINITE CASE W/WO reset,potential 
//local time distribution, mean with local time
//RR(6)t(10)dt(0.0001) = 1hr 20'
//for PDF RR=10^6, for mean RR=10^5


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
double time,t,mm,a1,a2,a3,PI,R1,R2,y1,eta,D,dt,r,r1,s,x_0,x_f,lambda,loc_time,total_loc,z_0,mfp,A[RR]={},bin,offset,mean,mean_theory=0,mean_theory_r=0;
init_genrand(4781);
FILE   *fp;


fp=fopen("pdf_semiinft_t20_r0_l-1.1","w");
//fp=fopen("t_mean_semiintf_r0_l-1.1","w");


D=0.5;         //diffusion rate
t=20.0;        //total observed time
dt=0.0001;     //timestep 0.0001
r=0.0;         //reset probability
s=1.0;         //reset position
lambda=-1.1;   //lambda positive:bounded potential
bin=0.001;     //bin size


//for(t=1.0;t<50;t+=1){

for(i=0;i<RR;i++){A[i]=0;}
      
//realization loop starts
for(realz=0;realz<RR;realz++){

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

//semi-infinite boundary
else if(x_f<0.0){break; }

                 }//timestep ends

total_loc = count*dt;

j    = total_loc/bin;
A[j] = A[j] + 1;

              }//realization ends


//local time distribution (PDF) and print
for(i=0;i<RR;i++)fprintf(fp,"%lf    %lf\n",i*bin/mfp,A[i]/bin/RR*mfp);


//theory mean
//mean_theory=1.0/D*(1.0-exp(D*t)*erfc(sqrt(D*t)));
//mean_theory_r=sqrt(r*1.0/D)*1.0/(1.0 + 1.0/(tanh(sqrt(1.0*r/D)*s)))*t;


//numerical mean (<mean>) of local time distribution and print
//offset=0; offset = A[0]/RR*bin/mfp;
//mean=0; for(i=0;i<RR;i++){mean += i*1.0*A[i]/RR*bin/mfp;}  fprintf(fp,"%lf    %lf\n",t,mean); }//t loop ends



clock_t toc = clock(); fprintf(fp,"#CPU RUNTIME: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
fclose(fp);
}//main ends
