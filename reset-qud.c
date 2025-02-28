//reset vs survival time with an absorbing boundary for quadratic potential
//RR(3)tmax(5)r(0.001-50-1.3)D(0.5) = 



# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <time.h>
# include "mt19937ar.c"

# define RR     1000
# define tmax   100000


int main()
{
clock_t tic = clock();

int    realz,count;
double t,mm,a1,a2,a3,PI,R1,R2,y1,eta,D,dt,r,r1,x_0,x_f,lambda,sv_time,total_sv;
init_genrand(4089);
FILE   *fp;

//fp=fopen("asd","w");
fp=fopen("rsvt-qud_D.5l-1.1","w");


D=0.5;        //diffusion rate
dt=0.001;     //timestep
//r=30.1;      //reset prob.
lambda=-1.1;   //lambda positive:bounded potential


//reset loop starts
for(r=0.001;r<50.0;r*=1.3){
count=0; total_sv=0;
x_0=1.0;

      

//realization loop starts
for(realz=0;realz<RR;realz++){
x_0=1.0;

//timestep loop starts
for(t=dt;t<tmax;t+=dt){
sv_time=0.0;  //survival time at each realization

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

if(a3<r1){x_f = 1.0; }                          //reset to initial position at 1
    else {x_f = x_0 - lambda*x_0*dt + eta; }    //diffusion with a potential

//final position update
x_0 = x_f;

//absorbing condition when absorbed going out of the timestep loop
if(x_f<0.0){sv_time = t;  count++;  
total_sv += sv_time;
break;     
           }

                 }//timestep ends

              }//realization ends

//printf("%lf   %lf\n",r,total_sv/count);
fprintf(fp,"%lf   %lf\n",r,total_sv/count);

       }//reset loop ends

clock_t toc = clock(); fprintf(fp,"#CPU RUNTIME: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
fclose(fp);
}//main ends
