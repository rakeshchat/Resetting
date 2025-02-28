//closed chain resetting with hardcore exclusion

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <time.h>
# include "mt19937ar.c"

# define L       500
# define tmax    500
# define tany    1000

int main()
{ 
clock_t tic = clock();

int    x,t,tt,s[L]={},s0[L]={},site[L]={},pt[L]={},i,a,b,c,d,e,f,max,min,rho1;
double n,tot,p,q,r;
init_genrand(4089);
FILE   *fp;

p=40.0; q=60.0; r=17.0;
tot = p+q+r;

fp=fopen("fileout1","w");

//doing all sites empty
for(x=0;x<L;x++){s[x]=0; s0[x]=0;}

s[230]=1; s[245]=2; s[259]=3; s[279]=4;

//saving initial configuration
for(x=0;x<L;x++){s0[x] = s[x];}


for(tt=0;tt<tany;tt++){

//initialization
for(x=0;x<L;x++){s[x] = s0[x];}


for(t=0;t<tmax;t++){

//MC cycle begins
for(x=0;x<L;x++){
//choosing site
a=genrand_real1()*L;
//choosing rate
n=genrand_real1();
if(n>0.0 && n<p/tot){ if(s[a]!=0 && s[(a+1)]==0){s[(a+1)]=s[a]; s[a]=0; } }

else if(n>p/tot && n<(p+q)/tot){ if(s[a]!=0 && s[(a-1)]==0){s[(a-1)]=s[a]; s[a]=0; } }

else if(n>(p+q)/tot && n<1.0){ if(s[a]!=0){
                                 for(x=0;x<L;x++){ if(s0[x]==s[a]){b=x;} }
                                 if(a>b){max=a; min=b;} else {max=b; min=a;}
                                 rho1=0; for(x=min+1;x<max;x++){rho1 += s[x];}
                                 if(s[b]==0 && rho1==0){s[b]=s[a]; s[a]=0;}
                                          } 
                             }
                }//MC ends

//for(x=0;x<L;x++){fprintf(fp,"%d",s[x]); } fprintf(fp,"\n");

                   }//tmax ends

for(x=0;x<L;x++){ if(s[x]==1){c=x;} } 
for(x=0;x<L;x++){ if(s[x]==2){d=x;} }
for(x=0;x<L;x++){ if(s[x]==3){e=x;} }
for(x=0;x<L;x++){ if(s[x]==4){f=x;} }

fprintf(fp,"%d  %d  %d  %d\n",c,d,e,f);

             }//tany ends

//clock_t toc = clock(); fprintf(fp,"#CPU RUNTIME: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
fclose(fp);
}//main ends
