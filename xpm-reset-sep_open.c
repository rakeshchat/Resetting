//open chain with one reset particle, hardcore exclusion, diffusion at boundaries

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <time.h>
# include "mt19937ar.c"

# define L       100
# define tmax    500
# define tany    1

int main()
{ 
clock_t tic = clock();

int    x,t,tt,s[L]={0},s0[L]={0},site[tmax]={0},i,a,b,c,max,min,rho1,M1;
double n,tot,p,q,r,fill,C1[tmax]={0};
init_genrand(4089);
FILE   *fp;

fill=0.2;
p=40.0; q=40.0; r=8.0;
tot = p+q+r;


fp=fopen("fileout2.xpm","w");

//print-------------------------
fprintf(fp,"/* XPM */");                  fprintf(fp,"\n");
fprintf(fp,"static char * MyArray[] = {");fprintf(fp,"\n");
fprintf(fp,"\"100 500 3 1\",");     fprintf(fp,"\n");
fprintf(fp,"\"0      c #FFFFFF\",");fprintf(fp,"\n");
fprintf(fp,"\"1      c #EE8167\",");fprintf(fp,"\n");
fprintf(fp,"\"2      c #1638E5\",");fprintf(fp,"\n");
//fprintf(fp,"\"3      c #8838FR\",");fprintf(fp,"\n");
//print-------------------------



//doing all sites empty
for(x=0;x<L;x++){s[x]=0; s0[x]=0;}
//filling normal particles
i=0; while (i<fill*L){x=genrand_real1()*L ; if(!s[x]){s[x]=1; i++ ; }  }
//one reset particle
M1=L/2; if(s[M1]==0 || s[M1!=0]){s[M1]=2; }

//saving initial configuration
for(x=0;x<L;x++){s0[x] = s[x];}


//initialization
for(x=0;x<L;x++){s[x] = s0[x];}


for(t=0;t<tmax;t++){

//MC cycle begins
for(x=0;x<L;x++){
//choosing site
a=genrand_real1()*L;
//choosing rate
n=genrand_real1();

switch(a){
case 0:
    switch((int)s[0])
        {
    case 0:
        if(n>0.0 && n<p/tot){s[0]=1; } break;
    case 1:
        if(n>p/tot && n<(p+q)/tot){s[0]=0; }
        else if(n>0.0 && n<p/tot && !s[1]){s[0]=0; s[1]=1; } break;
    case 2:
        if(n>(p+q)/tot && n<1.0){
                                 for(x=0;x<L;x++){ if(s0[x]==s[0]){b=x;} }
                                 rho1=0; for(x=1;x<b;x++){rho1 += s[x];}
                                 if(s[b]==0 && rho1==0){s[b]=s[0]; s[0]=0;}
                                } break;
        } break;
case L-1:
    switch((int)s[L-1])
        {
    case 0:
        if(n>p/tot && n<(p+q)/tot){s[L-1]=1; } break;
    case 1:
        if(n>0.0 && n<p/tot){s[L-1]=0; }
        else if(n>p/tot && n<(p+q)/tot && !s[L-2]){s[L-1]=0; s[L-2]=1; } break;
    case 2:
        if(n>(p+q)/tot && n<1.0){
                                 for(x=0;x<L;x++){ if(s0[x]==s[L-1]){b=x;} }
                                 rho1=0; for(x=b+1;x<L-1;x++){rho1 += s[x];}
                                 if(s[b]==0 && rho1==0){s[b]=s[L-1]; s[L-1]=0;}
                                } break;
        } break;
default:
        if(n>0.0 && n<p/tot){ if(s[a]!=0 && s[(a+1)]==0){s[(a+1)]=s[a]; s[a]=0; } }
        else if(n>p/tot && n<(p+q)/tot){ if(s[a]!=0 && s[(a-1)]==0){s[(a-1)]=s[a]; s[a]=0; } }
        else if(n>(p+q)/tot && n<1.0){ if(s[a]!=0){
                                 for(x=0;x<L;x++){ if(s0[x]==s[a]){b=x;} }
                                 if(a>b){max=a; min=b;} else {max=b; min=a;}
                                 rho1=0; for(x=min+1;x<max;x++){rho1 += s[x];}
                                 if(s[b]==0 && rho1==0){s[b]=s[a]; s[a]=0; }
                                                  }
                                     } break;

                             }//switch ends

                }//MC ends
//=======================

fprintf(fp,"\"");
for(x=0;x<L;x++){fprintf(fp,"%d",s[x]); } fprintf(fp,"\",");fprintf(fp,"\n");

                   }//tmax ends

//=======================

fprintf(fp,"};");

//clock_t toc = clock(); fprintf(fp,"#CPU RUNTIME: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
fclose(fp);
}//main ends
