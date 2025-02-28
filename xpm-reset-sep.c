//closed chain resetting with hardcore exclusion---xpm file generation

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <time.h>
# include "mt19937ar.c"

# define L       100
# define tmax    500
# define ens     1

int main()
{ 
//clock_t tic = clock();

int    x,t,s[L]={},s0[L]={},site[L]={},pt[L]={},i,a,b,max,min,rho1;
double n,tot,p,q,r;
init_genrand(4089);
FILE   *fp;

p=40.0; q=40.0; r=17.0;
tot = p+q+r;

fp=fopen("output2.xpm","w");

//print-------------------------
fprintf(fp,"/* XPM */");                  fprintf(fp,"\n");
fprintf(fp,"static char * MyArray[] = {");fprintf(fp,"\n");
fprintf(fp,"\"100 500 3 1\",");     fprintf(fp,"\n");
fprintf(fp,"\"0      c #FFFFFF\",");fprintf(fp,"\n");
fprintf(fp,"\"1      c #FF5733\",");fprintf(fp,"\n");
fprintf(fp,"\"2      c #1638E5\",");fprintf(fp,"\n");
//fprintf(fp,"\"3      c #8838FR\",");fprintf(fp,"\n");
//print-------------------------


//doing all sites empty
for(x=0;x<L;x++){s[x]=0; s0[x]=0;}
//randomly selecting NN number of sites and fill with tagged particles
//a=0; i=1; while (a<NN){b=genrand_real1()*L ; if(!s[b]){s[b]=i; a++ ; i++;}  }
//creating an array of tagged particles
//a=1; for(x=0;x<L;x++){if(s[x]!=0){pt[a] = s[x]; a++;  } }


s[51]=1; s[54]=2; //s[58]=3;
//saving initial configuration
for(x=0;x<L;x++){s0[x] = s[x];}


//printing lattice array
//for(x=0;x<L;x++)   {printf("%d",s[x]);     } printf("--lattice\n");
//printing tagged particles array
//for(x=1;x<NN+1;x++){printf(" %d",pt[x]);   } printf("--particles\n");
//printing occupied sites array
//for(x=0;x<NN;x++)  {printf(" %d",site[x]); } printf("--sites\n");



for(t=0;t<tmax;t++){

//MC cycle begins
for(x=0;x<L;x++){
//choosing site
a=genrand_real1()*L;
//choosing rate
n=genrand_real1();
if(n>0.0 && n<p/tot){ if(s[a]!=0 && s[(a+1)%L]==0){s[(a+1)%L]=s[a]; s[a]=0; } }

else if(n>p/tot && n<(p+q)/tot){ if(s[a]!=0 && s[(a-1+L)%L]==0){s[(a-1+L)%L]=s[a]; s[a]=0; } }

else if(n>(p+q)/tot && n<1.0){ if(s[a]!=0){
                                 for(x=0;x<L;x++){ if(s0[x]==s[a]){b=x; } }
                                 if(a>b){max=a; min=b;} else {max=b; min=a;}
                                 rho1=0; for(x=min+1;x<max;x++){rho1 += s[x];}
                                 if(s[b]==0 && rho1==0){s[b]=s[a]; s[a]=0;}
                                          } 
                             }
                }//MC ends

fprintf(fp,"\"");
for(x=0;x<L;x++){fprintf(fp,"%d",s[x]); } fprintf(fp,"\",");fprintf(fp,"\n");

                   }//tmax ends
fprintf(fp,"};");

//clock_t toc = clock(); fprintf(fp,"#CPU RUNTIME: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
fclose(fp);
}//main ends
