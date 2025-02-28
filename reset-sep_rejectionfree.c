//open chain resetting with hardcore exclusion

# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <time.h>
# include "mt19937ar.c"

# define L       10
# define NN      5
# define tmax    1
# define ens     1

int main()
{ clock_t tic = clock();

int    x,t,s[L]={},ss[L]={},site[L]={},pt[L]={},i,a,b,c,m;
double n,tot,p,q,r;
init_genrand(4089);
FILE   *fp;


p=20.0; q=20.0; r=5.0;
tot = p+q+r;
m=NN;

fp=fopen("outputfile-1","w");

//doing all sites empty
for(x=0;x<L;x++){s[x]=0; ss[x]=0;}
//randomly selecting NN number of sites and fill with tagged particles
a=0; i=1; while (a<NN){b=genrand_real1()*L ; if(!s[b]){s[b]=i; a++ ; i++;}  }
//creating an array of tagged particles
a=1; for(x=0;x<L;x++){if(s[x]!=0){pt[a] = s[x]; a++;  } }
//creating an array of occupied sites
a=0; for(x=0;x<L;x++){if(s[x]!=0){site[a] = x; a++;  } }

//saving initial configuration
for(x=0;x<L;x++){ss[x] = s[x];}

//printing lattice array
for(x=0;x<L;x++)   {printf("%d",s[x]);     } printf("--lattice\n");
//printing tagged particles array
for(x=1;x<NN+1;x++){printf(" %d",pt[x]);   } printf("--particles\n");
//printing occupied sites array
for(x=0;x<NN;x++)  {printf(" %d",site[x]); } printf("--sites\n");





//for(t=0;t<tmax;t++){

//MC cycle begins
for(x=0;x<NN;x++){
//a = randomly select from occupied site array
a=genrand_real1()*NN;
//printf("c=%d x=%d\n",c,x);
//b = site number
//b=site[a];
//c = which tagged particle is at that occupied site
//c=s[site[a]];
//printf("%d  %d  %d\n",a,b,c);


n=genrand_real1();
if(n>0.0 && n<p/tot)      { if(s[(site[a]+1)%L]==0){s[(site[a]+1)%L]=s[site[a]]; s[site[a]]=0; } }

else if(n>p/tot && n<(p+q)/tot){ if(s[(site[a]-1+L)%L]==0){s[(site[a]-1+L)%L]=s[site[a]]; s[site[a]]=0; } }

else if(n>(p+q)/tot && n<1.0)     { if(s[(site[a]+1)%L]==0){s[(site[a]+1)%L]=s[site[a]]; s[site[a]]=0; } }


//creating an array of tagged particles after update
//a=1; for(x=0;x<L;x++){if(s[x]!=0){pt[a] = s[x]; a++;  } }
//creating an array of occupied sites after update
//a=0; for(x=0;x<L;x++){if(s[x]!=0){site[a] = x; a++;  } }













                }//MC ends

//                   }//tmax ends


//for(x=0;x<L;x++){printf("%d",s[x]); } printf("----lattice\n");


//clock_t toc = clock(); fprintf(fp,"#CPU RUNTIME: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
fclose(fp);
}//main ends
