/* 
   Algoritmo de Metropolis para modelo de Ising 2D. 
   Implementação sem otimizações.
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define SEED 22171245

#define L 300
#define L2 (L*L)
#define N L2

#define temperatura 1.5  //Tc~2.27
#define tmax 1000

clock_t start, stop;

int right[N], left[N], up[N], down[N];

/*---------------------------------------------------------------------*/
/*----------P. Random Number Generator by Parisi & Rapuano-------------*/
/*---------------------------------------------------------------------*/
#define FNORM   (2.3283064365e-10)
#define RANDOM  ((ira[ip++] = ira[ip1++] + ira[ip2++]) ^ ira[ip3++])
#define FRANDOM (FNORM * RANDOM)
#define pm1 ((FRANDOM > 0.5) ? 1 : -1)

unsigned myrand, ira[256];
unsigned char ip, ip1, ip2, ip3;

unsigned rand4init(void) {
  unsigned long long y;
  
  y = (myrand*16807LL);
  myrand = (y&0x7fffffff) + (y>>31);	  	  	                  
  if (myrand&0x80000000)
    myrand = (myrand&0x7fffffff) + 1;
  return myrand;
}

void Init_Random(void) {
  unsigned i;
  
  ip=128;
  ip1=ip-24;
  ip2=ip-55;
  ip3=ip-61;
  
  for (i=ip3; i<ip; i++)
    ira[i] = rand4init();
}
/*---------------------------------------------------------------------*/
/*---------------------------------------------------------------------*/
/*---------------------------------------------------------------------*/

void vizinhos_ccp2d(void)
{
  int i;
  
  for (i=0; i<L2; i++){

    if (i % L==L-1)         /* ultima coluna */
      right[i] = i-L+1;
    else 
      right[i] = i+1; 
    
    if (i % L==0)           /* primeira coluna */
      left[i] = i+L-1;
    else 
      left[i] = i-1;
    
    if (i<L)                /* primeira linha */
     up[i] = L2-L+i; 
    else 
      up[i] = i - L;
    
    if (i>=L2-L)             /* ultima linha  */
      down[i] = (i % L);
    else 
      down[i] = i + L;
    
  }  
  
  return;
}

void oneSweep2D(int *s) {
  int i, site, soma;

  for (i = 0; i < N; i++) {
    site = FRANDOM * N;
    
    soma = s[site] * (s[right[site]] +  s[left[site]] +
		     s[up[site]] + s[down[site]]);

    if (soma <= 0)
      s[site] = -s[site];
    else if(FRANDOM < exp(-2.0 * soma / temperatura))
      s[site] = -s[site];
  }
  
  return;
}

void measures(int *s, int *ene, int *mag){
  int ix,iy,site;

  *ene=0;
  *mag=0;
  site=0;
  for (iy = 0; iy < L; iy++) {
    for (ix = 0; ix < L; ix++) {
      *ene -= s[site] * (s[right[site]] +  s[left[site]] +
			 s[up[site]] + s[down[site]]);
      *mag += s[site];
      site++;
    }
  }
  
  return;
}

int main(int argc, char *argv[]) {
  int i,t,ene,mag,ix,iy;
  int s[N];
  float norm = 1.0/N;
  
  //  inicializa PRNG
  myrand = SEED;
  Init_Random();

  // inicializa/calcula vizinhos
  vizinhos_ccp2d();

  // inicializa spins
  for(i=0;i<N;i++)
    s[i] = pm1;

  start = clock();

  for(t=0;t<tmax;t++){
    oneSweep2D(s);
    //measures(s,&ene,&mag);
    //printf("%i %f %f\n",t,ene*norm, mag*norm);
    // fflush(stdout);

    for(ix=0;ix<L;ix++)
      for(iy=0;iy<L;iy++)
    	printf("%d %d %d\n",ix,iy,s[iy+L*ix]);
    printf("#draw\n");


  }
  stop = clock();
  //  printf("%d %d %f\n",NmaxMCS, N, (float) (stop-start)/ ( CLOCKS_PER_SEC ));
  
  //printf("%f\n",x);

  return 0;
}
