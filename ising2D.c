/*********************************************************************** 
  Esse programa foi baseado no material do livro: 
  "Programmazione Scientifica",
  Pearson Education ed. (2006), by Barone, Marinari, Organtini and
  Ricci-Tersenghi. ISBN 8871922425
  http://chimera.roma1.infn.it/SP/PROGRAMMI/cap19/cap19_Ising2Dshort.c

  Para uma versão levemente mais optimizada, mas menos legível:
  http://chimera.roma1.infn.it/SP/PROGRAMMI/cap19/cap19_Ising2D.c
  
  Copyright (C) 2006 Federico Ricci-Tersenghi 
  (Federico.Ricci@roma1.infn.it)
  
  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.
***********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

//semente PRNG. 
#define SEED 12974236
//para nunca repetir mesma hist. termica: time(NULL)

#define L 100

#define L2 (L*L)
#define N L2
#define DIM 2
#define DIM2 (2 * DIM) //numero de vizinhos


#define temperatura 2.0 //Tc~2.27
#define tmax 10000


clock_t start, stop;

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

void Init_Prob(float *prob) {
  int i;

  if (temperatura > 0.0) {
    for (i = 0; i <= DIM2; i++)
      prob[i] = exp(-2.0 * i / temperatura);
  } else {
    for (i = 0; i <= DIM2; i++)
      prob[i] = 0.0;
  }
}

//atualização sequencial
void oneMCS(int *s, float *prob, int *pMag, int *pEner) {
  int site, ix, iy, soma;

  site = 0;
  for (iy = 0; iy < L; iy++) {
    for (ix = 0; ix < L; ix++, site++) {
      soma = s[site] * (s[ix + L * ((iy+L-1)%L)] +
		       s[((ix+L-1)%L) + L * iy] +
		       s[ix + L * ((iy+1)%L)] +
		       s[((ix+1)%L) + L * iy]);
      if (soma <= 0 || FRANDOM < prob[soma]) {
	s[site] = -s[site];
	*pMag += 2 * s[site];  //tire comentarios p/ medir a cada MCS
	*pEner += 2 * soma;
      }
    }
  }
  
  return;
}

//atualização aleatória
void oneMCSrand(int *s, float *prob, int *pMag, int *pEner) {
  int i, site, ix, iy, soma;

  for (i = 0; i < N; i++) {
    site = FRANDOM * N;
    ix = site % L;
    iy = (site - (site % L)) / L;
    
    soma = s[site] * (s[ix + L * ((iy+L-1)%L)] +
		      s[((ix+L-1)%L) + L * iy] +
		      s[ix + L * ((iy+1)%L)] +
		      s[((ix+1)%L) + L * iy]);
    if (soma <= 0 || FRANDOM < prob[soma]) {
      s[site] = -s[site];
      //*pMag += 2 * s[site];
      //*pEner += 2 * soma;
    }
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
      *ene -= s[site] * (s[ix + L * ((iy+L-1)%L)] +
			 s[((ix+L-1)%L) + L * iy] +
			 s[ix + L * ((iy+1)%L)] +
			 s[((ix+1)%L) + L * iy]);
      *mag += s[site];
      site++;
    }
  }
  
  return;
}

int main(int argc, char *argv[]) {
  int i,t,ene,mag,ix,iy;
  int s[N];
  float prob[DIM2+1];
  float norm = 1.0/N;
  
  //  inicializa PRNG
  myrand = SEED;
  Init_Random();

  Init_Prob(prob);

  // ---------- inicializa spins -----------------
  // c.i. desordenada: s[i] = pm1  ferro: s[i] = 1
  // ---------------------------------------------
  
  for(i=0;i<N;i++)
    s[i] = pm1;  //condição inicial desordenada
  
  measures(s,&ene,&mag);

  start = clock();

  for( t=0 ; t<tmax ; t++ ){
    oneMCS(s,prob,&mag,&ene);
    
    //para visualizar com DynamicLattice:
    /* for(iy=0;iy<L;iy++) */
    /*   for(ix=0;ix<L;ix++) */
    /* 	printf("%d %d %d\n",ix,iy,s[ix+L*iy]); */
    /* printf("#draw\n"); */

    //...medindo:
    printf("%i %f %f\n", t,  mag*norm ,ene*norm/2); 
  }
  
  stop = clock();
  //  printf("%d %d %f\n",tmax, N,(float)(stop-start)/(CLOCKS_PER_SEC));
  
  return 0;
}
