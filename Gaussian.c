#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int main(void){
  int iter,niter=100;
  int naccept;
  double step_size=0.5e0;
  double x,backup_x,dx;
  double action_init, action_fin;
  double metropolis;

  srand((unsigned)time(NULL)); 

  /*********************************/
  /* Set the initial configuration */
  /*********************************/      
  x=0e0;
  naccept=0;
  /*************/
  /* Main loop */
  /*************/
  for(iter=1;iter<niter+1;iter++){
    backup_x=x;    
    action_init=0.5e0*x*x;
    
    dx = (double)rand()/RAND_MAX;
    dx=(dx-0.5e0)*step_size*2e0;
    x=x+dx;
    
    action_fin=0.5e0*x*x;
    /*******************/
    /* Metropolis test */
    /*******************/
    metropolis = (double)rand()/RAND_MAX;    
    if(exp(action_init-action_fin) > metropolis)
      /* accept */
      naccept=naccept+1;
    else 
      /* reject */
      x=backup_x;
    /***************/
    /* data output */
    /***************/	
    printf("%f\n",x);}
}

