!##############################################################################
!######              Gaussian Integral by Metropolis                  #########
!##############################################################################
program Gaussian

  implicit none
  !---------------------------------
  integer iter, niter
  parameter(niter=10000)
  integer naccept
  double precision step_size
  parameter(step_size=4d0)
  double precision x, backup_x, dx,sumxx
  double precision action_init, action_fin
  double precision metropolis
  
  call random_seed()        

  !*************************************
  !*** Set the initial configuration ***
  !*************************************
  x=0d0
  sumxx=0d0!sum of x^2, to calculate <x^2>
  !****************************
  !*** Make the output file ***
  !****************************
  open(unit=10,status='REPLACE',file='output.txt',action='WRITE')
  write(10,*) "# iteration, x, <x^2>, acceptance rate"
  write(10,*)'#------------------------------------------------'
  !*****************
  !*** Main part ***
  !*****************  
  naccept=0 !counter for the number of acceptance
  do iter=1,niter    
     backup_x=x
     action_init=0.5d0*x**2d0

     call random_number(dx)
     dx=(dx-0.5d0)*step_size*2d0
     x=x+dx
     
     action_fin=0.5d0*x**2d0
     !***********************
     !*** Metropolis test ***
     !***********************
     call random_number(metropolis)
     if(dexp(action_init-action_fin) > metropolis)then
        !accept
        naccept=naccept+1
     else 
        !reject
        x=backup_x
     end if
     sumxx=sumxx+x*x
     !*******************
     !*** data output ***
     !*******************
     write(10,*)iter,x,sumxx/dble(iter),dble(naccept)/dble(iter)    
     write(*,*)iter,x,sumxx/dble(iter),dble(naccept)/dble(iter) 
 end do

  close(10)

end program Gaussian


