!############################################################################
!######              Matrix Integral by Metropolis                  #########
!######                  written by M. Hanada                       #########
!############################################################################
program phi4

  implicit none
  !---------------------------------
  integer nmat
  parameter(nmat=100)
  integer ninit
  parameter(ninit=0)!ninit=1 -> new config; ninit=0 -> old config
  integer iter,niter
  parameter(niter=10000)
  integer ntau
  parameter(ntau=20)
  double precision dtau
  parameter(dtau=0.005d0)
  integer naccept
  double complex phi(1:NMAT,1:NMAT),backup_phi(1:NMAT,1:NMAT)
  double precision ham_init,ham_fin,action,sum_action
  double precision tr_phi,tr_phi2
  double precision metropolis

  open(unit=10,status='REPLACE',file='matrix-HMC.txt',action='WRITE')
  !*************************************
  !*** Set the initial configuration ***
  !*************************************
  call pre_random
  if(ninit.EQ.1)then
     phi=(0d0,0d0)
  else if(ninit.EQ.0)then
     open(UNIT=22, File ='config.dat', STATUS = "OLD", ACTION = "READ")
     read(22,*) phi
     close(22)
  end if
  sum_action=0d0
  !*****************
  !*** Main part ***
  !*****************  
  naccept=0 !counter for the number of acceptance
  do iter=1,niter

     backup_phi=phi     
     call Molecular_Dynamics(nmat,phi,dtau,ntau,ham_init,ham_fin)
     !***********************
     !*** Metropolis test ***
     !***********************
     call random_number(metropolis)
     if(dexp(ham_init-ham_fin) > metropolis)then
        !accept
        naccept=naccept+1
     else 
        !reject
        phi=backup_phi
     end if
     !*******************
     !*** data output ***
     !*******************
     call calc_action(nmat,phi,action)
     sum_action=sum_action+action
     write(*,*)iter,action/dble(nmat*nmat),sum_action/dble(iter)/dble(nmat*nmat),dble(naccept)/dble(iter)

     write(10,*)iter,action/dble(nmat*nmat),sum_action/dble(iter)/dble(nmat*nmat),dble(naccept)/dble(iter)

  end do

  close(10)

  open(UNIT = 22, File = 'config.dat', STATUS = "REPLACE", ACTION = "WRITE")
  write(22,*) phi
  close(22)
     
  
  !  call Hamiltonian_conservation_test(nmat,phi)

  


end program Phi4
!###############################################
SUBROUTINE calc_action(nmat,phi,action)

  implicit none

  integer nmat
  double precision action
  double complex phi(1:NMAT,1:NMAT)
  double complex phi2(1:NMAT,1:NMAT)
  integer imat,jmat,kmat
  !*** phi2=phi*phi ***
  phi2=(0d0,0d0)
  do imat=1,nmat
     do jmat=1,nmat
        do kmat=1,nmat
           phi2(imat,jmat)=phi2(imat,jmat)+phi(imat,kmat)*phi(kmat,jmat)
        end do
     end do
  end do

  action=0d0
  !*** Tr phi^2 term ***
  do imat=1,nmat
     action=action+0.5d0*dble(phi2(imat,imat))
  end do
  !*** Tr phi^4 term ***
  do imat=1,nmat
     do jmat=1,nmat
        action=action+0.25d0*dble(phi2(imat,jmat)*phi2(jmat,imat))
     end do
  end do
  !*** overall normalization ***
  action=action*dble(nmat)

  return
  
END SUBROUTINE calc_action
!###############################################     
SUBROUTINE calc_hamiltonian(nmat,phi,P_phi,ham)

  implicit none

  integer nmat
  double precision action,ham
  double complex phi(1:NMAT,1:NMAT)
  double complex P_phi(1:NMAT,1:NMAT)
  integer imat,jmat

  call calc_action(nmat,phi,action)

  ham=action
  do imat=1,nmat
     do jmat=1,nmat
        ham=ham+0.5d0*dble(P_phi(imat,jmat)*P_phi(jmat,imat))
     end do
  end do
  
  return

END SUBROUTINE calc_hamiltonian
!###################################
subroutine Molecular_Dynamics(nmat,phi,dtau,ntau,ham_init,ham_fin)
  
  implicit none
  
  integer nmat
  integer ntau
  double precision dtau
  double precision r1,r2
  double precision ham_init,ham_fin
  double complex phi(1:NMAT,1:NMAT)
  double complex P_phi(1:NMAT,1:NMAT)
  double complex delh(1:NMAT,1:NMAT)
  integer imat,jmat,step
  !*** randomly generate auxiliary momenta ***
  do imat=1,nmat-1
     do jmat=imat+1,nmat
        call BoxMuller(r1,r2)
        P_phi(imat,jmat)=dcmplx(r1/dsqrt(2d0))+dcmplx(r2/dsqrt(2d0))*(0D0,1D0)
        P_phi(jmat,imat)=dcmplx(r1/dsqrt(2d0))-dcmplx(r2/dsqrt(2d0))*(0D0,1D0)
     end do
  end do
  do imat=1,nmat
     call BoxMuller(r1,r2)
     P_phi(imat,imat)=dcmplx(r1)
  end do
  !*** calculate Hamiltonian ***
  call calc_hamiltonian(nmat,phi,P_phi,ham_init)
  !*** first step of leap frog ***
  phi=phi+P_phi*dcmplx(0.5d0*dtau)
  !*** 2nd, ..., Ntau-th steps ***
  step=1
  do while (step.LT.ntau)
     step=step+1
     call calc_force(delh,phi,nmat)
     P_phi=P_phi-delh*dtau  
     phi=phi+P_phi*dcmplx(dtau)
  end do
  !*** last step of leap frog ***
  call calc_force(delh,phi,nmat)
  P_phi=P_phi-delh*dtau
  phi=phi+P_phi*dcmplx(0.5d0*dtau)
  !*** calculate Hamiltonian ***
  call calc_hamiltonian(nmat,phi,P_phi,ham_fin)
  
  return
  
END subroutine Molecular_Dynamics
!****************************************************************
SUBROUTINE BoxMuller(p,q)  
  implicit none 
  !***** output *****
  doubleprecision p,q
  !******************
  double precision r,s,Pi

  Pi=2d0*DASIN(1d0)
  !uniform random numbers between 0 and 1
  call random_number(r)
  call random_number(s)
  !Gaussian random numbers, 
  !with weights proportional to e^{-p^2/2} and e^{-q^2/2}
  p=dsqrt(-2d0*dlog(r))*DSIN(2d0*Pi*s)
  q=dsqrt(-2d0*dlog(r))*DCOS(2d0*Pi*s)

  return

END SUBROUTINE BoxMuller
!##############################
subroutine calc_force(delh,phi,nmat)

  implicit none

  integer nmat
  double complex phi(1:NMAT,1:NMAT),phi2(1:NMAT,1:NMAT),phi3(1:NMAT,1:NMAT)
  double complex delh(1:NMAT,1:NMAT)
  
  integer imat,jmat,kmat
  !*** phi2=phi*phi, phi3=phi*phi*phi ***  
  phi2=(0d0,0d0)
  phi3=(0d0,0d0)
  do imat=1,nmat
     do jmat=1,nmat
        do kmat=1,nmat
           phi2(imat,jmat)=phi2(imat,jmat)+phi(imat,kmat)*phi(kmat,jmat)
        end do
     end do
  end do
  do imat=1,nmat
     do jmat=1,nmat
        do kmat=1,nmat
           phi3(imat,jmat)=phi3(imat,jmat)+phi2(imat,kmat)*phi(kmat,jmat)
        end do
     end do
  end do
  !*** delh=dH/dphi *** 
  delh=phi+phi3
  delh=delh*dcmplx(nmat)

  return

END subroutine Calc_Force
!###################################
subroutine Hamiltonian_conservation_test(nmat,phi)
  
  implicit none
  
  integer nmat
  integer ntau
  double precision dtau
  double precision r1,r2
  double precision ham_init,ham_fin
  double complex phi(1:NMAT,1:NMAT),backup_phi(1:NMAT,1:NMAT)
  double complex P_phi(1:NMAT,1:NMAT),backup_P(1:NMAT,1:NMAT)
  double complex delh(1:NMAT,1:NMAT)
  
  integer imat,jmat,step
  
  do imat=1,nmat-1
     do jmat=imat+1,nmat
        call BoxMuller(r1,r2)
        P_phi(imat,jmat)=&
             dcmplx(r1/dsqrt(2d0))&
             &+dcmplx(r2/dsqrt(2d0))*(0D0,1D0)
        P_phi(jmat,imat)=&
             dcmplx(r1/dsqrt(2d0))&
             &-dcmplx(r2/dsqrt(2d0))*(0D0,1D0)
     end do
  end do
  do imat=1,nmat
     call BoxMuller(r1,r2)
     P_phi(imat,imat)=dcmplx(r1)
  end do

  backup_phi=phi
  backup_P=P_phi

  do ntau=3,1000
     dtau=0.1d0/dble(ntau)
     call calc_hamiltonian(nmat,phi,P_phi,ham_init)
     !*******************************
     !*** first step of leap frog ***                       
     !*******************************
     phi=phi+P_phi*dcmplx(0.5d0*dtau)
     
     step=1
     do while (step.LT.ntau)
        step=step+1
        call calc_force(delh,phi,nmat)
        P_phi=P_phi-delh*dtau  
        phi=phi+P_phi*dcmplx(dtau)
     end do
     !*****************
     !*** last step ***
     !*****************
     call calc_force(delh,phi,nmat)
     P_phi=P_phi-delh*dtau
     phi=phi+P_phi*dcmplx(0.5d0*dtau)
     
     call calc_hamiltonian(nmat,phi,P_phi,ham_fin)

     write(*,*)ntau,ham_fin-ham_init

     P_phi=backup_P
     phi=backup_phi

  end do
  
  return
  
END subroutine Hamiltonian_Conservation_Test

!*** Not my original. Got from somebody. ***
subroutine pre_random
  implicit none
  integer::seedsize,c
  integer,allocatable::seed(:)
 
  !In fortran90, seed is array.
  ! To get seedsize, use below.
  call random_seed(size=seedsize)

  !Allocate seed array.
  allocate(seed(1:seedsize))
 
  !Get system time.
  call system_clock(count=c)

  !Substitute "seed" using system time.
  seed=c

  !Set "seed" to produce random number obey to system time.
  call random_seed(put=seed)
 
  return
end subroutine pre_random
