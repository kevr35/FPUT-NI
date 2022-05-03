module toda_functions
implicit none
contains


  function eLa(u,c,dt,N) result(x)
    !Operator corresponding to the kinetic part of the Hamiltonian
    integer, parameter :: dp=selected_real_kind(15,307)
    integer, intent(in) :: N
    real(dp), intent(in) :: u(:,:), c, dt
    real(dp) :: x(N+2)
    integer :: i

    x(1) = 0.0_dp
    !$omp parallel do
    do i=2,N+1
      x(i) = u(1,i) + c*dt*u(2,i)
    end do
    !$omp end parallel do
    x(N+2) = 0.0_dp
  end function

  function eLb(u,dt,N,para) result(x)
    !Operator corresponding to the potential part of the Hamiltonian
    integer, parameter :: dp=selected_real_kind(15,307)
    integer, intent(in) :: N
    real(dp), intent(in) :: u(:,:), para, dt
    real(dp) :: d=0.5_dp, x(N+2)
    integer ::  i

    x(1) = 0.0_dp
    !$omp parallel do
    do i=2,N+1
      x(i) = u(2,i)+0.5_dp*d*dt*(-exp(2.0_dp*para*(u(1,i)-u(1,i-1)))+exp(2.0_dp*para*(u(1,i+1)-u(1,i))))/para
    end do
    !$omp end parallel do
    x(N+2) = 0.0_dp
  end function

  
end module

program toda
  use toda_functions
  use iso_fortran_env
  implicit none 
  integer, parameter :: dp=selected_real_kind(15,307)
  integer :: N, modes, m
  integer(kind=int64) :: i, j,l,k, it_max,jumps,spacing,it
  real(dp), parameter :: pi=4.0_dp*atan(1.0_dp),c1=(sqrt(3.0_dp)-1.0_dp)/(2.0_dp*sqrt(3.0_dp)), c2=1.0_dp/(sqrt(3.0_dp)),d1=0.5_dp
  real(dp) :: dt=0.1_dp, para, max_time, e, s, theta
  real(dp), allocatable :: u(:,:), a(:,:), nm_energy(:)
  !character(len=255) :: path, fname

  read *, N
  read *, e
  read *, max_time
  read *, m
  read *, spacing
  read *, theta
  !read *, path
  !read *, fname
  !write(*, '("Enter the value of alpha: ")', &
  !  advance='no')
  !read *, para
  para = 1.0_dp
  jumps = spacing/dt


  it_max = max_time/(spacing)
  !write(*, '("Enter the number of linear modes to be recorded: ")', &
  !  advance='no')
  !read *, modes
  modes = N
  allocate( u(2,N+2),a(2,N), nm_energy(modes))
  
  !Initial Conditions
  do i=1,(N+2)/2
    u(1,i) = sin(m*(i-1)*pi/(N+1))
    u(2,i) = sin(m*(i-1)*pi/(N+1))
  end do



  if (mod(m,2)==0) then !m is even

    if (mod(N,2)==0) then !N is even
      do i=(N+2)/2+1, N+2
        u(1,i) = -u(1,N+3-i)
        u(2,i) = -u(2,N+3-i)
      end do
    else !N is odd
      u(1,(N+2)/2+1) = 0.0_dp
      u(2,(N+2)/2+1) = 0.0_dp
      do i=(N+2)/2+2, N+2
        u(1,i) = -u(1,N+3-i)
        u(2,i) = -u(2,N+3-i)
      end do  
    end if



  else !m is odd

    if (mod(N,2)==0) then !N is even
      do i=(N+2)/2+1, N+2
        u(1,i) = u(1,N+3-i)
        u(2,i) = u(2,N+3-i)
      end do
    else !N is odd
      u(1,(N+2)/2+1) = (-1.0_dp)**((m+3)/(2))
      u(2,(N+2)/2+1) = (-1.0_dp)**((m+3)/(2))
      do i=(N+2)/2+2, N+2
        u(1,i) = u(1,N+3-i)
        u(2,i) = u(2,N+3-i)
      end do  
    end if

  end if
  
 !amplitudes of P and Q 
u(1,:)=cos(theta) * (1.0_dp/sin((m*pi/(2*N+2.0_dp))))*sqrt(e/(N+1.0_dp)) * u(1,:)
u(2,:)=sin(theta) * 2.0_dp*sqrt(e/(N+1.0_dp)) * u(2,:)
  
  !make folder
  !call system('mkdir '//trim(path)//'/ReissTemp') 
  
  !open file to write to
  !open(1,file= trim(path) // '/ReissTemp/' // fname,status="replace")
  print*, N, dt 
  print*, theta,m
  print*, para, e

  !Integration
  do it=1,it_max-1 !SABA2C
	    !Normal modes transformation
	do k=1,N
		a(:,k) = u(:,1) 
		do j=2,N+1
		a(:,k) = a(:,k) + sqrt(2.0_dp/(N+1.0_dp))*u(:,j)*sin(pi*(j-1.0_dp)*k/(N+1.0_dp))
		end do
	end do
	!Normal mode energy
	!$omp parallel do
	do k=1,modes
		nm_energy(k) = (0.5_dp)*(a(2,k))**2+2.0_dp*((sin(pi*k/(2.0_dp*N+2.0_dp)))**2)*(a(1,k))**2
	end do
	!$omp end parallel do
	!entropy
	s = (sum(nm_energy(:)/sum(nm_energy)*log(abs(nm_energy(:)/sum(nm_energy)))))/log(1.0_dp*N)+1
	!print
	print*, (it-1)*spacing, s
    !SABA2
    u(1,:) = eLa(u(:,:),c1,dt,N)
    u(2,:) = eLb(u(:,:),dt,N,para)
    u(1,:) = eLa(u(:,:),c2,dt,N)
    u(2,:) = eLb(u(:,:),dt,N,para)
    u(1,:) = eLa(u(:,:),c1,dt,N)
    do k=1,jumps-1 !SABA2C 
    !Second loops exists so not all iterations are stored and kept
      !SABA
      u(1,:) = eLa(u(:,:),c1,dt,N)
      u(2,:) = eLb(u(:,:),dt,N,para)
      u(1,:) = eLa(u(:,:),c2,dt,N)
      u(2,:) = eLb(u(:,:),dt,N,para)
      u(1,:) = eLa(u(:,:),c1,dt,N)
    end do
  end do
  
  !close(1)
  
  !move file to projectnb
  
  !call system('mv '//trim(path) // '/ReissTemp/' // fname//' /projectnb2/frgeeeph/FPUT_files/AlphaMeta/'//fname)


  !u(:,1,:) = u(:,it_max,:)
  !dt=-dt
  !t(1) = t(it_max)


end program toda