module functions
use iso_fortran_env
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
      x(i) = u(2,i)+d*dt*(u(1,i+1)+u(1,i-1)-2*u(1,i)+para*((u(1,i+1)-u(1,i))**2-(u(1,i)-u(1,i-1))**2))
    end do
    !$omp end parallel do
    x(N+2) = 0.0_dp

  end function

  function eLc(u,dt,N,para) result(x)
    !Operator corresponding to the corrector term C=[[A,B],B]
    integer, parameter :: dp=selected_real_kind(15,307)
    integer, intent(in) :: N
    real(dp), intent(in) :: u(:,:), para, dt
    integer :: i
    real(dp) ::  c, x(N+2)

    c=-(dt**3)*(2.0_dp-sqrt(3.0_dp))/48.0_dp

    x(1) = 0.0_dp
    i = 2
    x(i) = u(2,i)+c*(2.0_dp)*((2.0_dp*u(1,i+1)-u(1,i+2)-u(1,i)+para*((u(1,i+1)-u(1,i))**(2)-(u(1,i+2)-u(1,i+1))**(2))) &
               *( 1+para*(2)*(u(1,i+1)-u(1,i)) ) &
               + ( -2.0_dp*u(1,i)+u(1,i+1)+u(1,i-1)+para*( -(u(1,i)-u(1,i-1))**(2) + (u(1,i+1)-u(1,i))**(2) ) ) &
               *( 2.0_dp+para*(2)*( (u(1,i)-u(1,i-1)) + (u(1,i+1)-u(1,i)) ) ) ) 
    !$omp parallel do
    do i=3,N
      x(i) = u(2,i)+c*(2.0_dp)*((2.0_dp*u(1,i+1)-u(1,i+2)-u(1,i)+para*((u(1,i+1)-u(1,i))**(2)-(u(1,i+2)-u(1,i+1))**(2))) &
                 *( 1+para*(2)*(u(1,i+1)-u(1,i)) ) &
                 + ( -2.0_dp*u(1,i)+u(1,i+1)+u(1,i-1)+para*( -(u(1,i)-u(1,i-1))**(2) + (u(1,i+1)-u(1,i))**(2) ) ) &
                 *( 2.0_dp+para*(2)*( (u(1,i)-u(1,i-1)) + (u(1,i+1)-u(1,i)) ) ) & 
                 + ( 2.0_dp*u(1,i-1)-u(1,i)-u(1,i-2)+para*( (u(1,i-1)-u(1,i-2))**(2) - (u(1,i)-u(1,i-1))**(2) ) ) &
                 *( 1.0_dp+para*(2)*(u(1,i)-u(1,i-1)) ) ) 
    end do
    !$omp end parallel do
    i = N+1
    x(i) = u(2,i)+c*(2.0_dp)*(( -2.0_dp*u(1,i)+u(1,i+1)+u(1,i-1)+para*( -(u(1,i)-u(1,i-1))**(2) + (u(1,i+1)-u(1,i))**(2) ) ) &
               *( 2.0_dp+para*(2)*( (u(1,i)-u(1,i-1)) + (u(1,i+1)-u(1,i)) ) ) & 
               + ( 2.0_dp*u(1,i-1)-u(1,i)-u(1,i-2)+para*( (u(1,i-1)-u(1,i-2))**(2) - (u(1,i)-u(1,i-1))**(2) ) ) &
               *( 1.0_dp+para*(2)*(u(1,i)-u(1,i-1)) ) )   
    x(N+2) = 0.0_dp

  
  end function

  function kd(i,N) result(x)
    integer, parameter :: dp=selected_real_kind(15,307)
    integer, intent(in) :: N
	integer(kind=int64), intent(in) :: i
    integer :: x

    if (i==0) then
      x=1
    else if ( i==(2*N+2) ) then
      x=-1
    else
      x=0
    end if


  end function

end module

program fpu
  use functions
  use iso_fortran_env
  implicit none 
  integer, parameter :: dp=selected_real_kind(15,307)
  integer :: N, modes, m
  integer(kind=int64) :: i, j,l,k, it_max,jumps,spacing,it
  real(dp), parameter :: pi=4.0_dp*atan(1.0_dp),c1=(sqrt(3.0_dp)-1.0_dp)/(2.0_dp*sqrt(3.0_dp)), c2=1.0_dp/(sqrt(3.0_dp)),d1=0.5_dp
  real(dp) :: dt=0.1_dp, para, max_time, e, Amp, s
  real(dp), allocatable :: u(:,:), a(:,:), nm_energy(:)
  character(len=255) :: path, fname

  read *, N
  read *, e
  read *, max_time
  read *, m
  read *, spacing
  read *, path
  read *, fname
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
    u(2,i) = 0.0_dp
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
      u(2,(N+2)/2+1) = 0.0_dp
      do i=(N+2)/2+2, N+2
        u(1,i) = u(1,N+3-i)
        u(2,i) = u(2,N+3-i)
      end do  
    end if

  end if




  Amp=(1.0_dp/sin((m*pi/(2*N+2.0_dp))))*sqrt(e/(N+1.0_dp))


  u(1,:) = Amp*u(1,:)
  
  !make folder
  !call system('mkdir '//trim(path)//'/ReissTemp') 
  
  !open file to write to
  print*,path
  print*,fname
  print*, trim(path) // '/ReissTemp/' // fname
  open(1,file= trim(path) // '/ReissTemp/' // fname,status="replace")
  write(1,*) N, dt 
  write(1,*) Amp,m
  write(1,*) para, e

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
	s = (sum(nm_energy(:)/e*log(abs(nm_energy(:)/e))))/log(1.0_dp*CEILING(N/2.0))+1
	!print
	write(1,*) (it-1)*spacing, s
    !eLc
    u(2,:) = eLc(u(:,:),dt,N,para)
    !SABA2
    u(1,:) = eLa(u(:,:),c1,dt,N)
    u(2,:) = eLb(u(:,:),dt,N,para)
    u(1,:) = eLa(u(:,:),c2,dt,N)
    u(2,:) = eLb(u(:,:),dt,N,para)
    u(1,:) = eLa(u(:,:),c1,dt,N)
    !eLc
    u(2,:) = eLc(u(:,:),dt,N,para)
    do k=1,jumps-1 !SABA2C 
    !Second loops exists so not all iterations are stored and kept
      !eLc
      u(2,:) = eLc(u(:,:),dt,N,para)
      !SABA
      u(1,:) = eLa(u(:,:),c1,dt,N)
      u(2,:) = eLb(u(:,:),dt,N,para)
      u(1,:) = eLa(u(:,:),c2,dt,N)
      u(2,:) = eLb(u(:,:),dt,N,para)
      u(1,:) = eLa(u(:,:),c1,dt,N)
      !eLc
      u(2,:) = eLc(u(:,:),dt,N,para)
    end do
  end do
  
  close(1)
  
  !move file to projectnb
  
  call system('mv '//trim(path) // '/ReissTemp/' // fname//' /projectnb2/frgeeeph/FPUT_files/AlphaMeta/'//fname)


  !u(:,1,:) = u(:,it_max,:)
  !dt=-dt
  !t(1) = t(it_max)


end program fpu