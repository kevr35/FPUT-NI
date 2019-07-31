module functions
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
    integer, intent(in) :: i, N
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
  implicit none 
  integer, parameter :: dp=selected_real_kind(15,307)
  integer :: N, i, k,j, m=1, it_max, modes, jumps
  real(dp), parameter :: pi=4.0_dp*atan(1.0_dp),c1=(sqrt(3.0_dp)-1.0_dp)/(2.0_dp*sqrt(3.0_dp)), c2=1.0_dp/(sqrt(3.0_dp)),d1=0.5_dp
  real(dp) :: dt=0.1_dp, para, max_time, e, Amp, e1=0.0_dp, e2=0.0_dp
  real(dp), allocatable :: u(:,:,:), a(:,:,:), nm_energy(:,:), T(:), energy(:), p(:,:),ef(:)

  read *, N

  read *, e

  para = 1.0_dp
  jumps = 10/dt

  read *, max_time
  it_max = max_time/(dt*jumps)
  !write(*, '("Enter the number of linear modes to be recorded: ")', &
  !  advance='no')
  !read *, modes

  allocate( u(2,it_max,N+2), energy(it_max), t(it_max), p(it_max,N+2),ef(it_max)  )

  !Initial Conditions
  t(1) = 0.0_dp
  do i=1,(N+2)/2
    u(1,1,i) = sin(m*(i-1)*pi/(N+1))
    u(2,1,i) = 0.0_dp
  end do



  if (mod(m,2)==0) then !m is even

    if (mod(N,2)==0) then !N is even
      do i=(N+2)/2+1, N+2
        u(1,1,i) = -u(1,1,N+3-i)
        u(2,1,i) = -u(2,1,N+3-i)
      end do
    else !N is odd
      u(1,1,(N+2)/2+1) = 0.0_dp
      u(2,1,(N+2)/2+1) = 0.0_dp
      do i=(N+2)/2+2, N+2
        u(1,1,i) = -u(1,1,N+3-i)
        u(2,1,i) = -u(2,1,N+3-i)
      end do  
    end if



  else !m is odd

    if (mod(N,2)==0) then !N is even
      do i=(N+2)/2+1, N+2
        u(1,1,i) = u(1,1,N+3-i)
        u(2,1,i) = u(2,1,N+3-i)
      end do
    else !N is odd
      u(1,1,(N+2)/2+1) = (-1.0_dp)**((m+3)/(2))
      u(2,1,(N+2)/2+1) = 0.0_dp
      do i=(N+2)/2+2, N+2
        u(1,1,i) = u(1,1,N+3-i)
        u(2,1,i) = u(2,1,N+3-i)
      end do  
    end if

  end if




  u(1,1,:) = (1.0_dp/sin((m*pi/(2*N+2.0_dp))))*sqrt(e/(N+1.0_dp))*u(1,1,:)

  !Integration
  do i=1,it_max-1 !SABA2C
    t(i+1) = t(i) + dt
    !eLc
    u(1,i+1,:) = u(1,i,:)
    u(2,i+1,:) = eLc(u(:,i,:),dt,N,para)
    !SABA2
    u(1,i+1,:) = eLa(u(:,i+1,:),c1,dt,N)
    u(2,i+1,:) = eLb(u(:,i+1,:),dt,N,para)
    u(1,i+1,:) = eLa(u(:,i+1,:),c2,dt,N)
    u(2,i+1,:) = eLb(u(:,i+1,:),dt,N,para)
    u(1,i+1,:) = eLa(u(:,i+1,:),c1,dt,N)
    !eLc
    u(2,i+1,:) = eLc(u(:,i+1,:),dt,N,para)
    do k=1,jumps-1 !SABA2C 
    !Second loops exists so not all iterations are stored and kept
      t(i+1) = t(i+1) + dt
      !eLc
      u(2,i+1,:) = eLc(u(:,i+1,:),dt,N,para)
      !SABA
      u(1,i+1,:) = eLa(u(:,i+1,:),c1,dt,N)
      u(2,i+1,:) = eLb(u(:,i+1,:),dt,N,para)
      u(1,i+1,:) = eLa(u(:,i+1,:),c2,dt,N)
      u(2,i+1,:) = eLb(u(:,i+1,:),dt,N,para)
      u(1,i+1,:) = eLa(u(:,i+1,:),c1,dt,N)
      !eLc
      u(2,i+1,:) = eLc(u(:,i+1,:),dt,N,para)
    end do
  end do

  p = cshift(u(1,:,:),1,dim=2)
  energy = 0.5_dp*u(2,:,1)**2+0.5_dp*(p(:,1)-u(1,:,1))**2+(para/3.0_dp)*(p(:,1)-u(1,:,1))**(3)
  do i=2,N+1
      energy = energy+0.5_dp*u(2,:,i)**2+0.5_dp*(p(:,i)-u(1,:,i))**2+(para/3.0_dp)*(p(:,i)-u(1,:,i))**(3)
  end do

  do i=1,it_max
     ef(i) = abs(energy(i)-energy(1))/energy(1)
  end do

  print*, N, 3, para, e
  do i = 1,it_max
    print*, t(i), ef(i)
  end do



end program fpu
