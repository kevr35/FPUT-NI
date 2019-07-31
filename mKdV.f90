module mkdv_functions
implicit none
contains

  function d(phi,a,N) result(x)
    integer, parameter :: dp=selected_real_kind(15,307)
    real(dp), intent(in) :: phi(:), a
    integer, intent(in) :: N
    real(dp) :: x(2*N)

    x = (cshift(phi,1)-cshift(phi,-1)) / (2.0_dp*a)

  end function

  function mkdv(phi,a,dt,beta,N) result(x)
    integer, parameter :: dp=selected_real_kind(15,307)
    real(dp), intent(in) :: phi(:), a, beta, dt
    integer, intent(in) :: N
    integer ::  i
    real(dp) :: x(2*N)

    x = phi+dt*(-3.0_dp*beta*phi**(2)*d(phi,a,N) -(1.0_dp/12.0_dp) * d(d(d(phi,a,N),a,N),a,N) ) / (2.0_dp * (N+1)**3)

  end function

end module




program mkdv
  use mkdv_functions
  implicit none 
  integer, parameter :: dp=selected_real_kind(15,307)
  integer :: N, i, j, it_max, jumps
  real(dp), parameter :: pi=4.0_dp*atan(1.0_dp), L = 1.0_dp
  real(dp) :: dt=0.1_dp, k, beta, max_time, S, Amp, a
  real(dp), allocatable :: u(:,:), t(:)

  k = pi/L
  jumps = 5/dt  

  write(*, '("Enter the number of active particles on the lattice: ")', &
    advance='no')
  read *, N
  write(*, '("Enter the value of S: ")', &
    advance='no')
  read *, S

  beta = sign(1.0_dp, S)
  Amp = (2.0_dp/pi)*sqrt( abs(S)/abs(beta) )
  a = 1.0_dp/(N+1)**3

  if (-10.0_dp<=S<=4.0_dp) then
    max_time = ( 864.0_dp*S**2 - 5376.0_dp * pi**2 * S + 4096.0_dp*pi**4 ) / ( 405.0_dp * S**3 +4104 pi**2 S**2 - 4992 * pi**(4) * S + 2048*pi**6)
  elseif (S<-10.0_dp) then
    max_time = 0.3078_dp/sqrt(S)
  elseif (S>4.0_dp) then
    max_time = 0.5862_dp/sqrt(S)
  end if

  max_time = 1.25_dp*max_time
  it_max = max_time/(dt*jumps)

  allocate( u(it_max,2*N), t(it_max)  )

  !Initial Conditions
  t(1) = 0.0_dp
  do i=1,2*N+1
    u(1,i) = a*k*cos((i-1)*pi/(N+1))/2.0_dp
  end do

  !Integration
  do i=1,it_max-1 !SABA2C
    t(i+1) = t(i) + dt
    u(i+1,:) = mkdv( u(i,:) , a , dt , beta , N)


    do j=1,jumps-1 !SABA2C 
    !Second loops exists so not all iterations are stored and kept
    t(i+1) = t(i) + dt
    u(i+1,:) = mkdv( u(i,:) , a , dt , beta , N)
    end do
  end do

  print*,'integration completed'

  print*,'writing data to file...'
  open(1,file="mKdV",status="replace")
  write(1,*) N, beta, S
  do i = 1,it_max
    write(1,*) t(i), u(i,:)
  end do
  close(1)
  print*,'done'


end program mkdv