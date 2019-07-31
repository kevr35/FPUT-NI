module functions
implicit none
contains
  function rnd48(seed)

  integer, parameter :: dp=selected_real_kind(15,307)
  integer(dp), parameter :: seed_m=25214903917_dp, seed_a=11_dp
  real(dp), parameter :: twoneg48=0.35527136788005009e-14_dp
  real rnd48
  integer(dp) :: seed

  seed=iand(seed_m*seed+seed_a,281474976710655_dp)
  rnd48=twoneg48*seed

  end function rnd48

  function eLa(u,c,dt) result(x)
    !Operator corresponding to the kinetic part of the Hamiltonian
    integer, parameter :: dp=selected_real_kind(15,307)
    real(dp), intent(in) :: u(:,:), c, dt
    real(dp) :: x(2)

    x(:) = u(1,:) + c*dt*u(2,:)

  end function

  function eLb(u,dt) result(x)
    !Operator corresponding to the potential part of the Hamiltonian
    integer, parameter :: dp=selected_real_kind(15,307)
    real(dp), intent(in) :: u(:,:), dt
    real(dp) :: d=0.5_dp, x(2)

    x(1) = u(2,1) + d*dt*(-u(1,1)-2*u(1,1)*u(1,2))   
    x(2) = u(2,2) + d*dt*(-u(1,2)-(u(1,1))**2+(u(1,2))**2)   

  end function

  function eLc(u,dt) result(x)
    !Operator corresponding to the corrector term C=[[A,B],B]
    integer, parameter :: dp=selected_real_kind(15,307)
    real(dp), intent(in) :: u(:,:),  dt
    real(dp) ::  c, x(2)

    c=-(dt**3)*(2.0_dp-sqrt(3.0_dp))/48.0_dp

    x(1) = u(2,1) + c*(-2*(u(1,1)+2*u(1,1)*u(1,2))*(1+2*u(1,2))-2*(u(1,2)+(u(1,1))**2-(u(1,2))**2)*(2*u(1,1)))
    x(2) = u(2,2) + c*(-2*(u(1,1)+2*u(1,1)*u(1,2))*(2*u(1,1))-2*(u(1,2)+(u(1,1))**2-(u(1,2))**2)*(1-2*u(1,2)))

  
  end function

end module

program henonheiles
  use functions
  implicit none 
  integer, parameter :: dp=selected_real_kind(15,307)
  integer(dp) :: i,k,it_max,seed = 25676
  real(dp), parameter :: c1=(sqrt(3.0_dp)-1.0_dp)/(2.0_dp*sqrt(3.0_dp)), c2=1.0_dp/(sqrt(3.0_dp)),d1=0.5_dp
  real(dp) :: dt=0.1_dp, max_time, energy, p_x, x_x, y_y
  real(dp), allocatable :: u(:,:,:)

  write(*, '("Enter the total energy: ")', &
    advance='no')
  read *, energy


  max_time= 10000
  it_max = max_time/(dt)

  allocate( u(2,it_max,2) )
  open(1,file="henon-heiles",status="replace")
  write(1,*) energy
  !Initial Conditions
  do k = 1,100
    !Selected initial condtions ~uniform(-1.5,1.5)
    10 continue
    p_x = (rnd48(seed)-0.5_dp)*2.0_dp
    x_x = (rnd48(seed)-0.5_dp)*2.0_dp
    y_y = (rnd48(seed)-0.5_dp)*2.0_dp
    if (2.0_dp*energy-p_x*p_x-x_x*x_x-y_y*y_y-(2.0_dp*x_x*x_x*y_y)+((2.0_dp/3.0_dp)*y_y**3).lt.0) then 
      goto 10
    end if

    u(1,1,1) = x_x
    u(1,1,2) = y_y
    u(2,1,1) = p_x
    u(2,1,2) = (-1)**k*sqrt(2.0_dp*energy-p_x*p_x-x_x*x_x-y_y*y_y-(2.0_dp*x_x*x_x*y_y)+((2.0_dp/3.0_dp)*y_y**3))

    !Integration
    do i=1,it_max-1 !SABA2C
      u(1,i+1,:) = u(1,i,:)
      u(2,i+1,:) = eLc(u(:,i,:),dt)
      u(1,i+1,:) = eLa(u(:,i+1,:),c1,dt)
      u(2,i+1,:) = eLb(u(:,i+1,:),dt)
      u(1,i+1,:) = eLa(u(:,i+1,:),c2,dt)
      u(2,i+1,:) = eLb(u(:,i+1,:),dt)
      u(1,i+1,:) = eLa(u(:,i+1,:),c1,dt)
      u(2,i+1,:) = eLc(u(:,i+1,:),dt)
    end do
    !Write data
    do i = 1,it_max
      if (-0.001<u(1,i,1) .and. u(1,i,1)<0.001 .and. u(2,i,1)>0) then
        write(1,*) u(1,i,:), u(2,i,:)
      end if 
    end do
  end do
  close(1)
end program henonheiles