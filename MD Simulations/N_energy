program md_FPUT

  implicit none  !! Critical to begin with this statement.

  real*8 pi
  real*8 dt,total_time
  integer,parameter :: npart=10
  real*8,dimension(npart) :: mass
  real*8,dimension(npart,3) :: pos,vel,acc
  real*8 omega,x_0
  real*8 tim,pot,energy
  real*8 E_1stbond
  integer nsteps

  call setup_parameters       !! Step 1 - setting parameters
  call initial_conditions     !! Step 2 - initial conditions
  call evolve                 !! Step 3 - Evolution

  !! Main program finishes here.
  !! Defining functions and subroutine that can use the variables used so far without any definition
  contains
  !-----------------------------------------------------------------  

  !! Step 0: sets up value of important variables such as mass, potential parameters, dt, ...
  subroutine setup_parameters
    !! Note that the variables assigned here have been defined in the main program. They should not be re-defined here. The value set here will be persistent in the main code.
    implicit none

    pi=3.14159265359

    mass  = 21874.8d0             !! = 12.d0 * 1822.d0
    omega = 0.006825d0            !! = 1500.d0 * 4.55d-6   
    x_0   = 2.83d0                !! = 1.5 / 0.529
    
    dt = 0.1 * 2*pi/omega
    total_time = 20.0 * 2*pi/omega

    nsteps=nint(total_time/dt)  !## The number of time steps. Note use of int here.

  end subroutine setup_parameters
  !##############################################

  !! Step 2: subroutine to set initial conditions, to be run before start of the evolution
  subroutine initial_conditions
    implicit none
    integer i

    !! All particles at equilibrium, except first particle is shifted to left by 2 a.u.
    pos=0.d0
    do i=1,npart
      pos(i,1)=(i-1)*x_0
    enddo
    pos(1,1)=pos(1,1)-2.d0

    !! All initial velocities=0
    vel = 0.d0

    call compute_pot(pos,pot,acc)
    call compute_energy
    tim=0.0

  end subroutine initial_conditions
  !##############################################

  !! Step 3: Evolves MD
  subroutine evolve
    implicit none
    integer i

    open(10,file='coordinates.out')
    open(11,file="energy.out")
    open(13,file="vmd.xyz")
    do i=1,nsteps
      call write_output
      call write_vmd
      call evolve_1step
    end do
    close(10);close(11);close(13)

  end subroutine evolve
  !##############################################

  subroutine write_output
    implicit none

    write(10,*) tim,pos(1,1),vel(1,1)
    write(11,*) tim,energy,E_1stbond

  end subroutine write_output
  !##############################################

  !! Writes output compatible with VMD
  subroutine write_vmd
    implicit none
    integer i

    write(13,*) npart
    write(13,*)
    do i=1,npart
      write(13,*) "C",pos(i,:)
    enddo

  end subroutine write_vmd
  !##############################################

  !! Evolution by one timestep dt using Velocity-Verlet method
  subroutine evolve_1step
    !! Velocity-Verlet - https://en.wikipedia.org/wiki/Verlet_integration#Velocity_Verlet
    implicit none

    pos=pos+vel*dt+0.5*acc*dt*dt
    vel=vel+0.5*acc*dt
    call compute_pot(pos,pot,acc)         !## getting updated acceleration and potential
    vel=vel+0.5*acc*dt
    call compute_energy
    tim=tim+dt                       !## Updating current time

  end subroutine evolve_1step
  !##############################################

  !! Calculates potential (pot) and acceleration (acc) for the input position (x)
  subroutine compute_pot(x,pot,acc)
    implicit none !! ALWAYS START WITH THIS

    real*8,intent(in)::x(npart,3)  !! intent(in) tells that this is an input whose value cannot change within the subroutine
    real*8,intent(out)::pot,acc(npart,3) !! intent(out) tells that pot,acc are outputs of this subroutine and should be assigned a value.
    integer i

    pot=0.d0
    acc=0.d0
    do i=1,npart-1
      pot=pot+0.5*mass(i)*omega**2*(x(i+1,1)-x(i,1)-x_0)**2
      acc(i,1)=acc(i,1)-1.d0/mass(i) * (-mass(i)*omega**2*(x(i+1,1)-x(i,1)-x_0))
      acc(i+1,1)=acc(i+1,1)-1.d0/mass(i) * (mass(i)*omega**2*(x(i+1,1)-x(i,1)-x_0))
    enddo

  end subroutine compute_pot
  !##############################################

  !! Calcualtes total energy (energy) and energy of first bond (E_1stbond)
  !! compute_pot should be called before this to get the value of potential (pot)
  subroutine compute_energy
    implicit none
    integer i

    energy=pot
    do i=1,npart
      energy=energy+0.5*mass(i)*sum(vel(i,:)*vel(i,:))
    enddo
    E_1stbond = 0.5*mass(1)*(vel(2,1)-vel(1,1))**2 + 0.5*mass(1)*omega**2*(pos(2,1)-pos(1,1)-x_0)**2

  end subroutine compute_energy
  !##############################################

end program md_FPUT
