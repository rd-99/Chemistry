Module mod_verlet
    implicit none

    integer npart
    real*8, allocatable :: pos(:,:),vel(:,:),acc(:,:),mass(:)
    real*8 pot, energy, E_1stbond
    real*8 dt,total_time,curr_time
    real*8 dt_min,dt_max
    integer nsteps,nsteps_dt
    real*8 omega,x_0
    real*8 E_avg,Esq_avg,E_std_dev
    integer flag_write, flag_vmd

    contains
    !##############################################

    subroutine setup_parameters
        implicit none
        real*8 mass_atom

        !! reading parameters fron input file
        open(10,file="md_fput.inp")
        read(10,*) npart
        read(10,*) nsteps_dt
        read(10,*) dt_min
        read(10,*) dt_max
        read(10,*) total_time
        read(10,*) mass_atom
        read(10,*) omega
        read(10,*) x_0
        read(10,*) flag_write   !! if set to 1, writes energy as a function of time
        read(10,*) flag_vmd     !! if set to 1, writes vmd output (xyz format)
        close(10)

        !! allocating arrrays dependent on npart
        allocate(pos(npart,3),vel(npart,3),acc(npart,3),mass(npart))
        mass=mass_atom

    end subroutine setup_parameters
    !##############################################

    subroutine compute_std_dev_en
        !! Varies dt from dt_min to dt_max, runs a MD trajectory using velocity-Verlet and calculate the standard dev. of energy
        implicit none
        integer i

        open(100,file="std_en_dt.out")
        do i=1,nsteps_dt
            dt=dt_min+(i-1)*(dt_max-dt_min)/real(nsteps_dt-1)
            call initial_condition
            call evolve

            E_avg=E_avg/real(nsteps)
            Esq_avg=Esq_avg/real(nsteps)
            E_std_dev = dsqrt(Esq_avg-E_avg**2)
            write(100,*) dt, E_std_dev
        enddo
        close(100)

    end subroutine compute_std_dev_en
    !##############################################

    subroutine initial_condition
        implicit none
        integer i

        !! Setting initial position
        pos=0.d0
        do i=1,npart
          pos(i,1)=(i-1)*x_0
        enddo
        pos(1,1)=pos(1,1)-2.d0

        !! Setting initial velocity
        vel = 0.d0

        call compute_pot(pos,pot,acc)
        call compute_energy
        curr_time=0.0

        E_avg=0.d0
        Esq_avg=0.d0

    end subroutine initial_condition
    !##############################################

    subroutine evolve
        implicit none
        integer i

        nsteps=nint(total_time/dt)
        open(10,file="coordinates.out")
        open(11,file="energy.out")
        do i=1,nsteps
            call compute_energy
            if(flag_write==1) call write_output
            if(flag_vmd==1) call write_vmd
            call velocity_verlet
            curr_time=curr_time+dt
        enddo
        close(10); close(11)
    end subroutine evolve
    !##############################################

    subroutine velocity_verlet
        !!Velocity verlet algorithm
        implicit none

        pos=pos+vel*dt+0.5*acc*dt*dt
        vel=vel+0.5*acc*dt
        call compute_pot(pos,pot,acc)
        vel=vel+0.5*acc*dt

    end subroutine velocity_verlet
    !##############################################

    subroutine compute_pot(x,pot,acc)
        !! calculates potential (pot) and acceleration (acc) at the given positions (x)
        implicit none

        double precision,intent(in)::x(npart,3)
        double precision,intent(out)::pot,acc(npart,3)
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

    subroutine write_output
        implicit none

        write(10,*) curr_time,pos(1,1),vel(1,1)
        write(11,*) curr_time,energy

    end subroutine write_output
    !##############################################

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

    subroutine compute_energy
        !! calculate total energy, energy of 1st bond, average energy and average of energy squared
        implicit none
        integer i

        energy=pot
        do i=1,npart
          energy=energy+0.5*mass(i)*sum(vel(i,:)*vel(i,:))
        enddo
        E_1stbond = 0.5*mass(1)*(vel(2,1)-vel(1,1))**2 + 0.5*mass(1)*omega**2*(pos(2,1)-pos(1,1)-x_0)**2

        E_avg=E_avg+energy
        Esq_avg=Esq_avg+energy**2

    end subroutine compute_energy
    !##############################################

end module mod_verlet
