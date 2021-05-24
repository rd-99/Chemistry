Module mod_verlet
    implicit none

    integer npart
    real*8, allocatable :: pos(:,:),vel(:,:),acc(:,:),mass(:),mass_inv(:),pos_t0(:,:)
    real*8 pot, energy,r_c_pot,kb,pi,dr2
    real*8 total_time,curr_time,nstep_reset,dt
    integer nsteps
    real*8 epsilon,sigma,sigma6,eps4,eps48
    real*8 density,temperature,box_len,reset_temp_time
    real*8 E_avg,Esq_avg,E_std_dev,r_csq
    integer flag_write, flag_vmd

    contains
    !##############################################

    ! Step 1: setting up parameters and allocating required arrays
    ! Reads from the input file md_LJ_fluid.inp
    subroutine setup_parameters
        implicit none
        real*8 mass_atom

        open(10,file="md_LJ_fluid.inp")
        read(10,*) npart
        read(10,*) dt
        read(10,*) total_time
        read(10,*) mass_atom
        read(10,*) sigma
        read(10,*) epsilon
        read(10,*) density
        read(10,*) temperature
        read(10,*) flag_write
        read(10,*) flag_vmd
        read(10,*) reset_temp_time
        close(10)

        allocate(pos(npart,3),pos_t0(npart,3),vel(npart,3),acc(npart,3),mass(npart),mass_inv(npart))
        mass=mass_atom
        mass_inv= 1.d0/mass

        sigma6=sigma**6
        eps4=epsilon*4
        eps48=epsilon*48
        r_csq = (2.5d0 * sigma)**2
        r_c_pot = ( (1/2.5d0)**12 - (1/2.5d0)**6 )
        kb = 1.380649d0 * 10**(-23)
        pi = 3.142

    end subroutine setup_parameters
    !##############################################

    !! Step 2: setting initial conditions
    !! initial position set to a simple cubic lattice
    !! initial velocities=0
    subroutine initial_condition
        implicit none
        integer i,j,k
        integer npart_side,part_num
        real*8 lattice_len

        npart_side=nint(1.d0*npart**(1/3.d0))   !! no. of particles along each direction
        box_len=(npart/density)**(1/3.d0)       !! length of simulation box
        lattice_len = box_len/real(npart_side)  !! length of the cubic unit cell
        nstep_reset = nint(reset_temp_time/dt)
        write(*,*) "Initial conditions: placing particles in simple cubic lattice:"
        write(*,*) "total number of particles = ",npart
        write(*,*) "particles along each side = ",npart_side
        write(*,*) "box length = ",box_len
        write(*,*) "lattice length = ",lattice_len

        pos=0.d0
        call thermal_velocities(temperature)
        !!vel=  0.d0

        !! putting particles on a simple cubic lattice
        part_num=1
        do i=1,npart_side
            do j=1,npart_side
                do k=1,npart_side
                    pos(part_num,1)=(i-1)*lattice_len
                    pos(part_num,2)=(j-1)*lattice_len
                    pos(part_num,3)=(k-1)*lattice_len
                    part_num=part_num+1
                enddo
            enddo
        enddo

        !! Writing initial conditions in a vmd compatible xyz file
        open(13,file="init_cond.xyz")
        call write_vmd
        close(13)

        call compute_pot
        call compute_energy
        curr_time=0d0
        E_avg=0.d0
        Esq_avg=0.d0

    end subroutine initial_condition
    !##############################################

    !! Step 3: Evolution using velocity-Verlet
    !! Writes energy output if flag_write is set to 1 in the input file
    !! Writes vmd output if flag_vmd is set to 1 in the input file
    subroutine evolve
        implicit none
        integer i

        nsteps=nint(total_time/dt)
        open(10,file="coordinates.out")
        open(11,file="energy.out")
        open(13,file="vmd.xyz")
        open(14,file = "temperature.out")
        open(15,file = "dr2_vs_t.out")

        do i=1,int(nsteps*2/5)
            call compute_energy
            call get_inside_the_box
            call cal_tempreture
            if(mod(i,20)==1.and.flag_write==1) call write_output
            if(mod(i,100)==1.and.flag_vmd==1) call write_vmd
            if(mod(curr_time,nstep_reset) == 0) call thermal_velocities(0.71d0)
            call velocity_verlet
            curr_time=curr_time+dt
        enddo

        pos_t0 = pos(:,:)
        do i=int(nsteps*2/5),nsteps
            call compute_energy
            call get_inside_the_box
            call cal_tempreture
            if(mod(i,20)==1.and.flag_write==1) call write_output
            if(mod(i,5)==1.and.flag_write==1) call write_dr2_output
            if(mod(i,100)==1.and.flag_vmd==1) call write_vmd
            call velocity_verlet
            call cal_dr2
            curr_time=curr_time+dt
        enddo

        E_avg=E_avg/real(nsteps)
        Esq_avg=Esq_avg/real(nsteps)
        E_std_dev = dsqrt(Esq_avg-E_avg**2)
       !! write(11,*) "dt (a.u.) standard deviation of E (a.u)",dt, E_std_dev
        close(10); close(11); close(13)



    end subroutine evolve
    !##############################################

    !! velocity verlet algorithm
    !! also calls for computing potential and acceleration
    subroutine velocity_verlet
        implicit none

        pos=pos+vel*dt+0.5*acc*dt*dt
        vel=vel+0.5*acc*dt
        call compute_pot
        vel=vel+0.5*acc*dt

    end subroutine velocity_verlet
    !##############################################

    !! calculates pairwi'se LJ potential and acceleration
    subroutine compute_pot
        implicit none
        integer i,j
        real*8 rijsq,vec(3),V_LJ,dVLJ_drij
        real*8 dpot_dx(3)

        pot=0.d0
        acc=0.d0
        do i=1,npart-1
            do j=i+1,npart
                call distance(i,j,rijsq,vec)
                call LJ_potential(rijsq,V_LJ,dVLJ_drij)
                if (rijsq < r_csq) then
                    pot=pot+V_LJ
                    dpot_dx=dVLJ_drij*vec
                    acc(i,:)=acc(i,:)-mass_inv(i) * (dpot_dx)
                    acc(j,:)=acc(j,:)+mass_inv(j) * (dpot_dx)
                endif
            enddo
        enddo

    end subroutine compute_pot
    !##############################################

    !! calculates the distance squared between atoms i and j
    subroutine distance(i,j,rijsq,vec)
        !! input: i and j, that represent atoms i and j
        !! Output: vec(3) = pos(i,:)-pos(j,:) (vector from atom j to i)
        !! Output: rijsq: square of the distance between atoms i and j
        implicit none
        integer,intent(in) :: i,j
        real*8,intent(out) :: rijsq,vec(3)
        integer k
       !! call get_inside_the_box
        do k = 1,3
        vec(k) = pos(i,k)-pos(j,k)
        !!call get_inside_the_box
        if ( vec(k) > (0.5d0 * box_len) ) then
             vec(k) = vec(k) - box_len
        else if  (vec(k) < (-0.5d0 * box_len)) then
             vec(k) = vec(k) + box_len
        end if
        enddo
        rijsq = (sum(vec*vec))

    end subroutine distance
    !##############################################

    !! calculates LJ potential and 1/rij*dVLJ_dx
    !! Optimized
    subroutine LJ_potential(xsq,V_LJ,dVLJ_dx)
        !! Takes xsq=x_squares as the input
        !! Calculates V_LJ=LJ potential and dVLJ_dx=(1/x).derivative of potential
        implicit none
        real*8,intent(in)::xsq
        real*8,intent(out)::V_LJ,dVLJ_dx
        real*8 sig6_x6

        sig6_x6=sigma6/xsq**3

        V_LJ=eps4*sig6_x6*(sig6_x6-1.d0)
        dVLJ_dx=eps48/xsq*sig6_x6*(0.5-sig6_x6)

    end subroutine LJ_potential
    !##############################################

    !! writes first particles position and velocity and energy for the current time step
    subroutine write_output
        implicit none

        write(10,*) curr_time,pos(1,1),vel(1,1)
        write(11,*) curr_time,energy
        write(14,*) curr_time,temperature
        write(15,*) curr_time,dr2
    end subroutine write_output
    !##############################################

    subroutine write_dr2_output
        implicit none

        write(15,*) curr_time,dr2
    end subroutine write_dr2_output

    !! Writes vmd compatible file
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

    !! calculates energy, and average energy and energy squared for calculating standard deviation
    subroutine compute_energy
        implicit none
        integer i

        energy=pot
        do i=1,npart
          energy=energy+0.5*mass(i)*sum(vel(i,:)*vel(i,:))
        enddo

        E_avg=E_avg+energy
        Esq_avg=Esq_avg+energy**2

    end subroutine compute_energy
    !##############################################

    subroutine get_inside_the_box
        implicit none
        integer i,j
        do i=1,npart !! Loop over all particles
            do  j=1,3 !! Loop along x,y,z direction
                if (pos(i,j) > box_len) then
                     pos(i,j) = pos(i,j)-box_len
                else if (pos(i,j) < 0)then
                    pos(i,j) = pos(i,j)+box_len
                end if
            enddo
        enddo
    end subroutine get_inside_the_box

    subroutine gaussian_random_number(rnd)
    !! Using Box-Mueller transformation to generate a random number with gaussian distribution with mean=0 and sigma=1
    !! rnd*sig+x0 will give a gaussian random number with mean=x0 and standard deviation=sig,
        implicit none
        double precision,intent(out)::rnd
        double precision U1,U2
        call random_number(U1) !!! Calculates a random number over a uniform interval [0,1]
        call random_number(U2)
        rnd=sqrt(-2*log(U1))*cos(2*pi*U2)
    end subroutine gaussian_random_number


    subroutine thermal_velocities(temperature)
        implicit none
        double precision,intent(in) :: temperature
        double precision rnd
        integer i,j
        do i=1,npart
            do j=1,3
                call gaussian_random_number(rnd) !! Random number with gaussian distribution with mean=0, sigma=1
                vel(i,j)=rnd*sqrt(kb*temperature/mass(i)) !! Random number with gaussian distribution with mean=0, sigma=sqrt(kb*temperature/mass(i))
            enddo
        enddo
    end subroutine thermal_velocities


    subroutine cal_tempreture
        implicit none
        integer i,j
        do i=1,npart
            do j = 1,3
                temperature = temperature + (vel(i,j))**2
            end do
        end do
        temperature = temperature*mass(1)/(3*npart)
    end subroutine cal_tempreture



    subroutine cal_dr2
        implicit none
            dr2 = sum((pos - pos_t0)**2)
        dr2 = dr2/npart
    end subroutine cal_dr2

end module mod_verlet
