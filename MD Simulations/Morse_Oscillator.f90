program md
 implicit none
 !! variable declarations
    real*8 x,v,a,pot,mass,eps,pi,energy,sigma
    integer dt,total_time,no_of_steps,current_time,i
    x = 7
    v = 0
    eps = 0.000871
    mass = 145123
    sigma = 6.69
    dt = 75
    total_time = 1000000
    no_of_steps = total_time/dt
    pi = 3.14
    energy = 0
    call compute_pot(x,a,pot)
    
    
    open(10,file= "energy.out" )
    do i = 1,no_of_steps
        write(10,*) current_time,energy
        call evolve_by_1_step(x,v,a,energy)
        current_time = current_time + dt
        end do
        close(10)
 

    contains


    subroutine Evolve_by_1_step(x,v,a,energy)
        implicit none
        real*8,intent(inout) :: x,v,a,energy
    
        x = x + v*dt + 0.5*a*dt*dt
        v = v + 0.5*a*dt
        call compute_pot(x,a,pot)
        v = v + 0.5*a*dt
        energy = pot + 0.5*mass*v*v
    
        end subroutine Evolve_by_1_step
    
    subroutine calculate_pot(x,a,pot)
        implicit none
        real*8  sigma,eps
        real*8, intent(in)  :: x
        real*8, intent(out) :: a,pot
    
        a = 4*eps*(12*((sigma/12)/x**(13)) - 6*((sigma**6)/(x**7)))/mass
        pot = 4*eps*((sigma/x)**12 - (sigma/x)**6)
        end subroutine calculate_pot
    
end program md
    
