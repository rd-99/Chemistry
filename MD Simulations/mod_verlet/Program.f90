Program molecular_dynamics
    !! Modular molecular dynamics code in Fortran
    !! Output: a file called 'std_en_dt.out' containing 2 columns: time step (dt), standard deviation in energy (both in atomic units)
    !! Input read from md_fput.inp
    use mod_verlet
    implicit none

    call setup_parameters
    call compute_std_dev_en


end program molecular_dynamics
