PROGRAM darius

	USE SCALAR, only : store_wisdom, step
	USE VAR_MPI, only : num_procs, my_id, ierr
	USE FFTW
	
	
	IMPLICIT NONE
	
    INTEGER :: num_arguments, num_step
	INTEGER :: count(1:5), count_rate, count_max
  	CHARACTER(len = 32) :: arg
	
     
	num_arguments = command_argument_count()

	if(num_arguments .eq. 1) then

		CALL getarg(1, arg)
		read(arg, *) store_wisdom 

	else
	   write(0,'(a,a,i3,a)') 'Illegal number of arguments.'
	   stop 1
	end if
	
! 	step = 0
	num_step = 1
	
! 	CALL MPI_INIT(ierr)
! 	CALL fftw_mpi_init
! 	CALL MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)
! 	CALL MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
	
	CALL system_clock(count(1))

	if(my_id .eq. 0) then	
		CALL PART_INIT
	else
		CALL MPI_PART_INIT
	end if
	
	CALL system_clock(count(2))
	
! 	do step = 0, num_step-1
! 		write(*, '(a, I5, a)') "Step number ", step, " "
		CALL PROJECT_DENSITY
		
		CALL system_clock(count(3))
		
		CALL FIELDS
		
		CALL system_clock(count(4))
		
		CALL UPDATE_PARTICLES
		
		CALL system_clock(count(5), count_rate, count_max)
! 	end do
! 	CALL PROJECT_DENSITY

! 	CALL MPI_FINALIZE(ierr)

	if(my_id.eq.0) then
		write(0,'(a, f10.4, a)') 'Part. Init. = ', dble((count(2)-count(1)))/count_rate, ' s'
		write(0,'(a, f10.4, a)') 'Proj. Dens. = ', dble((count(3)-count(2)))/count_rate, ' s'
		write(0,'(a, f10.4, a)') 'Poten. FFT  = ', dble((count(4)-count(3)))/count_rate, ' s'
		write(0,'(a, f10.4, a)') 'Part. Updt. = ', dble((count(5)-count(4)))/count_rate, ' s'
	
		write(0,'(a, f10.4, a)') '  Total execution time = ', dble((count(5)-count(1)))/count_rate, ' s'
	end if



END PROGRAM
