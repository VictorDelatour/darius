PROGRAM darius

	USE SCALAR, only : store_wisdom, step
	USE VAR_MPI, only : num_procs, my_id, ierr
	USE FFTW
	
	
	IMPLICIT NONE
	
    INTEGER :: num_arguments, num_step
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
	
	CALL MPI_INIT(ierr)
	CALL fftw_mpi_init
	CALL MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)
	CALL MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
	
		
	if(my_id .eq. 0) then	
		CALL PART_INIT
	else
		CALL MPI_PART_INIT
	end if
	
! 	do step = 0, num_step-1
! 		write(*, '(a, I5, a)') "Step number ", step, " "
! 		CALL PROJECT_DENSITY
! 		CALL FIELDS
! 		CALL UPDATE_PARTICLES
! 	end do
! 	CALL PROJECT_DENSITY

	CALL MPI_FINALIZE(ierr)


END PROGRAM
