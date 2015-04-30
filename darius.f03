PROGRAM darius

	USE SCALAR, only : store_wisdom, step
	
	
	IMPLICIT NONE
	
    INTEGER :: num_arguments
  	CHARACTER(len = 32) :: arg
     
	num_arguments = command_argument_count()

	if(num_arguments .eq. 1) then

		CALL getarg(1, arg)
		read(arg, *) store_wisdom 

	else
	   write(0,'(a,a,i3,a)') 'Illegal number of arguments.'
	   stop 1
	end if
	
	step = 0
		
	CALL PART_INIT
	CALL PROJECT_DENSITY
	CALL FIELDS
	CALL UPDATE_PARTICLES
	CALL PROJECT_DENSITY

! 	CALL INIT_FIELDS 	! Fields are initialized
! 	CALL FIELDS 		! Fields are updated
! 	CALL WRITE_BOV
! 	CALL WRITE_FIELDS	! Fields are written to a data file for visualization

END PROGRAM
