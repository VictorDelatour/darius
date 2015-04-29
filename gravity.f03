PROGRAM gravity
	
	USE DIMENSION, only : nx, ny, nz  ! Contains dimension-related variables
	USE MATRIX, only : rho3d, phi3d	  ! Contains the variables rho3d and phi3d
	USE SCALAR, only : first
	
	IMPLICIT NONE
	
    INTEGER :: i
    INTEGER, EXTERNAL :: iargc
    CHARACTER(len = 32) :: arg
     

	if(iargc() .eq. 0) then
		nx = 64 ! Size of the cube 
		ny = nx
		nz = nx
	else if(iargc() .eq. 4) then
		CALL getarg(1, arg)
		read(arg, *) nx 
		
		CALL getarg(2, arg)
		read(arg, *) ny 
		
		CALL getarg(3, arg)
		read(arg, *) nz 
		
		CALL getarg(4, arg)
		read(arg, *) first ! First = is it the first time this simulation is run, i.e. should we store the wisdom? 1 = yes, 0 = no, wisdom already existent
		
		print*, nx, ny, nz, first

	else
	   write(0,'(a,a,i3,a)') 'Illegal number of arguments.'
	   stop 1
	end if
		

	
	ALLOCATE( rho3d(nx, ny, nz) )
	ALLOCATE( phi3d(nx, ny, nz) )
	

	CALL INIT_FIELDS 	! Fields are initialized
	CALL FIELDS 		! Fields are updated
	CALL WRITE_BOV
	!CALL WRITE_FIELDS	! Fields are written to a data file for visualization

END PROGRAM
