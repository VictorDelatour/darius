SUBROUTINE WRITE_FIELDS

	USE DIMENSION, only : nx, ny, nz  ! Contains dimension-related variables
	USE MATRIX, only : rho3d, phi3d	  ! Contains the variables rho3d and phi3d
	
	IMPLICIT NONE
	
	INTEGER :: ierror
	
	open(unit = 1, file = "rho3d.bin", status = "replace", action = "readwrite", iostat = ierror, &
	    &  form = "unformatted", access = "direct", recl = nx*ny*nz*8)
	
	if(ierror /= 0) then
		print*, "Failed to open 'rho3d.bin'"
		stop
	end if
	
	open(unit = 2, file = "phi3d.bin", status = "replace", action = "readwrite", iostat = ierror, &
	    &  form = "unformatted", access = "direct", recl = nx*ny*nz*8)

	
	if(ierror /= 0) then
		print*, "Failed to open 'phi3d.bin'"
		stop
	end if
	
	write(unit = 1, rec = 1) rho3d
	write(unit = 2, rec = 1) phi3d
	
	
	close(unit = 1)
	close(unit = 2)
	
END SUBROUTINE WRITE_FIELDS