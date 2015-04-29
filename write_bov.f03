SUBROUTINE WRITE_BOV
	
	USE DIMENSION, only : nx, ny, nz, low, len
	
	IMPLICIT NONE
	
	INTEGER :: rho_id, phi_id
	
	rho_id = 7
	phi_id = 8
	
	open(unit = rho_id, file = "density.bov", status = "replace", action = "readwrite")
	
	write(rho_id,*) "TIME: 1"
	write(rho_id,*) "DATA_FILE: rho3d_MPI.bin"
	write(rho_id,*) "DATA_SIZE: ", nx, ny, nz
	write(rho_id,*) "DATA_FORMAT: DOUBLE"
	write(rho_id,*) "VARIABLE: density"
	write(rho_id,*) "DATA_ENDIAN: LITTLE"
	write(rho_id,*) "CENTERING: zonal"
	write(rho_id,*) "BRICK_ORIGIN: ", low, low, low
	write(rho_id,*) "BRICK_SIZE: ", len, len, len
	
	close(unit = rho_id)
	
	
	open(unit = phi_id, file = "potential.bov", status = "replace", action = "readwrite")
	
	write(phi_id,*) "TIME: 1"
	write(phi_id,*) "DATA_FILE: phi3d_MPI.bin"
	write(phi_id,*) "DATA_SIZE: ", nx, ny, nz
	write(phi_id,*) "DATA_FORMAT: DOUBLE"
	write(phi_id,*) "VARIABLE: potential"
	write(phi_id,*) "DATA_ENDIAN: LITTLE"
	write(phi_id,*) "CENTERING: zonal"
	write(phi_id,*) "BRICK_ORIGIN: ", low, low, low
	write(phi_id,*) "BRICK_SIZE: ", len, len, len
	
	close(unit = phi_id)	
	
	
END SUBROUTINE WRITE_BOV