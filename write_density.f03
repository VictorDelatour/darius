SUBROUTINE WRITE_DENSITY(nx, ny, nz, density, step)
	
	IMPLICIT NONE
	
	INTEGER :: ierror, bov_id, density_id, step, origin
	INTEGER, INTENT(in) :: nx, ny, nz
	REAL*8, DIMENSION(nx, ny, nz), INTENT(in) :: density
	CHARACTER(LEN=128) :: bov_filename, dat_filename
	CHARACTER(LEN=128) :: step_string, nx_string, ny_string, nz_string, origin_string
	
	origin = 1
	
	write(step_string, "(I0)") step
	write(nx_string, "(I0)") nx
	write(ny_string, "(I0)") ny
	write(nz_string, "(I0)") nz
	write(origin_string, "(I0)") origin
	
	bov_filename = "./data/density_darius_" // TRIM(step_string) // ".bov"
	dat_filename = "./data/density_darius_" // TRIM(step_string) // ".dat"
	
	bov_id = 7
	density_id = 8
	
	open(unit = bov_id, file = TRIM(bov_filename), status = "replace", action = "readwrite")
	
	write(bov_id,*) "TIME: " //TRIM(step_string)
	write(bov_id,*) "DATA_FILE: density_darius_"//TRIM(step_string)//".dat"
	write(bov_id,*) "DATA_SIZE: "//TRIM(nx_string)//" "//TRIM(ny_string)//" "//TRIM(nz_string)
	write(bov_id,*) "DATA_FORMAT: DOUBLE"
	write(bov_id,*) "VARIABLE: density"
	write(bov_id,*) "DATA_ENDIAN: LITTLE"
	write(bov_id,*) "CENTERING: zonal"
	write(bov_id,*) "BRICK_ORIGIN: "//TRIM(origin_string)//" "//TRIM(origin_string)//" "//TRIM(origin_string)
	write(bov_id,*) "BRICK_SIZE: "//TRIM(nx_string)//" "//TRIM(ny_string)//" "//TRIM(nz_string)
	
	close(unit = bov_id)
	
	open(unit = density_id, file = TRIM(dat_filename), status = "replace", action = "readwrite", iostat = ierror, &
	    &  form = "unformatted", access = "direct", recl = nx*ny*nz*8)
		
	write(unit = density_id, rec = 1) density
	
	close(unit = density_id)
	
END SUBROUTINE