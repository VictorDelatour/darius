MODULE VAR_MPI
	
	USE, INTRINSIC :: iso_c_binding
	
	IMPLICIT NONE
	
	INTEGER :: num_procs, my_id, ierr
	
	INTEGER(C_INTPTR_T) :: Lx, Ly, Lz
	INTEGER(C_INTPTR_T) :: ii, jj, kk, dim3
	INTEGER(C_INTPTR_T) :: local_nx, local_ny, local_nz, alloc_local 
	INTEGER(C_INTPTR_T) :: local_x_offset, local_y_offset, local_z_offset
	
	
	
END MODULE VAR_MPI