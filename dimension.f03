MODULE DIMENSION
	
	USE, INTRINSIC :: iso_c_binding
	
	IMPLICIT NONE

	INTEGER :: nparticles
	INTEGER(C_INTPTR_T) :: nx, ny, nz
	REAL*8 :: low, upp, len
	
	!integer :: nx, ny, nz !dimensions of the cube
	
END MODULE DIMENSION