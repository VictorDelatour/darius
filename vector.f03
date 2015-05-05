MODULE VECTOR
	
	IMPLICIT NONE
	
	REAL*8, DIMENSION(:), ALLOCATABLE :: x, y, z, vx, vy, vz, mass
	REAL*8, DIMENSION(:), ALLOCATABLE :: part_list
	REAL*8, DIMENSION(:), ALLOCATABLE :: local_x, local_y, local_z, local_vx, local_vy, local_vz,local_mass
	
END MODULE VECTOR