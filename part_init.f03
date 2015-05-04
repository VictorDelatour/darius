SUBROUTINE PART_INIT
	
	USE DIMENSION, only: nx, ny, nz, nparticles
	USE VECTOR, only: x, y, z, vx, vy, vz, mass
	USE MATRIX, only: rho3d, phi3d
	USE VAR_MPI, only: num_procs, my_id
	
	IMPLICIT NONE
	
	INTEGER :: unit_part, unit_info, pos, min_level
	CHARACTER(LEN = 128) :: filename, buffer
	INTEGER, DIMENSION(num_procs) :: part_per_slab
	INTEGER, DIMENSION(:), ALLOCATABLE :: slab_position
	
	
	
	unit_info = 1
	unit_part = 10
	
	filename = './output_00003/info_00003.txt'


	open(unit = unit_info, file = filename, status = 'old', form = 'formatted')

	read(unit_info, '(A13, I11)')
	read(unit_info, '(A13, I11)')
	read(unit_info, '(A13, I11)') buffer, min_level

	close(unit_info)

	nx = 2**min_level
	ny = nx
	nz = nx
	
	ALLOCATE( rho3d(nx, ny, nz) )
	ALLOCATE( phi3d(nx, ny, nz) )
	

	filename = './output_00003/part_00003.out00001'

	open(unit = unit_part, file = filename, status = 'old', form = 'unformatted')
	read(unit_part)
	read(unit_part)
	read(unit_part) nparticles
	
	write(*,'(a i10 a)') 'There are', nparticles, ' particles'

	ALLOCATE( x(nparticles) )
	ALLOCATE( y(nparticles) )
	ALLOCATE( z(nparticles) )
	ALLOCATE( vx(nparticles) )
	ALLOCATE( vy(nparticles) )
	ALLOCATE( vz(nparticles) )
	ALLOCATE( mass(nparticles) )
	
	ALLOCATE( slab_position(nparticles) )

	do pos = 1,5
    	read(unit_part) !skip
	end do
	
	read(unit_part) x
	read(unit_part) y
	read(unit_part) z
	
	read(unit_part) vx
	read(unit_part) vy
	read(unit_part) vz

	read(unit_part) mass 
	
	slab_position = nint( z * (num_procs - 1) + 1  ) 
	
	x = x * (nx - 1) + 1.0
	y = y * (ny - 1) + 1.0
	z = z * (nz - 1) + 1.0

	close(unit_part)
	
	part_per_slab = 0
	
	do pos = 1, nparticles
		part_per_slab(slab_position(pos)) = part_per_slab(slab_position(pos)) + 1
	end do
	
	
	if (my_id .eq. 0) then
		write(*,*) part_per_slab
		write(*,*) sum(part_per_slab)
	end if
	
! 	do pos = 1, 10
! 		write(*,'(a F10.5 F10.5 F10.5 e15.6)') 'Position:', x(pos), y(pos), z(pos), mass(pos)
! 	end do
		
	
END SUBROUTINE PART_INIT