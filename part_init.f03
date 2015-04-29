SUBROUTINE PART_INIT
	
	USE DIMENSION, only: nx, ny, nz, nparticles
	USE VECTOR, only: x, y, z, vx, vy, vz, mass
	USE MATRIX, only: rho3d, phi3d
	
	IMPLICIT NONE
	
	INTEGER :: unit_part, unit_info, pos, min_level
	CHARACTER(LEN = 128) :: filename, buffer
	INTEGER :: wtf
	
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

	x = x * (nx - 1) + 1.0
	y = y * (ny - 1) + 1.0
	z = z * (nz - 1) + 1.0

	close(unit_part)
	
	do pos = 1, 10
		write(*,'(a F10.5 F10.5 F10.5 e15.6)') 'Position:', x(pos), y(pos), z(pos), mass(pos)
	end do
		
	
END SUBROUTINE PART_INIT