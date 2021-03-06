SUBROUTINE PART_INIT
	
	USE DIMENSION, only: nx, ny, nz, nparticles
	USE VECTOR, only: x, y, z, vx, vy, vz, mass, local_x, local_y, local_z, local_vx, local_vy, local_vz,local_mass
	USE MATRIX, only: rho3d, phi3d
	USE VAR_MPI, only: num_procs, my_id
	USE FFTW
	
	IMPLICIT NONE
	
	INTEGER :: unit_part, unit_info, pos, slab_pos, part_iter_pos, min_level, proc
	INTEGER :: ierror, send_request
	CHARACTER(LEN = 128) :: filename, buffer
	INTEGER, DIMENSION(num_procs) :: part_per_slab, part_iterator
	INTEGER, DIMENSION(:), ALLOCATABLE :: slab_position
	INTEGER, DIMENSION(:, :), ALLOCATABLE :: part_in_slab
	
	
	
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
	
	if (num_procs > nz) then
		write(0,('a')) "Cannot have more procs than there are z-slabs"
		stop
	end if
	
	! Don't need this anymore
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
	
	! Get
	slab_position = nint( z * (num_procs - 1) + 1  ) 
	
	x = x * (nx - 1) + 1.0
	y = y * (ny - 1) + 1.0
	z = z * (nz - 1) + 1.0

	close(unit_part)
	
	part_per_slab = 0

	
	do pos = 1, nparticles
		part_per_slab(slab_position(pos)) = part_per_slab(slab_position(pos)) + 1
	end do
	
	do proc = 1, num_procs-1
		CALL MPI_ISEND(part_per_slab(proc+1), 1, MPI_INT, proc, 1, MPI_COMM_WORLD, proc, ierror)
	end do
	
	
! 	if (my_id .eq. 0) then
! 		write(*,*) part_per_slab
! 		write(*,*) sum(part_per_slab)
! 		write(*,*) maxval(part_per_slab)
! 	end if


! This is cheating, particles will move and get out of slabs
	ALLOCATE( part_in_slab(num_procs, maxval(part_per_slab)))
	
	part_iterator = 1
	part_in_slab = -1
	
	do pos = 1, nparticles
		
		slab_pos = slab_position(pos)
		part_iter_pos = part_iterator(slab_pos)
		
		part_in_slab(slab_pos, part_iter_pos) = pos
		
		part_iterator(slab_pos) = part_iter_pos + 1
		
	end do
	
	do proc = 1, num_procs
		write(*,*) proc-1, part_per_slab(proc)
! 		write(*,*) part_in_slab(proc, part_iterator(proc))
! 		write(*,*) part_in_slab(proc, part_iterator(proc)-1)
	end do
	
	do proc = 1, num_procs-1
		write(*,'(a I0)') "Sending data to proc", proc
! 		CALL MPI_RECV(gradz(1, 1, 1 + i*nz/num_procs), alloc_local, MPI_DOUBLE, i, 4, MPI_COMM_WORLD, status, ierror)
! 		Consider using non-blocking communication
! ! 		CALL MPI_SEND(part_in_slab(proc, 1:(part_iterator(proc+1)-1)), part_iterator(proc+1)-1, MPI_DOUBLE, proc, 2, MPI_COMM_WORLD, ierror)
! 		CALL MPI_SEND( x(part_in_slab(proc+1, 1:(part_iterator(proc+1)-1))) , part_iterator(proc+1)-1, MPI_DOUBLE, proc, 2, MPI_COMM_WORLD, ierror)
! 		CALL MPI_SEND( y(part_in_slab(proc+1, 1:(part_iterator(proc+1)-1))) , part_iterator(proc+1)-1, MPI_DOUBLE, proc, 3, MPI_COMM_WORLD, ierror)
! 		CALL MPI_SEND( z(part_in_slab(proc+1, 1:(part_iterator(proc+1)-1))) , part_iterator(proc+1)-1, MPI_DOUBLE, proc, 4, MPI_COMM_WORLD, ierror)
! 		CALL MPI_SEND( vx(part_in_slab(proc+1, 1:(part_iterator(proc+1)-1))) , part_iterator(proc+1)-1, MPI_DOUBLE, proc, 5, MPI_COMM_WORLD, ierror)
! 		CALL MPI_SEND( vy(part_in_slab(proc+1, 1:(part_iterator(proc+1)-1))) , part_iterator(proc+1)-1, MPI_DOUBLE, proc, 6, MPI_COMM_WORLD, ierror)
! 		CALL MPI_SEND( vz(part_in_slab(proc+1, 1:(part_iterator(proc+1)-1))) , part_iterator(proc+1)-1, MPI_DOUBLE, proc, 7, MPI_COMM_WORLD, ierror)
! 		CALL MPI_SEND( mass(part_in_slab(proc+1, 1:(part_iterator(proc+1)-1))) , part_iterator(proc+1)-1, MPI_DOUBLE, proc, 8, MPI_COMM_WORLD, ierror)
	end do
	
! 	ALLOCATE( part_list(part_per_slab(1)))
	ALLOCATE( local_x(part_per_slab(1)) )
	ALLOCATE( local_y(part_per_slab(1)) )
	ALLOCATE( local_z(part_per_slab(1)) )
	ALLOCATE( local_vx(part_per_slab(1)) )
	ALLOCATE( local_vy(part_per_slab(1)) )
	ALLOCATE( local_vz(part_per_slab(1)) )
	ALLOCATE( local_mass(part_per_slab(1)) )
	
	local_x = x(part_in_slab(1, 1:(part_iterator(1)-1)))
	local_y = y(part_in_slab(1, 1:(part_iterator(1)-1)))
	local_z = z(part_in_slab(1, 1:(part_iterator(1)-1)))
	local_vx = vx(part_in_slab(1, 1:(part_iterator(1)-1)))
	local_vy = vy(part_in_slab(1, 1:(part_iterator(1)-1)))
	local_vz = vz(part_in_slab(1, 1:(part_iterator(1)-1)))
	local_mass = mass(part_in_slab(1, 1:(part_iterator(1)-1)))
	
	
! 	part_list = part_in_slab(proc, 1:(part_iterator(1)-1))
	
	
!
! 	do pos = 1, num_procs
! 		write(*,*) part_in_slab(pos, 1:10)
! 	end do
	
! 	do pos = 1, 10
! 		write(*,'(a F10.5 F10.5 F10.5 e15.6)') 'Position:', x(pos), y(pos), z(pos), mass(pos)
! 	end do
		
	
END SUBROUTINE PART_INIT