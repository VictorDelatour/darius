SUBROUTINE FIELDS

	USE DIMENSION, only : nx, ny, nz, low, upp, len  ! Contains dimension-related variables
	USE MATRIX, only : rho3d, phi3d, gradx, grady, gradz	  ! Contains the variables rho3d and phi3d
	USE SCALAR, only : i, j, k, gconst, at, ga, pi, store_wisdom, step
	USE FFTW
	USE VAR_MPI
	
	IMPLICIT NONE
	
	
	! Define values
	REAL*8 :: k2, invkk, L, beta, sigma, shift, sum, p
	REAL*8, DIMENSION(3) :: pos
	REAL*8, DIMENSION(nx) :: kq 
	REAL*8, DIMENSION(:,:,:), allocatable :: subarray
	REAL*8, DIMENSION(:,:,:), allocatable :: tempstorage, local_data
	REAL*8 :: temp
	
	INTEGER :: fh_init, fh_rho, fh_phi, ierror, import
	INTEGER :: status(MPI_STATUS_SIZE)
	INTEGER :: count(1:6), count_rate, count_max
	INTEGER, DIMENSION(3) :: sizes, subsizes, start
	INTEGER :: subarray_type
	INTEGER(KIND = MPI_OFFSET_KIND) :: offset
	INTEGER(C_INT) :: ret
	INTEGER, DIMENSION(1) :: send, recv

	TYPE(C_PTR) :: planFOR, planBACK, cdata, cgrad_x, cgrad_y, cgrad_z
	COMPLEX(C_DOUBLE_COMPLEX), pointer :: data(:,:,:), grad_x(:,:,:), grad_y(:,:,:), grad_z(:,:,:)
	COMPLEX(C_DOUBLE_COMPLEX), PARAMETER :: icomplex = (0.0,1.0)

	CALL system_clock(count(1))
	
	!MPI & FFTW initialization
	CALL MPI_INIT(ierr)
	CALL fftw_mpi_init
	CALL MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)
	CALL MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
	
	
		
	! Open file to store the density
	CALL MPI_FILE_OPEN(MPI_COMM_WORLD, "./data/rho3d_MPI.bin", &
	 IOR(MPI_MODE_CREATE, MPI_MODE_RDWR), MPI_INFO_NULL, fh_rho, ierr)
	! Open file to store the potential
	CALL MPI_FILE_OPEN(MPI_COMM_WORLD, "./data/phi3d_MPI.bin", &
	 IOR(MPI_MODE_CREATE, MPI_MODE_RDWR), MPI_INFO_NULL, fh_phi, ierr)

	
	! Allocate memory to z-slab for MPI
	alloc_local = fftw_mpi_local_size_3d(nx, ny, nz, MPI_COMM_WORLD, local_nz, local_z_offset)
	
	if(step .eq. 0) then
		ALLOCATE( gradx(nx, ny, nz) )
		ALLOCATE( grady(nx, ny, nz) )
		ALLOCATE( gradz(nx, ny, nz) )
		ALLOCATE( subarray(nx, ny, local_nz) )
	end if
	
	
	! Allocate memory
	cdata = fftw_alloc_complex(alloc_local)
	CALL c_f_pointer(cdata, data, [nx, ny, local_nz])
	
	
	cgrad_x = fftw_alloc_complex(alloc_local)
	cgrad_y = fftw_alloc_complex(alloc_local)
	cgrad_z = fftw_alloc_complex(alloc_local)

	CALL c_f_pointer(cgrad_x, grad_x, [nx, ny, local_nz])
	CALL c_f_pointer(cgrad_y, grad_y, [nx, ny, local_nz])
	CALL c_f_pointer(cgrad_z, grad_z, [nx, ny, local_nz])
	

	 
	! Prepare creation of subarray for parallel I/O
	sizes = (/ nx, ny, nz /) 
	subsizes = (/ nx, ny, local_nz /)
	start = (/ 0, 0, INT(local_z_offset) /)
	
	! Create subarray for parallel I/O
	CALL MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, start, MPI_ORDER_FORTRAN, MPI_DOUBLE, subarray_type, ierr)
	CALL MPI_TYPE_COMMIT(subarray_type, ierr)
	
	! Set view
	offset = 0 ! The start variable (^) takes care of the offset of the z-slab
	!CALL MPI_FILE_SET_VIEW(fh_rho, offset, MPI_DOUBLE, subarray_type, "native", MPI_INFO_NULL, ierr)
	CALL MPI_FILE_SET_VIEW(fh_rho, offset, MPI_DOUBLE, subarray_type, "native", MPI_INFO_NULL, ierr)
	CALL MPI_FILE_SET_VIEW(fh_init, offset, MPI_DOUBLE, subarray_type, "native", MPI_INFO_NULL, ierr)
	CALL MPI_FILE_SET_VIEW(fh_phi, offset, MPI_DOUBLE, subarray_type, "native", MPI_INFO_NULL, ierr)
	
	!Value of k^2
	! Unused
	do i = 1, nx
		kq(i) = ( 2*pi* (i-.5*(nx+1))/(.5*nx) ) ** 2 
	end do
	
	! Defined FFTW plans, using WISDOM
	! Should the plan be for local_nz, or for nz?
	
	
	
	if(store_wisdom /= 1) then
		if(my_id .eq. 0) then
			ret = FFTW_IMPORT_WISDOM_FROM_FILENAME(C_CHAR_'wisdom_FWD.dat' // C_NULL_CHAR)
		end if
		CALL FFTW_MPI_BROADCAST_WISDOM(MPI_COMM_WORLD)
	end if
	
	planFOR = fftw_mpi_plan_dft_3d(nx, ny, nz, data, data, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_PATIENT)
	
	if(store_wisdom .eq. 1) then
		CALL FFTW_MPI_GATHER_WISDOM(MPI_COMM_WORLD)
		if(my_id .eq. 0) then
			ret = FFTW_EXPORT_WISDOM_TO_FILENAME(C_CHAR_'wisdom_FWD.dat' // C_NULL_CHAR)
		end if
	end if
	

	if(store_wisdom /= 1) then
		if(my_id .eq. 0) then
			ret = FFTW_IMPORT_WISDOM_FROM_FILENAME(C_CHAR_'wisdom_BWD.dat' // C_NULL_CHAR)
		end if
		CALL FFTW_MPI_BROADCAST_WISDOM(MPI_COMM_WORLD)
	end if

	planBACK = fftw_mpi_plan_dft_3d(nx, ny, nz, data, data, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_PATIENT)
	
	if(store_wisdom .eq. 1) then
		CALL FFTW_MPI_GATHER_WISDOM(MPI_COMM_WORLD)
		if(my_id .eq. 0) then
			ret = FFTW_EXPORT_WISDOM_TO_FILENAME(C_CHAR_'wisdom_BWD.dat' // C_NULL_CHAR)
		end if
	end if
	

	
	low = 1
	upp = nx
	len = upp - low
	
	! Bad practice, only 1/nnum_procs will do this loop once, all of them will do it twice => loose time
		
	if( my_id .eq. 0) then
		
		do i = 1, num_procs-1
			! Assume wrongly that local size alloc_local is the same for each, which is wrong.
			CALL MPI_SEND(rho3d(1, 1, 1 + i*nz/num_procs), alloc_local, MPI_DOUBLE, i, 2, MPI_COMM_WORLD, ierror)
		end do
		 
		data(:, :, :) = rho3d(:, :, 1:(nz/num_procs))
		
	else
		
		if(step .eq. 0) then
			ALLOCATE( local_data(nx, ny, local_nz))
		end if
		
		CALL MPI_RECV(local_data(1,1,1), alloc_local, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, status, ierror)
		data = local_data ! Probably useless?
		
	end if
	

	
	
	
	if(my_id.eq.0) then
		print*, "DATA INITIALIZED"
	end if
	
! Put it inside the the if/else logical struct 
	
	subarray = real(data) ! Only useful for printing the data to the file
	
	if(my_id.eq.0) then
		print*, "WRITE RHO"
	end if
	
	
! 	We alreay do it in the density projection, but always good to have it?
	
	CALL MPI_FILE_WRITE_ALL(fh_rho, subarray, alloc_local, MPI_REAL8, MPI_STATUS_IGNORE)
	
	if(my_id.eq.0) then
		print*, "FILE CLOSED"
	end if
	
	CALL MPI_FILE_CLOSE(fh_rho, ierr)
	
	if(my_id.eq.0) then
		print*, "RHO WRITTEN"
	end if

	CALL system_clock(count(2)) ! Time for initialization
	
	
	if(my_id.eq.0) then
		print*, "FORWARD FFT"
	end if
	
	
	! Execute plan => Forward FFT
	CALL fftw_mpi_execute_dft(planFOR, data, data)
	
	if(my_id.eq.0) then
		print*, "COMPLETE..."
	end if
	
	
	CALL system_clock(count(3)) ! Time for the Forward FFT
	

	
	p = (2*pi*nx/(upp-low)) ** 2

	! calculate local gravitational field
	do k = 1, local_nz
		
		! Avoid doing two branches.
		if(2 * (k+local_z_offset) <= nz ) then
			pos(3) = k - 1
		else
			pos(3) = k - nz - 1
		end if
		
		
		do j = 1, ny
			
			if(2*j <= ny) then
				pos(2) = j - 1
			else 
				pos(2) = j - ny - 1
			end if
			
			
			do i = 1, nx
				
				if(2*i <= nx) then
					pos(1) = i - 1
				else
					pos(1) = i - nx - 1
				end if
				
				

 				k2 = p * ( pos(1)**2 + pos(2)**2 + pos(3)**2 ) / (nx*nx)
				if(k2 /= 0) then
					data(i, j, k) = data(i, j, k) * 4  * pi / k2 !* gconst
				else
					data(i, j, k) = 0.
				end if
				
				! Do the gradient directly here
				grad_x(i, j, k) = pos(1) * icomplex * data(i, j, k)
				grad_y(i, j, k) = pos(2) * icomplex * data(i, j, k)
				grad_z(i, j, k) = pos(3) * icomplex * data(i, j, k)

			end do
		end do
	end do
	
	if(my_id.eq.0) then
		print*, "GRAVITATIONAL FIELD COMPUTED"
	end if
	

	CALL system_clock(count(4)) ! Time for the update
	
	! Execute plan => Backward FFT
	CALL fftw_mpi_execute_dft(planBACK, data, data)	
	

	
	if(my_id.eq.0) then
		print*, "NORMALIZE PHI"
	end if
	
	! Normalize data
	subarray = real(data) / (nx * ny * nz )
	
	if(my_id.eq.0) then
		print*, "RETURN GRADIENT"
	end if
	
	! In place inverse FFT of the gradient
	! Don't need to write the gradient to file,
	! But you need to send it to the master process
	! Figure out how later
	CALL fftw_mpi_execute_dft(planBACK, grad_x, grad_x)
	CALL fftw_mpi_execute_dft(planBACK, grad_y, grad_y)
	CALL fftw_mpi_execute_dft(planBACK, grad_z, grad_z)
	
	CALL system_clock(count(5)) ! Time for Backward FFT
	
	
	if( my_id .eq. 0) then
		
		
		do i = 1, num_procs-1
			CALL MPI_RECV(gradx(1, 1, 1 + i*nz/num_procs), alloc_local, MPI_DOUBLE, i, 2, MPI_COMM_WORLD, status, ierror)
			CALL MPI_RECV(grady(1, 1, 1 + i*nz/num_procs), alloc_local, MPI_DOUBLE, i, 3, MPI_COMM_WORLD, status, ierror)
			CALL MPI_RECV(gradz(1, 1, 1 + i*nz/num_procs), alloc_local, MPI_DOUBLE, i, 4, MPI_COMM_WORLD, status, ierror)
		end do
		
	else
		CALL MPI_SEND(grad_x, alloc_local, MPI_DOUBLE, i, 2, MPI_COMM_WORLD, ierror)
		CALL MPI_SEND(grad_y, alloc_local, MPI_DOUBLE, i, 2, MPI_COMM_WORLD, ierror)
		CALL MPI_SEND(grad_z, alloc_local, MPI_DOUBLE, i, 2, MPI_COMM_WORLD, ierror)
	end if
	
	! So now you have the gradient stored in gradx, grady, gradz
	

	
	if(my_id.eq.0) then
		print*, "WRITE PHI"
	end if
	
	! Parallel write
	CALL MPI_FILE_WRITE_ALL(fh_phi, subarray, alloc_local, MPI_REAL8, MPI_STATUS_IGNORE)	
	CALL MPI_FILE_CLOSE(fh_phi, ierr)
	
	! Add part where every slave sends back its data to the master process
	
	CALL system_clock(count(6), count_rate, count_max)


	if(my_id.eq.0) then
		print*, "WRITTEN PHI"
	end if
	
	CALL system_clock(count(7), count_rate, count_max)
	
	
	if(my_id.eq.0) then
		write(0,'(a, f10.4, a)') '  Initialization = ', dble((count(2)-count(1)))/count_rate, ' s'
		write(0,'(a, f10.4, a)') '  Forward FFT = ', dble((count(3)-count(2)))/count_rate, ' s'
		write(0,'(a, f10.4, a)') '  Update = ', dble((count(4)-count(3)))/count_rate, ' s'
		write(0,'(a, f10.4, a)') '  Backward FFT = ', dble((count(5)-count(4)))/count_rate, ' s'
		write(0,'(a, f10.4, a)') '  Write = ', dble((count(6)-count(5)))/count_rate, ' s'
		
		write(0,'(a, f10.4, a)') '  Total execution time = ', dble((count(6)-count(1)))/count_rate, ' s'
		write(0,'(a, f10.4, a)') '  Total computation time = ', dble((count(5)-count(2)))/count_rate, ' s'
		
	end if
	
	
	CALL fftw_destroy_plan(planBACK)
 	CALL fftw_destroy_plan(planFOR)
	
	CALL MPI_FINALIZE(ierr)
	
END SUBROUTINE FIELDS
