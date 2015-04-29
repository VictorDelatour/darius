SUBROUTINE FIELDS

	USE DIMENSION, only : nx, ny, nz, low, upp, len  ! Contains dimension-related variables
	USE MATRIX, only : rho3d, phi3d	  ! Contains the variables rho3d and phi3d
	USE SCALAR, only : i, j, k, gconst, at, ga, pi, store_wisdom
	USE FFTW
	USE VAR_MPI
	
	IMPLICIT NONE
	
	
	! Define values
	REAL*8 :: k2, invkk, L, beta, sigma, shift, sum, p
	REAL*8, DIMENSION(3) :: center, pos, pos_2, pos_3 ! Center of the gaussian
	REAL*8, DIMENSION(nx) :: kq 
	REAL*8, DIMENSION(:,:,:), allocatable :: subarray
	REAL*4, DIMENSION(:,:,:), allocatable :: tempstorage, local_data
	REAL*8 :: temp
	
	INTEGER :: ierr, num_procs, my_id, fh_init, fh_rho, fh_phi, ierror, import
	INTEGER :: status(MPI_STATUS_SIZE)
	INTEGER :: count(1:6), count_rate, count_max
	INTEGER, DIMENSION(3) :: sizes, subsizes, start
	INTEGER :: subarray_type
	INTEGER(KIND = MPI_OFFSET_KIND) :: offset
	INTEGER(C_INT) :: ret
	INTEGER, DIMENSION(1) :: send, recv

	TYPE(C_PTR) :: planFOR, planBACK, cdata
	COMPLEX(C_DOUBLE_COMPLEX), pointer :: data(:,:,:)

	CALL system_clock(count(1))
	
	!MPI & FFTW initialization
	CALL MPI_INIT(ierr)
	CALL fftw_mpi_init
	CALL MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)
	CALL MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
	
	! Open file to store the density
	!CALL MPI_FILE_OPEN(MPI_COMM_WORLD, "cube128.bin", &
	 !MPI_MODE_RDONLY, MPI_INFO_NULL, fh_init, ierr)
	
	! Open file to store the density
	CALL MPI_FILE_OPEN(MPI_COMM_WORLD, "TO_PC/rho3d_MPI.bin", &
	 IOR(MPI_MODE_CREATE, MPI_MODE_RDWR), MPI_INFO_NULL, fh_rho, ierr)
	! Open file to store the potential
	CALL MPI_FILE_OPEN(MPI_COMM_WORLD, "TO_PC/phi3d_MPI.bin", &
	 IOR(MPI_MODE_CREATE, MPI_MODE_RDWR), MPI_INFO_NULL, fh_phi, ierr)

	
	! Allocate memory to z-slab for MPI
	alloc_local = fftw_mpi_local_size_3d(nx, ny, nz, MPI_COMM_WORLD, local_nz, local_z_offset)
	ALLOCATE( subarray(nx, ny, local_nz) )
	
	! Allocate memory
	cdata = fftw_alloc_complex(alloc_local)
	CALL c_f_pointer(cdata, data, [nx, ny, local_nz])
	print*, alloc_local
	
	 
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
	CALL MPI_FILE_SET_VIEW(fh_rho, offset, MPI_FLOAT, subarray_type, "native", MPI_INFO_NULL, ierr)
	CALL MPI_FILE_SET_VIEW(fh_init, offset, MPI_DOUBLE, subarray_type, "native", MPI_INFO_NULL, ierr)
	CALL MPI_FILE_SET_VIEW(fh_phi, offset, MPI_DOUBLE, subarray_type, "native", MPI_INFO_NULL, ierr)
	
	!Value of k^2
	do i = 1, nx
		kq(i) = ( 2*pi* (i-.5*(nx+1))/(.5*nx) ) ** 2 
	end do
	
	! Defined FFTW plans, using WISDOM
	
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
	
	
	low = -5 !previously, -1
	upp = 5 !previously, 1
	len = upp - low
	

	import = 1
	
	if(import .eq. 1) then
		
		if( my_id .eq. 0) then
			! Read only on unit 1
			
			ALLOCATE( tempstorage(nx, ny, nz) ) ! Allocate memory only on process 0
			
			open(unit = 1, file = "cube128.bin", action = "read", iostat = ierror, &
				 & form = "unformatted", access = "direct", recl = nx*ny*nz*4 )

		

			read(unit = 1, rec = 1) tempstorage
			print*, nx, ny, nz
			
! 			do i = 1, 10
! 				print*, my_id, tempstorage(i,1,1)
! 			end do
			

			do i = 1, num_procs-1
				send = my_id + i
! 				
! 				CALL MPI_SEND(send, 1, MPI_INTEGER, i, 1, MPI_COMM_WORLD, ierror)
				CALL MPI_SEND(tempstorage(1, 1, 1 + i*nz/num_procs), alloc_local, &
					 			& MPI_FLOAT, i, 2, MPI_COMM_WORLD, ierror)

			end do
			 
			data(:, :, :) = tempstorage(:, :, 1:(nz/num_procs))
	
! 			do i = 1, 3
! 				print*, my_id, data(i,1,1)
! 			end do
			
		else
! 			CALL MPI_RECV(recv, 1,MPI_INTEGER, 0, 1, MPI_COMM_WORLD, status, ierror)
			
			ALLOCATE( local_data(nx, ny, local_nz))
			
			CALL MPI_RECV(local_data(1,1,1), alloc_local, MPI_FLOAT, 0, 2, MPI_COMM_WORLD, status, ierror)
!
! 			do i = 1, 3
! 				print*, my_id,  local_data(i,1,1)
! 			end do
			
			data = local_data
			
! 			do i = 1, 3
! 				print*, my_id, data(i,1,1)
! 			end do
			
			
			
		end if
		
	else
		
		! Initialize data with three exponentials
		! Center of the gaussian at the middle of the box	
	
		beta = 1
		sigma = .2
		center = .5*( 1 + (/nx, ny, nz/) )
		shift = nx/4

		do i = 1, nx
	! 		pos(1) = len*(i - center(1))/nx
			pos(1) = len*(i - center(1) - shift)/nx
			pos_2(1) = len*(i - center(1) + shift)/nx
			pos_3(1) = len*(i - center(1) + shift)/nx
		
			do j = 1, ny
	! 			pos(2) = len*(j - center(2))/ny
				pos(2) = len*(j - center(2) - shift)/ny
				pos_2(2) = len*(j - center(2) + shift)/ny
				pos_3(2) = len*(j - center(2) - shift)/ny
			
				do k = 1, local_nz
	! 				pos(3) = len*(k+local_z_offset - center(3))/nz				
					pos(3) = len*(k + local_z_offset - center(3) - shift)/nz
					pos_2(3) = len*(k + local_z_offset - center(3) + shift)/nz
					pos_3(3) = len*(k + local_z_offset - center(3) - shift)/nz
	! 				sum = pos(1)**2 + pos(2)**2 + pos(3)**2
	! 				data(i, j, k) = -(4 * sum - 6) / (4*pi) * exp(-sum)
					data(i, j, k) = beta * exp( -( pos(1)**2 + pos(2)**2 + pos(3)**2 ) / (sigma**2) )
					data(i, j, k) = data(i, j, k) + beta * exp( -( pos_2(1)**2 + pos_2(2)**2 + pos_2(3)**2 ) / (sigma**2) )
					data(i, j, k) = data(i, j, k) + beta * exp( -( pos_3(1)**2 + pos_3(2)**2 + pos_3(3)**2 ) / (sigma**2) )
				
					!data(i, j, k) = ( pos(1) + pos(2) + pos(3) )/3
				
				end do
			end do
		end do
		
	end if 

	print*, my_id, "DATA INITIALIZED"
	
! 	if(my_id.eq.0) then
! 		print*, "DATA INITIALIZED"
! 	end if
	
	subarray = real(data)
	
	print*, my_id, "WRITE RHO"
	
! 	if(my_id.eq.0) then
! 		print*, "WRITE RHO"
! 	end if
	
	CALL MPI_FILE_WRITE_ALL(fh_rho, subarray, alloc_local, MPI_REAL8, MPI_STATUS_IGNORE)	
	CALL MPI_FILE_CLOSE(fh_rho, ierr)
	
	
	
	print*, my_id, "RHO WRITTEN"
	
! 	if(my_id.eq.0) then
! 		print*, "RHO WRITTEN"
! 	end if

	CALL system_clock(count(2)) ! Time for initialization
	
	print*, my_id, "FORWARD FFT"
	
! 	if(my_id.eq.0) then
! 		print*, "FORWARD FFT"
! 	end if
	
	
	! Execute plan => Forward FFT
	CALL fftw_mpi_execute_dft(planFOR, data, data)
	
	print*, my_id, "COMPLETE..."
	
! 	if(my_id.eq.0) then
! 		print*, "COMPLETE..."
! 	end if
	
	
	CALL system_clock(count(3)) ! Time for the Forward FFT
	

	
	p = (2*pi*nx/(upp-low)) ** 2

	! calculate local gravitational field
	do k = 1, local_nz
		
		if( 2 * (k+local_z_offset) > nz ) then
			pos(3) = k - nz - 1
		else
			pos(3) = k - 1
		end if
		
		do j = 1, ny
			
			if( 2*j > ny ) then
				pos(2) = j - ny - 1
			else
				pos(2) = j - 1
			end if
			
			do i = 1, nx
				
				if( 2*i > nx ) then
					pos(1) = i - nx - 1
				else
					pos(1) = i - 1
				end if

 				k2 = p * ( pos(1)**2 + pos(2)**2 + pos(3)**2 ) / (nx*nx)
				if(k2 /= 0) then
					data(i, j, k) = data(i, j, k) * 4  * pi / k2 !* gconst
				else
					data(i, j, k) = 0.
				end if

			end do
		end do
	end do
	
	print*, my_id, "GRAVITATIONAL FIELD COMPUTED"
	
! 	if(my_id.eq.0) then
! 		print*, "GRAVITATIONAL FIELD COMPUTED"
! 	end if
	

	CALL system_clock(count(4)) ! Time for the update
	
	! Execute plan => Backward FFT
	CALL fftw_mpi_execute_dft(planBACK, data, data)	
	
	CALL system_clock(count(5)) ! Time for Backward FFT
	
	print*, my_id, "NORMALIZE PHI"
	
! 	if(my_id.eq.0) then
! 		print*, "NORMALIZE PHI"
! 	end if
	
	
	! Normalize data
	subarray = real(data) / (nx * ny * nz )
	
	print*, my_id, "WRITE PHI"
	
! 	if(my_id.eq.0) then
! 		print*, "WRITE PHI"
! 	end if
	
	! Parallel write
	CALL MPI_FILE_WRITE_ALL(fh_phi, subarray, alloc_local, MPI_REAL8, MPI_STATUS_IGNORE)	
	CALL MPI_FILE_CLOSE(fh_phi, ierr)

	print*, my_id, "WRITTEN PHI"	

! 	if(my_id.eq.0) then
! 		print*, "WRITTEN PHI"
! 	end if
	
	CALL system_clock(count(6), count_rate, count_max)
	
	

	if(my_id.eq.0) then
		write(0,'(a, f10.4, a)') 'Initialization = ', dble((count(2)-count(1)))/count_rate, ' s'
		write(0,'(a, f10.4, a)') 'Forward FFT = ', dble((count(3)-count(2)))/count_rate, ' s'
		write(0,'(a, f10.4, a)') 'Update = ', dble((count(4)-count(3)))/count_rate, ' s'
		write(0,'(a, f10.4, a)') 'Backward FFT = ', dble((count(5)-count(4)))/count_rate, ' s'
		write(0,'(a, f10.4, a)') 'Write = ', dble((count(6)-count(5)))/count_rate, ' s'
		
		write(0,'(a, f10.4, a)') 'Total execution time = ', dble((count(6)-count(1)))/count_rate, ' s'
		write(0,'(a, f10.4, a)') 'Total computation time = ', dble((count(5)-count(2)))/count_rate, ' s'
		
	end if
	
	
	
	CALL fftw_destroy_plan(planBACK)
 	CALL fftw_destroy_plan(planFOR)
	
	
	CALL MPI_FINALIZE(ierr)
	
	
	
END SUBROUTINE FIELDS
