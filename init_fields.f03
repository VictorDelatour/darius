SUBROUTINE INIT_FIELDS
	
	USE DIMENSION, only : nx, ny, nz  ! Contains dimension-related variables
	USE MATRIX, only : rho3d, phi3d	  ! Contains the variables rho3d and phi3d
	USE SCALAR, only : i, j, k, pi
	
	IMPLICIT NONE
	
	!!! SIMPLY CALL ADD_EXPONENTIAL
	!!! Give (0,0,0) as the center you want
	!!! Data should be given as
	
	REAL*8, DIMENSION(3) :: center, pos ! Center of the gaussian
	REAL*8 :: beta, sigma 				! Factor in front of the gaussian and width of the gaussian
	
	
	! Initialize variables
	! Center of the gaussian at the middle of the box
	center(1) = .5 * (1+nx) ! Check with visualisation if this is correct
	center(2) = .5 * (1+ny)
	center(3) = .5 * (1+nz)
	
	beta = 1
	sigma = .5
	
	! 
	do i = 1, nx
		
		pos(1) = 2*(i - center(1))/nx
		
		do j = 1, ny
			
			pos(2) = 2*(j - center(2))/ny
			
			do k = 1, nz
				
				pos(3) = 2*(k - center(3))/nz
				rho3d(i, j, k) = beta * exp( -( pos(1)**2 + pos(2)**2 + pos(3)**2 ) / (sigma**2) )
				!rho3d(i, j, k) = sin(4*pi*i/nx)
				
			end do 
		end do 
	end do 
	
END SUBROUTINE INIT_FIELDS