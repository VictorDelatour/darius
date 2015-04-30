SUBROUTINE UPDATE_PARTICLES
	
	USE DIMENSION, ONLY: nx, ny, nz, nparticles
	USE VECTOR, ONLY: x, y, z, vx, vy, vz, mass
	USE MATRIX, ONLY: rho3d, phi3d, gradx, grady, gradz
	USE SCALAR, ONLY: step, deltat
	
	IMPLICIT NONE
		
	REAL*8 :: didx, didy, didz, didxp,  didyp,  didzp
	REAL*8 :: d1, d2, d3, d4, d5, d6, d7, d8
	REAL*8 :: force_x, force_y, force_z
	REAL*8 :: x_part, y_part, z_part
	REAL*8 :: vx_part, vy_part, vz_part
	
	INTEGER :: particle
	INTEGER :: idx, idy, idz, idxp, idyp, idzp
	INTEGER :: ierror
	
	deltat = 1.0
	
	write(*,*) "Updating particles"
	
	do particle = 1, nparticles
		
		! Local variable, we're going to access it a lot!
		x_part = x(particle)
		y_part = y(particle)
		z_part = z(particle)
		vx_part = vx(particle)
		vy_part = vy(particle)
		vz_part = vz(particle)
		
		
		idx = int( x_part + 0.5 )
		idy = int( y_part + 0.5 )
		idz = int( z_part + 0.5 )
		
		if(idx .eq. 0) idx = nx
		if(idy .eq. 0) idy = ny
		if(idz .eq. 0) idz = nz
		
		idxp = idx + 1
		idyp = idy + 1
		idzp = idz + 1
		
		if(idx .eq. nx) idxp = 1
		if(idy .eq. ny) idyp = 1
		if(idz .eq. nz) idzp = 1
		
		! Computes the distance from the middle point of the cell?
		didx = x_part - float(idx) + 0.5
		didy = y_part - float(idy) + 0.5
		didz = z_part - float(idz) + 0.5
		
		didxp = 1.0 - didx
		didyp = 1.0 - didy
		didzp = 1.0 - didz
				
		! Computes the weights for trilinear interpolation
		d1 = didx 	* didy 	* didz
		d2 = didxp 	* didy 	* didz
		d3 = didx 	* didyp * didz
		d4 = didxp 	* didyp * didz
		d5 = didx 	* didy 	* didzp
		d6 = didxp 	* didy 	* didzp
		d7 = didx 	* didyp * didzp
		d8 = didxp 	* didyp * didzp
		
		force_x = 0
		force_y = 0
		force_z = 0
		
		force_x = force_x + d1 * gradx(idx, idy, idz)+ d2 * gradx(idxp, idy, idz) + d3 * gradx(idx, idyp, idz) + d4 * gradx(idxp, idyp, idz)
		force_x = force_x + d5 * gradx(idx, idy, idzp) + d6 * gradx(idxp, idy, idzp)+ d7 * gradx(idx, idyp, idzp) + d8 * gradx(idxp, idyp, idzp)
		
		force_y = force_y + d1 * grady(idx, idy, idz)+ d2 * grady(idxp, idy, idz) + d3 * grady(idx, idyp, idz) + d4 * grady(idxp, idyp, idz)
		force_y = force_y + d5 * grady(idx, idy, idzp) + d6 * grady(idxp, idy, idzp)+ d7 * grady(idx, idyp, idzp) + d8 * grady(idxp, idyp, idzp)
		
		force_z = force_z + d1 * gradz(idx, idy, idz)+ d2 * gradz(idxp, idy, idz) + d3 * gradz(idx, idyp, idz) + d4 * gradz(idxp, idyp, idz)
		force_z = force_z + d5 * gradz(idx, idy, idzp) + d6 * gradz(idxp, idy, idzp)+ d7 * gradz(idx, idyp, idzp) + d8 * gradz(idxp, idyp, idzp)
		
		! Leap-frog scheme
		
		vx_part = vx_part + force_x * deltat
		vy_part = vy_part + force_y * deltat
		vz_part = vz_part + force_z * deltat
		
		x_part = x_part + vx_part * deltat
		y_part = y_part + vy_part * deltat
		z_part = z_part + vz_part * deltat
		
		! Take care of the boundary conditions
		
		if( x_part < 1) then
			x_part = x_part + (nx-1)
		else if( x_part > nx) then
			x_part = x_part - (nx-1)
		end if
		
		if( y_part < 1) then
			y_part = y_part + (ny-1)
		else if( y_part > ny) then
			y_part = y_part - (ny-1)
		end if
		
		if( z_part < 1) then
			z_part = z_part + (nz-1)
		else if( z_part > nz) then
			z_part = z_part - (nz-1)
		end if
		
		vx(particle) = vx_part
		vy(particle) = vy_part
		vz(particle) = vz_part
		
		x(particle) = x_part
		y(particle) = y_part
		z(particle) = z_part
		
			
	end do
	
	write(*,*) "Done"
	
	step = step + 1
	
	
END SUBROUTINE UPDATE_PARTICLES