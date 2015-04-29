SUBROUTINE PROJECT_DENSITY
	
	USE DIMENSION, only: nx, ny, nz, nparticles
	USE VECTOR, only: x, y, z, vx, vy, vz, mass
	USE MATRIX, only: rho3d
	USE SCALAR, only: step
	
	IMPLICIT NONE
		
	REAL*8 :: didx, didy, didz, didxp,  didyp,  didzp
	REAL*8 :: d1, d2, d3, d4, d5, d6, d7, d8
	
	INTEGER :: particle
	INTEGER :: idx, idy, idz, idxp, idyp, idzp
	INTEGER :: ierror
	
	do particle = 1, nparticles
		
		idx = int( x(particle) + 0.5 )
		idy = int( y(particle) + 0.5 )
		idz = int( z(particle) + 0.5 )
		
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
		didx = x(particle) - float(idx) + 0.5
		didy = y(particle) - float(idy) + 0.5
		didz = z(particle) - float(idz) + 0.5
		
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
		
		
		rho3d(idx,  idy,  idz)  = rho3d(idx,  idy,  idz)  + d1 * mass(particle)
		rho3d(idxp, idy,  idz)  = rho3d(idxp, idy,  idz)  + d2 * mass(particle)
		rho3d(idx,  idyp, idz)  = rho3d(idx,  idyp, idz)  + d3 * mass(particle)
		rho3d(idxp, idyp, idz)  = rho3d(idxp, idyp, idz)  + d4 * mass(particle)
		rho3d(idx,  idy,  idzp) = rho3d(idx,  idy,  idzp) + d5 * mass(particle)
		rho3d(idxp, idy,  idzp) = rho3d(idxp, idy,  idzp) + d6 * mass(particle)
		rho3d(idx,  idyp, idzp) = rho3d(idx,  idyp, idzp) + d7 * mass(particle)
		rho3d(idxp, idyp, idzp) = rho3d(idxp, idyp, idzp) + d8 * mass(particle)
		
	end do
	
	CALL WRITE_DENSITY(nx, ny, nz, rho3d, step)
	
END SUBROUTINE PROJECT_DENSITY