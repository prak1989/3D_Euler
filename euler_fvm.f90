module Variables
implicit none
contains


!************************************* Compute volume ****************************************!
	subroutine Volume(B,H,A_iph, A_jph, A_kph,Cell_Volume)
	implicit none
		double precision, dimension(1:3) :: B, H, d1, d2, A_iph, A_jph, A_kph
		double precision :: Cell_Volume, Vol, d4

		Vol = 0.0d0
		d2(1) = B(1) - H(1)
		d2(2) = B(2) - H(2)
		d2(3) = B(3) - H(3)

		d1(1) = A_iph(1)
		d1(2) = A_iph(2)
		d1(3) = A_iph(3)

		call dot_prod(d1,d2,d4)
		Vol = Vol + d4

		d1(1) = A_jph(1)
		d1(2) = A_jph(2)
		d1(3) = A_jph(3)

		call dot_prod(d1,d2,d4)
		Vol = Vol + d4

		d1(1) = A_kph(1)
		d1(2) = A_kph(2)
		d1(3) = A_kph(3)

		call dot_prod(d1,d2,d4)
		Vol = Vol + d4

		Cell_Volume = (1.0d0/3.0d0)*Vol 
	end subroutine Volume



!************************************ Face centroids ****************************************!
	subroutine Face_centroids(x,y,z,Fc_iph,Fc_jph,Fc_kph,Fc_imh,Fc_jmh,Fc_kmh)
	implicit none
		integer :: i,j,k
		integer,parameter :: Nx = 10, Ny = 10, Nz = 10
		double precision, dimension(0:Nx+2,0:Ny+2,0:Nz+2) :: x, y, z
		double precision, dimension(0:Nx+2,0:Ny+2,0:Nz+2,1:3) :: Fc_iph,Fc_jph,Fc_kph,Fc_imh,Fc_jmh,Fc_kmh

		do i = 0,Nx+1
			do j = 0,Ny+1
				do k = 0, Nz+1
					Fc_iph(i,j,k,1) = 0.25*(x(i+1,j,k+1) + x(i+1,j+1,k+1) + x(i+1,j+1,k) + x(i+1,j,k))
					Fc_iph(i,j,k,2) = 0.25*(y(i+1,j,k+1) + y(i+1,j+1,k+1) + y(i+1,j+1,k) + y(i+1,j,k))
					Fc_iph(i,j,k,3) = 0.25*(z(i+1,j,k+1) + z(i+1,j+1,k+1) + z(i+1,j+1,k) + z(i+1,j,k))

					Fc_jph(i,j,k,1) = 0.25*(x(i+1,j+1,k+1) + x(i+1,j+1,k) + x(i,j+1,k) + x(i,j+1,k+1))
					Fc_jph(i,j,k,2) = 0.25*(y(i+1,j+1,k+1) + y(i+1,j+1,k) + y(i,j+1,k) + y(i,j+1,k+1))
					Fc_jph(i,j,k,3) = 0.25*(z(i+1,j+1,k+1) + z(i+1,j+1,k) + z(i,j+1,k) + z(i,j+1,k+1))
	
					Fc_kph(i,j,k,1) = 0.25*(x(i+1,j+1,k+1) + x(i+1,j,k+1) + x(i,j,k+1) + x(i,j+1,k+1))
					Fc_kph(i,j,k,2) = 0.25*(y(i+1,j+1,k+1) + y(i+1,j,k+1) + y(i,j,k+1) + y(i,j+1,k+1))
					Fc_kph(i,j,k,3) = 0.25*(z(i+1,j+1,k+1) + z(i+1,j,k+1) + z(i,j,k+1) + z(i,j+1,k+1))
			
					Fc_imh(i,j,k,1) = 0.25*(x(i,j,k+1) + x(i,j+1,k+1) + x(i,j+1,k) + x(i,j,k))
					Fc_imh(i,j,k,2) = 0.25*(y(i,j,k+1) + y(i,j+1,k+1) + y(i,j+1,k) + y(i,j,k))
					Fc_imh(i,j,k,3) = 0.25*(z(i,j,k+1) + z(i,j+1,k+1) + z(i,j+1,k) + z(i,j,k))

					Fc_jmh(i,j,k,1) = 0.25*(x(i,j,k+1) + x(i,j,k) + x(i+1,j,k) + x(i+1,j,k+1))
					Fc_jmh(i,j,k,2) = 0.25*(y(i,j,k+1) + y(i,j,k) + y(i+1,j,k) + y(i+1,j,k+1))
					Fc_jmh(i,j,k,3) = 0.25*(z(i,j,k+1) + z(i,j,k) + z(i+1,j,k) + z(i+1,j,k+1))
			
					Fc_kmh(i,j,k,1) = 0.25*(x(i,j,k) + x(i+1,j,k) + x(i+1,j+1,k) + x(i,j+1,k))
					Fc_kmh(i,j,k,2) = 0.25*(y(i,j,k) + y(i+1,j,k) + y(i+1,j+1,k) + y(i,j+1,k))
					Fc_kmh(i,j,k,3) = 0.25*(z(i,j,k) + z(i+1,j,k) + z(i+1,j+1,k) + z(i,j+1,k))
				end do
			end do
		end do
	end subroutine Face_centroids
	


!******************************* Cross Product *************************************************!
	subroutine cross_prod(d1,d2,d3)
	implicit none
		double precision, dimension(1:3) :: d1, d2, d3	
		d3(1) = ((d1(2)*d2(3))- (d2(2)*d1(3)))
		d3(2) = -((d1(1)*d2(3))- (d2(1)*d1(3)))
		d3(3) = ((d1(1)*d2(2))- (d2(1)*d1(2)))
		
	end subroutine cross_prod



!******************************** Dot Product ***************************************************!
	subroutine dot_prod(d1,d2,d4)
	implicit none
		double precision, dimension(1:3) :: d1, d2
		double precision :: d4 
		d4 = d1(1)*d2(1) + d1(2)*d2(2) + d1(3)*d2(3)

	end subroutine dot_prod
end module Variables




program Euler
	use Variables
	implicit none	
	integer :: i,j,k
	double precision :: dx, dy, dz
	double precision, dimension(1:3) :: d1,d2,d3
	integer,parameter :: Nx = 10, Ny = 10, Nz = 10
	double precision, dimension(1:3) :: A,B,C,D,E,F,G,H
	double precision, dimension(1:nx, 1:ny, 1:nz) :: vol
	double precision, dimension(0:Nx+2,0:Ny+2,0:Nz+2) :: x, y, z
	double precision, dimension(0:Nx+2,0:Ny+2,0:Nz+2) :: xc,yc,zc
	double precision, dimension(0:Nx+2,0:Ny+2,0:Nz+2,1:3) :: A_iph,A_imh,A_jph,A_jmh,A_kph,A_kmh
	double precision, dimension(0:Nx+2,0:Ny+2,0:Nz+2,1:3) :: n_iph,n_imh,n_jph,n_jmh,n_kph,n_kmh
	double precision, dimension(0:Nx+2,0:Ny+2,0:Nz+2,1:3) :: Fc_iph,Fc_jph,Fc_kph,Fc_imh,Fc_jmh,Fc_kmh  
	double precision :: surf_mag_iph, surf_mag_imh, surf_mag_jph, surf_mag_jmh, surf_mag_kph, surf_mag_kmh, d4, Cell_Volume
	double precision, parameter :: gamma = 1.40d0, XUpper = 1.0d0, XLower = 0.0d0, YUpper = 0.50d0, YLower = 0.0d0, ZUpper = 0.50d0, ZLower = 0.0d0
	

	dx = (XUpper - XLower)/Nx
	dy = (YUpper - YLower)/Ny
	dz = (ZUpper - ZLower)/Nz

	!******************* Grid generation (node points) ******************!
	
	do k = 0,Nz+2
		do j = 0,Ny+2
			do i = 0, Nx+2
				x(i,j,k) = (i-1)*dx
				y(i,j,k) = (j-1)*dy
				z(i,j,k) = (k-1)*dz
			end do
		end do
	end do


	!************************ Center points *****************************!
	
	do k = 0,Nz+1
		do j = 0,Ny+1
			do i = 0, Nx+1
				xc(i,j,k) = 0.50d0*(x(i+1,j,k) + x(i,j,k))
				yc(i,j,k) = 0.50d0*(y(i,j+1,k) + y(i,j,k))
				zc(i,j,k) = 0.50d0*(z(i,j,k+1) + z(i,j,k))
			end do
		end do
	end do

	
	call Face_centroids(x,y,z, Fc_iph, Fc_jph, Fc_kph, Fc_imh, Fc_jmh, Fc_kmh)


	!********************* Defining the vertices ************************!
	do k = 0,Nz+1
		do j = 0,Ny+1
			do i = 0, Nx+1
				A(1) = x(i+1,j,k+1)
				A(2) = y(i+1,j,k+1)
				A(3) = z(i+1,j,k+1)

				B(1) = x(i+1,j+1,k+1)
				B(2) = y(i+1,j+1,k+1)
				B(3) = z(i+1,j+1,k+1)

				C(1) = x(i+1,j+1,k)
				C(2) = y(i+1,j+1,k)
				C(3) = z(i+1,j+1,k)

				D(1) = x(i+1,j,k)
				D(2) = y(i+1,j,k)
				D(3) = z(i+1,j,k)

				E(1) = x(i,j,k+1)
				E(2) = y(i,j,k+1)
				E(3) = z(i,j,k+1)

				F(1) = x(i,j+1,k+1)
				F(2) = y(i,j+1,k+1)
				F(3) = z(i,j+1,k+1)

				G(1) = x(i,j+1,k)
				G(2) = y(i,j+1,k)
				G(3) = z(i,j+1,k)

				H(1) = x(i,j,k)
				H(2) = y(i,j,k)
				H(3) = z(i,j,k)



	!***************** Area and unit normal for the 6 faces **********************!
			!********** Area and unit normal @ i+1/2 *************!
			 	d1(1) = C(1) - A(1)
				d1(2) = C(2) - A(2)
				d1(3) = C(3) - A(3)
	
				d2(1) = B(1) - D(1)
				d2(2) = B(2) - D(2)
				d2(3) = B(3) - D(3)

				call cross_prod(d1,d2,d3)
				A_iph(i,j,k,1) = 0.50d0*d3(1)
				A_iph(i,j,k,2) = 0.50d0*d3(2)
				A_iph(i,j,k,3) = 0.50d0*d3(3)

				surf_mag_iph = (A_iph(i,j,k,1)**2 + A_iph(i,j,k,2)**2 + A_iph(i,j,k,3)**2)**0.50d0 !*** magnitude of Area vector ***!
			
				n_iph(i,j,k,1) = A_iph(i,j,k,1)/surf_mag_iph
				n_iph(i,j,k,2) = A_iph(i,j,k,2)/surf_mag_iph
				n_iph(i,j,k,3) = A_iph(i,j,k,3)/surf_mag_iph

			!************ Area and unit normal @ i-1/2 ************!
				d1(1) = G(1) - E(1)
				d1(2) = G(2) - E(2)
				d1(3) = G(3) - E(3)
	
				d2(1) = F(1) - H(1)
				d2(2) = F(2) - H(2)
				d2(3) = F(3) - H(3)

				call cross_prod(d1,d2,d3)
				A_imh(i,j,k,1) = 0.50d0*d3(1)
				A_imh(i,j,k,2) = 0.50d0*d3(2)
				A_imh(i,j,k,3) = 0.50d0*d3(3)

				surf_mag_imh = (A_imh(i,j,k,1)**2 + A_imh(i,j,k,2)**2 + A_imh(i,j,k,3)**2)**0.50d0 !*** magnitude of Area vector ***!
			
				n_imh(i,j,k,1) = A_imh(i,j,k,1)/surf_mag_imh
				n_imh(i,j,k,2) = A_imh(i,j,k,2)/surf_mag_imh
				n_imh(i,j,k,3) = A_imh(i,j,k,3)/surf_mag_imh

			!************* Area and unit normal @ j+1/2 ************!
			 	d1(1) = F(1) - C(1)
				d1(2) = F(2) - C(2)
				d1(3) = F(3) - C(3)
	
				d2(1) = B(1) - G(1)
				d2(2) = B(2) - G(2)
				d2(3) = B(3) - G(3)

				call cross_prod(d1,d2,d3)
				A_jph(i,j,k,1) = 0.50d0*d3(1)
				A_jph(i,j,k,2) = 0.50d0*d3(2)
				A_jph(i,j,k,3) = 0.50d0*d3(3)

				surf_mag_jph = (A_jph(i,j,k,1)**2 + A_jph(i,j,k,2)**2 + A_jph(i,j,k,3)**2)**0.50d0 !*** magnitude of Area vector ***!
			
				n_jph(i,j,k,1) = A_jph(i,j,k,1)/surf_mag_jph
				n_jph(i,j,k,2) = A_jph(i,j,k,2)/surf_mag_jph
				n_jph(i,j,k,3) = A_jph(i,j,k,3)/surf_mag_jph
			
			!************** Area and unit normal @ j-1/2 ************!
				d1(1) = E(1) - D(1)
				d1(2) = E(2) - D(2)
				d1(3) = E(3) - D(3)
	
				d2(1) = A(1) - H(1)
				d2(2) = A(2) - H(2)
				d2(3) = A(3) - H(3)

				call cross_prod(d1,d2,d3)
				A_jmh(i,j,k,1) = 0.50d0*d3(1)
				A_jmh(i,j,k,2) = 0.50d0*d3(2)
				A_jmh(i,j,k,3) = 0.50d0*d3(3)

				surf_mag_jmh = (A_jmh(i,j,k,1)**2 + A_jmh(i,j,k,2)**2 + A_jmh(i,j,k,3)**2)**0.50d0 !*** magnitude of Area vector ***!
			
				n_jmh(i,j,k,1) = A_jmh(i,j,k,1)/surf_mag_jmh
				n_jmh(i,j,k,2) = A_jmh(i,j,k,2)/surf_mag_jmh
				n_jmh(i,j,k,3) = A_jmh(i,j,k,3)/surf_mag_jmh	
				
			!************** Area and unit normal @ k+1/2 *************!
			 	d1(1) = B(1) - E(1)
				d1(2) = B(2) - E(2)
				d1(3) = B(3) - E(3)
	
				d2(1) = F(1) - A(1)
				d2(2) = F(2) - A(2)
				d2(3) = F(3) - A(3)

				call cross_prod(d1,d2,d3)
				A_kph(i,j,k,1) = 0.50d0*d3(1)
				A_kph(i,j,k,2) = 0.50d0*d3(2)
				A_kph(i,j,k,3) = 0.50d0*d3(3)

				surf_mag_kph = (A_kph(i,j,k,1)**2 + A_kph(i,j,k,2)**2 + A_kph(i,j,k,3)**2)**0.50d0 !*** magnitude of Area vector ***!
			
				n_kph(i,j,k,1) = A_kph(i,j,k,1)/surf_mag_kph
				n_kph(i,j,k,2) = A_kph(i,j,k,2)/surf_mag_kph
				n_kph(i,j,k,3) = A_kph(i,j,k,3)/surf_mag_kph

			!************** Area and unit normal @ k-1/2 **************!
				d1(1) = C(1) - H(1)
				d1(2) = C(2) - H(2)
				d1(3) = C(3) - H(3)
	
				d2(1) = G(1) - D(1)
				d2(2) = G(2) - D(2)
				d2(3) = G(3) - D(3)

				call cross_prod(d1,d2,d3)
				A_kmh(i,j,k,1) = 0.50d0*d3(1)
				A_kmh(i,j,k,2) = 0.50d0*d3(2)
				A_kmh(i,j,k,3) = 0.50d0*d3(3)

				surf_mag_kmh = (A_kmh(i,j,k,1)**2 + A_kmh(i,j,k,2)**2 + A_kmh(i,j,k,3)**2)**0.50d0 !*** magnitude of Area vector ***!
			
				n_kmh(i,j,k,1) = A_kmh(i,j,k,1)/surf_mag_kmh
				n_kmh(i,j,k,2) = A_kmh(i,j,k,2)/surf_mag_kmh
				n_kmh(i,j,k,3) = A_kmh(i,j,k,3)/surf_mag_kmh

				call Volume(B,H,A_iph(i,j,k,1:3), A_jph(i,j,k,1:3), A_kph(i,j,k,1:3),Cell_Volume)
				vol(i,j,k) = Cell_Volume
			end do
		end do
	end do


	
	open(unit = 45, file="volume.dat")	
	do k = 1, nz
		do j = 1, ny
			do i = 1, nx
				write(45,'(1E30.18)') vol(i,j,k)
			end do
		end do
	end do
	close(45)

end program Euler






