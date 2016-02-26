!#####################################################################################################################################!
!#####################################################################################################################################!
!***********************************  3D - FINITE VOLUME SOLVER FOR EULER EQUATION    ************************************************!
!***********************************       AUTHOR: PRAKASH (AE15D015)                 ************************************************!
!***********************************       DATE  : 20 FEBRUARY, 2016                  ************************************************!
!***********************************       LAST MODIFIED: 22 FEBRUARY, 2016           ************************************************!
!#####################################################################################################################################!
!#####################################################################################################################################!

module Geometry
implicit none
		integer,parameter :: Nx = 4, Ny = 4, Nz = 4
		double precision, parameter :: gamma = 1.40d0, del_t = 1E-3
		double precision, parameter :: XUpper = 1.0d0, XLower = 0.0d0, YUpper = 1.0d0, YLower = 0.0d0, ZUpper = 1.0d0, ZLower = 0.0d0 
contains
!---------------------------------------------------------------------------------------------------------------------------!
!********************************************** Compute volume *************************************************************!
!---------------------------------------------------------------------------------------------------------------------------!

	subroutine Volume(B,H,A_iph, A_jph, A_kph,Cell_Volume)
	implicit none
		double precision, dimension(1:3) :: B, H, d1, d2, A_iph, A_jph, A_kph
		double precision :: Cell_Volume, Vol, Dot_Soln

		Vol = 0.0d0
		d2(1) = B(1) - H(1)
		d2(2) = B(2) - H(2)
		d2(3) = B(3) - H(3)

		d1(1) = A_iph(1)
		d1(2) = A_iph(2)
		d1(3) = A_iph(3)

		call dot_prod(d1,d2,Dot_Soln)
		Vol = Vol + Dot_Soln

		d1(1) = A_jph(1)
		d1(2) = A_jph(2)
		d1(3) = A_jph(3)

		call dot_prod(d1,d2,Dot_Soln)
		Vol = Vol + Dot_Soln

		d1(1) = A_kph(1)
		d1(2) = A_kph(2)
		d1(3) = A_kph(3)

		call dot_prod(d1,d2,Dot_Soln)
		Vol = Vol + Dot_Soln

		Cell_Volume = (1.0d0/3.0d0)*Vol 
	end subroutine Volume
!############################################################################################################################!



!----------------------------------------------------------------------------------------------------------------------------!
!**************************************************** Face centroids ********************************************************!
!----------------------------------------------------------------------------------------------------------------------------!

	subroutine Face_centroids(x,y,z,Fc_iph,Fc_jph,Fc_kph,Fc_imh,Fc_jmh,Fc_kmh)
	implicit none
		integer :: i,j,k
		double precision, dimension(0:Nx+2,0:Ny+2,0:Nz+2) :: x, y, z
		double precision, dimension(0:Nx+2,0:Ny+2,0:Nz+2,1:3) :: Fc_iph,Fc_jph,Fc_kph,Fc_imh,Fc_jmh,Fc_kmh

		do k = 0,Nz+1
			do j = 0,Ny+1
				do i = 0, Nx+1
					Fc_iph(i,j,k,1) = 0.250d0*(x(i+1,j,k+1) + x(i+1,j+1,k+1) + x(i+1,j+1,k) + x(i+1,j,k))
					Fc_iph(i,j,k,2) = 0.250d0*(y(i+1,j,k+1) + y(i+1,j+1,k+1) + y(i+1,j+1,k) + y(i+1,j,k))
					Fc_iph(i,j,k,3) = 0.250d0*(z(i+1,j,k+1) + z(i+1,j+1,k+1) + z(i+1,j+1,k) + z(i+1,j,k))

					Fc_jph(i,j,k,1) = 0.250d0*(x(i+1,j+1,k+1) + x(i+1,j+1,k) + x(i,j+1,k) + x(i,j+1,k+1))
					Fc_jph(i,j,k,2) = 0.250d0*(y(i+1,j+1,k+1) + y(i+1,j+1,k) + y(i,j+1,k) + y(i,j+1,k+1))
					Fc_jph(i,j,k,3) = 0.250d0*(z(i+1,j+1,k+1) + z(i+1,j+1,k) + z(i,j+1,k) + z(i,j+1,k+1))
	
					Fc_kph(i,j,k,1) = 0.250d0*(x(i+1,j+1,k+1) + x(i+1,j,k+1) + x(i,j,k+1) + x(i,j+1,k+1))
					Fc_kph(i,j,k,2) = 0.250d0*(y(i+1,j+1,k+1) + y(i+1,j,k+1) + y(i,j,k+1) + y(i,j+1,k+1))
					Fc_kph(i,j,k,3) = 0.250d0*(z(i+1,j+1,k+1) + z(i+1,j,k+1) + z(i,j,k+1) + z(i,j+1,k+1))
			
					Fc_imh(i,j,k,1) = 0.250d0*(x(i,j,k+1) + x(i,j+1,k+1) + x(i,j+1,k) + x(i,j,k))
					Fc_imh(i,j,k,2) = 0.250d0*(y(i,j,k+1) + y(i,j+1,k+1) + y(i,j+1,k) + y(i,j,k))
					Fc_imh(i,j,k,3) = 0.250d0*(z(i,j,k+1) + z(i,j+1,k+1) + z(i,j+1,k) + z(i,j,k))

					Fc_jmh(i,j,k,1) = 0.250d0*(x(i,j,k+1) + x(i,j,k) + x(i+1,j,k) + x(i+1,j,k+1))
					Fc_jmh(i,j,k,2) = 0.250d0*(y(i,j,k+1) + y(i,j,k) + y(i+1,j,k) + y(i+1,j,k+1))
					Fc_jmh(i,j,k,3) = 0.250d0*(z(i,j,k+1) + z(i,j,k) + z(i+1,j,k) + z(i+1,j,k+1))
			
					Fc_kmh(i,j,k,1) = 0.250d0*(x(i,j,k) + x(i+1,j,k) + x(i+1,j+1,k) + x(i,j+1,k))
					Fc_kmh(i,j,k,2) = 0.250d0*(y(i,j,k) + y(i+1,j,k) + y(i+1,j+1,k) + y(i,j+1,k))
					Fc_kmh(i,j,k,3) = 0.250d0*(z(i,j,k) + z(i+1,j,k) + z(i+1,j+1,k) + z(i,j+1,k))
				end do
			end do
		end do
	end subroutine Face_centroids
!############################################################################################################################!



!----------------------------------------------------------------------------------------------------------------------------!
!*********************************************** Cell centroids *************************************************************!
!----------------------------------------------------------------------------------------------------------------------------!
	subroutine Cell_centroid(A,B,C,D,E,F,G,H,Cell_center)
	implicit none
		double precision, dimension(1:3) :: A,B,C,D,E,F,G,H,Cell_center
		Cell_center(1) = 0.1250d0*(A(1) + B(1) + C(1) + D(1) + E(1) + F(1) + G(1) + H(1))	
		Cell_center(2) = 0.1250d0*(A(2) + B(2) + C(2) + D(2) + E(2) + F(2) + G(2) + H(2)) 
		Cell_center(3) = 0.1250d0*(A(3) + B(3) + C(3) + D(3) + E(3) + F(3) + G(3) + H(3))

	end subroutine Cell_centroid

!----------------------------------------------------------------------------------------------------------------------------!
!*********************************************** Cross Product **************************************************************!
!----------------------------------------------------------------------------------------------------------------------------!
	subroutine cross_prod(d1,d2,Cross_Soln)
	implicit none
		double precision, dimension(1:3) :: d1, d2, Cross_Soln	
		Cross_Soln(1) = ((d1(2)*d2(3))- (d2(2)*d1(3)))
		Cross_Soln(2) = -((d1(1)*d2(3))- (d2(1)*d1(3)))
		Cross_Soln(3) = ((d1(1)*d2(2))- (d2(1)*d1(2)))
	
	open(unit = 65, file="Cross_prod.dat")		
	write(65,'(3E30.18)')  Cross_Soln(1), Cross_Soln(2), Cross_Soln(3) 
	write (65,*)
	end subroutine cross_prod
!############################################################################################################################!




!----------------------------------------------------------------------------------------------------------------------------!
!**************************************************** Dot Product ***********************************************************!
!----------------------------------------------------------------------------------------------------------------------------!
	subroutine dot_prod(d1,d2,Dot_Soln)
	implicit none
		double precision, dimension(1:3) :: d1, d2
		double precision :: Dot_Soln 
		Dot_Soln = d1(1)*d2(1) + d1(2)*d2(2) + d1(3)*d2(3)

	open(unit = 75, file="Dot_prod.dat")
	write(75,'(1E30.18)')  Dot_Soln
	write(75,*)
	end subroutine dot_prod
end module Geometry
!############################################################################################################################!



!----------------------------------------------------------------------------------------------------------------------------!
!*************************************************** Main loop **************************************************************!
!----------------------------------------------------------------------------------------------------------------------------!

program Euler
	use Geometry
	implicit none	
	integer :: i,j,k,l
	double precision, dimension(1:3) :: d1,d2,Cross_Soln	
	double precision :: dx, dy, dz,Dot_Soln, Cell_Volume
	double precision, dimension(1:Nx, 1:Ny, 1:Nz) :: vol
	double precision :: IC_RHO, IC_UX, IC_UY, IC_UZ, IC_PR
	double precision, dimension(0:Nx+2,0:Ny+2,0:Nz+2) :: x, y, z
	double precision, dimension(0:Nx+2,0:Ny+2,0:Nz+2) :: xc,yc,zc
	double precision, dimension(1:3) :: A,B,C,D,E,F,G,H,Cell_center
	double precision, dimension(0:Nx+2,0:Ny+2,0:Nz+2) :: rho,ux,uy,uz,p	
	double precision :: BCTOP_RHO, BCTOP_UX, BCTOP_UY, BCTOP_UZ, BCTOP_PR
	double precision :: BCLEFT_RHO, BCLEFT_UX, BCLEFT_UY, BCLEFT_UZ, BCLEFT_PR
	double precision :: BCBACK_RHO, BCBACK_UX, BCBACK_UY, BCBACK_UZ, BCBACK_PR
	double precision :: BCRIGHT_RHO, BCRIGHT_UX, BCRIGHT_UY, BCRIGHT_UZ, BCRIGHT_PR
	double precision :: BCFRONT_RHO, BCFRONT_UX, BCFRONT_UY, BCFRONT_UZ, BCFRONT_PR
	double precision, dimension(0:Nx+2,0:Ny+2,0:Nz+2,1:5) :: U, Flux_x, Flux_y, Flux_z 
	double precision :: BCBOTTOM_RHO, BCBOTTOM_UX, BCBOTTOM_UY, BCBOTTOM_UZ, BCBOTTOM_PR, time
	double precision, dimension(0:Nx+2,0:Ny+2,0:Nz+2,1:3) :: n_iph,n_imh,n_jph,n_jmh,n_kph,n_kmh
	double precision, dimension(0:Nx+2,0:Ny+2,0:Nz+2,1:3) :: Fc_iph,Fc_jph,Fc_kph,Fc_imh,Fc_jmh,Fc_kmh  
	double precision, dimension(0:Nx+2,0:Ny+2,0:Nz+2,1:3) :: A_iph,A_imh,A_jph,A_jmh,A_kph,A_kmh,Vol_Centroid
	double precision, dimension(0:Nx+2,0:Ny+2,0:Nz+2,1:5) :: Flux_x_iph, Flux_x_jph, Flux_x_kph, Flux_x_imh, Flux_x_jmh, Flux_x_kmh
	double precision, dimension(0:Nx+2,0:Ny+2,0:Nz+2,1:5) :: Flux_y_iph, Flux_y_jph, Flux_y_kph, Flux_y_imh, Flux_y_jmh, Flux_y_kmh
	double precision, dimension(0:Nx+2,0:Ny+2,0:Nz+2,1:5) :: Flux_z_iph, Flux_z_jph, Flux_z_kph, Flux_z_imh, Flux_z_jmh, Flux_z_kmh
	double precision, dimension(0:Nx+1, 0:Ny+1, 0:Nz+1) :: surf_mag_iph, surf_mag_imh, surf_mag_jph, surf_mag_jmh, surf_mag_kph, surf_mag_kmh	

	
!----------------------------------------------------------------------------------------------------------------------------!
!********************************************* Grid generation (node points) ************************************************!
!----------------------------------------------------------------------------------------------------------------------------!	
	dx = (XUpper - XLower)/Nx
	dy = (YUpper - YLower)/Ny
	dz = (ZUpper - ZLower)/Nz


	do k = 0,Nz+2
		do j = 0,Ny+2
			do i = 0, Nx+2
				x(i,j,k) = (i-1)*dx
				y(i,j,k) = (j-1)*dy
				z(i,j,k) = (k-1)*dz
			end do
		end do
	end do
!############################################################################################################################!




!----------------------------------------------------------------------------------------------------------------------------!
!************************************************** Center points ***********************************************************!
!----------------------------------------------------------------------------------------------------------------------------!	
	do k = 0,Nz+1
		do j = 0,Ny+1
			do i = 0, Nx+1
				xc(i,j,k) = 0.50d0*(x(i+1,j,k) + x(i,j,k))
				yc(i,j,k) = 0.50d0*(y(i,j+1,k) + y(i,j,k))
				zc(i,j,k) = 0.50d0*(z(i,j,k+1) + z(i,j,k))
			end do
		end do
	end do
!############################################################################################################################!
	
	
	
!----------------------------------------------------------------------------------------------------------------------------!
!********************************************* Defining the vertices ********************************************************!
!----------------------------------------------------------------------------------------------------------------------------!
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
!----------------------------------------------------------------------------------------------------------------------------!
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Area and unit normal for the 6 faces <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
!----------------------------------------------------------------------------------------------------------------------------!
!************************************* Area and unit normal @ i+1/2 *********************************************************!
!----------------------------------------------------------------------------------------------------------------------------!
			 	d1(1) = C(1) - A(1)
				d1(2) = C(2) - A(2)
				d1(3) = C(3) - A(3)
	
				d2(1) = B(1) - D(1)
				d2(2) = B(2) - D(2)
				d2(3) = B(3) - D(3)

				call cross_prod(d1,d2,Cross_Soln)
				A_iph(i,j,k,1) = 0.50d0*Cross_Soln(1)
				A_iph(i,j,k,2) = 0.50d0*Cross_Soln(2)
				A_iph(i,j,k,3) = 0.50d0*Cross_Soln(3)

				surf_mag_iph(i,j,k) = (A_iph(i,j,k,1)**2 + A_iph(i,j,k,2)**2 + A_iph(i,j,k,3)**2)**0.50d0 !*** magnitude of Area vector ***!
			
				n_iph(i,j,k,1) = A_iph(i,j,k,1)/surf_mag_iph(i,j,k)
				n_iph(i,j,k,2) = A_iph(i,j,k,2)/surf_mag_iph(i,j,k)
				n_iph(i,j,k,3) = A_iph(i,j,k,3)/surf_mag_iph(i,j,k)

!----------------------------------------------------------------------------------------------------------------------------!
!************************************** Area and unit normal @ i-1/2 ********************************************************!
!----------------------------------------------------------------------------------------------------------------------------!
				d1(1) = G(1) - E(1)
				d1(2) = G(2) - E(2)
				d1(3) = G(3) - E(3)
	
				d2(1) = F(1) - H(1)
				d2(2) = F(2) - H(2)
				d2(3) = F(3) - H(3)

				call cross_prod(d1,d2,Cross_Soln)
				A_imh(i,j,k,1) = 0.50d0*Cross_Soln(1)
				A_imh(i,j,k,2) = 0.50d0*Cross_Soln(2)
				A_imh(i,j,k,3) = 0.50d0*Cross_Soln(3)

				surf_mag_imh(i,j,k) = (A_imh(i,j,k,1)**2 + A_imh(i,j,k,2)**2 + A_imh(i,j,k,3)**2)**0.50d0 !*** magnitude of Area vector ***!
			
				n_imh(i,j,k,1) = A_imh(i,j,k,1)/surf_mag_imh(i,j,k)
				n_imh(i,j,k,2) = A_imh(i,j,k,2)/surf_mag_imh(i,j,k)
				n_imh(i,j,k,3) = A_imh(i,j,k,3)/surf_mag_imh(i,j,k)
				
!----------------------------------------------------------------------------------------------------------------------------!
!************************************* Area and unit normal @ j+1/2 *********************************************************!
!----------------------------------------------------------------------------------------------------------------------------!
			 	d1(1) = F(1) - C(1)
				d1(2) = F(2) - C(2)
				d1(3) = F(3) - C(3)
	
				d2(1) = B(1) - G(1)
				d2(2) = B(2) - G(2)
				d2(3) = B(3) - G(3)

				call cross_prod(d1,d2,Cross_Soln)
				A_jph(i,j,k,1) = 0.50d0*Cross_Soln(1)
				A_jph(i,j,k,2) = 0.50d0*Cross_Soln(2)
				A_jph(i,j,k,3) = 0.50d0*Cross_Soln(3)

				surf_mag_jph(i,j,k) = (A_jph(i,j,k,1)**2 + A_jph(i,j,k,2)**2 + A_jph(i,j,k,3)**2)**0.50d0 !*** magnitude of Area vector ***!
			
				n_jph(i,j,k,1) = A_jph(i,j,k,1)/surf_mag_jph(i,j,k)
				n_jph(i,j,k,2) = A_jph(i,j,k,2)/surf_mag_jph(i,j,k)
				n_jph(i,j,k,3) = A_jph(i,j,k,3)/surf_mag_jph(i,j,k)
				
!----------------------------------------------------------------------------------------------------------------------------!
!************************************* Area and unit normal @ j-1/2 *********************************************************!
!----------------------------------------------------------------------------------------------------------------------------!
				d1(1) = E(1) - D(1)
				d1(2) = E(2) - D(2)
				d1(3) = E(3) - D(3)
	
				d2(1) = A(1) - H(1)
				d2(2) = A(2) - H(2)
				d2(3) = A(3) - H(3)

				call cross_prod(d1,d2,Cross_Soln)
				A_jmh(i,j,k,1) = 0.50d0*Cross_Soln(1)
				A_jmh(i,j,k,2) = 0.50d0*Cross_Soln(2)
				A_jmh(i,j,k,3) = 0.50d0*Cross_Soln(3)

				surf_mag_jmh(i,j,k) = (A_jmh(i,j,k,1)**2 + A_jmh(i,j,k,2)**2 + A_jmh(i,j,k,3)**2)**0.50d0 !*** magnitude of Area vector ***!
			
				n_jmh(i,j,k,1) = A_jmh(i,j,k,1)/surf_mag_jmh(i,j,k)
				n_jmh(i,j,k,2) = A_jmh(i,j,k,2)/surf_mag_jmh(i,j,k)
				n_jmh(i,j,k,3) = A_jmh(i,j,k,3)/surf_mag_jmh(i,j,k)

!----------------------------------------------------------------------------------------------------------------------------!
!************************************ Area and unit normal @ k+1/2 **********************************************************!
!----------------------------------------------------------------------------------------------------------------------------!
			 	d1(1) = B(1) - E(1)
				d1(2) = B(2) - E(2)
				d1(3) = B(3) - E(3)
	
				d2(1) = F(1) - A(1)
				d2(2) = F(2) - A(2)
				d2(3) = F(3) - A(3)

				call cross_prod(d1,d2,Cross_Soln)
				A_kph(i,j,k,1) = 0.50d0*Cross_Soln(1)
				A_kph(i,j,k,2) = 0.50d0*Cross_Soln(2)
				A_kph(i,j,k,3) = 0.50d0*Cross_Soln(3)

				surf_mag_kph(i,j,k) = (A_kph(i,j,k,1)**2 + A_kph(i,j,k,2)**2 + A_kph(i,j,k,3)**2)**0.50d0 !*** magnitude of Area vector ***!
			
				n_kph(i,j,k,1) = A_kph(i,j,k,1)/surf_mag_kph(i,j,k)
				n_kph(i,j,k,2) = A_kph(i,j,k,2)/surf_mag_kph(i,j,k)
				n_kph(i,j,k,3) = A_kph(i,j,k,3)/surf_mag_kph(i,j,k)
			
!----------------------------------------------------------------------------------------------------------------------------!
!*********************************** Area and unit normal @ k-1/2 ***********************************************************!
!----------------------------------------------------------------------------------------------------------------------------!
				d1(1) = C(1) - H(1)
				d1(2) = C(2) - H(2)
				d1(3) = C(3) - H(3)
	
				d2(1) = G(1) - D(1)
				d2(2) = G(2) - D(2)
				d2(3) = G(3) - D(3)

				call cross_prod(d1,d2,Cross_Soln)
				A_kmh(i,j,k,1) = 0.50d0*Cross_Soln(1)
				A_kmh(i,j,k,2) = 0.50d0*Cross_Soln(2)
				A_kmh(i,j,k,3) = 0.50d0*Cross_Soln(3)

				surf_mag_kmh(i,j,k) = (A_kmh(i,j,k,1)**2 + A_kmh(i,j,k,2)**2 + A_kmh(i,j,k,3)**2)**0.50d0 !*** magnitude of Area vector ***!
			
				n_kmh(i,j,k,1) = A_kmh(i,j,k,1)/surf_mag_kmh(i,j,k)
				n_kmh(i,j,k,2) = A_kmh(i,j,k,2)/surf_mag_kmh(i,j,k)
				n_kmh(i,j,k,3) = A_kmh(i,j,k,3)/surf_mag_kmh(i,j,k)

!##############################################################################################################################!
				
				

				call Cell_centroid(A,B,C,D,E,F,G,H,Cell_center) !************** Calculation of cell centroid ******************!
				open(unit = 35, file="Volume_centroid.dat")
				do l = 1,3
					Vol_Centroid(i,j,k,l) = Cell_center(l)
					write(35,'(3I2,1E25.10)') k,j,i,Vol_Centroid(i,j,k,l)
				end do
				write(35,*)

				
				call Volume(B,H,A_iph(i,j,k,1:3), A_jph(i,j,k,1:3), A_kph(i,j,k,1:3),Cell_Volume) !************** Calculation of cell volume ******************!
				vol(i,j,k) = Cell_Volume
				
				
				call Face_centroids(x,y,z, Fc_iph, Fc_jph, Fc_kph, Fc_imh, Fc_jmh, Fc_kmh) !************** Calculation of face centroid ******************!
	
!--------------------------------------------------------------------------------------------------------------------------------!
!*********************************** Writing the output to a file ***************************************************************!
!--------------------------------------------------------------------------------------------------------------------------------!
				open(unit = 45, file="volume.dat")
				write(45,'(3I2,1E25.10)') k,j,i,vol(i,j,k)
				

				open(unit = 55, file="face_centroid.dat")
				write(55,'(3I2,6E25.10)')k,j,i,Fc_iph(i,j,k,l), Fc_jph(i,j,k,l), Fc_kph(i,j,k,l), Fc_imh(i,j,k,l), Fc_jmh(i,j,k,l), Fc_kmh(i,j,k,l)
				

				open(unit = 85, file="Area_unit_normal.dat")
				write(85,'(3I2,7E25.10)') k,j,i,A_iph(i,j,k,1), A_iph(i,j,k,2), A_iph(i,j,k,3), surf_mag_iph(i,j,k), n_iph(i,j,k,1), n_iph(i,j,k,2), n_iph(i,j,k,3)
				write(85,'(3I2,7E25.10)') k,j,i,A_imh(i,j,k,1), A_imh(i,j,k,2), A_imh(i,j,k,3), surf_mag_imh(i,j,k), n_imh(i,j,k,1), n_imh(i,j,k,2), n_imh(i,j,k,3)
				write(85,'(3I2,7E25.10)') k,j,i,A_jph(i,j,k,1), A_jph(i,j,k,2), A_jph(i,j,k,3), surf_mag_jph(i,j,k), n_jph(i,j,k,1), n_jph(i,j,k,2), n_jph(i,j,k,3)
				write(85,'(3I2,7E25.10)') k,j,i,A_jmh(i,j,k,1), A_jmh(i,j,k,2), A_jmh(i,j,k,3), surf_mag_jmh(i,j,k), n_jmh(i,j,k,1), n_jmh(i,j,k,2), n_jmh(i,j,k,3)
				write(85,'(3I2,7E25.10)') k,j,i,A_kph(i,j,k,1), A_kph(i,j,k,2), A_kph(i,j,k,3), surf_mag_kph(i,j,k), n_kph(i,j,k,1), n_kph(i,j,k,2), n_kph(i,j,k,3)
				write(85,'(3I2,7E25.10)') k,j,i,A_kmh(i,j,k,1), A_kmh(i,j,k,2), A_kmh(i,j,k,3), surf_mag_kmh(i,j,k), n_kmh(i,j,k,1), n_kmh(i,j,k,2), n_kmh(i,j,k,3)
				write(85,*)
				write(85,*)
				write(85,*) 'For the next cell volume'
				
			end do
		end do
	end do
	close(35)
	close(45)
	close(55)
	close(65)
	close(75)
	close(85)
!#################################################################################################################################!


				
!---------------------------------------------------------------------------------------------------------------------------------!	
!********************************************** Initial conditions ***************************************************************!
!---------------------------------------------------------------------------------------------------------------------------------!	
	do i = 1,Nx
		do j = 1,Ny
			do k = 1, Nz
				rho(i,j,k) = IC_RHO 
				ux(i,j,k)  = IC_UX 
				uy(i,j,k)  = IC_UY 
				uz(i,j,k)  = IC_UZ 
				p(i,j,k)   = IC_PR 
			end do
		end do
	end do


	do i = 1,Nx
		do j = 1,Ny
			do k = 1, Nz 
				U(i,j,k,1) = rho(i,j,k)
				U(i,j,k,2) = rho(i,j,k) * ux(i,j,k)
				U(i,j,k,3) = rho(i,j,k) * uy(i,j,k)
				U(i,j,k,4) = rho(i,j,k) * uz(i,j,k)
				U(i,j,k,5) = (p(i,j,k)/(gamma-1.0d0)) + ((rho(i,j,k)*(ux(i,j,k)**2 + uy(i,j,k)**2 + uz(i,j,k)**2))/2.0d0)
			end do
		end do
	end do

!############################################## End of Initial conditions ###########################################################!



!----------------------------------------------------------------------------------------------------------------------------------!
!********************************************** Boundary conditions ***************************************************************!
!----------------------------------------------------------------------------------------------------------------------------------!		
	time = 0.0d0	
	do while(time <= 1.0) !********* Time loop *********!
		do k = 0,Nz+1
			do j = 0,Ny+1
				rho(0,j,k) = BCLEFT_RHO
				ux(0,j,k) = BCLEFT_UX 
				uy(0,j,k) = BCLEFT_UY 
				uz(0,j,k) = BCLEFT_UZ 
				p(0,j,k) = BCLEFT_PR 
		
				rho(Nx+1,j,k) = BCRIGHT_RHO 
				ux(Nx+1,j,k) = BCRIGHT_UX 
				uy(Nx+1,j,k) = BCRIGHT_UY 
				uz(Nx+1,j,k) = BCRIGHT_UZ 
				p(Nx+1,j,k) = BCRIGHT_PR 
			end do
		end do


		do k = 0,Nz+1
			do i = 0,Nx+1
				rho(i,0,k) = BCBOTTOM_RHO
				ux(i,0,k) = BCBOTTOM_UX 
				uy(i,0,k) = BCBOTTOM_UY 
				uz(i,0,k) = BCBOTTOM_UZ 
				p(i,0,k) = BCBOTTOM_PR 
		
				rho(i,Ny+1,k) = BCTOP_RHO
				ux(i,Ny+1,k) = BCTOP_UX 
				uy(i,Ny+1,k) = BCTOP_UY 
				uz(i,Ny+1,k) = BCTOP_UZ 
				p(i,Ny+1,k) = BCTOP_PR
			end do
		end do


		do j = 0,Ny+1
			do i = 0,Nx+1
				rho(i,j,0) = BCBACK_RHO 
				ux(i,j,0) = BCBACK_UX 
				uy(i,j,0) = BCBACK_UY 
				uz(i,j,0) = BCBACK_UZ 
				p(i,j,0) = BCBACK_PR 
		
				rho(i,j,Nz+1) = BCFRONT_RHO
				ux(i,j,Nz+1) = BCFRONT_UX 
				uy(i,j,Nz+1) = BCFRONT_UY 
				uz(i,j,Nz+1) = BCFRONT_UZ 
				p(i,j,Nz+1) = BCFRONT_PR 
			end do
		end do
!########################################### End of Boundary conditions #############################################################!	



!-----------------------------------------------------------------------------------------------------------------------------------!
!*********************************** Calculation of fluxes including the boundaries ************************************************!	
!-----------------------------------------------------------------------------------------------------------------------------------!
		do k = 0,Nz+1
			do j = 0,Ny+1
				do i = 0,Nx+1
					Flux_x(i,j,k,1) = rho(i,j,k) * ux(i,j,k)
					Flux_x(i,j,k,2) = (rho(i,j,k) * ux(i,j,k)**2) + p(i,j,k)
					Flux_x(i,j,k,3) = rho(i,j,k) * ux(i,j,k) * uy(i,j,k)
					Flux_x(i,j,k,4) = rho(i,j,k) * ux(i,j,k) * uz(i,j,k)
					Flux_x(i,j,k,5) = (((p(i,j,k)/(gamma-1.0d0)) + ((rho(i,j,k)*(ux(i,j,k)**2 + uy(i,j,k)**2 + uz(i,j,k)**2))/2.0d0)) + p(i,j,k)) * ux(i,j,k)

					Flux_y(i,j,k,1) = rho(i,j,k) * uy(i,j,k)
					Flux_y(i,j,k,2) = rho(i,j,k) * uy(i,j,k) * ux(i,j,k)
					Flux_y(i,j,k,3) = (rho(i,j,k) * uy(i,j,k)**2) + p(i,j,k)          
					Flux_y(i,j,k,4) = rho(i,j,k) * uy(i,j,k) * uz(i,j,k)
					Flux_y(i,j,k,5) = (((p(i,j,k)/(gamma-1.0d0)) + ((rho(i,j,k)*(ux(i,j,k)**2 + uy(i,j,k)**2 + uz(i,j,k)**2))/2.0d0)) + p(i,j,k)) * uy(i,j,k)

					Flux_z(i,j,k,1) = rho(i,j,k) * uy(i,j,k)
					Flux_z(i,j,k,2) = rho(i,j,k) * uz(i,j,k) * ux(i,j,k)
					Flux_z(i,j,k,3) = rho(i,j,k) * uz(i,j,k) * uy(i,j,k)          
					Flux_z(i,j,k,4) = (rho(i,j,k) * uz(i,j,k)**2) + p(i,j,k)  
					Flux_z(i,j,k,5) = (((p(i,j,k)/(gamma-1.0d0)) + ((rho(i,j,k)*(ux(i,j,k)**2 + uy(i,j,k)**2 + uz(i,j,k)**2))/2.0d0)) + p(i,j,k)) * uz(i,j,k)
				end do
			end do
		end do
!#####################################################################################################################################!



!-------------------------------------------------------------------------------------------------------------------------------------!
!************************************************ Fluxes at the interface ************************************************************!
!-------------------------------------------------------------------------------------------------------------------------------------!
		do l = 1,5
			do k = 1,Nz
				do j = 1,Ny
					do i = 1,Nx
						Flux_x_iph(i,j,k,l) = 0.50d0*(Flux_x(i,j,k,l) + Flux_x(i+1,j,k,l))
						Flux_x_imh(i,j,k,l) = 0.50d0*(Flux_x(i,j,k,l) + Flux_x(i-1,j,k,l))
						Flux_x_jph(i,j,k,l) = 0.50d0*(Flux_x(i,j,k,l) + Flux_x(i,j+1,k,l))
						Flux_x_jmh(i,j,k,l) = 0.50d0*(Flux_x(i,j,k,l) + Flux_x(i,j-1,k,l))
						Flux_x_kph(i,j,k,l) = 0.50d0*(Flux_x(i,j,k,l) + Flux_x(i,j,k+1,l))
						Flux_x_kmh(i,j,k,l) = 0.50d0*(Flux_x(i,j,k,l) + Flux_x(i,j,k-1,l))	

						Flux_y_iph(i,j,k,l) = 0.50d0*(Flux_y(i,j,k,l) + Flux_y(i+1,j,k,l))
						Flux_y_imh(i,j,k,l) = 0.50d0*(Flux_y(i,j,k,l) + Flux_y(i-1,j,k,l))
						Flux_y_jph(i,j,k,l) = 0.50d0*(Flux_y(i,j,k,l) + Flux_y(i,j+1,k,l))
						Flux_y_jmh(i,j,k,l) = 0.50d0*(Flux_y(i,j,k,l) + Flux_y(i,j-1,k,l))
						Flux_y_kph(i,j,k,l) = 0.50d0*(Flux_y(i,j,k,l) + Flux_y(i,j,k+1,l))
						Flux_y_kmh(i,j,k,l) = 0.50d0*(Flux_y(i,j,k,l) + Flux_y(i,j,k-1,l))

						Flux_z_iph(i,j,k,l) = 0.50d0*(Flux_z(i,j,k,l) + Flux_z(i+1,j,k,l))
						Flux_z_imh(i,j,k,l) = 0.50d0*(Flux_z(i,j,k,l) + Flux_z(i-1,j,k,l))
						Flux_z_jph(i,j,k,l) = 0.50d0*(Flux_z(i,j,k,l) + Flux_z(i,j+1,k,l))
						Flux_z_jmh(i,j,k,l) = 0.50d0*(Flux_z(i,j,k,l) + Flux_z(i,j-1,k,l))
						Flux_z_kph(i,j,k,l) = 0.50d0*(Flux_z(i,j,k,l) + Flux_z(i,j,k+1,l))
						Flux_z_kmh(i,j,k,l) = 0.50d0*(Flux_z(i,j,k,l) + Flux_z(i,j,k-1,l))	
					end do
				end do
			end do
		end do
!######################################################################################################################################!



!--------------------------------------------------------------------------------------------------------------------------------------!
!************************************************* Updating next time step ************************************************************!
!--------------------------------------------------------------------------------------------------------------------------------------!

		do l = 1,5
			do k = 1,Nz
				do j = 1,Ny
					do i = 1,Nx
						U(i,j,k,l) = U(i,j,k,l) - ((del_t/Vol(i,j,k))*(((Flux_x_iph(i,j,k,l)*n_iph(i,j,k,1) + Flux_y_iph(i,j,k,l)*n_iph(i,j,k,2) + Flux_z_iph(i,j,k,l)*n_iph(i,j,k,3))*surf_mag_iph(i,j,k)) - ((Flux_x_imh(i,j,k,l)*n_imh(i,j,k,1) + Flux_y_imh(i,j,k,l)*n_imh(i,j,k,2) + Flux_z_imh(i,j,k,l)*n_imh(i,j,k,3))*surf_mag_imh(i,j,k)) + ((Flux_x_jph(i,j,k,l)*n_jph(i,j,k,1) + Flux_y_jph(i,j,k,l)*n_jph(i,j,k,2) + Flux_z_jph(i,j,k,l)*n_jph(i,j,k,3))*surf_mag_jph(i,j,k)) - ((Flux_x_jmh(i,j,k,l)*n_jmh(i,j,k,1) + Flux_y_jmh(i,j,k,l)*n_jmh(i,j,k,2) + Flux_z_jmh(i,j,k,l)*n_jmh(i,j,k,3))*surf_mag_jmh(i,j,k)) + ((Flux_x_kph(i,j,k,l)*n_kph(i,j,k,1) + Flux_y_kph(i,j,k,l)*n_kph(i,j,k,2) + Flux_z_kph(i,j,k,l)*n_kph(i,j,k,3))*surf_mag_kph(i,j,k)) - ((Flux_x_kmh(i,j,k,l)*n_kmh(i,j,k,1) + Flux_y_kmh(i,j,k,l)*n_kmh(i,j,k,2) + Flux_z_kmh(i,j,k,l)*n_kmh(i,j,k,3))*surf_mag_kmh(i,j,k))))
					end do
				end do
			end do
		end do
!######################################################################################################################################!



!--------------------------------------------------------------------------------------------------------------------------------------!
!*************************************** Retrieving the primitive variables ***********************************************************!
!--------------------------------------------------------------------------------------------------------------------------------------!
		do k = 1,Nz
			do j = 1,Ny
				do i = 1, Nx 
					rho(i,j,k) = U(i,j,k,1) 
					ux(i,j,k)  = U(i,j,k,2)/U(i,j,k,1)  
					uy(i,j,k)  = U(i,j,k,3)/U(i,j,k,1) 
					uz(i,j,k)  = U(i,j,k,4)/U(i,j,k,1)
					p(i,j,k)   = (gamma-1.0d0)*(U(i,j,k,5) - ((rho(i,j,k)*(ux(i,j,k)**2 + uy(i,j,k)**2 + uz(i,j,k)**2))/2.0d0))  
				end do
			end do
		end do
!#######################################################################################################################################!
	
		time = time + del_t 
	
	end do !********* End of while loop ***********!

end program Euler




