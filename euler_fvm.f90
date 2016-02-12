module variables
	integer :: i,j,k	
	real*4 :: dx, dy, dz
	integer :: Nx = 10, Ny = 10, Nz = 10
	real*4 :: g = 1.4, XUpper, XLower, YUpper, YLower, ZUpper, ZLower, d4 
	real*4, dimension(0:Nx+2,0:Ny+2,0:Nz+2) :: ds
	real*4, dimension(0:Nx+2,0:Ny+2,0:Nz+2) :: x, y, z
	real*4, dimension(0:Nx+2,0:Ny+2,0:Nz+2) :: xc,yc,zc
	real*4, dimension(0:Nx+2,0:Ny+2,0:Nz+2) :: rho,u,v,w,pr,c	
	real*4, dimension(0:Nx+2,0:Ny+2,0:Nz+2) :: u1,u2,u3,u4,u5
 	real*4, dimension(0:Nx+2,0:Ny+2,0:Nz+2) :: f1,f2,f3,f4,f5
	real*4, dimension(0:Nx+2,0:Ny+2,0:Nz+2) :: g1,g2,g3,g4,g5
	real*4, dimension(0:Nx+2,0:Ny+2,0:Nz+2) :: h1,h2,h3,h4,h5
	real*4, dimension(1:3) :: A,B,C,D,E,F,G,H,d1,d2,d3
	real*4, dimension(0:Nx+2,0:Ny+2,0:Nz+2,1:3) :: A_iph,A_imh,A_jph,A_jmh,A_kph,A_kmh
	real*4, dimension(0:Nx+2,0:Ny+2,0:Nz+2,1:3) :: n_iph,n_imh,n_jph,n_jmh,n_kph,n_kmh
	real*4 :: A_IPH, A_IMH, A_JPH, A_JMH, A_KPH, A_KMH
	!real*8, dimension(0:Nx+2,0:Ny+2,0:Nz+2) :: s1,s2,s3,s4,s5,s6
	!real*8, dimension(0:Nx+2,0:Ny+2,0:Nz+2) :: n1,n2,n3,n4,n5,n6
end module variables



module var_declare
	use variables
	implicit none
	u1(i,j,k) = rho(i,j,k)
	u2(i,j,k) = rho(i,j,k) * u(i,j,k)
	u3(i,j,k) = rho(i,j,k) * v(i,j,k)
	u4(i,j,k) = rho(i,j,k) * w(i,j,k)
	u5(i,j,k) = (p(i,j,k)/(g-1.0)) + ((rho(i,j,k)*(u(i,j,k)**2 + v(i,j,k)**2 + w(i,j,k)**2))/2.0)

	f1(i,j,k) = rho(i,j,k) * u(i,j,k)
	f2(i,j,k) = (rho(i,j,k) * u(i,j,k)**2) + p(i,j,k))
	f3(i,j,k) = rho(i,j,k) * u(i,j,k) * v(i,j,k)
	f4(i,j,k) = rho(i,j,k) * u(i,j,k) * w(i,j,k)
	f5(i,j,k) = (u5(i,j,k) + p(i,j,k)) * u(i,j,k)

	g1(i,j,k) = rho(i,j,k) * v(i,j,k)
	g2(i,j,k) = rho(i,j,k) * v(i,j,k) * u(i,j,k)
	g3(i,j,k) = (rho(i,j,k) * v(i,j,k)**2) + p(i,j,k))          
	g4(i,j,k) = rho(i,j,k) * v(i,j,k) * w(i,j,k)
	g5(i,j,k) = (u5(i,j,k) + p(i,j,k)) * v(i,j,k)

	h1(i,j,k) = rho(i,j,k) * v(i,j,k)
	h2(i,j,k) = rho(i,j,k) * w(i,j,k) * u(i,j,k)
	h3(i,j,k) = rho(i,j,k) * w(i,j,k) * v(i,j,k)          
	h4(i,j,k) = (rho(i,j,k) * w(i,j,k)**2) + p(i,j,k))  
	h5(i,j,k) = (u5(i,j,k) + p(i,j,k)) * w(i,j,k)

	!******************** del_x, del_y, del_z **********************!
	
	dx = (XUpper - XLower)/Nx
	dy = (YUpper - YLower)/Ny
	dz = (ZUpper - ZLower)/Nz

	
end module var_declare



program euler
	use variables
	use var_declare
	implicit none	
	
	!******************* Grid generation (node points) ******************!
	
	do i = 0,Nx+2
		do j = 0,Ny+2
			do k = 0, Nz+2
				x(i,j,k) = (i-1)*dx
				y(i,j,k) = (j-1)*dy
				z(i,j,k) = (k-1)*dz
			end do
		end do
	end do


	!************************ Center points *****************************!
	
	do i = 0,Nx+1
		do j = 0,Ny+1
			do k = 0, Nz+1
				xc(i,j,k) = 0.5*(x(i+1,j,k) + x(i,j,k))
				yc(i,j,k) = 0.5*(y(i,j+1,k) + y(i,j,k))
				zc(i,j,k) = 0.5*(z(i,j,k+1) + z(i,j,k))
			end do
		end do
	end do

	
	!********************* Defining the vertices ************************!
	do i = 0,Nx+1
		do j = 0,Ny+1
			do k = 0, Nz+1
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
			!******* Area and unit normal @ i+1/2 *******!
			 	d1(1) = C(1) - A(1)
				d1(2) = C(2) - A(2)
				d1(3) = C(3) - A(3)
	
				d2(1) = B(1) - D(1)
				d2(2) = B(2) - D(2)
				d2(3) = B(3) - D(3)

				call cross_prod(d1,d2)
				A_iph(i,j,k,1) = 0.5*d3(1)
				A_iph(i,j,k,2) = 0.5*d3(2)
				A_iph(i,j,k,3) = 0.5*d3(3)

				A_IPH = (A_iph(i,j,k,1)**2 + A_iph(i,j,k,2)**2 + A_iph(i,j,k,3)**2)**0.5 !***magnitude***!
			
				n_iph(i,j,k,1) = A_iph(i,j,k,1)/A_IPH
				n_iph(i,j,k,2) = A_iph(i,j,k,2)/A_IPH
				n_iph(i,j,k,3) = A_iph(i,j,k,3)/A_IPH

			!******* Area and unit normal @ i-1/2 *******!
				d1(1) = G(1) - E(1)
				d1(2) = G(2) - E(2)
				d1(3) = G(3) - E(3)
	
				d2(1) = F(1) - H(1)
				d2(2) = F(2) - H(2)
				d2(3) = F(3) - H(3)

				call cross_prod(d1,d2)
				A_imh(i,j,k,1) = 0.5*d3(1)
				A_imh(i,j,k,2) = 0.5*d3(2)
				A_imh(i,j,k,3) = 0.5*d3(3)

				A_IMH = (A_imh(i,j,k,1)**2 + A_imh(i,j,k,2)**2 + A_imh(i,j,k,3)**2)**0.5 !***magnitude***!
			
				n_imh(i,j,k,1) = A_imh(i,j,k,1)/A_IMH
				n_imh(i,j,k,2) = A_imh(i,j,k,2)/A_IMH
				n_imh(i,j,k,3) = A_imh(i,j,k,3)/A_IMH

			!******* Area and unit normal @ j+1/2 *******!
			 	d1(1) = F(1) - C(1)
				d1(2) = F(2) - C(2)
				d1(3) = F(3) - C(3)
	
				d2(1) = B(1) - G(1)
				d2(2) = B(2) - G(2)
				d2(3) = B(3) - G(3)

				call cross_prod(d1,d2)
				A_jph(i,j,k,1) = 0.5*d3(1)
				A_jph(i,j,k,2) = 0.5*d3(2)
				A_jph(i,j,k,3) = 0.5*d3(3)

				A_JPH = (A_jph(i,j,k,1)**2 + A_jph(i,j,k,2)**2 + A_jph(i,j,k,3)**2)**0.5 !***magnitude***!
			
				n_jph(i,j,k,1) = A_jph(i,j,k,1)/A_JPH
				n_jph(i,j,k,2) = A_jph(i,j,k,2)/A_JPH
				n_jph(i,j,k,3) = A_jph(i,j,k,3)/A_JPH
			
			!******* Area and unit normal @ j-1/2 *******!
				d1(1) = E(1) - D(1)
				d1(2) = E(2) - D(2)
				d1(3) = E(3) - D(3)
	
				d2(1) = A(1) - H(1)
				d2(2) = A(2) - H(2)
				d2(3) = A(3) - H(3)

				call cross_prod(d1,d2)
				A_jmh(i,j,k,1) = 0.5*d3(1)
				A_jmh(i,j,k,2) = 0.5*d3(2)
				A_jmh(i,j,k,3) = 0.5*d3(3)

				A_JMH = (A_jmh(i,j,k,1)**2 + A_jmh(i,j,k,2)**2 + A_jmh(i,j,k,3)**2)**0.5 !***magnitude***!
			
				n_jmh(i,j,k,1) = A_jmh(i,j,k,1)/A_JMH
				n_jmh(i,j,k,2) = A_jmh(i,j,k,2)/A_JMH
				n_jmh(i,j,k,3) = A_jmh(i,j,k,3)/A_JMH	
				
			!******* Area and unit normal @ k+1/2 *******!
			 	d1(1) = B(1) - E(1)
				d1(2) = B(2) - E(2)
				d1(3) = B(3) - E(3)
	
				d2(1) = F(1) - A(1)
				d2(2) = F(2) - A(2)
				d2(3) = F(3) - A(3)

				call cross_prod(d1,d2)
				A_kph(i,j,k,1) = 0.5*d3(1)
				A_kph(i,j,k,2) = 0.5*d3(2)
				A_kph(i,j,k,3) = 0.5*d3(3)

				A_KPH = (A_kph(i,j,k,1)**2 + A_kph(i,j,k,2)**2 + A_kph(i,j,k,3)**2)**0.5 !***magnitude***!
			
				n_kph(i,j,k,1) = A_kph(i,j,k,1)/A_KPH
				n_kph(i,j,k,2) = A_kph(i,j,k,2)/A_KPH
				n_kph(i,j,k,3) = A_kph(i,j,k,3)/A_KPH

			!******* Area and unit normal @ k-1/2 *******!
				d1(1) = C(1) - H(1)
				d1(2) = C(2) - H(2)
				d1(3) = C(3) - H(3)
	
				d2(1) = G(1) - D(1)
				d2(2) = G(2) - D(2)
				d2(3) = G(3) - D(3)

				call cross_prod(d1,d2)
				A_kmh(i,j,k,1) = 0.5*d3(1)
				A_kmh(i,j,k,2) = 0.5*d3(2)
				A_kmh(i,j,k,3) = 0.5*d3(3)

				A_KMH = (A_kmh(i,j,k,1)**2 + A_kmh(i,j,k,2)**2 + A_kmh(i,j,k,3)**2)**0.5 !***magnitude***!
			
				n_kmh(i,j,k,1) = A_kmh(i,j,k,1)/A_KMH
				n_kmh(i,j,k,2) = A_kmh(i,j,k,2)/A_KMH
				n_kmh(i,j,k,3) = A_kmh(i,j,k,3)/A_KMH
			end do
		end do
	end do

end program euler







subroutine cross_prod(d1,d2)
	use variables
	use var_declare
	implicit none
		
	d3(1) = ((d1(2)*d2(3))- (d2(2)*d1(3)))
	d3(2) = -((d1(1)*d2(3))- (d2(1)*d1(3)))
	d3(3) = ((d1(1)*d2(2))- (d2(1)*d1(2)))
		
end subroutine cross_prod




