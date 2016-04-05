module timeStep

	use modSource

	implicit none

contains

	subroutine forwardsweep(this,xp0,vp0,qs,ms,rho_back,source)
		type(PM1D), intent(inout) :: this
		real(mp), dimension(this%n), intent(in) :: xp0, vp0
		real(mp), intent(in) :: qs(this%n), ms(this%n), rho_back
		integer :: k
		interface
			subroutine source(pm,k,str)
				use modPM1D
				type(PM1D), intent(inout) :: pm
				integer, intent(in) :: k
				character(len=*), intent(in) :: str
			end subroutine
		end interface

		!initial condition
		call setPlasma(this%p,xp0,vp0,qs,ms)
		call setMesh(this%m,rho_back)

		!Time stepping
		call halfStep(this)
		do k=1,this%nt
			call updatePlasma(this,source,k)
		end do
	end subroutine

	subroutine halfStep(this)
		type(PM1D), intent(inout) :: this
		integer :: i, j
		real(mp) :: rhs(this%ng-1)
		real(mp) :: phi1(this%ng-1)
		real(mp) :: dt, L
		integer :: N,Ng
		dt = this%dt
		L = this%L
		N = this%N
		Ng = this%ng

		call assignMatrix(this%a,this%m,this%p%xp)

		!charge assignment
		call chargeAssign(this%a,this%p,this%m)

		!field equation rhs
		rhs = -this%m%rho(1:Ng-1)/this%eps0

		!solve field equation
		call CG_K(phi1,rhs,this%m%dx)				!Conjugate-Gradient Iteration
		this%m%phi(1:Ng-1) = phi1
		this%m%phi(Ng) = 0.0_mp

		!Electric field : -D*phi
		this%m%E = - multiplyD(this%m%phi,this%m%dx)

		!Force assignment : mat'*E
		call forceAssign(this%a,this%p,this%m)

		!Half time step advancement in velocity
		this%p%vp = this%p%vp + dt/2.0_mp*this%p%qs/this%p%ms*this%p%Ep
	end subroutine

	subroutine updatePlasma(this,source,k)
		type(PM1D), intent(inout) :: this
		integer, intent(in) :: k
		real(mp) :: rhs(this%ng-1), phi1(this%ng-1)
		real(mp) :: dt, L, B
		integer :: N, Ng
		interface
			subroutine source(pm,k,str)
				use modPM1D
				type(PM1D), intent(inout) :: pm
				integer, intent(in) :: k
				character(len=*), intent(in) :: str
			end subroutine
		end interface
		dt = this%dt
		L = this%L
		N = this%n
		B = this%B0
		Ng = this%ng

		call recordPlasma(this%r, this%p, this%m ,k)									!record for n=0~(Nt-1)

		call source(this,k,'xp')
		this%p%xp = this%p%xp + dt*this%p%vp

		call assignMatrix(this%a,this%m,this%p%xp)

		!charge assignment
		call chargeAssign(this%a,this%p,this%m)

		!field equation rhs
		rhs = -this%m%rho(1:Ng-1)/this%eps0
		!solve field equation
		call CG_K(phi1,rhs,this%m%dx)				!Conjugate-Gradient Iteration
		this%m%phi(1:Ng-1) = phi1
		this%m%phi(Ng) = 0.0_mp

		!Electric field : -D*phi
		this%m%E = - multiplyD(this%m%phi,this%m%dx)

		!Force assignment : mat'*E
		call forceAssign(this%a, this%p, this%m)

		!Half time step advancement in velocity
		this%p%vp = this%p%vp + dt*this%p%qs/this%p%ms*this%p%Ep
	end subroutine
!
!	subroutine QOI(this,J)
!		type(plasma), intent(in) :: this
!		real(mp), intent(out) :: J
!		J = 1.0_mp/Ng/(Nt-Ni)*SUM(this%Edata(:,Ni+1:Nt)**2)
!	end subroutine
!
!	subroutine QOItwo(this,A,B,J)
!	type(plasma), intent(in) :: this
!	integer, intent(in) :: A, B
!	real(mp), intent(out) :: J
!	J = 1.0_mp/Ng/(B-A)*SUM(this%Edata(:,A+1:B)**2)		!omitted 1/N/T for the sake of machine precision
!	end subroutine
!

!
!	subroutine adjoint(this,delJdelA)
!		type(plasma), intent(inout) :: this
!		real(mp), intent(out) :: delJdelA
!		real(mp) :: dts
!		real(mp) :: xps(size(this%xp)), vps(size(this%vp)), Es(size(this%E)), phis(size(this%phi))
!		real(mp) :: dEs(size(this%E))
!		real(mp) :: rhs(size(this%E)), phis1(size(this%phi)-1)
!		integer :: i,k, nk
!		real(mp) :: coeff(this%n)
!		dts = -dt
!		xps = 0.0_mp
!		vps = 0.0_mp
!
!		do k=1,this%nt
!			!vps update
!			nk = this%nt+1-k
!!			dvps = merge(2.0_mp/N/(Nt-Ni)/dt*this%vpdata(:,this%nt+1-k),0.0_mp,nk>=Ni+1)
!			if (k==this%nt) then
!!				vps = 1.0_mp/2.0_mp*( vps + dts*( -xps + dvps ) )
!				vps = 1.0_mp/2.0_mp*( vps + dts*( -xps ) )
!			else
!!				vps = vps + dts*( -xps + dvps )
!				vps = vps + dts*( -xps )
!			end if
!			!assignment
!			call assignMatrix(this,this%xpdata(:,nk))
!
!			!Es update : (qe/me/dx)*mat*vps
!			call vpsAssign(this,Es,vps)
!			dEs = merge( -2.0_mp/Ng/(Nt-Ni)/dt/dx*this%Edata(:,this%nt+1-k),0.0_mp,nk>=Ni+1)
!			Es = Es + dEs
!
!			!phis update
!			!rhs = D*Es
!			rhs = multiplyD(Es,dx)
!
!			!K*phis = D*Es
!			call CG_K(phis1,rhs(1:Ng-1),dx)				!Conjugate-Gradient Iteration
!			phis(1:Ng-1) = phis1
!			phis(Ng) = 0.0_mp
!
!			!xps update
!			if(nk==this%ni) then
!				xps = xps - dts*xps*(B*L/N*2.0_mp*pi*mode/L)*COS(2.0_mp*pi*mode*this%xpdata(:,nk)/L)
!			end if
!			xps = xps + xpsUpdate(this,vps,phis,dts,k)
!
!			call recordAdjoint(this,nk,xps)					!recorde xps from n=(Nt-1) to n=0
!			if( nk<this%ni ) then
!				exit
!			end if
!		end do
!
!		coeff = - L/N*SIN( 2.0_mp*pi*mode*this%xpdata(:,this%ni)/L )
!		delJdelA = 0.0_mp
!		do i = 1,N
!			delJdelA = delJdelA + dt*coeff(i)*this%xpsdata(i,this%ni+1)
!		end do
!	end subroutine

end module