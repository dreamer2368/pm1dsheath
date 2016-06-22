module timeStep

	use modSource
	use modTarget
	use modBC
	use modRecord
	use ArMCC

	implicit none

contains

	subroutine forwardsweep(this,r,target_input,source)
		type(PM1D), intent(inout) :: this
		type(recordData), intent(inout) :: r
		integer :: k
		interface
			subroutine target_input(pm,k,str)
				use modPM1D
				type(PM1D), intent(inout) :: pm
				integer, intent(in) :: k
				character(len=*), intent(in) :: str
			end subroutine
		end interface
		interface
			subroutine source(pm)
				use modPM1D
				type(PM1D), intent(inout) :: pm
			end subroutine
		end interface

		!Time stepping
		call halfStep(this)
		do k=1,this%nt
			call updatePlasma(this,r,target_input,source,k)
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

		call applyBC(this)
		do i=1,this%n
			call assignMatrix(this%a(i),this%m,this%p(i)%xp)
		end do
		call adjustGrid(this)

		!charge assignment
		call chargeAssign(this%a,this%p,this%m)

		call solveMesh(this%m,this%eps0)

		!Electric field : -D*phi
		this%m%E = - multiplyD(this%m%phi,this%m%dx,this%m%BCindex)

		!Force assignment : mat'*E
		do i=1, this%n
			call forceAssign(this%a(i),this%p(i),this%m)
		end do

		!Half time step advancement in velocity
		do i=1, this%n
			call accelSpecies(this%p(i),dt/2.0_mp)
		end do
	end subroutine

	subroutine updatePlasma(this,r,target_input,source,k)
		type(PM1D), intent(inout) :: this
		type(recordData), intent(inout) :: r
		integer, intent(in) :: k
		real(mp) :: rhs(this%ng-1), phi1(this%ng-1)
		real(mp) :: dt, L
		integer :: N, Ng, i
		interface
			subroutine target_input(pm,k,str)
				use modPM1D
				type(PM1D), intent(inout) :: pm
				integer, intent(in) :: k
				character(len=*), intent(in) :: str
			end subroutine
		end interface
		interface
			subroutine source(pm)
				use modPM1D
				type(PM1D), intent(inout) :: pm
			end subroutine
		end interface
		dt = this%dt
		L = this%L
		N = this%n
		Ng = this%ng
		call recordPlasma(r, this, k)									!record for n=0~(Nt-1)

		call target_input(this,k,'xp')

		call source(this)

		do i=1,this%n
			call moveSpecies(this%p(i),dt)
		end do

		call applyBC(this)
		do i=1, this%n
			call assignMatrix(this%a(i),this%m,this%p(i)%xp)
		end do
		call adjustGrid(this)

		!charge assignment
		call chargeAssign(this%a,this%p,this%m)

		call solveMesh(this%m,this%eps0)

		!Electric field : -D*phi
		this%m%E = - multiplyD(this%m%phi,this%m%dx,this%m%BCindex)

		!Force assignment : mat'*E
		do i=1, this%n
			call forceAssign(this%a(i), this%p(i), this%m)
		end do

		do i=1, this%n
			call accelSpecies(this%p(i),dt)
		end do
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