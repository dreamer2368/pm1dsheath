module timeStep

	use declaration
	use assignFunctions
	use MatrixVector

	implicit none

contains

	subroutine halfStep(this)
		type(plasma), intent(inout) :: this
		integer :: i, j
		real(mp) :: rhs(Ng-1)
		real(mp) :: phi1(Ng-1)

		call assignMatrix(this,this%xp)

		!charge assignment
		call chargeAssign(this)

		!field equation rhs
		rhs = -this%rho(1:Ng-1)/eps0
		!solve field equation
		call CG_K(phi1,rhs,dx)				!Conjugate-Gradient Iteration
		this%phi(1:Ng-1) = phi1
		this%phi(Ng) = 0.0_mp

		!Electric field : -D*phi
		this%E = - multiplyD(this%phi,dx)

		!Force assignment : mat'*E
		call forceAssign(this)

		!Half time step advancement in velocity
		this%vp = this%vp + dt/2.0_mp*qe/me*this%Ep
	end subroutine

	subroutine updatePlasma(this)
		type(plasma), intent(inout) :: this
		integer :: i, j, k
		real(mp) :: rhs(Ng-1)
		real(mp) :: phi1(Ng-1)
		this%time = 0.0_mp
		do k = 1,this%nt
			call recordPlasma(this,k)									!record for n=0~(Nt-1)
			if( k<this%nt ) then
				this%time(k+1) = this%time(k) + dt
			end if
			if(k==this%ni) then
				this%xp = this%xp + dt*B*L/N*SIN(4.0_mp*pi*this%xp/L*mode)		!xp_(Ni-1) : k=Ni, xp_Ni : k=(Ni+1)
			end if
			this%xp = this%xp + dt*this%vp
			call assignMatrix(this,this%xp)

			!Sometimes time stepping continues even though particles are assigned way outside the boundary. Even though it continues without error message, it's wrong.
			if( MAXVAL(this%gr)>this%ng .or. MAXVAL(this%gl)>this%ng ) then
				print *, 'n= ', k, MAXVAL(this%gr), MAXVAL(this%gl)
				print *, 'Boundary handling is failed. particle is way outside BC. stopped time stepping.'
				stop
			end if
			if( MINVAL(this%gl)<1 .or. MINVAL(this%gr)<1 ) then
				print *, 'n= ', k, MINVAL(this%gr), MINVAL(this%gl)
				print *, 'Boundary handling is failed. particle is way outside BC. stopped time stepping.'
				stop
			end if

			!charge assignment
			call chargeAssign(this)
			!field equation rhs
			rhs = -this%rho(1:Ng-1)/eps0
			!solve field equation
			call CG_K(phi1,rhs,dx)				!Conjugate-Gradient Iteration
			this%phi(1:Ng-1) = phi1
			this%phi(Ng) = 0.0_mp
			!Electric field : -D*phi
			this%E = - multiplyD(this%phi,dx)
			!Force assignment : mat'*E
			call forceAssign(this)
			!Half time step advancement in velocity
			this%vp = this%vp + dt*qe/me*this%Ep
		end do
	end subroutine

	subroutine QOI(this,J)
		type(plasma), intent(in) :: this
		real(mp), intent(out) :: J
		J = 1.0_mp/N/(Nt-Ni)*SUM(this%vpdata(:,Ni+1:Nt)**2)		!omitted 1/N/T for the sake of machine precision
	end subroutine

	subroutine QOItwo(this,A,B,J)
	type(plasma), intent(in) :: this
	integer, intent(in) :: A, B
	real(mp), intent(out) :: J
	J = 1.0_mp/N/(B-A)*SUM(this%vpdata(:,A+1:B)**2)		!omitted 1/N/T for the sake of machine precision
	end subroutine

	subroutine forwardsweep(this)
		type(plasma), intent(inout) :: this

		call buildPlasma(this,Nt,Ni,N,Ng,xp0,vp0)
		!Time stepping
		call halfStep(this)
		call updatePlasma(this)
		!export the data - when needed
		call printPlasma(this,L,dt)
	end subroutine

	subroutine adjoint(this,delJdelA)
		type(plasma), intent(inout) :: this
		real(mp), intent(out) :: delJdelA
		real(mp) :: dts
		real(mp) :: xps(size(this%xp)), vps(size(this%vp)), Es(size(this%E)), phis(size(this%phi))
		real(mp) :: dvps(size(this%vp))
		real(mp) :: rhs(size(this%E)), phis1(size(this%phi)-1)
		integer :: i,k, nk
		real(mp) :: coeff(this%n)
		dts = -dt
		xps = 0.0_mp
		vps = 0.0_mp
		do k=1,this%nt
			!vps update
			nk = this%nt+1-k
			dvps = merge(2.0_mp/N/(Nt-Ni)/dt*this%vpdata(:,this%nt+1-k),0.0_mp,nk>=Ni+1)
			if (k==this%nt) then
				vps = 1.0_mp/2.0_mp*( vps + dts*( -xps + dvps ) )	!omitted 1/N/T(see QOI function)
			else
				vps = vps + dts*( -xps + dvps )
			end if
			!assignment
			call assignMatrix(this,this%xpdata(:,nk))

			!Es update : (qe/me/dx)*mat*vps
			call vpsAssign(this,Es,vps)

			!phis update
			!rhs = D*Es
			rhs = multiplyD(Es,dx)

			!K*phis = D*Es
			call CG_K(phis1,rhs(1:Ng-1),dx)				!Conjugate-Gradient Iteration
			phis(1:Ng-1) = phis1
			phis(Ng) = 0.0_mp

			!xps update
			if(nk==this%ni) then
				xps = xps - dts*xps*(B*L/N*4.0_mp*pi*mode/L)*COS(4.0_mp*pi*mode*this%xpdata(:,nk)/L)
			end if
			xps = xps + xpsUpdate(this,vps,phis,dts,k)

			call recordAdjoint(this,nk,xps)					!recorde xps from n=(Nt-1) to n=0
			if( nk<this%ni ) then
				exit
			end if
		end do
		coeff = - L/N*SIN( 4.0_mp*pi*mode*this%xpdata(:,this%ni)/L )
		delJdelA = 0.0_mp
		do i = 1,N
			delJdelA = delJdelA + dt*coeff(i)*this%xpsdata(i,this%ni+1)
		end do
	end subroutine

end module