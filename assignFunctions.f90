module assignFunctions

	use declaration

	implicit none

contains

	subroutine assignMatrix(this,xp)	!apply BC and create assignment matrix
		type(plasma), intent(inout) :: this
		real(mp), intent(inout) :: xp(this%n)
		integer :: i
		integer :: g		!nearest grid point
		integer :: gl		!left grid point
		integer :: gr		!right grid point
		real(mp) :: fracl		!fraction for left grid point
		real(mp) :: frac0		!fraction for nearest grid point
		real(mp) :: fracr		!fraction for right grid point
		real(mp) :: h			!distance from nearest grid point

		!apply BC
		do i=1,this%n
			if( xp(i)<0 ) then
				xp(i) = xp(i) + L
			elseif( xp(i)>=L ) then
				xp(i) = xp(i) - L
			end if
		end do

		!assignment matrix
		do i=1,this%n
			g = FLOOR(xp(i)/dx + 0.5_mp)+1
			gl = g-1
			gr = g+1
			h = xp(i)/dx - g + 1.0_mp

			if (g<1) then
				g = g + this%ng
			elseif (g>this%ng) then
				g = g - this%ng
			end if
			if (gl<1) then
				gl = gl + this%ng
			elseif (gl>this%ng) then
				gl = gl - this%ng
			end if
			if (gr<1) then
				gr = gr + this%ng
			elseif (gr>this%ng) then
				gr = gr - this%ng
			end if

			frac0 = 0.75_mp - h*h
			fracl = 0.5_mp*(0.5_mp-h)*(0.5_mp-h)
			fracr = 0.5_mp*(0.5_mp+h)*(0.5_mp+h)

			this%g(i) = g
			this%gl(i) = gl
			this%gr(i) = gr
			this%h(i) = h
			this%f0(i) = frac0
			this%fl(i) = fracl
			this%fr(i) = fracr
		end do
	end subroutine

	subroutine chargeAssign(this)
		type(plasma), intent(inout) :: this
		integer :: i

		this%rho = 0.0_mp
		do i=1,this%n
			this%rho( this%gl(i) ) = this%rho( this%gl(i) ) + qe/dx*this%fl(i)
			this%rho( this%g(i) ) = this%rho( this%g(i) ) + qe/dx*this%f0(i)
			this%rho( this%gr(i) ) = this%rho( this%gr(i) ) + qe/dx*this%fr(i)
		end do
		this%rho = this%rho + rho_back
	end subroutine

	subroutine forceAssign(this)
		type(plasma), intent(inout) :: this
		integer :: i

		this%Ep = 0.0_mp
		do i=1,this%n
			this%Ep(i) = this%fl(i)*this%E( this%gl(i) ) + this%fr(i)*this%E( this%gr(i) ) + this%f0(i)*this%E( this%g(i) )
		end do
	end subroutine

	subroutine vpsAssign(this,Es,vps)
		type(plasma), intent(in) :: this
		real(mp), intent(in) :: vps(:)
		real(mp), intent(out) :: Es( size(this%E) )
		integer :: i

		Es = 0.0_mp
		do i=1,this%n
			Es( this%gl(i) ) = Es( this%gl(i) ) + this%fl(i)*qe/me/dx*vps(i)
			Es( this%g(i) ) = Es( this%g(i) ) + this%f0(i)*qe/me/dx*vps(i)
			Es( this%gr(i) ) = Es( this%gr(i) ) + this%fr(i)*qe/me/dx*vps(i)
		end do
	end subroutine

	function xpsUpdate(this,vps,phis,dts,k) result(dxps)
		type(plasma), intent(in) :: this
		real(mp), intent(in) :: vps(:)
		real(mp), intent(in) :: phis(:)
		real(mp), intent(in) :: dts
		integer, intent(in) :: k					!current time step
		real(mp) :: dxps( size(vps) )
		integer :: i

		dxps = 0.0_mp
		do i=1,this%n
!			dxps(i) = dts*( -qe/me/dx*vps(i)*( this%Edata(this%gr(i),this%nt+1-k) - this%Edata(this%gl(i),this%nt+1-k) ) )
!			dxps(i) = dxps(i) + dts*qe/dx/eps0*( phis(this%gr(i)) - phis(this%gl(i)) )
			dxps(i) = dts*(-qe)/me/dx*vps(i)*this%Edata(this%gl(i),this%nt+1-k)*( -0.5_mp+this%h(i) )
			dxps(i) = dxps(i) + dts*(-qe)/me/dx*vps(i)*this%Edata(this%g(i),this%nt+1-k)*( -2.0_mp*this%h(i) )
			dxps(i) = dxps(i) + dts*(-qe)/me/dx*vps(i)*this%Edata(this%gr(i),this%nt+1-k)*( +0.5_mp+this%h(i) )
			dxps(i) = dxps(i) + dts*qe/dx/eps0*phis(this%gl(i))*( -0.5_mp+this%h(i) )
			dxps(i) = dxps(i) + dts*qe/dx/eps0*phis(this%g(i))*( -2.0_mp*this%h(i) )
			dxps(i) = dxps(i) + dts*qe/dx/eps0*phis(this%gr(i))*( +0.5_mp+this%h(i) )
		end do
	end function

end module