module modAssign

	use modPlasma
	use modMesh

	implicit none

	type pmAssign
		integer :: n
		integer :: ng

		integer, allocatable :: g(:,:)
		real(mp), allocatable :: frac(:,:)
		real(mp), allocatable :: h(:)
	end type

contains

	subroutine buildAssign(this,n,ng)
		type(pmAssign), intent(out) :: this
		integer, intent(in) :: n
		integer, intent(in) :: ng

		this%n = n
		this%ng = ng

		allocate(this%g(n,3))
		allocate(this%frac(n,3))
		allocate(this%h(n))
	end subroutine

	subroutine destroyAssign(this)
		type(pmAssign), intent(inout) :: this

		deallocate(this%g)
		deallocate(this%frac)
		deallocate(this%h)
	end subroutine

	subroutine assignMatrix(this,m,xp)	!apply BC and create assignment matrix
		type(pmAssign), intent(inout) :: this
		type(mesh), intent(inout) :: m
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
				xp(i) = xp(i) + m%L
			elseif( xp(i)>=m%L ) then
				xp(i) = xp(i) - m%L
			end if
		end do

		if( MINVAL(xp)<0 .or. MAXVAL(xp)>m%L ) then
			print *, MINVAL(xp), MAXVAL(xp)
			print *, 'Boundary handling is failed. particle is way outside BC. stopped time stepping.'
			stop
		end if

		!assignment matrix
		do i=1,this%n
			g = FLOOR(xp(i)/m%dx + 0.5_mp)+1
			gl = g-1
			gr = g+1
			h = xp(i)/m%dx - g + 1.0_mp

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

			this%g(i,:) = (/gl,g,gr/)
			this%h(i) = h
			this%frac(i,:) = (/ fracl, frac0, fracr /)
		end do
	end subroutine

	subroutine chargeAssign(this,p,m)
		type(pmAssign), intent(inout) :: this
		type(plasma), intent(inout) :: p
		type(mesh), intent(inout) :: m
		integer :: i

		m%rho = 0.0_mp
		do i=1,this%n
			m%rho( this%g(i,:) ) = m%rho( this%g(i,:) ) + p%qs/m%dx*this%frac(i,:)
		end do
		m%rho = m%rho + m%rho_back
	end subroutine

	subroutine forceAssign(this,p,m)
		type(pmAssign), intent(inout) :: this
		type(plasma), intent(inout) :: p
		type(mesh), intent(inout) :: m
		integer :: i

		p%Ep = 0.0_mp
		do i=1,this%n
			p%Ep(i) = sum( this%frac(i,:)*m%E(this%g(i,:)) )
		end do
	end subroutine

end module