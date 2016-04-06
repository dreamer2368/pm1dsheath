module modBC

	use modPM1D

	implicit none

contains

	subroutine applyBC(pm)
		type(PM1D), intent(inout) :: pm

		select case(pm%BCindex)
			case(0)	!periodic
				call applyBC_periodic(pm%a,pm%p,pm%m)
		end select
	end subroutine

	subroutine adjustGrid(pm)
		type(PM1D), intent(inout) :: pm

		select case(pm%BCindex)
			case(0)	!periodic
				call adjustGrid_periodic(pm%a)
		end select
	end subroutine

	subroutine applyBC_periodic(a,p,m)
		type(pmassign), intent(in) :: a
		type(plasma), intent(inout) :: p
		type(mesh), intent(inout) :: m
		integer :: i

		!apply BC
		do i=1,p%n
			if( p%xp(i)<0 ) then
				p%xp(i) = p%xp(i) + m%L
			elseif( p%xp(i)>=m%L ) then
				p%xp(i) = p%xp(i) - m%L
			end if
		end do

		if( MINVAL(p%xp)<0 .or. MAXVAL(p%xp)>m%L ) then
			print *, MINVAL(p%xp), MAXVAL(p%xp)
			print *, 'Boundary handling is failed. particle is way outside BC. stopped time stepping.'
			stop
		end if
	end subroutine

	subroutine adjustGrid_periodic(a)
		type(pmassign), intent(inout) :: a
		integer :: i

		do i=1,a%n
			if( a%g(i,1)<1 ) then
				a%g(i,1) = a%g(i,1) + a%ng
			elseif( a%g(i,1)>a%ng ) then
				a%g(i,1) = a%g(i,1) - a%ng
			end if
			if( a%g(i,2)<1 ) then
				a%g(i,2) = a%g(i,2) + a%ng
			elseif( a%g(i,2)>a%ng ) then
				a%g(i,2) = a%g(i,2) - a%ng
			end if
		end do
	end subroutine

end module