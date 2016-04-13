module modBC

	use modPM1D

	implicit none

contains

	subroutine applyBC(pm)
		type(PM1D), intent(inout) :: pm
		integer :: i

		select case(pm%BCindex)
			case(0)	!periodic
				do i=1,pm%n
					call applyBC_periodic(pm%p(i),pm%m)
				end do
		end select
	end subroutine

	subroutine adjustGrid(pm)
		type(PM1D), intent(inout) :: pm
		integer :: i

		select case(pm%BCindex)
			case(0)	!periodic
				do i=1,pm%n
					call adjustGrid_periodic(pm%a(i))
				end do
		end select
	end subroutine

	subroutine applyBC_periodic(p,m)
		type(species), intent(inout) :: p
		type(mesh), intent(inout) :: m
		integer :: i

		!apply BC
		do i=1,p%np
			if( p%xp(i)<0 ) then
				p%xp(i) = p%xp(i) + m%L
			elseif( p%xp(i)>=m%L ) then
				p%xp(i) = p%xp(i) - m%L
			end if
		end do
	end subroutine

	subroutine adjustGrid_periodic(a)
		type(pmassign), intent(inout) :: a
		integer :: i

		do i=1,a%np
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

		if( MINVAL(a%g(:,1))<1 .or. MAXVAL(a%g(:,2))>a%ng ) then
			print *, MINVAL(a%g(:,1)), MAXVAL(a%g(:,2))
			print *, 'Boundary handling is failed. particle is way outside BC. stopped time stepping.'
			stop
		end if
	end subroutine

end module