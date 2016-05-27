module modBC

	use modPM1D
	use random

	implicit none

contains

!==============================particle BC=====================================

	subroutine applyBC(pm)
		type(PM1D), intent(inout) :: pm
		integer :: i

		select case(pm%pBCindex)
			case(0)	!periodic
				do i=1,pm%n
					call applyBC_periodic(pm%p(i),pm%m)
				end do
			case(1) !absorbing-absorbing
				do i=1,pm%n
					call applyBC_absorbing(pm%p(i),pm%m)
				end do
			case(2) !refluxing-absorbing
				do i=1,pm%n
					call applyBC_refluxing_absorbing(pm%p(i),pm%m,pm%dt,pm%A0(i))
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

	subroutine applyBC_absorbing(p,m)
		type(species), intent(inout) :: p
		type(mesh), intent(inout) :: m
		integer :: i, np1
		real(mp), allocatable :: vec(:)

		np1 = p%np
		i = 1
		!apply BC
		do while( i .le. np1 )
			if( p%xp(i).le.0.0_mp ) then
				p%xp(i) = p%xp(np1)					!replacement with the last particle
				p%vp(i) = p%vp(np1)
				p%Ep(i) = p%Ep(np1)
				m%rho_back(1) = m%rho_back(1) + p%spwt*p%qs
				np1 = np1-1
				i = i-1
			elseif( p%xp(i).ge.m%L ) then
				p%xp(i) = p%xp(np1)
				p%vp(i) = p%vp(np1)
				p%Ep(i) = p%Ep(np1)
				m%rho_back(m%ng) = m%rho_back(m%ng) + p%spwt*p%qs
				np1 = np1-1
				i = i-1
			end if
			i = i+1
		end do
		p%np = np1

		allocate(vec(np1))

		vec = p%xp(1:np1)
		deallocate(p%xp)
		allocate(p%xp(np1))
		p%xp = vec

		vec = p%vp(1:np1)
		deallocate(p%vp)
		allocate(p%vp(np1))
		p%vp = vec

		vec = p%Ep(1:np1)
		deallocate(p%Ep)
		allocate(p%Ep(np1))
		p%Ep = vec

		deallocate(vec)
	end subroutine

	subroutine applyBC_refluxing_absorbing(p,m,dt,vT)			!refluxing at the left plane, absorbing at the right plane
		type(species), intent(inout) :: p
		type(mesh), intent(inout) :: m
		real(mp), intent(in) :: dt, vT
		real(mp) :: temp(1)
		integer :: i, np1
		real(mp), allocatable :: vec(:)

		!apply refluxing BC
		do i=1,p%np
			if( p%xp(i).le.0.0_mp ) then
				temp = abs(vT*randn(1))
				p%vp(i) = temp(1)
				call RANDOM_NUMBER(temp)
				p%xp(i) = temp(1)*dt*p%vp(i)
			end if
		end do

		np1 = p%np
		i = 1
		!apply absorbing BC
		do while( i .le. np1 )
			if( p%xp(i).ge.m%L ) then
				p%xp(i) = p%xp(np1)
				p%vp(i) = p%vp(np1)
				p%Ep(i) = p%Ep(np1)
				m%rho_back(m%ng) = m%rho_back(m%ng) + p%spwt*p%qs
				np1 = np1-1
				i = i-1
			end if
			i = i+1
		end do
		p%np = np1

		allocate(vec(np1))

		vec = p%xp(1:np1)
		deallocate(p%xp)
		allocate(p%xp(np1))
		p%xp = vec

		vec = p%vp(1:np1)
		deallocate(p%vp)
		allocate(p%vp(np1))
		p%vp = vec

		vec = p%Ep(1:np1)
		deallocate(p%Ep)
		allocate(p%Ep(np1))
		p%Ep = vec

		deallocate(vec)
	end subroutine

!===========================grid adjustment for BC=============================

	subroutine adjustGrid(pm)
		type(PM1D), intent(inout) :: pm
		integer :: i

		select case(pm%mBCindex)
			case(0)	!periodic
				do i=1,pm%n
					call adjustGrid_periodic(pm%a(i))
				end do
			case(1)	!Dirichlet-Dirichlet
				do i=1,pm%n
					call adjustGrid_absorbing(pm%a(i),pm%p(i))
				end do
			case(2)	!Dirichlet-Neumann
				do i=1,pm%n
					call adjustGrid_absorbing(pm%a(i),pm%p(i))
				end do
		end select
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

	subroutine adjustGrid_absorbing(a,p)
		type(pmassign), intent(inout) :: a
		type(species), intent(in) :: p
		integer :: i

		if( MINVAL(a%g(:,1))<1 .or. MAXVAL(a%g(:,2))>a%ng ) then
			print *, MINVAL(a%g(:,1)), p%xp(minloc(a%g(:,1))), MAXVAL(a%g(:,2)), p%xp(maxloc(a%g(:,2)))
			print *, 'Boundary handling is failed. particle is way outside BC. stopped time stepping.'
			stop
		end if

		!adjustment for boundary : charge/(dx/2)
		do i=1,a%np
			if( a%g(i,1).eq.1 ) then
				a%frac(i,1) = a%frac(i,1)*2.0_mp
			end if
			if( a%g(i,2).eq.a%ng ) then
				a%frac(i,2) = a%frac(i,2)*2.0_mp
			end if
		end do
	end subroutine

end module