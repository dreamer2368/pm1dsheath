module modSpecies

	use constants

	implicit none

	type species
		integer :: np
		real(mp), allocatable :: xp(:)
		real(mp), allocatable :: vp(:)
		real(mp), allocatable :: Ep(:)

		real(mp) :: ms, qs
	end type

contains

	subroutine buildSpecies(this,qs,ms)
		type(species), intent(out) :: this
		real(mp), intent(in) :: ms, qs

		this%ms = ms
		this%qs = qs
	end subroutine

	subroutine setSpecies(this,np0,xp0,vp0)
		type(species), intent(inout) :: this
		integer, intent(in) :: np0
		real(mp), intent(in) :: xp0(np0)
		real(mp), intent(in) :: vp0(np0)

		this%np = np0
		allocate(this%xp(np0))
		allocate(this%vp(np0))
		allocate(this%Ep(np0))

		this%xp = xp0
		this%vp = vp0
	end subroutine

	subroutine destroySpecies(this)
		type(species), intent(inout) :: this

		deallocate(this%xp)
		deallocate(this%vp)
		deallocate(this%Ep)
	end subroutine

	subroutine moveSpecies(this,dt)
		type(species), intent(inout) :: this
		real(mp), intent(in) :: dt

		this%xp = this%xp + dt*this%vp
	end subroutine

	subroutine accelSpecies(this,dt)
		type(species), intent(inout) :: this
		real(mp), intent(in) :: dt

		this%vp = this%vp + dt*this%qs/this%ms*this%Ep
	end subroutine

end module