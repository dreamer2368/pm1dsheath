module modPlasma

	use constants

	implicit none

	type plasma
		integer :: n
		real(mp), allocatable :: xp(:)
		real(mp), allocatable :: vp(:)
		real(mp), allocatable :: Ep(:)

		real(mp), allocatable :: ms(:), qs(:)
	end type

contains

	subroutine buildPlasma(this,n)
		type(plasma), intent(out) :: this
		integer, intent(in) :: n

		this%n = n

		allocate(this%xp(n))
		allocate(this%vp(n))
		allocate(this%Ep(n))

		allocate(this%qs(n))
		allocate(this%ms(n))
	end subroutine

	subroutine setPlasma(this,xp0,vp0,qs,ms)
		type(plasma), intent(inout) :: this
		real(mp), intent(in) :: xp0(this%n)
		real(mp), intent(in) :: vp0(this%n)
		real(mp), intent(in) :: qs(this%n), ms(this%n)

		this%xp = xp0
		this%vp = vp0

		this%qs = qs
		this%ms = ms
	end subroutine

	subroutine destroyPlasma(this)
		type(plasma), intent(inout) :: this

		deallocate(this%xp)
		deallocate(this%vp)
		deallocate(this%Ep)
		deallocate(this%qs)
		deallocate(this%ms)
	end subroutine

end module