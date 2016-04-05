module modSource

	use modPM1D

	implicit none

contains
!==============Wave perturbation on position===============================

	subroutine IC_wave(this,k,str)
		type(PM1D), intent(inout) :: this
		integer, intent(in) :: k
		character(len=*), intent(in) :: str

		SELECT CASE (str)
			CASE('xp')
				if(k==this%ni) then
					this%p%xp = this%p%xp + this%dt*this%B0*this%L/this%n*SIN(4.0_mp*pi*this%p%xp/this%L)		!xp_(Ni-1) : k=Ni, xp_Ni : k=(Ni+1)
				end if
		END SELECT
	end subroutine

end module