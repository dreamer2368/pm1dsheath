module modTarget

	use modPM1D

	implicit none

contains
!==============Default=====================================================

	subroutine Null_input(this,k,str)
		type(PM1D), intent(inout) :: this
		integer, intent(in) :: k
		character(len=*), intent(in) :: str
	end subroutine

!==============Wave perturbation on position===============================

	subroutine IC_wave(this,k,str)
		type(PM1D), intent(inout) :: this
		integer, intent(in) :: k
		character(len=*), intent(in) :: str

		SELECT CASE (str)
			CASE('xp')
				if(k==this%ni) then
					this%p(1)%xp = this%p(1)%xp + this%dt*this%A0(2)*this%L/this%p(1)%np*SIN(4.0_mp*pi*this%p(1)%xp/this%L)		!xp_(Ni-1) : k=Ni, xp_Ni : k=(Ni+1)
				end if
		END SELECT
	end subroutine

end module