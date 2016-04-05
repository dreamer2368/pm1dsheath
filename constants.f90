module constants

	implicit none

	integer, parameter :: mp = SELECTED_REAL_KIND(15,307)
	real(mp), parameter :: pi = 4.0_mp*ATAN(1.0_mp)
	complex(mp), parameter :: eye = (0.0_mp,1.0_mp)

contains

end module