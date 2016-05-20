MODULE random

	use constants

	implicit none

contains

	function randn(N) result(x)
		integer, intent(in) :: N
		integer :: i,nseed,clock
		integer, allocatable :: seed(:)
		real(mp) :: r1(N), r2(N), x(N)

		call RANDOM_SEED(size=nseed)
		allocate(seed(nseed))
		call SYSTEM_CLOCK(COUNT=clock)
		seed = clock + 37*(/ ( i-1, i=1,nseed ) /)
		call RANDOM_SEED(put=seed)
		deallocate(seed)

		call RANDOM_NUMBER(r1)
		call RANDOM_NUMBER(r2)
		x = sqrt( -2.0_mp*log(r1) )*cos( 2.0_mp*pi*r2 )
	end function

END MODULE random