MODULE random

	use constants

	implicit none

contains

	function randn(N) result(x)
		integer, intent(in) :: N
		real(mp) :: r1(N), r2(N), x(N)

		call RANDOM_NUMBER(r1)
		call RANDOM_NUMBER(r2)
		x = sqrt( -2.0_mp*log(r1) )*cos( 2.0_mp*pi*r2 )
	end function

	function randr(N) result(x)
		integer, intent(in) :: N
		real(mp) :: r1(N/2), r2(N-N/2), x(N)

		call RANDOM_NUMBER(r1)
		call RANDOM_NUMBER(r2)
		x(1:N/2) = sqrt( -2.0_mp*log(r1) )
		x(N/2+1:N) = -sqrt( -2.0_mp*log(r2) )
	end function

END MODULE random