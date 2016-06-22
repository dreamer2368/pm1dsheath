module ArMCC

	use modPM1D
	use random
	implicit none

contains

!=======================================================
!  e + Ar -> e + Ar  Elastic      (with Ramseur minimum)
!=======================================================
	function asigma1(energy) result(sig1)
		real(mp), intent(in) :: energy
		real(mp) :: sig1

		if( energy < 0.2_mp ) then
			sig1 = 1.0_mp/( 10.0_mp**(19.0_mp+energy/0.11_mp) )
		else
			sig1 = 9.070000E-19*(energy**1.55_mp)*( (energy+70.0_mp)**1.10_mp )/( (14.0_mp+energy)**3.25_mp )
		end if
	end function

!=======================================================
!  e + Ar -> e + Ar  Excitation
!=======================================================
	function asigma2(energy) result(sig2)
		real(mp), intent(in) :: energy
		real(mp) :: sig2

		if( energy > 12.0_mp ) then
			sig2 = ( 3.85116E-19*log(energy/3.4015_mp) - 4.85227E-19 )/energy
		else
			sig2 = 0.0_mp
		end if
	end function

!=======================================================
!  e + Ar -> e + e + Ar+  Ion.
!=======================================================
	function asigma3(energy) result(sig3)
		real(mp), intent(in) :: energy
		real(mp) :: sig3

		if( energy > 15.76_mp ) then
			sig3 = 1.3596E-18/energy*log( (energy+120.0_mp/energy)/15.76_mp )*	&
					( atan( (energy**2 - 9.76_mp*energy + 2.4_mp)/(20.6_mp*energy + 206.0_mp) ) +	&
						atan( (2.0_mp*energy - 80.0_mp)/(10.3_mp*energy+103.0_mp) ) )
		else
			sig3 = 0.0_mp
		end if
	end function

!=======================================================
!  Ar + Ar+ -> Ar+ + Ar  Charge X
!=======================================================
	function asigma4(energy) result(sig4)
		real(mp), intent(in) :: energy
		real(mp) :: sig4

		if( energy > 4.0_mp ) then
			sig4 = 2.0E-19 + 5.5E-19/sqrt(energy)
		else
			sig4 = -2.95E-19*sqrt(energy) + 10.65E-19
		end if
	end function

!=======================================================
!  Ar + Ar+ -> Ar + Ar+   Scat.
!=======================================================
	function asigma5(energy) result(sig5)
		real(mp), intent(in) :: energy
		real(mp) :: sig5

		if( energy > 4.0_mp ) then
			sig5 = 1.8E-19 + 4.0E-19/sqrt(energy)
		else
			sig5 = -2.0E-19*sqrt(energy) + 7.8E-19
		end if
	end function

end module