program main

	use convergence

	implicit none

!	real(mp) :: Bdata(1001)
!	Bdata = (/ ( 0.0002_mp*(i-501), i=1,1001 ) /)

	! print to screen
	print *, 'calling program main'

	call setup					!A,B and QOI are initialized here
!	B = B0 + B0*0.1_mp**(7.0_mp)
	call particle_initialize
	call forwardsweep(langmuir)
!	fDA = (/ ( 9.0_mp-i, i=1,30) /)
!	fDA = EXP(fDA)
!	print *, fDA
!	call error(fDA(18),e)
!	call errorConvergence(fDA,ek)
!	call partialTimeInterval(fDA,e)
!	call energyHistory(langmuir)
!	call Jchaos(Bdata)
	! print to screen
	print *, 'program main...done.'

contains

	! You can add custom subroutines/functions here later, if you want

end program
