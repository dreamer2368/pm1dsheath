program main

	use init
	use timeStep

	implicit none

	! print to screen
	print *, 'calling program main'

	call testforward(64,10000,0.2_mp,20.0_mp,60.0_mp)

!	call setup					!A,B and QOI are initialized here
!	B = B0 + B0*EXP(-12.0_mp)
!	call particle_initialize
!	call forwardsweep(langmuir)

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

	subroutine testforward(Ng,N,dt,Ti,Tf)
		type(PM1D) :: twostream
		integer, intent(in) :: Ng, N
		real(mp), intent(in) :: dt,Ti,Tf
		real(mp) :: xp0(N), vp0(N), qs(N), ms(N), rho_back

		call buildPM1D(twostream,Tf,Ti,Ng,N,dir='test')

		call particle_initialize(twostream,0.2_mp,0.0_mp,1,xp0,vp0,qs,ms,rho_back)

		call forwardsweep(twostream,xp0,vp0,qs,ms,rho_back,IC_wave)

		call printPlasma(twostream%r)

		call destroyPM1D(twostream)
	end subroutine

end program
