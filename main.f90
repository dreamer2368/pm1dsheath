program main

	use init
	use timeStep

	implicit none

	! print to screen
	print *, 'calling program main'

	call testforward(64,10000,1,0.2_mp,20.0_mp,60.0_mp)
!	call test_DST

	! print to screen
	print *, 'program main...done.'

contains

	! You can add custom subroutines/functions here later, if you want

	subroutine testforward(Ng,N,order,dt,Ti,Tf)
		type(PM1D) :: twostream
		integer, intent(in) :: Ng, N, order
		real(mp), intent(in) :: dt,Ti,Tf
		real(mp) :: xp0(N), vp0(N), qs(N), ms(N), rho_back

		call buildPM1D(twostream,Tf,Ti,Ng,N,order,dir='test')

		call particle_initialize(twostream,0.2_mp,0.0_mp,1,xp0,vp0,qs,ms,rho_back)

		call forwardsweep(twostream,xp0,vp0,qs,ms,rho_back,IC_wave)

		call printPlasma(twostream%r)

		call destroyPM1D(twostream)
	end subroutine

	subroutine test_DST
		type(mesh) :: this
		integer, parameter :: Ng = 64
		real(mp) :: L = 2.0_mp, rhs(Ng)

		call buildMesh(this,L,Ng)

		rhs = 1.0_mp
		call DSTPoisson(this%phi,rhs,this%W)

		open(unit=301,file='data/phi.bin',status='replace',form='unformatted',access='stream')
		write(301) this%phi
		close(301)

		call destroyMesh(this)
	end subroutine

end program
