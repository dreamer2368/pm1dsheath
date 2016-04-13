program main

	use init
	use timeStep

	implicit none

	! print to screen
	print *, 'calling program main'

	call testforward(64,5000,1,0.2_mp,20.0_mp,60.0_mp)
!	call test_DST

	! print to screen
	print *, 'program main...done.'

contains

	! You can add custom subroutines/functions here later, if you want

	subroutine testforward(Ng,Np,order,dt,Ti,Tf)
		type(PM1D) :: twostream
		type(recordData) :: r
		integer, intent(in) :: Ng, Np, order
		real(mp), intent(in) :: dt,Ti,Tf
		integer :: N=2

		call buildPM1D(twostream,Tf,Ti,Ng,N,0,order)
		call buildRecord(r,twostream%nt,N,twostream%L,Ng,'species_test')

		call twostream_initialize(twostream,Np,0.2_mp,0.0_mp,1)

		call forwardsweep(twostream,r,IC_wave)

		call printPlasma(r)

		call destroyRecord(r)
		call destroyPM1D(twostream)
	end subroutine

	subroutine test_DST
		type(mesh) :: this
		integer, parameter :: Ng = 64
		real(mp) :: L = 2.0_mp, rhs(Ng)

		call buildMesh(this,L,Ng,0)

		rhs = 1.0_mp
		call DSTPoisson(this%phi,rhs,this%W)

		open(unit=301,file='data/phi.bin',status='replace',form='unformatted',access='stream')
		write(301) this%phi
		close(301)

		call destroyMesh(this)
	end subroutine

end program
