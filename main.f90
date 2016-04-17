program main

	use init
	use timeStep

	implicit none

	! print to screen
	print *, 'calling program main'

	call testsheath
!	call testforward(64,5000,1,0.2_mp,20.0_mp,60.0_mp)
!	call test_DST

	! print to screen
	print *, 'program main...done.'

contains

	! You can add custom subroutines/functions here later, if you want

	subroutine testsheath
		type(PM1D) :: sheath
		type(recordData) :: r
		integer :: Ng=128, Ne=100000, Ni=50000, order=1
		real(mp) :: n0 = 1.00000E-14
		real(mp) :: eps = 8.85418782E-12, Kb = 1.38065E-23
		real(mp) :: dt=5E-11,Time_i=5.0E-7,Tf=1.0E-6
		real(mp) :: Te = 1.0_mp, Ti = 0.1_mp
		real(mp) :: me = 9.10938215E-31, qe = 1.602176565E-19
		real(mp) :: AMU = 1.660538921E-27
		real(mp) :: L = 0.04_mp, spwt

		call buildPM1D(sheath,Tf,Time_i,Ng,2,1,order,dt=dt,L=L)
		call buildRecord(r,sheath%nt,2,sheath%L,Ng,'sheath1D',100)

		spwt = n0*L/Ne
		call buildSpecies(sheath%p(1),-qe,me,spwt)
		spwt = n0*L/Ni
		call buildSpecies(sheath%p(2),qe,8*AMU,spwt)

		call sheath_initialize(sheath,Ne,Ni,Te,Ti,Kb)

		call forwardsweep(sheath,r,Null_input)

		call printPlasma(r)

		call destroyRecord(r)
		call destroyPM1D(sheath)
	end subroutine

	subroutine testforward(Ng,Np,order,dt,Ti,Tf)
		type(PM1D) :: twostream
		type(recordData) :: r
		integer, intent(in) :: Ng, Np, order
		real(mp), intent(in) :: dt,Ti,Tf
		integer :: N=2

		call buildPM1D(twostream,Tf,Ti,Ng,N,0,order)
		call buildRecord(r,twostream%nt,N,twostream%L,Ng,'species_test',1)

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
