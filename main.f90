program main

	use init
	use timeStep

	implicit none

	! print to screen
	print *, 'calling program main'

!	call testsheath
!	call testforward(64,5000,1,0.2_mp,20.0_mp,60.0_mp)
!	call test_DST
!	call test_Poisson_DN
!	call test_twostream
!	call test_absorbing_boundary
!	call test_refluxing_boundary
!	call test_source
!	call Procassini
!	call test_randn
	call cross_section

	! print to screen
	print *, 'program main...done.'

contains

	! You can add custom subroutines/functions here later, if you want

	subroutine cross_section
		integer, parameter :: N=10000
		real(mp), dimension(N) :: energy, sig1, sig2, sig3, sig4, sig5
		integer :: i

		energy = exp( log(10.0_mp)*( (/ (i,i=1,N) /)/(0.2_mp*N) - 2.0_mp ) )
		do i=1,N
			sig1(i) = asigma1(energy(i))
			sig2(i) = asigma2(energy(i))
			sig3(i) = asigma3(energy(i))
			sig4(i) = asigma4(energy(i))
			sig5(i) = asigma5(energy(i))
		end do

		call system('mkdir -p data/cross_section')
		open(unit=301,file='data/cross_section/sig1.bin',status='replace',form='unformatted',access='stream')
		open(unit=302,file='data/cross_section/sig2.bin',status='replace',form='unformatted',access='stream')
		open(unit=303,file='data/cross_section/sig3.bin',status='replace',form='unformatted',access='stream')
		open(unit=304,file='data/cross_section/sig4.bin',status='replace',form='unformatted',access='stream')
		open(unit=305,file='data/cross_section/sig5.bin',status='replace',form='unformatted',access='stream')
		open(unit=306,file='data/cross_section/energy.bin',status='replace',form='unformatted',access='stream')
		write(301) sig1
		write(302) sig2
		write(303) sig3
		write(304) sig4
		write(305) sig5
		write(306) energy
		close(301)
		close(302)
		close(303)
		close(304)
		close(305)
		close(306)
	end subroutine

	subroutine Procassini
		type(PM1D) :: sheath
		type(recordData) :: r
		real(mp), parameter :: Kb = 1.38065E-23, EV_TO_K = 11604.52_mp, eps = 8.85418782E-12
		real(mp), parameter :: Te = 50.0_mp*EV_TO_K, tau = 100.0_mp
		real(mp), parameter :: me = 9.10938215E-31, qe = 1.602176565E-19, mu = 1836
		real(mp), parameter :: n0 = 2.00000000E14
		integer, parameter :: Ne = 80000, Ni = 80000
		real(mp) :: mi, Ti, wp0, lambda0, dt, dx, L
		real(mp) :: ve0, vi0, Time_f
		real(mp) :: A(4)
		integer :: i

		mi = mu*me
		Ti = Te/tau
		wp0 = sqrt(n0*qe*qe/me/eps)
		lambda0 = sqrt(eps*Kb*Te/n0/qe/qe)
		L = 20.0_mp*lambda0

		print *, 'L = ',L,', lambda0 = ',lambda0,' e = lambda/L = ',lambda0/L

		dt = 0.1_mp/wp0
		dx = 0.2_mp*lambda0
!		dt = 0.5_mp*dx/(lambda0*wp0)

		ve0 = sqrt(Kb*Te/me)
		vi0 = sqrt(Kb*Ti/mi)
		Time_f = 1.0_mp*L/vi0

		A = (/ ve0, vi0, 0.2_mp, 1.0_mp*Ni /)
		call buildPM1D(sheath,Time_f,0.0_mp,ceiling(L/dx),2,pBC=2,mBC=2,order=1,A=A,L=L,dt=dt,eps=eps)
		sheath%wp = wp0
		call buildRecord(r,sheath%nt,2,sheath%L,sheath%ng,'Rayleigh',20)

		call buildSpecies(sheath%p(1),-qe,me,n0*L/Ne)
		call buildSpecies(sheath%p(2),qe,mi,n0*L/Ni)

		call sheath_initialize(sheath,Ne,Ni,Te,Ti,Kb)
		call forwardsweep(sheath,r,Null_input,PartialUniform_Rayleigh2)

		call printPlasma(r)

		call destroyRecord(r)
		call destroyPM1D(sheath)
	end subroutine

	subroutine test_source
		type(PM1D) :: source
		type(recordData) :: r
		integer :: i, N=1
		real(mp) :: rho_back(64)

		call buildPM1D(source,2.0_mp,4.0_mp,64,2,pBC=1,mBC=2,order=1,A=(/0.1_mp, 0.1_mp, 1.0_mp, 10000.0_mp /))
		call buildRecord(r,source%nt,1,source%L,64,'test_source',1)

		call buildSpecies(source%p(1),1.0_mp,1.0_mp,1.0_mp)
		call buildSpecies(source%p(1),-1.0_mp,1.0_mp,1.0_mp)
		call setSpecies(source%p(1),1,0.5_mp*source%L*(/ 1.0_mp /),0.0_mp*(/ 1.0_mp /))
		call setSpecies(source%p(2),1,0.5_mp*source%L*(/ 1.0_mp /),0.0_mp*(/ 1.0_mp /))
		do i=1,N
			call PartialUniform_Rayleigh(source)
		end do

		rho_back = 0.0_mp
		call setMesh(source%m,rho_back)

		call recordPlasma(r,source,1)
		call printPlasma(r)

		call destroyRecord(r)
		call destroyPM1D(source)
	end subroutine

	subroutine test_refluxing_boundary
		type(PM1D) :: reflux
		type(recordData) :: r
		integer, parameter :: Ng=64, N=10000, order=1
		real(mp) :: Ti=20, Tf = 40
		real(mp) :: xp0(N), vp0(N), rho_back(Ng), qe, me
		integer :: i

		call buildPM1D(reflux,Tf,Ti,Ng,1,pBC=2,mBC=2,order=order,A=(/ 1.0_mp, 1.0_mp /))
		call buildRecord(r,reflux%nt,1,reflux%L,Ng,'test_reflux',1)

		xp0 = -0.5_mp*reflux%L
		vp0 = 0.0_mp
		rho_back = 0.0_mp
		qe = -(0.1_mp)**2/(N/reflux%L)
		me = -qe
		rho_back(Ng) = -qe
		call buildSpecies(reflux%p(1),qe,me,1.0_mp)
		call setSpecies(reflux%p(1),N,xp0,vp0)
		call setMesh(reflux%m,rho_back)

		call applyBC(reflux)
		call recordPlasma(r,reflux,1)
		call printPlasma(r)

		call destroyRecord(r)
		call destroyPM1D(reflux)
	end subroutine

	subroutine test_absorbing_boundary
		type(PM1D) :: absorb
		type(recordData) :: r
		integer, parameter :: Ng=64, N=1, order=1
		real(mp) :: Ti=20, Tf = 40
		real(mp) :: xp0(N), vp0(N), rho_back(Ng), qe, me
		integer :: i

		call buildPM1D(absorb,Tf,Ti,Ng,1,pBC=1,mBC=2,order=order)
		call buildRecord(r,absorb%nt,1,absorb%L,Ng,'test_absorb',1)

		xp0(1) = absorb%L + 0.4_mp*absorb%m%dx
		vp0 = 0.0_mp
		rho_back = 0.0_mp
		qe = -(0.1_mp)**2/(N/absorb%L)
		me = -qe
		rho_back(Ng) = -qe
		call buildSpecies(absorb%p(1),qe,me,1.0_mp)
		call setSpecies(absorb%p(1),N,xp0,vp0)
		call setMesh(absorb%m,rho_back)

		call applyBC(absorb)
		call assignMatrix(absorb%a(1),absorb%m,absorb%p(1)%xp)
		call adjustGrid(absorb)
		print *, 'xp = ',absorb%p(1)%xp,', gL = ',absorb%a(1)%g(1,1),', gR = ',absorb%a(1)%g(1,2),	&
			', fL = ',absorb%a(1)%frac(1,1),', fR = ',absorb%a(1)%frac(1,2)

		call destroyRecord(r)
		call destroyPM1D(absorb)
	end subroutine

	subroutine test_twostream
		type(PM1D) :: twostream
		type(recordData) :: r
		integer :: Ng=64, Ne=10000, order=1
		real(mp) :: Ti=40,Tf=80

		call buildPM1D(twostream,Tf,Ti,Ng,1,pBC=0,mBC=0,order=order,dt=0.2_mp)
		call buildRecord(r,twostream%nt,1,twostream%L,Ng,'twostream1D',1)

		call twostream_initialize(twostream,Ne,0.2_mp,0.0_mp,1)

		call forwardsweep(twostream,r,Null_input,Null_source)

		call printPlasma(r)

		call destroyRecord(r)
		call destroyPM1D(twostream)
	end subroutine

	subroutine testsheath
		type(PM1D) :: sheath
		type(recordData) :: r
		integer :: Ng=400, Ne=200000, Ni=100000, order=1
		real(mp) :: n0 = 1.00000E14
		real(mp) :: eps = 8.85418782E-12, Kb = 1.38065E-23
		real(mp) :: dt=5E-11,Time_i=5.0E-9,Tf=1.5E-7
		real(mp) :: Te = 1.0_mp, Ti = 0.1_mp, EV_TO_K = 11604.52_mp		! 1eV in Kelvin, QE/K
		real(mp) :: me = 9.10938215E-31, qe = 1.602176565E-19
		real(mp) :: AMU = 1.660538921E-27
		real(mp) :: L = 0.04_mp, spwt

		call buildPM1D(sheath,Tf,Time_i,Ng,2,pBC=1,mBC=1,order=order,dt=dt,L=L,eps=eps)
		call buildRecord(r,sheath%nt,2,sheath%L,Ng,'sheath1D',25)

		spwt = n0*L/Ne
		call buildSpecies(sheath%p(1),-qe,me,spwt)
		spwt = n0*L/Ni
		call buildSpecies(sheath%p(2),qe,8*AMU,spwt)

		call sheath_initialize(sheath,Ne,Ni,Te*EV_TO_K,Ti*EV_TO_K,Kb)

		call forwardsweep(sheath,r,Null_input,Null_source)

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

		call buildPM1D(twostream,Tf,Ti,Ng,N,pBC=0,mBC=0,order=order,A=(/1.0_mp,0.0_mp/))
		call buildRecord(r,twostream%nt,N,twostream%L,Ng,'species_test',1)

		call twostream_initialize(twostream,Np,0.2_mp,0.0_mp,1)

		call forwardsweep(twostream,r,IC_wave,Null_source)

		call printPlasma(r)

		call destroyRecord(r)
		call destroyPM1D(twostream)
	end subroutine

	subroutine test_randn
		real(mp), allocatable :: test(:)
		integer :: N = 100000
		integer :: i,nseed,clock
		integer, allocatable :: seed(:)
		real(mp) :: temp(1)

		call RANDOM_SEED(size=nseed)
		allocate(seed(nseed))
		call SYSTEM_CLOCK(COUNT=clock)
		seed = 37*(/ ( i-1, i=1,nseed ) /)
		call RANDOM_SEED(put=seed)
		deallocate(seed)

		allocate(test(N))
!		test = randn(N)
		do i=1,N
			temp = randn(1)
			test(i) = temp(1)
		end do
		open(unit=301,file='data/randn.bin',status='replace',form='unformatted',access='stream')
		write(301) test
		close(301)
		deallocate(test)
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

	subroutine test_Poisson_DD
		type(mesh) :: test
		integer :: Ng(3), i, j
		real(mp) :: L, eps
		real(mp), allocatable :: rho(:), phi_sol(:), xg(:)

		L = 1.0_mp
		eps = 1.0_mp
		Ng = (/ 2**5, 2**6, 2**7 /) + 1
		do i=1,3
			allocate(rho(Ng(i)))
			allocate(xg(Ng(i)))
			allocate(phi_sol(Ng(i)))

			call buildMesh(test,L,Ng(i),1)
			xg = (/ (j,j=0,Ng(i)-1) /)*test%dx
			rho = -4.0_mp*pi*pi/L/L*sin(2.0_mp*pi*xg/L)
			test%rho = rho
			call solveMesh(test,eps)
			phi_sol = -sin(2.0_mp*pi*xg/L)
			print *, 'Ng: ',Ng(i),', err: ',maxval(abs(test%phi-phi_sol))

			deallocate(rho)
			deallocate(xg)
			deallocate(phi_sol)
		end do
	end subroutine

	subroutine test_Poisson_DN
		integer :: Ng(3), i, j
		real(mp) :: L, eps, dx, bc
		real(mp), allocatable :: rhs1(:), rhs(:), phi1(:), phi(:), phi_sol(:), xg(:), test(:), co1(:), co2(:), co3(:)
		character(len=100) :: Nstr
		real(mp) :: start,finish

		L = 1.0_mp
		bc = 3.0_mp
		Ng = (/ 2**5, 2**6, 2**7 /) + 1
		do i=1,3
			allocate(rhs(Ng(i)))
			allocate(xg(Ng(i)))
			allocate(phi(Ng(i)))
			allocate(phi_sol(Ng(i)))

			allocate(rhs1(Ng(i)-1))
			allocate(phi1(Ng(i)-1))

			allocate(co1(Ng(i)-1))
			allocate(co2(Ng(i)-1))
			allocate(co3(Ng(i)-1))
			dx = L/(Ng(i)-1)
			co2 = -2.0_mp/dx/dx
			co1 = 1.0_mp/dx/dx
			co3 = 1.0_mp/dx/dx
			co2(Ng(i)-1) = 1.0_mp/dx
			co1(Ng(i)-1) = -1.0_mp/dx

			xg = (/ (j,j=0,Ng(i)-1) /)*dx
			rhs = cos(xg)
			rhs1(1:Ng(i)-2) = rhs(2:Ng(i)-1)
			rhs1(Ng(i)-1) = bc

			call cpu_time(start)
			call solve_tridiag(co1,co2,co3,rhs1,phi1,Ng(i)-1)						!faster
			call cpu_time(finish)
			print *, 'Time: ',finish-start

			phi(2:Ng(i)) = phi1
			phi(1) = 0.0_mp
			phi_sol = 1.0_mp-cos(xg) + (bc-sin(L-0.5_mp*dx))*xg
			print *, 'Ng: ',Ng(i),', err: ',maxval(abs(phi-phi_sol))
			write(Nstr,*) Ng(i)
			open(unit=305,file='data/testphi_'//trim(adjustl(Nstr))//'.bin',status='replace',form='unformatted',access='stream')
			write(305) phi
			close(305)
			open(unit=306,file='data/phisol_'//trim(adjustl(Nstr))//'.bin',status='replace',form='unformatted',access='stream')
			write(306) phi_sol
			close(306)

			deallocate(rhs)
			deallocate(xg)
			deallocate(phi)
			deallocate(phi_sol)
			deallocate(rhs1)
			deallocate(phi1)
			deallocate(co1)
			deallocate(co2)
			deallocate(co3)
		end do
	end subroutine

end program
