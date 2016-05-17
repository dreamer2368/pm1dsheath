program main

	use init
	use timeStep

	implicit none

	! print to screen
	print *, 'calling program main'

!	call testsheath
!	call testforward(64,5000,1,0.2_mp,20.0_mp,60.0_mp)
!	call test_DST
	call test_Poisson_DN

	! print to screen
	print *, 'program main...done.'

contains

	! You can add custom subroutines/functions here later, if you want

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

		call buildPM1D(sheath,Tf,Time_i,Ng,2,1,order,dt=dt,L=L,eps=eps)
		call buildRecord(r,sheath%nt,2,sheath%L,Ng,'sheath1D',25)

		spwt = n0*L/Ne
		call buildSpecies(sheath%p(1),-qe,me,spwt)
		spwt = n0*L/Ni
		call buildSpecies(sheath%p(2),qe,8*AMU,spwt)

		call sheath_initialize(sheath,Ne,Ni,Te*EV_TO_K,Ti*EV_TO_K,Kb)

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
