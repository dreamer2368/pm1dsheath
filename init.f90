module init

	use modPM1D
	use random
	implicit none

contains

	subroutine twostream_initialize(this,Np,v0,vT,mode)		!generate initial distribution
		type(PM1D), intent(inout) :: this
		integer, intent(in) :: Np, mode
		real(mp), intent(in) :: v0, vT
		real(mp) :: xp0(Np), vp0(Np),rho_back
		real(mp) :: L,qs,ms
		integer :: i,j,N,pm(Np)
		L = this%L
		N = this%n

		qs = -this%wp*this%wp/(this%n*Np/L)
		ms = -qs
		rho_back = -qs*this%n*Np/L

		!velocity distribution initialize
		vp0 = vT*randn(Np)
		pm = (/ ( i, i=1,Np ) /)
		pm = 1 - 2*MOD(pm,2)
		vp0 = vp0 + pm*v0

		do i=1,this%n
			!spatial distribution initialize
			xp0 = (/ ( j*L/Np, j=0,Np-1 ) /) + 0.5*(i-1)*L/Np
			xp0 = xp0 + this%A0*L/Np*SIN( 2.0_mp*pi*xp0/L*mode )

			call buildSpecies(this%p(i),qs,ms,1.0_mp)
			call setSpecies(this%p(i),Np,xp0,vp0)
		end do
		call setMesh(this%m,rho_back*(/ ( 1, i=1,this%m%ng) /))
	end subroutine

	subroutine sheath_initialize(this,Ne,Ni,Te,Ti,Kb)
		type(PM1D), intent(inout) :: this
		integer, intent(in) :: Ne, Ni
		real(mp), intent(in) :: Te, Ti, Kb
		real(mp) :: Vth_e, Vth_i
		real(mp) :: xpe(Ne), vpe(Ne), xpi(Ni), vpi(Ni)
		integer :: i,nseed,clock
		integer, allocatable :: seed(:)

		call RANDOM_SEED(size=nseed)
		allocate(seed(nseed))
		call SYSTEM_CLOCK(COUNT=clock)
		seed = clock + 127*(/ ( i, i=1,nseed ) /)
		call RANDOM_SEED(put=seed)
		call RANDOM_NUMBER(xpe)
		call RANDOM_NUMBER(xpi)
		xpe = xpe*this%L
		xpi = xpi*this%L

		Vth_e = sqrt(2.0_mp*Kb*Te/this%p(1)%ms)
		Vth_i = sqrt(2.0_mp*Kb*Ti/this%p(2)%ms)
		print *, 'Vth_e: ',Vth_e,', Vth_i: ',Vth_i
		vpe = Vth_e/sqrt(2.0_mp)*randn(Ne)
		vpi = Vth_i/sqrt(2.0_mp)*randn(Ni)

		call setSpecies(this%p(1),Ne,xpe,vpe)
		call setSpecies(this%p(2),Ni,xpi,vpi)

		call setMesh(this%m, (/ (0.0_mp, i=1,this%m%ng) /))
	end subroutine

end module