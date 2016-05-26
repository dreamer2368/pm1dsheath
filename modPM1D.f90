module modPM1D

	use modSpecies
	use modMesh
	use modAssign

	implicit none

	type PM1D
		integer :: nt, ni, n, ng, pBCindex, mBCindex
		real(mp) :: eps0, wp
		real(mp) :: dt, L
		real(mp), allocatable :: A0(:)

		type(species), allocatable :: p(:)
		type(mesh) :: m
		type(pmAssign), allocatable :: a(:)
	end type

contains

	subroutine buildPM1D(this,Tf,Ti,Ng,N,pBC,mBC,order,dt,L,A,eps)
		type(PM1D), intent(out) :: this
		real(mp), intent(in) :: Tf,Ti
		integer, intent(in) :: Ng, N, pBC, mBC, order
		real(mp), intent(in), optional :: dt, A(:), L, eps
		real(mp) :: L0
		integer :: i,nseed
		integer, allocatable :: seed(:)
		if( present(dt) ) then
			this%dt = dt
		else
			this%dt = 0.2_mp
		end if
		if( present(A) ) then
			allocate(this%A0(size(A)))
			this%A0 = A
		else
			allocate(this%A0(1))
			this%A0 = 1.0_mp
		end if
		if( present(L) ) then
			this%L = L
		else
			this%L = 2*pi/( sqrt(3.0_mp)/2.0_mp/sqrt(2.0_mp)/0.2_mp )
		end if
		if( present(eps) ) then
			this%eps0 = eps
		else
			this%eps0 = 1.0_mp
		end if
		this%n = N
		this%ng = Ng
		this%pBCindex = pBC
		this%mBCindex = mBC
		this%nt = CEILING(Tf/this%dt)
		this%dt = Tf/this%nt
		this%ni = FLOOR(Ti/this%dt) + 1
		print *, 'Plasma is created'
		print *, 'L = (',this%L,')'
		print *, 'Ng = (',this%ng,'), N = ',this%n
		print *, 'Particle BC : ', this%pBCindex
		print *, 'Mesh BC : ', this%mBCindex
		print *, 'A = ',this%A0
		print *, 'Ni = ',this%ni,', Nt = ',this%nt,', dt = ',this%dt

		this%wp = 1.0_mp
		allocate(this%p(N))
		call buildMesh(this%m,this%L,Ng,this%mBCindex)
		allocate(this%a(N))
		do i=1,N
			call buildAssign(this%a(i),Ng,order)
		end do

		call RANDOM_SEED(size=nseed)
		allocate(seed(nseed))
		seed = 3433*(/ ( i-1, i=1,nseed ) /)
		call RANDOM_SEED(put=seed)
		deallocate(seed)
	end subroutine

	subroutine destroyPM1D(this)
		type(PM1D), intent(inout) :: this
		integer :: i

		deallocate(this%A0)
		do i=1,this%n
			call destroySpecies(this%p(i))
			call destroyAssign(this%a(i))
		end do
		deallocate(this%p)
		call destroyMesh(this%m)
		deallocate(this%a)
	end subroutine

end module
