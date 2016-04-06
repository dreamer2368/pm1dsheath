module modPM1D

	use modPlasma
	use modMesh
	use modAssign
	use modRecord

	implicit none

	type PM1D
		integer :: nt, ni, n, ng, BCindex
		real(mp) :: L, eps0, wp
		real(mp) :: dt, A0, B0

		type(plasma) :: p
		type(mesh) :: m
		type(recordData) :: r
		type(pmAssign) :: a
	end type

contains

	subroutine buildPM1D(this,Tf,Ti,Ng,N,BC,order,dt,L,A,B,dir)
		type(PM1D), intent(out) :: this
		real(mp), intent(in) :: Tf,Ti
		integer, intent(in) :: Ng, N, BC, order
		real(mp), intent(in), optional :: dt, A, B, L
		character(len=*), intent(in), optional :: dir
		if( present(dt) ) then
			this%dt = dt
		else
			this%dt = 0.2_mp
		end if
		if( present(A) ) then
			this%A0 = A
		else
			this%A0 = 1.0_mp
		end if
		if( present(B) ) then
			this%B0 = B
		else
			this%B0 = 0.0_mp
		end if
		if( present(L) ) then
			this%L = L
		else
			this%L = 2*pi/( sqrt(3.0_mp)/2.0_mp/sqrt(2.0_mp)/0.2_mp )
		end if
		this%ng = Ng
		this%n = N
		this%BCindex = BC
		this%nt = CEILING(Tf/this%dt)
		this%dt = Tf/this%nt
		this%ni = FLOOR(Ti/this%dt) + 1
		print *, 'Plasma is created'
		print *, 'L = (',this%L,')'
		print *, 'Ng = (',Ng,'), N = ',N
		print *, 'BC : ', this%BCindex
		print *, 'A = ',this%A0
		print *, 'Ni = ',this%ni,', Nt = ',this%nt,', dt = ',this%dt

		this%eps0 = 1.0_mp
		this%wp = 1.0_mp
		call buildPlasma(this%p,N)
		call buildMesh(this%m,this%L,Ng)
		call buildAssign(this%a,N,Ng,order)
		if( present(dir) ) then
			call buildRecord(this%r,this%nt,this%n,this%L,this%ng,dir)
		else
			call buildRecord(this%r,this%nt,this%n,this%L,this%ng)
		end if
	end subroutine

	subroutine destroyPM1D(this)
		type(PM1D), intent(inout) :: this

		call destroyPlasma(this%p)
		call destroyMesh(this%m)
		call destroyAssign(this%a)
		call destroyRecord(this%r)
	end subroutine

end module
