module modMesh

	use MatrixVector

	implicit none

	type mesh
		integer :: ng
		real(mp) :: L, dx, rho_back

		real(mp), allocatable :: E(:)
		real(mp), allocatable :: rho(:)
		real(mp), allocatable :: phi(:)

		complex(mp), allocatable :: W(:)
	end type

contains

	subroutine buildMesh(this,L,ng)
		type(mesh), intent(out) :: this
		integer, intent(in) :: ng
		real(mp), intent(in) :: L

		this%L = L
		this%ng = ng
		this%dx = L/ng
		this%rho_back = 0.0_mp

		allocate(this%E(ng))
		allocate(this%phi(ng))
		allocate(this%rho(ng))

		allocate(this%W(ng))

		call DSTPoisson_setup(this%ng,this%L,this%W)
	end subroutine

	subroutine setMesh(this,rho_back)
		type(mesh), intent(inout) :: this
		real(mp), intent(in) :: rho_back

		this%rho_back = rho_back
	end subroutine

	subroutine destroyMesh(this)
		type(mesh), intent(inout) :: this

		deallocate(this%E)
		deallocate(this%rho)
		deallocate(this%phi)
		deallocate(this%W)
	end subroutine

end module