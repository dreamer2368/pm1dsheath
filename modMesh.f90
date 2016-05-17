module modMesh

	use MatrixVector

	implicit none

	type mesh
		integer :: ng, BCindex
		real(mp) :: L, dx

		real(mp), allocatable :: E(:)
		real(mp), allocatable :: rho(:)
		real(mp), allocatable :: rho_back(:)				!1D sheath: surface charge
		real(mp), allocatable :: phi(:)

		complex(mp), allocatable :: W(:)
	end type

contains

	subroutine buildMesh(this,L,ng,BC)
		type(mesh), intent(out) :: this
		integer, intent(in) :: ng, BC
		real(mp), intent(in) :: L

		this%L = L
		this%ng = ng
		this%BCindex = BC
		select case (this%BCindex)
			case(0)						!periodic
				this%dx = L/ng
			case(1)						!Dirichlet-Dirichlet
				this%dx = L/(ng-1)
		end select

		allocate(this%E(ng))
		allocate(this%phi(ng))
		allocate(this%rho(ng))
		allocate(this%rho_back(ng))
		this%rho_back = 0.0_mp

		allocate(this%W(ng))

		call DSTPoisson_setup(this%ng,this%L,this%W)
	end subroutine

	subroutine setMesh(this,rho_back)
		type(mesh), intent(inout) :: this
		real(mp), intent(in) :: rho_back(this%ng)

		this%rho_back = rho_back
	end subroutine

	subroutine destroyMesh(this)
		type(mesh), intent(inout) :: this

		deallocate(this%E)
		deallocate(this%rho)
		deallocate(this%rho_back)
		deallocate(this%phi)
		deallocate(this%W)
	end subroutine

!===========Mesh Solver===============================================

	subroutine solveMesh(this,eps)
		type(mesh), intent(inout) :: this
		real(mp), intent(in) :: eps

		select case(this%BCindex)
			case(0)
				call solveMesh_periodic(this,eps)
			case(1)
				call solveMesh_D_D(this,eps)
		end select
	end subroutine

	subroutine solveMesh_periodic(this,eps)
		type(mesh), intent(inout) :: this
		real(mp), intent(in) :: eps
		real(mp), dimension(this%ng-1) :: rhs, phi1

		rhs = -( this%rho(1:this%ng-1) + this%rho_back(1:this%ng-1) )/eps
		call CG_K(multiplyK,phi1,rhs,this%dx)
		this%phi(1:this%ng-1) = phi1
		this%phi(this%ng) = 0.0_mp
	end subroutine

	subroutine solveMesh_D_D(this,eps)
		type(mesh), intent(inout) :: this
		real(mp), intent(in) :: eps
		real(mp), dimension(this%ng-2) :: rhs, phi1

		rhs = -( this%rho(2:this%ng-1) + this%rho_back(1:this%ng-1) )/eps
		call CG_K(multiplyK,phi1,rhs,this%dx)
		this%phi(2:this%ng-1) = phi1
		this%phi(1) = 0.0_mp
		this%phi(this%ng) = 0.0_mp
	end subroutine

	subroutine solveMesh_D_N(this,eps)						!D(i=1), N(i=N)
		type(mesh), intent(inout) :: this
		real(mp), intent(in) :: eps
		real(mp), dimension(this%ng-1) :: rhs, phi1

		rhs = -this%rho(2:this%ng)/eps
		call TTA(K_DN,phi1,rhs,this%dx)
		this%phi(2:this%ng-1) = phi1
		this%phi(1) = 0.0_mp
		this%phi(this%ng) = 0.0_mp
	end subroutine

end module