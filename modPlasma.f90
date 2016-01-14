module modPlasma

	implicit none

	integer, parameter :: mp = SELECTED_REAL_KIND(15)

	type plasma
		integer :: nt
		integer :: ni
		integer :: n
		integer :: ng
		real(mp), allocatable :: time(:)
		real(mp), allocatable :: xp(:)
		real(mp), allocatable :: vp(:)
		real(mp), allocatable :: E(:)
		real(mp), allocatable :: rho(:)
		real(mp), allocatable :: phi(:)

!		real(mp), allocatable :: mat(:,:)
		integer, allocatable :: g(:)
		integer, allocatable :: gl(:)
		integer, allocatable :: gr(:)
		real(mp), allocatable :: fl(:)
		real(mp), allocatable :: f0(:)
		real(mp), allocatable :: fr(:)
		real(mp), allocatable :: h(:)

		real(mp), allocatable :: Ep(:)
		real(mp), allocatable :: xpdata(:,:)
		real(mp), allocatable :: vpdata(:,:)
		real(mp), allocatable :: Edata(:,:)

		real(mp), allocatable :: xpsdata(:,:)
	end type

contains

	subroutine buildPlasma(this,nt,ni,n,ng,xp0,vp0)

		type(plasma), intent(out) :: this
		integer, intent(in) :: nt
		integer, intent(in) :: ni
		integer, intent(in) :: n
		integer, intent(in) :: ng
		real(mp), intent(in) :: xp0(:)
		real(mp), intent(in) :: vp0(:)

		this%nt = nt
		this%ni = ni
		this%n = n
		this%ng = ng

		allocate(this%time(nt))
		allocate(this%xp(n))
		allocate(this%vp(n))
		allocate(this%E(ng))
		allocate(this%phi(ng))
		allocate(this%rho(ng))

!		allocate(this%mat(ng,n))
		allocate(this%g(n))
		allocate(this%gl(n))
		allocate(this%gr(n))
		allocate(this%fl(n))
		allocate(this%f0(n))
		allocate(this%fr(n))
		allocate(this%h(n))

		allocate(this%Ep(n))

		this%xp = xp0
		this%vp = vp0
		this%g = 0

		allocate(this%xpdata(n,nt))
		allocate(this%vpdata(n,nt))
		allocate(this%Edata(ng,nt))

		allocate(this%xpsdata(n,nt))

	end subroutine

	subroutine destroyPlasma(this)
		type(plasma), intent(inout) :: this
		deallocate(this%time)
		deallocate(this%xp)
		deallocate(this%vp)
		deallocate(this%E)
		deallocate(this%rho)
		deallocate(this%phi)

!		deallocate(this%mat)
		deallocate(this%g)
		deallocate(this%gl)
		deallocate(this%gr)
		deallocate(this%f0)
		deallocate(this%fl)
		deallocate(this%fr)
		deallocate(this%h)

		deallocate(this%Ep)
		deallocate(this%xpdata)
		deallocate(this%vpdata)
		deallocate(this%Edata)

		deallocate(this%xpsdata)
	end subroutine

	subroutine recordPlasma(this,k)

		type(plasma), intent(inout) :: this
		integer, intent(in) :: k

		this%xpdata(:,k) = this%xp
		this%vpdata(:,k) = this%vp
		this%Edata(:,k) = this%E

	end subroutine

	subroutine printPlasma(this,L,dt)
		type(plasma), intent(in) :: this
		real(mp), intent(in) :: L, dt

		open(unit=1,file='data/record.out',status='replace')
		write(1,*) this%nt
		write(1,*) this%n
		write(1,*) this%ng
		write(1,*) L
		write(1,*) dt
		close(1)

		open(unit=1,file='data/xp.bin',status='replace',form='unformatted',access='stream')
		open(unit=2,file='data/vp.bin',status='replace',form='unformatted',access='stream')
		open(unit=3,file='data/E.bin',status='replace',form='unformatted',access='stream')
		write(1) this%xpdata
		write(2) this%vpdata
		write(3) this%Edata
		close(1)
		close(2)
		close(3)
	end subroutine

	subroutine recordAdjoint(this,k,xps)
		type(plasma), intent(inout) :: this
		integer, intent(in) :: k
		real(mp), intent(in) :: xps(:)

		this%xpsdata(:,k) = xps
	end subroutine

end module