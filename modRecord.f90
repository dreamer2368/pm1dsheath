module modRecord

	use modPlasma
	use modMesh

	implicit none

	type recordData
		integer :: nt, n, ng
		real(mp) :: L, dx
		character(len=:), allocatable :: dir

		real(mp), allocatable :: xpdata(:,:)
		real(mp), allocatable :: vpdata(:,:)
		real(mp), allocatable :: xpsdata(:,:)
		real(mp), allocatable :: vpsdata(:,:)

		real(mp), allocatable :: Epdata(:,:)
		real(mp), allocatable :: Edata(:,:)
		real(mp), allocatable :: rhodata(:,:)
		real(mp), allocatable :: PE(:), KE(:)
	end type

contains

	subroutine buildRecord(this,nt,n,L,ng,input_dir)
		type(recordData), intent(out) :: this
		integer, intent(in) :: nt, n, ng
		real(mp), intent(in) :: L
		character(len=*), intent(in), optional :: input_dir

		this%nt = nt
		this%n = n
		this%L = L
		this%ng = ng

		allocate(this%xpdata(n,nt))
		allocate(this%vpdata(n,nt))
		allocate(this%xpsdata(n,nt))
		allocate(this%vpsdata(n,nt))

		allocate(this%Epdata(n,nt))
		allocate(this%Edata(ng,nt))
		allocate(this%rhodata(ng,nt))

		allocate(this%PE(nt))
		allocate(this%KE(nt))

		if( present(input_dir) ) then
			allocate(character(len=len(input_dir)) :: this%dir)
			this%dir = input_dir
		else
			allocate(character(len=0) :: this%dir)
			this%dir = ''
		end if

		call system('mkdir -p data/'//this%dir)
	end subroutine

	subroutine destroyRecord(this)
		type(recordData), intent(inout) :: this

		deallocate(this%xpdata)
		deallocate(this%vpdata)
		deallocate(this%xpsdata)
		deallocate(this%vpsdata)

		deallocate(this%Epdata)
		deallocate(this%Edata)
		deallocate(this%rhodata)

		deallocate(this%PE)
		deallocate(this%KE)

		deallocate(this%dir)
	end subroutine

	subroutine recordPlasma(this,p,m,k)
		type(recordData), intent(inout) :: this
		type(plasma), intent(in) :: p
		type(mesh), intent(in) :: m
		integer, intent(in) :: k

		this%xpdata(:,k) = p%xp
		this%vpdata(:,k) = p%vp
		this%Epdata(:,k) = p%Ep
		this%Edata(:,k) = m%E
		this%rhodata(:,k) = m%rho
		this%PE(k) = 0.5_mp*SUM(m%E**2)*m%dx
		this%KE(k) = 0.5_mp*SUM(p%ms*(p%vp**2))
	end subroutine

	subroutine printPlasma(this,str)
		type(recordData), intent(in) :: this
		character(len=*), intent(in), optional :: str
		integer :: i

		if( present(str) ) then
			open(unit=300,file='data/'//this%dir//'/record'//str,status='replace')
			open(unit=301,file='data/'//this%dir//'/xp'//str//'.bin',status='replace',form='unformatted',access='stream')
			open(unit=302,file='data/'//this%dir//'/vp'//str//'.bin',status='replace',form='unformatted',access='stream')
			open(unit=303,file='data/'//this%dir//'/E'//str//'.bin',status='replace',form='unformatted',access='stream')
			open(unit=304,file='data/'//this%dir//'/rho'//str//'.bin',status='replace',form='unformatted',access='stream')
			open(unit=305,file='data/'//this%dir//'/PE'//str//'.bin',status='replace',form='unformatted',access='stream')
			open(unit=306,file='data/'//this%dir//'/KE'//str//'.bin',status='replace',form='unformatted',access='stream')
		else
			open(unit=300,file='data/'//this%dir//'/record',status='replace')
			open(unit=301,file='data/'//this%dir//'/xp.bin',status='replace',form='unformatted',access='stream')
			open(unit=302,file='data/'//this%dir//'/vp.bin',status='replace',form='unformatted',access='stream')
			open(unit=303,file='data/'//this%dir//'/E.bin',status='replace',form='unformatted',access='stream')
			open(unit=304,file='data/'//this%dir//'/rho.bin',status='replace',form='unformatted',access='stream')
			open(unit=305,file='data/'//this%dir//'/PE.bin',status='replace',form='unformatted',access='stream')
			open(unit=306,file='data/'//this%dir//'/KE.bin',status='replace',form='unformatted',access='stream')
		end if

		write(300,*) this%n, this%ng, this%nt, this%L
		close(300)

		do i = 1,this%nt
			write(301) this%xpdata(:,i)
			write(302) this%vpdata(:,i)

!			write(303) this%Epdata(:,i)
			write(303) this%Edata(:,i)
			write(304) this%rhodata(:,i)
			write(305) this%PE(i)
			write(306) this%KE(i)
		end do
		close(301)
		close(302)
		close(303)
		close(304)
		close(305)
		close(306)
	end subroutine

end module
