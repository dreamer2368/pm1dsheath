module modRecord

	use modPM1D

	implicit none

	type array
		real(mp), allocatable :: vec(:)
	end type

	type recordData
		integer :: nt, n, ng, mod
		integer*4, allocatable :: np(:,:)
		real(mp) :: L, dx
		character(len=:), allocatable :: dir

		type(array), allocatable :: xpdata(:,:)
		type(array), allocatable :: vpdata(:,:)
		type(array), allocatable :: xpsdata(:,:)
		type(array), allocatable :: vpsdata(:,:)
		type(array), allocatable :: Epdata(:,:)

		real(mp), allocatable :: Edata(:,:)
		real(mp), allocatable :: rhodata(:,:)
		real(mp), allocatable :: PE(:), KE(:,:)
	end type

contains

	subroutine buildRecord(this,nt,n,L,ng,input_dir,mod)
		type(recordData), intent(out) :: this
		integer, intent(in) :: nt, n, ng, mod
		real(mp), intent(in) :: L
		character(len=*), intent(in), optional :: input_dir

		this%nt = nt
		this%n = n
		this%L = L
		this%ng = ng
		this%mod = mod

		allocate(this%np(n,nt))
		allocate(this%xpdata(n,nt))
		allocate(this%vpdata(n,nt))
		allocate(this%xpsdata(n,nt))
		allocate(this%vpsdata(n,nt))
		allocate(this%Epdata(n,nt))

		allocate(this%Edata(ng,nt))
		allocate(this%rhodata(ng,nt))

		allocate(this%PE(nt))
		allocate(this%KE(n,nt))

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
		integer :: i,j

		do i=1,this%n
			do j=1, this%nt
				if( allocated(this%xpdata(i,j)%vec) ) then
					deallocate(this%xpdata(i,j)%vec)
				end if
				if( allocated(this%vpdata(i,j)%vec) ) then
					deallocate(this%vpdata(i,j)%vec)
				end if
				if( allocated(this%Epdata(i,j)%vec) ) then
					deallocate(this%Epdata(i,j)%vec)
				end if
!				deallocate(this%xpsdata(i,j)%vec)
!				deallocate(this%vpsdata(i,j)%vec)
			end do
		end do

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

	subroutine recordPlasma(this,pm,k)
		type(recordData), intent(inout) :: this
		type(PM1D), intent(in) :: pm
		integer, intent(in) :: k					!k : time step
		integer :: n								!n : species

		do n=1,pm%n
			this%np(n,k) = pm%p(n)%np
			allocate(this%xpdata(n,k)%vec(pm%p(n)%np))
			allocate(this%vpdata(n,k)%vec(pm%p(n)%np))
			allocate(this%Epdata(n,k)%vec(pm%p(n)%np))

			this%xpdata(n,k)%vec = pm%p(n)%xp
			this%vpdata(n,k)%vec = pm%p(n)%vp
			this%Epdata(n,k)%vec = pm%p(n)%Ep
			this%KE(n,k) = 0.5_mp*SUM(pm%p(n)%ms*(pm%p(n)%vp**2))
		end do

		this%Edata(:,k) = pm%m%E
		this%rhodata(:,k) = pm%m%rho
		this%PE(k) = 0.5_mp*SUM(pm%m%E**2)*pm%m%dx
	end subroutine

	subroutine printPlasma(this,str)
		type(recordData), intent(in) :: this
		character(len=*), intent(in), optional :: str
		character(len=100) :: s
		integer :: i,j

		if( present(str) ) then
			open(unit=300,file='data/'//this%dir//'/record'//str,status='replace')
			open(unit=301,file='data/'//this%dir//'/E'//str//'.bin',status='replace',form='unformatted',access='stream')
			open(unit=302,file='data/'//this%dir//'/rho'//str//'.bin',status='replace',form='unformatted',access='stream')
			open(unit=303,file='data/'//this%dir//'/PE'//str//'.bin',status='replace',form='unformatted',access='stream')
			open(unit=304,file='data/'//this%dir//'/Np'//str//'.bin',status='replace',form='unformatted',access='stream')
			do i=1,this%n
				write(s,*) i
				open(unit=305+3*i,file='data/'//this%dir//'/xp_'//trim(s)//'_'//str//'.bin',status='replace',form='unformatted',access='stream')
				open(unit=306+3*i,file='data/'//this%dir//'/vp_'//trim(s)//'_'//str//'.bin',status='replace',form='unformatted',access='stream')
				open(unit=307+3*i,file='data/'//this%dir//'/KE_'//trim(s)//'_'//str//'.bin',status='replace',form='unformatted',access='stream')
			end do
		else
			open(unit=300,file='data/'//this%dir//'/record',status='replace')
			open(unit=301,file='data/'//this%dir//'/E.bin',status='replace',form='unformatted',access='stream')
			open(unit=302,file='data/'//this%dir//'/rho.bin',status='replace',form='unformatted',access='stream')
			open(unit=303,file='data/'//this%dir//'/PE.bin',status='replace',form='unformatted',access='stream')
			open(unit=304,file='data/'//this%dir//'/Np.bin',status='replace',form='unformatted',access='stream')
			do i=1,this%n
				write(s,*) i
				open(unit=305+3*i,file='data/'//this%dir//'/xp_'//trim(adjustl(s))//'.bin',status='replace',form='unformatted',access='stream')
				open(unit=306+3*i,file='data/'//this%dir//'/vp_'//trim(adjustl(s))//'.bin',status='replace',form='unformatted',access='stream')
				open(unit=307+3*i,file='data/'//this%dir//'/KE_'//trim(adjustl(s))//'.bin',status='replace',form='unformatted',access='stream')
			end do
		end if

		write(300,*) this%n, this%ng, this%nt, this%L, this%mod
		close(300)

		do i = 0,this%nt/this%mod-1
			write(301) this%Edata(:,i*this%mod+1)
			write(302) this%rhodata(:,i*this%mod+1)
			write(303) this%PE(i*this%mod+1)
			write(304) this%np(:,i*this%mod+1)

			do j=1,this%n
				write(305+3*j) this%xpdata(j,i*this%mod+1)%vec
				write(306+3*j) this%vpdata(j,i*this%mod+1)%vec
				write(307+3*j) this%KE(j,i*this%mod+1)
			end do
		end do
		close(301)
		close(302)
		close(303)
		close(304)
		do i=1,this%n
			close(305+3*i)
			close(306+3*i)
			close(307+3*i)
		end do
	end subroutine

end module
