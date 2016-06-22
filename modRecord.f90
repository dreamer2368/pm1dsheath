module modRecord

	use modPM1D

	implicit none

	type recordData
		integer :: nt, n, ng, mod
		real(mp) :: L, dx
		character(len=:), allocatable :: dir

		integer*4, allocatable :: np(:,:)
		real(mp), allocatable :: phidata(:,:)
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
		integer :: nr
		nr = (nt-1)/mod+1

		this%nt = nt
		this%n = n
		this%L = L
		this%ng = ng
		this%mod = mod

		allocate(this%np(n,nr))

		allocate(this%phidata(ng,nr))
		allocate(this%Edata(ng,nr))
		allocate(this%rhodata(ng,nr))

		allocate(this%PE(nr))
		allocate(this%KE(n,nr))

		if( present(input_dir) ) then
			allocate(character(len=len(input_dir)) :: this%dir)
			this%dir = input_dir
		else
			allocate(character(len=0) :: this%dir)
			this%dir = ''
		end if

		call system('mkdir -p data/'//this%dir//'/xp')
		call system('mkdir -p data/'//this%dir//'/vp')

		call system('rm data/'//this%dir//'/xp/*.*')
		call system('rm data/'//this%dir//'/vp/*.*')
	end subroutine

	subroutine destroyRecord(this)
		type(recordData), intent(inout) :: this
		integer :: i,j

		deallocate(this%np)
		deallocate(this%phidata)
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
		character(len=100) :: nstr, kstr
		real(mp) :: qe = 1.602176565E-19

		if( (this%mod.eq.1) .or. (mod(k,this%mod).eq.1) ) then
			do n=1,pm%n
				write(nstr,*) n
				if( this%mod .eq. 1 ) then
					write(kstr,*) k
				else
					write(kstr,*) k/this%mod + 1
				end if
				open(unit=305,file='data/'//this%dir//'/xp/'//trim(adjustl(kstr))//'_'	&
					//trim(adjustl(nstr))//'.bin',status='replace',form='unformatted',access='stream')
				open(unit=306,file='data/'//this%dir//'/vp/'//trim(adjustl(kstr))//'_'	&
					//trim(adjustl(nstr))//'.bin',status='replace',form='unformatted',access='stream')
				write(305) pm%p(n)%xp
				write(306) pm%p(n)%vp
				close(305)
				close(306)
				this%np(n,k/this%mod+1) = pm%p(n)%np
				this%KE(n,k/this%mod+1) = 0.5_mp*SUM(pm%p(n)%ms*pm%p(n)%spwt*(pm%p(n)%vp**2))
			end do

			this%phidata(:,k/this%mod+1) = pm%m%phi
			this%Edata(:,k/this%mod+1) = pm%m%E
			this%rhodata(:,k/this%mod+1) = pm%m%rho
			this%PE(k/this%mod+1) = 0.5_mp*SUM(pm%m%E**2)*pm%m%dx
			print *, '============= ',k,'-th Time Step ================='
			do n=1,pm%n
				print *, 'Species(',n,'): ',pm%p(n)%np, ', KE: ', 0.5_mp*pm%p(n)%ms*sum((pm%p(n)%vp)**2)
			end do
		end if
	end subroutine

	subroutine printPlasma(this)
		type(recordData), intent(in) :: this
		character(len=100) :: s
		integer :: i,j

		open(unit=300,file='data/'//this%dir//'/record',status='replace')
		open(unit=301,file='data/'//this%dir//'/E.bin',status='replace',form='unformatted',access='stream')
		open(unit=302,file='data/'//this%dir//'/rho.bin',status='replace',form='unformatted',access='stream')
		open(unit=303,file='data/'//this%dir//'/PE.bin',status='replace',form='unformatted',access='stream')
		open(unit=304,file='data/'//this%dir//'/Np.bin',status='replace',form='unformatted',access='stream')
		open(unit=305,file='data/'//this%dir//'/phi.bin',status='replace',form='unformatted',access='stream')
		do i=1,this%n
			write(s,*) i
			open(unit=307+i,file='data/'//this%dir//'/KE_'//trim(adjustl(s))//'.bin',status='replace',form='unformatted',access='stream')
		end do

		write(300,*) this%n, this%ng, this%nt, this%L, this%mod
		close(300)

		do i = 1,this%nt/this%mod
			write(301) this%Edata(:,i)
			write(302) this%rhodata(:,i)
			write(303) this%PE(i)
			write(304) this%np(:,i)
			write(305) this%phidata(:,i)
			do j=1,this%n
				write(307+j) this%KE(j,i)
			end do
		end do

		close(301)
		close(302)
		close(303)
		close(304)
		close(305)
		do i=1,this%n
			close(307+i)
		end do
	end subroutine

end module
