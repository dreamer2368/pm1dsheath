module convergence

use declaration
use init
use timeStep
use assignFunctions

contains

	subroutine error(fDA,ek)

		real(mp), intent(in) :: fDA
		real(mp), intent(out) :: ek
		!Quantity of Interest
		real(mp) :: J0, J1
		real(mp) :: DJDA		!finite approximation of sensitivity
		real(mp) :: delJdelA	!adjoint sensitivity

		call setup					!A,B and QOI are initialized here
		call particle_initialize
		call forwardsweep(langmuir)

		!quantity of interest
		call QOI(langmuir,J0)

		!adjoint
		call adjoint(langmuir,delJdelA)

		!finite approximation
		call destroyPlasma(langmuir)
		B = B0 + B0*fDA									!DEBUG : USE B0 instead of B. Also, consider different amplitudes of B.
		call particle_initialize
		call forwardsweep(langmuir)
		call QOI(langmuir,J1)

		DJDA = (J1-J0)/(B0*fDA)							!DEBUG : USE B0 instead of B
		ek = ABS( ABS(delJdelA-DJDA)/delJdelA )

		call destroyPlasma(langmuir)
		print *, 'J0=',J0,', J1=',J1,', delJdelA=',delJdelA,', DJDA=',DJDA,', ek=',ek
	end subroutine

	subroutine errorConvergence(f,e)
		real(mp), intent(in) :: f(:)
		real(mp), intent(inout) :: e(:)
		!Quantity of Interest
		real(mp) :: J0, J1
		real(mp) :: DJDA		!finite approximation of sensitivity
		real(mp) :: delJdelA	!adjoint sensitivity
		integer :: i

		open(unit=1,file='data/fDA.out',status='replace')
		open(unit=2,file='data/ek.out',status='replace')
		do i=1,size(f)
			call error(f(i),e(i))
			write(1,*) f(i)
			write(2,*) e(i)
		end do
		close(1)
		close(2)
	end subroutine

	subroutine Jchaos(f)
		real(mp), intent(in) :: f(:)										!state space mesh
		real(mp) :: period(8), Jdata(size(f),size(period))						!J
		integer :: i, j, Nperiod(size(period))

		period = (/ 1.0_mp/2.0_mp/pi, 0.5_mp, 1.0_mp, 2.0_mp, 4.0_mp, 6.0_mp, 8.0_mp, 10.0_mp /)

		open(unit=1,file='data/Jdata.out',status='replace')
		call setup					!A is initialized here
		do i=1, size(f)
			B = B0 + B0*f(i)
			print *, 'B=', B
			call particle_initialize
			call forwardsweep(langmuir)
			do j=1, size(period)
				call QOItwo(langmuir,Ni,CEILING((Ti+period(j)*Tp)/dt),Jdata(i,j))
				print *, 'Ni=',Ni,', Nf=', CEILING((Ti+period(j)*Tp)/dt),', J(',i,',',j,')=',Jdata(i,j)
			end do
			write(1,*) Jdata(i,:)
			call destroyPlasma(langmuir)
		end do
		close(1)
	end subroutine

	subroutine partialTimeInterval(f,e)
		real(mp), intent(in) :: f(:)
		real(mp), intent(out) :: e(:,:)
		real(mp) :: freq(8), Tf(size(freq)), J0(size(freq)), delJdelA(size(freq))
		real(mp) :: J1(size(f),size(freq)), DJDA(size(f),size(freq))
		integer :: i,k

		freq = (/ 1.0_mp/2.0_mp/pi, 0.5_mp, 1.0_mp, 2.0_mp, 4.0_mp, 6.0_mp, 8.0_mp, 10.0_mp /)

		open(unit=1,file='data/fDA.out',status='replace')
		write(1,*) f(:)
		close(1)
		open(unit=1,file='data/Tf.out',status='replace')
		write(1,*) freq(:)*Tp
		close(1)

		open(unit=1,file='data/ek_twostream.out',status='replace')
		do k=1,size(freq)
			T = Ti + freq(k)*Tp
			do i=1,size(f)
			call error(f(i),e(k,i))
			end do
			write(1,*) e(k,:)
		end do
		close(1)
	end subroutine

	subroutine energyHistory(this)
		type(plasma), intent(inout) :: this
		real(mp) :: Kinetic, Potential, Total
		integer :: i

		call setup					!A is initialized here
		call particle_initialize
		call buildPlasma(this,Nt,Ni,N,Ng,xp0,vp0)

		!!Time stepping
		call halfStep(this)
		call updatePlasma(this)
		open(unit=1,file='data/energy.out',status='replace')
		write(1,*) 'time	kinetic	potential'
		do i=1,this%nt
			write(1,*) this%time(i),'	',0.5_mp*me*SUM(this%vpdata(:,i)**2),'	',0.5_mp*eps0*dx*SUM(this%Edata(:,i)**2)
		end do
		close(1)
	end subroutine

end module