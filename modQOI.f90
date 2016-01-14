module modQOI

	use modPlasma

	implicit none

	type QOI
		integer :: Ni
		integer :: Nf
		real(mp) :: value
	end type

contains

	subroutine buildQOI(this,Ti,Tf,dt)
		real(mp), intent(in) :: Ti, Tf, dt
		type(QOI), intent(inout) :: this
		this%Ni = FLOOR(Ti/dt) + 1
		this%Nf = CEILING(Tf/dt)
		print *,'QOI_Ni=', this%Ni,', QOI_Nf=',this%Nf
	end subroutine

	subroutine EvalQOI(plasmaA,this)
		type(plasma), intent(in) :: plasmaA
		type(QOI), intent(inout) :: this
		if (this%Nf>plasmaA%nt) then
			print *,'================================================'
			print *,'=========FAULT : Nf > total time step==========='
			print *,'================================================'
		end if
		this%value = 1.0_mp/plasmaA%n/(this%Nf-this%Ni)*SUM( plasmaA%vpdata(:,this%Ni+1:this%Nf)**2 )
	end subroutine

	function AdjQOI(plasmaA,this,dt,k) result(dJ)
		type(plasma), intent(in) :: plasmaA
		type(QOI), intent(in) :: this
		integer, intent(in) :: k
		real(mp), intent(in) :: dt
		integer :: nk
		real(mp), allocatable :: dJ(:)

		allocate(dJ(plasmaA%n))

		nk = this%Nf+1-k
		if( nk>=this%Ni+1 ) then
			dJ = 2.0_mp/plasmaA%n/(this%Nf-this%Ni)/dt*plasmaA%vpdata(:,nk)
		else
			dJ = 0.0_mp
		end if
	end function

end module