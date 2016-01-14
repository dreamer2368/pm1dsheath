module declaration

	use modPlasma

	implicit none

	! parameters setup
	integer, parameter :: N = 10000
	real(mp), parameter :: pi = 4.0_mp*ATAN(1.0_mp)
	real(mp), parameter :: L = 2*pi/( sqrt(3.0_mp)/2.0_mp/sqrt(2.0_mp)/0.2_mp )
	real(mp), parameter :: eps0 = 1.0_mp
	real(mp), parameter :: wp = 1.0_mp
	real(mp), parameter :: qe = -wp*wp/(N/L)
	real(mp), parameter :: me = -qe
	real(mp), parameter :: rho_back = -qe*N/L

	real(mp) :: Ti = 40.0_mp
	real(mp) :: Tp = 2.0_mp*pi/wp
!	real(mp) :: T = 40.0_mp
!	real(mp) :: T = 40.0_mp + 10.0_mp*(2.0_mp*pi/wp)
	real(mp) :: T = 20.0_mp*(2.0_mp*pi/wp)

	!for loop index
	integer :: i

	! time step parameter
	real(mp) :: dt
	integer :: Nt, Ni

	!initial spatial distribution
	real(mp) :: xp0(N)
	real(mp), parameter :: A0 = 1.0_mp, B0 = 0.0_mp
	real(mp) :: A = A0
	real(mp) :: B = B0
	integer, parameter :: mode = 1

	!initial velocity distribution
	real(mp), parameter :: vT = 0.0_mp
	real(mp), parameter :: v0 = 0.0_mp
	real(mp) :: vp0(N)
	integer :: pm(N)

	!grid and operators setup
	integer, parameter :: Ng = 64
	real(mp) :: dx
	real(mp) :: xg(Ng)
	real(mp) :: rhs(Ng-1)

	!time-stepping particle & field variable are composed of type fluid
	type(plasma) :: langmuir

	!error convergence variable
	real(mp) :: fDA(30)
	real(mp) :: ek(30)
!	real(mp) :: e
	real(mp) :: e(8,30)

!	character(32) :: filename ! You can make this longer than 32 characters, if needed

end module