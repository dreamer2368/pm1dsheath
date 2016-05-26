module modSource

	use modPM1D
	use random

	implicit none

contains

	subroutine Null_source(pm)
		type(PM1D), intent(inout) :: pm
	end subroutine

!======= [Spatial]_[Velocity] source distribution ==================

	!Uniform on x in [(0.5-A(1))*L, (0.5+A(1))*L]
	!Rayleigh on v with \sigma = A(2),A(3) (electron,ion)
	!Keep A(4) number of ion
	subroutine PartialUniform_Rayleigh(pm)
		type(PM1D), intent(inout) :: pm
		integer :: Nadd, newN, i
		real(mp), allocatable :: xp_add(:), vp_add(:,:), temp_x(:), temp_v(:)

		Nadd = floor(pm%A0(4)-pm%p(2)%np)
		allocate(xp_add(Nadd))
		allocate(vp_add(2,Nadd))

		call RANDOM_NUMBER(xp_add)
		xp_add = xp_add*2.0_mp*pm%A0(1)*pm%L + (0.5_mp - pm%A0(1))*pm%L
		vp_add(1,:) = pm%A0(2)*randr(Nadd)
		vp_add(2,:) = pm%A0(3)*randr(Nadd)

		do i=1,pm%n
			newN = pm%p(i)%np+Nadd
			allocate(temp_x(newN))
			allocate(temp_v(newN))
			temp_x(1:pm%p(i)%np) = pm%p(i)%xp
			temp_x(pm%p(i)%np+1:newN) = xp_add
			temp_v(1:pm%p(i)%np) = pm%p(i)%vp
			temp_v(pm%p(i)%np+1:newN) = vp_add(i,:)

			call destroySpecies(pm%p(i))
			call setSpecies(pm%p(i),newN,temp_x,temp_v)

			deallocate(temp_x)
			deallocate(temp_v)
		end do
	end subroutine

	!Uniform on x in [0, A(3)]
	!Rayleigh on v with \sigma = A(1),A(2) (electron,ion)
	!Keep A(4) number of ion
	subroutine PartialUniform_Rayleigh2(pm)
		type(PM1D), intent(inout) :: pm
		integer :: Nadd, newN, i
		real(mp), allocatable :: xp_add(:), vp_add(:,:), temp_x(:), temp_v(:)

		Nadd = floor(pm%A0(4)-pm%p(2)%np)
		allocate(xp_add(Nadd))
		allocate(vp_add(2,Nadd))

		call RANDOM_NUMBER(xp_add)
		xp_add = xp_add*pm%A0(3)*pm%L
		vp_add(1,:) = pm%A0(1)*randr(Nadd)
		vp_add(2,:) = pm%A0(2)*randr(Nadd)

		do i=1,pm%n
			newN = pm%p(i)%np+Nadd
			allocate(temp_x(newN))
			allocate(temp_v(newN))
			temp_x(1:pm%p(i)%np) = pm%p(i)%xp
			temp_x(pm%p(i)%np+1:newN) = xp_add
			temp_v(1:pm%p(i)%np) = pm%p(i)%vp
			temp_v(pm%p(i)%np+1:newN) = vp_add(i,:)

			call destroySpecies(pm%p(i))
			call setSpecies(pm%p(i),newN,temp_x,temp_v)

			deallocate(temp_x)
			deallocate(temp_v)
		end do
	end subroutine

	!Uniform on x in [(0.5-A(1))*L, (0.5+A(1))*L]
	!Maxwellian on v with \sigma = A(2),A(3) (electron,ion)
	!Keep A(4) number of ion
	subroutine PartialUniform_Maxwellian(pm)
		type(PM1D), intent(inout) :: pm
		integer :: Nadd, newN, i
		real(mp), allocatable :: xp_add(:), vp_add(:,:), temp_x(:), temp_v(:)

		Nadd = floor(pm%A0(4)-pm%p(2)%np)
		allocate(xp_add(Nadd))
		allocate(vp_add(2,Nadd))

		call RANDOM_NUMBER(xp_add)
		xp_add = xp_add*2.0_mp*pm%A0(1)*pm%L + (0.5_mp - pm%A0(1))*pm%L
		vp_add(1,:) = pm%A0(2)*randn(Nadd)
		vp_add(2,:) = pm%A0(3)*randn(Nadd)

		do i=1,pm%n
			newN = pm%p(i)%np+Nadd
			allocate(temp_x(newN))
			allocate(temp_v(newN))
			temp_x(1:pm%p(i)%np) = pm%p(i)%xp
			temp_x(pm%p(i)%np+1:newN) = xp_add
			temp_v(1:pm%p(i)%np) = pm%p(i)%vp
			temp_v(pm%p(i)%np+1:newN) = vp_add(i,:)

			call destroySpecies(pm%p(i))
			call setSpecies(pm%p(i),newN,temp_x,temp_v)

			deallocate(temp_x)
			deallocate(temp_v)
		end do
	end subroutine

	!Uniform on x in [0, A(3)]
	!Maxwellian on v with \sigma = A(1),A(2) (electron,ion)
	!Keep A(4) number of ion
	subroutine PartialUniform_Maxwellian2(pm)
		type(PM1D), intent(inout) :: pm
		integer :: Nadd, newN, i
		real(mp), allocatable :: xp_add(:), vp_add(:,:), temp_x(:), temp_v(:)

		Nadd = floor(pm%A0(4)-pm%p(2)%np)
		allocate(xp_add(Nadd))
		allocate(vp_add(2,Nadd))

		call RANDOM_NUMBER(xp_add)
		xp_add = xp_add*pm%A0(3)*pm%L
		vp_add(1,:) = pm%A0(1)*randn(Nadd)
		vp_add(2,:) = pm%A0(2)*randn(Nadd)

		do i=1,pm%n
			newN = pm%p(i)%np+Nadd
			allocate(temp_x(newN))
			allocate(temp_v(newN))
			temp_x(1:pm%p(i)%np) = pm%p(i)%xp
			temp_x(pm%p(i)%np+1:newN) = xp_add
			temp_v(1:pm%p(i)%np) = pm%p(i)%vp
			temp_v(pm%p(i)%np+1:newN) = vp_add(i,:)

			call destroySpecies(pm%p(i))
			call setSpecies(pm%p(i),newN,temp_x,temp_v)

			deallocate(temp_x)
			deallocate(temp_v)
		end do
	end subroutine

end module