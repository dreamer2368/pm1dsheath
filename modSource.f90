module modSource

	use modPM1D
	use random

	implicit none

contains

!======= [Spatial]_[Velocity] source distribution ==================

	!A(4) number of particle per time step, Uniform on x\in[(0.5-A(1))*L, (0.5+A(1))*L], Maxwellian on v with \sigma = A(2),A(3) (electron,ion)
	subroutine PartialUniform_Maxwellian(pm)
		type(PM1D), intent(inout) :: pm
		integer :: i,nseed,clock
		integer, allocatable :: seed(:)
		integer :: Nadd, newN
		real(mp), allocatable :: xp_add(:), vp_add(:,:), temp_x(:), temp_v(:)

		Nadd = floor(pm%A0(4))
		allocate(xp_add(Nadd))
		allocate(vp_add(2,Nadd))

		call RANDOM_SEED(size=nseed)
		allocate(seed(nseed))
		call SYSTEM_CLOCK(COUNT=clock)
		seed = clock + 59*(/ ( i-1, i=1,nseed ) /)
		call RANDOM_SEED(put=seed)
		deallocate(seed)

		call RANDOM_NUMBER(xp_add)
		xp_add = xp_add*2.0_mp*pm%A0(1)*pm%L + 0.5_mp*pm%L
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

	subroutine Null_source(pm)
		type(PM1D), intent(inout) :: pm
	end subroutine

end module