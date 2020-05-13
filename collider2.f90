subroutine collider(u1,u2,m1,m2,mu)

	use phys

	implicit none

	double precision, intent(in) :: m1, m2, mu
!	double precision, intent(inout), dimension(3) :: u1, u2
	double precision alpha, theta, phi, gamma, randnum, v1, v2, Cg, Sg, Vg
	double precision, dimension(3) :: vp, r, vcm, u1, u2, v1v, v2v, vrel1
	double precision, dimension(3,3) :: Rrg

!	write(*,*) 0.5*m1*sum(u1(:)**2) + 0.5*m2*sum(u2(:)**2)

        theta = acos(mu)

          !* Calculate pair's relative speed

!	vcm(:) = m1/(m1+m2)*u1(:)		! center-of-mass velocity vector, old
	vcm(:) = (m1*u1+m2*u2)/(m1+m2)		! center-of-mass velocity vector, new

!	v1 = sqrt(sum(u1(:)**2))*m2/(m1+m2)	! particle 1 c-o-m scalar velocity, old

	v1v(:) = u1(:) - vcm(:)			! particle 1 c-o-m velocity, new
	v2v(:) = u2(:) - vcm(:)			! particle 2 c-o-m velocity, new
	v1 = sqrt(sum(v1v(:)**2))		! particle 1 c-o-m scalar velocity, new
	v2 = sqrt(sum(v2v(:)**2))		! particle 2 c-o-m scalar velocity, new

	r(:) = u1(:)/sqrt(sum(u1(:)**2))	! unit vector parallel to particle 1 velocity

!	write(*,*) r

!	construct rotation matrix
!	alpha = atan(u1(2)/u1(1))			! old
	alpha = atan2(u1(2),u1(1))			! new

!	phi = atan(u1(3)/sqrt(u1(1)**2+u1(2)**2))       ! old
	phi = atan2(u1(3),sqrt(u1(1)**2+u1(2)**2))      ! new

!	write(*,*) 'alpha', alpha
!	write(*,*) 'phi', phi

	call random_number(randnum)
	gamma = 2d0*pi*randnum			! Collision angle gamma

	vp(1) = v1*cos(alpha)*cos(phi-theta)
	vp(2) = v1*sin(alpha)*cos(phi-theta)
	vp(3) = v1*sin(phi-theta)

!	write(*,*) 'vp', vp(:)

	Cg = cos(gamma)
	Sg = sin(gamma)
	Vg = (1d0-Cg)

	Rrg(1,1) = r(1)*r(1)*Vg+Cg
	Rrg(2,1) = r(1)*r(2)*Vg+r(3)*Sg
	Rrg(3,1) = r(1)*r(3)*Vg-r(2)*Sg

	Rrg(1,2) = r(1)*r(2)*Vg-r(3)*Sg
	Rrg(2,2) = r(2)*r(2)*Vg+Cg
	Rrg(3,2) = r(2)*r(3)*Vg+r(1)*Sg

	Rrg(1,3) = r(1)*r(3)*Vg+r(2)*Sg
	Rrg(2,3) = r(2)*r(3)*Vg-r(1)*Sg
	Rrg(3,3) = r(3)*r(3)*Vg+Cg

!	write(*,*) Rrg

	vrel1 = matmul(Rrg,vp)

!	*** update post-collision velocities ***
	u1(:) = vcm(:) + vrel1
	u2(:) = vcm(:) - m1/m2*vrel1

!	write(*,*) 0.5*m1*sum(u1(:)**2) + 0.5*m2*sum(u2(:)**2)
!	write(*,*)

	return
end subroutine collider
