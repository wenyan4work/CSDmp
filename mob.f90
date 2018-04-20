! optimized version of mob.f90
! to be tested

module mob
use prutil

contains
subroutine grmobmxcalc(ntotal,p_pos,grmobmx)
implicit none
integer, intent(in) :: ntotal
real(kind=cp), dimension(ntotal,3), intent(in) :: p_pos
real(kind=SP), dimension(11*ntotal,11*ntotal), intent(out) :: grmobmx

integer :: i,j,k,l,alpha,beta
real(kind=cp), allocatable, dimension(:,:,:,:) :: a
real(kind=cp), allocatable, dimension(:,:,:,:) :: b
real(kind=cp), allocatable, dimension(:,:,:,:) :: g
real(kind=cp), allocatable, dimension(:,:,:,:) :: c
real(kind=cp), allocatable, dimension(:,:,:,:) :: h
real(kind=cp), allocatable, dimension(:,:,:,:) :: m

real(kind=cp), dimension(3,3,3) :: gtensor
real(kind=cp), dimension(3,3,3) :: htensor
real(kind=cp), dimension(3,3,3,3) :: mtensor
real(kind=cp), dimension(3,3,3,3) :: mdx,mdy,mdz
real(kind=cp) :: r1,r2,r3,rinv,xa,ya,yb,xc,yc,xg,yg,yh,xm,ym,zm
real(kind=cp), dimension(3) :: e

grmobmx=0.0_cp
gtensor=0.0_cp
htensor=0.0_cp
mtensor=0.0_cp
mdx=0.0_cp
mdy=0.0_cp
mdz=0.0_cp
allocate(a(ntotal,ntotal,3,3))
allocate(b(ntotal,ntotal,3,3))
allocate(c(ntotal,ntotal,3,3))
allocate(g(ntotal,ntotal,5,3))
allocate(h(ntotal,ntotal,5,3))
allocate(m(ntotal,ntotal,5,5))

a=0.0_cp
b=0.0_cp
g=0.0_cp
c=0.0_cp
h=0.0_cp
m=0.0_cp

! self terms
self: do alpha=1,ntotal
        beta=alpha
        e(1)=0
        e(2)=0
        e(3)=0
        rinv=0 !rinv=1/r
! calculate the mob functions
        xa=1.0_cp
        ya=1.0_cp
        yb=0.0_cp
        xc=0.75_cp
        yc=xc
        xg=0.0_cp
        yg=0.0_cp
        yh=0.0_cp
        xm=0.9_cp
        ym=xm
        zm=xm
        ! simplified expression with e1=e2=e3=0 
        forall(i=1:3,j=1:3)
        a(alpha,beta,i,j)=ya*(kd(i,j))
        b(alpha,beta,i,j)=0.0_cp
        c(alpha,beta,i,j)=yc*(kd(i,j))
        end forall
        
        forall(i=1:3,j=1:3,k=1:3)
        gtensor(i,j,k)=0.0_cp
        htensor(i,j,k)=0.0_cp
        end forall
       
        forall(i=1:3,j=1:3,k=1:3,l=1:3)
        mdx(i,j,k,l)=1.5_cp*xm*(-kd(i,j)/3.0_cp)*(-kd(k,l)/3.0_cp)

        mdy(i,j,k,l)=0.0_cp
        
        mdz(i,j,k,l)=0.5_cp*zm*(kd(i,k)*kd(j,l)+kd(j,k)*kd(i,l) &
        & -kd(i,j)*kd(k,l))
        
        mtensor(i,j,k,l)=mdx(i,j,k,l)+mdy(i,j,k,l)+mdz(i,j,k,l)
        end forall
! (U,Omega,-E) = - GrMob * (F,T,S) fluid force on particle
! (F,T,S) = -GrRe * (U,Omega,-E) invert
! E vector conversion
! EV1=E11-E33, EV2=2E12, EV3=2E13, EV4=2E23, EV5=E22-E33
! in the leftside, it should be -E, so the elements are -E11+E33, etc
        g(alpha,beta,1,1)=gtensor(1,1,1)-gtensor(3,3,1)
        g(alpha,beta,1,2)=gtensor(1,1,2)-gtensor(3,3,2)
        g(alpha,beta,1,3)=gtensor(1,1,3)-gtensor(3,3,3)
        g(alpha,beta,2,1)=2.0_cp*gtensor(1,2,1)
        g(alpha,beta,2,2)=2.0_cp*gtensor(1,2,2)
        g(alpha,beta,2,3)=2.0_cp*gtensor(1,2,3)
        g(alpha,beta,3,1)=2.0_cp*gtensor(1,3,1)
        g(alpha,beta,3,2)=2.0_cp*gtensor(1,3,2)
        g(alpha,beta,3,3)=2.0_cp*gtensor(1,3,3)
        g(alpha,beta,4,1)=2.0_cp*gtensor(2,3,1)
        g(alpha,beta,4,2)=2.0_cp*gtensor(2,3,2)
        g(alpha,beta,4,3)=2.0_cp*gtensor(2,3,3)
        g(alpha,beta,5,1)=gtensor(2,2,1)-gtensor(3,3,1)
        g(alpha,beta,5,2)=gtensor(2,2,2)-gtensor(3,3,2)
        g(alpha,beta,5,3)=gtensor(2,2,3)-gtensor(3,3,3)
        
        
        h(alpha,beta,1,1)=htensor(1,1,1)-htensor(3,3,1)
        h(alpha,beta,1,2)=htensor(1,1,2)-htensor(3,3,2)
        h(alpha,beta,1,3)=htensor(1,1,3)-htensor(3,3,3)
        h(alpha,beta,2,1)=2.0_cp*htensor(1,2,1)
        h(alpha,beta,2,2)=2.0_cp*htensor(1,2,2)
        h(alpha,beta,2,3)=2.0_cp*htensor(1,2,3)
        h(alpha,beta,3,1)=2.0_cp*htensor(1,3,1)
        h(alpha,beta,3,2)=2.0_cp*htensor(1,3,2)
        h(alpha,beta,3,3)=2.0_cp*htensor(1,3,3)
        h(alpha,beta,4,1)=2.0_cp*htensor(2,3,1)
        h(alpha,beta,4,2)=2.0_cp*htensor(2,3,2)
        h(alpha,beta,4,3)=2.0_cp*htensor(2,3,3)
        h(alpha,beta,5,1)=htensor(2,2,1)-htensor(3,3,1)
        h(alpha,beta,5,2)=htensor(2,2,2)-htensor(3,3,2)
        h(alpha,beta,5,3)=htensor(2,2,3)-htensor(3,3,3)
        
        ! line 1 of m 5*5
        m(alpha,beta,1,1) = &
        & mtensor(1,1,1,1)-mtensor(1,1,3,3)-mtensor(3,3,1,1)+mtensor(3,3,3,3) 
        m(alpha,beta,1,2) = &
        & mtensor(1,1,1,2)+mtensor(1,1,2,1)-mtensor(3,3,1,2)-mtensor(3,3,2,1) 
        m(alpha,beta,1,3) = &
        & mtensor(1,1,1,3)+mtensor(1,1,3,1)-mtensor(3,3,1,3)-mtensor(3,3,3,1) 
        m(alpha,beta,1,4) = &
        & mtensor(1,1,2,3)+mtensor(1,1,3,2)-mtensor(3,3,2,3)-mtensor(3,3,3,2) 
        m(alpha,beta,1,5) = &
        & mtensor(1,1,2,2)-mtensor(1,1,3,3)-mtensor(3,3,2,2)+mtensor(3,3,3,3) 
        
        ! line 2 of m 5*5
        m(alpha,beta,2,1) = 2.0_cp*(mtensor(1,2,1,1)-mtensor(1,2,3,3))
        m(alpha,beta,2,2) = 2.0_cp*(mtensor(1,2,1,2)+mtensor(1,2,2,1))
        m(alpha,beta,2,3) = 2.0_cp*(mtensor(1,2,1,3)+mtensor(1,2,3,1))
        m(alpha,beta,2,4) = 2.0_cp*(mtensor(1,2,2,3)+mtensor(1,2,3,2))
        m(alpha,beta,2,5) = 2.0_cp*(mtensor(1,2,2,2)-mtensor(1,2,3,3))
        
        ! line 3 of m 5*5
        m(alpha,beta,3,1) = 2.0_cp*(mtensor(1,3,1,1)-mtensor(1,3,3,3))
        m(alpha,beta,3,2) = 2.0_cp*(mtensor(1,3,1,2)+mtensor(1,3,2,1))
        m(alpha,beta,3,3) = 2.0_cp*(mtensor(1,3,1,3)+mtensor(1,3,3,1))
        m(alpha,beta,3,4) = 2.0_cp*(mtensor(1,3,2,3)+mtensor(1,3,3,2))
        m(alpha,beta,3,5) = 2.0_cp*(mtensor(1,3,2,2)-mtensor(1,3,3,3))
        
        ! line 4 of m 5*5
        m(alpha,beta,4,1) = 2.0_cp*(mtensor(2,3,1,1)-mtensor(2,3,3,3))
        m(alpha,beta,4,2) = 2.0_cp*(mtensor(2,3,1,2)+mtensor(2,3,2,1))
        m(alpha,beta,4,3) = 2.0_cp*(mtensor(2,3,1,3)+mtensor(2,3,3,1))
        m(alpha,beta,4,4) = 2.0_cp*(mtensor(2,3,2,3)+mtensor(2,3,3,2))
        m(alpha,beta,4,5) = 2.0_cp*(mtensor(2,3,2,2)-mtensor(2,3,3,3))
        
        ! line 5 of m 5*5
        m(alpha,beta,5,1) = &
        & mtensor(2,2,1,1)-mtensor(2,2,3,3)-mtensor(3,3,1,1)+mtensor(3,3,3,3) 
        m(alpha,beta,5,2) = &
        & mtensor(2,2,1,2)+mtensor(2,2,2,1)-mtensor(3,3,1,2)-mtensor(3,3,2,1) 
        m(alpha,beta,5,3) = &
        & mtensor(2,2,1,3)+mtensor(2,2,3,1)-mtensor(3,3,1,3)-mtensor(3,3,3,1) 
        m(alpha,beta,5,4) = &
        & mtensor(2,2,2,3)+mtensor(2,2,3,2)-mtensor(3,3,2,3)-mtensor(3,3,3,2) 
        m(alpha,beta,5,5) = &
        & mtensor(2,2,2,2)-mtensor(2,2,3,3)-mtensor(3,3,2,2)+mtensor(3,3,3,3) 
        
end do self




! pairwise terms

pairwisealpha: do alpha=1,ntotal
        pairwisebeta: do beta=alpha+1,ntotal

        r1=p_pos(beta,1)-p_pos(alpha,1) 
        r2=p_pos(beta,2)-p_pos(alpha,2) 
        r3=p_pos(beta,3)-p_pos(alpha,3) 
        rinv=1/sqrt(r1*r1+r2*r2+r3*r3)
        e(1)=r1*rinv
        e(2)=r2*rinv
        e(3)=r3*rinv
!calculate the mob functions
! all functions are 12 pair. x12a, y12a, etc...
        xa = 1.5_cp*rinv-rinv**3
        ya = 0.75_cp*rinv+0.5*rinv**3
        yb = -0.75_cp*rinv**2
        xc = 0.75_cp*rinv**3
        yc = -0.375_cp*rinv**3
        xg = 2.25_cp*rinv**2-3.6_cp*rinv**4
        yg = 1.2_cp*rinv**4
        yh = -1.125_cp*rinv**3
        xm = -4.5_cp*rinv**3+10.8_cp*rinv**5
        ym = 2.25_cp*rinv**3-7.2_cp*rinv**5
        zm = 1.8_cp*rinv**5
        
        forall(i=1:3,j=1:3)
        a(alpha,beta,i,j)=xa*e(i)*e(j)+ya*(kd(i,j)-e(i)*e(j))
        b(alpha,beta,i,j)=yb*(per(i,j,1)*e(1)+per(i,j,2)*e(2)+per(i,j,3)*e(3))
        c(alpha,beta,i,j)=xc*e(i)*e(j)+yc*(kd(i,j)-e(i)*e(j))
        end forall
        
        forall(i=1:3,j=1:3)
        a(beta,alpha,i,j)=a(alpha,beta,i,j)
        b(beta,alpha,i,j)=-b(alpha,beta,i,j)
        c(beta,alpha,i,j)=c(alpha,beta,i,j)
        end forall

        forall(i=1:3,j=1:3,k=1:3)
        gtensor(i,j,k)=xg*(e(i)*e(j)-kd(i,j)/3.0_cp)*e(k)+yg*(e(i)*kd(j,k) &
        & +e(j)*kd(i,k)-2.0_cp*e(i)*e(j)*e(k))
       
        htensor(i,j,k)=yh*(e(i)*(per(j,k,1)*e(1)+per(j,k,2)*e(2)+per(j,k,3) &
        & *e(3))+e(j)*(per(i,k,1)*e(1)+per(i,k,2)*e(2)+per(i,k,3)*e(3)))
       
        end forall
       
        forall(i=1:3,j=1:3,k=1:3,l=1:3)
        mdx(i,j,k,l)=1.5_cp*xm*(e(i)*e(j)-kd(i,j)/3.0_cp)*(e(k)*e(l) &
        & -kd(k,l)/3.0_cp)

        mdy(i,j,k,l)=0.5_cp*ym*(e(i)*kd(j,l)*e(k)+e(j)*kd(i,l)*e(k) &
        & +e(i)*kd(j,k)*e(l)+e(j)*kd(i,k)*e(l)-4.0_cp*e(i)*e(j)*e(k)*e(l))
        
        mdz(i,j,k,l)=0.5_cp*zm*(kd(i,k)*kd(j,l)+kd(j,k)*kd(i,l) &
        & -kd(i,j)*kd(k,l)+e(i)*e(j)*kd(k,l)+kd(i,j)*e(k)*e(l) &
        & +e(i)*e(j)*e(k)*e(l)-e(i)*kd(j,l)*e(k)-e(j)*kd(i,l)*e(k) &
        & -e(i)*kd(j,k)*e(l)-e(j)*kd(i,k)*e(l))
        
        mtensor(i,j,k,l)=mdx(i,j,k,l)+mdy(i,j,k,l)+mdz(i,j,k,l)
        end forall
! (U,Omega,-E) = - GrMob * (F,T,S) fluid force on particle
! (F,T,S) = -GrRe * (U,Omega,-E) invert
! E vector conversion
! EV1=E11-E33, EV2=2E12, EV3=2E13, EV4=2E23, EV5=E22-E33
! in the leftside, it should be -E, so the elements are -E11+E33, etc
        g(alpha,beta,1,1)=gtensor(1,1,1)-gtensor(3,3,1)
        g(alpha,beta,1,2)=gtensor(1,1,2)-gtensor(3,3,2)
        g(alpha,beta,1,3)=gtensor(1,1,3)-gtensor(3,3,3)
        g(alpha,beta,2,1)=2.0_cp*gtensor(1,2,1)
        g(alpha,beta,2,2)=2.0_cp*gtensor(1,2,2)
        g(alpha,beta,2,3)=2.0_cp*gtensor(1,2,3)
        g(alpha,beta,3,1)=2.0_cp*gtensor(1,3,1)
        g(alpha,beta,3,2)=2.0_cp*gtensor(1,3,2)
        g(alpha,beta,3,3)=2.0_cp*gtensor(1,3,3)
        g(alpha,beta,4,1)=2.0_cp*gtensor(2,3,1)
        g(alpha,beta,4,2)=2.0_cp*gtensor(2,3,2)
        g(alpha,beta,4,3)=2.0_cp*gtensor(2,3,3)
        g(alpha,beta,5,1)=gtensor(2,2,1)-gtensor(3,3,1)
        g(alpha,beta,5,2)=gtensor(2,2,2)-gtensor(3,3,2)
        g(alpha,beta,5,3)=gtensor(2,2,3)-gtensor(3,3,3)
        
        
        h(alpha,beta,1,1)=htensor(1,1,1)-htensor(3,3,1)
        h(alpha,beta,1,2)=htensor(1,1,2)-htensor(3,3,2)
        h(alpha,beta,1,3)=htensor(1,1,3)-htensor(3,3,3)
        h(alpha,beta,2,1)=2.0_cp*htensor(1,2,1)
        h(alpha,beta,2,2)=2.0_cp*htensor(1,2,2)
        h(alpha,beta,2,3)=2.0_cp*htensor(1,2,3)
        h(alpha,beta,3,1)=2.0_cp*htensor(1,3,1)
        h(alpha,beta,3,2)=2.0_cp*htensor(1,3,2)
        h(alpha,beta,3,3)=2.0_cp*htensor(1,3,3)
        h(alpha,beta,4,1)=2.0_cp*htensor(2,3,1)
        h(alpha,beta,4,2)=2.0_cp*htensor(2,3,2)
        h(alpha,beta,4,3)=2.0_cp*htensor(2,3,3)
        h(alpha,beta,5,1)=htensor(2,2,1)-htensor(3,3,1)
        h(alpha,beta,5,2)=htensor(2,2,2)-htensor(3,3,2)
        h(alpha,beta,5,3)=htensor(2,2,3)-htensor(3,3,3)
        
        ! line 1 of m 5*5
        m(alpha,beta,1,1) = &
        & mtensor(1,1,1,1)-mtensor(1,1,3,3)-mtensor(3,3,1,1)+mtensor(3,3,3,3) 
        m(alpha,beta,1,2) = &
        & mtensor(1,1,1,2)+mtensor(1,1,2,1)-mtensor(3,3,1,2)-mtensor(3,3,2,1) 
        m(alpha,beta,1,3) = &
        & mtensor(1,1,1,3)+mtensor(1,1,3,1)-mtensor(3,3,1,3)-mtensor(3,3,3,1) 
        m(alpha,beta,1,4) = &
        & mtensor(1,1,2,3)+mtensor(1,1,3,2)-mtensor(3,3,2,3)-mtensor(3,3,3,2) 
        m(alpha,beta,1,5) = &
        & mtensor(1,1,2,2)-mtensor(1,1,3,3)-mtensor(3,3,2,2)+mtensor(3,3,3,3) 
        
        ! line 2 of m 5*5
        m(alpha,beta,2,1) = 2.0_cp*(mtensor(1,2,1,1)-mtensor(1,2,3,3))
        m(alpha,beta,2,2) = 2.0_cp*(mtensor(1,2,1,2)+mtensor(1,2,2,1))
        m(alpha,beta,2,3) = 2.0_cp*(mtensor(1,2,1,3)+mtensor(1,2,3,1))
        m(alpha,beta,2,4) = 2.0_cp*(mtensor(1,2,2,3)+mtensor(1,2,3,2))
        m(alpha,beta,2,5) = 2.0_cp*(mtensor(1,2,2,2)-mtensor(1,2,3,3))
        
        ! line 3 of m 5*5
        m(alpha,beta,3,1) = 2.0_cp*(mtensor(1,3,1,1)-mtensor(1,3,3,3))
        m(alpha,beta,3,2) = 2.0_cp*(mtensor(1,3,1,2)+mtensor(1,3,2,1))
        m(alpha,beta,3,3) = 2.0_cp*(mtensor(1,3,1,3)+mtensor(1,3,3,1))
        m(alpha,beta,3,4) = 2.0_cp*(mtensor(1,3,2,3)+mtensor(1,3,3,2))
        m(alpha,beta,3,5) = 2.0_cp*(mtensor(1,3,2,2)-mtensor(1,3,3,3))
        
        ! line 4 of m 5*5
        m(alpha,beta,4,1) = 2.0_cp*(mtensor(2,3,1,1)-mtensor(2,3,3,3))
        m(alpha,beta,4,2) = 2.0_cp*(mtensor(2,3,1,2)+mtensor(2,3,2,1))
        m(alpha,beta,4,3) = 2.0_cp*(mtensor(2,3,1,3)+mtensor(2,3,3,1))
        m(alpha,beta,4,4) = 2.0_cp*(mtensor(2,3,2,3)+mtensor(2,3,3,2))
        m(alpha,beta,4,5) = 2.0_cp*(mtensor(2,3,2,2)-mtensor(2,3,3,3))
        
        ! line 5 of m 5*5
        m(alpha,beta,5,1) = &
        & mtensor(2,2,1,1)-mtensor(2,2,3,3)-mtensor(3,3,1,1)+mtensor(3,3,3,3) 
        m(alpha,beta,5,2) = &
        & mtensor(2,2,1,2)+mtensor(2,2,2,1)-mtensor(3,3,1,2)-mtensor(3,3,2,1) 
        m(alpha,beta,5,3) = &
        & mtensor(2,2,1,3)+mtensor(2,2,3,1)-mtensor(3,3,1,3)-mtensor(3,3,3,1) 
        m(alpha,beta,5,4) = &
        & mtensor(2,2,2,3)+mtensor(2,2,3,2)-mtensor(3,3,2,3)-mtensor(3,3,3,2) 
        m(alpha,beta,5,5) = &
        & mtensor(2,2,2,2)-mtensor(2,2,3,3)-mtensor(3,3,2,2)+mtensor(3,3,3,3) 
        
        forall(i=1:5,j=1:3)
        g(beta,alpha,i,j)=-g(alpha,beta,i,j)
        h(beta,alpha,i,j)=h(alpha,beta,i,j)
        end forall

        forall(i=1:5,j=1:5)
        m(beta,alpha,i,j)=m(alpha,beta,i,j)
        end forall

        end do pairwisebeta
end do pairwisealpha








do alpha=1,ntotal
do beta=1,ntotal
        do i=1,3
        do j=1,3
        !fill a,b,c
        grmobmx(3*(alpha-1)+i,3*(beta-1)+j)=a(alpha,beta,i,j)
        grmobmx(3*ntotal+3*(alpha-1)+i,3*(beta-1)+j)=b(alpha,beta,i,j)
        grmobmx(3*ntotal+3*(alpha-1)+i,3*ntotal+3*(beta-1)+j)=c(alpha,beta,i,j)
        !fill btilda by symmetry
        grmobmx(3*(beta-1)+j,3*ntotal+3*(alpha-1)+i)=b(alpha,beta,i,j)
        end do
        end do
        do i=1,5
        do j=1,3
        !fill g,h
        grmobmx(6*ntotal+5*(alpha-1)+i,3*(beta-1)+j)=g(alpha,beta,i,j)
        grmobmx(6*ntotal+5*(alpha-1)+i,3*ntotal+3*(beta-1)+j)=h(alpha,beta,i,j)
        !fill gtilda,htilda by symmetry
        grmobmx(3*(beta-1)+j,6*ntotal+5*(alpha-1)+i)=g(alpha,beta,i,j)
        grmobmx(3*ntotal+3*(beta-1)+j,6*ntotal+5*(alpha-1)+i)=h(alpha,beta,i,j)
        end do
        end do
        do i=1,5
        do j=1,5
        !fill in m
        grmobmx(6*ntotal+5*(alpha-1)+i,6*ntotal+5*(beta-1)+j)=m(alpha,beta,i,j)
        end do
        end do
end do
end do

!open(unit=101,file='grmobmx',status='replace')
!write(101,*) grmobmx
!close(unit=101)



deallocate(a)
deallocate(b)
deallocate(c)
deallocate(g)
deallocate(h)
deallocate(m)
end subroutine grmobmxcalc





end module mob
