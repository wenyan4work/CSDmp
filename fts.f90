module fts
! use the intel imkl blas & lapack lib for fortran 95
use prutil
use f95_precision

contains 
subroutine remxcalc(ntotal,grmobmx,rfu,rfe,rsu,rse)
implicit none
integer, intent(in) :: ntotal
real(kind=SP), dimension(11*ntotal,11*ntotal), intent(in) :: grmobmx
real(kind=SP), dimension(6*ntotal,6*ntotal), intent(inout) :: rfu
real(kind=SP), dimension(6*ntotal,5*ntotal), intent(inout) :: rfe
real(kind=SP), dimension(5*ntotal,6*ntotal), intent(inout) :: rsu
real(kind=SP), dimension(5*ntotal,5*ntotal), intent(inout) :: rse
real(kind=SP), dimension(:,:), allocatable :: invmobmx
integer :: n  ! dimension of the 3 grand matrices
integer :: i,j
real(kind=cp) :: sycheck
n=11*ntotal
allocate(invmobmx(n,n))

! invert the grand mobility matrix
call invmat_chol_sp(n,grmobmx,invmobmx)

forall(i=1:6*ntotal,j=1:6*ntotal)
rfu(i,j)=rfu(i,j)+invmobmx(i,j)
end forall

forall(i=1:6*ntotal,j=1:5*ntotal)
rfe(i,j)=rfe(i,j)+invmobmx(i,6*ntotal+j)
end forall

forall(i=1:5*ntotal,j=1:6*ntotal)
rsu(i,j)=rsu(i,j)+invmobmx(6*ntotal+i,j)
end forall

forall(i=1:5*ntotal,j=1:5*ntotal)
rse(i,j)=rse(i,j)+invmobmx(6*ntotal+i,6*ntotal+j)
end forall
end subroutine remxcalc

subroutine hydrostep(ntotal,p_pos,p_u,p_o,p_fe,u_bg,omega_bg,einf_bg &
&,rfu,rfe,timestepsize)
implicit none
integer, intent(in) :: ntotal
real(kind=cp), dimension(ntotal,3), intent(inout) :: p_pos
real(kind=cp), dimension(ntotal,3), intent(out) :: p_u
real(kind=cp), dimension(ntotal,3), intent(out) :: p_o
real(kind=cp), dimension(ntotal,3), intent(in) :: p_fe
real(kind=cp), dimension(3), intent(in) :: u_bg
real(kind=cp), dimension(3), intent(in) :: omega_bg
real(kind=cp), dimension(3,3), intent(in) :: einf_bg
real(kind=SP), dimension(6*ntotal,6*ntotal), intent(in) :: rfu
real(kind=SP), dimension(6*ntotal,5*ntotal), intent(in) :: rfe
real(kind=cp), intent(in) :: timestepsize

real(kind=SP), allocatable, dimension(:,:) :: invrfu
real(kind=SP), allocatable, dimension(:) :: uo_sp
real(kind=cp), allocatable, dimension(:) :: uo_bg
real(kind=SP), allocatable, dimension(:) :: einf_sp
real(kind=SP), allocatable, dimension(:) :: fh_sp
real(kind=cp), dimension(3) :: u_bg_local

integer :: i,j

allocate(invrfu(6*ntotal,6*ntotal))
allocate(uo_sp(6*ntotal))
allocate(uo_bg(6*ntotal))
allocate(fh_sp(6*ntotal))
allocate(einf_sp(5*ntotal))
invrfu=0.0_SP
uo_sp=0.0_SP
uo_bg=0.0_cp
fh_sp=0.0_SP
einf_sp=0.0_SP
!construct the uo_bg, einf vector
do i=1,ntotal
        do j=1,3
        ! U_bg_local=U_bg + Omega_bg * x + Einf_bg . x
        u_bg_local(j)=u_bg(j)+(per(j,1,2)*omega_bg(1)*p_pos(i,2) &
        & +per(j,2,1)*omega_bg(2)*p_pos(i,1)+per(j,1,3)*omega_bg(1)*p_pos(i,3) &
        & +per(j,3,1)*omega_bg(3)*p_pos(i,1)+per(j,2,3)*omega_bg(2)*p_pos(i,3) &
        & +per(j,3,2)*omega_bg(3)*p_pos(i,2)) &
        &+ (einf_bg(j,1)*p_pos(i,1)+einf_bg(j,2)*p_pos(i,2)+einf_bg(j,3)*p_pos(i,3))
        ! fill in the U-Uinf, Omega-Omegainf vector
        uo_bg(3*(i-1)+j)=u_bg_local(j)
        uo_bg(3*ntotal+3*(i-1)+j)=omega_bg(j)
        end do
        ! fill in the -Einf Vector
        ! E vector conversion
        ! EV1=E11-E33, EV2=2E12, EV3=2E13, EV4=2E23, EV5=E22-E33
        einf_sp(5*(i-1)+1) = -einf_bg(3,3)+einf_bg(1,1)
        einf_sp(5*(i-1)+2) = 2.0_cp*einf_bg(1,2)
        einf_sp(5*(i-1)+3) = 2.0_cp*einf_bg(1,3)
        einf_sp(5*(i-1)+4) = 2.0_cp*einf_bg(2,3)
        einf_sp(5*(i-1)+5) = -einf_bg(3,3)+einf_bg(2,2)
end do


!matrix multiplication. Pay attention to the minus sign 
!U = rfu^-1(Fe+RfeEinf)+U_inf
!delta X = U * delta t

!fh=matmul(rfe,einf)
!forall(i=1:ntotal,j=1:3)
!fh(3*(i-1)+j)=fh(3*(i-1)+j)+p_fe(i,j)
!end forall
!uo=matmul(invrfu,fh)
!uo=uo+uo_bg

call invmat_chol_sp(6*ntotal,rfu,invrfu)
call mvmul_blas_sp(6*ntotal,5*ntotal,rfe,einf_sp,fh_sp)
forall(i=1:ntotal,j=1:3)
fh_sp(3*(i-1)+j)=fh_sp(3*(i-1)+j)+p_fe(i,j)
end forall
call mvmul_blas_sp(6*ntotal,6*ntotal,invrfu,fh_sp,uo_sp)
uo_sp=uo_sp+uo_bg

forall(i=1:ntotal,j=1:3)
p_u(i,j)=uo_sp(3*(i-1)+j)
p_o(i,j)=uo_sp(3*ntotal+3*(i-1)+j)
end forall
p_pos=p_pos+p_u*timestepsize

end subroutine hydrostep

subroutine ftscalc(ntotal,p_pos,p_force,p_torque,p_stresslet,p_vel,p_omega, & 
&u_bg,omega_bg,einf_bg,grremx)
implicit none
integer, intent(in) :: ntotal
real(kind=cp), dimension(ntotal,3), intent(in) :: p_pos
real(kind=cp), dimension(ntotal,3), intent(out) :: p_force
real(kind=cp), dimension(ntotal,3), intent(out) :: p_torque
real(kind=cp), dimension(ntotal,5), intent(out) :: p_stresslet
real(kind=cp), dimension(ntotal,3), intent(in) :: p_vel
real(kind=cp), dimension(ntotal,3), intent(in) :: p_omega
real(kind=cp), dimension(3), intent(in) :: u_bg
real(kind=cp), dimension(3), intent(in) :: omega_bg
real(kind=cp), dimension(3,3), intent(in) :: einf_bg
real(kind=cp), dimension(11*ntotal,11*ntotal), intent(in) :: grremx

real(kind=cp), allocatable, dimension(:) :: uoe
real(kind=cp), allocatable, dimension(:) :: fts
real(kind=cp), dimension(3) :: u_bg_local

integer :: i,j

allocate(uoe(11*ntotal))
allocate(fts(11*ntotal))

!construct the u omega e vector
do i=1,ntotal
        do j=1,3
        ! U_bg_local=U_bg + Omega_bg * x + Einf_bg . x
        u_bg_local(j)=u_bg(j)+(per(j,1,2)*omega_bg(1)*p_pos(i,2) &
        & +per(j,2,1)*omega_bg(2)*p_pos(i,1)+per(j,1,3)*omega_bg(1)*p_pos(i,3) &
        & +per(j,3,1)*omega_bg(3)*p_pos(i,1)+per(j,2,3)*omega_bg(2)*p_pos(i,3) &
        & +per(j,3,2)*omega_bg(3)*p_pos(i,2)) &
        &+ (einf_bg(j,1)*p_pos(i,1)+einf_bg(j,2)*p_pos(i,2)+einf_bg(j,3)*p_pos(i,3))
        ! fill in the U-Uinf, Omega-Omegainf vector
        uoe(3*(i-1)+j)=p_vel(i,j)-u_bg_local(j)
        uoe(3*ntotal+3*(i-1)+j)=p_omega(i,j)-omega_bg(j)
        end do
        ! fill in the -Einf Vector
! E vector conversion
! EV1=E11-E33, EV2=2E12, EV3=2E13, EV4=2E23, EV5=E22-E33
! in the leftside, it should be -E, so the elements are -E11+E33, etc
        uoe(6*ntotal+5*(i-1)+1) = einf_bg(3,3)-einf_bg(1,1)
        uoe(6*ntotal+5*(i-1)+2) = -2.0_cp*einf_bg(1,2)
        uoe(6*ntotal+5*(i-1)+3) = -2.0_cp*einf_bg(1,3)
        uoe(6*ntotal+5*(i-1)+4) = -2.0_cp*einf_bg(2,3)
        uoe(6*ntotal+5*(i-1)+5) = einf_bg(3,3)-einf_bg(2,2)
        
end do
!matrix multiplication. Pay attention to the minus sign 
!(F,T,S)=- grremx * (U-U_bg,Omega-Omega_bg,-Einf_bg)
fts = - matmul(grremx,uoe)
!call mvmul_blas(11*ntotal,grremx,uoe,fts) ! step1 : multiplication
!forall(i=1:11*ntotal) ! step2 : the negative sign
!fts(i)=-fts(i)
!end forall

do i=1,ntotal
        do j=1,3
        p_force(i,j)=fts(3*(i-1)+j)
        p_torque(i,j)=fts(3*ntotal+3*(i-1)+j)
        end do
        ! sv1=s11 sv2=s12=s21 sv3=s13=s31 sv4=s23=s32 sv5=s22
        p_stresslet(i,1)=fts(6*ntotal+5*(i-1)+1)
        p_stresslet(i,2)=fts(6*ntotal+5*(i-1)+2)
        p_stresslet(i,3)=fts(6*ntotal+5*(i-1)+3)
        p_stresslet(i,4)=fts(6*ntotal+5*(i-1)+4)
        p_stresslet(i,5)=fts(6*ntotal+5*(i-1)+5)
end do

end subroutine ftscalc

subroutine contactcheck(ntotal,p_pos_local,contactflag)
implicit none
integer,intent(in) :: ntotal
real(kind=cp),dimension(ntotal,3),intent(in) :: p_pos_local
integer, intent(out) :: contactflag
integer :: i,j
real(kind=cp) :: r,r1,r2,r3
contactflag=0

do i=1,ntotal
do j=i+1,ntotal
        r1=p_pos_local(j,1)-p_pos_local(i,1)
        r2=p_pos_local(j,2)-p_pos_local(i,2)
        r3=p_pos_local(j,3)-p_pos_local(i,3)
        r=sqrt(r1*r1+r2*r2+r3*r3)
!        write(*,*) 'r= ',r,i,j
        if (r.lt.2.0_cp) then
        contactflag=1
        exit
        end if
end do
end do

end subroutine contactcheck



end module fts
