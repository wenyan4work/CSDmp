module prutil     !precision and utils

integer, parameter :: SP=4  ! single precision
integer, parameter :: DP=8  ! double precision

integer, parameter :: cp=DP ! current precision used in computation

contains 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Matrix inversion with cholesky
!           Single Precision
! mxout = invert(mxin), n*n matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine invmat_chol_sp(n,mxin,mxout)
implicit none
integer, intent(in) :: n
integer :: info ! error info of LAPACK Call
real(kind=SP), intent(in), dimension(n,n) :: mxin
real(kind=SP), intent(out), dimension(n,n) :: mxout
integer :: i,j

forall(i=1:n,j=1:n)
mxout(i,j)=mxin(i,j)
end forall
! use cholesky decomposition to calculate the inverse of symmetric positive
! definitive matrix

!call potrf(mxout,'U',info)
call SPOTRF('U', n, mxout, n, info)
if (info.ne.0) then
write(*,*) 'info potrf= ', info
stop
end if

!call potri(mxout,'U',info)
call SPOTRI('U', n, mxout, n, info)
if (info.ne.0) then
write(*,*) 'info potri= ', info
stop
end if

do i=1,n
        do j=1,i
                mxout(i,j)=mxout(j,i)
        end do
end do


end subroutine invmat_chol_sp


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Matrix vector multiplication
!           Single Precision
! mxout = invert(mxin), n*n matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine mvmul_blas_sp(m,n,mxin,vecin,vecout)
! vecout=mxin*vecin
implicit none
integer, intent(in) :: m,n
real(kind=SP), dimension(m,n), intent(in) :: mxin
real(kind=SP), dimension(n), intent(in) :: vecin
real(kind=SP), dimension(m), intent(out) :: vecout
! prepare vecout for calculation
vecout=0.0_SP
!call gemv(mxin,vecin,vecout)
call SGEMV('N',m,n,1.0,mxin,m,vecin,1,0.0,vecout,1)

end subroutine mvmul_blas_sp


! Kronecker delta
pure function kd(i,j)
integer, intent(in) :: i,j
integer :: kd
kd=0
if(i.eq.j) then
kd = 1
else 
kd = 0
end if
end function kd

! permutation  symbol
pure function per(i,j,k)
integer,intent(in) :: i,j,k
integer :: per
per=0
if((i.eq.1).and.(j.eq.2).and.(k.eq.3)) then
per = 1
else if((i.eq.2).and.(j.eq.3).and.(k.eq.1)) then
per = 1
else if((i.eq.3).and.(j.eq.1).and.(k.eq.2)) then
per = 1
else if((i.eq.1).and.(j.eq.3).and.(k.eq.2)) then
per = -1
else if((i.eq.2).and.(j.eq.1).and.(k.eq.3)) then
per = -1
else if((i.eq.3).and.(j.eq.2).and.(k.eq.1)) then
per = -1
else
        per=0
end if
end function per

end module prutil
