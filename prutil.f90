module prutil     !precision and utils
use f95_precision  ! provide the DP and SP kind of real
use lapack95    !interface to mkl
use blas95 

integer, parameter :: cp=DP  ! or DP

contains 
! mxout = invert(mxin), n*n matrix
subroutine invmat_chol_sp(n,mxin,mxout)
implicit none
integer, intent(in) :: n
integer :: info ! error info of LAPACK Call
real(kind=SP), intent(in), dimension(n,n) :: mxin
real(kind=SP), intent(out), dimension(n,n) :: mxout
integer, allocatable, dimension(:) :: ipiv
integer :: i,j

allocate(ipiv(n))

forall(i=1:n,j=1:n)
mxout(i,j)=mxin(i,j)
end forall
! use cholesky decomposition to calculate the inverse of symmetric positive
! definitive matrix

call potrf(mxout,'U',info)
if (info.ne.0) then
write(*,*) 'info potrf= ', info
!open(unit=101,file='errormxin',status='replace')
!do i=1,n
!write(101,*) mxin(i,:)
!write(101,*) '----------------------------------------------------'
!end do
!close(unit=101)
stop
end if

call potri(mxout,'U',info)
if (info.ne.0) then
write(*,*) 'info potri= ', info
stop
end if

do i=1,n
        do j=1,i
                mxout(i,j)=mxout(j,i)
        end do
end do

deallocate(ipiv)

end subroutine invmat_chol_sp

subroutine invmat_gen(n,mxin,mxout)
implicit none
integer, intent(in) :: n
integer :: info ! error info of LAPACK Call
real(kind=cp), intent(in), dimension(n,n) :: mxin
real(kind=cp), intent(out), dimension(n,n) :: mxout
integer, allocatable, dimension(:) :: ipiv
integer :: i,j

allocate(ipiv(n))

forall(i=1:n,j=1:n)
mxout(i,j)=mxin(i,j)
end forall
! use cholesky decomposition to calculate the inverse of symmetric positive
! definitive matrix

call getrf(mxout,ipiv,info)
if (info.ne.0) then
write(*,*) 'info getrf=',info
stop
end if
call getri(mxout,ipiv,info)
if (info.ne.0) then
write(*,*) 'info getri=',info
end if
deallocate(ipiv)

end subroutine invmat_gen

subroutine mvmul_blas_sp(m,n,mxin,vecin,vecout)
! vecout=mxin*vecin
implicit none
integer, intent(in) :: m,n
real(kind=SP), dimension(m,n), intent(in) :: mxin
real(kind=SP), dimension(n), intent(in) :: vecin
real(kind=SP), dimension(m), intent(out) :: vecout
! prepare vecout for calculation
vecout=0.0_SP
call gemv(mxin,vecin,vecout)

end subroutine mvmul_blas_sp



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
