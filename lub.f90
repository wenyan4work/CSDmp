module lub
use prutil
use f95_precision
use in_out

contains 
! calculate the lubrication contribution to rfu rfe rsu rse
! rfu is symmetric
! rfe = transpose of rsu
! rse is symmetric
subroutine lubmxcalc(ntotal,p_pos,rfu,rfe,rsu,rse)
implicit none
integer, intent(in) :: ntotal
real(kind=cp), dimension(ntotal,3), intent(in) :: p_pos
real(kind=SP), dimension(6*ntotal,6*ntotal), intent(out) :: rfu
real(kind=SP), dimension(6*ntotal,5*ntotal), intent(out) :: rfe
real(kind=SP), dimension(5*ntotal,6*ntotal), intent(out) :: rsu
real(kind=SP), dimension(5*ntotal,5*ntotal), intent(out) :: rse

real(kind=cp) :: r,r1,r2,r3
real(kind=cp) :: x11a,x12a,y11a,y12a,y11b,y12b,x11c,x12c,y11c,y12c
real(kind=cp) :: x11g,x12g,y11g,y12g,y11h,y12h,xm,ym,zm
real(kind=cp), dimension(3) :: e
real(kind=cp), dimension(22,22) :: lub2bmx
real(kind=cp), dimension(2,2,3,3) :: alub
real(kind=cp), dimension(2,2,3,3) :: blub
real(kind=cp), dimension(2,2,3,3) :: club
real(kind=cp), dimension(2,2,3,3,3) :: glubtensor
real(kind=cp), dimension(2,2,3,3,3) :: hlubtensor
real(kind=cp), dimension(3,3,3,3) :: mlubtensor
real(kind=cp), dimension(5,5) :: mmx


integer :: alpha,beta
integer :: i,j,k,l

rfu=0.0_SP
rfe=0.0_SP
rsu=0.0_SP
rse=0.0_SP


do alpha=1,ntotal
        do beta=alpha+1,ntotal
        r1=p_pos(beta,1)-p_pos(alpha,1)
        r2=p_pos(beta,2)-p_pos(alpha,2)
        r3=p_pos(beta,3)-p_pos(alpha,3)
        r=sqrt(r1*r1+r2*r2+r3*r3)
        lub2bmx=0.0_cp
        alub=0.0_cp
        blub=0.0_cp
        club=0.0_cp
        glubtensor=0.0_cp
        hlubtensor=0.0_cp
        mlubtensor=0.0_cp
        
        if(r.lt.2.0_cp) then
                write(*,*) 'Contact Error, r='
                write(*,*) r,alpha,beta
                stop
        else if(r.gt.4.0_cp) then
                lub2bmx=0.0_cp                
        else
        ! initialize the lub2bmx matrix
        ! build the 2body lubrication resistance matrix lub2bmx(22*22) 
        call lubcalc(r,x11a,x12a,y11a,y12a,y11b,y12b,x11c,x12c,y11c,y12c,&
                        & x11g,x12g,y11g,y12g,y11h,y12h,xm,ym,zm)
        
        e(1)=r1/r
        e(2)=r2/r
        e(3)=r3/r
        
        forall(i=1:3,j=1:3)
        alub(1,1,i,j)=x11a*e(i)*e(j)+y11a*(kd(i,j)-e(i)*e(j))
        alub(1,2,i,j)=x12a*e(i)*e(j)+y12a*(kd(i,j)-e(i)*e(j))
        blub(1,1,i,j)=y11b*(per(i,j,1)*e(1)+per(i,j,2)*e(2)+per(i,j,3)*e(3))
        blub(1,2,i,j)=y12b*(per(i,j,1)*e(1)+per(i,j,2)*e(2)+per(i,j,3)*e(3))
        club(1,1,i,j)=x11c*e(i)*e(j)+y11c*(kd(i,j)-e(i)*e(j))
        club(1,2,i,j)=x12c*e(i)*e(j)+y12c*(kd(i,j)-e(i)*e(j))
        end forall

        forall(i=1:3,j=1:3,k=1:3)
        glubtensor(1,1,i,j,k)=x11g*(e(i)*e(j)-kd(i,j)/3.0_cp)*e(k)& 
        &+y11g*(e(i)*kd(j,k)+e(j)*kd(i,k)-2*e(i)*e(j)*e(k))
        glubtensor(1,2,i,j,k)=x12g*(e(i)*e(j)-kd(i,j)/3.0_cp)*e(k)& 
        &+y12g*(e(i)*kd(j,k)+e(j)*kd(i,k)-2*e(i)*e(j)*e(k))
        hlubtensor(1,1,i,j,k)=y11h*((per(i,k,1)*e(1)+per(i,k,2)*e(2)+&
        &per(i,k,3)*e(3))*e(j) + (per(j,k,1)*e(1)+per(j,k,2)*e(2)+&
        &per(j,k,3)*e(3))*e(i))
        hlubtensor(1,2,i,j,k)=y12h*((per(i,k,1)*e(1)+per(i,k,2)*e(2)+&
        &per(i,k,3)*e(3))*e(j) + (per(j,k,1)*e(1)+per(j,k,2)*e(2)+&
        &per(j,k,3)*e(3))*e(i))
        end forall 
        
        forall(i=1:3,j=1:3,k=1:3,l=1:3)
        mlubtensor(i,j,k,l)=xm*(3.0_cp*(e(i)*e(j)-kd(i,j)/3.0_cp)*(e(k)*e(l)&
        & -kd(k,l)/3.0_cp)/2.0_cp)+ym*(e(i)*kd(j,l)*e(k)+e(j)*kd(i,l)*e(k)&
        & +e(i)*kd(j,k)*e(l)+e(j)*kd(i,k)*e(l)-4.0_cp*e(i)*e(j)*e(k)*e(l))&
        & /2.0_cp+zm*(kd(i,k)*kd(j,l)+kd(j,k)*kd(i,l)-kd(i,j)*kd(k,l)& 
        & +e(i)*e(j)*kd(k,l)+kd(i,j)*e(k)*e(l)-e(i)*kd(j,l)*e(k)&
        & -e(j)*kd(i,l)*e(k)-e(i)*kd(j,k)*e(l)-e(j)*kd(i,k)*e(l)&
        & +e(i)*e(j)*e(k)*e(l))/2.0_cp
        end forall
      
        forall(i=1:3,j=1:3)
        alub(2,2,i,j)=alub(1,1,i,j) !x22a=x11a,y22a=y11a
        alub(2,1,i,j)=alub(1,2,i,j) !x12a=x21a,y12a=y21a
        blub(2,2,i,j)=-blub(1,1,i,j) !y11b=-y22b
        blub(2,1,i,j)=-blub(1,2,i,j) !y12b=-y21b
        club(2,2,i,j)=club(1,1,i,j) !x11c=x22c
        club(2,1,i,j)=club(1,2,i,j) !x21c=x12c
        end forall

        forall(i=1:3,j=1:3,k=1:3)
        glubtensor(2,2,i,j,k)=-glubtensor(1,1,i,j,k) ! x22g=-x11g, y22g=-y11g
        glubtensor(2,1,i,j,k)=-glubtensor(1,2,i,j,k) ! x21g=-x12g, y21g=-y12g
        hlubtensor(2,2,i,j,k)=hlubtensor(1,1,i,j,k) !  y22h=y11h 
        hlubtensor(2,1,i,j,k)=hlubtensor(1,2,i,j,k) ! y21h=y12h
        end forall
        
        forall(i=1:3,j=1:3)
                ! fill in abc
                lub2bmx(i,j)=alub(1,1,i,j)
                lub2bmx(i+3,j)=alub(2,1,i,j)
                lub2bmx(i,j+3)=alub(1,2,i,j)
                lub2bmx(i+3,j+3)=alub(2,2,i,j)
                lub2bmx(6+i,6+j)=club(1,1,i,j)
                lub2bmx(9+i,6+j)=club(2,1,i,j)
                lub2bmx(6+i,9+j)=club(1,2,i,j)
                lub2bmx(9+i,9+j)=club(2,2,i,j)
                lub2bmx(6+i,j)=blub(1,1,i,j)
                lub2bmx(6+i,3+j)=blub(1,2,i,j)
                lub2bmx(9+i,j)=blub(2,1,i,j)
                lub2bmx(9+i,3+j)=blub(2,2,i,j)
                ! fill in btilda by symmetry 
                lub2bmx(j,6+i)=blub(1,1,i,j)
                lub2bmx(3+j,6+i)=blub(1,2,i,j)
                lub2bmx(j,9+i)=blub(2,1,i,j)
                lub2bmx(3+j,9+i)=blub(2,2,i,j)
        end forall

        forall(j=1:3)
                ! fill in g & h
                ! ev1=e11-e33 ev2=2e12 ev3=2e13 ev4=2e23 ev5=e22-e33
                ! sv1=s11 sv2=s12=s21 sv3=s13=s31 sv4=s23=s32 sv5=s22
                lub2bmx(13,j)=glubtensor(1,1,1,1,j)
                lub2bmx(14,j)=glubtensor(1,1,1,2,j)
                lub2bmx(15,j)=glubtensor(1,1,1,3,j)
                lub2bmx(16,j)=glubtensor(1,1,2,3,j)
                lub2bmx(17,j)=glubtensor(1,1,2,2,j)
                lub2bmx(18,j)=glubtensor(2,1,1,1,j)
                lub2bmx(19,j)=glubtensor(2,1,1,2,j)
                lub2bmx(20,j)=glubtensor(2,1,1,3,j)
                lub2bmx(21,j)=glubtensor(2,1,2,3,j)
                lub2bmx(22,j)=glubtensor(2,1,2,2,j)
                lub2bmx(13,3+j)=glubtensor(1,2,1,1,j)
                lub2bmx(14,3+j)=glubtensor(1,2,1,2,j)
                lub2bmx(15,3+j)=glubtensor(1,2,1,3,j)
                lub2bmx(16,3+j)=glubtensor(1,2,2,3,j)
                lub2bmx(17,3+j)=glubtensor(1,2,2,2,j)
                lub2bmx(18,3+j)=glubtensor(2,2,1,1,j)
                lub2bmx(19,3+j)=glubtensor(2,2,1,2,j)
                lub2bmx(20,3+j)=glubtensor(2,2,1,3,j)
                lub2bmx(21,3+j)=glubtensor(2,2,2,3,j)
                lub2bmx(22,3+j)=glubtensor(2,2,2,2,j)
                lub2bmx(j,13)=glubtensor(1,1,1,1,j)
                lub2bmx(j,14)=glubtensor(1,1,1,2,j)
                lub2bmx(j,15)=glubtensor(1,1,1,3,j)
                lub2bmx(j,16)=glubtensor(1,1,2,3,j)
                lub2bmx(j,17)=glubtensor(1,1,2,2,j)
                lub2bmx(j,18)=glubtensor(2,1,1,1,j)
                lub2bmx(j,19)=glubtensor(2,1,1,2,j)
                lub2bmx(j,20)=glubtensor(2,1,1,3,j)
                lub2bmx(j,21)=glubtensor(2,1,2,3,j)
                lub2bmx(j,22)=glubtensor(2,1,2,2,j)
                lub2bmx(3+j,13)=glubtensor(1,2,1,1,j)
                lub2bmx(3+j,14)=glubtensor(1,2,1,2,j)
                lub2bmx(3+j,15)=glubtensor(1,2,1,3,j)
                lub2bmx(3+j,16)=glubtensor(1,2,2,3,j)
                lub2bmx(3+j,17)=glubtensor(1,2,2,2,j)
                lub2bmx(3+j,18)=glubtensor(2,2,1,1,j)
                lub2bmx(3+j,19)=glubtensor(2,2,1,2,j)
                lub2bmx(3+j,20)=glubtensor(2,2,1,3,j)
                lub2bmx(3+j,21)=glubtensor(2,2,2,3,j)
                lub2bmx(3+j,22)=glubtensor(2,2,2,2,j)

                lub2bmx(13,6+j)=hlubtensor(1,1,1,1,j)
                lub2bmx(14,6+j)=hlubtensor(1,1,1,2,j)
                lub2bmx(15,6+j)=hlubtensor(1,1,1,3,j)
                lub2bmx(16,6+j)=hlubtensor(1,1,2,3,j)
                lub2bmx(17,6+j)=hlubtensor(1,1,2,2,j)
                lub2bmx(18,6+j)=hlubtensor(2,1,1,1,j)
                lub2bmx(19,6+j)=hlubtensor(2,1,1,2,j)
                lub2bmx(20,6+j)=hlubtensor(2,1,1,3,j)
                lub2bmx(21,6+j)=hlubtensor(2,1,2,3,j)
                lub2bmx(22,6+j)=hlubtensor(2,1,2,2,j)
                lub2bmx(13,9+j)=hlubtensor(1,2,1,1,j)
                lub2bmx(14,9+j)=hlubtensor(1,2,1,2,j)
                lub2bmx(15,9+j)=hlubtensor(1,2,1,3,j)
                lub2bmx(16,9+j)=hlubtensor(1,2,2,3,j)
                lub2bmx(17,9+j)=hlubtensor(1,2,2,2,j)
                lub2bmx(18,9+j)=hlubtensor(2,2,1,1,j)
                lub2bmx(19,9+j)=hlubtensor(2,2,1,2,j)
                lub2bmx(20,9+j)=hlubtensor(2,2,1,3,j)
                lub2bmx(21,9+j)=hlubtensor(2,2,2,3,j)
                lub2bmx(22,9+j)=hlubtensor(2,2,2,2,j)
                lub2bmx(6+j,13)=hlubtensor(1,1,1,1,j)
                lub2bmx(6+j,14)=hlubtensor(1,1,1,2,j)
                lub2bmx(6+j,15)=hlubtensor(1,1,1,3,j)
                lub2bmx(6+j,16)=hlubtensor(1,1,2,3,j)
                lub2bmx(6+j,17)=hlubtensor(1,1,2,2,j)
                lub2bmx(6+j,18)=hlubtensor(2,1,1,1,j)
                lub2bmx(6+j,19)=hlubtensor(2,1,1,2,j)
                lub2bmx(6+j,20)=hlubtensor(2,1,1,3,j)
                lub2bmx(6+j,21)=hlubtensor(2,1,2,3,j)
                lub2bmx(6+j,22)=hlubtensor(2,1,2,2,j)
                lub2bmx(9+j,13)=hlubtensor(1,2,1,1,j)
                lub2bmx(9+j,14)=hlubtensor(1,2,1,2,j)
                lub2bmx(9+j,15)=hlubtensor(1,2,1,3,j)
                lub2bmx(9+j,16)=hlubtensor(1,2,2,3,j)
                lub2bmx(9+j,17)=hlubtensor(1,2,2,2,j)
                lub2bmx(9+j,18)=hlubtensor(2,2,1,1,j)
                lub2bmx(9+j,19)=hlubtensor(2,2,1,2,j)
                lub2bmx(9+j,20)=hlubtensor(2,2,1,3,j)
                lub2bmx(9+j,21)=hlubtensor(2,2,2,3,j)
                lub2bmx(9+j,22)=hlubtensor(2,2,2,2,j)

                ! the rewritted matrix is symmetric
                
        end forall

                !fill in m11
                mmx(1,1)=(2.0_cp*mlubtensor(1,1,1,1) & 
                & -mlubtensor(1,1,2,2)-mlubtensor(1,1,3,3))/3.0_cp
                mmx(1,2)=(mlubtensor(1,1,1,2)+mlubtensor(1,1,2,1))/2.0_cp
                mmx(1,3)=(mlubtensor(1,1,1,3)+mlubtensor(1,1,3,1))/2.0_cp
                mmx(1,4)=(mlubtensor(1,1,2,3)+mlubtensor(1,1,3,2))/2.0_cp
                mmx(1,5)=(2.0_cp*mlubtensor(1,1,2,2) &
                & -mlubtensor(1,1,1,1)-mlubtensor(1,1,3,3))/3.0_cp

                mmx(2,1)=(2.0_cp*mlubtensor(1,2,1,1) & 
                & -mlubtensor(1,2,2,2)-mlubtensor(1,2,3,3))/3.0_cp
                mmx(2,2)=(mlubtensor(1,2,1,2)+mlubtensor(1,2,2,1))/2.0_cp
                mmx(2,3)=(mlubtensor(1,2,1,3)+mlubtensor(1,2,3,1))/2.0_cp
                mmx(2,4)=(mlubtensor(1,2,2,3)+mlubtensor(1,2,3,2))/2.0_cp
                mmx(2,5)=(2.0_cp*mlubtensor(1,2,2,2) &
                & -mlubtensor(1,2,1,1)-mlubtensor(1,2,3,3))/3.0_cp
                
                mmx(3,1)=(2.0_cp*mlubtensor(1,3,1,1) & 
                & -mlubtensor(1,3,2,2)-mlubtensor(1,3,3,3))/3.0_cp
                mmx(3,2)=(mlubtensor(1,3,1,2)+mlubtensor(1,3,2,1))/2.0_cp
                mmx(3,3)=(mlubtensor(1,3,1,3)+mlubtensor(1,3,3,1))/2.0_cp
                mmx(3,4)=(mlubtensor(1,3,2,3)+mlubtensor(1,3,3,2))/2.0_cp
                mmx(3,5)=(2.0_cp*mlubtensor(1,3,2,2) &
                & -mlubtensor(1,3,1,1)-mlubtensor(1,3,3,3))/3.0_cp
                
                mmx(4,1)=(2.0_cp*mlubtensor(2,3,1,1) & 
                & -mlubtensor(2,3,2,2)-mlubtensor(2,3,3,3))/3.0_cp
                mmx(4,2)=(mlubtensor(2,3,1,2)+mlubtensor(2,3,2,1))/2.0_cp
                mmx(4,3)=(mlubtensor(2,3,1,3)+mlubtensor(2,3,3,1))/2.0_cp
                mmx(4,4)=(mlubtensor(2,3,2,3)+mlubtensor(2,3,3,2))/2.0_cp
                mmx(4,5)=(2.0_cp*mlubtensor(2,3,2,2) &
                & -mlubtensor(2,3,1,1)-mlubtensor(2,3,3,3))/3.0_cp

                mmx(5,1)=(2.0_cp*mlubtensor(2,2,1,1) & 
                & -mlubtensor(2,2,2,2)-mlubtensor(2,2,3,3))/3.0_cp
                mmx(5,2)=(mlubtensor(2,2,1,2)+mlubtensor(2,2,2,1))/2.0_cp
                mmx(5,3)=(mlubtensor(2,2,1,3)+mlubtensor(2,2,3,1))/2.0_cp
                mmx(5,4)=(mlubtensor(2,2,2,3)+mlubtensor(2,2,3,2))/2.0_cp
                mmx(5,5)=(2.0_cp*mlubtensor(2,2,2,2) &
                & -mlubtensor(2,2,1,1)-mlubtensor(2,2,3,3))/3.0_cp
               
                !replicate m11 to m12,m21,m22
                forall(i=1:5,j=1:5)
                        lub2bmx(12+i,12+j)=mmx(i,j)
                        lub2bmx(17+i,12+j)=mmx(i,j)
                        lub2bmx(12+i,17+j)=mmx(i,j)
                        lub2bmx(17+i,17+j)=mmx(i,j)
                end forall
        end if

        ! add the lub2bmx to rfu , rsu , rse
        ! step1: add a,b,btilda,c to rfu
        forall(i=1:3,j=1:3)
        ! a
        rfu(3*(alpha-1)+i,3*(alpha-1)+j)=rfu(3*(alpha-1)+i,3*(alpha-1)+j)&
        & +lub2bmx(i,j) ! a alpha alpha
        rfu(3*(alpha-1)+i,3*(beta-1)+j)=rfu(3*(alpha-1)+i,3*(beta-1)+j)&
        & +lub2bmx(i,3+j) ! a alpha beta
        rfu(3*(beta-1)+i,3*(alpha-1)+j)=rfu(3*(beta-1)+i,3*(alpha-1)+j)&
        & +lub2bmx(3+i,j) ! a beta alpha
        rfu(3*(beta-1)+i,3*(beta-1)+j)=rfu(3*(beta-1)+i,3*(beta-1)+j)&
        & +lub2bmx(3+i,3+j) ! a beta beta
        ! b
        rfu(3*ntotal+3*(alpha-1)+i,3*(alpha-1)+j)=rfu(3*ntotal+3*(alpha-1)&
        & +i,3*(alpha-1)+j)+lub2bmx(6+i,j) !balphaalpha
        rfu(3*ntotal+3*(alpha-1)+i,3*(beta-1)+j)=rfu(3*ntotal+3*(alpha-1)&
        & +i,3*(beta-1)+j)+lub2bmx(6+i,3+j) !balphabeta
        rfu(3*ntotal+3*(beta-1)+i,3*(alpha-1)+j)=rfu(3*ntotal+3*(beta-1)&
        & +i,3*(alpha-1)+j)+lub2bmx(9+i,j) !b beta alpha
        rfu(3*ntotal+3*(beta-1)+i,3*(beta-1)+j)=rfu(3*ntotal+3*(beta-1)&
        & +i,3*(beta-1)+j)+lub2bmx(9+i,3+j) !b beta beta
        ! btilda by symmetry
        rfu(3*(alpha-1)+j,3*ntotal+3*(alpha-1)+i)=rfu(3*(alpha-1)+j,&
        & 3*ntotal+3*(alpha-1)+i)+lub2bmx(6+i,j) !balphaalpha
        rfu(3*(beta-1)+j,3*ntotal+3*(alpha-1)+i)=rfu(3*(beta-1)+j,&
        & 3*ntotal+3*(alpha-1)+i)+lub2bmx(6+i,3+j) !balphabeta
        rfu(3*(alpha-1)+j,3*ntotal+3*(beta-1)+i)=rfu(3*(alpha-1)+j,&
        & 3*ntotal+3*(beta-1)+i)+lub2bmx(9+i,j) !b beta alpha
        rfu(3*(beta-1)+j,3*ntotal+3*(beta-1)+i)=rfu(3*(beta-1)+j,&
        & 3*ntotal+3*(beta-1)+i)+lub2bmx(9+i,3+j) !b beta beta
        ! c 
        rfu(3*ntotal+3*(alpha-1)+i,3*ntotal+3*(alpha-1)+j)=rfu(3*ntotal&
        & +3*(alpha-1)+i,3*ntotal+3*(alpha-1)+j)+lub2bmx(6+i,6+j) !c a a
        rfu(3*ntotal+3*(alpha-1)+i,3*ntotal+3*(beta-1)+j)=rfu(3*ntotal&
        & +3*(alpha-1)+i,3*ntotal+3*(beta-1)+j)+lub2bmx(6+i,9+j) !c a b
        rfu(3*ntotal+3*(beta-1)+i,3*ntotal+3*(alpha-1)+j)=rfu(3*ntotal&
        & +3*(beta-1)+i,3*ntotal+3*(alpha-1)+j)+lub2bmx(9+i,6+j) !c b a
        rfu(3*ntotal+3*(beta-1)+i,3*ntotal+3*(beta-1)+j)=rfu(3*ntotal&
        & +3*(beta-1)+i,3*ntotal+3*(beta-1)+j)+lub2bmx(9+i,9+j) !c b b
        end forall 
        
        !step2: add g & h to rsu
        forall(i=1:5,j=1:3)
        ! g
        rsu(5*(alpha-1)+i,3*(alpha-1)+j)=rsu(5*(alpha-1)+i,3*(alpha-1)+j)&
        & +lub2bmx(12+i,j) ! g a a
        rsu(5*(alpha-1)+i,3*(beta-1)+j)=rsu(5*(alpha-1)+i,3*(beta-1)+j)&
        & +lub2bmx(12+i,3+j) ! g a b
        rsu(5*(beta-1)+i,3*(alpha-1)+j)=rsu(5*(beta-1)+i,3*(alpha-1)+j)&
        & +lub2bmx(17+i,j) ! g b a
        rsu(5*(beta-1)+i,3*(beta-1)+j)=rsu(5*(beta-1)+i,3*(beta-1)+j)&
        & +lub2bmx(17+i,3+j) ! g b b
        ! h
        rsu(5*(alpha-1)+i,3*(alpha-1)+j+3*ntotal)=&
        & rsu(5*(alpha-1)+i,3*(alpha-1)+j+3*ntotal)+lub2bmx(12+i,6+j) ! h a a
        rsu(5*(alpha-1)+i,3*(beta-1)+j+3*ntotal)=&
        & rsu(5*(alpha-1)+i,3*(beta-1)+j+3*ntotal)+lub2bmx(12+i,9+j) ! h a b
        rsu(5*(beta-1)+i,3*(alpha-1)+j+3*ntotal)=&
        & rsu(5*(beta-1)+i,3*(alpha-1)+j+3*ntotal)+lub2bmx(17+i,6+j) ! h b a
        rsu(5*(beta-1)+i,3*(beta-1)+j+3*ntotal)=&
        & rsu(5*(beta-1)+i,3*(beta-1)+j+3*ntotal)+lub2bmx(17+i,9+j) ! h b b
        end forall

        !step3: add m to rse
        forall(i=1:5,j=1:5)
        rse(5*(alpha-1)+i,5*(alpha-1)+j)=rse(5*(alpha-1)+i,5*(alpha-1)+j)&
        & +lub2bmx(12+i,12+j) ! m a a
        rse(5*(alpha-1)+i,5*(beta-1)+j)=rse(5*(alpha-1)+i,5*(beta-1)+j)&
        & +lub2bmx(12+i,17+j) ! m a b
        rse(5*(beta-1)+i,5*(alpha-1)+j)=rse(5*(beta-1)+i,5*(alpha-1)+j)&
        & +lub2bmx(17+i,12+j) ! m b a
        rse(5*(beta-1)+i,5*(beta-1)+j)=rse(5*(beta-1)+i,5*(beta-1)+j)&
        & +lub2bmx(17+i,17+j) ! m b b
        end forall
        
        !step4: fill rfe by symmetry with rsu
        forall(i=1:6*ntotal,j=1:5*ntotal)
        rfe(i,j)=rsu(j,i)
        endforall
        
        ! output rfu rfe rsu rse for debugging
        !open(unit=2001,file='rfu',status='replace')
        !write(*,*) rfu
        !close(unit=2001)
        !open(unit=2002,file='rfe',status='replace')
        !write(*,*) rfe
        !close(unit=2002)
        !open(unit=2003,file='rsu',status='replace')
        !write(*,*) rsu
        !close(unit=2003)
        !open(unit=2004,file='rse',status='replace')
        !write(*,*) rse
        !close(unit=2004)

        end do
end do

end subroutine lubmxcalc

subroutine lubcalc(r,x11a,x12a,y11a,y12a,y11b,y12b,x11c,x12c,y11c,y12c,&
 & x11g,x12g,y11g,y12g,y11h,y12h,xm,ym,zm)
implicit none
real(kind=cp), intent(in) :: r
real(kind=cp), intent(out) :: x11a
real(kind=cp), intent(out) :: x12a
real(kind=cp), intent(out) :: y11a
real(kind=cp), intent(out) :: y12a
real(kind=cp), intent(out) :: y11b
real(kind=cp), intent(out) :: y12b
real(kind=cp), intent(out) :: x11c
real(kind=cp), intent(out) :: x12c
real(kind=cp), intent(out) :: y11c
real(kind=cp), intent(out) :: y12c
real(kind=cp), intent(out) :: x11g
real(kind=cp), intent(out) :: x12g
real(kind=cp), intent(out) :: y11g
real(kind=cp), intent(out) :: y12g
real(kind=cp), intent(out) :: y11h
real(kind=cp), intent(out) :: y12h
real(kind=cp), intent(out) :: xm
real(kind=cp), intent(out) :: ym
real(kind=cp), intent(out) :: zm

real(kind=cp) :: xi
real(kind=cp) :: invxi
real(kind=cp) :: loginvxi
real(kind=cp) :: xiloginvxi

if (r.lt.2.1_cp) then
        xi=r-2.0_cp
        invxi=1.0_cp/xi
        loginvxi=log(invxi)
        xiloginvxi=xi*loginvxi
        
        x11a=-1.23041_cp+0.25_cp*invxi+1.8918_cp*xi+9.0_cp*loginvxi/40.0_cp & 
        &+3.0_cp*xiloginvxi/112.0_cp
        x12a=-x11a+0.00312_cp-0.0011_cp*xi

        y11a=-0.39394_cp+0.95665_cp*xi+loginvxi/6.0_cp
        y12a=-y11a+0.004636_cp-0.007049_cp*xi

        y11b=0.408286_cp-0.84055_cp*xi-loginvxi/6.0_cp-xiloginvxi/12.0_cp

        x11c=0.0479_cp+0.12494_cp*xi-xiloginvxi/6.0_cp
        x12c=-x11c+0.016869_cp+0.049536_cp*xi

        y11c=-0.605434_cp+0.939139_cp*xi+4.0_cp*loginvxi/15.0_cp & 
        &+94.0_cp*xiloginvxi/375.0_cp
        !y12c corrected, compare to the old code. 62/375
        !y12c tabulated data is the same
        y12c=loginvxi/15.0_cp+62.0_cp*xiloginvxi/375.0_cp-0.219032_cp&
        &+0.332843_cp*xi
        y12b=-0.405978_cp+0.833042_cp*xi+loginvxi/6.0_cp+xiloginvxi/12.0_cp

        x11g=-1.16897_cp+1.47882_cp*xi+invxi/4.0_cp & 
        &+9.0_cp*loginvxi/40.0_cp+39.0_cp*xiloginvxi/280.0_cp
        x12g=-x11g+0.01_cp-0.00167_cp*xi

        y11g=-0.2041_cp+0.44226_cp*xi+loginvxi/12.0_cp+xiloginvxi/24.0_cp
        y12g=-y11g+0.012265_cp-0.02757_cp*xi

        y11h=-0.143777_cp+0.264207_cp*xi+loginvxi/30.0_cp & 
        &+137.0_cp*xiloginvxi/1500.0_cp
        y12h=-0.298166_cp+0.534123_cp*xi+2.0_cp*loginvxi/15.0_cp &
        &+113.0_cp*xiloginvxi/1500.0_cp
        
        ! pay attention to m terms. here, in the equation,
        ! xm=(x11m+x12m)/2, ym=(y11m+y12m)/2, zm=(z11m+z12m)/2
        ! however, in the tabulated data, the data are 
        ! (x11m+x12m), (y11m+y12m), (z11m+z12m)
        ! so, in the tabcalc function, 
        ! xm ym zm tab values should be divided by 2
        ! this gives the same results as all the old codes

        xm=invxi/6.0_cp+3.0_cp*loginvxi/20.0_cp+47.0_cp*xiloginvxi/280.0_cp&
        & -0.740815_cp+0.706802_cp*xi
        ! the ym equation contains an extra term (xiloginvxi)
        ! compared to the old code, so the O(1), O(xi) terms are modified
        ym=-0.193518_cp+0.553883_cp*xi+loginvxi/12.0_cp&
        & -29.0_cp*xiloginvxi/500.0_cp
        ! the zm equation & tabulation are not modified
        zm = 0.00645755_cp - 0.021142_cp * xi
else
        call lubtabcalc(r,x11a,x12a,y11a,y12a,y11b,y12b,x11c,x12c,y11c,y12c,&
 & x11g,x12g,y11g,y12g,y11h,y12h,xm,ym,zm)
end if

end subroutine lubcalc

subroutine lubtabcalc(r,x11a,x12a,y11a,y12a,y11b,y12b,x11c,x12c,y11c,y12c,&
 & x11g,x12g,y11g,y12g,y11h,y12h,xm,ym,zm)
implicit none
 
real(kind=cp), intent(in) :: r
real(kind=cp), intent(out) :: x11a
real(kind=cp), intent(out) :: x12a
real(kind=cp), intent(out) :: y11a
real(kind=cp), intent(out) :: y12a
real(kind=cp), intent(out) :: y11b
real(kind=cp), intent(out) :: y12b
real(kind=cp), intent(out) :: x11c
real(kind=cp), intent(out) :: x12c
real(kind=cp), intent(out) :: y11c
real(kind=cp), intent(out) :: y12c
real(kind=cp), intent(out) :: x11g
real(kind=cp), intent(out) :: x12g
real(kind=cp), intent(out) :: y11g
real(kind=cp), intent(out) :: y12g
real(kind=cp), intent(out) :: y11h
real(kind=cp), intent(out) :: y12h
real(kind=cp), intent(out) :: xm
real(kind=cp), intent(out) :: ym
real(kind=cp), intent(out) :: zm

integer, parameter :: ntabc=39
real(kind=cp),dimension(ntabc) , save :: trabc
real(kind=cp),dimension(ntabc) , save:: tx11a
real(kind=cp),dimension(ntabc) , save:: tx12a
real(kind=cp),dimension(ntabc) , save:: ty11a
real(kind=cp),dimension(ntabc) , save:: ty12a
real(kind=cp),dimension(ntabc) , save:: ty11b
real(kind=cp),dimension(ntabc) , save:: ty12b
real(kind=cp),dimension(ntabc) , save:: tx11c
real(kind=cp),dimension(ntabc) , save:: tx12c
real(kind=cp),dimension(ntabc) , save:: ty11c
real(kind=cp),dimension(ntabc) , save:: ty12c

integer, parameter :: ntgh=47
real(kind=cp),dimension(ntgh) , save:: trgh
real(kind=cp),dimension(ntgh) , save:: tx11g
real(kind=cp),dimension(ntgh) , save:: tx12g
real(kind=cp),dimension(ntgh) , save:: ty11g
real(kind=cp),dimension(ntgh) , save:: ty12g
real(kind=cp),dimension(ntgh) , save:: ty11h
real(kind=cp),dimension(ntgh) , save:: ty12h

integer, parameter :: ntm=47
real(kind=cp),dimension(ntm) , save:: trm
real(kind=cp),dimension(ntm) , save:: txm
real(kind=cp),dimension(ntm) , save:: tym
real(kind=cp),dimension(ntm) , save:: tzm

integer :: i

integer, save :: readflag=0

if (readflag.eq.0) then

call readlubabc(ntabc,trabc,tx11a,tx12a,ty11a,ty12a,ty11b,ty12b,&
                &tx11c,tx12c,ty11c,ty12c)
call readlubgh(ntgh,trgh,tx11g,tx12g,ty11g,ty12g,ty11h,ty12h)
call readlubm(ntm,trm,txm,tym,tzm)

readflag = 1
end if

x11a=0
x12a=0
y11a=0
y12a=0
y11b=0
y12b=0
x11c=0
x12c=0
y11c=0
y12c=0
x11g=0
x12g=0
y11g=0
y12g=0
y11h=0
y12h=0
xm=0
ym=0
zm=0

do i=1,ntabc-1
        if((trabc(i).le.r).and.(trabc(i+1).gt.r)) then
        x11a=tx11a(i)+(r-trabc(i))*(tx11a(i+1)-tx11a(i))/(trabc(i+1)-trabc(i))
        x12a=tx12a(i)+(r-trabc(i))*(tx12a(i+1)-tx12a(i))/(trabc(i+1)-trabc(i))
        y11a=ty11a(i)+(r-trabc(i))*(ty11a(i+1)-ty11a(i))/(trabc(i+1)-trabc(i))
        y12a=ty12a(i)+(r-trabc(i))*(ty12a(i+1)-ty12a(i))/(trabc(i+1)-trabc(i))
        y11b=ty11b(i)+(r-trabc(i))*(ty11b(i+1)-ty11b(i))/(trabc(i+1)-trabc(i))
        y12b=ty12b(i)+(r-trabc(i))*(ty12b(i+1)-ty12b(i))/(trabc(i+1)-trabc(i))
        x11c=tx11c(i)+(r-trabc(i))*(tx11c(i+1)-tx11c(i))/(trabc(i+1)-trabc(i))
        x12c=tx12c(i)+(r-trabc(i))*(tx12c(i+1)-tx12c(i))/(trabc(i+1)-trabc(i))
        y11c=ty11c(i)+(r-trabc(i))*(ty11c(i+1)-ty11c(i))/(trabc(i+1)-trabc(i))
        y12c=ty12c(i)+(r-trabc(i))*(ty12c(i+1)-ty12c(i))/(trabc(i+1)-trabc(i))
        exit
        end if
end do

do i=1,ntgh-1
        if((trgh(i).le.r).and.(trgh(i+1).gt.r)) then
        x11g=tx11g(i)+(r-trgh(i))*(tx11g(i+1)-tx11g(i))/(trgh(i+1)-trgh(i))
        x12g=tx12g(i)+(r-trgh(i))*(tx12g(i+1)-tx12g(i))/(trgh(i+1)-trgh(i))
        y11g=ty11g(i)+(r-trgh(i))*(ty11g(i+1)-ty11g(i))/(trgh(i+1)-trgh(i))
        y12g=ty12g(i)+(r-trgh(i))*(ty12g(i+1)-ty12g(i))/(trgh(i+1)-trgh(i))
        y11h=ty11h(i)+(r-trgh(i))*(ty11h(i+1)-ty11h(i))/(trgh(i+1)-trgh(i))
        y12h=ty12h(i)+(r-trgh(i))*(ty12h(i+1)-ty12h(i))/(trgh(i+1)-trgh(i))
        exit
        end if
end do

do i=1,ntm-1
        if((trm(i).le.r).and.(trm(i+1).gt.r)) then
        xm=(txm(i)+(r-trm(i))*(txm(i+1)-txm(i))/(trm(i+1)-trm(i)))/2.0_cp
        ym=(tym(i)+(r-trm(i))*(tym(i+1)-tym(i))/(trm(i+1)-trm(i)))/2.0_cp
        zm=(tzm(i)+(r-trm(i))*(tzm(i+1)-tzm(i))/(trm(i+1)-trm(i)))/2.0_cp
        exit
        end if
end do

end subroutine lubtabcalc


end module lub
