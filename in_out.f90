module in_out
use prutil

contains
subroutine readconf(timestepsize,totaltime,interval,u_bg,omega_bg,einf_bg)
implicit none
real(kind=cp), intent(out) :: timestepsize
real(kind=cp), intent(out) :: totaltime
integer, intent(out) :: interval
real(kind=cp), dimension(3), intent(out) :: u_bg
real(kind=cp), dimension(3), intent(out) :: omega_bg
real(kind=cp), dimension(3,3), intent(out) :: einf_bg

character(len=80) :: line

open(unit=2,file='conf.in',status='old')
read(2,'(a80)') line
read(2,*) timestepsize
read(2,'(a80)') line
read(2,*) totaltime
read(2,'(a80)') line
read(2,*) interval
read(2,'(a80)') line
read(2,*) u_bg(1),u_bg(2),u_bg(3)
read(2,'(a80)') line
read(2,*) omega_bg(1),omega_bg(2),omega_bg(3)
read(2,'(a80)') line
read(2,*) einf_bg(1,1),einf_bg(1,2),einf_bg(1,3),einf_bg(2,2),einf_bg(2,3)
einf_bg(2,1)=einf_bg(1,2)
einf_bg(3,1)=einf_bg(1,3)
einf_bg(3,2)=einf_bg(2,3)
einf_bg(3,3)=0-einf_bg(1,1)-einf_bg(2,2)
close(unit=2)

end subroutine readconf

subroutine readparticlein(ntotal,p_pos,p_fe)
implicit none
integer, intent(in):: ntotal
real(kind=cp),dimension(ntotal,3), intent(out) :: p_pos
real(kind=cp),dimension(ntotal,3), intent(out) :: p_fe
character(len=8) :: pname
character(len=80) :: line
integer :: i,j

open(unit=3,file='particles.in',status='old')
read(3,'(a80)') line
read(3,'(a80)') line
! initialize
p_pos=0
p_fe=0
! read
do i=1,ntotal
read(3,'(A8,3ES23.16)') pname,p_pos(i,1),p_pos(i,2),p_pos(i,3)
read(3,'(A8,3ES23.16)') pname,p_fe(i,1),p_fe(i,2),p_fe(i,3)
read(3,'(a80)') line
end do

close(unit=3)
end subroutine readparticlein

subroutine writeposxyz_indiv(ntotal, time, p_pos)
implicit none
integer, intent(in) :: ntotal
real(kind=cp), intent(in) :: time
real(kind=cp), dimension(ntotal,3), intent(in) :: p_pos

character(len=8) :: pname='A  '
integer, save :: replaceflag=1
integer :: i,j
if(replaceflag.eq.1) then
        open(unit=101,file='out1.xyz',status='replace')
        open(unit=102,file='out2.xyz',status='replace')
        open(unit=103,file='out3.xyz',status='replace')
        write(101,*) 'time      x       y       z'
        write(102,*) 'time      x       y       z'
        write(103,*) 'time      x       y       z'
else
        open(unit=101,file='out1.xyz',status='old',position='append')
        open(unit=102,file='out2.xyz',status='old',position='append')
        open(unit=103,file='out3.xyz',status='old',position='append')
end if

replaceflag= 0 

write(101,'(4F11.4)') time, p_pos(1,1),p_pos(1,2),p_pos(1,3)
write(102,'(4F11.4)') time, p_pos(2,1),p_pos(2,2),p_pos(2,3)
write(103,'(4F11.4)') time, p_pos(3,1),p_pos(3,2),p_pos(3,3)

close(unit=101)
close(unit=102)
close(unit=103)

end subroutine writeposxyz_indiv



subroutine writeposxyz(ntotal, time, p_pos)
implicit none
integer, intent(in) :: ntotal
real(kind=cp), intent(in) :: time
real(kind=cp), dimension(ntotal,3), intent(in) :: p_pos

character(len=8) :: pname='A  '
integer, save :: replaceflag=1
integer :: i,j
if(replaceflag.eq.1) then
        open(unit=101,file='out.xyz',status='replace')
else
        open(unit=101,file='out.xyz',status='old',position='append')
end if

replaceflag= 0 

write(101,*) ntotal
write(101,*) time, 'x y z'
do i=1,ntotal
write(101,'(a8,3ES23.16)') pname, p_pos(i,1),p_pos(i,2),p_pos(i,3)
end do

close(unit=101)

end subroutine writeposxyz

subroutine writeparticledata(ntotal,time,p_pos,p_u,p_o,p_sh)
implicit none
integer, intent(in) :: ntotal
real(kind=cp), intent(in) :: time
real(kind=cp), dimension(ntotal,3), intent(in) :: p_pos
real(kind=cp), dimension(ntotal,3), intent(in) :: p_u
real(kind=cp), dimension(ntotal,3), intent(in) :: p_o
real(kind=cp), dimension(ntotal,5), intent(in) :: p_sh

integer :: i,j
integer, save :: replaceflag=1

if(replaceflag.eq.1) then
        open(unit=101,file='particles.out.data',status='replace')
else
        open(unit=101,file='particles.out.data',status='old',position='append')
end if
replaceflag = 0

write(101,*) ntotal
write(101,*) 'time:', time, 'x y z'

do i=1,ntotal
        write(101,*) i,'pos', p_pos(i,1),p_pos(i,2),p_pos(i,3)
        write(101,*) i,'vel', p_u(i,1),p_u(i,2),p_u(i,3)
        write(101,*) i,'omega', p_o(i,1),p_o(i,2),p_o(i,3)
        write(101,*) i,'hydro-stresslet',p_sh(i,1),p_sh(i,2),p_sh(i,3)
end do
write(101,*) '--------------------------------------------------------'
close(unit=101)

end subroutine writeparticledata

subroutine readlubabc(ntabc,trabc,tx11a,tx12a,ty11a,ty12a,ty11b,ty12b,&
&tx11c,tx12c,ty11c,ty12c)
implicit none
integer, intent(in) :: ntabc
real(kind=cp),dimension(ntabc),intent(out) :: trabc
real(kind=cp),dimension(ntabc),intent(out) :: tx11a
real(kind=cp),dimension(ntabc),intent(out) :: tx12a
real(kind=cp),dimension(ntabc),intent(out) :: ty11a
real(kind=cp),dimension(ntabc),intent(out) :: ty12a
real(kind=cp),dimension(ntabc),intent(out) :: ty11b
real(kind=cp),dimension(ntabc),intent(out) :: ty12b
real(kind=cp),dimension(ntabc),intent(out) :: tx11c
real(kind=cp),dimension(ntabc),intent(out) :: tx12c
real(kind=cp),dimension(ntabc),intent(out) :: ty11c
real(kind=cp),dimension(ntabc),intent(out) :: ty12c

integer :: i
character(len=80) :: line
open(unit = 1001,file = 'lubdat/r2babc.dat',status = 'old')
read(1001,'(a80)') line

do i = 1, ntabc
read( 1001, * ) trabc(i),tx11a(i),tx12a(i),ty11a(i),ty12a(i),&
        & ty11b(i),ty12b(i),tx11c(i),tx12c(i),ty11c(i),ty12c(i)
end do

close(unit=1001)
end subroutine readlubabc

subroutine readlubgh(ntgh,trgh,tx11g,tx12g,ty11g,ty12g,ty11h,ty12h)
implicit none
integer, intent(in) :: ntgh
real(kind=cp),dimension(ntgh),intent(out) ::trgh
real(kind=cp),dimension(ntgh),intent(out) ::tx11g
real(kind=cp),dimension(ntgh),intent(out) ::tx12g
real(kind=cp),dimension(ntgh),intent(out) ::ty11g
real(kind=cp),dimension(ntgh),intent(out) ::ty12g
real(kind=cp),dimension(ntgh),intent(out) ::ty11h
real(kind=cp),dimension(ntgh),intent(out) ::ty12h
integer :: i 
character(len=80) :: line
open( unit = 1002, file = 'lubdat/r2bgh.dat', status = 'old' )
read(1002, '(a80)') line
do i = 1, ntgh
read( 1002, * ) trgh(i),tx11g(i),tx12g(i),ty11g(i),ty12g(i),ty11h(i),ty12h(i)
end do

close(unit=1002)
end subroutine readlubgh

subroutine readlubm(ntm,trm,txm,tym,tzm)
implicit none
integer, intent(in) :: ntm
real(kind=cp),dimension(ntm),intent(out) ::trm
real(kind=cp),dimension(ntm),intent(out) ::txm
real(kind=cp),dimension(ntm),intent(out) ::tym
real(kind=cp),dimension(ntm),intent(out) ::tzm
integer :: i
character(len=80) :: line
open( unit = 1003, file = 'lubdat/r2bm.dat', status = 'old' )
read(1003,'(a80)') line
do i=1,ntm
read(1003,*) trm(i),txm(i),tym(i),tzm(i)
end do

close(unit=1003)

end subroutine readlubm

end module in_out 
