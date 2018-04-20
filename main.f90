program sd
use f95_precision  ! provide the dp and sp kind of real
use mob ! subroutine grmobmxcalc
use lub ! subroutine lubcalc
use fts ! subroutine ftscalc
use prutil  ! select precision
use in_out ! input and output

implicit none
! particle pointer
integer :: i=0
integer :: j=0
integer :: ntotal=0 ! total number of particles

! the properties of the background flow
real(kind=cp), dimension(3) :: u_bg=(/0,0,0/)
real(kind=cp), dimension(3) :: omega_bg=(/0,0,0/)
real(kind=cp), dimension(3,3) :: einf_bg

! the properties of each particle
real(kind=cp), allocatable, dimension(:,:) :: p_pos             ! n*3
real(kind=cp), allocatable, dimension(:,:) :: p_pos_local       ! n*3
real(kind=cp), allocatable, dimension(:,:) :: p_u             ! n*3
real(kind=cp), allocatable, dimension(:,:) :: p_o           ! n*3
! ev1=e11-e33 ev2=2e12 ev3=2e13 ev4=2e23 ev5=e22-e33
real(kind=cp), allocatable, dimension(:,:) :: p_fe       ! n*3
real(kind=cp), allocatable, dimension(:,:) :: p_te      ! n*3
real(kind=cp), allocatable, dimension(:,:) :: p_sh       ! n*5
! sv1=s11 sv2=s12=s21 sv3=s13=s31 sv4=s23=s32 sv5=s22

! the grand mobility, added lubrication and grand resistance matrix
real(kind=SP), allocatable, dimension(:,:) :: grmobmx           !11n*11n
real(kind=SP), allocatable, dimension(:,:) :: rfu             !6n*6n
real(kind=SP), allocatable, dimension(:,:) :: rfe             !6n*5n
real(kind=SP), allocatable, dimension(:,:) :: rsu             !5n*6n
real(kind=SP), allocatable, dimension(:,:) :: rse             !5n*5n
! the time iteration info
real(kind=cp) :: timestepsize=0.0_cp
real(kind=cp) :: totaltime=0.0_cp
real(kind=cp) :: time=0.0_cp
integer :: intervalcount=0
integer :: interval=0

! other stuff
character(len=80) :: line
character(len=8) :: pname
!non-dimensionalization process:

! initialize

! read bgflow\time\timestepsize
call readconf(timestepsize,totaltime,interval,u_bg,omega_bg,einf_bg)
intervalcount=interval

! read the particle initial position and velocity
open(unit=3,file='particles.in',status='old')
read(3,'(i5)') ntotal
close(unit=3)
! allocate the particle properties
allocate(p_pos(ntotal,3))
allocate(p_pos_local(ntotal,3))
allocate(p_u(ntotal,3))
allocate(p_o(ntotal,3))
allocate(p_fe(ntotal,3))
allocate(p_te(ntotal,3))
allocate(p_sh(ntotal,5))
! initialize the external force & torque
p_te=0.0_cp
p_fe=0.0_cp
p_sh=0.0_cp

call readparticlein(ntotal,p_pos,p_fe)

! allocate & initialize the matrix
allocate(grmobmx(11*ntotal,11*ntotal))
allocate(rfu(6*ntotal,6*ntotal))
allocate(rfe(6*ntotal,5*ntotal))
allocate(rsu(5*ntotal,6*ntotal))
allocate(rse(5*ntotal,5*ntotal))
grmobmx=0.0_SP
rfu=0.0_SP
rfe=0.0_SP
rsu=0.0_SP
rse=0.0_SP

! end of initialization
 
! begin iteration

do
        ! calculate the farfield grand mobility matrix
        call grmobmxcalc(ntotal,p_pos,grmobmx)

        ! calculate the lub contribution to rfu,rfe,rsu,rse
        call lubmxcalc(ntotal,p_pos,rfu,rfe,rsu,rse)

        ! calculate the resistance sub matrices rfu,rfe,rsu,rse
        call remxcalc(ntotal,grmobmx,rfu,rfe,rsu,rse)

        ! update u & pos
        ! U = Rfu^-1(Fe+Rfe.EV)+Uinf
        ! pos = pos + U * timestepsize
        call hydrostep(ntotal,p_pos,p_u,p_o,p_fe,u_bg,omega_bg,einf_bg&
        & ,rfu,rfe,timestepsize)
        
        ! output the data
        if(intervalcount.eq.interval) then
        call writeparticledata(ntotal,time,p_pos,p_u,p_o,p_sh)
        call writeposxyz(ntotal, time, p_pos)
        intervalcount=0
        end if
        
        intervalcount=intervalcount+1
        
        ! stop condition for the whole simulation
        time = time+timestepsize
        if (time.gt.totaltime) then
                exit
        end if
end do


! end of execution, deallocate all the arrays
deallocate(p_pos)
deallocate(p_u)
deallocate(p_o)
deallocate(p_fe)
deallocate(p_te)
deallocate(p_sh)
deallocate(grmobmx)
deallocate(rfu)
deallocate(rfe)
deallocate(rsu)
deallocate(rse)

stop
end program sd
