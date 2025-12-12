program rad_rho
implicit none
character*256 :: frho,fout,fctl
real*8 :: cc(3)
logical :: iqst
type rho_data
  integer*4 :: n1,n2,n3,nnodes_tmp
  real*8 :: AL(3,3)
  real*8,allocatable,dimension(:,:,:) :: rho
end type rho_data
type(rho_data) :: dat
real*8,allocatable :: lr(:), hr(:)
call getarg(1,frho)
call getarg(2,fctl)
write(fout,141) trim(frho),".withr"
141 format(A,A)

! default set 1/2 as the radial center
cc=[1,1,1]/2.d0

! ensure
write(6,*) "input: ", trim(frho)
write(6,*) "outpt: ", trim(fout)
write(6,*) "radcc: ", cc

iqst=check_file(frho)

call read_rho(frho,dat)

call getrad(dat,cc,fctl,lr,hr)

call write_rad(fout,lr,hr)

contains
subroutine write_rad(fout,lr,hr)
implicit none
integer*4 :: new_unit,nsiz,i
character(len=*),intent(in) :: fout
real*8, intent(in) :: lr(:),hr(:)
nsiz=size(lr)
open(newunit=new_unit,file=trim(fout))
rewind(new_unit)
write(new_unit,161) "r(ang)","rho(e/Bhor3)"
161 format(A24,4x,A24)
do i=1,nsiz
write(new_unit,151) lr(i),hr(i)
end do
151 format(E24.16,4x,E24.16)
close(new_unit)
end subroutine write_rad

subroutine getrad(dat,cc,fctl,lr,hr)
implicit none
type(rho_data),intent(in) :: dat
real*8,intent(in) :: cc(3)
character(len=*),intent(in) :: fctl
real*8,allocatable,intent(out) :: lr(:),hr(:)
integer*4 :: n1,n2,n3,nsiz,nseg,i,j,k,na
real*8 :: AL(3,3),dr,rmax,rr(3),fa
integer*4,allocatable :: cn(:)
n1=dat%n1
n2=dat%n2
n3=dat%n3
AL=dat%AL
select case(trim(fctl))
  case ("diag")
  rmax=0.5d0*dsqrt(AL(1,1)**2+AL(2,2)**2+AL(3,3)**2)
  case ("edge")
  rmax=0.5d0*min(AL(1,1),AL(2,2),AL(3,3))
  case DEFAULT
  write(6,*) "wrong control char : ", trim(fctl)
end select
! take the digonal line as the line mesh
!
! collected points owed to the mesh point
! **](***](***](***](***](**
! +----+----+----+----+----+
! ^ mesh grid point
! 1    2    3    4 .. nseg nsiz
nsiz=min(n1,n2,n3)+1
nseg=nsiz-1
dr=rmax/nseg
allocate(lr(nsiz),hr(nsiz),cn(nsiz))
lr=[((i-1)*dr,i=1,nsiz)]
hr=0.d0
cn=0
do k=1,n3
rr(3)=dabs(dist_min((k-1.d0)/n3-cc(3)))*AL(3,3)
do j=1,n2
rr(2)=dabs(dist_min((j-1.d0)/n2-cc(2)))*AL(2,2)
do i=1,n1
rr(1)=dabs(dist_min((i-1.d0)/n1-cc(1)))*AL(1,1)
fa=norm2(rr)/dr
na=floor(fa)
fa=fa-na
if(fa.le.0.5d0) then
  na=na+1
else
  na=na+2
end if
if(na.lt.1.or.na.gt.nsiz) then
  ! write(6,*) "na broken, stop"
  ! write(6,*) "na = ", na
  ! write(6,*) "nsiz=", nsiz
  ! stop
  cycle
end if
cn(na)=cn(na)+1
hr(na)=hr(na)+dat%rho(i,j,k)
end do
end do
end do
hr=hr/cn
end subroutine getrad

function dist_min(x) result(y)
        implicit none
        real*8,intent(in) :: x
        real*8 :: y
        y=mod(mod(x,1.d0)+1.5d0,1.d0)-0.5d0
end function dist_min


function check_file(frho) result(iqst)
implicit none
character(len=*) :: frho
logical :: iqst
inquire(file=trim(frho),exist=iqst)
if(.not.iqst) then
  write(6,*) trim(frho), " not exist, stop"
  stop
end if
end function check_file
subroutine write_rho(file,p)
implicit none
character(len=*) :: file
type(rho_data) :: p
integer*4 :: n1,n2,n3,nnodes_tmp,nr_n,iread,i,j,k,ii,jj
real*8 :: AL(3,3)
real*8,allocatable :: vr0(:,:,:),vr_tmp(:)
nnodes_tmp=p%nnodes_tmp
n1=p%n1
n2=p%n2
n3=p%n3
AL=p%AL     
allocate(vr0(n1,n2,n3))
vr0=p%rho
open(13,file=trim(adjustl(file)),form="unformatted")
rewind(13)
write(13) n1,n2,n3,nnodes_tmp
write(13) AL
nr_n=(n1*n2*n3)/nnodes_tmp
allocate(vr_tmp(nr_n))
do iread=1,nnodes_tmp
   do ii=1,nr_n
      jj=ii+(iread-1)*nr_n
      i=(jj-1)/(n2*n3)+1
      j=(jj-1-(i-1)*n2*n3)/n3+1
      k=jj-(i-1)*n2*n3-(j-1)*n3
      vr_tmp(ii)=vr0(i,j,k)
   enddo
   if(iread.eq.1) then
   write(6,*) "vr_tmp first 9 "
   write(6,'(3(E14.7,1x))') vr_tmp(1:3)
   write(6,'(3(E14.7,1x))') vr_tmp(4:6)
   write(6,'(3(E14.7,1x))') vr_tmp(7:9)
   endif
   write(13) vr_tmp
enddo
close(13)
deallocate(vr0)
deallocate(vr_tmp)
end subroutine write_rho

subroutine read_rho(file,p)
implicit none
character(len=*) :: file
type(rho_data) :: p
integer*4 :: n1,n2,n3,nnodes_tmp,nr_n,iread,i,j,k,ii,jj
real*8 :: AL(3,3)
real*8,allocatable :: vr0(:,:,:),vr_tmp(:)
open(11,file=trim(adjustl(file)),form="unformatted")
rewind(11)
read(11) n1,n2,n3,nnodes_tmp
read(11) AL
nr_n=(n1*n2*n3)/nnodes_tmp
allocate(vr0(n1,n2,n3))
allocate(vr_tmp(nr_n))
do iread=1,nnodes_tmp
   read(11) vr_tmp
   if(iread.eq.1) then
   write(6,*) "read file1 ..."
   write(6,*) "vr_tmp first 9 "
   write(6,'(3(E14.7,1x))') vr_tmp(1:3)
   write(6,'(3(E14.7,1x))') vr_tmp(4:6)
   write(6,'(3(E14.7,1x))') vr_tmp(7:9)
  endif
  do ii=1,nr_n
      jj=ii+(iread-1)*nr_n
      i=(jj-1)/(n2*n3)+1
      j=(jj-1-(i-1)*n2*n3)/n3+1
      k=jj-(i-1)*n2*n3-(j-1)*n3
      vr0(i,j,k)=vr_tmp(ii)
   enddo
enddo
close(11)
deallocate(vr_tmp)
p%n1=n1
p%n2=n2
p%n3=n3
p%nnodes_tmp=nnodes_tmp
p%AL=AL
allocate(p%rho(n1,n2,n3))
p%rho=vr0
deallocate(vr0)
end subroutine read_rho


end program rad_rho
