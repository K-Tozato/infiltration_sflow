! \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_
module coord
  integer :: nx, ny, nr
  double precision ::  dx, dy 
  double precision, allocatable :: xx(:), yy(:), gl(:,:)
end module

module time
  integer :: nti, tstep, tt
  double precision :: dt, dti, ww(2), wout
end module

module inp_data
  integer :: fl 
  double precision, allocatable :: k0(:,:), thi(:,:), ths(:,:), psi(:,:), nn(:,:)
  double precision, allocatable :: rain(:,:), qin(:,:), ax(:,:), ay(:,:), slp(:,:), raint(:), nrain(:)
  integer, allocatable :: rid(:,:)
end module

module solve_data
  double precision, allocatable :: hs(:,:), hsn(:,:), hmx(:,:), zz(:,:), zzn(:,:), fn(:,:)
end module
! \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_

!------------------------------------------------------------------------
program gwater
use coord, only : nx, ny, nr, dx, dy
use time
use inp_data, only : k0, raint, rain, nrain
use solve_data, only : hs
implicit none
double precision :: dts
integer :: val(8), i
character(10) :: ct(3)

open(11,file='./raindata/time_condition.txt',status='old')
read(11,*) nti, dti
close(11)
tt = 0 ; tstep = 0

call datain
call angle
call readrain

ww(2) = 0.d0 ; wout = 0.d0
dt = 0.01d0 ; dts = 0.d0
do
  tt = tt + 1
  if(tt == 1) then
    tstep = tstep + 1
    call date_and_time(ct(1),ct(2),ct(3),val)
    write(*,'(a,i4,a,i3,a,i2,a,i2)') "time step : ", tstep, " |" , val(5), ":", val(6), ":", val(7)
    write(*,'(a,f5.1)') "raindata [mm/h] : ", rain(tstep,1)*1000.d0*3600.d0
    write(*,'(a,f8.4)') "time     [h]    : ", dts/3600.d0
    do i = 1, nr
      ww(2) = ww(2) + dble(nrain(i)) * (rain(tstep,i)+rain(tstep+1,i))*dti*0.5d0
    end do
  end if  
  
  dts = dts + dt
  
  do i = 1, nr
    raint(i) = rain(tstep,i) + (rain(tstep+1,i)-rain(tstep,i))/dti * (dts - dble(tstep-1)*dti)
  end do

  call flow
  if(dts >= dble(tstep)*dti) then
    write(*,'(a,f7.4)') "dt : ", dt
    write(*,'(a,2f12.3)') "ww : ", ww(1), ww(2)
    write(*,'(a,f12.3,a)') "error : ", (ww(1) - ww(2))/ww(2)*100.d0, "%"
    write(*,*) "- - - - - - - - - - - - - - - - - "
    tt = 0
    call data_out
    if(tstep == nti) exit
  end if
end do

end program 
!------------------------------------------------------------------------

! \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_
subroutine datain
use coord
use time, only : dt, ww, nti
use solve_data
use inp_data
implicit none
integer :: i, j, k
character(100) :: fname, fname1
open(11,file='./input/num_node.txt',status='old')
read(11,*) nx
read(11,*) ny
read(11,*) fl
close(11)

allocate(xx(nx), yy(ny), gl(nx,ny))
allocate(hs(nx-1,ny-1), hsn(nx-1,ny-1), hmx(nx-1,ny-1), zz(nx-1,ny-1), zzn(nx-1,ny-1), fn(nx-1,ny-1))
allocate(ax(nx,ny-1), ay(nx-1,ny), qin(nx-1,ny-1), slp(nx-1,ny-1))
allocate(k0(nx,ny),psi(nx,ny),thi(nx,ny),ths(nx,ny),nn(nx,ny))

open(11,file='./input/coordinate.txt',status='old')
open(12,file='./input/parameter_infil.txt',status='old')
i = 0 ; j = 1
do k = 1, nx*ny
  i = i + 1
  read(11,*) xx(i), yy(j), gl(i,j)
  read(12,*) k0(i,j), psi(i,j), thi(i,j), ths(i,j), nn(i,j)
  if(i==nx) then
    i = 0
    j = j + 1
  endif
end do
close(11)
close(12)
dx = xx(2) - xx(1)
dy = yy(2) - yy(1)

hs(:,:) = 0.d0
zz(:,:) = 0.d0 ; zzn(:,:) = 0.d0

end subroutine
! \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_
subroutine readrain
use coord, only : nr, nx, ny, xx, yy
use time, only :  nti
use inp_data, only : rain, raint, rid, nrain
implicit none
integer :: i, j, k
double precision, allocatable :: xr(:), yr(:)
double precision :: dmin, dd
character(100) :: fname, fname1

nr = 0
fname = './raindata/coordinate_rain.txt'
open(11,file=fname,status='old')
do 
  read(11,*, iostat=i) 
  if(i < 0) exit
  nr = nr + 1
end do
close(11)

allocate(xr(nr), yr(nr))
allocate(raint(nr), rain(nti+1,nr), rid(nx,ny), nrain(nr))

fname = './raindata/coordinate_rain.txt'
open(11,file=fname,status='old')
do i = 1,nr
  read(11,*) xr(i), yr(i)
end do
close(11)

fname1 = './raindata/rain.txt'
open(11,file=fname1,status='old')
do i = 1,nti+1
  read(11,*) (rain(i,j), j=1,nr)
end do
close(11)

! calculate rid and nrain
nrain(:) = 0
if(nr==1) then
  nrain = nx*ny
  rid(:,:) = 1
end if

if(nr > 1) then
  do j = 1,ny
    do i = 1,nx
      dmin = maxval(xx) - minval(xx)
      do k = 1,nr
        dd = dsqrt( (xx(i)-xr(k))**2.d0 + (yy(j)-yr(k))**2.d0 )
        if(dd < dmin) then
          dmin = dd ; rid(i,j) = k
        end if
      end do
      nrain(rid(i,j)) = nrain(rid(i,j)) + 1
    end do
  end do
end if

end subroutine
! \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_
subroutine angle
use coord, only : nx, ny, dx, dy, gl, xx, yy
use inp_data, only : ax, ay, slp
implicit none
integer :: i, j
double precision :: a1, a2, aa
double precision, allocatable :: xx1(:), yy1(:), gl1(:,:)

allocate(xx1(nx-1),yy1(ny-1),gl1(nx-1,ny-1))

do j = 1,(ny-1)
  do i = 1,(nx-1)
    xx1(i) = 0.5d0 * (xx(i+1) + xx(i))
    yy1(j) = 0.5d0 * (yy(j+1) + yy(j))
    gl1(i,j) = (gl(i,j) + gl(i+1,j) + gl(i,j+1) + gl(i+1,j+1)) * 0.25d0
    a1 = 0.5d0 * ( (gl(i+1,j)+gl(i+1,j+1)) - (gl(i,j)+gl(i,j+1)) ) / dx
    a2 = 0.5d0 * ( (gl(i,j+1)+gl(i+1,j+1)) - (gl(i,j)+gl(i+1,j)) ) / dy
    aa = dsqrt(a1**2.d0 + a2**2.d0)
    slp(i,j) = datan(aa)
    if(a1==0.d0 .and. a2==0.d0) slp(i,j) = 0.d0
  end do 
end do

deallocate(xx,yy,gl)
nx = nx - 1 ; ny = ny - 1
allocate(xx(nx),yy(ny),gl(nx,ny))
xx = xx1 ; yy = yy1 ; gl = gl1
deallocate(xx1,yy1,gl1)

do j = 1,ny
  do i = 1,nx
    if(i >= 2) then
      ax(i,j) = (gl(i,j)-gl(i-1,j))/dx
    end if
    if(j >= 2) then
      ay(i,j) = (gl(i,j)-gl(i,j-1))/dy
    end if
  end do 
end do
ax(1,:) = ax(2,:) ; ay(:,1) = ay(:,2)
ax(nx+1,:) = ax(nx,:) ; ay(:,ny+1) = ay(:,ny)
end subroutine
! \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_
subroutine flow
use coord, only : nx, ny, dx, dy, gl
use time, only : dt, ww, wout
use solve_data, only : hs, hsn, hmx, fn, zz, zzn 
use inp_data, only : k0, psi, thi, ths, nn, fl, qin, raint, slp, ax, ay, rid
implicit none
integer :: i, j, k, num(2)
double precision :: h1, h2, grad(2), hh, ff, dt2, dt3, ang
double precision :: fnn, f1, fo, vv(2), rr
double precision :: qsx(2), qsy(2), qgx(2), qgy(2)
double precision, parameter :: gg = 9.81d0

qsx = 0.0d0 ; qsy = 0.0d0 ; qgx = 0.0d0 ; qgy = 0.0d0

ang = 0.1d0
!$OMP parallel
!$OMP do private(ff,fnn,fo,f1,vv,h1,h2,hh,grad,qsx,qsy,rr) reduction(+:wout,hsn)
do j = 1,ny
  do i = 1, nx
    ! Green Ampt  - - - - - - - - - - - - - - - - - - -
    rr = raint(rid(i,j))
    if(fn(i,j) .ne. 0) then
      ff = k0(i,j)*(dcos(slp(i,j))+psi(i,j)*(ths(i,j)-thi(i,j))/fn(i,j))
    else
      ff = k0(i,j)*100.d0
    end if
    qin(i,j) = rr - min(rr+hs(i,j)/dt,ff)

    if(ff > rr+hs(i,j)/dt) then
      fnn = fn(i,j) + (rr+hs(i,j)/dt)*dt
      fn(i,j) = fnn
    else
      fo = fn(i,j)
      f1 = fn(i,j)
      do k = 1,100000
        vv(1) = f1 - fo - k0(i,j)*dcos(slp(i,j))*dt - psi(i,j)*(ths(i,j)-thi(i,j))/dcos(slp(i,j)) &
              *dlog((f1*dcos(slp(i,j))+psi(i,j)*(ths(i,j)-thi(i,j)))/(fo*dcos(slp(i,j))+psi(i,j)*(ths(i,j)-thi(i,j))))
        vv(2) = 1 - psi(i,j)*(ths(i,j)-thi(i,j))/(f1*dcos(slp(i,j))+psi(i,j)*(ths(i,j)-thi(i,j)))
        fnn = f1 - vv(1)/vv(2)
        if((fnn-f1)/f1 < 1.0d-6) exit
        f1 = fnn
      end do
      fn(i,j) = fnn
    end if
    zzn(i,j) = fnn/(ths(i,j)-thi(i,j))/dcos(slp(i,j))
    ! - - - - - - - - - - - - - - - - - - - - - - - -      
    
    if(fl==1) wout = wout + qin(i,j)*dt

    if(fl==2) then
      !surface flow - - - - - - 
      hh = max(hs(i,j), 0.d0)    
      h1 = 0.d0 ; h2 = 0.d0 ; qsx(:) = 0.d0
      if(i==1) then
        grad(1) = ax(2,j)
        grad(2) = (gl(i+1,j)+hs(i+1,j) -gl(i,j)-hh)/dx
        h1 = hh
        if(hh*grad(2)<0.d0 .or. hs(i+1,j)*grad(2)>0.d0) h2 = 0.5d0*(hh + hs(i+1,j))
      else if(i==nx) then
        grad(1) = (gl(i,j)+hh -gl(i-1,j)-hs(i-1,j))/dx
        grad(2) = ax(nx,j)
        if(hh*grad(1)>0.d0 .or. hs(i-1,j)*grad(1)<0.d0) h1 = 0.5d0*(hh + hs(i-1,j))
        h2 = hh
      else 
        grad(1) = (gl(i,j)+hh -gl(i-1,j)-hs(i-1,j))/dx
        grad(2) = (gl(i+1,j)+hs(i+1,j) -gl(i,j)-hh)/dx
        if(hh*grad(1)>0.d0 .or. hs(i-1,j)*grad(1)<0.d0) h1 = 0.5d0*(hh + hs(i-1,j))
        if(hh*grad(2)<0.d0 .or. hs(i+1,j)*grad(2)>0.d0) h2 = 0.5d0*(hh + hs(i+1,j))
      end if
   
      qsx(1) = - (h1**(5.d0/3.d0))*dsign(1.d0,grad(1))*dsqrt(dabs(grad(1)))/nn(i,j)  
      qsx(2) = - (h2**(5.d0/3.d0))*dsign(1.d0,grad(2))*dsqrt(dabs(grad(2)))/nn(i,j)
  
      h1 = 0.d0 ; h2 = 0.d0 ; qsy(:) = 0.d0
      if(j==1) then
        grad(1) = ay(i,2)
        grad(2) = (gl(i,j+1)+hs(i,j+1) -gl(i,j)-hh)/dy
        h1 = hh
        if(hh*grad(2)<0.d0 .or. hs(i,j+1)*grad(2)>0.d0) h2 = 0.5d0*(hh + hs(i,j+1))
      else if(j==ny) then
        grad(1) = (gl(i,j)+hh -gl(i,j-1)-hs(i,j-1))/dy
        grad(2) = ay(i,ny)
        if(hh*grad(1)>0.d0 .or. hs(i,j-1)*grad(1)<0.d0) h1 = 0.5d0*(hh + hs(i,j-1))
        h2 = hh
      else 
        grad(1) = (gl(i,j)+hh -gl(i,j-1)-hs(i,j-1))/dy
        grad(2) = (gl(i,j+1)+hs(i,j+1) -gl(i,j)-hh)/dy
        if(hh*grad(1)>0.d0 .or. hs(i,j-1)*grad(1)<0.d0) h1 = 0.5d0*(hh + hs(i,j-1))
        if(hh*grad(2)<0.d0 .or. hs(i,j+1)*grad(2)>0.d0) h2 = 0.5d0*(hh + hs(i,j+1))
      end if
 
      qsy(1) = - (h1**(5.d0/3.d0))*dsign(1.d0,grad(1))*dsqrt(dabs(grad(1)))/nn(i,j)
      qsy(2) = - (h2**(5.d0/3.d0))*dsign(1.d0,grad(2))*dsqrt(dabs(grad(2)))/nn(i,j)
   
      if(i == 1 ) qsx(1) = min(qsx(1),0.d0)
      if(i == nx) qsx(2) = max(qsx(2),0.d0)
      if(j == 1 ) qsy(1) = min(qsy(1),0.d0)
      if(j == ny) qsy(2) = max(qsy(2),0.d0)

      f1 = - dt/dx*(qsx(2)-qsx(1)) - dt/dy*(qsy(2)-qsy(1)) + qin(i,j)*dt
     
      ff = - min(qsx(1),0.d0)/dx + max(qsx(2),0.d0)/dx - min(qsy(1),0.d0)/dy + max(qsy(2),0.d0)/dy
      if(hh < ff*dt) then
        ff = hh/(ff*dt)
        if(qsx(1)<0.d0) then
          if(i>1) hsn(i-1,j) = hsn(i-1,j) + qsx(1)*dt/dx*(1-ff)
          qsx(1) = qsx(1)*ff  
        end if
        if(qsx(2)>0.d0) then
          if(i<nx) hsn(i+1,j) = hsn(i+1,j) - qsx(2)*dt/dx*(1-ff)
          qsx(2) = qsx(2)*ff
        end if
        if(qsy(1)<0.d0) then
          if(j>1) hsn(i,j-1) = hsn(i,j-1) + qsy(1)*dt/dy*(1-ff)
          qsy(1) = qsy(1)*ff
        end if
        if(qsy(2)>0.d0) then 
          if(j<ny) hsn(i,j+1) = hsn(i,j+1) - qsy(2)*dt/dy*(1-ff)
          qsy(2) = qsy(2)*ff
        end if
        f1 = - dt/dx*(qsx(2)-qsx(1)) - dt/dy*(qsy(2)-qsy(1)) + qin(i,j)*dt
      end if

      if(i == 1 ) wout = wout - min(qsx(1),0.d0)*dt/dx
      if(i == nx) wout = wout + max(qsx(2),0.d0)*dt/dx
      if(j == 1 ) wout = wout - min(qsy(1),0.d0)*dt/dy
      if(j == ny) wout = wout + max(qsy(2),0.d0)*dt/dy
    
      hsn(i,j) = hsn(i,j) + f1
    end if
    
  end do
end do
!$OMP end do  


ww(1) = wout ; dt2 = 10.d0 ; dt3 = 10.d0
!$OMP do private(h1, h2, grad) reduction(+:ww) reduction(MIN:dt2,dt3)
do j = 1,ny
  do i = 1,nx
    if(hsn(i,j)<0.d0) hsn(i,j) = 0.d0
    hmx(i,j) = max(hsn(i,j),hmx(i,j))    

    if(i<nx)then
      if(hsn(i,j)*hsn(i+1,j)>1.d-16) then
        h1 = hsn(i,j) ; h2 = hsn(i+1,j)
        grad(1) = dsqrt(dabs(h2-h1))/dx
        if(dabs(ax(i+1,j))>ang)  dt2 = min(dt2, (dx**2.d0)*0.5d0*grad(1)*nn(i,j)/((h1+h2)*0.5d0)**(5.d0/3.d0))
        if(dabs(ax(i+1,j))<=ang) dt3 = min(dt3, 0.7d0*dx/dsqrt(gg*(h1+h2)*0.5d0))
      end if
    end if
    if(j<ny) then 
      if(hsn(i,j)*hsn(i,j+1)>1.d-16) then
        h1 = hsn(i,j) ; h2 = hsn(i,j+1)
        grad(1) = dsqrt(dabs(h2-h1))/dy
        if(dabs(ay(i,j+1))>ang)  dt2 = min(dt2, (dy**2.d0)*0.5d0*grad(1)*nn(i,j)/((h1+h2)*0.5d0)**(5.d0/3.d0))
        if(dabs(ay(i,j+1))<=ang) dt3 = min(dt3, 0.7d0*dy/dsqrt(gg*(h1+h2)*0.5d0))
      end if
    end if

    ww(1) = ww(1) + hsn(i,j) + zzn(i,j)*(ths(i,j)-thi(i,j))*dcos(slp(i,j))
  end do
end do
!$OMP end do  
!$OMP end parallel

hs = hsn ; zz = zzn
dt = max(min(dt2, dt3), 1.0d-6)
if(dt2 < 0.1d0-6) write(*,*) "dt is very small !!!!", dt2, dt3
end subroutine
! \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_
subroutine data_out
use coord, only : xx, yy, nx, ny
use time, only : tstep, nti
use inp_data, only : fl
use solve_data, only : hs, hmx, zzn
implicit none
integer :: i, j, k
character(100) :: fname, fname1

fname1 = './output/water-depth-level_'
write(fname,'(a,i4.4,a4)') trim(adjustl(fname1)),tstep,'.txt'
open(12, file=fname, status='replace')
do j = 1,ny
  do i = 1,nx
    write(12,'(2f10.2,2f8.5)') xx(i), yy(j), zzn(i,j), max(hs(i,j),0.d0)
  end do
end do
close(12)

if(tstep==nti .and. fl==2) then
  open(12, file="./output/max_s_water_level.txt", status='replace')
  do i = 1,nx
    write(12,*) (real(hmx(i,j)), j=1,ny)
  end do
  close(12)
end if

end subroutine
! \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_ \_




