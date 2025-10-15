module gpu_sorting
contains


attributes(global) subroutine keys(nn,xyz,key,newID,xmin,ymin,zmin,h,nzx,nzy,nzz,xmax,ymax,zmax)
Implicit none

  Integer,value,intent(in) :: nn
  Integer,value,intent(in) :: nzx,nzy,nzz
  Real(8),value,intent(in) :: xmin,ymin,zmin,h
  Real(8),value,intent(in) :: xmax,ymax,zmax
  Real(8)                  :: xyz(3,nn)

  Integer :: i
  Integer :: xNow,yNow,zNow
  Integer :: newID(nn)

  Real(8) :: key(nn)
  Real(8) :: dist,distMax
  Real(8) :: xMm,yMm
  Real(8) :: xyzmin(3),xyzmax(3)

  i=(blockIdx%x-1)*blockDim%x+threadIdx%x
  if(i>nn) goto 999

  key(i)=0

  xNow=int((xyz(1,i)-xmin)/(2*h))+1
  yNow=int((xyz(2,i)-ymin)/(2*h))+1
  zNow=int((xyz(3,i)-zmin)/(2*h))+1

  if(xNow>nzx) xNow=nzx
  if(yNow>nzy) yNow=nzy
  if(zNow>nzz) zNow=nzz
  
  xyzmin(1)=(xNow-1)*2d0*h+xmin
  xyzmin(2)=(yNow-1)*2d0*h+ymin
  xyzmin(3)=(zNow-1)*2d0*h+zmin
  xyzmax(1)=(xNow)*2d0*h+xmin
  xyzmax(2)=(yNow)*2d0*h+ymin
  xyzmax(3)=(zNow)*2d0*h+zmin
  
  xMm=1+(xyzmax(1)-xyzmin(1))*19
  yMm=1+(xyzmax(2)-xyzmin(2))*29
  
  distMax=xMm*yMm*(xyzmax(3)-xyzmin(3))+xMm*(xyzmax(2)-xyzmin(2))+(xyzmax(1)-xyzmin(1))+1
  dist=xMm*yMm*(xyz(3,i)-xyzmin(3))+xMm*(xyz(2,i)-xyzmin(2))+(xyz(1,i)-xyzmin(1))
  
  key(i)=real(nzx*nzy*(zNow-1)+nzx*(yNow-1)+xNow,8)
  key(i)=key(i)+dist/(distMax+1)
  
  newID(i)=i

999 continue
end subroutine keys


attributes(global) subroutine upDumMatR(arr1,arr2,n_ini,n,mat,matDum)
Implicit none

  Integer,value :: n_ini,n,arr1,arr2
  Real(8) :: mat(arr1,arr2)
  Real(8) :: matDum(arr1*arr2)
  Integer :: i,j

  i=n_ini+(blockIdx%x-1)*blockDim%x+threadIdx%x
  if(i>n_ini.and.i<=n) then
    do j=1,arr1
      matDum(arr1*(i-1)+j)=mat(j,i)
    enddo
  endif

end subroutine upDumMatR

attributes(global) subroutine upDumMatI(arr1,arr2,n_ini,n,mat,matDum)
Implicit none

  Integer,value :: n_ini,n,arr1,arr2
  Integer :: mat(arr1,arr2)
  Integer :: matDum(arr1*arr2)
  Integer :: i,j

  i=n_ini+(blockIdx%x-1)*blockDim%x+threadIdx%x
  if(i>n_ini.and.i<=n) then
    do j=1,arr1
      matDum(arr1*(i-1)+j)=mat(j,i)
    enddo
  endif

end subroutine upDumMatI

attributes(global) subroutine upDumVecR(arr,n_ini,n,vec,vecDum)
Implicit none

  Integer,value :: n_ini,n,arr
  Real(8) :: vec(arr)
  Real(8) :: vecDum(arr)
  Integer :: i

  i=n_ini+(blockIdx%x-1)*blockDim%x+threadIdx%x
  if(i>n_ini.and.i<=n) then
    vecDum(i)=vec(i)
  endif

end subroutine upDumVecR

attributes(global) subroutine upDumVecI(arr,n_ini,n,vec,vecDum)
Implicit none

  Integer,value :: n_ini,n,arr
  Integer :: vec(arr)
  Integer :: vecDum(arr)
  Integer :: i

  i=n_ini+(blockIdx%x-1)*blockDim%x+threadIdx%x
  if(i>n_ini.and.i<=n) then
    vecDum(i)=vec(i)
  endif

end subroutine upDumVecI


attributes(global) subroutine orderingMatR(arr1,arr2,n_ini,n,newID,mat,matDum)
Implicit none

  Integer,value :: n_ini,n,arr1,arr2
  Integer :: newID(arr2)
  Real(8) :: mat(arr1,arr2)
  Real(8) :: matDum(arr1*arr2)
  Integer :: i,j

  i=n_ini+(blockIdx%x-1)*blockDim%x+threadIdx%x
  if(i>n_ini.and.i<=n) then
    do j=1,arr1
      mat(j,i)=matDum(arr1*(newID(i)-1)+j)
    enddo
  endif

end subroutine orderingMatR

attributes(global) subroutine orderingMatI(arr1,arr2,n_ini,n,newID,mat,matDum)
Implicit none

  Integer,value :: n_ini,n,arr1,arr2
  Integer :: newID(arr2)
  Integer :: mat(arr1,arr2)
  Integer :: matDum(arr1*arr2)
  Integer :: i,j

  i=n_ini+(blockIdx%x-1)*blockDim%x+threadIdx%x
  if(i>n_ini.and.i<=n) then
    do j=1,arr1
      mat(j,i)=matDum(arr1*(newID(i)-1)+j)
    enddo
  endif

end subroutine orderingMatI

attributes(global) subroutine orderingVecR(arr,n_ini,n,newID,vec,vecDum)
Implicit none

  Integer,value :: n_ini,n,arr
  Integer :: newID(arr)
  Real(8) :: vec(arr)
  Real(8) :: vecDum(arr)
  Integer :: i

  i=n_ini+(blockIdx%x-1)*blockDim%x+threadIdx%x
  if(i>n_ini.and.i<=n) then
    vec(i)=vecDum(newID(i))
  endif

end subroutine orderingVecR

attributes(global) subroutine orderingVecI(arr,n_ini,n,newID,vec,vecDum)
Implicit none

  Integer,value :: n_ini,n,arr
  Integer :: newID(arr)
  Integer :: vec(arr)
  Integer :: vecDum(arr)
  Integer :: i

  i=n_ini+(blockIdx%x-1)*blockDim%x+threadIdx%x
  if(i>n_ini.and.i<=n) then
    vec(i)=vecDum(newID(i))
  endif

end subroutine orderingVecI


end module gpu_sorting



subroutine sorting
use cudafor
use thrust

use properties
use parameters
use gpu_properties

use gpu_sorting

Implicit none

  Integer :: p,istat


  call keys<<<bl_fluid,Threads>>>(nf,xyz_f_d,key_d,&
                                  newID_d,xmin,ymin,zmin,h, &
                                  nzx,nzy,nzz,xmax,ymax,zmax)

  istat=cudaDeviceSynchronize
  call sort_by_key_double(key_d,nf,newID_d)
  istat=cudaDeviceSynchronize


  !Transfering xyz
  call upDumMatR<<<bl_fluid,Threads>>>(3,nf,0,nf,xyz_f_d,a_d)
  call orderingMatR<<<bl_fluid,Threads>>>(3,nf,0,nf,newID_d,xyz_f_d,a_d)

  !Transfering vel
  call upDumMatR<<<bl_fluid,Threads>>>(3,nf,0,nf,vel_f_d,a_d)
  call orderingMatR<<<bl_fluid,Threads>>>(3,nf,0,nf,newID_d,vel_f_d,a_d)

  !Transfering pres
  call upDumVecR<<<bl_fluid,Threads>>>(nf,0,nf,pres_f_d,a_d)
  call orderingVecR<<<bl_fluid,Threads>>>(nf,0,nf,newID_d,pres_f_d,a_d)

return
end subroutine sorting




subroutine sorting_wall
use cudafor
use thrust

use properties
use parameters
use gpu_properties

use gpu_sorting

Implicit none

  Integer :: p,istat


  call keys<<<bl_wall,Threads>>>(nwall,xyz_w_d,key_d,&
                                 newID_d,xmin,ymin,zmin,h, &
                                 nzx,nzy,nzz,xmax,ymax,zmax)

  istat=cudaDeviceSynchronize
  call sort_by_key_double(key_d,nwall,newID_d)
  istat=cudaDeviceSynchronize


  !Transfering xyz
  call upDumMatR   <<<bl_wall,Threads>>>(3,nwall,0,nwall,xyz_w_d,a_d)
  call orderingMatR<<<bl_wall,Threads>>>(3,nwall,0,nwall,newID_d,xyz_w_d,a_d)

  !Transfering vel
  call upDumMatR   <<<bl_wall,Threads>>>(3,nwall,0,nwall,vel_w_d,a_d)
  call orderingMatR<<<bl_wall,Threads>>>(3,nwall,0,nwall,newID_d,vel_w_d,a_d)

  !Transfering xyzbound
  call upDumMatR   <<<bl_wall,Threads>>>(3,nwall,0,nwall,xyzbound_w_d,a_d)
  call orderingMatR<<<bl_wall,Threads>>>(3,nwall,0,nwall,newID_d,xyzbound_w_d,a_d)


return
end subroutine sorting_wall


































