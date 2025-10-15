module gpu_tracking
contains


attributes(global) subroutine partzone_zero(arr1,arr2,ptz_i,ptz_f, &
                                            nztotal,partZone)
Implicit none

  Integer,value :: arr1,arr2,nztotal
  Integer,value :: ptz_i,ptz_f

  Integer :: partZone(arr1,arr2)

  Integer :: i,j

  i=(blockIdx%x-1)*blockDim%x+threadIdx%x
  if(i>nztotal) goto 999

  do j=ptz_i,ptz_f
    partZone(i,j)=0
  enddo

999 continue
end subroutine partzone_zero



attributes(global) subroutine z_ini(n,n_ini,np,ptz_i, &
                                    nztotal,xyz,partZone, &
                                    nzx,nzy,nzz,xmin,ymin,zmin,h)
Implicit none

  Integer,value :: n,n_ini,np
  Integer,value :: ptz_i
  Integer,value :: nztotal
  Integer,value :: nzx,nzy,nzz
  Real(8),value :: xmin,ymin,zmin,h

  Real(8) :: xyz(3,np)
  Integer :: partZone(nztotal,4)

  Integer :: i
  Integer :: xZ,yZ,zZ,zoneI
  Integer :: istat

  i=n_ini+(blockIdx%x-1)*blockDim%x+threadIdx%x
  if(i>n) goto 999

  xZ=int((xyz(1,i)-xmin)/(2*h))+1
  yZ=int((xyz(2,i)-ymin)/(2*h))+1
  zZ=int((xyz(3,i)-zmin)/(2*h))+1
  zoneI=nzx*nzy*(zZ-1)+nzx*(yZ-1)+xZ
  if(zoneI<=0.or.zoneI>nztotal) goto 999

  partZone(zoneI,ptz_i)=2*np
  partZone(zoneI,ptz_i+1)=-1

999 continue
end subroutine z_ini


attributes(global) subroutine zones(n,n_ini,np,ptz_i, &
                                    nztotal,xyz,partZone, &
                                    nzx,nzy,nzz,xmin,ymin,zmin,h)
Implicit none

  Integer,value :: n,n_ini,np
  Integer,value :: ptz_i
  Integer,value :: nztotal
  Integer,value :: nzx,nzy,nzz
  Real(8),value :: xmin,ymin,zmin,h

  Real(8) :: xyz(3,np)
  Integer :: partZone(nztotal,4)

  Integer :: i
  Integer :: xZ,yZ,zZ,zoneI
  Integer :: istat

  i=n_ini+(blockIdx%x-1)*blockDim%x+threadIdx%x
  if(i>n) goto 999

  xZ=int((xyz(1,i)-xmin)/(2*h))+1
  yZ=int((xyz(2,i)-ymin)/(2*h))+1
  zZ=int((xyz(3,i)-zmin)/(2*h))+1
  zoneI=nzx*nzy*(zZ-1)+nzx*(yZ-1)+xZ
  if(zoneI<=0.or.zoneI>nztotal) goto 999

  istat=atomicmin(partZone(zoneI,ptz_i),i)
  istat=atomicmax(partZone(zoneI,ptz_i+1),i)

999 continue
end subroutine zones


end module gpu_tracking


subroutine tracking_wall
use cudafor

use properties
use parameters
use gpu_properties

use gpu_tracking

Implicit none

  Integer :: istat

  call partzone_zero<<<bl_zones,Threads>>>(nztotal,4,3,4,&
                                           nztotal,partZone_d)

  call z_ini<<<bl_wall,Threads>>>(nwall,0,nwall,3, &
                                  nztotal,xyz_w_d,&
                                  partZone_d, &
                                  nzx,nzy,nzz,xmin,ymin,zmin,h)
  call zones<<<bl_wall,Threads>>>(nwall,0,nwall,3, &
                                  nztotal,xyz_w_d,&
                                  partZone_d, &
                                  nzx,nzy,nzz,xmin,ymin,zmin,h)


return
end subroutine tracking_wall


subroutine tracking_fluid
use cudafor

use properties
use parameters
use gpu_properties

use gpu_tracking

Implicit none

  Integer :: istat

  call partzone_zero<<<bl_zones,Threads>>>(nztotal,4,1,2,&
                                           nztotal,partZone_d)

  call z_ini<<<bl_fluid,Threads>>>(nf,0,nf,1, &
                                   nztotal,xyz_f_d,&
                                   partZone_d, &
                                   nzx,nzy,nzz,xmin,ymin,zmin,h)
  call zones<<<bl_fluid,Threads>>>(nf,0,nf,1, &
                                   nztotal,xyz_f_d,&
                                   partZone_d, &
                                   nzx,nzy,nzz,xmin,ymin,zmin,h)

return
end subroutine tracking_fluid











