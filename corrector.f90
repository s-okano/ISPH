module gpu_corrector

contains

attributes(global) subroutine pres_up(n,arr,pres,a,bound)
Implicit none

  Integer,value :: n,arr
  Real(8) :: pres(arr),a(arr)
  Integer :: bound(arr)

  Integer :: i
  
  i=(blockIdx%x-1)*blockDim%x+threadIdx%x
  if(i>n) goto 999

  if(bound(i)==1) then
    pres(i)=a(i)
    if(a(i)<0d0) pres(i)=0d0
  else
    pres(i)=0d0
  endif


999 continue
return
endsubroutine pres_up


attributes(global) subroutine corrector(nf,xyz_f,vel,pres_f,bound,&
                                        nwall,xyz_w,xyzb,pres_w,&
                                        mass,h,wBeta,rho0,dt, &
                                        partZone,xmin,ymin,zmin,nzx,nzy,nzz,nztotal)
Implicit none

  Integer,value :: nf,nwall,nztotal
  Integer,value :: nzx,nzy,nzz
  Real(8),value :: mass,h,wBeta,rho0,dt
  Real(8),value :: xmin,ymin,zmin

  Real(8) :: pres_f(nf),xyz_f(3,nf),vel(3,nf)
  Real(8) :: xyzb(3,nwall),xyz_w(3,nwall),pres_w(nwall)

  Real(8) :: dist,R
  Real(8) :: wFunc_0,wFunc_x,wFunc_y,wFunc_z
  Real(8) :: nablaPx,nablaPy,nablaPz

  Integer :: partZone(nztotal,4),bound(nf)

  Integer :: zi,zf,yi,yf,xi,xf,x,y,z,loop
  Integer :: xNow,yNow,zNow,zoneJ,zoneI
  Integer :: i,j,k

  Integer :: istat

  i=(blockIdx%x-1)*blockDim%x+threadIdx%x
  if(i>nf) goto 999
  if(bound(i)==7) goto 999

  xNow=int((xyz_f(1,i)-xmin)/(2*h))+1
  yNow=int((xyz_f(2,i)-ymin)/(2*h))+1
  zNow=int((xyz_f(3,i)-zmin)/(2*h))+1

  nablaPx=0d0; nablaPy=0d0; nablaPz=0d0

  zi=zNow-1
  zf=zNow+1
  if(zNow<=1) zi=1
  if(zNow>=nzz) zf=nzz
  yi=yNow-1
  yf=yNow+1
  if(yNow<=1) yi=1
  if(yNow>=nzy) yf=nzy
  xi=xNow-1
  xf=xNow+1
  if(xNow<=1) xi=1
  if(xNow>=nzx) xf=nzx

  do z=zi,zf
    do y=yi,yf
      do x=xi,xf
        zoneJ=nzx*nzy*(z-1)+nzx*(y-1)+x
        loop=1
        if(partZone(zoneJ,2*loop-1)>0) then
          do j=partZone(zoneJ,2*loop-1),partZone(zoneJ,2*loop)
            if(j/=i) then
              dist=(xyz_f(1,i)-xyz_f(1,j))**2 &
                  +(xyz_f(2,i)-xyz_f(2,j))**2 &
                  +(xyz_f(3,i)-xyz_f(3,j))**2
              dist=sqrt(dist)
              if(dist<2d0*h) then
                if(dist<0.1d0*h) dist=0.1d0*h
                R=dist/h
                if(R<1d0) wFunc_0=-3d0*R+(9d0/4d0)*(R**2)
                if(R>=1d0) wFunc_0=(-3d0/4d0)*((2d0-R)**2)
                wFunc_0=wFunc_0*wBeta/h
                wFunc_x=wFunc_0*(xyz_f(1,i)-xyz_f(1,j))/dist
                wFunc_y=wFunc_0*(xyz_f(2,i)-xyz_f(2,j))/dist
                wFunc_z=wFunc_0*(xyz_f(3,i)-xyz_f(3,j))/dist
                nablaPx=nablaPx+(mass/rho0)*(pres_f(j)+pres_f(i))*wFunc_x
                nablaPy=nablaPy+(mass/rho0)*(pres_f(j)+pres_f(i))*wFunc_y
                nablaPz=nablaPz+(mass/rho0)*(pres_f(j)+pres_f(i))*wFunc_z
              endif
            endif
          enddo
        endif
        loop=2
        if(partZone(zoneJ,2*loop-1)>0) then
          do j=partZone(zoneJ,2*loop-1),partZone(zoneJ,2*loop)
            dist=(xyz_f(1,i)-xyz_w(1,j))**2 &
                +(xyz_f(2,i)-xyz_w(2,j))**2 &
                +(xyz_f(3,i)-xyz_w(3,j))**2
            dist=sqrt(dist)
            if(dist<2d0*h) then
              if(dist<0.1d0*h) dist=0.1d0*h
              R=dist/h
              if(R<1d0) wFunc_0=-3d0*R+(9d0/4d0)*(R**2)
              if(R>=1d0) wFunc_0=(-3d0/4d0)*((2d0-R)**2)
              wFunc_0=wFunc_0*wBeta/h
              wFunc_x=wFunc_0*(xyz_f(1,i)-xyz_w(1,j))/dist
              wFunc_y=wFunc_0*(xyz_f(2,i)-xyz_w(2,j))/dist
              wFunc_z=wFunc_0*(xyz_f(3,i)-xyz_w(3,j))/dist
              nablaPx=nablaPx+(mass/rho0)*(pres_w(j)+pres_f(i))*wFunc_x
              nablaPy=nablaPy+(mass/rho0)*(pres_w(j)+pres_f(i))*wFunc_y
              nablaPz=nablaPz+(mass/rho0)*(pres_w(j)+pres_f(i))*wFunc_z
            endif
          enddo
        endif
      enddo
    enddo
  enddo

  vel(1,i)=vel(1,i)-dt*nablaPx/rho0 
  vel(2,i)=vel(2,i)-dt*nablaPy/rho0 
  vel(3,i)=vel(3,i)-dt*nablaPz/rho0


999 continue
return
end subroutine corrector





attributes(global) subroutine PPEexp(nf,xyz_f,vel,rho,pres_f,a,bound,&
                                     nwall,xyz_w,pres_w,vel_w,&
                                     mass,h,wBeta,rho0,eta2,dt,alfa, &
                                     partZone,xmin,ymin,zmin,nzx,nzy,nzz,nztotal)
Implicit none

  Integer,value :: nf,nwall,maxne,nztotal
  Integer,value :: nzx,nzy,nzz
  Real(8),value :: xmin,ymin,zmin
  Real(8),value :: mass,h,wBeta,rho0,eta2,dt,alfa

  Real(8) :: xyz_f(3,nf),vel(3,nf),pres_f(nf),rho(nf)
  Real(8) :: xyz_w(3,nwall),pres_w(nwall),vel_w(3,nwall)
  Real(8) :: a(nf)

  Real(8) :: dist,wFunc_0,wFunc_x,wFunc_y,wFunc_z,R
  Real(8) :: Aij,nablaU,SA,vx,vy,vz,rxW,bb

  Integer :: bound(nf)
  Integer :: partZone(nztotal,4)

  Integer :: i,j,k
  Integer :: zi,zf,yi,yf,xi,xf,x,y,z,loop
  Integer :: xNow,yNow,zNow,zoneJ,zoneI

  i=(blockIdx%x-1)*blockDim%x+threadIdx%x
  if(i>nf) goto 999
  if(bound(i)==7) goto 999

  xNow=int((xyz_f(1,i)-xmin)/(2*h))+1
  yNow=int((xyz_f(2,i)-ymin)/(2*h))+1
  zNow=int((xyz_f(3,i)-zmin)/(2*h))+1

  zi=zNow-1
  zf=zNow+1
  if(zNow<=1) zi=1
  if(zNow>=nzz) zf=nzz
  yi=yNow-1
  yf=yNow+1
  if(yNow<=1) yi=1
  if(yNow>=nzy) yf=nzy
  xi=xNow-1
  xf=xNow+1
  if(xNow<=1) xi=1
  if(xNow>=nzx) xf=nzx
 
  a(i)=0d0
  bb=0d0
  if(bound(i)==1) then
    SA=0d0; nablaU=0d0

    do z=zi,zf
      do y=yi,yf
        do x=xi,xf
          zoneJ=nzx*nzy*(z-1)+nzx*(y-1)+x
          loop=1
          if(partZone(zoneJ,2*loop-1)>0) then
            do j=partZone(zoneJ,2*loop-1),partZone(zoneJ,2*loop)
              if(j/=i) then
                dist=(xyz_f(1,i)-xyz_f(1,j))**2 &
                    +(xyz_f(2,i)-xyz_f(2,j))**2 &
                    +(xyz_f(3,i)-xyz_f(3,j))**2
                dist=sqrt(dist)
                if(dist<2d0*h) then
                  if(dist<0.1d0*h) dist=0.1d0*h
                  R=dist/h
                  if(R<1d0) wFunc_0=-3d0*R+(9d0/4d0)*(R**2)
                  if(R>=1d0) wFunc_0=(-3d0/4d0)*((2d0-R)**2)
                  wFunc_0=wFunc_0*wBeta/h
                  wFunc_x=wFunc_0*(xyz_f(1,i)-xyz_f(1,j))/dist
                  wFunc_y=wFunc_0*(xyz_f(2,i)-xyz_f(2,j))/dist
                  wFunc_z=wFunc_0*(xyz_f(3,i)-xyz_f(3,j))/dist
                  nablaU=nablaU+(mass/rho0)*((vel(1,j)-vel(1,i))*wFunc_x &
                                            +(vel(2,j)-vel(2,i))*wFunc_y &
                                            +(vel(3,j)-vel(3,i))*wFunc_z)
                  rxW=(xyz_f(1,i)-xyz_f(1,j))*wFunc_x &
                     +(xyz_f(2,i)-xyz_f(2,j))*wFunc_y &
                     +(xyz_f(3,i)-xyz_f(3,j))*wFunc_z
                  Aij=2d0*mass*(rxW/rho0)/(dist**2+eta2)
                  SA=SA+Aij
                  if(bound(j)==1) then
                    bb=bb+Aij*pres_f(j)
                  endif
                endif
              endif
            enddo
          endif
          loop=2
          if(partZone(zoneJ,2*loop-1)>0) then
            do j=partZone(zoneJ,2*loop-1),partZone(zoneJ,2*loop)
              if(j/=i) then
                dist=(xyz_f(1,i)-xyz_w(1,j))**2 &
                    +(xyz_f(2,i)-xyz_w(2,j))**2 &
                    +(xyz_f(3,i)-xyz_w(3,j))**2
                dist=sqrt(dist)
                if(dist<2d0*h) then
                  if(dist<0.1d0*h) dist=0.1d0*h
                  R=dist/h
                  if(R<1d0) wFunc_0=-3d0*R+(9d0/4d0)*(R**2)
                  if(R>=1d0) wFunc_0=(-3d0/4d0)*((2d0-R)**2)
                  wFunc_0=wFunc_0*wBeta/h
                  wFunc_x=wFunc_0*(xyz_f(1,i)-xyz_w(1,j))/dist
                  wFunc_y=wFunc_0*(xyz_f(2,i)-xyz_w(2,j))/dist
                  wFunc_z=wFunc_0*(xyz_f(3,i)-xyz_w(3,j))/dist
                  nablaU=nablaU+(mass/rho0)*((vel_w(1,j)-vel(1,i))*wFunc_x &
                                            +(vel_w(2,j)-vel(2,i))*wFunc_y &
                                            +(vel_w(3,j)-vel(3,i))*wFunc_z)
                  rxW=(xyz_f(1,i)-xyz_w(1,j))*wFunc_x &
                     +(xyz_f(2,i)-xyz_w(2,j))*wFunc_y &
                     +(xyz_f(3,i)-xyz_w(3,j))*wFunc_z
                  Aij=2d0*mass*(rxW/rho0)/(dist**2+eta2)
                  SA=SA+Aij
                  bb=bb+Aij*pres_w(j)
                endif
              endif
            enddo
          endif
        enddo
      enddo
    enddo

    bb=bb+alfa*((rho0-rho(i))/(dt**2))
    bb=bb+(rho0*nablaU/dt)
    a(i)=(bb/SA)
    !!a(i)=0.5d0*(bb/SA)+0.5d0*pres_f(i)
  else
    a(i)=0d0
    pres_f(i)=0d0
  endif

 
999 continue
return
end subroutine PPEexp




attributes(global) subroutine defMatPres(nf,xyz_f,vel,rho,pres_f,&
                                         bound,M,a,b,mCol,mRow, &
                                         nwall,xyz_w,pres_w,vel_w, &
                                         mass,h,wBeta,rho0,eta2,dt,alfa,maxne, &
                                         partZone,xmin,ymin,zmin,nzx,nzy,nzz,nztotal)
Implicit none

  Integer,value :: nf,nwall,maxne,nztotal
  Integer,value :: nzx,nzy,nzz
  Real(8),value :: mass,h,wBeta,rho0,eta2,dt,alfa
  Real(8),value :: xmin,ymin,zmin

  Real(8) :: xyz_f(3,nf),vel(3,nf)
  Real(8) :: pres_f(nf),rho(nf)
  Real(8) :: xyz_w(3,nwall),pres_w(nwall),vel_w(3,nwall)
  Real(8) :: M(maxne*nf),a(nf),b(nf)
  Integer :: mCol(maxne*nf),mRow(nf+1)

  Integer :: bound(nf)
  Integer :: partZone(nztotal,4)

  Integer :: zi,zf,yi,yf,xi,xf,x,y,z,loop
  Integer :: xNow,yNow,zNow,zoneJ
  
  Integer :: i,j,k
  Integer :: diag 
  Real(8) :: dist,wFunc_0,wFunc_x,wFunc_y,wFunc_z,R
  Real(8) :: SA,Aij,nablaU,rxW

  Real(8) :: beta
 
  i=(blockIdx%x-1)*blockDim%x+threadIdx%x
  if(i>nf) goto 999
  if(bound(i)==7) goto 999

  xNow=int((xyz_f(1,i)-xmin)/(2*h))+1
  yNow=int((xyz_f(2,i)-ymin)/(2*h))+1
  zNow=int((xyz_f(3,i)-zmin)/(2*h))+1
  beta=0.d0
  
  a(i)=0d0
  b(i)=0d0
  k=0   !COURTER

  zi=zNow-1
  zf=zNow+1
  if(zNow<=1) zi=1
  if(zNow>=nzz) zf=nzz
  yi=yNow-1
  yf=yNow+1
  if(yNow<=1) yi=1
  if(yNow>=nzy) yf=nzy
  xi=xNow-1
  xf=xNow+1
  if(xNow<=1) xi=1
  if(xNow>=nzx) xf=nzx
  if(bound(i)==1) then   !ONLY INNER PARTICLES
    nablaU=0d0; SA=0d0
    do z=zi,zf
      do y=yi,yf
        do x=xi,xf
          zoneJ=nzx*nzy*(z-1)+nzx*(y-1)+x
          loop=1       !WATER
          if(partZone(zoneJ,2*loop-1)>0) then
            do j=partZone(zoneJ,2*loop-1),partZone(zoneJ,2*loop)
              dist=(xyz_f(1,i)-xyz_f(1,j))**2 &
                  +(xyz_f(2,i)-xyz_f(2,j))**2 &
                  +(xyz_f(3,i)-xyz_f(3,j))**2
              dist=sqrt(dist)
              if(dist<2d0*h) then
                k=k+1
                mCol(mRow(i)+k-1)=j
                if(j==i) then
                  diag=k
                else
                  if(dist<0.1d0*h) dist=0.1d0*h
                  R=dist/h
                  if(R<1d0) wFunc_0=-3d0*R+(9d0/4d0)*(R**2)
                  if(R>=1d0) wFunc_0=(-3d0/4d0)*((2d0-R)**2)
                  wFunc_0=wFunc_0*wBeta/h
                  wFunc_x=wFunc_0*(xyz_f(1,i)-xyz_f(1,j))/dist
                  wFunc_y=wFunc_0*(xyz_f(2,i)-xyz_f(2,j))/dist
                  wFunc_z=wFunc_0*(xyz_f(3,i)-xyz_f(3,j))/dist
                  nablaU=nablaU+(mass/rho0)*((vel(1,j)-vel(1,i))*wFunc_x &
                                            +(vel(2,j)-vel(2,i))*wFunc_y &
                                            +(vel(3,j)-vel(3,i))*wFunc_z)
                  rxW=(xyz_f(1,i)-xyz_f(1,j))*wFunc_x &
                     +(xyz_f(2,i)-xyz_f(2,j))*wFunc_y &
                     +(xyz_f(3,i)-xyz_f(3,j))*wFunc_z
                  Aij=2d0*mass*(rxW/rho0)/(dist**2+eta2)
                  SA=SA+Aij

                  a(i)=a(i)+Aij*pres_f(j)

                  if(bound(j)==1) then
                    M(mRow(i)+k-1)=-Aij
                  endif
                endif
              endif
            enddo
          endif

          loop=2          !WALL SEARCH
          if(partZone(zoneJ,2*loop-1)>0) then
            do j=partZone(zoneJ,2*loop-1),partZone(zoneJ,2*loop)
              dist=(xyz_f(1,i)-xyz_w(1,j))**2 &
                  +(xyz_f(2,i)-xyz_w(2,j))**2 &
                  +(xyz_f(3,i)-xyz_w(3,j))**2
              dist=sqrt(dist)
              if(dist<2d0*h) then
                if(dist<0.1d0*h) dist=0.1d0*h
                R=dist/h
                if(R<1d0) wFunc_0=-3d0*R+(9d0/4d0)*(R**2)
                if(R>=1d0) wFunc_0=(-3d0/4d0)*((2d0-R)**2)
                wFunc_0=wFunc_0*wBeta/h
                wFunc_x=wFunc_0*(xyz_f(1,i)-xyz_w(1,j))/dist
                wFunc_y=wFunc_0*(xyz_f(2,i)-xyz_w(2,j))/dist
                wFunc_z=wFunc_0*(xyz_f(3,i)-xyz_w(3,j))/dist
                nablaU=nablaU+(mass/rho0)*((vel_w(1,j)-vel(1,i))*wFunc_x &
                                          +(vel_w(2,j)-vel(2,i))*wFunc_y &
                                          +(vel_w(3,j)-vel(3,i))*wFunc_z)
                rxW=(xyz_f(1,i)-xyz_w(1,j))*wFunc_x &
                   +(xyz_f(2,i)-xyz_w(2,j))*wFunc_y &
                   +(xyz_f(3,i)-xyz_w(3,j))*wFunc_z
                Aij=2d0*mass*(rxW/rho0)/(dist**2+eta2)
                SA=SA+Aij*beta

                a(i)=a(i)+Aij*pres_w(j)

                b(i)=b(i)+Aij*pres_w(j)
              endif
            enddo
          endif
        enddo
      enddo
    enddo

    !a(i)=pres_f(i)
    b(i)=b(i)+alfa*((rho0-rho(i))/(dt**2))   !DENSITY
    b(i)=b(i)+(rho0*nablaU/dt)               !VELOCITY 

    a(i)=a(i)+alfa*((rho0-rho(i))/(dt**2))   !DENSITY
    a(i)=a(i)+(rho0*nablaU/dt)               !VELOCITY 
    a(i)=a(i)/SA

    M(mRow(i)+diag-1)=SA
  else
    do z=zi,zf
      do y=yi,yf
        do x=xi,xf
          zoneJ=nzx*nzy*(z-1)+nzx*(y-1)+x
          loop=1
          if(partZone(zoneJ,2*loop-1)>0) then
            do j=partZone(zoneJ,2*loop-1),partZone(zoneJ,2*loop)
              dist=(xyz_f(1,i)-xyz_f(1,j))**2 &
                  +(xyz_f(2,i)-xyz_f(2,j))**2 &
                  +(xyz_f(3,i)-xyz_f(3,j))**2
              dist=sqrt(dist)
              if(dist<2d0*h) then
                k=k+1
                mCol(mRow(i)+k-1)=j
                M(mRow(i)+k-1)=0d0
              endif
            enddo
          endif
        enddo
      enddo
    enddo

  endif 
 
999 continue
return
end subroutine defMatPres




attributes(global) subroutine x_up(nf,xyz_f,vel,bound,&
                                   nwall,xyz_w,xyzb,&
                                   h,dt, &
                                   partZone,xmin,ymin,zmin,nzx,nzy,nzz,nztotal)
Implicit none

  Integer,value :: nf,nwall
  Integer,value :: nzx,nzy,nzz,nztotal 
  Real(8),value :: dt,h
  Real(8),value :: xmin,ymin,zmin

  Real(8) :: xyz_f(3,nf),vel(3,nf)
  Real(8) :: xyzb(3,nwall),xyz_w(3,nwall)

  Integer :: i,j,k
  Integer :: partZone(nztotal,4),bound(nf)
  Integer :: zi,zf,yi,yf,xi,xf,x,y,z,loop
  Integer :: xNow,yNow,zNow,zoneJ,zoneI

  Real(8) :: dist,nx,ny,nz,norm
  Real(8) :: dist_min,VxN,nn_b(3)

  i=(blockIdx%x-1)*blockDim%x+threadIdx%x
  if(i>nf) goto 999
  if(bound(i)==7) goto 999

  xNow=int((xyz_f(1,i)-xmin)/(2*h))+1
  yNow=int((xyz_f(2,i)-ymin)/(2*h))+1
  zNow=int((xyz_f(3,i)-zmin)/(2*h))+1

  k=0
  dist_min=1.1d0*h/1.2d0

  zi=zNow-1
  zf=zNow+1
  if(zNow<=1) zi=1
  if(zNow>=nzz) zf=nzz
  yi=yNow-1
  yf=yNow+1
  if(yNow<=1) yi=1
  if(yNow>=nzy) yf=nzy
  xi=xNow-1
  xf=xNow+1
  if(xNow<=1) xi=1
  if(xNow>=nzx) xf=nzx
  do z=zi,zf
    do y=yi,yf
      do x=xi,xf
        zoneJ=nzx*nzy*(z-1)+nzx*(y-1)+x
        loop=2
        if(partZone(zoneJ,2*loop-1)>0) then
          do j=partZone(zoneJ,2*loop-1),partZone(zoneJ,2*loop)
            dist=(xyz_f(1,i)-xyz_w(1,j))**2 &
                +(xyz_f(2,i)-xyz_w(2,j))**2 &
                +(xyz_f(3,i)-xyz_w(3,j))**2
            dist=sqrt(dist)
            if(dist<1d0*h) then
              norm=sqrt(xyzb(1,j)**2+xyzb(2,j)**2+xyzb(3,j)**2)
              nx=xyzb(1,j)/norm
              ny=xyzb(2,j)/norm
              nz=xyzb(3,j)/norm
              VxN=vel(1,i)*nx+vel(2,i)*ny+vel(3,i)*nz
              if(dist<dist_min.and.VxN<0d0) then
                dist_min=dist
                k=j
                nn_b(1)=nx
                nn_b(2)=ny
                nn_b(3)=nz
              endif
            endif
          enddo
        endif
      enddo
    enddo
  enddo

  if(k/=0) then
    VxN=vel(1,i)*nn_b(1)+vel(2,i)*nn_b(2)+vel(3,i)*nn_b(3)
    vel(1,i)=vel(1,i)-VxN*nn_b(1)
    vel(2,i)=vel(2,i)-VxN*nn_b(2)
    vel(3,i)=vel(3,i)-VxN*nn_b(3)
  endif

  xyz_f(1,i)=xyz_f(1,i)+vel(1,i)*dt 
  xyz_f(2,i)=xyz_f(2,i)+vel(2,i)*dt
  xyz_f(3,i)=xyz_f(3,i)+vel(3,i)*dt

999 continue
return
end subroutine x_up







attributes(global) subroutine virt1(nf,xyz_f,vel,visc,veddy, &
                                    nwall,xyz_w,xyzb,vel_w,pres_w,&
                                    mass,h,wBeta,rho0,eta2,grav, &
                                    partZone,xmin,ymin,zmin,nzx,nzy,nzz,nztotal,phi,pi)
Implicit none

  Integer,value :: nf,nwall
  Integer,value :: nzx,nzy,nzz,nztotal
  Real(8),value :: mass,h,wBeta,rho0,eta2
  Real(8),value :: xmin,ymin,zmin,phi,pi
  Real(8) :: grav(3)

  Real(8) :: xyz_f(3,nf),vel(3,nf),veddy(nf),visc(nf)
  Real(8) :: xyz_w(3,nwall),xyzb(3,nwall),vel_w(3,nwall),pres_w(nwall)

  Real(8) :: dist,R,wFunc
  Real(8) :: ved,sumW,temp
  Real(8) :: vx,vy,vz

  Real(8) :: vi,rxW
  Real(8) :: wFunc_0,wFunc_x,wFunc_y,wFunc_z
  Real(8) :: LapUx,LapUy,LapUz

  Real(8) :: vnn
  Real(8) :: nx,ny,nz,res,const

  Integer :: partZone(nztotal,4)
  Integer :: zi,zf,yi,yf,xi,xf,x,y,z,loop
  Integer :: xNow,yNow,zNow,zoneJ,zoneI

  Integer :: i,j,k
  Integer :: cont

  i=(blockIdx%x-1)*blockDim%x+threadIdx%x
  if(i>nwall) goto 999

  xNow=int((xyz_w(1,i)+xyzb(1,i)-xmin)/(2*h))+1
  yNow=int((xyz_w(2,i)+xyzb(2,i)-ymin)/(2*h))+1
  zNow=int((xyz_w(3,i)+xyzb(3,i)-zmin)/(2*h))+1

  vx=0d0; vy=0d0; vz=0d0; sumW=0d0
  cont=0

  zi=zNow-1
  zf=zNow+1
  if(zNow<=1) zi=1
  if(zNow>=nzz) zf=nzz
  yi=yNow-1
  yf=yNow+1
  if(yNow<=1) yi=1
  if(yNow>=nzy) yf=nzy
  xi=xNow-1
  xf=xNow+1
  if(xNow<=1) xi=1
  if(xNow>=nzx) xf=nzx
  do z=zi,zf
    do y=yi,yf
      do x=xi,xf
        zoneJ=nzx*nzy*(z-1)+nzx*(y-1)+x
        loop=1
        if(partZone(zoneJ,2*loop-1)>0) then
          do j=partZone(zoneJ,2*loop-1),partZone(zoneJ,2*loop)
            dist=(xyz_w(1,i)+xyzb(1,i)-xyz_f(1,j))**2 &
                +(xyz_w(2,i)+xyzb(2,i)-xyz_f(2,j))**2 &
                +(xyz_w(3,i)+xyzb(3,i)-xyz_f(3,j))**2
            dist=sqrt(dist)
            if(dist<2d0*h) then
              if(dist<0.1d0*h) dist=0.1d0*h
              R=dist/h
              if(R<1d0) wFunc=1d0-(3d0/2d0)*(R**2)+(3e0/4e0)*(R**3)
              if(R>=1d0) wFunc=(1d0/4d0)*((2d0-R)**3)
              wFunc=wFunc*wBeta
              temp=(mass/rho0)*wFunc
              sumW=sumW+temp
              vx=vx+temp*vel(1,j)
              vy=vy+temp*vel(2,j)
              vz=vz+temp*vel(3,j)
              cont=cont+1
            endif
          enddo
        endif
      enddo
    enddo
  enddo
  
  if(cont>=1) then
    if(abs(sumW)<=1d-5*h) then
      vx=0d0; vy=0d0; vz=0d0
    else
      vx=vx/sumW
      vy=vy/sumW
      vz=vz/sumW
    endif
    LapUx=0d0; LapUy=0d0; LapUz=0d0

    do z=zi,zf
      do y=yi,yf
        do x=xi,xf
          zoneJ=nzx*nzy*(z-1)+nzx*(y-1)+x
          loop=1
            if(partZone(zoneJ,2*loop-1)>0) then
              do j=partZone(zoneJ,2*loop-1),partZone(zoneJ,2*loop)
                dist=(xyz_w(1,i)+xyzb(1,i)-xyz_f(1,j))**2 &
                    +(xyz_w(2,i)+xyzb(2,i)-xyz_f(2,j))**2 &
                    +(xyz_w(3,i)+xyzb(3,i)-xyz_f(3,j))**2
                dist=sqrt(dist)
                if(dist<2d0*h) then
                  if(dist<0.1d0*h) dist=0.1d0*h
                  R=dist/h
                  if(R<1d0) wFunc_0=-3d0*R+(9d0/4d0)*(R**2)
                  if(R>=1d0) wFunc_0=(-3d0/4d0)*((2d0-R)**2)
                  wFunc_0=wFunc_0*wBeta/h
                  wFunc_x=wFunc_0*(xyz_w(1,i)+xyzb(1,i)-xyz_f(1,j))/dist
                  wFunc_y=wFunc_0*(xyz_w(2,i)+xyzb(2,i)-xyz_f(2,j))/dist
                  wFunc_z=wFunc_0*(xyz_w(3,i)+xyzb(3,i)-xyz_f(3,j))/dist
                  vi=(2d0*(visc(j)+veddy(j)))/(rho0)
                  rxW=(xyz_w(1,i)+xyzb(1,i)-xyz_f(1,j))*wFunc_x+ &
                      (xyz_w(2,i)+xyzb(2,i)-xyz_f(2,j))*wFunc_y+ &
                      (xyz_w(3,i)+xyzb(3,i)-xyz_f(3,j))*wFunc_z
                  LapUx=LapUx+(mass*vi*rxW/(dist**2+eta2))*(vx-vel(1,j))
                  LapUy=LapUy+(mass*vi*rxW/(dist**2+eta2))*(vy-vel(2,j))
                  LapUz=LapUz+(mass*vi*rxW/(dist**2+eta2))*(vz-vel(3,j))
                endif
              enddo
            endif
          !enddo
        enddo
      enddo
    enddo

    pres_w(i)=-rho0*(xyzb(1,i)*(LapUx+grav(1))+ &
                     xyzb(2,i)*(LapUy+grav(2))+ &
                     xyzb(3,i)*(LapUz+grav(3)))
    if(pres_w(i)<0d0) pres_w(i)=0d0   !avoid wall suction

    vel_w(1,i)=0d0
    vel_w(2,i)=0d0
    vel_w(3,i)=0d0

  else
    pres_w(i)=0d0
    vel_w(1,i)=0d0
    vel_w(2,i)=0d0
    vel_w(3,i)=0d0
  endif

999 continue
end subroutine virt1


attributes(global) subroutine pres_virt(nf,xyz_f,pres_f,&
                                        nwall,xyz_w,xyzb,pres_w,&
                                        mass,h,wBeta,rho0, &
                                        partZone,xmin,ymin,zmin,nzx,nzy,nzz,nztotal)
Implicit none
  
  Integer,value :: nf,nwall
  Integer,value :: nzx,nzy,nzz,nztotal 
  Real(8),value :: mass,h,wBeta,rho0
  Real(8),value :: xmin,ymin,zmin

  Real(8) :: xyz_f(3,nf),pres_f(nf)
  Real(8) :: xyz_w(3,nwall),xyzb(3,nwall),pres_w(nwall)

  Real(8) :: dist,wFunc,pre,sumW,R

  Integer :: partZone(nztotal,4)
  Integer :: zi,zf,yi,yf,xi,xf,x,y,z,loop
  Integer :: xNow,yNow,zNow,zoneJ,zoneI
  Integer :: i,j,k
  Integer :: cont

  i=(blockIdx%x-1)*blockDim%x+threadIdx%x
  if(i>nwall) goto 999

  xNow=int((xyz_w(1,i)+xyzb(1,i)-xmin)/(2*h))+1
  yNow=int((xyz_w(2,i)+xyzb(2,i)-ymin)/(2*h))+1
  zNow=int((xyz_w(3,i)+xyzb(3,i)-zmin)/(2*h))+1

  cont=0
  pre=0d0; sumW=0d0

  zi=zNow-1
  zf=zNow+1
  if(zNow<=1) zi=1
  if(zNow>=nzz) zf=nzz
  yi=yNow-1
  yf=yNow+1
  if(yNow<=1) yi=1
  if(yNow>=nzy) yf=nzy
  xi=xNow-1
  xf=xNow+1
  if(xNow<=1) xi=1
  if(xNow>=nzx) xf=nzx
  do z=zi,zf
    do y=yi,yf
      do x=xi,xf
        zoneJ=nzx*nzy*(z-1)+nzx*(y-1)+x
        loop=1
          if(partZone(zoneJ,2*loop-1)>0) then
            do j=partZone(zoneJ,2*loop-1),partZone(zoneJ,2*loop)
              dist=(xyz_w(1,i)+xyzb(1,i)-xyz_f(1,j))**2 &
                  +(xyz_w(2,i)+xyzb(2,i)-xyz_f(2,j))**2 &
                  +(xyz_w(3,i)+xyzb(3,i)-xyz_f(3,j))**2
              dist=sqrt(dist)
              if(dist<2d0*h) then
                R=dist/h
                if(R<1d0) wFunc=1d0-(3d0/2d0)*(R**2)+(3d0/4d0)*(R**3)
                if(R>=1d0) wFunc=(1d0/4d0)*((2d0-R)**3)
                wFunc=wFunc*wBeta
                sumW=sumW+mass*wFunc/rho0
                pre=pre+mass*pres_f(j)*wFunc/rho0
                cont=cont+1
              endif
            enddo
          endif
        !enddo
      enddo
    enddo
  enddo

  if(cont>=1) then
    if(abs(sumW)>1d-5*h) pres_w(i)=pres_w(i)+pre/sumW
    if(abs(sumW)<=1d-5*h) pres_w(i)=0d0
  else
    pres_w(i)=0e0
  endif


999 continue
return
endsubroutine pres_virt



end module gpu_corrector








