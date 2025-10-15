
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
  Integer :: mCol(maxne*3*nf),mRow(nf+1)

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
  beta=1d-2
  
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
                b(i)=b(i)+Aij*pres_w(j)
              endif
            enddo
          endif
        enddo
      enddo
    enddo

    a(i)=pres_f(i)
    b(i)=b(i)+alfa*((rho0-rho(i))/(dt**2))   !DENSITY
    b(i)=b(i)+(rho0*nablaU/dt)               !VELOCITY 
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




end module
