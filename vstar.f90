module gpu_vstar

contains
attributes(global) subroutine Init(nf,xyz_f,vel,veddy,visc,rho,bound,pres,&
                                   nwall,xyz_w,vel_w,&
                                   mass,h,wBeta,rho0,eta2,ne,&
                                   visc0,viscMax,coes,phi,&
                                   partZone,xmin,ymin,zmin,nzx,nzy,nzz,nztotal)
Implicit none

  Integer,value :: nf,nwall
  Integer,value :: nzx,nzy,nzz,nztotal 
  Real(8),value :: visc0,viscMax,coes,phi
  Real(8),value :: mass,h,wBeta,rho0,eta2
  Real(8),value :: xmin,ymin,zmin

  Real(8) :: xyz_f(3,nf),vel(3,nf),veddy(nf),visc(nf),rho(nf)
  Real(8) :: pres(nf)
  Real(8) :: xyz_w(3,nwall),vel_w(3,nwall)

  Real(8) :: dist,R,wFunc
  Real(8) :: wFunc_0,wFunc_x,wFunc_y,wFunc_z
  Real(8) :: den,vmod2

  Real(8) :: cs
  Real(8) :: delta

  Integer :: bound(nf)
  Integer :: partZone(nztotal,4)
  Integer :: ne(nf+1)

  Integer :: i,j,k
  Integer :: zi,zf,yi,yf,xi,xf,x,y,z,loop
  Integer :: xNow,yNow,zNow,zoneJ,zoneI
  Integer :: cont

  Real(8) :: SROJ

  i=(blockIdx%x-1)*blockDim%x+threadIdx%x
  if(i>nf) goto 999  

  ne(i)=0
  cont=0

  xNow=int((xyz_f(1,i)-xmin)/(2*h))+1
  yNow=int((xyz_f(2,i)-ymin)/(2*h))+1
  zNow=int((xyz_f(3,i)-zmin)/(2*h))+1

  if(xNow<=0.or.xNow>nzx.or.&
     yNow<=0.or.yNow>nzy.or.&
     zNow<=0.or.zNow>nzz) then
    bound(i)=7
    goto 999
  endif

  rho(i)=0d0; SROJ=0d0

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
            dist=(xyz_f(1,i)-xyz_f(1,j))**2 &
                +(xyz_f(2,i)-xyz_f(2,j))**2 &
                +(xyz_f(3,i)-xyz_f(3,j))**2
            dist=sqrt(dist)
            if(dist<2e0*h) then 
              R=dist/h
              if(R<1d0) wFunc=1d0-(3d0/2d0)*(R**2)+(3d0/4d0)*(R**3)
              if(R>=1d0) wFunc=(1d0/4d0)*((2d0-R)**3)
              wFunc=wFunc*wBeta
              rho(i)=rho(i)+mass*wFunc
              cont=cont+1
              ne(i)=ne(i)+1

              if(j/=i) then
                if(dist<0.1d0*h) dist=0.1d0*h
                R=dist/h
                if(R<1d0) wFunc_0=-3d0*R+(9d0/4d0)*(R**2)
                if(R>=1d0) wFunc_0=(-3d0/4d0)*((2d0-R)**2)
                wFunc_0=wFunc_0*wBeta/h
                wFunc_x=wFunc_0*(xyz_f(1,i)-xyz_f(1,j))/dist
                wFunc_y=wFunc_0*(xyz_f(2,i)-xyz_f(2,j))/dist
                wFunc_z=wFunc_0*(xyz_f(3,i)-xyz_f(3,j))/dist
                vmod2=(vel(1,j)-vel(1,i))**2+ &
                      (vel(2,j)-vel(2,i))**2+ &
                      (vel(3,j)-vel(3,i))**2
                den=2d0/(rho0)
                SROJ=SROJ+mass*den*(vmod2/(dist**2+eta2)) &
                                  *((xyz_f(1,j)-xyz_f(1,i))*wFunc_x &
                                   +(xyz_f(2,j)-xyz_f(2,i))*wFunc_y &
                                   +(xyz_f(3,j)-xyz_f(3,i))*wFunc_z)
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
              R=dist/h
              if(R<1d0) wFunc=1d0-(3d0/2d0)*(R**2)+(3d0/4d0)*(R**3)
              if(R>=1d0) wFunc=(1d0/4d0)*((2d0-R)**3)
              wFunc=wFunc*wBeta
              rho(i)=rho(i)+mass*wFunc
              cont=cont+1
              if(dist<0.1e0*h) dist=0.1e0*h
              R=dist/h

              if(R<1e0) wFunc_0=-3e0*R+(9e0/4e0)*(R**2e0)
              if(R>=1e0) wFunc_0=(-3e0/4e0)*((2e0-R)**2e0)
              wFunc_0=wFunc_0*wBeta/h
              wFunc_x=wFunc_0*(xyz_f(1,i)-xyz_w(1,j))/dist
              wFunc_y=wFunc_0*(xyz_f(2,i)-xyz_w(2,j))/dist
              wFunc_z=wFunc_0*(xyz_f(3,i)-xyz_w(3,j))/dist
              vmod2=(vel(1,i)-vel_w(1,j))**2+ &
                    (vel(2,i)-vel_w(2,j))**2+ &
                    (vel(3,i)-vel_w(3,j))**2
              den=2d0/(rho0)
              SROJ=SROJ+mass*den*(vmod2/(dist**2+eta2)) &
                                *((xyz_w(1,j)-xyz_f(1,i))*wFunc_x &
                                 +(xyz_w(2,j)-xyz_f(2,i))*wFunc_y &
                                 +(xyz_w(3,j)-xyz_f(3,i))*wFunc_z)
            endif
          enddo
        endif
      enddo
    enddo
  enddo


  if(cont>48) bound(i)=1
  if(cont<=48.and.cont>=5) bound(i)=2
  if(ne(i)<=5) bound(i)=3
  if(cont<=5) bound(i)=3

  cs=1d-1
  delta=h

  SROJ=sqrt(0.5d0*SROJ)
  veddy(i)=SROJ*(cs*delta)**2

!!  !!for Non-Newtonian fluid
!!  if(SROJ<0.001d0*h) then
!!    visc(i)=viscMax
!!  else
!!    visc(i)=visc0+((coes+pres(i)*tan(phi))/(rho0*SROJ))   !Cohesion and friction angle
!!!    visc(i)=visc0+((yield)/(rho0*SROJ))  !Fixed yield stress
!!
!!    if(visc(i)<=visc0) visc(i)=visc0
!!    if(visc(i)>=viscMax) visc(i)=viscMax
!!  endif

999 continue
end subroutine Init 



attributes(global) subroutine vel_up(n,arr,vel,a,bound,grav,dt)
Implicit none

  Integer,value :: n,arr
  Real(8),value :: dt
  Real(8) :: grav(3)
  Real(8) :: vel(3,arr),a(3*arr)
  Integer :: bound(arr)

  Integer :: i

  i=(blockIdx%x-1)*blockDim%x+threadIdx%x
  if(i>n) goto 999

  vel(1,i)=a(3*(i-1)+1)
  vel(2,i)=a(3*(i-1)+2)
  vel(3,i)=a(3*(i-1)+3)

999 continue
end subroutine vel_up



attributes(global) subroutine vstarEXP(nf,xyz_f,vel,veddy,visc,a,bound,&
                                       nwall,xyz_w,vel_w,&
                                       mass,h,wBeta,rho0,eta2,dt,grav, &
                                       partZone,xmin,ymin,zmin,nzx,nzy,nzz,nztotal)
Implicit none

  Integer,value :: nf,nwall
  Integer,value :: nzx,nzy,nzz,nztotal
  Real(8),value :: mass,h,wBeta,rho0,eta2,dt
  Real(8),value :: xmin,ymin,zmin
  Real(8) :: grav(3)

  Real(8) :: xyz_f(3,nf),vel(3,nf),veddy(nf),visc(nf),a(3*nf)
  Real(8) :: xyz_w(3,nwall),vel_w(3,nwall)

  Real(8) :: dist,R
  Real(8) :: wFunc_0,wFunc_x,wFunc_y,wFunc_z
  Real(8) :: vi,rxW,Bij,SB_x,SB_y,SB_z

  Integer :: partZone(nztotal,4)
  Integer :: bound(nf)

  Integer :: zi,zf,yi,yf,xi,xf,x,y,z,loop
  Integer :: xNow,yNow,zNow,zoneJ
  Integer :: i,j,k

  i=(blockIdx%x-1)*blockDim%x+threadIdx%x  !Matrix row
  if(i>nf) goto 999  
  if(bound(i)==7) then
    a(3*(i-1)+1)=vel(1,i)+dt*(grav(1))
    a(3*(i-1)+2)=vel(2,i)+dt*(grav(2))
    a(3*(i-1)+3)=vel(3,i)+dt*(grav(3))
    goto 999
  endif

  SB_x=0d0; SB_y=0d0; SB_z=0d0

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
                R=dist/h
                if(R<1d0) wFunc_0=-3d0*R+(9d0/4d0)*(R**2)
                if(R>=1d0) wFunc_0=(-3d0/4d0)*((2d0-R)**2)
                wFunc_0=wFunc_0*wBeta/h
                wFunc_x=wFunc_0*(xyz_f(1,i)-xyz_f(1,j))/dist
                wFunc_y=wFunc_0*(xyz_f(2,i)-xyz_f(2,j))/dist
                wFunc_z=wFunc_0*(xyz_f(3,i)-xyz_f(3,j))/dist
                vi=((visc(i)+veddy(i))+(visc(j)+veddy(j)))/(rho0)
                rxW=(xyz_f(1,i)-xyz_f(1,j))*wFunc_x+ &
                    (xyz_f(2,i)-xyz_f(2,j))*wFunc_y+ &
                    (xyz_f(3,i)-xyz_f(3,j))*wFunc_z
                Bij=mass*vi*rxW/(dist**2+eta2)
                SB_x=SB_x+Bij*(vel(1,i)-vel(1,j))
                SB_y=SB_y+Bij*(vel(2,i)-vel(2,j))
                SB_z=SB_z+Bij*(vel(3,i)-vel(3,j))
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
            if(dist<2e0*h) then
              R=dist/h
              if(R<1d0) wFunc_0=-3d0*R+(9d0/4d0)*(R**2)
              if(R>=1d0) wFunc_0=(-3d0/4d0)*((2d0-R)**2)
              wFunc_0=wFunc_0*wBeta/h
              wFunc_x=wFunc_0*(xyz_f(1,i)-xyz_w(1,j))/dist
              wFunc_y=wFunc_0*(xyz_f(2,i)-xyz_w(2,j))/dist
              wFunc_z=wFunc_0*(xyz_f(3,i)-xyz_w(3,j))/dist
              vi=2d0*(visc(i)+veddy(i))/(rho0)
              rxW=(xyz_f(1,i)-xyz_w(1,j))*wFunc_x+ &
                  (xyz_f(2,i)-xyz_w(2,j))*wFunc_y+ &
                  (xyz_f(3,i)-xyz_w(3,j))*wFunc_z
              Bij=mass*vi*rxW/(dist**2+eta2)
              SB_x=SB_x+Bij*(vel(1,i)-vel_w(1,j))
              SB_y=SB_y+Bij*(vel(2,i)-vel_w(2,j))
              SB_z=SB_z+Bij*(vel(3,i)-vel_w(3,j))
            endif
          enddo
        endif
      enddo
    enddo
  enddo

  a(3*(i-1)+1)=vel(1,i)+dt*(SB_x+grav(1))
  a(3*(i-1)+2)=vel(2,i)+dt*(SB_y+grav(2))
  a(3*(i-1)+3)=vel(3,i)+dt*(SB_z+grav(3))

999 continue
return
end subroutine vstarEXP



end module gpu_vstar


