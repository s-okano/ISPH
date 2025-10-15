subroutine input
use cudafor
use properties
use parameters
use gpu_properties

use thrust
use gpu_sorting

  Implicit none
  Integer :: i,j,k,l,nul
  Integer :: thr,istat
  Integer :: dum_i
  Real(8) :: dist,xb,yb,zb
  Real(8) :: dum_r(10)


  open(unit=700, file='data_input.txt', access='sequential')
  read(700,*) dev,kwrite,iter,time_int
  read(700,*) d,dt,alfa
  read(700,*) rho0,visc0,viscMax
  read(700,*) coes,phi,yield
  close(700)
  istat=cudaSetDevice(dev)

  open(unit=300, file='xyzWater.txt', access='sequential')
  open(unit=400, file='xyzWall.txt', access='sequential')

  read(300,*) nf
  read(400,*) nwall

  write(*,*) "# of fluid particles: ", nf
  write(*,*) "# of wall  particles: ", nwall
  write(*,*) "Total # of particles: ", nf+nwall

  !CPU============================
  allocate(      xyz_f(3,nf))
  allocate(      vel_f(3,nf))
  allocate(     pres_f(nf))
  allocate(      rho_f(nf))
  allocate(     visc_f(nf))
  allocate(    bound_f(nf))
  allocate(    veddy_f(nf))

  allocate(      xyz_w(3,nwall))
  allocate(      vel_w(3,nwall))
  allocate( xyzbound_w(3,nwall))
  allocate(     pres_w(nwall))
  !CPU============================

  xyz_f =0d0
  vel_f =0d0
  pres_f=0d0
  rho_f =0d0
  visc_f=0d0

  pres_w=0d0

  h=1.2d0*d           ! m 
  mass=(d**3)*rho0  ! kg
  eta2=1d-4*h*h 



  wBeta=1d0/(pi*(h**3))
  wBeta=wBeta*1000d0/1000.80954839d0


  phi=phi*pi/180e0

  n=max(nf,nwall)
  write(*,*) "  Inputing CPU data..."
  do i=1,n
    if(i<=nf) then
      read(300,*) nul,dum_i,dum_r(1),dum_r(2),dum_r(3),dum_r(4),dum_r(5),dum_r(6),&
                  dum_r(7),dum_r(8),dum_r(9),dum_r(10)
  
      xyz_f(1,i)=dum_r(1)
      xyz_f(2,i)=dum_r(2)
      xyz_f(3,i)=dum_r(3)
      vel_f(1,i)=dum_r(4)
      vel_f(2,i)=dum_r(5)
      vel_f(3,i)=dum_r(6)

      if(i==1) then
        xmin=dum_r(1)
        xmax=dum_r(1)
        ymin=dum_r(2)
        ymax=dum_r(2)
        zmin=dum_r(3)
        zmax=dum_r(3)
      endif

      if(dum_r(1)<xmin) xmin=dum_r(1)
      if(dum_r(1)>xmax) xmax=dum_r(1)
      if(dum_r(2)<ymin) ymin=dum_r(2)
      if(dum_r(2)>ymax) ymax=dum_r(2)
      if(dum_r(3)<zmin) zmin=dum_r(3)
      if(dum_r(3)>zmax) zmax=dum_r(3)

    endif
    if(i<=nwall) then
      read(400,*) nul,dum_i,dum_r(1),dum_r(2),dum_r(3),dum_r(4),dum_r(5),dum_r(6),&
                  dum_r(7),dum_r(8),dum_r(9),dum_r(10)

      xyz_w(1,i)=dum_r(1)
      xyz_w(2,i)=dum_r(2)
      xyz_w(3,i)=dum_r(3)
      vel_w(1,i)=dum_r(4)
      vel_w(2,i)=dum_r(5)
      vel_w(3,i)=dum_r(6)

      xyzbound_w(1,i)=-dum_r(7)*dum_r(8)
      xyzbound_w(2,i)=-dum_r(7)*dum_r(9)
      xyzbound_w(3,i)=-dum_r(7)*dum_r(10)

      if(dum_r(1)<xmin) xmin=dum_r(1)
      if(dum_r(1)>xmax) xmax=dum_r(1)
      if(dum_r(2)<ymin) ymin=dum_r(2)
      if(dum_r(2)>ymax) ymax=dum_r(2)
      if(dum_r(3)<zmin) zmin=dum_r(3)
      if(dum_r(3)>zmax) zmax=dum_r(3)
    endif

  enddo
  write(*,*) "  End of Inputing CPU data"
  close(300)
  close(400)

  write(*,*) " "

  Threads=128
  if(mod(nf,Threads)==0) bl_fluid=int(nf/Threads)
  if(mod(nf,Threads)/=0) bl_fluid=int(nf/Threads)+1

  if(mod(nwall,Threads)==0) bl_wall=int((nwall)/Threads)
  if(mod(nwall,Threads)/=0) bl_wall=int((nwall)/Threads)+1
  write(*,*) "Blocks: ", bl_fluid,bl_wall
  write(*,*) "Threads: ", Threads

  xmin=xmin-5*d
  xmax=xmax+5*d
  ymin=ymin-5*d
  ymax=ymax+500*d
  zmin=zmin-5*d
  zmax=zmax+5*d

  write(*,"(A4,F10.5,A1,F10.5,A1)") "x= {",xmin,",",xmax,"}"
  write(*,"(A4,F10.5,A1,F10.5,A1)") "y= {",ymin,",",ymax,"}"
  write(*,"(A4,F10.5,A1,F10.5,A1)") "z= {",zmin,",",zmax,"}"

  nzx=int((xmax-xmin)/(2*h))+1
  nzy=int((ymax-ymin)/(2*h))+1
  nzz=int((zmax-zmin)/(2*h))+1
  nztotal=nzx*nzy*nzz


  if(mod(nztotal,Threads)==0) bl_zones=int(nztotal/Threads)
  if(mod(nztotal,Threads)/=0) bl_zones=int(nztotal/Threads)+1

  write(*,*) "  GPU allocation..."
  !GPU============================
  allocate(      xyz_f_d(3,nf))
  allocate(      vel_f_d(3,nf))
  allocate(     pres_f_d(nf))
  allocate(      rho_f_d(nf))
  allocate(    veddy_f_d(nf))
  allocate(     visc_f_d(nf))
  allocate(    bound_f_d(nf))

  allocate(      xyz_w_d(3,nwall))
  allocate(      vel_w_d(3,nwall))
  allocate( xyzbound_w_d(3,nwall))
  allocate(     pres_w_d(nwall))

  allocate(ne_d(nf+1))
  ne_d=0
  if(time_int==2) then
    allocate(mRow_d(nf+1))
  
    allocate(M_d(maxne*nf))
    allocate(b_d(nf))
    allocate(mCol_d(maxne*nf))
  
    allocate(r_d(nf))
    allocate(p_d(nf))
    allocate(q_d(nf))

    mRow_d=0
    M_d=0
    b_d=0
    mCol_d=0
    r_d=0
    p_d=0
    q_d=0

  endif

  !GPU============================
  allocate(partZone_d(nztotal,4))

  allocate(key_d(n))
  allocate(newID_d(n))
  allocate(a_d(3*n))
               
  allocate(grav_d(3))

  write(*,*) "  End of GPU allocation."

  istat=cudaDeviceSynchronize
  xyz_w_d=xyz_w 
  xyzbound_w_d=xyzbound_w
  vel_w_d=vel_w

  xyz_f_d=xyz_f
  vel_f_d=vel_f
  istat=cudaDeviceSynchronize

  !Initialization
  pres_f_d=0
  rho_f_d=0
  visc_f_d=visc0
  veddy_f_d=0
  bound_f_d=1

  pres_w_d=0

  partZone_d=0

  grav_d(1)=grav_x
  grav_d(2)=grav_y
  grav_d(3)=grav_z

  !GPU INPUT=======================================================

return
endsubroutine



subroutine init_random_seed()

Integer :: i,n,clock
Integer,dimension(:),allocatable :: seed

  call random_seed(size=n)
  allocate(seed(n))

  call system_clock(count=clock)

  seed=clock+37*(/(i-1,i=1,n)/)
  call random_seed(put=seed)

  deallocate(seed)

end subroutine















