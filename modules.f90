module properties
  Implicit none

  Integer :: n,nf,nwall
  Real    :: t
  Integer :: Threads,Blocks,dev
  Integer :: bl_fluid,bl_wall,bl_zones
  Integer :: kwrite,iter

  Real(8),allocatable ::      xyz_f(:,:)
  Real(8),allocatable ::      vel_f(:,:)
  Real(8),allocatable ::     pres_f(:)
  Real(8),allocatable ::      rho_f(:)
  Real(8),allocatable ::    veddy_f(:)
  Real(8),allocatable ::     visc_f(:)
  Integer,allocatable ::    bound_f(:)

  Real(8),allocatable ::      xyz_w(:,:)
  Real(8),allocatable ::      vel_w(:,:)
  Real(8),allocatable :: xyzbound_w(:,:)
  Real(8),allocatable ::     pres_w(:)
  Real(8),allocatable ::      rho_w(:)
  Integer,allocatable ::    bound_w(:)

  Real(8),save        :: xmin,ymin,zmin
  Real(8),save        :: xmax,ymax,zmax
  Integer             :: nzx,nzy,nzz

  Integer :: time_int ! =1: fully-explicit
                      ! =2: semi-implicit

  Integer :: maxne=150

end module properties

module parameters

  Real(8),parameter :: grav_x=0d0          ! m/s2
  Real(8),parameter :: grav_y=-9.8d0       ! m/s2
  Real(8),parameter :: grav_z=0d0          ! m/s2
  Real(8),parameter :: pi=4d0*atan(1d0)

  Real(8) :: d
  Real(8) :: h
  Real(8) :: rho0
  Real(8) :: mass
  Real(8) :: eta2
  Real(8) :: visc0,viscMax         ! m2/s 
  Real(8) :: dt                    ! s
  Real(8) :: alfa
  Real(8) :: wBeta
  Real(8) :: yield     ! Pa
  Real(8) :: phi,coes

  Integer :: nztotal

end module parameters

module gpu_properties

  Integer,allocatable,device ::    bound_f_d(:)
  Real(8),allocatable,device ::      xyz_f_d(:,:)
  Real(8),allocatable,device ::      vel_f_d(:,:)
  Real(8),allocatable,device ::     pres_f_d(:)
  Real(8),allocatable,device ::      rho_f_d(:)
  Real(8),allocatable,device ::    veddy_f_d(:)
  Real(8),allocatable,device ::     visc_f_d(:)

  Real(8),allocatable,device ::      xyz_w_d(:,:)
  Real(8),allocatable,device ::      vel_w_d(:,:)
  Real(8),allocatable,device :: xyzbound_w_d(:,:)
  Real(8),allocatable,device ::     pres_w_d(:)
  Real(8),allocatable,device ::      rho_w_d(:)

  Real(8),allocatable,device :: M_d(:),a_d(:),b_d(:)
  Integer,allocatable,device :: mCol_d(:)
  Integer,allocatable,device :: mRow_d(:)
  Integer,allocatable,device :: ne_d(:)

  Real(8),allocatable,device :: r_d(:),p_d(:),q_d(:)

  Real(8),allocatable,device :: grav_d(:)
  Integer,allocatable,device :: partZone_d(:,:)

  Real(8),allocatable,device :: key_d(:)
  Integer,allocatable,device :: newID_d(:)

end module gpu_properties


