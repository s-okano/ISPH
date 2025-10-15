program main
!Libraries
use cudafor  !Cuda fortran
use thrust   !Thrust: sorting, prefixSum

!Modules to define the GPU kernels
use gpu_vstar
use gpu_corrector

!Modules to define the most important variables
use properties      !CPU variables
use parameters      !CPU parameters (constants)
use gpu_properties  !GPU variables

Implicit none

  Integer :: i,j,k,maxiter,max_i
  Integer :: istat,status
  Integer :: dum1,dum2,step
  Integer :: aaaa
  
  Real :: time

  type(cudaEvent) :: start_e,stop_e

  Real(8) :: conv
  Real(8) :: pre


  write(*,*) " "
  write(*,*) '==============================================='
  write(*,*) '   =============3D SPH on GPU=============='
  write(*,*) '==============================================='

  call input

  istat=cudaDeviceSynchronize
  write(*,*)istat

  call output()
  conv=1d-12

  max_i=kwrite*iter

  write(*,*) " "
  write(*,*) " "
  write(*,*) "Main Properties"
  write(*,*) " "
  write(*,"(A21,I15)"     ) "Device:              ",dev
  write(*,"(A21,I15)"     ) "kwrite:              ",kwrite
  write(*,"(A21,f15.3,A7)") "Max simulation time: ",dt*max_i,"s"
  write(*,"(A21,f15.9,A7)") "Particle diameter:   ",d,"m"
  write(*,"(A21,f15.9   )") "Alfa:                ",alfa
  write(*,"(A21,f15.9,A7)") "dt:                  ",dt,"s"
  write(*,"(A21,f15.9,A7)") "Rho:                 ",rho0,"kg/m3"
  write(*,"(A21,f15.9,A7)") "Viscosity 0:         ",visc0,"m2/s"
  write(*,"(A21,f15.9,A7)") "Viscosity max:       ",viscMax,"m2/s"
!  write(*,"(A21,f15.9,A7)") "Cohesion:            ",coes,"Pa"
!  write(*,"(A21,f15.9,A7)") "Friction angle:      ",phi,"rad"
!  write(*,"(A21,f15.9,A7)") "Yield stress:        ",yield,"Pa"

  istat=cudaDeviceSynchronize
  istat=cudaEventCreate(start_e)
  istat=cudaEventCreate(stop_e)
  istat=cudaEventRecord(start_e,0)
  istat=cudaDeviceSynchronize

  call sorting_wall
  call tracking_wall

  step=0

  do k=1,max_i
    write(*,*) '.'
    write(*,*) '===============New iteration==============='
    write(*,"(a11,I10)") 'Particles= ', n
    write(*,"(a5,I10)") 'step= ', k
    write(*,"(a5,f20.9)") 't= ', t

    write(*,*) 'Tracking...'
    call sorting
    call tracking_fluid
    write(*,*) 'End tracking.'


    write(*,*) 'Vstar...'
    call Init<<<bl_fluid,Threads>>>(nf,xyz_f_d,vel_f_d,veddy_f_d,visc_f_d,rho_f_d,bound_f_d,&
                                    pres_f_d,&
                                    nwall,xyz_w_d,vel_w_d,&
                                    mass,h,wBeta,rho0,eta2,&
                                    ne_d,&
                                    visc0,viscMax,coes,phi,&
                                    partZone_d,xmin,ymin,zmin,nzx,nzy,nzz,nztotal)

    call vstarEXP<<<bl_fluid,Threads>>>(nf,xyz_f_d,vel_f_d,veddy_f_d,visc_f_d,a_d,bound_f_d,&
                                        nwall,xyz_w_d,vel_w_d,&
                                        mass,h,wBeta,rho0,eta2,dt,grav_d,&
                                        partZone_d,xmin,ymin,zmin,nzx,nzy,nzz,nztotal)

    call vel_up<<<bl_fluid,Threads>>>(nf,nf,vel_f_d,a_d,bound_f_d,grav_d,dt)

    write(*,*) 'End of Vstar.'


    write(*,*) 'Corrector...'

    if(time_int==2) then
      istat=cudaDeviceSynchronize
      call prefixSum(ne_d,mRow_d,nf+1)

      istat=cudaDeviceSynchronize
      aaaa=mRow_d(nf+1)
      istat=cudaDeviceSynchronize
      aaaa=aaaa-1

      call virt1<<<bl_wall,Threads>>>(nf,xyz_f_d,vel_f_d,visc_f_d,veddy_f_d,&
                                      nwall,xyz_w_d,xyzbound_w_d,vel_w_d,pres_w_d,&
                                      mass,h,wBeta,rho0,eta2,grav_d,&
                                      partZone_d,xmin,ymin,zmin,nzx,nzy,nzz,nztotal,phi,pi)

      call defMatPres<<<bl_fluid,Threads>>>(nf,xyz_f_d,vel_f_d,rho_f_d,pres_f_d,&
                                            bound_f_d,M_d,a_d,b_d,mCol_d,mRow_d, &
                                            nwall,xyz_w_d,pres_w_d,vel_w_d, &
                                            mass,h,wBeta,rho0,eta2,dt,alfa,maxne, &
                                            partZone_d,xmin,ymin,zmin,nzx,nzy,nzz,nztotal)

      call iccg(conv,20,nf,aaaa,bl_fluid)
                                                           
!      call iccg(conv,20,nf,maxne,bl_fluid)

      call pres_up<<<bl_fluid,Threads>>>(nf,nf,pres_f_d,a_d,bound_f_d)

      call pres_virt<<<bl_wall,Threads>>>(nf,xyz_f_d,pres_f_d,&
                                          nwall,xyz_w_d,xyzbound_w_d,pres_w_d,&
                                          mass,h,wBeta,rho0,&
                                          partZone_d,xmin,ymin,zmin,nzx,nzy,nzz,nztotal)

    elseif(time_int==1) then
      call PPEexp<<<bl_fluid,Threads>>>(nf,xyz_f_d,vel_f_d,rho_f_d,pres_f_d,a_d,bound_f_d,&
                                        nwall,xyz_w_d,pres_w_d,vel_w_d,&
                                        mass,h,wBeta,rho0,eta2,dt,alfa,&
                                        partZone_d,xmin,ymin,zmin,nzx,nzy,nzz,nztotal)

      call pres_up<<<bl_fluid,Threads>>>(nf,nf,pres_f_d,a_d,bound_f_d)

      call virt1<<<bl_wall,Threads>>>(nf,xyz_f_d,vel_f_d,visc_f_d,veddy_f_d,&
                                      nwall,xyz_w_d,xyzbound_w_d,vel_w_d,pres_w_d,&
                                      mass,h,wBeta,rho0,eta2,grav_d,&
                                      partZone_d,xmin,ymin,zmin,nzx,nzy,nzz,nztotal,phi,pi)

      call pres_virt<<<bl_wall,Threads>>>(nf,xyz_f_d,pres_f_d,&
                                          nwall,xyz_w_d,xyzbound_w_d,pres_w_d,&
                                          mass,h,wBeta,rho0,&
                                          partZone_d,xmin,ymin,zmin,nzx,nzy,nzz,nztotal)
    endif

    call corrector<<<bl_fluid,Threads>>>(nf,xyz_f_d,vel_f_d,pres_f_d,bound_f_d,&
                                         nwall,xyz_w_d,xyzbound_w_d,pres_w_d,&
                                         mass,h,wBeta,rho0,dt,&
                                         partZone_d,xmin,ymin,zmin,nzx,nzy,nzz,nztotal)

    call x_up<<<bl_fluid,Threads>>>(nf,xyz_f_d,vel_f_d,bound_f_d,&
                                    nwall,xyz_w_d,xyzbound_w_d,h,dt,&
                                    partZone_d,xmin,ymin,zmin,nzx,nzy,nzz,nztotal)
    write(*,*) 'End of the corrector.'


    t=t+dt

    if(mod(k,kwrite)==0) then
      write(*,*) 'Output...'
      call output()
      write(*,*) 'End of output.'
    endif


  enddo

  istat=cudaDeviceSynchronize
  istat=cudaEventRecord(stop_e,0)
  istat=cudaEventSynchronize(stop_e)
  istat=cudaEventElapsedTime(time,start_e,stop_e)
  istat=cudaDeviceSynchronize
  write(*,*) "time: ",time,"ms"


  !!write(*,*) 'Output...'
  !!call output()
  !!write(*,*) 'End of output.'


  write(*,*) '==============================================='
  write(*,*) '   =============3D SPH on GPU=============='
  write(*,*) '==============================================='
  write(*,*) 'Program successfully finished'
  write(*,*) ' '

endprogram

