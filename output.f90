subroutine  output()
  use cudafor
  use properties
  use parameters
  use gpu_properties
  use vtk_module

  implicit none
  integer,save :: stepNo = 0
  character(len=256) :: name
  integer :: len_name,istat

  Real :: vector_single(3,n)
  Real :: scalar_single(n)

  istat=cudaDeviceSynchronize
 xyz_f(:,1:nf)=xyz_f_d(:,1:nf)
 vel_f(:,1:nf)=vel_f_d(:,1:nf)
  pres_f(1:nf)= pres_f_d(1:nf)
   rho_f(1:nf)=  rho_f_d(1:nf)
  visc_f(1:nf)= visc_f_d(1:nf)
 veddy_f(1:nf)=veddy_f_d(1:nf)
 bound_f(1:nf)=bound_f_d(1:nf)
  istat=cudaDeviceSynchronize

  call vtk_initialize()
  vector_single(:,1:nf)=real(xyz_f(:,1:nf))
  call vtk_setPoints(nf,vector_single)

  vector_single(:,1:nf)=real(vel_f(:,1:nf))
  name = 'Velocity (m/s)'
  call vtk_setPointVectorValue(nf,vector_single, len_trim(name), trim(name))

  scalar_single(1:nf)=real(pres_f(1:nf))
  name = 'Pressure (Pa)'
  call vtk_setPointScalarValue(nf,scalar_single, len_trim(name), trim(name))

  scalar_single(1:nf)=real(visc_f(1:nf))
  name = 'Viscosity (m2/s)'
  call vtk_setPointScalarValue(nf,scalar_single, len_trim(name), trim(name))

  scalar_single(1:nf)=real(veddy_f(1:nf))
  name = 'Eddy viscosity (m2/s)'
  call vtk_setPointScalarValue(nf,scalar_single, len_trim(name), trim(name))

  scalar_single(1:nf)=real(rho_f(1:nf))
  name = 'Density (kg/m3)'
  call vtk_setPointScalarValue(nf,scalar_single, len_trim(name), trim(name))

  name = 'Boundary condition'
  call vtk_setPointScalarValue(nf,bound_f, len_trim(name), trim(name))

  write(name,'("resultWater_",I5.5,".vtu")') stepNo
  call vtk_output(len_trim(name), trim(name))
  call vtk_finalize()



  istat=cudaDeviceSynchronize
   xyz_w=     xyz_w_d
   vel_w=     vel_w_d
  pres_w=    pres_w_d
  istat=cudaDeviceSynchronize

  vector_single(:,1:nwall)=real(xyz_w(:,1:nwall))
  call vtk_initialize()
  call vtk_setPoints(nwall,vector_single)

  vector_single(:,1:nwall)=real(vel_w(:,1:nwall))
  name = 'Velocity (m/s)'
  call vtk_setPointVectorValue(nwall,vector_single, len_trim(name), trim(name))

  scalar_single(1:nwall)=real(pres_w( 1:nwall))
  name = 'Pressure (Pa)'
  call vtk_setPointScalarValue(nwall,scalar_single, len_trim(name), trim(name))

  write(name,'("resultWall_",I5.5,".vtu")') stepNo
  call vtk_output(len_trim(name), trim(name))
  call vtk_finalize()



  stepNo = stepNo +1

end subroutine output
