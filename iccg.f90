module gpu_iccg

contains
attributes(global) subroutine finalSum(partial,total)
Implicit none

  Real(8) :: partial(*)
  Real(8) :: total
  Real(8),shared :: psum(256)
  Integer :: t,inext

  t=threadIdx%x

  psum(t)=partial(t)
  call syncthreads()

  inext=blockDim%x/2
  do while(inext>=1)
    if(t<=inext) psum(t)=psum(t)+psum(t+inext)
    inext=inext/2
    call syncthreads()
  enddo

  if(t==1) total=psum(1)

end subroutine finalSum


attributes(global) subroutine partialSum(vec1,vec2,partial,n)
Implicit none

  Real(8) :: vec1(*),vec2(*)
  Real(8) :: partial(*)
  Real(8),shared :: psum(512)
  Integer, value :: n,blockSize
  Integer :: t,loc_t,inext

  Integer :: i,elemPt

  t=threadIdx%x+(blockIdx%x-1)*blockDim%x
  loc_t=threadIdx%x
  psum(loc_t)=0d0

  elemPt=(n/gridDim%x)/blockDim%x+1

  do i=(t-1)*elemPt+1,t*elemPt
    if(i<=n) psum(loc_t)=psum(loc_t)+vec1(i)*vec2(i)
  enddo

  call syncthreads()

  inext=blockDim%x/2
  do while(inext>=1)
    if(loc_t<=inext) psum(loc_t)=psum(loc_t)+psum(loc_t+inext)
    inext=inext/2
    call syncthreads()
  enddo
  if(loc_t==1) partial(blockIdx%x)=psum(1)

999 continue
end subroutine partialSum


attributes(global) subroutine p_update(n,beta,r,p)
Implicit none

  Integer,value :: n
  Real(8),value :: beta
  Real(8) :: r(n),p(n)

  Integer :: i

  i=(blockIdx%x-1)*blockDim%x+threadIdx%x
  if(i>n) goto 999  

  p(i)=r(i)+beta*p(i)

999 continue
return
end subroutine p_update


attributes(global) subroutine a_r_update(n,alfa,a,r,p,q)
Implicit none

  Integer,value :: n
  Real(8) :: a(n)
  Real(8) :: r(n),p(n),q(n)

  Real(8),value :: alfa
  Integer :: i

  i=(blockIdx%x-1)*blockDim%x+threadIdx%x
  if(i>n) goto 999  
  
  a(i)=a(i)+alfa*p(i)
  r(i)=r(i)-alfa*q(i)

999 continue
return
end subroutine a_r_update





attributes(global) subroutine s_updateBI(n,s,r,Ap,alfa)
Implicit none

  Integer,value :: n
  Real(8),value :: alfa
  Real(8) :: s(n),r(n),Ap(n)

  Integer :: i

  i=(blockIdx%x-1)*blockDim%x+threadIdx%x
  if(i>n) goto 999  

  s(i)=r(i)-alfa*Ap(i)

999 continue
return
end subroutine s_updateBI


attributes(global) subroutine p_updateBI(n,p,r,Ap,beta,omega)
Implicit none

  Integer,value :: n
  Real(8),value :: beta,omega
  Real(8) :: p(n),r(n),Ap(n)

  Integer :: i

  i=(blockIdx%x-1)*blockDim%x+threadIdx%x
  if(i>n) goto 999  

  p(i)=r(i)+beta*(p(i)-omega*Ap(i))

999 continue
return
end subroutine p_updateBI




attributes(global) subroutine a_r_updateBI(n,a,r,p,s,As,alfa,omega)
Implicit none

  Integer,value :: n
  Real(8),value :: alfa,omega
  Real(8) :: a(n),r(n),p(n),s(n),As(n)

  Integer :: i

  i=(blockIdx%x-1)*blockDim%x+threadIdx%x
  if(i>n) goto 999

  a(i)=a(i)+alfa*p(i)+omega*s(i)
  r(i)=s(i)-omega*As(i)

999 continue
return
end subroutine a_r_updateBI















end module gpu_iccg


subroutine iccg(conv,mo,nnn,non_zero,bl_iccg)
use cudafor
use cusparse
use cublas

use properties
use parameters
use gpu_properties

use gpu_iccg

Implicit none

  Real(8) :: conv
  Integer :: i,p,k,nnn,non_zero,bl_iccg
  Integer :: mo,istat

  Real(8),device :: temp
  Real(8) :: beta,alfa0,pq,rho_0,rho_1,err1

  type(cusparseHandle) :: handle
  Integer :: status
  Real(8) :: alpha2,beta2

  Real(8),device :: part(512) !!!!!! FOR THE partialSum,finalSum subroutines


  type(cusparseSpMatDescr) :: matM
  type(cusparseDnVecDescr) :: vecA,vecR,vecP,vecQ
  Integer(8) :: bufferSize

  Integer(8),save :: bufferSize_alloc = 0
  Integer(1),pointer,device :: buffer(:)

  ! initalize CUSPARSE and matrix descriptor
  status=cusparseCreate(handle)
  if(status/=CUSPARSE_STATUS_SUCCESS) write(*,*) 'cusparseCreate error: ', status


  err1=1d3
  alpha2=-1d0
  beta2=1d0
  write(*,*) "  ICCG method"
  write(*,*) "        Step   Error"

  istat=cudaDeviceSynchronize
  istat=cudaMemcpy(r_d(1:nnn),b_d(1:nnn),nnn)
  istat=cudaDeviceSynchronize

  status=cusparseCreateCsr(matM,nnn,nnn,non_zero, &  !mat(cusparse), rows, cols, non zero
                           mRow_d,mCol_d,M_d,  &
                           CUSPARSE_INDEX_32I,CUSPARSE_INDEX_32I, &
                           CUSPARSE_INDEX_BASE_ONE,CUDA_R_64F)

  status=cusparseCreateDnVec(vecA,nnn,a_d,CUDA_R_64F)  ! in

  status=cusparseCreateDnVec(vecR,nnn,r_d,CUDA_R_64F)   ! in

  status=cusparseSpMV_buffersize(handle,CUSPARSE_OPERATION_NON_TRANSPOSE, & ! in   op
                                 alpha2,matM,vecA,  & ! in                         alph A, X
                                 beta2,vecR,  & ! in/out                           beta Y
                                 CUDA_R_64F,CUSPARSE_CSRMV_ALG1,bufferSize)   ! in
  
  istat=cudaDeviceSynchronize



  if(allocated(buffer)) deallocate(buffer)
  if( bufferSize == 0 ) then
    nullify(buffer)
  else
    allocate(buffer(bufferSize))
  endif
  

  !! Y = alph op(A) * X + beta Y
  status=cusparseSpMV(handle,CUSPARSE_OPERATION_NON_TRANSPOSE, & ! in   op
                      alpha2,matM,vecA,  & ! in                         alph A, X
                      beta2,vecR,  & ! in/out                           beta Y
                      CUDA_R_64F,CUSPARSE_CSRMV_ALG1,buffer)   ! in
  

  istat=cudaDeviceSynchronize
  istat=cudaMemcpy(p_d(1:nnn),r_d(1:nnn),nnn)
  istat=cudaDeviceSynchronize


  status=cusparseCreateDnVec(vecP,  &  ! out
                             nnn,  &  ! in
                             p_d,  &  ! in
                             CUDA_R_64F)   ! in

  status=cusparseCreateDnVec(vecQ,  &  ! out
                             nnn,  &  ! in
                             q_d,  &  ! in
                             CUDA_R_64F)   ! in
  istat=cudaDeviceSynchronize
  

  alpha2=1d0
  beta2=0d0
  k=1
  do while(err1>conv)  !=========================================


    temp=0d0
    istat=cudaDeviceSynchronize
    call partialSum<<<256,512,512*8>>>(r_d,r_d,part,nnn)
    call finalSum<<<1,256,256*8>>>(part,temp)
    istat=cudaDeviceSynchronize

    rho_1=temp

    if(rho_1<conv) then
      err1=rho_1
      goto 123
    endif



    istat=cudaDeviceSynchronize
    if(k>1) then
      beta=rho_1/rho_0
      call p_update<<<bl_iccg,Threads>>>(nnn,beta,r_d,p_d)
    endif


    !! Y = alph op(A) * X + beta Y
    status=cusparseSpMV(handle,CUSPARSE_OPERATION_NON_TRANSPOSE, & ! in   op
                        alpha2,matM,vecP,  & ! in                         alph A, X
                        beta2,vecQ,  & ! in/out                           beta Y
                        CUDA_R_64F,CUSPARSE_CSRMV_ALG1,buffer)   ! in


    istat=cudaDeviceSynchronize
    err1=rho_1

    if(mod(k,mo)==0) write(*,*) k,err1
    if(err1>=0) then
    else
      write(*,*) k,err1," Solution diverged =(!"
!      call output()
      stop
    endif
    if(k>5000) then
      write(*,*) k,err1," Not converging! =("
!      call output()
      stop
    endif

    rho_0=rho_1
    k=k+1


    temp=0d0
    istat=cudaDeviceSynchronize
    call partialSum<<<256,512,512*8>>>(p_d,q_d,part,nnn)
    call finalSum<<<1,256,256*8>>>(part,temp)
    istat=cudaDeviceSynchronize

    pq=temp
    istat=cudaDeviceSynchronize
    alfa0=rho_1/pq

    call a_r_update<<<bl_iccg,Threads>>>(nnn,alfa0,&
                                      a_d,r_d,p_d,q_d)
    
  enddo !=======================================================
  123 continue

  write(*,*) k,err1," Converged!"

return
end subroutine iccg


















subroutine biccg_stab(conv,mo,nnn,non_zero,bl_bi)!,steps)
use cudafor
use cusparse

use properties
use parameters
use gpu_properties

use gpu_iccg

Implicit none

  Integer :: mo,nnn,non_zero,bl_bi
  Real(8) :: conv
  Integer :: k,steps
  Real(8) :: err1,err0
  Real(8),device :: temp
  Integer :: istat
  
  Integer :: i,j
 
  Real(8),device,allocatable :: Ap_d(:),rr0_d(:)
  Real(8),device,allocatable :: s_d(:),dumR_d(:),As_d(:)
  Real(8),device,allocatable :: part(:) !!!!!! FOR THE partialSum,finalSum subroutines
  Integer(1),pointer,device :: buffer(:)

  Real(8) :: alfa_d,beta_d,omega_d,dum1,dum2,old_r_rr

  type(cusparseHandle) :: handle
  Integer :: status
  Real(8) :: alpha2,beta2

  type(cusparseSpMatDescr) :: matM
  type(cusparseDnVecDescr) :: vecA,vecR,vecAp,vecAs,vecdumR,vecP,vecS
  Integer(8) :: bufferSize


  Real(8),allocatable :: rand(:)


  allocate(Ap_d(nnn),rr0_d(nnn))
  allocate(s_d(nnn),dumR_d(nnn),As_d(nnn),part(512))
  allocate(rand(nnn))

        


  ! initalize CUSPARSE and matrix descriptor
  status=cusparseCreate(handle)
  if(status/=CUSPARSE_STATUS_SUCCESS) write(*,*) 'cusparseCreate error: ', status


  istat=cudaDeviceSynchronize
  alpha2=-1d0
  beta2=1d0
  istat=cudaMemcpy(r_d(1:nnn),b_d(1:nnn),nnn)
  istat=cudaDeviceSynchronize


  status=cusparseCreateCsr(matM,nnn,nnn,non_zero, &  !mat(cusparse), rows, cols, non zero
                           mRow_d,mCol_d,M_d,  &
                           CUSPARSE_INDEX_32I,CUSPARSE_INDEX_32I, &
                           CUSPARSE_INDEX_BASE_ONE,CUDA_R_64F)

  status=cusparseCreateDnVec(vecA,nnn,a_d,CUDA_R_64F)  ! in

  status=cusparseCreateDnVec(vecR,nnn,r_d,CUDA_R_64F)   ! in

  status=cusparseSpMV_buffersize(handle,CUSPARSE_OPERATION_NON_TRANSPOSE, & ! in   op
                                 alpha2,matM,vecA,  & ! in                         alph A, X
                                 beta2,vecR,  & ! in/out                           beta Y
                                 CUDA_R_64F,CUSPARSE_CSRMV_ALG1,bufferSize)   ! in
  istat=cudaDeviceSynchronize

  if(allocated(buffer)) deallocate(buffer)
  if( bufferSize == 0 ) then
    bufferSize=1
          !nullify(buffer)
  !else
  endif
  allocate(buffer(bufferSize))



  write(*,*) "  BICCG STAB method"
  write(*,*) "        Step   Error"


  !! Y = alph op(A) * X + beta Y
  status=cusparseSpMV(handle,CUSPARSE_OPERATION_NON_TRANSPOSE, & ! in   op
                      alpha2,matM,vecA,  & ! in                         alph A, X
                      beta2,vecR,  & ! in/out                           beta Y
                      CUDA_R_64F,CUSPARSE_CSRMV_ALG1,buffer)   ! in
  istat=cudaDeviceSynchronize


  istat=cudaMemcpy(p_d(1:nnn),r_d(1:nnn),nnn)


!  istat=cudaMemcpy(rr0_d(1:nnn),r_d(1:nnn),nnn)

  !!!rr0_d is arbitrary!!!!

  call init_random_seed()
  call random_number(rand)
  rr0_d(1:nnn)=rand(1:nnn)




  Ap_d=0d0
  As_d=0d0
  s_d=0d0
  dumR_d=0d0
  temp=0d0
  istat=cudaDeviceSynchronize
  call partialSum<<<256,512,512*8>>>(r_d,rr0_d,part,nnn)
  call finalSum<<<1,256,256*8>>>(part,temp)
  istat=cudaDeviceSynchronize
  old_r_rr=temp



  status=cusparseCreateDnVec(vecAp,nnn,Ap_d,CUDA_R_64F)   

  status=cusparseCreateDnVec(vecAs,nnn,As_d,CUDA_R_64F)   

  status=cusparseCreateDnVec(vecdumR,nnn,dumR_d,CUDA_R_64F)  

  status=cusparseCreateDnVec(vecP,nnn,p_d,CUDA_R_64F)  

  status=cusparseCreateDnVec(vecS,nnn,s_d,CUDA_R_64F)  



  
  temp=0e0
  istat=cudaDeviceSynchronize
  call partialSum<<<256,512,512*8>>>(r_d,r_d,part,nnn)
  call finalSum<<<1,256,256*8>>>(part,temp)
  istat=cudaDeviceSynchronize
  err1=temp
  err0=err1

  k=0
  !write(*,"(A8,I5,ES14.5)") "        ",k,err1
  write(*,*) k,err1


  k=1
  do while(err1>conv)

    alpha2=1d0
    beta2=0d0
    istat=cudaDeviceSynchronize

    !! Y = alph op(A) * X + beta Y
    status=cusparseSpMV(handle,CUSPARSE_OPERATION_NON_TRANSPOSE, & ! in   op
                        alpha2,matM,vecP,  & ! in                         alph A, X
                        beta2,vecAp,  & ! in/out                           beta Y
                        CUDA_R_64F,CUSPARSE_CSRMV_ALG1,buffer)   ! in

    istat=cudaDeviceSynchronize


    temp=0e0
    istat=cudaDeviceSynchronize
    call partialSum<<<256,512,512*8>>>(Ap_d,rr0_d,part,nnn)
    call finalSum<<<1,256,256*8>>>(part,temp)
    istat=cudaDeviceSynchronize

    dum2=temp
    alfa_d=old_r_rr/dum2


    istat=cudaDeviceSynchronize

    call s_updateBI<<<bl_bi,Threads>>>(nnn,s_d,r_d,Ap_d,alfa_d)


    alpha2=1d0
    beta2=0d0
    istat=cudaDeviceSynchronize

    !! Y = alph op(A) * X + beta Y
    status=cusparseSpMV(handle,CUSPARSE_OPERATION_NON_TRANSPOSE, & ! in   op
                        alpha2,matM,vecS,  & ! in                         alph A, X
                        beta2,vecAs,  & ! in/out                           beta Y
                        CUDA_R_64F,CUSPARSE_CSRMV_ALG1,buffer)   ! in

    temp=0e0
    istat=cudaDeviceSynchronize
    call partialSum<<<256,512,512*8>>>(As_d,s_d,part,nnn)
    call finalSum<<<1,256,256*8>>>(part,temp)
    istat=cudaDeviceSynchronize
    dum1=temp
    temp=0e0
    istat=cudaDeviceSynchronize
    call partialSum<<<256,512,512*8>>>(As_d,As_d,part,nnn)
    call finalSum<<<1,256,256*8>>>(part,temp)
    istat=cudaDeviceSynchronize
    dum2=temp
    omega_d=dum1/dum2


    istat=cudaDeviceSynchronize
    call a_r_updateBI<<<bl_bi,Threads>>>(nnn,a_d,r_d,p_d,s_d,As_d,alfa_d,omega_d)

    istat=cudaDeviceSynchronize
    temp=0e0
    istat=cudaDeviceSynchronize
    call partialSum<<<256,512,512*8>>>(r_d,rr0_d,part,nnn)
    call finalSum<<<1,256,256*8>>>(part,temp)
    istat=cudaDeviceSynchronize
    dum1=temp
    istat=cudaDeviceSynchronize
    beta_d=(dum1/old_r_rr)*(alfa_d/omega_d)
    old_r_rr=dum1

    call p_updateBI<<<bl_bi,Threads>>>(nnn,p_d,r_d,Ap_d,beta_d,omega_d)




!!!!verify
    istat=cudaDeviceSynchronize
    alpha2=-1d0
    beta2=1d0
    istat=cudaMemcpy(dumR_d(1:nnn),b_d(1:nnn),nnn)
    istat=cudaDeviceSynchronize

    !! Y = alph op(A) * X + beta Y
    status=cusparseSpMV(handle,CUSPARSE_OPERATION_NON_TRANSPOSE, & ! in   op
                        alpha2,matM,vecA,  & ! in                         alph A, X
                        beta2,vecdumR,  & ! in/out                           beta Y
                        CUDA_R_64F,CUSPARSE_CSRMV_ALG1,buffer)   ! in

    istat=cudaDeviceSynchronize
    temp=0e0
    istat=cudaDeviceSynchronize
    call partialSum<<<256,512,512*8>>>(dumR_d,dumR_d,part,nnn)
    call finalSum<<<1,256,256*8>>>(part,temp)
    istat=cudaDeviceSynchronize

    err1=temp
!!!!verify


!!    istat=cudaDeviceSynchronize
!!    temp=0e0
!!    istat=cudaDeviceSynchronize
!!    call partialSum<<<256,512,512*8>>>(r_d,r_d,part,nnn)
!!    call finalSum<<<1,256,256*8>>>(part,temp)
!!    istat=cudaDeviceSynchronize
!!    err1=temp



    if(mod(k,mo)==0) then
       
       if(abs(err1-err0)/err0<1d-7) then
         write(*,*) k,err1," Not converging! =("
         steps=10001
         goto 123
       else
         write(*,*) k,err1
       endif
       err0=err1
    endif

    if(err1>0e0) then
    else
      write(*,"(A8,I5,ES14.5,A22)") "        ",k,err1," Solution diverged =(!"

      write(*,*) k,err1
      istat=cudaDeviceSynchronize
      a_d=0d0
      istat=cudaDeviceSynchronize
      steps=10001
      goto 123
      !call output
      !stop
    endif
    if(k>5000) then
      write(*,"(A8,I5,ES14.5,A19)") "        ",k,err1," Not converging! =("
      steps=10001
      if(err1<=1d-2) goto 123
      !call output
      stop
    endif

    k=k+1


  enddo


  write(*,*) k,err1," Converged!"
  steps=k

123 continue


  if(allocated(buffer)) deallocate(buffer)

  status=cusparseDestroySpMat(matM)
  status=cusparseDestroyDnVec(vecA)
  status=cusparseDestroyDnVec(vecR)

  status=cusparseDestroyDnVec(vecAp)   
  status=cusparseDestroyDnVec(vecAs)   
  status=cusparseDestroyDnVec(vecdumR)  
  status=cusparseDestroyDnVec(vecP)  
  status=cusparseDestroyDnVec(vecS)  

  status=cusparseDestroy(handle)  

  deallocate(Ap_d,rr0_d)
  deallocate(s_d,dumR_d,As_d)
  deallocate(rand)

end subroutine biccg_stab











subroutine CGNR(conv,mo,nnn,non_zero,bl_bi,steps)
use cudafor
use cusparse

use properties
use parameters
use gpu_properties

use gpu_iccg

Implicit none

  Integer :: mo,nnn,non_zero,bl_bi
  Real(8) :: conv
  Integer :: k,steps
  Real(8) :: err1,err0
  Real(8),device :: temp
  Integer :: istat
  
  Integer :: i,j
 
  Real(8),device,allocatable :: w_d(:)
  Real(8),device,allocatable :: z_d(:)
  Real(8),device,allocatable :: dumR_d(:)
  Real(8),device,allocatable :: part(:) !!!!!! FOR THE partialSum,finalSum subroutines

  Real(8) :: alfa_d,beta_d,omega_d,dum1,dum2

  type(cusparseHandle) :: handle
  Integer :: status
  Real(8) :: alpha2,beta2

  type(cusparseSpMatDescr) :: matM
  type(cusparseDnVecDescr) :: vecA,vecR,vecW,vecZ,vecP
  type(cusparseDnVecDescr) :: vecdumR
  Integer(8) :: bufferSize

  Integer(1),pointer,device :: buffer(:)


  allocate(w_d(nnn),z_d(nnn),dumR_d(nnn),part(512))


  ! initalize CUSPARSE and matrix descriptor
  status=cusparseCreate(handle)
  if(status/=CUSPARSE_STATUS_SUCCESS) write(*,*) 'cusparseCreate error: ', status


  istat=cudaDeviceSynchronize
  istat=cudaMemcpy(r_d(1:nnn),b_d(1:nnn),nnn)
  w_d=0d0
  z_d=0d0
  istat=cudaDeviceSynchronize


  status=cusparseCreateCsr(matM,nnn,nnn,non_zero, &  !mat(cusparse), rows, cols, non zero
                           mRow_d,mCol_d,M_d,  &
                           CUSPARSE_INDEX_32I,CUSPARSE_INDEX_32I, &
                           CUSPARSE_INDEX_BASE_ONE,CUDA_R_64F)

  status=cusparseCreateDnVec(vecA,nnn,a_d,CUDA_R_64F)  ! in

  status=cusparseCreateDnVec(vecR,nnn,r_d,CUDA_R_64F)   ! in

  status=cusparseCreateDnVec(vecW,nnn,w_d,CUDA_R_64F)   

  status=cusparseCreateDnVec(vecZ,nnn,z_d,CUDA_R_64F)   

  status=cusparseCreateDnVec(vecP,nnn,p_d,CUDA_R_64F)  

  status=cusparseCreateDnVec(vecdumR,nnn,dumR_d,CUDA_R_64F)  

  status=cusparseSpMV_buffersize(handle,CUSPARSE_OPERATION_NON_TRANSPOSE, & ! in   op
                                 alpha2,matM,vecA,  & ! in                         alph A, X
                                 beta2,vecR,  & ! in/out                           beta Y
                                 CUDA_R_64F,CUSPARSE_CSRMV_ALG1,bufferSize)   ! in
  istat=cudaDeviceSynchronize
  if(allocated(buffer)) deallocate(buffer)
  if( bufferSize == 0 ) then
    bufferSize=1
  endif
    !  nullify(buffer)
  !else
    allocate(buffer(bufferSize))
  !endif




  write(*,*) "  CGNR method"
  write(*,*) "        Step   Error"


  !! Y = alph op(A) * X + beta Y
  alpha2=-1d0
  beta2=1d0
  status=cusparseSpMV(handle,CUSPARSE_OPERATION_NON_TRANSPOSE, & ! in   op
                      alpha2,matM,vecA,  & ! in                         alph A, X
                      beta2,vecR,  & ! in/out                           beta Y
                      CUDA_R_64F,CUSPARSE_CSRMV_ALG1,buffer)   ! in

  alpha2=1d0
  beta2=0d0
  !! Y = alph op(A) * X + beta Y
  status=cusparseSpMV(handle,CUSPARSE_OPERATION_TRANSPOSE, & ! in   op
                      alpha2,matM,vecR,  & ! in                         alph A, X
                      beta2,vecZ,  & ! in/out                           beta Y
                      CUDA_R_64F,CUSPARSE_CSRMV_ALG1,buffer)   ! in
            
  istat=cudaDeviceSynchronize

  istat=cudaMemcpy(p_d(1:nnn),z_d(1:nnn),nnn)

  


  
  temp=0e0
  istat=cudaDeviceSynchronize
  call partialSum<<<256,512,512*8>>>(r_d,r_d,part,nnn)
  call finalSum<<<1,256,256*8>>>(part,temp)
  istat=cudaDeviceSynchronize
  err1=temp
  err0=err1

  k=0
  !write(*,"(A8,I5,ES14.5)") "        ",k,err1
  write(*,*) k,err1


  k=1
  do while(err1>conv)

    alpha2=1d0
    beta2=0d0
    istat=cudaDeviceSynchronize

    !! Y = alph op(A) * X + beta Y
    status=cusparseSpMV(handle,CUSPARSE_OPERATION_NON_TRANSPOSE, & ! in   op
                        alpha2,matM,vecP,  & ! in                         alph A, X
                        beta2,vecW,  & ! in/out                           beta Y
                        CUDA_R_64F,CUSPARSE_CSRMV_ALG1,buffer)   ! in
    istat=cudaDeviceSynchronize



    temp=0e0
    istat=cudaDeviceSynchronize
    call partialSum<<<256,512,512*8>>>(z_d,z_d,part,nnn)
    call finalSum<<<1,256,256*8>>>(part,temp)
    istat=cudaDeviceSynchronize
    dum1=temp

    temp=0e0
    istat=cudaDeviceSynchronize
    call partialSum<<<256,512,512*8>>>(w_d,w_d,part,nnn)
    call finalSum<<<1,256,256*8>>>(part,temp)
    istat=cudaDeviceSynchronize
    dum2=temp

    istat=cudaDeviceSynchronize
    alfa_d=dum1/dum2
    istat=cudaDeviceSynchronize
    !!! dum1 = zi**2


    call a_r_update<<<bl_bi,Threads>>>(nnn,alfa_d,a_d,r_d,p_d,w_d)



    !! Y = alph op(A) * X + beta Y
    status=cusparseSpMV(handle,CUSPARSE_OPERATION_TRANSPOSE, & ! in   op
                        alpha2,matM,vecR,  & ! in                         alph A, X
                        beta2,vecZ,  & ! in/out                           beta Y
                        CUDA_R_64F,CUSPARSE_CSRMV_ALG1,buffer)   ! in

    temp=0e0
    istat=cudaDeviceSynchronize
    call partialSum<<<256,512,512*8>>>(z_d,z_d,part,nnn)
    call finalSum<<<1,256,256*8>>>(part,temp)
    istat=cudaDeviceSynchronize
    dum2=temp

    beta_d=dum2/dum1


    call p_update<<<bl_bi,Threads>>>(nnn,beta_d,z_d,p_d)




!!!!verify
    istat=cudaDeviceSynchronize
    alpha2=-1d0
    beta2=1d0
    istat=cudaMemcpy(dumR_d(1:nnn),b_d(1:nnn),nnn)
    istat=cudaDeviceSynchronize

    !! Y = alph op(A) * X + beta Y
    status=cusparseSpMV(handle,CUSPARSE_OPERATION_NON_TRANSPOSE, & ! in   op
                        alpha2,matM,vecA,  & ! in                         alph A, X
                        beta2,vecdumR,  & ! in/out                           beta Y
                        CUDA_R_64F,CUSPARSE_CSRMV_ALG1,buffer)   ! in

    istat=cudaDeviceSynchronize
    temp=0e0
    istat=cudaDeviceSynchronize
    call partialSum<<<256,512,512*8>>>(dumR_d,dumR_d,part,nnn)
    call finalSum<<<1,256,256*8>>>(part,temp)
    istat=cudaDeviceSynchronize

    err1=temp
!!!!verify


!    istat=cudaDeviceSynchronize
!    temp=0e0
!    istat=cudaDeviceSynchronize
!    call partialSum<<<256,512,512*8>>>(r_d,r_d,part,nnn)
!    call finalSum<<<1,256,256*8>>>(part,temp)
!    istat=cudaDeviceSynchronize
!    err1=temp



    if(mod(k,mo)==0) then
       if(abs(err1-err0)/err0<1d-7) then
         write(*,*) k,err1," Not converging! =("
         steps=50001
         goto 123
       else
         write(*,*) k,err1
       endif
       err0=err1
    endif

    if(err1>0e0) then
    else
      write(*,"(A8,I5,ES14.5,A22)") "        ",k,err1," Solution diverged =(!"
      istat=cudaDeviceSynchronize
      a_d=0d0
      istat=cudaDeviceSynchronize
      steps=50001
      goto 123

      !write(*,*) k,err1
      !call output
      !stop
    endif
    if(k>50000) then
      write(*,"(A8,I5,ES14.5,A19)") "        ",k,err1," Not converging! =("
      !call output
       !stop
      steps=50001
      goto 123
    endif

    k=k+1


  enddo

  steps=k

  write(*,*) k,err1," Converged!"

123 continue

  if(allocated(buffer)) deallocate(buffer)

  status=cusparseDestroySpMat(matM)
  status=cusparseDestroyDnVec(vecA)
  status=cusparseDestroyDnVec(vecR)

  status=cusparseDestroyDnVec(vecW)   
  status=cusparseDestroyDnVec(vecZ)   
  status=cusparseDestroyDnVec(vecP)  
  status=cusparseDestroyDnVec(vecdumR)  

  status=cusparseDestroy(handle)  

  deallocate(w_d,z_d,dumR_d)

end subroutine CGNR

















subroutine CGNE(conv,mo,nnn,non_zero,bl_bi,steps)
use cudafor
use cusparse

use properties
use parameters
use gpu_properties

use gpu_iccg

Implicit none

  Integer :: mo,nnn,non_zero,bl_bi
  Real(8) :: conv
  Integer :: k,steps
  Real(8) :: err1,err0
  Real(8),device :: temp
  Integer :: istat
  
  Integer :: i,j
 

  Real(8),device,allocatable :: w_d(:)
  Real(8),device,allocatable :: z_d(:)
  Real(8),device,allocatable :: dumR_d(:)
  Real(8),device,allocatable :: part(:) !!!!!! FOR THE partialSum,finalSum subroutines

  Real(8) :: alfa_d,beta_d,omega_d,dum1,dum2

  type(cusparseHandle) :: handle
  Integer :: status
  Real(8) :: alpha2,beta2



  type(cusparseSpMatDescr) :: matM
  type(cusparseDnVecDescr) :: vecA,vecR,vecW,vecZ,vecP
  type(cusparseDnVecDescr) :: vecdumR
  Integer(8) :: bufferSize

  Integer(1),pointer,device :: buffer(:)




  allocate(w_d(nnn),z_d(nnn),dumR_d(nnn),part(512))


  ! initalize CUSPARSE and matrix descriptor
  status=cusparseCreate(handle)
  if(status/=CUSPARSE_STATUS_SUCCESS) write(*,*) 'cusparseCreate error: ', status


  istat=cudaDeviceSynchronize
  istat=cudaMemcpy(r_d(1:nnn),b_d(1:nnn),nnn)
  w_d=0d0
  z_d=0d0
  istat=cudaDeviceSynchronize


  status=cusparseCreateCsr(matM,nnn,nnn,non_zero, &  !mat(cusparse), rows, cols, non zero
                           mRow_d,mCol_d,M_d,  &
                           CUSPARSE_INDEX_32I,CUSPARSE_INDEX_32I, &
                           CUSPARSE_INDEX_BASE_ONE,CUDA_R_64F)

  status=cusparseCreateDnVec(vecA,nnn,a_d,CUDA_R_64F)  ! in

  status=cusparseCreateDnVec(vecR,nnn,r_d,CUDA_R_64F)   ! in

  status=cusparseCreateDnVec(vecW,nnn,w_d,CUDA_R_64F)   

  status=cusparseCreateDnVec(vecZ,nnn,z_d,CUDA_R_64F)   

  status=cusparseCreateDnVec(vecP,nnn,p_d,CUDA_R_64F)  

  status=cusparseCreateDnVec(vecdumR,nnn,dumR_d,CUDA_R_64F)  

  status=cusparseSpMV_buffersize(handle,CUSPARSE_OPERATION_NON_TRANSPOSE, & ! in   op
                                 alpha2,matM,vecA,  & ! in                         alph A, X
                                 beta2,vecR,  & ! in/out                           beta Y
                                 CUDA_R_64F,CUSPARSE_CSRMV_ALG1,bufferSize)   ! in
  istat=cudaDeviceSynchronize
  if(allocated(buffer)) deallocate(buffer)
  if( bufferSize == 0 ) then
    bufferSize=1
  endif
    !  nullify(buffer)
  !else
    allocate(buffer(bufferSize))
  !endif




  write(*,*) "  CGNE method"
  write(*,*) "        Step   Error"


  !! Y = alph op(A) * X + beta Y
  alpha2=-1d0
  beta2=1d0
  status=cusparseSpMV(handle,CUSPARSE_OPERATION_NON_TRANSPOSE, & ! in   op
                      alpha2,matM,vecA,  & ! in                         alph A, X
                      beta2,vecR,  & ! in/out                           beta Y
                      CUDA_R_64F,CUSPARSE_CSRMV_ALG1,buffer)   ! in

  alpha2=1d0
  beta2=0d0
  !! Y = alph op(A) * X + beta Y
  status=cusparseSpMV(handle,CUSPARSE_OPERATION_TRANSPOSE, & ! in   op
                      alpha2,matM,vecR,  & ! in                         alph A, X
                      beta2,vecP,  & ! in/out                           beta Y
                      CUDA_R_64F,CUSPARSE_CSRMV_ALG1,buffer)   ! in
            
  istat=cudaDeviceSynchronize

!!  istat=cudaMemcpy(p_d(1:nnn),z_d(1:nnn),nnn)

  


  
  temp=0e0
  istat=cudaDeviceSynchronize
  call partialSum<<<256,512,512*8>>>(r_d,r_d,part,nnn)
  call finalSum<<<1,256,256*8>>>(part,temp)
  istat=cudaDeviceSynchronize
  err1=temp
  err0=err1

  k=0
  !write(*,"(A8,I5,ES14.5)") "        ",k,err1
  write(*,*) k,err1


  k=1
  do while(err1>conv)

    alpha2=1d0
    beta2=0d0
    istat=cudaDeviceSynchronize

    !! Y = alph op(A) * X + beta Y
    status=cusparseSpMV(handle,CUSPARSE_OPERATION_NON_TRANSPOSE, & ! in   op
                        alpha2,matM,vecP,  & ! in                         alph A, X
                        beta2,vecW,  & ! in/out                           beta Y
                        CUDA_R_64F,CUSPARSE_CSRMV_ALG1,buffer)   ! in
    istat=cudaDeviceSynchronize



    temp=0e0
    istat=cudaDeviceSynchronize
    call partialSum<<<256,512,512*8>>>(r_d,r_d,part,nnn)
    call finalSum<<<1,256,256*8>>>(part,temp)
    istat=cudaDeviceSynchronize
    dum1=temp

    temp=0e0
    istat=cudaDeviceSynchronize
    call partialSum<<<256,512,512*8>>>(p_d,p_d,part,nnn)
    call finalSum<<<1,256,256*8>>>(part,temp)
    istat=cudaDeviceSynchronize
    dum2=temp

    istat=cudaDeviceSynchronize
    alfa_d=dum1/dum2
    istat=cudaDeviceSynchronize
    !!! dum1 = ri**2


    call a_r_update<<<bl_bi,Threads>>>(nnn,alfa_d,a_d,r_d,p_d,w_d)



    !! Y = alph op(A) * X + beta Y
    status=cusparseSpMV(handle,CUSPARSE_OPERATION_TRANSPOSE, & ! in   op
                        alpha2,matM,vecR,  & ! in                         alph A, X
                        beta2,vecZ,  & ! in/out                           beta Y
                        CUDA_R_64F,CUSPARSE_CSRMV_ALG1,buffer)   ! in

    temp=0e0
    istat=cudaDeviceSynchronize
    call partialSum<<<256,512,512*8>>>(r_d,r_d,part,nnn)
    call finalSum<<<1,256,256*8>>>(part,temp)
    istat=cudaDeviceSynchronize
    dum2=temp

    beta_d=dum2/dum1


    call p_update<<<bl_bi,Threads>>>(nnn,beta_d,z_d,p_d)




!!!!verify
    istat=cudaDeviceSynchronize
    alpha2=-1d0
    beta2=1d0
    istat=cudaMemcpy(dumR_d(1:nnn),b_d(1:nnn),nnn)
    istat=cudaDeviceSynchronize

    !! Y = alph op(A) * X + beta Y
    status=cusparseSpMV(handle,CUSPARSE_OPERATION_NON_TRANSPOSE, & ! in   op
                        alpha2,matM,vecA,  & ! in                         alph A, X
                        beta2,vecdumR,  & ! in/out                           beta Y
                        CUDA_R_64F,CUSPARSE_CSRMV_ALG1,buffer)   ! in

    istat=cudaDeviceSynchronize
    temp=0e0
    istat=cudaDeviceSynchronize
    call partialSum<<<256,512,512*8>>>(dumR_d,dumR_d,part,nnn)
    call finalSum<<<1,256,256*8>>>(part,temp)
    istat=cudaDeviceSynchronize

    err1=temp
!!!!verify


!    istat=cudaDeviceSynchronize
!    temp=0e0
!    istat=cudaDeviceSynchronize
!    call partialSum<<<256,512,512*8>>>(r_d,r_d,part,nnn)
!    call finalSum<<<1,256,256*8>>>(part,temp)
!    istat=cudaDeviceSynchronize
!    err1=temp



    if(mod(k,mo)==0) then
       if(abs(err1-err0)/err0<1d-7) then
         write(*,*) k,err1," Not converging! =("
         steps=50001
         goto 123
       else
         write(*,*) k,err1
       endif
       err0=err1
    endif

    if(err1>0e0) then
    else
      write(*,"(A8,I5,ES14.5,A22)") "        ",k,err1," Solution diverged =(!"
      istat=cudaDeviceSynchronize
      a_d=0d0
      istat=cudaDeviceSynchronize
      steps=50001
      goto 123

      !write(*,*) k,err1
      !call output
      !stop
    endif
    if(k>50000) then
      write(*,"(A8,I5,ES14.5,A19)") "        ",k,err1," Not converging! =("
      steps=50001
      goto 123
      !call output
      !stop
    endif

    k=k+1


  enddo

  steps=k
  !write(*,"(A8,I5,ES14.5,A11)") "        ",k,err1," Converged!"
  write(*,*) k,err1," Converged!"





123 continue


  if(allocated(buffer)) deallocate(buffer)

  status=cusparseDestroySpMat(matM)
  status=cusparseDestroyDnVec(vecA)
  status=cusparseDestroyDnVec(vecR)

  status=cusparseDestroyDnVec(vecW)   
  status=cusparseDestroyDnVec(vecZ)   
  status=cusparseDestroyDnVec(vecP)  
  status=cusparseDestroyDnVec(vecdumR)  

  status=cusparseDestroy(handle)  


  deallocate(w_d,z_d,dumR_d)

end subroutine CGNE













subroutine solver(conv,mo)
!Libraries
use cudafor  !Cuda fortran
use thrust   !Thrust: sorting, prefixSum

!Modules to define the most important variables
use properties      !CPU variables
use parameters      !CPU parameters (constants)
use gpu_properties  !GPU variables

Implicit none

  Integer :: mo
  Real(8) :: conv
  real(8),allocatable :: rand(:)

  Integer :: j,istat,dum1
  Integer :: iccg_now




  istat=cudaDeviceSynchronize
!  dum1=mRow_d(n+1)
  dum1=mRow_d(nf+1)
  istat=cudaDeviceSynchronize
  dum1=dum1-1


  678 continue



  call iccg(conv,20,nf,dum1,bl_fluid)

!  call iccg(conv,20,n,dum1,bl)



!  call biccg_stab(conv,mo,n,dum1,bl_fluid,iccg_now)
!
!  if(iccg_now>10000) then
!
!    j=0
!    234 continue
!
!    call CGNE(conv,mo,n,dum1,bl,iccg_now)
!
!    if(iccg_now>50000) then
!  
!      345 continue
!
!      call CGNR(conv,mo,n,dum1,bl,iccg_now)
!
!      if(iccg_now>50000) then
!
!        allocate(rand(n))
!
!        call init_random_seed()
!        call random_number(rand)
!        istat=cudaDeviceSynchronize
!        a_d(1:n)=rand(1:n)
!        istat=cudaDeviceSynchronize
!
!        deallocate(rand)
!
!
!        j=j+1
!        if(j==1) goto 234
!        if(j==2) goto 345
!        if(j==3) then
!
!          conv=conv*10d0 
!          if(conv<1.1d-6) then
!            goto 678
!          else
!            print *, conv
!            call output
!            stop
!          endif
!
!        endif
!      endif
!
!
!    endif
!
!  endif



        
endsubroutine




