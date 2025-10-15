module thrust
  interface thrustsort 

    subroutine sort_int(input,N) bind(C,name="sort_int_wrapper") 
    use iso_c_binding 
      integer(c_int),device:: input(*) 
      integer(c_int),value:: N 
    end subroutine 

    subroutine sort_double(input,N) bind(C,name="sort_double_wrapper") 
    use iso_c_binding 
      real(c_double),device:: input(*) 
      integer(c_int),value:: N 
    end subroutine 

    subroutine sort_by_key_int(key,N,input) bind(C,name="sort_by_key_int_wrapper") 
    use iso_c_binding 
      integer(c_int),device:: key(*)
      integer(c_int),device:: input(*) 
      integer(c_int),value:: N 
    end subroutine 

    subroutine sort_by_key_double(key,N,input) bind(C,name="sort_by_key_double_wrapper") 
    use iso_c_binding 
      real(c_double),device:: key(*)
      Integer(c_int),device:: input(*) 
      integer(c_int),value:: N 
    end subroutine 

    subroutine prefixSum(input,output,N) bind(C,name="prefixSum_wrapper")
    use iso_c_binding
      integer(c_int),device :: input(*)
      integer(c_int),device :: output(*)
      integer(c_int),value  :: N
    end subroutine

  end interface 

end module thrust

