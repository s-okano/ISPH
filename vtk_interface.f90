module vtk_module
      implicit none
      interface vtk_initialize
              subroutine initialize() bind(C,name="initialize")
                      use iso_c_binding
              end subroutine initialize
      end interface vtk_initialize

      interface vtk_setPoints
              subroutine setPoints(n, coord) bind(C,name="setPoints")
                      use iso_c_binding
                      integer(c_int),value :: n
                      real(c_float) :: coord(:,:)
              end subroutine setPoints
      end interface vtk_setPoints

      interface vtk_setPointScalarValue
              subroutine setPointFloatScalarValue(n, value, len_name, name) bind(C,name="setPointFloatScalarValue")
                      use iso_c_binding
                      integer(c_int),value :: n
                      real(c_float) :: value(:)
                      integer(c_int),value :: len_name
                      character(*,c_char) :: name
              end subroutine setPointFloatScalarValue

              subroutine setPointIntScalarValue(n, value, len_name, name) bind(C,name="setPointIntScalarValue")
                use, intrinsic :: iso_c_binding
                integer(c_int),value :: n
                integer(c_int) :: value(:)
                integer(c_int),value :: len_name
                character(*,kind=c_char),intent(in) :: name
              end subroutine setPointIntScalarValue
      end interface vtk_setPointScalarValue

      interface vtk_setPointVectorValue
              subroutine setPointFloatVectorValue(n, value, len_name, name) bind(C,name="setPointFloatVectorValue")
                      use iso_c_binding
                      integer(c_int),value :: n
                      real(c_float) :: value(:,:)
                      integer(c_int),value :: len_name
                      character(c_char),intent(in) :: name
              end subroutine setPointFloatVectorValue
      end interface vtk_setPointVectorValue

      interface vtk_output
              subroutine output(len_filename, filename) bind(C,name="output")
                      use iso_c_binding
                      integer(c_int),value :: len_filename
                      character(c_char),intent(in) :: filename
              end subroutine output
      end interface vtk_output

      interface vtk_finalize
              subroutine finalize() bind(C,name="finalize")
                      use iso_c_binding
              end subroutine finalize
      end interface vtk_finalize
end module vtk_module
