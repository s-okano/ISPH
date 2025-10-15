TARGET  = GPU_SPH.exe
FC      = pgfortran
CC      = g++

###c
CCU	= nvcc
C_OPT	= -std=c++14

###C++
CpDIR   = -I/home/daniel/VTK/include/vtk-9.1
CFLAGS  = -static -std=gnu++14 $(CpDIR)

###Fortran
FDIR	= -L/home/daniel/VTK/lib64
FDIR2   = -L/usr/local/cuda-12.2/lib64 -lcusparse
FLIBS   = -lvtkCommonDataModel-9.1 -lvtkCommonCore-9.1 -lvtkIOXML-9.1 -lvtkCommonMath-9.1 -lvtkCommonMisc-9.1 -lvtkCommonSystem-9.1 -lvtkCommonTransforms-9.1 -lvtksys-9.1 -lvtkIOGeometry-9.1 -lvtkIOXMLParser-9.1 -lvtkIOCore-9.1 -lvtkCommonExecutionModel-9.1

FFLAGS  = -cpp $(DOPT) -cuda -gpu=cuda12.2 -cuda=charstring $(FDIR) $(FLIBS) -lstdc++ #-tp=haswell -Bstatic_pgi -mp


OBJS  = \
	csort.o\
	modules.o \
	thrust_module.o\
	sorting.o \
	tracking.o \
	iccg.o \
	corrector.o \
	vstar.o \
	output.o \
	input.o \
	vtk_interface_f90.o \
	vtk_interface_cxx.o \
	main.o

# --------------------------------------

.SUFFIXES: .o .mod .f90 .cxx .cu

.f90.o:
	$(FC) -c $(FFLAGS) $<
.f90.mod:
	$(FC) -c $(FFLAGS) $<

# --------------------------------------

$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) -o $(TARGET) $(OBJS) $(FDIR2)
	@echo make done.

# --------------------------------------

csort.o: csort.cu
	$(CCU) $(C_OPT) -c csort.cu -o csort.o

modules.o gpu_parameters.mod  gpu_properties.mod  parameters.mod  properties.mod: modules.f90
	$(FC) -c $(FFLAGS) modules.f90

thrust_module.o thrust.mod: thrust_module.f90
	$(FC) -c $(FFLAGS) thrust_module.f90

corrector.o gpu_corrector.mod: corrector.f90
	$(FC) -c $(FFLAGS) corrector.f90

vstar.o gpu_vstar.mod: vstar.f90
	$(FC) -c $(FFLAGS) vstar.f90

iccg.o gpu_iccg.mod: iccg.f90 gpu_properties.mod parameters.mod properties.mod
	$(FC) -c $(FFLAGS) iccg.f90 $(FDIR2)

tracking.o gpu_tracking.mod: tracking.f90 gpu_properties.mod parameters.mod properties.mod
	$(FC) -c $(FFLAGS) tracking.f90

sorting.o gpu_sorting.mod: sorting.f90 gpu_properties.mod parameters.mod properties.mod
	$(FC) -c $(FFLAGS) sorting.f90

input.o: input.f90 properties.mod parameters.mod gpu_properties.mod
	$(FC) -c $(FFLAGS) input.f90

output.o: output.f90 properties.mod parameters.mod vtk_module.mod vtk_interface_cxx.o vtk_interface_f90.o
	$(FC) -c $(FFLAGS) output.f90

vtk_interface_cxx.o: vtk_interface.cxx
	$(CC) $(CFLAGS) -c vtk_interface.cxx -o vtk_interface_cxx.o

vtk_interface_f90.o vtk_module.mod: vtk_interface.f90
	$(FC) $(FFLAGS) -c vtk_interface.f90 -o vtk_interface_f90.o

main.o: main.f90 gpu_vstar.mod gpu_corrector.mod gpu_properties.mod properties.mod parameters.mod thrust.mod
	$(FC) -c $(FFLAGS) main.f90 $(FDIR2)

# --------------------------------------

.PHONY: all clean
all: $(TARGET)

.PHONY: clean
clean:
	rm -f $(TARGET) *.o *.mod
