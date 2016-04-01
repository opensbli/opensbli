#
# The following environment variables should be predefined:
#
# OPS_INSTALL_PATH
# OPS_COMPILER (gnu,intel,etc)
#

#
# set paths for header files and libraries
#
OPS_INSTALL_PATH := $(OPS_INSTALL_PATH)/c
OPS_INC   = -I$(OPS_INSTALL_PATH)/include
OPS_LIB   = -L$(OPS_INSTALL_PATH)/lib
CUDA_LIB = -L$(CUDA_INSTALL_PATH)/lib64
CUDA_INC = -I$(CUDA_INSTALL_PATH)/include

OPENCL_INC = -I$(OPENCL_INSTALL_PATH)/include
OPENCL_LIB = -L$(OPENCL_INSTALL_PATH)/lib64 -lOpenCL


ifeq ($(OPS_COMPILER),gnu)
  CPP   = g++
ifdef DEBUG
  CPPFLAGS  = -O0 -g -DUNIX -Wall
else
  CPPFLAGS  = -O3 -fPIC -DUNIX -Wall
endif
  OMPFLAGS  = -fopenmp
  MPICPP  = $(MPI_INSTALL_PATH)/bin/mpiCC
  MPIFLAGS  = $(CCFLAGS)
else
ifeq ($(OPS_COMPILER),intel)
  CPP   = icpc
ifdef DEBUG
  CCFLAGS = -O0 -g -no-prec-div $(HDF5_LIB) -openmp -fp-model strict -fp-model source -prec-div -prec-sqrt -DMPICH_IGNORE_CXX_SEEK #-DOPS_DEBUG
else
  #CCFLAGS  = -O3 -ipo -no-prec-div -fp-model strict -fp-model source -prec-div -prec-sqrt -vec-report2 -xSSE4.2 -parallel #-DCOMM_PERF #-DDEBUG
  CCFLAGS = -O3 -ipo -no-prec-div -restrict -fno-alias -fp-model strict -fp-model source -prec-div -prec-sqrt -DMPICH_IGNORE_CXX_SEEK#-vec-report
endif
  CPPFLAGS  = $(CCFLAGS)
  OMPFLAGS  = $(HDF5_LIB) -openmp $(HDF5_LIB) -openmp-report2
  MPICPP  = $(MPI_INSTALL_PATH)/bin/mpicxx
  MPIFLAGS  = $(CCFLAGS)
else
ifeq ($(OPS_COMPILER),mic)
  CPP   = icpc
ifdef DEBUG
  CCFLAGS = -O0 -g -mmic -DMPICH_IGNORE_CXX_SEEK#-vec-report
else
  #CCFLAGS  = -O3 -ipo -no-prec-div -restrict -fno-alias -fp-model strict -fp-model source -prec-div -prec-sqrt -DMPICH_IGNORE_CXX_SEEK#-vec-report
  #CCFLAGS  = -O3 -ipo -no-prec-div -restrict -fno-alias -mmic -vec-report -finline -inline-forceinline -DSIMD_VEC=32 -DMPICH_IGNORE_CXX_SEEK#-vec-report
  CCFLAGS = -O3 -ipo -no-prec-div -restrict -fno-alias -mmic -vec-report -DSIMD_VEC=32 -DMPICH_IGNORE_CXX_SEEK -finline#-vec-report
endif
  CPPFLAGS  = $(CCFLAGS)
  OMPFLAGS  = $(HDF5_LIB) -openmp $(HDF5_LIB) -openmp-report2
  MPICPP  = $(MPI_INSTALL_PATH)/bin/mpicxx
  MPIFLAGS  = $(CCFLAGS)
  OPENCL_LIB = -L$(OPENCL_INSTALL_PATH)/lib64 -lintelocl -lOpenCL #for offload opencl execution on the phi
else
ifeq ($(OPS_COMPILER),pgi)
  CPP     = pgCC
  CUDA_ALIGN_FLAG := -D__x86_64 -D__align__\(n\)=__attribute__\(\(aligned\(n\)\)\) -D__location__\(a\)=__annotate__\(a\) -DCUDARTAPI=
ifdef DEBUG
  CCFLAGS   = -O2 -Minline -Kieee -Minform=inform -Minfo=all
else
  CCFLAGS   = -O3 -Kieee -Minline -Minform=severe -Minfo=all
endif
  CPPFLAGS  = $(CCFLAGS)
  OMPFLAGS  = -mp
  MPICPP  = $(MPI_INSTALL_PATH)/bin/mpicxx
  MPICC   = $(MPI_INSTALL_PATH)/bin/mpicc
  MPIFLAGS  = $(CPPFLAGS)
else
ifeq ($(OPS_COMPILER),cray)
  CPP       = CC
    CUDA_ALIGN_FLAG := -D__x86_64 -D__align__\(n\)=__attribute__\(\(aligned\(n\)\)\) -D__location__\(a\)=__annotate__\(a\) -DCUDARTAPI=
ifdef DEBUG
  CCFLAGS   = -O2
else
  CCFLAGS   = -O3
endif
  CPPFLAGS  = $(CCFLAGS)
  OMPFLAGS  =
  MPICPP    = CC
  MPICC     = cc
  MPIFLAGS  = $(CPPFLAGS)
else
print:
	@echo "unrecognised value for OPS_COMPILER"
endif
endif
endif
endif
endif

ifdef HDF5_INSTALL_PATH
  HDF5_INC 	  	:= -DCHECKPOINTING -I$(HDF5_INSTALL_PATH)/include
  HDF5_LIB 	  	:= -L$(HDF5_INSTALL_PATH)/lib -lhdf5_hl -lhdf5 -lz
  CC = $(MPICC)
  CPP = $(MPICPP)
  OPS_INC += $(HDF5_INC)
endif

# flags for nvcc
# set NV_ARCH to select the correct one
ifndef NV_ARCH
  MESSAGE=select an NVIDA device to compile in CUDA, e.g. make NV_ARCH=KEPLER
  NV_ARCH=Kepler
endif
ifeq ($(NV_ARCH),Fermi)
  CODE_GEN_CUDA=-gencode arch=compute_20,code=sm_21
else
ifeq ($(NV_ARCH),Kepler)
  CODE_GEN_CUDA=-gencode arch=compute_35,code=sm_35
endif
endif


ifdef DEBUG
  NVCC_FLAGS = -O0 -g -G
else
  NVCC_FLAGS = -O3
endif


#
# master to make all versions
#

TARGETS = advection_diffusion_2d_dev_seq advection_diffusion_2d_dev_mpi advection_diffusion_2d_mpi advection_diffusion_2d_mpi_openmp advection_diffusion_2d_openmp advection_diffusion_2d_seq advection_diffusion_2d_cuda advection_diffusion_2d_mpi_cuda advection_diffusion_2d_opencl advection_diffusion_2d_mpi_opencl

ifeq ($(OPS_COMPILER),pgi)
TARGETS += advection_diffusion_2d_openacc advection_diffusion_2d_mpi_openacc
else
TARGETS += openacc_error
endif



all: clean $(TARGETS)

advection_diffusion_2d_dev_seq: Makefile advection_diffusion_2d.cpp advection_diffusion_2d_block_0_kernel.h $(OPS_INSTALL_PATH)/lib/libops_seq.a
	$(CPP) $(CPPFLAGS) $(OPS_INC) $(OPS_LIB) advection_diffusion_2d.cpp -lops_seq $(HDF5_LIB) -o advection_diffusion_2d_dev_seq


advection_diffusion_2d_dev_mpi: Makefile advection_diffusion_2d.cpp advection_diffusion_2d_block_0_kernel.h $(OPS_INSTALL_PATH)/lib/libops_mpi.a
	$(MPICPP) $(MPIFLAGS) -DOPS_MPI $(OPS_INC) $(OPS_LIB) advection_diffusion_2d.cpp -lops_mpi $(HDF5_LIB) -o advection_diffusion_2d_dev_mpi

#
# mpi version
#

SEQ_KERNELS = advection_diffusion_2d_seq_kernels.cpp

advection_diffusion_2d_mpi: Makefile advection_diffusion_2d_ops.cpp advection_diffusion_2d_block_0_kernel.h $(OPS_INSTALL_PATH)/lib/libops_mpi.a
	$(MPICPP) $(MPIFLAGS) -DOPS_MPI $(OPS_INC) $(OPS_LIB) advection_diffusion_2d_ops.cpp -I. ./MPI/$(SEQ_KERNELS) -lops_mpi $(HDF5_LIB) -o advection_diffusion_2d_mpi

#
# OpenMP version
#
OMP_KERNELS = advection_diffusion_2d_omp_kernels.cpp

advection_diffusion_2d_mpi_openmp: Makefile advection_diffusion_2d_ops.cpp advection_diffusion_2d_block_0_kernel.h $(OPS_INSTALL_PATH)/lib/libops_mpi.a
	$(MPICPP) $(OMPFLAGS) $(MPIFLAGS) -DOPS_MPI $(OPS_INC) $(OPS_LIB) advection_diffusion_2d_ops.cpp -I. ./MPI_OpenMP/$(OMP_KERNELS) -lops_mpi $(HDF5_LIB) -o advection_diffusion_2d_mpi_openmp


advection_diffusion_2d_openmp: Makefile advection_diffusion_2d_ops.cpp advection_diffusion_2d_block_0_kernel.h $(OPS_INSTALL_PATH)/lib/libops_seq.a
	$(CPP) $(OMPFLAGS) $(MPIFLAGS) $(OPS_INC) $(OPS_LIB) advection_diffusion_2d_ops.cpp -I. ./MPI_OpenMP/$(OMP_KERNELS) -lops_seq $(HDF5_LIB) -o advection_diffusion_2d_openmp

#
# Sequential version
#
advection_diffusion_2d_seq: clean Makefile advection_diffusion_2d_ops.cpp advection_diffusion_2d_block_0_kernel.h $(OPS_INSTALL_PATH)/lib/libops_seq.a
	$(CPP) $(MPIFLAGS) $(OPS_INC) $(OPS_LIB) advection_diffusion_2d_ops.cpp -I. ./MPI/$(SEQ_KERNELS) -lops_seq $(HDF5_LIB) -o advection_diffusion_2d_seq

#
# CUDA version
#

advection_diffusion_2d_cuda: Makefile ./CUDA/advection_diffusion_2d_kernels_cu.o advection_diffusion_2d_ops.cpp advection_diffusion_2d_block_0_kernel.h $(OPS_INSTALL_PATH)/lib/libops_cuda.a
	$(CPP) $(OMPFLAGS) $(CPPFLAGS) $(CUDA_INC) $(OPS_INC) $(OPS_LIB) $(CUDA_LIB) advection_diffusion_2d_ops.cpp ./CUDA/advection_diffusion_2d_kernels_cu.o -lcudart -lops_cuda $(HDF5_LIB) -o advection_diffusion_2d_cuda

advection_diffusion_2d_mpi_cuda: Makefile ./CUDA/advection_diffusion_2d_kernels_mpi_cu.o advection_diffusion_2d_ops.cpp advection_diffusion_2d_block_0_kernel.h $(OPS_INSTALL_PATH)/lib/libops_mpi_cuda.a
	$(MPICPP) $(OMPFLAGS) $(CPPFLAGS) -DOPS_MPI $(CUDA_INC) $(OPS_INC) $(OPS_LIB) $(CUDA_LIB) advection_diffusion_2d_ops.cpp ./CUDA/advection_diffusion_2d_kernels_mpi_cu.o -lcudart -lops_mpi_cuda $(HDF5_LIB) -o advection_diffusion_2d_mpi_cuda

./CUDA/advection_diffusion_2d_kernels_cu.o: ./CUDA/advection_diffusion_2d_kernels.cu Makefile
	nvcc $(VAR) $(INC) $(NVCC_FLAGS) $(CODE_GEN_CUDA) $(OPS_INC) -I. -c -o ./CUDA/advection_diffusion_2d_kernels_cu.o ./CUDA/advection_diffusion_2d_kernels.cu

./CUDA/advection_diffusion_2d_kernels_mpi_cu.o: ./CUDA/advection_diffusion_2d_kernels.cu Makefile
	nvcc $(VAR) $(INC) $(NVCC_FLAGS) -DOPS_MPI $(CODE_GEN_CUDA) $(OPS_INC) -I. -I$(MPI_INSTALL_PATH)/include \
	-DMPICH_IGNORE_CXX_SEEK -c -o ./CUDA/advection_diffusion_2d_kernels_mpi_cu.o ./CUDA/advection_diffusion_2d_kernels.cu

openacc_obj_list = $(shell find OpenACC/ -name "*.cpp" | sed s/\\.cpp/\\.o/g | tr "\n" " ")
openacc_c_obj_list = $(shell find OpenACC/ -name "*_c.c" | sed s/\\.c/\\.o/g | tr "\n" " ")

./OpenACC/advection_diffusion_2d_kernels.o: ./OpenACC/advection_diffusion_2d_kernels.cpp Makefile
	$(MPICPP) -acc -ta=tesla:cc35 -ta=tesla:maxregcount:64 $(MPIFLAGS) $(CUDA_ALIGN_FLAG) $(OPS_INC) $(CUDA_INC) \
	-I. OpenACC/advection_diffusion_2d_kernels.cpp -c -o OpenACC/advection_diffusion_2d_kernels.o

./OpenACC/%.o: ./OpenACC/%.cpp ./OpenACC/%_c.c Makefile
	@echo "Compiling kernel" $*
	$(MPICC) -acc -ta=tesla:cc35 -ta=tesla:maxregcount:64 $(MPIFLAGS) $(CUDA_ALIGN_FLAG) $(OPS_INC) $(CUDA_INC) \
        -I. OpenACC/$*_c.c -c $(HDF5_LIB) -o OpenACC/$*_c.o
	$(MPICPP) -acc -ta=tesla:cc35 -ta=tesla:maxregcount:64 $(MPIFLAGS) $(CUDA_ALIGN_FLAG) $(OPS_INC) $(CUDA_INC) \
        -I. OpenACC/$*.cpp -c -o OpenACC/$*.o

advection_diffusion_2d_openacc: $(openacc_obj_list) advection_diffusion_2d_ops.cpp advection_diffusion_2d_block_0_kernel.h Makefile $(OPS_INSTALL_PATH)/lib/libops_cuda.a
	$(MPICPP) -acc -ta=tesla:cc35 $(MPIFLAGS) $(OPS_INC) $(OPS_LIB) $(CUDA_INC) $(CUDA_LIB) \
    advection_diffusion_2d_ops.cpp -I. $(openacc_obj_list) $(openacc_c_obj_list) -lcudart -lops_cuda $(HDF5_LIB) -o advection_diffusion_2d_openacc

advection_diffusion_2d_mpi_openacc: $(openacc_obj_list) advection_diffusion_2d_ops.cpp advection_diffusion_2d_block_0_kernel.h Makefile $(OPS_INSTALL_PATH)/lib/libops_mpi_cuda.a
	$(MPICPP) -acc -ta=tesla:cc35 $(MPIFLAGS) $(OPS_INC) $(OPS_LIB) $(CUDA_INC) $(CUDA_LIB) -DOPS_MPI \
    advection_diffusion_2d_ops.cpp -I. $(openacc_obj_list) $(openacc_c_obj_list) -lcudart -lops_mpi_cuda $(HDF5_LIB) -o advection_diffusion_2d_mpi_openacc

advection_diffusion_2d_opencl: ./OpenCL/advection_diffusion_2d_opencl_kernels.o advection_diffusion_2d_ops.cpp advection_diffusion_2d_block_0_kernel.h Makefile $(OPS_INSTALL_PATH)/lib/libops_opencl.a
	$(CPP) $(MPIFLAGS) $(OPS_INC) $(OPS_LIB) $(OPENCL_LIB) advection_diffusion_2d_ops.cpp ./OpenCL/advection_diffusion_2d_opencl_kernels.o -I. -lops_opencl  $(HDF5_LIB) -o advection_diffusion_2d_opencl

./OpenCL/advection_diffusion_2d_opencl_kernels.o: ./OpenCL/advection_diffusion_2d_opencl_kernels.cpp Makefile
	$(CPP) $(CPPFLAGS) $(VAR) $(OPENCL_INC) $(INC) $(OPS_INC) -I. \
                -c -o ./OpenCL/advection_diffusion_2d_opencl_kernels.o ./OpenCL/advection_diffusion_2d_opencl_kernels.cpp

advection_diffusion_2d_mpi_opencl: ./OpenCL/advection_diffusion_2d_mpi_opencl_kernels.o advection_diffusion_2d_ops.cpp advection_diffusion_2d_block_0_kernel.h Makefile $(OPS_INSTALL_PATH)/lib/libops_mpi_opencl.a
	$(MPICPP) $(MPIFLAGS) $(OPS_INC) $(OPS_LIB) $(OPENCL_LIB) advection_diffusion_2d_ops.cpp ./OpenCL/advection_diffusion_2d_mpi_opencl_kernels.o -I. -lops_mpi_opencl  $(HDF5_LIB) -o advection_diffusion_2d_mpi_opencl

./OpenCL/advection_diffusion_2d_mpi_opencl_kernels.o: ./OpenCL/advection_diffusion_2d_opencl_kernels.cpp Makefile
	$(MPICPP) $(CPPFLAGS) $(VAR) $(OPENCL_INC) $(INC) $(OPS_INC) -I. \
                -c -o ./OpenCL/advection_diffusion_2d_mpi_opencl_kernels.o ./OpenCL/advection_diffusion_2d_opencl_kernels.cpp

#
#error
#
openacc_error:
	@echo "advection_diffusion_2d_mpi_openacc cannot be built with $(OPS_COMPILER) compilers .. requires an OpenACC compiler (e.g. pgi)"


#
# cleanup
#

clean:
	rm -f advection_diffusion_2d_dev_seq advection_diffusion_2d_dev_mpi advection_diffusion_2d_mpi advection_diffusion_2d_seq advection_diffusion_2d_openmp advection_diffusion_2d_mpi_openmp advection_diffusion_2d_cuda advection_diffusion_2d_mpi_cuda advection_diffusion_2d_openacc advection_diffusion_2d_mpi_openacc ./CUDA/*.o ./OpenACC/*.o *.o advection_diffusion_2d_opencl advection_diffusion_2d_mpi_opencl ./OpenCL/*.o *.o
	rm -f *.h5
