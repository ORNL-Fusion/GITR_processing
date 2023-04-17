#!/bin/bash
source ../environment_files/env.fusiont5.sh

cmake -DTHRUST_INCLUDE_DIR=/home/tqd/Code/thrust/thrust \
    -DNETCDF_CXX_INCLUDE_DIR=/home/tqd/Code/spack/opt/spack/linux-ubuntu18.04-x86_64/gcc-7.4.0/netcdf-cxx4-4.3.0-rokguzynggbzgzwzkx27v477ldmiz4ux/include \
    -DNETCDF_CXX_LIBRARY=/home/tqd/Code/spack/opt/spack/linux-ubuntu18.04-x86_64/gcc-7.4.0/netcdf-cxx4-4.3.0-rokguzynggbzgzwzkx27v477ldmiz4ux/lib/libnetcdf_c++4.so \
    -DBoost_DIR=/home/tqd/Code/boost/boostBuild \
    -DBoost_INCLUDE_DIR=/home/tqd/Code/boost/boostBuild/include \
    -DLIBCONFIGPP_LIBRARY=$LIBCONFIGLIB \
    -DMPI_C_LIBRARIES=/cm/shared/apps/mpich/ge/gcc/64/3.2.1/lib \
    -DUSE_CUDA=1 \
    -DUSEMPI=0 \
    -DUSE_MPI=0 \
    -DUSE_BOOST=1 \
    -DUSEIONIZATION=0 \
    -DUSERECOMBINATION=0 \
    -DUSEPERPDIFFUSION=0 \
    -DUSECOULOMBCOLLISIONS=1 \
    -DUSEFRICTION=1 \
    -DUSEANGLESCATTERING=1 \
    -DUSEHEATING=1 \
    -DUSETHERMALFORCE=0 \
    -DUSESURFACEMODEL=0 \
    -DUSESHEATHEFIELD=1 \
    -DBIASED_SURFACE=0 \
    -DUSE_SURFACE_POTENTIAL=0 \
    -DUSEPRESHEATHEFIELD=0 \
    -DBFIELD_INTERP=0 \
    -DLC_INTERP=0 \
    -DGENERATE_LC=0 \
    -DEFIELD_INTERP=0 \
    -DPRESHEATH_INTERP=0 \
    -DDENSITY_INTERP=0 \
    -DTEMP_INTERP=0 \
    -DFLOWV_INTERP=0 \
    -DGRADT_INTERP=0 \
    -DODEINT=0 \
    -DFIXEDSEEDS=1 \
    -DPARTICLESEEDS=1 \
    -DGEOM_TRACE=0 \
    -DGEOM_HASH=0 \
    -DGEOM_HASH_SHEATH=0 \
    -DPARTICLE_TRACKS=0 \
    -DPARTICLE_SOURCE_SPACE=0 \
    -DPARTICLE_SOURCE_ENERGY=0 \
    -DPARTICLE_SOURCE_ANGLE=0 \
    -DPARTICLE_SOURCE_FILE=1 \
    -DSPECTROSCOPY=0 \
    -DUSE3DTETGEOM=0 \
    -DUSECYLSYMM=0 \
    -DUSEFIELDALIGNEDVALUES=0 \
    -DFLUX_EA=1 \
    -DFORCE_EVAL=0 \
    -DUSE_SORT=0 \
    -DCHECK_COMPATIBILITY=1 \
    ..
    #-DMANUAL_FLAGS={$USECUDA_ $USEMPI_ $USEBOOST_ $USEIONIZATION_} \
    #-DMANUAL_FLAGS="-DUSEIONIZATION=1" \
#~/cmake/cmake-3.7.0-rc1-Linux-x86_64/bin/cmake -DCMAKE_BUILD_TYPE=Debug -DTHRUST_INCLUDE_DIR=$CUDA_PATH/include -DNETCDF_DIR=$NETCDF -DNETCDF_CXX_ROOT=$NETCDFCXX4 -DLIBCONFIGPP_LIBRARY=$LIBCONFIGDIR/$LIBCONFIGLIB  ..
