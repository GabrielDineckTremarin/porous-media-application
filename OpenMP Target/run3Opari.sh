#rm *.mod* # deleta arquivos módulos
rm *.out # deleta executável
rm *.o
sh cleaning.sh  # apaga arquivos que não serão utilizados

rm *.opari.*
rm *.F90
rm *.o

export PATH=/opt/scorep/bin/:${PATH}
export SCOREP_MEMORY_RECORDING=true
export OMP_NUM_THREADS=16
export NO_STOP_MESSAGE=yes
export SCOREP_ENABLE_PROFILING=false
export SCOREP_ENABLE_TRACING=true
export SCOREP_TOTAL_MEMORY=150MB
export SCOREP_EXPERIMENT_DIRECTORY=scorep_51_TARGET_16T
export SCOREP_FILTERING_FILE=scorep.filt
#export SCOREP_METRIC_PERF_SEP=:
#export SCOREP_METRIC_PAPI=PAPI_FP_OPS #, PAPI_HW_INT 
#export SCOREP_METRIC_PERF=cycles,page-faults,LLC-load-misses
#export SCOREP_METRIC_RUSAGE=ru_stime:ru_majflt
#export SCOREP_METRIC_RUSAGE=ru_ixrss #ru_utime,ru_stime

######################################################
#  Compiler                                          #
######################################################

#opari2 --nm /opt/scorep/bin/opari2-config 
opari2 comum.f90
opari2 properties_CH4.f90
opari2 initial.f90
opari2 nonsymetric_mesh.f90
opari2 comp_mean.f90
opari2 boundary.f90
opari2 main.f90
opari2 equations.f90
opari2 convergence.f90
opari2 transient.f90
opari2 output.f90
opari2 probe.f90
opari2 flametip.f90
opari2 properties_CH4.f90


/opt/nvidia/hpc_sdk/Linux_x86_64/21.2/compilers/bin/pgf90 -fopenmp -O2 -g -c \
          comum.mod.F90\
 properties_CH4.mod.F90\
       initial.mod.F90 \
nonsymetric_mesh.mod.F90 \
     comp_mean.mod.F90 \
      boundary.mod.F90 \
         main.mod.F90  \
     equations.mod.F90 \
   convergence.mod.F90 \
     transient.mod.F90 \
        output.mod.F90 \
         probe.mod.F90 \
      flametip.mod.F90 #\
# -o cylinder_solver.out


/opt/scorep/bin/scorep-pgfortran -fopenmp -O2 -g  *.mod.o -o cylinder_solver.out

./cylinder_solver.out

#rm *.opari.*
#rm *.input.*
#rm *.prep*
#rm *mod*
