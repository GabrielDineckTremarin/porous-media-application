rm *.mod # deleta arquivos módulos
rm *.out # deleta executável
rm *.o
sh cleaning.sh  # apaga arquivos que não serão utilizados

export PATH=/opt/scorep/bin/:${PATH}

#export OMP_NUM_THREADS=32
export NO_STOP_MESSAGE=yes
export SCOREP_ENABLE_PROFILING=true
export SCOREP_ENABLE_TRACING=true
#export SCOREP_TOTAL_MEMORY=150MB
export SCOREP_EXPERIMENT_DIRECTORY=scorep_200_32T_2
export SCOREP_FILTERING_FILE=scorep.filt
#export SCOREP_METRIC_PERF_SEP=:

scorep --thread=omp /opt/nvidia/hpc_sdk/Linux_x86_64/23.1/compilers/bin/pgf90 -O3 -fopenmp\
         comum.f90 \
properties_CH4.f90 \
       initial.f90 \
nonsymetric_mesh.f90 \
     comp_mean.f90 \
      boundary.f90 \
         main.f90  \
     equations.f90 \
   convergence.f90 \
     transient.f90 \
        output.f90 \
         probe.f90 \
      flametip.f90 \
-o cylinder_solver.out

#./cylinder_solver.out