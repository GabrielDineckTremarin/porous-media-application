rm *.mod # deleta arquivos módulos
rm *.out # deleta executável
sh cleaning.sh  # apaga arquivos que não serão utilizados

#scorep --compiler --openacc  /opt/nvidia/hpc_sdk/Linux_x86_64/21.2/compilers/bin/pgf90  -O3 -fast -acc -Minfo=all \
##/opt/scorep/bin/scorep-pgfortran -02 -fast -acc --openacc\

#export PATH=/opt/scorep/bin/:${PATH}

#export SCOREP_OPENACC_ENABLE=yes
#export OMP_NUM_THREADS=1
#export NO_STOP_MESSAGE=yes
#export SCOREP_ENABLE_PROFILING=false
#export SCOREP_ENABLE_TRACING=true
#export SCOREP_TOTAL_MEMORY=150MB
#export SCOREP_EXPERIMENT_DIRECTORY=scorep_51_ACC
#export SCOREP_FILTERING_FILE=scorep.filt
#export ACC_PROFLIB=/opt/scorep/lib/libscorep_adapter_openacc_event.so
#export ACC_PROFLIB=$SCOREP_LIB/libscorep_adapter_openacc_event.so


#scorep /opt/nvidia/hpc_sdk/Linux_x86_64/23.3/compilers/bin/pgf90  -mp -fast -acc -g -O3\
/opt/nvidia/hpc_sdk/Linux_x86_64/23.3/compilers/bin/pgf90 -Minfo=all -mp -fast -acc -g -O3\
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


# sh pos.sh

