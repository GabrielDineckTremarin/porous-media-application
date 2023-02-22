#rm *.mod* # deleta arquivos módulos
rm *.out # deleta executável
rm *.o
sh cleaning.sh  # apaga arquivos que não serão utilizados

rm *.opari.*
rm *.F90
rm *.o

export PATH=/opt/scorep/bin/:${PATH}
export SCOREP_MEMORY_RECORDING=true
export OMP_NUM_THREADS=1
export NO_STOP_MESSAGE=yes
export SCOREP_ENABLE_PROFILING=true
export SCOREP_ENABLE_TRACING=true
export SCOREP_TOTAL_MEMORY=150MB
export SCOREP_EXPERIMENT_DIRECTORY=scorep_51_ACC
#export SCOREP_FILTERING_FILE=scorep.filt
#export ACC_PROFLIB=$SCOREP_LIB/libscorep_adapter_openacc_event.la
export ACC_PROFLIB=/opt/nvidia/hpc_sdk/Linux_x86_64/21.2/compilers/lib/libaccprof.so
export SCOREP_OPENACC_ENABLE=yes

######################################################
#  Compiler                                          #
######################################################

/opt/nvidia/hpc_sdk/Linux_x86_64/21.2/compilers/bin/pgf90 -acc -O2 -g -c \
          comum.f90\
 properties_CH4.f90\
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
      flametip.f90 #\
# -o cylinder_solver.out
#scorep --openacc /opt/nvidia/hpc_sdk/Linux_x86_64/21.2/compilers/bin/pgf90 -fast -acc \
#scorep --openacc /opt/scorep/bin/scorep-pgfortran -fast -acc *.mod -o cylinder_solver.out #\
scorep --openacc  /opt/nvidia/hpc_sdk/Linux_x86_64/21.2/compilers/bin/pgf90 -L/opt/scorep/lib -lscorep_adapter_openacc_event -lscorep_adapter_openacc_mgmt -fast -acc *.o -o cylinder_solver.out #\
#         comum.f90 \
#properties_CH4.f90 \
#       initial.f90 \
#nonsymetric_mesh.f90 \
#     comp_mean.f90 \
#      boundary.f90 \
#         main.f90  \
#     equations.f90 \
#   convergence.f90 \
#     transient.f90 \
#        output.f90 \
#         probe.f90 \
#      flametip.f90 \
#-o cylinder_solver.out


#scorep /opt/nvidia/hpc_sdk/Linux_x86_64/21.2/compilers/bin/pgf90 -acc -O2 -g  -o cylinder_solver.out

./cylinder_solver.out

#rm *.opari.*
#rm *.input.*
#rm *.prep*
#rm *mod*
