rm *.mod # deleta arquivos módulos
rm *.out # deleta executável
sh cleaning.sh  # apaga arquivos que não serão utilizados

if [ ! -d "Resultados" ]; then
	mkdir Resultados
else
	rm Resultados/*.txt
fi

export NO_STOP_MESSAGE=yes

#/opt/nvidia/hpc_sdk/Linux_x86_64/22.11/compilers/bin/pgf90 -O3 -fopenmp -ta=tesla:cc50 -Minfo=all \
/opt/nvidia/hpc_sdk/Linux_x86_64/23.1/compilers/bin/pgf90 -O3 -mp=gpu -gpu=cc70 -Minfo=all \
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

#for i in $(seq 1 7);
#do
#       (time ./cylinder_solver.out) 2>> Resultados/exec_51_target.txt
#done

# sh pos.sh

