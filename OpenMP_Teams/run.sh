rm *.mod # deleta arquivos módulos
rm *.out # deleta executável
sh cleaning.sh  # apaga arquivos que não serão utilizados

if [ ! -d "Resultados" ]; then
	mkdir Resultados
else
	rm Resultados/*.txt
fi

export NO_STOP_MESSAGE=yes

#/opt/nvidia/hpc_sdk/Linux_x86_64/21.2/compilers/bin/pgf90 -O3 -fopenmp \
/opt/nvidia/hpc_sdk/Linux_x86_64/23.1/compilers/bin/pgf90 -O3 -mp=multicore -Minfo=all \
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

#export OMP_NUM_THREADS=32

#time ./cylinder_solver.out

#export OMP_NUM_THREADS=32
#for i in $(seq 1 30);
#do
#    (time ./cylinder_solver.out) 2>> Resultados/exec_teams_200_32.txt
#done

#export OMP_NUM_THREADS=16
#for i in $(seq 1 30);
#do
#        (time ./cylinder_solver.out) 2>> Resultados/exec_200_16.txt
#done

#export OMP_NUM_THREADS=8
#for i in $(seq 1 30);
#do
#       (time ./cylinder_solver.out) 2>> Resultados/exec_200_8.txt
#done

#export OMP_NUM_THREADS=4
#for i in $(seq 1 30);
#do
#        (time ./cylinder_solver.out) 2>> Resultados/exec_200_4.txt
#done

#export OMP_NUM_THREADS=2
#for i in $(seq 1 30);
#do
#       (time ./cylinder_solver.out) 2>> Resultados/exec_200_2.txt
#done

#export OMP_NUM_THREADS=1
#for i in $(seq 1 1);
#do
#        (time ./cylinder_solver.out) 2>> Resultados/exec_51_1.txt
#done


#(time ./cylinder_solver.out) 2> Resultados/exec_sequencial_200.txt
#0_2_thread_51_malha.txt 
#(time ./cylinder_solver.out) 2> Resultados/exec_1_2_thread_51_malha.txt 
#(time ./cylinder_solver.out) 2> Resultados/exec_2_2_thread_51_malha.txt 
#(time ./cylinder_solver.out) 2> Resultados/exec_3_2_thread_51_malha.txt 
#(time ./cylinder_solver.out) 2> Resultados/exec_4_2_thread_51_malha.txt 




# sh pos.sh

