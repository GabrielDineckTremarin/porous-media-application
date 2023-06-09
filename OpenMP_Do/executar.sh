for i in $(seq 1 30)
do

time ./cylinder_solver.out >> saida100x124.txt
export OMP_NUM_THREADS=32
time ./cylinder_solver.out >> saida100x124.txt
export OMP_NUM_THREADS=16
time ./cylinder_solver.out >> saida100x124.txt
export OMP_NUM_THREADS=8
time ./cylinder_solver.out >> saida100x124.txt
export OMP_NUM_THREADS=4
time ./cylinder_solver.out >> saida100x124.txt
export OMP_NUM_THREADS=2
time ./cylinder_solver.out >> saida100x124.txt
export OMP_NUM_THREADS=1
time ./cylinder_solver.out >> saida100x124.txt

done
