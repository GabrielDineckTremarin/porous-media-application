#TAMANHO DA MALHA 
MESH_SIZE="100x124"        
                    
#NOME DO DIRETÓRIO QUE VAI SER GERADO
FILENAME="Resultados_O2_$MESH_SIZE"

#VARIÁVEIS COM OS NOMES DOS ARQUIVOS TXT GERADOS NAS SAIDAS
RES_32T="$FILENAME/32T_$MESH_SIZE.txt"
RES_16T="$FILENAME/16T_$MESH_SIZE.txt"
RES_8T="$FILENAME/8T_$MESH_SIZE.txt"
RES_4T="$FILENAME/4T_$MESH_SIZE.txt"
RES_2T="$FILENAME/2T_$MESH_SIZE.txt"
RES_1T="$FILENAME/1T_$MESH_SIZE.txt"


if [ ! -d "$FILENAME" ]; then #SE O DIRETÓRIO NÃO EXISTIR
        mkdir "$FILENAME"
else
        rm "$FILENAME"/*.txt # APAGAR TODOS OS ARQUIVOS .txt DESSE DIRETÓRIO
fi

for i in $(seq 1 30)
do

#echo "ITERAÇÃO ${i}" | tee -a "$RES_32T" "$RES_16T" "$RES_8T" "$RES_4T" "$RES_2T" "$RES_1T"

export OMP_NUM_THREADS=32
./cylinder_solver.out >> "$RES_32T"

export OMP_NUM_THREADS=16
./cylinder_solver.out >> "$RES_16T"

export OMP_NUM_THREADS=8
./cylinder_solver.out >> "$RES_8T"

export OMP_NUM_THREADS=4
./cylinder_solver.out >> "$RES_4T"

export OMP_NUM_THREADS=2
./cylinder_solver.out >> "$RES_2T"

export OMP_NUM_THREADS=1
./cylinder_solver.out >> "$RES_1T"

#echo " " | tee -a "$RES_32T" "$RES_16T" "$RES_8T" "$RES_4T" "$RES_2T" "$RES_1T"

done
