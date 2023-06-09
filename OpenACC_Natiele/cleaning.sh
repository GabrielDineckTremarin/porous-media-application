######### rotina para limpar o código #########

#apaga parte transiente
rm -rfv transient/data/*.dat 
rm -rfv transient/*.png

#apaga os restarts transientes
rm -rfv data/results/PTZH/*.dat 
rm -rfv data/results/V/*.dat
rm -rfv data/results/U/*.dat

#rm -rfv data/restart/*.dat

#apaga graficos gerados pelo pos_graphics
rm -rfv output/*.png

#apaga dados gerados pelo pos_process.f90
rm -rfv data/*.dat
rm -rfv data/error.dat
rm -rfv data/flametip.dat
rm -rfv data/probe.dat
rm -rfv data/probe_point.dat
rm -rfv data/*.vtk
rm *.log
rm *.gif
