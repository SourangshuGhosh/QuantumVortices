#sh clear
#echo " $(make total_clean): Object, mod file etc. Flushed out."
#echo "Make to compile the parallel code."
#make ./ns3d.exe
#(time mpiexec -n 4 ./ns3d.exe <parameter.txt) > result.out 2> error.out < /dev/null &
#cd ../src
#sh clear.sh
#make ns3d.exe
#cd ../Run
#cp ../src/ns3d.exe .
#cp data/ux/ux43p* restart/ux/
#cp data/uy/uy43p* restart/uy/
#cp data/uz/uz43p* restart/uz/
#cp data/ranf_seeds/ranf43p* restart/ranf_seeds/
mpiexec -n 4 ./GP2D.exe <parameter.txt > result.out 2> error.out < /dev/null &
#echo "Done successfully."
