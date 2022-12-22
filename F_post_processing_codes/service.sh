#vi input.txt

for i in 0 1 2 3 4
do
cd active_ssdd_${i}
cp ../input.txt .
cat "input.txt" | gmx_mpi  trjconv -f md_0_6.xtc -s md_0_6.tpr -n index.ndx -pbc mol -ur compact -center -o md_0_6_center.xtc
cat "input.txt" |  gmx_mpi  trjconv -f md_0_3.xtc -s md_0_3.tpr -n index.ndx -pbc mol -ur compact -center -o md_0_3_center.xtc
cat "input.txt" | gmx_mpi  trjconv -f md_0_4.xtc -s md_0_4.tpr -n index.ndx -pbc mol -ur compact -center -o md_0_4_center.xtc
cat "input.txt" | gmx_mpi  trjconv -f md_0_5.xtc -s md_0_5.tpr -n index.ndx -pbc mol -ur compact -center -o md_0_5_center.xtc
rm -r trajectories_${i}
rm trajectories_${i}.tar.gz
mkdir trajectories_${i}
cd trajectories_${i}
cp ../md_0_6_center.xtc .
cp ../md_0_3_center.xtc .
cp ../md_0_4_center.xtc .
cp ../md_0_5_center.xtc .
cd ../
tar -czvf trajectories_${i}.tar.gz  trajectories_${i}
cd ../
cp active_ssdd_${i}/trajectories_${i}.tar.gz traj/
done

for i in 5 6 7 8 9 
do
cd inactive_ssdd_${i}
cp ../input.txt .
cat "input.txt" | gmx_mpi  trjconv -f md_0_6.xtc -s md_0_6.tpr -n index.ndx -pbc mol -ur compact -center -o md_0_6_center.xtc
cat "input.txt" |  gmx_mpi  trjconv -f md_0_3.xtc -s md_0_3.tpr -n index.ndx -pbc mol -ur compact -center -o md_0_3_center.xtc
cat "input.txt" | gmx_mpi  trjconv -f md_0_4.xtc -s md_0_4.tpr -n index.ndx -pbc mol -ur compact -center -o md_0_4_center.xtc
cat "input.txt" | gmx_mpi  trjconv -f md_0_5.xtc -s md_0_5.tpr -n index.ndx -pbc mol -ur compact -center -o md_0_5_center.xtc
rm -r trajectories_${i}
rm trajectories_${i}.tar.gz
mkdir trajectories_${i}
cd trajectories_${i}
cp ../md_0_6_center.xtc .
cp ../md_0_3_center.xtc .
cp ../md_0_4_center.xtc .
cp ../md_0_5_center.xtc .
cd ../
tar -czvf trajectories_${i}.tar.gz  trajectories_${i}
cd ../
cp inactive_ssdd_${i}/trajectories_${i}.tar.gz traj/
done
