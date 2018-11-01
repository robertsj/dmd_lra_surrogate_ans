# Acquire the .p files
read -p "Enter username: " name
name=${name:-$USER}
# 22x22 2-g diffusion data
scp -P 22222 $name@eigendoit.mne.ksu.edu:/opt/data/LRA/diffusion2x2_ref_with_mesh_temps.p

