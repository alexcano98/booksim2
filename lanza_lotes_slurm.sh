rm -rf flat_fly/*
rm -rf hyperx/*
#sh scripts/launcher_slurm.sh runfiles/flat_fly/flatfly_topology_2vcs_36x36x36 params_script/ 1 #flat_fly

sh scripts/launcher_slurm_seq.sh runfiles/hyperx/hyperx_topology_2vcs_4x4x4 params_script/ 1 #flat_fly
