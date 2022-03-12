sh scripts/launcher_local.sh runfiles/hyperx/hyperx_topology_2vcs params_script/ $1
sh scripts/parser_get_stats.sh results/hyperx/hyperx_topology_2vcs_sim
python3 scripts/parser_get_plots.py results/hyperx/hyperx_topology_2vcs_sim/parsed/ hyperx
