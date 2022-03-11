sh scripts/launcher_local.sh runfiles/hyperx/hyperx_topology_2vcs params_script/ $0
sh parser_get_stats.sh results/hyperx/hyperx_topology_2vcs
python3 scripts/parser_get_plots.py results/hyperx/hyperx_topology_2vcs/parsed/
