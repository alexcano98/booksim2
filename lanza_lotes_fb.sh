sh scripts/launcher_local.sh runfiles/flat_fly/flatfly_topology_2vcs params_script/ $1
sh scripts/parser_get_stats.sh results/flat_fly/flatfly_topology_2vcs_sim
python3 scripts/parser_get_plots.py results/flat_fly/flatfly_topology_2vcs_sim/parsed/ flatfly
