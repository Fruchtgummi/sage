FILES="mac_lane/gauss_valuation.py mac_lane/trivial_valuation.py mac_lane/padic_valuation.py mac_lane/valuation.py mac_lane/value_group.py mac_lane/valuation_space.py mac_lane/limit_valuation.py mac_lane/function_field_valuation.py"; while true; do tput reset && (/usr/bin/sage -tp --warn-long .5 --optional=sage,standalone --initial $=FILES && /usr/bin/sage -tp --optional=sage,standalone --long $=FILES && /usr/bin/sage -coverage $=FILES) | less; done
