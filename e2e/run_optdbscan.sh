julia extract_alloc_cycle_split.jl -s 0 -e 10 -m configAll.ext &


julia optimize_fulldbscan_oc-split.jl -s 10 -e 20 -o votlist-mcmcfull &
julia optimize_fulldbscan_oc-split.jl -s 20 -e 30 -o votlist-mcmcfull &
julia optimize_fulldbscan_oc-split.jl -s 30 -e 40 -o votlist-mcmcfull &
julia optimize_fulldbscan_oc-split.jl -s 40 -e 50 -o votlist-mcmcfull &
julia optimize_fulldbscan_oc-split.jl -s 50 -e 60 -o votlist-mcmcfull &
julia optimize_fulldbscan_oc-split.jl -s 60 -e 70 -o votlist-mcmcfull &
julia optimize_fulldbscan_oc-split.jl -s 70 -e 80 -o votlist-mcmcfull &
julia optimize_fulldbscan_oc-split.jl -s 80 -e 90 -o votlist-mcmcfull &
julia optimize_fulldbscan_oc-split.jl -s 90 -e 100 -o votlist-mcmcfull &
