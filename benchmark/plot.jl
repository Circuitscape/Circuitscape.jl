using Plots
using Measures
using StatPlots

plotly()

timings = [ 537.810477        106.3954835        89.60318678
			4711.366189        1217.901604        543.0623558
			9486.912186        2337.547827        1124.275759]

labels = ["Python" "Julia" "Julia-CHOLMOD"]
# categories = repeat(str, inner = 3)
groups = repeat(["1m", "6m", "12m"], outer = 3)

groupedbar(groups, timings, bar_position = :dodge, bar_width=0.7, 
		   title = "Circuitscape Benchmarks - 16 processes", 
		   top_margin = 5mm, left_margin = 5mm,
		   xlabel = "Sizes", ylabel = "Time (sec)", framestyle = :box, 
		   label = labels)
