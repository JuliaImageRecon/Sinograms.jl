# exts.jl

# methods defined in the extension(s)
export sino_plot_rays, sino_geom_plot!
export ct_geom_plot2!, ct_geom_plot3

sino_plot_rays(::Nothing) = throw("Run `using Plots` first")
sino_geom_plot!(::Nothing) = throw("Run `using Plots` first")

ct_geom_plot2!(::Nothing) = throw("Run `using Plots` first")
ct_geom_plot3(::Nothing) = throw("Run `using Plots` first")
