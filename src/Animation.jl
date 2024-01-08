import Pkg
Pkg.activate(".")
using Plots
using DelimitedFiles


time_step = 10.
t_max = 1000.0

include("actions.jl")

init_lattice = readdlm("Data/test/init_lattice_config_w_bound.txt")
handles = Int64.(readdlm("Data/test/action_handles.txt"))
final_lattice = readdlm("Data/test/final_lattice_config_w_bound.txt")
boundaries = readdlm("Data/test/lattice_boundaries.txt")
times = readdlm("Data/test/transition_times.txt")

bounds = Dict()
for i in 1:4
    bounds[boundaries[i,1]] = boundaries[i,2]
end

l = lattice(init_lattice[2:end-1, 2:end-1], init_lattice, bounds)


function plot_lattice(l::lattice,t::Float64)
    x_dim = size(l.state,1)
    y_dim = size(l.state,2)
    ratio = (x_dim+1)/(y_dim+1)

    plt = plot(seriestype=:scatter,size=(800,800*ratio),xlims=(0.5,y_dim+0.5),ylims=(0.5,x_dim+0.5),legend=false, grid=false, ticks=false, framestyle = :box)
    title!("t = $t")

    ms = 800/(2*y_dim)
    for i in 1:x_dim, j in 1:y_dim
        if l(i,j) == 0
            scatter!([j],[x_dim - i + 1],color=:white,markersize=ms,markerstrokecolor=:white)
        elseif l(i,j) == 1
            scatter!([j],[x_dim - i + 1],color=:red,markersize=ms)
        elseif l(i,j) == 2
            scatter!([j],[x_dim - i + 1],color=:blue,markersize=ms)
        end
    end

    plt
end


curr_step = 1
anim = @animate for (i,t) in enumerate(0:time_step:t_max)
    while times[curr_step] < t
        act = handle_to_action(handles[curr_step], l)
        act(l)
        global curr_step += 1
    end
    plot_lattice(l,t)
    if (i-1) % 10 == 0
        println("Progress: t = $t, t_max = $t_max")
    end
end

gif(anim, "/scratch/d/Daniel.Pals/Masterthesis/Coding/Lattice_model/Animations/test2.gif", fps=15)






