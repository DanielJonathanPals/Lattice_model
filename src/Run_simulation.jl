using Pkg
Pkg.activate(".")
using DelimitedFiles

include("Next_reaction_methods.jl")
include("unpack_ARGS.jl")


lattice_params, model_params, simulation_params = unpack_ARGS(ARGS)
l_size = lattice_params["lattice_size"]
l_init = lattice_params["lattice_init"]
l_upper = lattice_params["upper_bound"]
l_lower = lattice_params["lower_bound"]
l_left = lattice_params["left_bound"]
l_right = lattice_params["right_bound"]


println("\nInitializing lattice...")

l = lattice(lattice_size=l_size, 
            initialization=l_init,  
            upper_boundary = l_upper, 
            lower_boundary = l_lower, 
            left_boundary = l_left, 
            right_boundary = l_right)

println("Saving initial lattice configuration...")

isdir("Data") || mkdir("Data")
isdir("Data/" * simulation_params["name"]) || mkdir("Data/" * simulation_params["name"])

open("Data/" * simulation_params["name"] * "/init_lattice_config_w_bound.txt", "w") do io
    writedlm(io, l.state_with_boundaries)
end

println("\nInitializing heap...")

heap = init_heap(l, model_params; actions = simulation_params["allowed_actions"], type_of_diffusion = simulation_params["type_of_diffusion"])

println("\nInitializing dependency graph...")

dep_graph = init_dependency_graph(l)

println("\nRunning simulation...")

trans_count = 0
t = 0
times = Vector{Float64}(undef, simulation_params["max_transitions"])
actions = Vector{Int64}(undef, simulation_params["max_transitions"])

for i in 1:simulation_params["max_transitions"]
    global trans_count += 1
    if i != 1 && times[i-1] > simulation_params["t_max"]
        break
    end

    t = heap.nodes[1].time
    times[i] = t
    top_handle = heap.nodes[1].handle
    actions[i] = top_handle
    act = handle_to_action(heap.nodes[1].handle, l)
    act(l)

    dep_handles = filter(!iszero,dep_graph[top_handle,:])
    old_rates = [heap.nodes[heap.node_map[handle]].trans_rate for handle in dep_handles]
    new_rates = [get_trans_rate(handle, l, model_params; type_of_diffusion=simulation_params["type_of_diffusion"]) for handle in dep_handles]

    # update the drawn node
    new_top_rate = get_trans_rate(top_handle, l, model_params; type_of_diffusion=simulation_params["type_of_diffusion"])
    t_new = t + randexp() / new_top_rate
    update!(heap, top_handle, t_new, new_top_rate)

    # update the dependent nodes
    for (j,handle) in enumerate(dep_handles)
        if handle == 0 
            break
        elseif handle != top_handle
            heap_idx = heap.node_map[handle]
            t_new = old_rates[j]/new_rates[j] * (heap.nodes[heap_idx].time - t) + t
            update!(heap, handle, t_new, new_rates[j])
        end
    end
end

println("\nSimulation finished. Saving data...")

open("Data/" * simulation_params["name"] * "/transition_times.txt", "w") do io
    writedlm(io, times[1:trans_count])
end

open("Data/" * simulation_params["name"] * "/action_handles.txt", "w") do io
    writedlm(io, actions[1:trans_count])
end

open("Data/" * simulation_params["name"] * "/final_lattice_config_w_bound.txt", "w") do io
    writedlm(io, l.state_with_boundaries)
end