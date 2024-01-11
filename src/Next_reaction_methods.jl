include("Heap.jl")


# Initialize nodes. If `actions` only contains "reaction" then the argument `type_of_diffusion` is ignored.
function init_nodes(l::lattice, params::Dict = default_params; actions::Vector{String}=["reaction", "diffusion"], type_of_diffusion::String="exclusion")
    l_react = []
    l_diff = []
    if "reaction" in actions
        l_react = collect(1:length(l.state) * 3)
    end
    if "diffusion" in actions
        l_diff = collect(length(l.state) * 3 +1:length(l.state) * 5)
    end

    handles = vcat(l_react, l_diff)
    nodes = Vector{MutableBinaryHeapNode}(undef, length(handles))
    times = zeros(length(handles))
    for (i,handle) in enumerate(handles)
        rate = get_trans_rate(handle, l, params; type_of_diffusion=type_of_diffusion)
        times[i] = randexp() / rate
        nodes[i] = MutableBinaryHeapNode(times[i], handle, rate)
    end
    return nodes
end


function init_heap(l::lattice, params::Dict = default_params; actions::Array{String}=["reaction", "diffusion"], type_of_diffusion::String="exclusion")
    nodes = init_nodes(l, params; actions=actions, type_of_diffusion=type_of_diffusion)
    heap = _make_mutable_binary_heap(nodes)
    return heap
end


function find_dependent_handles(handle::Int64, l::lattice; actions::Array{String}=["reaction", "diffusion"])
    action = handle_to_action(handle, l)
    handles = []
    x = action.x_coord
    y = action.y_coord

    if action isa reaction

        # first update other reactions
        has_nbr, x_coord, y_coord = neighbour_coords(l, x, y)
        for c in 0:2
            push!(handles, action_to_handle(reaction(x, y, c), l))
            for i in 1:4
                if has_nbr[i]
                    push!(handles, action_to_handle(reaction(x_coord[i],y_coord[i], c), l))
                end
            end
        end

        # update diffusion reactions
        if "diffusion" in actions
            has_nbr, x_coord, y_coord = neighbour_coords(l, x, y)
            has_nbr[1] && push!(handles, action_to_handle(diffusion(x, y, "above", "exclusion"), l))            # the Type of diffusion is irrelevant here
            has_nbr[2] && push!(handles, action_to_handle(diffusion(x, y, "right", "exclusion"), l))
            has_nbr[3] && push!(handles, action_to_handle(diffusion(x_coord[3], y_coord[3], "above", "exclusion"), l))
            has_nbr[4] && push!(handles, action_to_handle(diffusion(x_coord[4], y_coord[4], "right", "exclusion"), l))
        end
    
    elseif action isa diffusion

        # update other diffusion reactions
        has_nbr, x_coord, y_coord = neighbour_coords(l, x, y)
        has_nbr[1] && push!(handles, action_to_handle(diffusion(x, y, "above", "exclusion"), l))            # the Type of diffusion is irrelevant here
        has_nbr[2] && push!(handles, action_to_handle(diffusion(x, y, "right", "exclusion"), l))
        has_nbr[3] && push!(handles, action_to_handle(diffusion(x_coord[3], y_coord[3], "above", "exclusion"), l))
        has_nbr[4] && push!(handles, action_to_handle(diffusion(x_coord[4], y_coord[4], "right", "exclusion"), l))

        action.swap_with == "above" && (x2 = x_coord[1] ; y2 = y_coord[1])
        action.swap_with == "right" && (x2 = x_coord[2] ; y2 = y_coord[2])
        has_nbr2, x_coord2, y_coord2 = neighbour_coords(l, x2, y2)

        has_nbr2[1] && push!(handles, action_to_handle(diffusion(x2, y2, "above", "exclusion"), l))
        has_nbr2[2] && push!(handles, action_to_handle(diffusion(x2, y2, "right", "exclusion"), l))
        (action.swap_with == "above" && has_nbr2[4]) && push!(handles, action_to_handle(diffusion(x_coord2[4], y_coord2[4], "right", "exclusion"), l))
        (action.swap_with == "right" && has_nbr2[3]) && push!(handles, action_to_handle(diffusion(x_coord2[3], y_coord2[3], "above", "exclusion"), l))

        # update reactions
        if "reaction" in actions
            for c in 0:2
                for i in 1:4
                    if has_nbr[i]
                        push!(handles, action_to_handle(reaction(x_coord[i],y_coord[i], c), l))
                    end
                    if has_nbr2[i]
                        push!(handles, action_to_handle(reaction(x_coord2[i],y_coord2[i], c), l))
                    end
                end
            end
        end
    end
    return handles
end


# returns an array where the i-th row contains the handles of the actions that depend on the i-th handle
function init_dependency_graph(l::lattice)
    graph = zeros(Int64, length(l.state) * 5, 31)                   # Each handle can have at most 31 dependent handles
    for handle in 1:length(l.state) * 5
        dependent_handles = find_dependent_handles(handle, l)
        graph[handle, :] = vcat(dependent_handles, zeros(Int64, 31-length(dependent_handles)))
    end
    return graph
end