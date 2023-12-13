include("actions.jl")
include("Heap.jl")


# Initialize nodes. If `actions` only contains "reaction" then the argument `type_of_diffusion` is ignored.
function init_nodes(l::lattice, params::Dict = default_params; actions::Array{String}=["reaction", "diffusion"], type_of_diffusion::String="exclusion")
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
        nodes[i] = MutableBinaryHeapNode(times[i], handle)
    end
    return nodes
end


# Note that this only yields the correct results when applied BEFORE the action is applied to the lattice.
function find_dependent_handles(handle::Int64, l::lattice; actions::Array{String}=["reaction", "diffusion"])
    action = handle_to_action(handle, l)
    handles = []
    x = action.x_coord
    y = action.y_coord

    if action isa reaction

        # first update other reactions
        for c in 0:2
            handles.push(action_to_handle(reaction(x, y, c), l))
        end
        if l(x,y) == 1 || action.new_state == 1             # if the reaction is a binding or unbinding reaction
            for (x_neighbor,y_neighbor) in [(x-1,y), (x,y+1), (x+1,y), (x,y-1)]
                if l(x_neighbor,x_neighbor) == 1
                    handles.push(action_to_handle(reaction(x_neighbor, y_neighbor, 0), l))
                    handles.push(action_to_handle(reaction(x_neighbor, y_neighbor, 2), l))
                elseif l(x_neighbor,x_neighbor) == 2
                    handles.push(action_to_handle(reaction(x_neighbor, y_neighbor, 1), l))
                end
            end
        end

        # update diffusion reactions
        if "diffusion" in actions
            handles.push(action_to_handle(diffusion(x, y, "above", "exclusion"), l))            # the Type of diffusion is irrelevant here
            handles.push(action_to_handle(diffusion(x, y, "right", "exclusion"), l))
            handles.push(action_to_handle(diffusion(x+1, y, "above", "exclusion"), l))
            handles.push(action_to_handle(diffusion(x, y-1, "right", "exclusion"), l))
        end
    
    elseif action isa diffusion

        # update other diffusion reactions
        handles.push(action_to_handle(diffusion(x, y, "above", "exclusion"), l))
        handles.push(action_to_handle(diffusion(x, y, "right", "exclusion"), l))
        handles.push(action_to_handle(diffusion(x+1, y, "above", "exclusion"), l))
        handles.push(action_to_handle(diffusion(x, y-1, "right", "exclusion"), l))

        if action.swap_with == "above"
            handles.push(action_to_handle(diffusion(x-1, y, "above", "exclusion"), l))
            handles.push(action_to_handle(diffusion(x-1, y, "right", "exclusion"), l))
            handles.push(action_to_handle(diffusion(x-1, y-1, "right", "exclusion"), l))
        elseif action.swap_with == "right"
            handles.push(action_to_handle(diffusion(x, y+1, "above", "exclusion"), l))
            handles.push(action_to_handle(diffusion(x, y+1, "right", "exclusion"), l))
            handles.push(action_to_handle(diffusion(x+1, y+1, "above", "exclusion"), l))
        end

        # update reactions
        if "reaction" in actions
            for c in 0:2
                handles.push(action_to_handle(reaction(x, y, c), l))
            end
            
        end
    end
end