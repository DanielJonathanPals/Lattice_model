include("lattice.jl")


abstract type action end

struct reaction <: action
    x_coord::Int64
    y_coord::Int64
    new_state::Int64
end

struct diffusion <: action
    x_coord::Int64
    y_coord::Int64
    swap_with::String               # "above" or "right"
    type_of_diffusion::String       # "exclusion", "exchange"
end


# Applies the action `a` to the lattice `l`
function (a::action)(l::lattice)
    if a isa reaction
        update_lattice!(l, a.x_coord, a.y_coord, a.new_state)
    elseif a isa diffusion
        if a.swap_with == "above"
            state1 = l(a.x_coord, a.y_coord)
            state2 = l(a.x_coord-1, a.y_coord)
            update_lattice!(l, a.x_coord, a.y_coord, state2)
            update_lattice!(l, a.x_coord-1, a.y_coord, state1)
        elseif a.swap_with == "right"
            state1 = l(a.x_coord, a.y_coord)
            state2 = l(a.x_coord, a.y_coord+1)
            update_lattice!(l, a.x_coord, a.y_coord, state2)
            update_lattice!(l, a.x_coord, a.y_coord+1, state1)
        end
    end
end


# For the datastructure it is easier to encode the actions into handles which essentially are integers.
function action_to_handle(a::action, l::lattice)
    if a isa reaction
        handle = a.new_state * length(l.state) + (a.x_coord-1)*size(l.state,2) + a.y_coord
    elseif a isa diffusion
        if a.swap_with == "above"
            handle = 3 * length(l.state) + (a.x_coord-1)*size(l.state,2) + a.y_coord
        elseif a.swap_with == "right"
            handle = 4 * length(l.state) + (a.x_coord-1)*size(l.state,2) + a.y_coord
        end
    end
    return Int64(handle)
end


# This function decodes the handle `handle` into the corresponding action for the lattice `l`. Sinc the handels cannot distinguish between different types of diffusion, the type of
# diffusion has to be specified in the argument `type_of_diffusion`. If the handle corresponds to a reaction, the argument `type_of_diffusion` is ignored.
function handle_to_action(handle::Int64, l::lattice; type_of_diffusion::String="exclusion")
    if handle <= 3 * length(l.state)
        new_state = (handle-1) ÷ length(l.state)
        x_coord = ((handle-1) % length(l.state)) ÷ size(l.state,2) + 1
        y_coord = ((handle-1) % length(l.state)) % size(l.state,2) + 1
        return reaction(x_coord, y_coord, new_state)
    else
        x_coord = ((handle-1) % length(l.state)) ÷ size(l.state,2) + 1
        y_coord = ((handle-1) % length(l.state)) % size(l.state,2) + 1
        if (handle - 1) ÷ length(l.state) == 3
            return diffusion(x_coord, y_coord, "above", type_of_diffusion)
        elseif (handle - 1) ÷ length(l.state) == 4
            return diffusion(x_coord, y_coord, "right", type_of_diffusion)
        end
    end
end


default_params = Dict("epsilon" => -2.95, "k_IB" => (df, dμ, u) -> exp(-df-dμ-u), "dμ" => 1.0, "z_I" => 2.0, "z_B" => 1.0, "D" => 1.0, "γ" => 1.0)

# returns the transition rate corresponding to the action `a` for the lattice `l`
function get_trans_rate(a::reaction, l::lattice, params::Dict = default_params)
    old_state = l(a.x_coord, a.y_coord)
    u = loc_energy(l,a.x_coord, a.y_coord, epsilon=params["epsilon"])
    if old_state == 0
        if a.new_state == 1
            return params["D"] * params["z_B"]
        elseif a.new_state == 2
            return params["D"] * params["z_I"]
        elseif a.new_state == 0
            return 0.0
        end
    elseif old_state == 1
        if a.new_state == 2
            return params["k_IB"](log(params["z_I"]/params["z_B"]),params["dμ"],u) * exp(log(params["z_I"]/params["z_B"]) + params["dμ"] + u)
        elseif a.new_state == 0
            return params["D"] * exp(u)
        elseif a.new_state == 1
            return 0.0
        end
    elseif old_state == 2
        if a.new_state == 1
            return params["k_IB"](log(params["z_I"]/params["z_B"]),params["dμ"],u)
        elseif a.new_state == 0
            return params["D"]
        elseif a.new_state == 2
            return 0.0
        end
    end
end

function get_trans_rate(a::diffusion, l::lattice, params::Dict = default_params)
    if a.swap_with == "above"
        if (a.x_coord == 1) && (l.boundaries["lower_boundary"] != "periodic")
            return 0.0
        elseif a.type_of_diffusion == "exclusion"
            if (l(a.x_coord-1, a.y_coord) == 0) ⊻ (l(a.x_coord, a.y_coord) == 0)       # Diffusion occurs if one of the two sites is empty and the other is occupied
                return params["γ"]
            else
                return 0.0
            end
        elseif a.type_of_diffusion == "exchange"
            if l(a.x_coord-1, a.y_coord) != l(a.x_coord, a.y_coord)                     # Diffusion occurs if the two sites are occupied by different particles species (empty is also possible)
                return params["γ"]
            else
                return 0.0
            end
        end
    elseif a.swap_with == "right"
        if (a.y_coord == size(l.state,2)) && (l.boundaries["left_boundary"] != "periodic")
            return 0.0
        elseif a.type_of_diffusion == "exclusion"
            if (l(a.x_coord, a.y_coord+1) == 0) ⊻ (l(a.x_coord, a.y_coord) == 0)        # Diffusion occurs if one of the two sites is empty and the other is occupied
                return params["γ"]
            else
                return 0.0
            end
        elseif a.type_of_diffusion == "exchange"
            if l(a.x_coord, a.y_coord+1) != l(a.x_coord, a.y_coord)                     # Diffusion occurs if the two sites are occupied by different particles species (empty is also possible)
                return params["γ"]
            else
                return 0.0
            end
        end
    end
end

function get_trans_rate(handle::Int64, l::lattice, params::Dict = default_params; type_of_diffusion::String="exclusion")
    action = handle_to_action(handle, l; type_of_diffusion = type_of_diffusion)
    return get_trans_rate(action, l, params)
end