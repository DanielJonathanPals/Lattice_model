module Lattice_model

export lattice
export neighbours
export has_neighbour
export neighbour_coords
export update_lattice!
export loc_energy

export action
export reaction
export diffusion
export handle_to_action
export action_to_handle
export get_trans_rate

export animate

export init_nodes
export init_heap
export find_dependent_handles
export init_dependency_graph


include("actions.jl")
include("Animation.jl")
include("Next_reaction_methods.jl")

end