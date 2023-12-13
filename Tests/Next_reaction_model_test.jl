using Test

include("../src/Next_reaction_model.jl")

Random.seed!(1234)
l = lattice(lattice_size=(2, 2), initialization="empty", upper_boundary = "periodic", lower_boundary = "periodic", left_boundary = "bounding", right_boundary = "empty")
params = Dict("epsilon" => -log(2.), "k_IB" => (df, dμ, u) -> exp(-df-dμ-u), "dμ" => log(2.), "z_I" => 2., "z_B" => 1.0, "D" => 1.0, "γ" => 1.0)
update_lattice!(l, 1,1,2)
update_lattice!(l, 2,1,1)

@testset "lattice" begin
    init_nodes(l, params)
end