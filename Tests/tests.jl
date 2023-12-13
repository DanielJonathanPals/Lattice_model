using Test

@testset "Tests" begin
    include("lattice_test.jl")
    include("actions_test.jl")
    include("heap_test.jl")
end