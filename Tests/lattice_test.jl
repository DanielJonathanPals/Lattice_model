using Test

include("../src/lattice.jl")

l = lattice(lattice_size=(2, 2), initialization="empty", upper_boundary = "periodic", lower_boundary = "periodic", left_boundary = "bounding", right_boundary = "empty")

@testset "lattice" begin
    @test l.state == [0 0; 0 0]
    @test l.state_with_boundaries == [0 0 0 0; 1 0 0 0; 1 0 0 0; 0 0 0 0]
    update_lattice!(l, 1,1,2)
    update_lattice!(l, 2,1,1)
    @test l.state == [2 0; 1 0]
    @test l.state_with_boundaries == [0 1 0 0; 1 2 0 0; 1 1 0 0; 0 2 0 0]
    @test loc_energy(l,epsilon=1.0) == [3 0; 1 1]
    @test l(2,2) == 0
    @test l(2,1) == 1
    @test neighbours(l,2,1) == [2,0,2,1]
    @test loc_energy(l,2,1,epsilon=1.0) == 1
    @test loc_energy(l,1,1,epsilon=1.0) == 3
    @test neighbour_coords(l,1,1)[1] == [true, true, true, false]
    @test neighbour_coords(l,1,1)[2] == [2,1,2,0]
    @test neighbour_coords(l,1,1)[3] == [1,2,1,0]
    @test neighbour_coords(l,2,1)[1] == [true, true, true, false]
    @test neighbour_coords(l,2,1)[2] == [1,2,1,0]
    @test neighbour_coords(l,2,1)[3] == [1,2,1,0]
    @test has_neighbour(l,1,1,"above") == (true,(2,1))
    @test has_neighbour(l,1,1,"right") == (true,(1,2))
    @test has_neighbour(l,1,1,"lower") == (true,(2,1))
    @test has_neighbour(l,1,1,"left") == (false,(0,0))
end