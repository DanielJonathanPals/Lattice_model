using Test

include("../src/actions.jl")

l = lattice(lattice_size=(2, 2), initialization="empty", upper_boundary = "periodic", lower_boundary = "periodic", left_boundary = "bounding", right_boundary = "empty")
params = Dict("epsilon" => -log(2.), "k_IB" => (df, dμ, u) -> exp(-df-dμ-u), "dμ" => log(2.), "z_I" => 2., "z_B" => 1.0, "D" => 1.0, "γ" => 1.0)
a1 = reaction(1,1,2)
a2 = reaction(2,1,1)
a3 = diffusion(1,1,"above","exclusion")
a4 = diffusion(2,1,"right","exclusion")
a5 = diffusion(2,2,"above","exchange")
a6 = diffusion(2,2,"right","exchange")
actions = [a1,a2,a3,a4,a5,a6]

@testset "actions" begin
    a1(l)
    @test l.state == [2 0; 0 0]
    a2(l)
    @test l.state == [2 0; 1 0]
    a3(l)
    @test l.state == [1 0; 2 0]
    a4(l)
    @test l.state == [1 0; 0 2]
    handles = [action_to_handle(a, l) for a in actions]
    @test handles == [9, 7, 13, 19, 16, 20]
    @test [handle_to_action(h, l) for h in handles[1:4]] == actions[1:4]
    @test [handle_to_action(h, l; type_of_diffusion = "exchange") for h in handles[5:6]] == actions[5:6]
    @test get_trans_rate(a1, l, params) == 1.0
    @test get_trans_rate(handles[1], l, params) == 1.0
    @test get_trans_rate(a2, l, params) == 1.0
    @test get_trans_rate(handles[2], l, params) == 1.0
    @test get_trans_rate(a3, l, params) == 1.0
    @test get_trans_rate(handles[3], l, params) == 1.0
    @test get_trans_rate(a4, l, params) == 1.0
    @test get_trans_rate(handles[4], l, params) == 1.0
    a2(l)
    @test get_trans_rate(a2, l, params) == 0.0
    @test get_trans_rate(handles[2], l, params) == 0.0
    @test get_trans_rate(a3, l, params) == 0.0
    @test get_trans_rate(handles[3], l, params) == 0.0
    @test get_trans_rate(handles[3], l, params, type_of_diffusion = "exchange") == 0.0
    @test get_trans_rate(a4, l, params) == 0.0
    @test get_trans_rate(handles[4], l, params) == 0.0
    @test get_trans_rate(handles[4], l, params, type_of_diffusion = "exchange") == 1.0
    @test get_trans_rate(a5, l, params) == 1.0
    @test get_trans_rate(handles[5], l, params) == 1.0
    @test get_trans_rate(handles[4], l, params, type_of_diffusion = "exchange") == 1.0
    @test get_trans_rate(a6, l, params) == 0.0
    @test get_trans_rate(handles[6], l, params) == 0.0
    @test get_trans_rate(handles[6], l, params, type_of_diffusion = "exchange") == 0.0
end