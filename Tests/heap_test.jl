using Test

include("../src/Heap.jl")

times = collect(10:-1:1)
handles = collect(11:20)
nodes = [MutableBinaryHeapNode(t, h) for (t, h) in zip(times, handles)]
heap = _make_mutable_binary_heap(nodes)

@testset "Heap" begin
    @test heap.node_map == [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 6, 10, 9, 7, 3, 4, 5, 2, 1]
    @test [node.time for node in heap.nodes] == [1.0, 2.0, 5.0, 4.0, 3.0, 9.0, 6.0, 10.0, 7.0, 8.0]
    @test [node.handle for node in heap.nodes] == [20, 19, 16, 17, 18, 12, 15, 11, 14, 13]
    update!(heap, 15, 0.5)
    @test [node.time for node in heap.nodes] == [0.5, 2.0, 1.0, 4.0, 3.0, 9.0, 5.0, 10.0, 7.0, 8.0]
    @test [node.handle for node in heap.nodes] == [15, 19, 20, 17, 18, 12, 16, 11, 14, 13]
    update!(heap, 19, Inf)
    @test [node.time for node in heap.nodes] == [0.5, 3.0, 1.0, 4.0, 8.0, 9.0, 5.0, 10.0, 7.0, Inf]
    @test [node.handle for node in heap.nodes] == [15, 18, 20, 17, 13, 12, 16, 11, 14, 19]
end 