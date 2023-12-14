Int64(str::String) = Int64(parse(Float64, str))
Float64(str::String) = parse(Float64, str)


function unpack_ARGS(args)
    lattice_params = Dict("lattice_size" => Int64.((args[1], args[2])),
                            "lattice_init" => args[3],
                            "upper_bound"  => args[4],
                            "lower_bound"  => args[5],
                            "left_bound"   => args[6],
                            "right_bound"  => args[7])
    model_params = Dict("epsilon" => Float64(args[8]), 
                        "k_IB" => eval(Meta.parse(args[9])), 
                        "dμ" => Float64(args[10]), 
                        "z_I" => Float64(args[11]), 
                        "z_B" => Float64(args[12]), 
                        "D" => Float64(args[13]), 
                        "γ" => Float64(args[14]))

    simulation_params = Dict("t_max" => Float64(args[15]), 
                            "max_transitions" => Int64(args[16]),
                            "allowed_actions" => [String(strip(action,[' '])) for action in split(args[17], ",")],
                            "type_of_diffusion" => args[18],
                            "name" => args[19])


    return lattice_params, model_params, simulation_params
end
