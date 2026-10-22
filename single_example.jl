include("utilities.jl")
using DiscreteMarkovChains
using LinearAlgebra
using CairoMakie
using IntervalRootFinding


min_val = 1
max_val = 5
max_den = 200

numbers = generate_rationals(min_val, max_val, max_den)
τ_values, ϕ_values = find_n_branch_maps(numbers[1:1000], 3)

τ_special = []
ϕ_special = []

no_samples = 100
for i in 1:100
    index = sample(collect(1:length(τ_values)))

    τ = τ_values[index]
    ϕ = ϕ_values[index]

    map_parameters = compute_map(τ, ϕ)   
    map_parameters["no_samples"] = 1e5
    map_parameters["no_attractor_samples"] = 1e4

    temp_periods = zeros(no_samples)

    for k in 1:no_samples
        map_parameters["x₀"] = rationalize(rand(1)[1])
    
        # Run the map.
        x, S = run_map(map_parameters)

        map_data = Dict(
                            "x" => x,
                            "S" => S,
                            "c_points" => map_parameters["c_points"],
                            "partitions" => map_parameters["partitions"],
                            "vector_fields" => map_parameters["vector_fields"],
                            "no_branches" => map_parameters["no_branches"]
        )

        dynamics_data = compute_dynamics(map_data)

        temp_periods[k] = dynamics_data["period"]
    end
    if length(keys(counter(temp_periods))) > 1
        println("τ: ", τ)
        println("ϕ: ", ϕ)
    end

end

upper_bounds  = Int.(sum.(compute_upper_bounds(Z₀_intervals, map_partitions, maximum([gcd_1, gcd_2]))))
period_bound = compute_actual_bound(map_parameters, upper_bounds)


###########################################################################################

τ = 85//84
ϕ = 178//165

map_parameters = compute_map(τ, ϕ)   


Z₀_intervals =  compute_Z₀_intervals(map_parameters)
map_partitions = map(x->convert_to_interval(x), eachrow(map_parameters["partitions"]))
Z₀_gaps = compute_Z₀_gaps(Z₀_intervals, map_partitions)

map_parameters["no_samples"] = 1e5
map_parameters["no_attractor_samples"] = 1e4

Z₀_intervals =  compute_Z₀_intervals(map_parameters)
map_partitions = map(x->convert_to_interval(x), eachrow(map_parameters["partitions"]))
Z₀_gaps = compute_Z₀_gaps(Z₀_intervals, map_partitions)

upper_bounds  = Int.(sum.(compute_upper_bounds(Z₀_intervals, map_partitions, maximum([gcd_1, gcd_2]))))
solutions = compute_all_solutions(map_parameters, upper_bounds)


gcd_1 = gcd(map_parameters["vector_fields"])
gcd_2 = gcd(diam.(Z₀_gaps))



map_parameters["x₀"] = rationalize(rand(1)[1])

x, S = run_map(map_parameters)
map_data = Dict(
    "x" => x,
    "S" => S,
    "c_points" => map_parameters["c_points"],
    "partitions" => map_parameters["partitions"],
    "vector_fields" => map_parameters["vector_fields"],
    "no_branches" => map_parameters["no_branches"]
)

dynamics_data = compute_dynamics(map_data)

fig,ax = plot_map(map_parameters)
plot_intervals(fig, ax, Z₀_intervals, :blue, 0.5)
plot_trajectory(x, 200, :orange)

fig


vector_fields = map_parameters["vector_fields"]
parts = map_parameters["partitions"][:, 1]
a = Z₀_gaps[3]
a = compute_interval_image(a , parts, vector_fields)
plot_intervals(fig, ax, a, :green, 0.5)
fig





plot_trajectory(x, 200, :green)
#plot_intervals(fig, ax, Z₀_gaps, :yellow, 0.5)
fig

save("plots/traj_with_z_star_4.png", fig; pt_per_unit=1)








end

function compute_Z₀_gaps(τ_range, ϕ_range)


    no_τ_points = length(τ_range)
    no_ϕ_points = length(ϕ_range)

    no_gaps = zeros(no_τ_points, no_ϕ_points)
    denominators = zeros(no_τ_points, no_ϕ_points)

    @Threads.threads for i in range(1, no_τ_points)
        
        map_parameters = compute_map.(τ_range[i], ϕ_range)
        no_gaps[i, :] = compute_no_Z₀_gaps.(map_parameters)
        denominators[i, :] = map(x-> denominator(gcd(x["vector_fields"])), map_parameters)
        #denominators[i, :] = denominator.(ϕ_range)

    end

    return no_gaps, denominators
end

function compute_no_Z₀_gaps(map_parameters)
    """
    This function computes the upper bound on the period.
    """
    v = map_parameters["vector_fields"]

    if 0//1 ∈ v || map_parameters["no_branches"] == 2 
        Z₀_gaps = 1.0
    else
        Z₀_intervals =  compute_Z₀_intervals(map_parameters)
        map_partitions = map(x->convert_to_interval(x), eachrow(map_parameters["partitions"]))
        Z₀_gaps = compute_Z₀_gaps_count(Z₀_intervals, map_partitions)
    end

    return length(Z₀_gaps)
end

function compute_period_ub_new(map_parameters)
    """
    This function computes the upper bound on the period.
    """
    v = map_parameters["vector_fields"]

    if 0//1 ∈ v 
        period_bound = 1.0
    elseif map_parameters["no_branches"] == 2 
        v = abs.(v ./ abs.(maximum(v)))
        if sum(v) == 2
            winding_no = 1//1
        else
            winding_no = v[findall(x-> x != 1//1, v)][1]
        end
        period_bound = numerator(winding_no) + denominator(winding_no)
    else
       
        Z₀_intervals =  compute_Z₀_intervals(map_parameters)
        map_partitions = map(x->convert_to_interval(x), eachrow(map_parameters["partitions"]))
        Z₀_gaps = compute_Z₀_gaps(Z₀_intervals, map_partitions)


        gcd_1 = gcd(map_parameters["vector_fields"])
        gcd_2 = gcd(diam.(Z₀_gaps))

        upper_bounds  = Int.(sum.(compute_upper_bounds(Z₀_intervals, map_partitions, maximum([gcd_1, gcd_2]))))
        period_bound = compute_actual_bound(map_parameters, upper_bounds)
    end

    return period_bound

end

function simulate_periods(τ_range, ϕ_range, no_samples)


    no_τ_points = length(τ_range)
    no_ϕ_points = length(ϕ_range)

    periods = zeros(no_τ_points, no_ϕ_points)

    @Threads.threads for i in range(1, no_τ_points)

        τ = τ_range[i]


        for j in range(1, no_ϕ_points)
            
            ϕ = ϕ_range[j]

            # Compute the map.
            map_parameters = compute_map(τ, ϕ)   
            map_parameters["no_samples"] = 1e5
            map_parameters["no_attractor_samples"] = 1e4
    
            temp_periods = zeros(no_samples)

            for k in 1:no_samples
                map_parameters["x₀"] = rationalize(rand(1)[1])
            
                # Run the map.
                x, S = run_map(map_parameters)

                map_data = Dict(
                                    "x" => x,
                                    "S" => S,
                                    "c_points" => map_parameters["c_points"],
                                    "partitions" => map_parameters["partitions"],
                                    "vector_fields" => map_parameters["vector_fields"],
                                    "no_branches" => map_parameters["no_branches"]
                )

                dynamics_data = compute_dynamics(map_data)

                temp_periods[k] = dynamics_data["period"]
            end
            periods[i, j] = maximum(temp_periods)
        end
    end

    return periods
end



function compare_bounds(τ_list, ϕ_list, no_samples, no_initial_points)
    
    computed_upper_bounds = zeros(no_samples)
    sim_upper_bounds = zeros(no_samples)
    τ_values = zeros(Rational{Int64}, no_samples)
    ϕ_values = zeros(Rational{Int64}, no_samples)

    for i in 1:no_samples

        index = sample(range(1, length(τ_list)))
        τ = τ_list[index]
        ϕ = ϕ_list[index]

        τ_values[i] = τ
        ϕ_values[i] = ϕ        
        
        map_parameters = compute_map(τ, ϕ)   
        map_parameters["no_samples"] = 1e5
        map_parameters["no_attractor_samples"] = 1e4
  
        periods = zeros(no_initial_points)
        for j in 1:no_initial_points

            map_parameters["x₀"] = rationalize(rand(1)[1])
          
            # Run the map.
            x, S = run_map(map_parameters)
           
            map_data = Dict(
                "x" => x,
                "S" => S,
                "c_points" => map_parameters["c_points"],
                "partitions" => map_parameters["partitions"],
                "vector_fields" => map_parameters["vector_fields"],
                "no_branches" => map_parameters["no_branches"]
            )

            dynamics_data = compute_dynamics(map_data)
            periods[j] = dynamics_data["period"]
        end


        g_cd = gcd(map_parameters["vector_fields"])
        Z₀_intervals =  compute_Z₀_intervals(map_parameters)
        map_partitions = map(x->convert_to_interval(x), eachrow(map_parameters["partitions"]))
        Z₀_gaps, Z₀_gaps_count = compute_Z₀_gaps(Z₀_intervals, map_partitions, g_cd)
        gcd_ratio = gcd(diam.(Z₀_gaps))/g_cd

        if gcd_ratio <= 1//1
            upper_bounds  = Int.(sum.(compute_upper_bounds(Z₀_intervals, map_partitions, g_cd)))
            computed_bound = compute_actual_bound(map_parameters, upper_bounds)
        elseif length(keys(counter(diam.(Z₀_gaps)))) == 2 && gcd_ratio - floor(gcd_ratio) == 0.0
            computed_bound = Int(collect(values(counter(filter(x -> x>0, diam.(Z₀_gaps)))))[1])
        else 
            upper_bounds = Int.(sum.(compute_upper_bounds(Z₀_intervals, map_partitions, g_cd*gcd_ratio)))
            computed_bound = compute_actual_bound(map_parameters, upper_bounds)
        end
        
        computed_upper_bounds[i] = computed_bound
        sim_upper_bounds[i] = maximum(periods) 

    end

    return computed_upper_bounds, sim_upper_bounds,τ_values,ϕ_values
end







################# Step 1 #################


