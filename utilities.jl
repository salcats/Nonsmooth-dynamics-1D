using StatsBase
using IntervalArithmetic
using IntervalUnionArithmetic
using Combinatorics
using DynamicalSystems
using DataStructures
using CairoMakie
using DataFrames
# Fix the bounding interval I.

const I = intervalUnion(interval(0//1, 1//1))
const map_y_tick_label = Dict(
                                2=> ["c₁"],
                                3=> ["c₂"],
                                4=> ["c₃", "c₂"],
                                5=> ["c₂", "c₄"],
                                6=> ["c₃", "c₂", "c₄"]
)


function fᵣ(x, k, ϕ)
    """
    Vector field to apply in Πᵣ.
    """
    return x - k*ϕ
end

function fₗ(x, j)
    """
    Vector field to apply in Πₗ.
    """
    return x + j
end 


function Dᵣ_check(x)
    """
    Check if the trajectory is currently in Dᵣ.
    """
    return 0 < x < 1  
end


function Dₗ_check(x, ϕ, τ)
    """
    Check if the trajectory is currently in Dₗ.
    """
    return -τ*(ϕ+1)+1 < x < -τ*(ϕ+1)+1+ϕ
end


function convert_to_dataset(array)
    """
    Converts array to dataset for use with the Dynamical Systems package. 
    """
    return DelayEmbeddings.Dataset(array)
end


function construct_K(τ, ϕ)
    """
    Determine the fᵣ iterates.
    """
        k = (τ-1)*((1/ϕ) + 1)
    
    if ceil(k) == k
        k₁ = ceil(k) + 1
    else
        k₁ = ceil(k)
    end
    C = τ/(k₁ + 1 - τ)
    
    if ϕ >= C
        K = [k₁]
    else 
        K = [k₁, k₁ + 1]
    end
       
    return K 
end

function construct_J(τ, ϕ)
    """
    Determine the fₗ iterates.
    """
    j = (τ-1)*(ϕ + 1)

    if ceil(j) == j
        j₁ = ceil(j) + 1
    else 
        j₁ = ceil(j)
    end

    j = τ*(ϕ + 1)-1
    
    if ceil(j) == j
        jₘ = ceil(j) + 1
    else 
        jₘ = ceil(j)
    end

    J = []
    for m in range(0, jₘ - j₁)
        push!(J, j₁ + m)
    end

    return J
end

#### TODO - Come back to this.
function find_nearest_index(a,x)
    idx = searchsortedfirst(a,x)
    if (idx==1); return idx; end
    if (idx>length(a)); return length(a); end
    if (a[idx]==x); return idx; end
    if (abs(a[idx]-x) < abs(a[idx-1]-x))
       return idx
    else
       return idx-1
    end
end

#### TODO - Come back to this.
function compute_dynamics(map_data)
    """
    This function determines the asymptotic behaviour of the map and 
    characterises whether the attractor is cyclic etc.
    """
    x = map_data["x"]
    S = map_data["S"]
    c_points = map_data["c_points"]
    parts = map_data["partitions"]
    vector_fields = map_data["vector_fields"]
    no_branches = map_data["no_branches"]

    # Determine the cycle count.
    period = length(keys(counter(x)))

    # Determine the symbolic itinerary.
    branch_symbols = determine_branch_symbols(vector_fields)

    Z_intervals, Zₖ_sequence = compute_Zₖ_intervals(parts, vector_fields, no_branches)

    attractor_symbols =  S[end-period+1:end]

    attractor_itinerary = join(map(x->branch_symbols[x], attractor_symbols), "-") 
    
    attractor_symbol_count = Dict()

    for (index, Sᵢ) in enumerate(branch_symbols)
        attractor_symbol_count[Sᵢ] = counter(attractor_symbols)[index]
    end

    data = Dict(
                        "period" => period,
                        "attractor_itinerary" => attractor_itinerary,
                        "attractor_symbol_count" => attractor_symbol_count,
                        "Zₖ_sequence" => Zₖ_sequence,
                        "Z_intervals" => Z_intervals,
    )
    return data
end



function compute_Z₀_iterates(map_parameters)
    """
    This function computes the iterates of Z₀ that do not fall inside Z₂.
    """

    vector_fields = map_parameters["vector_fields"]
    parts = map_parameters["partitions"][:, 1]
    no_samples = map_parameters["no_samples"]
    no_attractor_samples = map_parameters["no_attractor_samples"]
    no_branches = map_parameters["no_branches"]
    
    Z_intervals, Zₖ_sequence = compute_Zₖ_intervals(map_parameters["partitions"], vector_fields, no_branches)

    Z₀ = Z_intervals["Z₀"][1]
    Z₂ = Z_intervals["Z₂"][1]

    x₀ = Z₀[1]
    y₀ = Z₀[2]

    x = Vector{Rational}()
    xₙ = x₀
    xₙ₊₁ = x₀

    y = Vector{Rational}()
    yₙ = y₀
    yₙ₊₁ = y₀

    push!(x, xₙ₊₁)
    push!(y, yₙ₊₁)

    for i in range(2, 20)

        s_x_i = searchsortedlast(parts, xₙ)
        s_y_i = searchsortedlast(parts, yₙ)
        
        xₙ₊₁ = xₙ + vector_fields[s_x_i]
        yₙ₊₁ = yₙ + vector_fields[s_y_i]
        
        if determine_set_membership(xₙ₊₁, Z₂)  && !determine_set_membership(yₙ₊₁, Z₂)  
            
            xₙ₊₁ = Z₂[2]
            
            push!(x, xₙ₊₁)
            push!(y, yₙ₊₁)

            xₙ = xₙ₊₁
            yₙ = yₙ₊₁

        elseif !determine_set_membership(xₙ₊₁, Z₂)  && determine_set_membership(yₙ₊₁, Z₂)      
            
            yₙ₊₁ = Z₂[1]

            push!(x, xₙ₊₁)
            push!(y, yₙ₊₁)

            xₙ = xₙ₊₁
            yₙ = yₙ₊₁

        elseif !determine_set_membership(xₙ₊₁, Z₂)  && !determine_set_membership(yₙ₊₁, Z₂)  
            
            push!(x, xₙ₊₁)
            push!(y, yₙ₊₁)

            xₙ = xₙ₊₁
            yₙ = yₙ₊₁ 
        
        else 
            println("absorbed")
            
        end
    end

    return hcat(x, y)
end

#### THIS NEEDS FIXING WITH INTERVAL ARITHMETIC INTERSECTIONS ETC 
function compute_Zₖ_intervals(parts, vector_fields, no_branches)

    inverse_branches = Vector{Interval{Rational{Int64}}}()

    for i in 1:no_branches
        push!(inverse_branches, Interval(parts[i, 1] + vector_fields[i], parts[i, 2] + vector_fields[i]))
    end

    sort!(inverse_branches)

    Z_interval_list = Vector{Interval{Rational{Int64}}}()
    Zₖ_sequence = Vector{Int64}()

    Z_intervals = Dict(
        "Z₀" => Vector{Interval{Rational{Int64}}}(),
        "Z₁" => Vector{Interval{Rational{Int64}}}(),
        "Z₂" => Vector{Interval{Rational{Int64}}}(),

    )
    total_union = compute_sequential_union(inverse_branches)

    Z₀_intervals = filter!(x -> diam(x) > 0, (I \ total_union).v)

    for inter in Z₀_intervals
        push!(Z_intervals["Z₀"], inter[1])
        push!(Z_interval_list, inter[1])
        push!(Zₖ_sequence, 0)
    end

    Z₂_intervals = map(x -> compute_sequential_intersection(x), combinations(inverse_branches, 2))
    filter!(x -> !isempty(x.v) , Z₂_intervals)
    filter!(x -> diam(x[1]) > 0 , Z₂_intervals)

    for inter in Z₂_intervals
        push!(Z_intervals["Z₂"], inter[1])
        push!(Z_interval_list, inter[1])
        push!(Zₖ_sequence, 2)
    end

    if !isempty(Z₀_intervals)
        Z₁_intervals = I \ compute_sequential_union(Z₀_intervals)
    end

    if !isempty(Z₂_intervals)
        Z₁_intervals = Z₁_intervals \ compute_sequential_union(Z₂_intervals)
    end

    if isempty(Z₀_intervals) && isempty(Z₂_intervals)
        for branch in inverse_branches
            push!(Z_intervals["Z₁"], branch)
            push!(Z_interval_list, branch)
            push!(Zₖ_sequence, 1)
        end
    else
        filter!(x -> diam(x) > 0, Z₁_intervals.v)
        for inter in Z₁_intervals.v
            push!(Z_intervals["Z₁"], inter[1])
            push!(Z_interval_list, inter[1])
            push!(Zₖ_sequence, 1)
        end
    end
    Zₖ_sequence = Zₖ_sequence[sortperm(Z_interval_list)]
    
    return Z_intervals, convert_Zₖ_sequence_to_string(Zₖ_sequence)
end


function run_map(map_parameters)
    
    x₀ = map_parameters["x₀"]
    vector_fields = map_parameters["vector_fields"]
    part_start_points = map_parameters["partitions"][:, 1]
    no_samples = map_parameters["no_samples"]
    no_attractor_samples = map_parameters["no_attractor_samples"]

    x = Vector{Rational}()
    xₙ = x₀
    xₙ₊₁ = x₀
    
    attractor_starting_point = no_samples - no_attractor_samples

    for i in range(2, no_samples)
        
        Sᵢ = searchsortedlast(part_start_points, xₙ)
        
        xₙ₊₁ = xₙ + vector_fields[Sᵢ]
        
        if i > attractor_starting_point
            push!(x, xₙ₊₁)
        end

        xₙ = xₙ₊₁
    end

    S = map(x-> searchsortedlast(part_start_points, x), x)
   
    return x, S

end


function compute_band_length(x_sorted, start_index, end_index)
    return x_sorted[end_index] - x_sorted[start_index]
end


function determine_μ(f_components, partition, ϕ)
    """ 
    Determine μ for a given partition.
    """
    μ  =  f_components[partition, 1] - f_components[partition, 2]*ϕ

    return μ[1]
end


function determine_partition(x, parts)
    """
    Determine partition for a given x.
    """
    return findall(x-> x==1, parts[:, 1] .<= x .<= parts[:, 2])

end
   
function compute_map_topology(τ, ϕ, k₁, k₂, j₁, jₘ)

    condition_1 = k₁ == k₂
    #k₁*ϕ > (1+τ*(ϕ+1))
    condition_2 = k₁*ϕ ∈ range(j₁, jₘ)
    condition_3 = ϕ <= j₁/k₁
    condition_4 = ϕ > (j₁/k₁)
    condition_5 = ϕ >= jₘ/k₂
    condition_6 = ϕ < jₘ/k₂


    if condition_1 
        if condition_2
            map_top = 1 
        else 
            map_top = 2
        end
    else
        if condition_3 && condition_5
            map_top = 3
        elseif condition_4 && condition_5
            map_top = 4
        elseif condition_3 && condition_6
            map_top = 5
        else
            map_top = 6
        end
    end

    return map_top
end
   
function compute_discontinuity_points(τ, ϕ, k₁, k₂, j₁, jₘ)

    # Construct dictionary.
    c_points = Dict(
                            "c₁" => Rational(0.0),
                            "c₂" => Rational(0.0),
                            "c₃" => Rational(0.0),
                            "c₄" => Rational(0.0)
    )
    c_points["c₁"] = k₁*ϕ - floor(k₁*ϕ)
    c_points["c₂"] = -τ*(ϕ+1) + k₁*ϕ
    c_points["c₃"] = k₁*ϕ - j₁
    c_points["c₄"] = k₂*ϕ - (jₘ - 1)

    return c_points
end

function compute_map(τ, ϕ)
    """ 
    This functions computes the map. 
    """
    # Compute the kᵢ and jᵢ functions.
    k₁ = compute_k₁(τ, ϕ)
    k₂ = compute_k₂(τ, ϕ)
    j₁ = compute_j₁(τ, ϕ)
    jₘ = compute_jₘ(τ, ϕ)


    map_topology = compute_map_topology(τ, ϕ, k₁, k₂, j₁, jₘ)
    c_points_dict = compute_discontinuity_points(τ, ϕ, k₁, k₂, j₁, jₘ)

    if map_topology == 1
        map_parameters = compute_1_map(τ, ϕ, k₁, k₂, j₁, jₘ, c_points_dict)
    elseif map_topology == 2
        map_parameters = compute_2_map(τ, ϕ, k₁, k₂, j₁, jₘ, c_points_dict)
    elseif map_topology == 3
        map_parameters = compute_3_map(τ, ϕ, k₁, k₂, j₁, jₘ, c_points_dict)
    elseif map_topology == 4
        map_parameters = compute_4_map(τ, ϕ, k₁, k₂, j₁, jₘ, c_points_dict)
    elseif map_topology == 5
        map_parameters = compute_5_map(τ, ϕ, k₁, k₂, j₁, jₘ, c_points_dict)
    elseif map_topology == 6
        map_parameters = compute_6_map(τ, ϕ, k₁, k₂, j₁, jₘ, c_points_dict)
    end
        
    return map_parameters
    
end



function compute_1_map(τ, ϕ, k₁, k₂, j₁, jₘ, c_points_dict)

    parts = zeros(Rational, 1, 2)
    parts[1, :] = [0//1, 1//1]

    map_parameters = Dict(
        "c_points" => [],
        "partitions" => parts, 
        "vector_fields" => [0//1],
        "map_topology" => 1,
        "k_seq" => [k₁*ϕ], 
        "j_seq" => [k₁*ϕ],
        "no_branches" => 1
    )

    return map_parameters

end

function compute_2_map(τ, ϕ, k₁, k₂, j₁, jₘ, c_points_dict)

    c₁ = c_points_dict["c₁"]
    parts = zeros(Rational, 2, 2)
    parts[1, :] = [0//1, c₁]
    parts[2, :] = [c₁, 1//1]

    map_parameters = Dict(
        "c_points" => [c₁],
        "partitions" => parts, 
        "vector_fields" => [1//1 - c₁, -c₁],
        "map_topology" => 2,
        "k_seq" => [k₁, k₁], 
        "j_seq" => [floor(k₁*ϕ), floor(k₁*ϕ) - 1//1],
        "no_branches" => 2
    )

    return map_parameters

end

function compute_3_map(τ, ϕ, k₁, k₂, j₁, jₘ, c_points_dict)

    c₂ = c_points_dict["c₂"]
    parts =  zeros(Rational,2, 2)
    parts[1, :] = [0//1, c₂]
    parts[2, :] = [c₂, 1//1]

    map_parameters = Dict(
        "c_points" => [c₂],
        "partitions" => parts, 
        "vector_fields" => [-k₁*ϕ + j₁, -k₂*ϕ + jₘ],
        "map_topology" => 3,
        "k_seq" => [k₁, k₂], 
        "j_seq" => [j₁, jₘ],
        "no_branches" => 2
    )

    return map_parameters
end

function compute_4_map(τ, ϕ, k₁, k₂, j₁, jₘ, c_points_dict)

    c₂ = c_points_dict["c₂"]
    c₃ = c_points_dict["c₃"]

    parts = zeros(Rational,3, 2)
    parts[1, :] = [0//1, c₃]
    parts[2, :] = [c₃, c₂]
    parts[3, :] = [c₂, 1//1]

    map_parameters = Dict(
        "c_points" => [c₂, c₃],
        "partitions" => parts, 
        "vector_fields" => [
                                -k₁*ϕ + j₁ + 1//1,
                                -k₁*ϕ + j₁,
                                -k₂*ϕ + jₘ,
        ],
        "map_topology" => 4,
        "k_seq" => [k₁, k₁, k₂], 
        "j_seq" => [j₁+1//1, j₁, jₘ],
        "no_branches" => 3
    )


    return map_parameters
end 

function compute_5_map(τ, ϕ, k₁, k₂, j₁, jₘ, c_points_dict)
    
    c₂ = c_points_dict["c₂"]
    c₄ = c_points_dict["c₄"]

    parts =  zeros(Rational,3, 2)
    parts[1, :] = [0//1, c₂]
    parts[2, :] = [c₂, c₄]
    parts[3, :] = [c₄, 1//1]


    map_parameters = Dict(
        "c_points" => [c₂, c₄],
        "partitions" => parts, 
        "vector_fields" => [
                            -k₁*ϕ + j₁,
                            -k₂*ϕ + jₘ,
                            -k₂*ϕ + jₘ - 1//1,
        ],
        "map_topology" => 5,
        "k_seq" => [k₁, k₂, k₂], 
        "j_seq" => [j₁, jₘ, jₘ-1//1],
        "no_branches" => 3
    )

    return map_parameters
end

function compute_6_map(τ, ϕ, k₁, k₂, j₁, jₘ, c_points_dict)

    c₂ = c_points_dict["c₂"]
    c₃ = c_points_dict["c₃"]
    c₄ = c_points_dict["c₄"]

    parts =  zeros(Rational,4, 2)
    parts[1, :] = [0, c₃]
    parts[2, :] = [c₃, c₂]
    parts[3, :] = [c₂, c₄]
    parts[4, :] = [c₄, 1]


    map_parameters = Dict(
        "c_points" => [c₂, c₃, c₄],
        "partitions" => parts, 
        "vector_fields" => [
                        -k₁*ϕ + j₁+1//1,
                        -k₁*ϕ + j₁,
                        -k₂*ϕ + jₘ,
                        -k₂*ϕ + jₘ-1//1,
        ],
        "map_topology" => 6,
        "k_seq" => [k₁, k₁, k₂, k₂], 
        "j_seq" => [j₁+1//1, j₁, jₘ, jₘ-1//1],
        "no_branches" => 4
    )

    return map_parameters
end



function σ(x, σₓ)
    """
    Generic discontinuity function.
    """
    return x > σₓ
end

function run_simulation(x₀, ϕ, τ, ϵₜ, no_samples)
    """
    This function runs the 1D discrete with delay system simulation.
    """

    # Initialise the array.
    x = zeros(no_samples)

    # Set the simulation parameters. 
    σₓ = -(τ-1)

    # Set the initial history conditions. 
    x[1:τ] = x₀

    # Initialise fᵣ counter. 
    fᵣ_counter = 0 

    # Run the simulation.
    for i in range(τ+1, no_samples)

        Πᵣ = σ(x[i-τ], σₓ) 

        if Πᵣ==true
            x[i] = fᵣ(x[i-1], 1, ϕ) 
            fᵣ_counter += 1
        else
            x[i] = fₗ(x[i-1], 1) 
        end

    end

    # Compute the proportion of time spent following fᵣ.
    λ = fᵣ_counter/(no_samples-τ)

    # Filter out the points that land in Dᵣ.
    x_Dᵣ = filter(x_i -> Dᵣ_check(x_i), x)

    return x_Dᵣ, λ
end


function generate_map_ticks(map_topology, c_points) 

    c_points = sort(c_points)
    x_ticks=[0.0]
    x_tick_labels = ["0"]
    x_tick_label_alignment=[(:right, :top)]
    
    for (index, c) in enumerate(c_points)
        push!(x_ticks, c)
        push!(x_tick_labels, map_y_tick_label[map_topology][index])
        if index % 2 == 1
            push!(x_tick_label_alignment, (:left, :top))
        else
            push!(x_tick_label_alignment, (:right, :top))
        end
    end

    push!(x_ticks, 1.0)
    push!(x_tick_labels, "1")
    push!(x_tick_label_alignment, (:left, :top))

    return x_ticks, x_tick_labels, x_tick_label_alignment

end


#### TODO - Can I add labels for the discontinuity points here?
function plot_map(map_parameters)
    """
    This function plots the map.
    """
    # Initialise the figure.
    fig = CairoMakie.Figure()
    fontsize_theme = CairoMakie.Theme(fontsize = 30)
    CairoMakie.set_theme!(fontsize_theme)

    # Extract the midpoint of the intervals for plotting.
    parts = map_parameters["partitions"]
    vector_fields = map_parameters["vector_fields"]
    c_points = map_parameters["c_points"]
    map_topology = map_parameters["map_topology"]
    # Plot the diagonal.
    x_points = CairoMakie.LinRange(0, 1, 100)

    #### TODO - Define a little helper function to produce the desired x ticks 
    #### using the map topology and the c points 

    x_ticks, x_tick_labels, x_tick_label_alignment = generate_map_ticks(map_topology, c_points)


    ax = CairoMakie.Axis(
                fig[1, 1],
                xlabel="xₙ",
                ylabel="xₙ₊₁",
                yticks=([0, 1],["0", "1"]),
                xticks=(x_ticks, x_tick_labels),
                xticklabelalign=x_tick_label_alignment,
    )
    CairoMakie.lines!(ax, x_points, x_points, color=:gray, alpha=0.2)

    CairoMakie.lines!(ax, [0, 0], [0, 1], color=:gray, linestyle=:dash, alpha=0.2)
    CairoMakie.lines!(ax, [1, 0], [1, 1], color=:gray, linestyle=:dash, alpha=0.2)
    CairoMakie.lines!(ax, [1, 1], [1, 0], color=:gray, linestyle=:dash, alpha=0.2)
    CairoMakie.lines!(ax, [1, 0], [0, 0], color=:gray, linestyle=:dash, alpha=0.2)

    # Plot the map lines.
    for i in range(1, size(parts)[1])

        x_n = CairoMakie.LinRange(parts[i,1], parts[i,2], 100)
        x_n_1 = x_n .+ vector_fields[i]
        CairoMakie.lines!(x_n, x_n_1, color=:black, linewidth=4)
    end


    # Plot the discontinuity lines.
    for i in range(1, size(parts)[1]-1)

        x_1 = parts[i,2]
        y_1 = x_1 .+ vector_fields[i]

        x_2 = parts[i+1,1]
        y_2 = x_2 .+ vector_fields[i+1]
    

        CairoMakie.lines!([x_1, x_2], [y_1, y_2], color=:red, linestyle=:dash, linewidth=3)

    end

    #### Change y-axis ticks to just 0 and 1 


    return fig, ax
end




function compute_Q_det(τ::Float64, ϕ::Float64)
    
    no_samples = Int(1e5)
    no_attractor_samples = Int(1e4)

    # Set the initial condition.
    x₀ = interval(rand(1)[1])

    # Compute the map.
    c_points, parts, vector_fields, map_topology, k_seq, j_seq = compute_map(τ, ϕ)

    # Get the partition start points.
    part_start_points = parts[:, 1] 
            
    map_parameters = Dict(  
                            "x₀" => x₀,
                            "vector_fields" => vector_fields,
                            "part_start_points" => part_start_points,
                            "no_samples" => no_samples,
                            "no_attractor_samples" => no_attractor_samples
    )

    # Run the map.
    x, S = run_map_quick(map_parameters)

    map_data = Dict(
                    "x" => x,
                    "S" => S,
                    "c_points" => c_points,
                    "parts" => parts,
                    "vector_fields" => vector_fields
    )

    #data = compute_dynamics(map_data)
    seq = compute_kj_sequence(k_seq, j_seq, S)

    Q, state_count = compute_empirical_transition_matrix(seq)

    chain = DiscreteMarkovChain(Q)

    ϵ = 1e-10
    mixing_time_ub = (1/(1-eigvals(Q)[1]))*log(1/(ϵ*minimum(stationary_distribution(chain))))


    return mixing_time_ub
    #return data["cycle_count"]


end

function compute_kj_sequence(k_seq, j_seq, S)

    seq = []

    for f ∈ S

        for k in range(1, k_seq[f])
            push!(seq, 1)
        end

        for j in range(1, j_seq[f])
            push!(seq, 2)
        end
    end

    return Int.(seq)
end



function run_loop(τ_range, ϕ_range, n)

    mixing_time = zeros(n, n)

    for i in range(1, n)
        for j in range(1, n)
            mixing_time[i, j]  = compute_Q_det(τ_range[i], ϕ_range[j])
        end
    end
        return mixing_time
end


"""
function add_bifurcation_lines(fig, ax, τ_min, τ_max, ϕ_min, ϕ_max)

    τ_range = LinRange(Float64(τ_min), Float64(τ_max), 100000)
    ϕ_range = LinRange(Float64(ϕ_min), Float64(ϕ_max), 100000)
    
    τ_points_b, ϕ_points_b = compute_bisections(τ_range, ϕ_range)
    τ_points_z_0, ϕ_points_z_0, τ_points_z_0_d, ϕ_points_z_0_d = compute_zeros(floor(τ_min), ceil(τ_max), floor(ϕ_min), ceil(ϕ_max))
        
    
    τ_points_z = vcat(τ_points_z_0, τ_points_z_0_d)
    ϕ_points_z = vcat(ϕ_points_z_0, ϕ_points_z_0_d)

    τ_indices = τ_min .< τ_points_z .< τ_max
    ϕ_indices = ϕ_min .< ϕ_points_z .< ϕ_max

    if !isempty(τ_indices)
        CairoMakie.scatter!(
                    ax,
                    τ_points_z[τ_indices .+ ϕ_indices .== 2],
                    ϕ_points_z[τ_indices .+ ϕ_indices .== 2],
                    color=:black,
                    markersize=2,
                    transparency=false,
        )
    end
  

    for key in keys(τ_points_b) 
        if !isempty(τ_points_b[key])

            τ_points = τ_points_b[key]
            ϕ_points = ϕ_points_b[key]
                    
            indices = compute_k₁_float.(τ_points, ϕ_points).*ϕ_points .<= 1 .+ τ_points.*(ϕ_points .+ 1)
                    
            τ_indices = τ_min .< τ_points .< τ_max
            ϕ_indices = ϕ_min .< ϕ_points .< ϕ_max

            
            if !isempty(indices)
                CairoMakie.scatter!(
                    ax,
                    τ_points[τ_indices .+ ϕ_indices .+ indices .== 3],
                    ϕ_points[τ_indices .+ ϕ_indices .+ indices .== 3],
                    color=:black,
                    markersize=2,
                    transparency=false,
                )
            end
        end
    end
    return fig, ax
end
"""
function compute_map_set(ϕ_range, τ_range)
    """
    This function computes the set of possible maps by 
    computing all permissible pairs of (k, j). 
    """
    k_possible = []

    # Compute possible values for k.
    for ϕ in ϕ_range
        k_possible = push!(k_possible, ceil.((0 .+τ_range.*ϕ)/ϕ))
        k_possible = push!(k_possible, ceil.((1 .+τ_range.*ϕ)/ϕ))
    end

    k_possible = Float64.(keys(counter(vcat(k_possible...))))

    kj_dict = Dict(k => [] for k in k_possible)

    # For each k, compute possible values of j
    for k ∈ k_possible

        j_possible = []
    
        for ϕ in ϕ_range
            j_possible = push!(j_possible, ceil.(k*ϕ .+ τ_range .- 0))
            j_possible = push!(j_possible, ceil.(k*ϕ .+ τ_range .- 1))
        end
        j_possible = Float64.(keys(counter(vcat(j_possible...))))
        kj_dict[k] = j_possible
    end

    # Determine total number of unique maps.
    no_maps = sum(map(x-> length(kj_dict[x]), k_possible))

    # Initialise data.
    map_set = zeros(no_maps, 2)

    inc = 1
    for k in sort(k_possible)
        for j in kj_dict[k]
            map_set[inc, :] = [k, j]
            inc += 1
        end
    end
    
    return map_set

end



#### TODO - This needs updating
function overlay_inverse_maps(parts, vector_fields)

    # Extract the midpoints of the intervals for plotting.
    parts = mid.(parts)
    vector_fields = mid.(vector_fields)

    # Plot
    for i in range(1, size(parts)[1])
        x_n = LinRange(parts[i, 1] + vector_fields[i], parts[i, 2] + vector_fields[i], 100)
        x_n_1 = x_n .- vector_fields[i]
        lines!(x_n, x_n_1, color=:orange, linewidth=4)
    end
end



function plot_trajectory(x, no_points, colour)
    """
    Plot the trajectory for a given number of points.
    """
    # Extract the midpoints for plotting.
    x = Float64.(x)

    for i in range(1, no_points)
        x_p = x[end-(no_points+3) + i]
        y_p = x[end-(no_points+3) + i+1]
        f_y_p = x[end-(no_points+3) + i+2]

        CairoMakie.lines!([x_p, y_p], [y_p, y_p],  color=colour, alpha=0.01)
        CairoMakie.lines!([y_p, y_p], [y_p, f_y_p], color=colour, alpha=0.01)
    end

end

function calkin_wilf(n)
    # Computes the first n terms in the calkin-wilf rational number tree.
    cw = zeros(Rational, n + 1)
    for i in 2:n + 1
        t = Int(floor(cw[i - 1])) * 2 - cw[i - 1] + 1
        cw[i] = 1 // t
    end
    return cw[2:end]
end


function farey(n::Int)
    rst = SortedSet{Rational}()
    for den in 1:n, num in 1:den-1
        push!(rst, Rational(num, den))
    end
    return collect(rst)
end


function generate_maps(τ_range, ϕ_range, ϵₜ)

    x_max = 1

    τ_dict = Dict(
                    "0" => [],
                    "1a" => [],
                    "1b" => [],
                    "2a" => [],
                    "2b" => [],
                    "3" => []
    )

    ϕ_dict = Dict(
                    "0" => [],
                    "1a" => [],
                    "1b" => [],
                    "2a" => [],
                    "2b" => [],
                    "3" => []
    )



    for τ in τ_range
        for ϕ in ϕ_range
        
            x_discontinuities, label = compute_discontinuities(τ, ϕ)

            push!(τ_dict[label], τ)
            push!(ϕ_dict[label], Float64(ϕ))
        end
    end

    return τ_dict, ϕ_dict
end




function check_bands_cyclicity(x, x_sorted, boundary_indices, cycle_count, no_attractor_samples)
    """
    This function checks if the band structure is cyclic or not.
    """
    # Find boundary points.
    boundary_points = x_sorted[boundary_indices]
        
    # Compute band itinerary.
    band_itinerary = map(x-> searchsortedlast(boundary_points, x), x)
        
    # Run through each band and check if it is cyclic.
    band_cyclicity = 0
    for band_number in range(1, cycle_count)
        band_cyclicity += check_band_cyclicity(band_itinerary, cycle_count, band_number, no_attractor_samples-cycle_count)
    end

    return band_cyclicity==cycle_count
end

function check_band_cyclicity(band_itinerary, cycle_count, band_number, no_samples)
    return length(unique(band_itinerary[band_number+cycle_count:cycle_count:no_samples]))==1
end

function compute_map_type_masks(τ_range, ϕ_range)
    """
    This function computes a binary mask over the (τ, ϕ)-parameter space for each 
    map type.
    """
    no_τ = length(τ_range)
    no_ϕ = length(ϕ_range)

    map_type_masks = Dict(    
                            "0" => zeros(no_ϕ, no_τ),
                            "1a" => zeros(no_ϕ, no_τ),
                            "1b" => zeros(no_ϕ, no_τ),
                            "2a" => zeros(no_ϕ, no_τ),
                            "2b" => zeros(no_ϕ, no_τ),
                            "3" => zeros(no_ϕ, no_τ)
    )

    ϕ_array = zeros(no_ϕ, no_τ)
    τ_array = zeros(no_ϕ, no_τ)

    for ϕ_i in eachindex(ϕ_range)
        for τ_i in eachindex(τ_range)
            
            ϕ = ϕ_range[ϕ_i]
            τ = τ_range[τ_i]
            
            ϕ_array[ϕ_i, τ_i] = ϕ
            τ_array[ϕ_i, τ_i] = τ

            K = construct_K(τ, ϕ)
            J = construct_J(τ, ϕ)
            
            ϕ_inv = 1/ϕ
            
            C = Float64(τ/(K[1] + 1 - τ))


            if ϕ >= C
                map_type_masks["0"][ϕ_i, τ_i] = K[1]*ϕ ∈ J
                map_type_masks["1a"][ϕ_i, τ_i] = sum(K[1]*ϕ-1 .< J .< K[1]*ϕ) >=1
            else
                map_type_masks["1b"][ϕ_i, τ_i] = J[end]/K[2] <= ϕ <= J[1]/K[1]
                map_type_masks["2a"][ϕ_i, τ_i] = (J[1]/K[1] < ϕ) && (J[end]/K[2] <= ϕ)
                map_type_masks["2b"][ϕ_i, τ_i] =  (ϕ < J[end]/K[2]) && (ϕ <= J[1]/K[1])
                map_type_masks["3"][ϕ_i, τ_i] = J[1]/K[1] < ϕ < J[end]/K[2]
            end
        end

    end

    return map_type_masks, ϕ_array, τ_array

end


function subscriptnumber(i::Int)
    if i < 0
        c = [Char(0x208B)]
    else
        c = []
    end
    for d in reverse(digits(abs(i)))
        push!(c, Char(0x2080+d))
    end
    return join(c)
end

function superscriptnumber(i::Int)
    if i < 0
        c = [Char(0x207B)]
    else
        c = []
    end
    for d in reverse(digits(abs(i)))
        if d == 0 push!(c, Char(0x2070)) end
        if d == 1 push!(c, Char(0x00B9)) end
        if d == 2 push!(c, Char(0x00B2)) end
        if d == 3 push!(c, Char(0x00B3)) end
        if d > 3 push!(c, Char(0x2070+d)) end
    end
    return join(c)
end



function determine_branch_symbols(vector_fields)
    """
    Determine the alphabet.
    """

    L_counter = 0
    R_counter = 0

    branch_symbols = Vector{String}()
    for f in vector_fields

        if f > 0 
            push!(branch_symbols, "L"*subscriptnumber(L_counter + 1))
            L_counter += 1
        else
            push!(branch_symbols, "R"*subscriptnumber(R_counter + 1))
            R_counter +=1
        end
    
    end

    return branch_symbols
end

function convert_Zₖ_sequence_to_string(Zₖ_sequence)

    Zₖ_string = ""
    for k in Zₖ_sequence
        Zₖ_string = Zₖ_string*"Z"*subscriptnumber(k)*"-"
    end

    return String(collect(Zₖ_string)[1:end-1])
end

function compute_Z_index(Z, interval, indices)
    return collect(keys(Z))[map(k -> interval ∈ Z[k], indices)]
end


function is_member_open(x, set_bounds)
    return set_bounds[1] < x < set_bounds[2]
end

function is_member_closed(x, set_bounds)
    return set_bounds[1] <= x <= set_bounds[2]
end
function compute_inverse_branches(parts, vector_fields, no_branches)
    """
    This function computes the inverse branches of the map.
    """

    no_branches = size(parts)[1]
    inverse_branches = zeros(Rational,no_branches+1, 2)
    
    for i in range(1, no_branches)
    
        inverse_branches[i, :] = [parts[i, 1] + vector_fields[i], parts[i, 2] + vector_fields[i]]
        
    end

    inverse_branches[1:no_branches, :] = sort(inverse_branches[1:no_branches, :], dims=1)

    inverse_branches[end, :] = [1//1, 1//1]

    return inverse_branches
end

function compute_Zₖ(inverse_branches, k)
    """
    This function computes the intervals in I which have k preimages. 
    """

    # Compute unique pairs of intervals
    interval_iterator = combinations(inverse_branches, k)

    # Compute the sequential intersections. 
    Zₖ = map(x -> compute_sequential_intersection(x), interval_iterator)
    
    # Filter out empty intervals
    filter!(x->!isempty(x), Zₖ)

    return Zₖ

end


function compute_all_Zₖ(inverse_branches, I)

    # Convert to interval union.
    inverse_branches = intervalUnion.(inverse_branches)
    
    # Compute number of inverse branches.
    no_inv_branches = length(inverse_branches)

    # Construct empty dict.
    Z = Dict()

    total_union = compute_sequential_union(inverse_branches)

    # Compute the maximum Zₖ.
    if isempty(compute_Zₖ(inverse_branches, no_inv_branches))
        Z[no_inv_branches] = ∅
    else
        Z[no_inv_branches] = remove_eps_less_intervals(compute_Zₖ(inverse_branches, no_inv_branches)[1])
    end

    for k in reverse(1:no_inv_branches-1)
        Zₖ = compute_Zₖ(inverse_branches, k)
        Z[k] = remove_empty_sets(map(inter -> compute_recursive_set_diff(inter, Z, k, no_inv_branches), Zₖ))
    end
    
    if sum(map(x-> diam(x), setdiff(I, total_union).v)) > eps(1.0)
        Z₀ = setdiff(I, total_union)
        Z[0] = remove_empty_sets([map(interval -> compute_recursive_set_diff(interval, Z, 0, no_inv_branches), Z₀)])
    end

    filter!(p -> !isempty(p.second), Z)

    Z = Dict(keys(Z).=> [convert_union_to_intervals(Z[i]) for i in keys(Z)])

    return Z
end


function remove_empty_sets(intervals)
    filter!(x->!isempty(x), intervals)
end

function convert_union_to_intervals(interval_unions)
    
    intervals = Vector{Interval{Float64}}()

    for interval_union in interval_unions
        no_intervals = length(interval_union.v)
        for i in range(1, no_intervals)
            push!(intervals, popat!(interval_union.v, 1))
        end
    end
    return intervals
end

function remove_eps_less_intervals(interval_union)
    """
    This function removes any intervals that have diameter less than 
    machine precision ϵ.
    """
    intervalUnion(deleteat!(interval_union.v, (i for i in eachindex(interval_union.v) if diam((interval_union.v[i] ∩ I)[1] )<eps(1.0) )))
end

function compute_recursive_set_diff(inter, Z, k, no_inv_branches)

    for i in reverse(k+1:no_inv_branches)
        
        previous_union = compute_sequential_union(Z[i])
        inter = remove_eps_less_intervals(isempty(previous_union) ? inter : setdiff(inter, previous_union))
    end

    return inter
end

function compute_sequential_union(intervals)
    """
    This function computes the sequential union of n intervals.
    """
    sequential_union = intervals[1]
    for i in 2:length(intervals)
        sequential_union = sequential_union ∪ intervals[i]
    end
    
    return sequential_union
end


function compute_sequential_intersection(intervals)
    """
    This function computes the sequential intersection of n intervals.
    """
    sequential_intersection = intervalUnion(interval(-1//1, 1//1))
    for interval in intervals
        sequential_intersection = sequential_intersection ∩ interval
    end

    return sequential_intersection
end





function approximate_ceil(x)
    γ = 1
    return 1 + x - 0.5 - atan(-γ*sin(2*π*x)/(1-γ*cos(2*π*x)))/π
end





function construct_range(min_value, max_value, no_points)

    value_range = max_value - min_value
    no_range_points = Int(floor(no_points/value_range))
    array = []
    for i in range(min_value, max_value)
        array = vcat(array, LinRange(i, i+1, no_range_points))
    end

    return array
end



function plot_heatmap(τ_array, ϕ_array, colour_map, map_types, no_ticks)
        
    τ_range = τ_array[1,: ]
    ϕ_range = ϕ_array[:, 1]
    fig, ax, hm = GLMakie.heatmap(
                                map_types', 
                                colormap=colour_map,     
                                transparency=false,
                                axis=(
                                    xticks=(compute_ticks_and_tick_labels(τ_range, no_ticks, 0)), 
                                    yticks=(compute_ticks_and_tick_labels(ϕ_range, no_ticks, 1)),
                                    xlabel="τ",
                                    ylabel="ϕ",
                                    xlabelsize=30,
                                    ylabelsize=30
                                ),                  
    )        
    cbar = Colorbar(fig[1, 2], hm, label="Map type")
    #cbar.ticks = ([1, 2, 3, 4, 5, 6], ["0", "1a", "1b", "2a", "2b", "3"])

    fig
    return fig, ax, hm

end



function compute_ticks_and_tick_labels(array, no_ticks, dp)

    tick_values = LinRange(array[1], array[end], no_ticks)
    tick_labels = string.(round.(tick_values, digits=dp))
    
    return tick_values, tick_labels

end



function map_type_plot(ax, map_type_arrays, map_type_label, colour)
    """
    This function plots the map types.
    """
    
    τ_array = map_type_arrays["τ_array"]
    ϕ_array = map_type_arrays["ϕ_array"]

    if map_type_label == "1a"

        map_type_array = map_type_arrays["map_type_"*map_type_label]
        
        indices = findall(vcat(map_type_array...).==1)
        τ_points = vcat(τ_array...)[indices]
        ϕ_points = vcat(ϕ_array...)[indices]
        
        scatter!(
            ax,
            τ_points,
            ϕ_points,
            markersize=2,
            color=colour,
            transparency=false,
            label=map_type_label,
        )
    elseif map_type_label=="0"

        τ_points = map_type_arrays["map_type_"*map_type_label*"_τ"]
        ϕ_points = map_type_arrays["map_type_"*map_type_label*"_ϕ"]

        scatter!(
            τ_points,
            ϕ_points,
            markersize=2,
            color=colour,
            label=map_type_label
        )


    else

        map_type_array = map_type_arrays["map_type_"*map_type_label]
        
        indices = findall(vcat(map_type_array...).==1)
        τ_points = vcat(τ_array...)[indices]
        ϕ_points = vcat(ϕ_array...)[indices]

        scatter!(
            τ_points,
            ϕ_points,
            markersize=2,
            color=colour,
            label=map_type_label
        )
    end
end


function compute_continuous_map_types(τ_min, τ_max, ϕ_min, ϕ_max, no_τ, no_ϕ)
    """
    This function computes all the parameter space points with a given map type
    within a range of ϕ and τ.
    """
    τ_range = construct_range(τ_min, τ_max, no_τ)
    ϕ_range = construct_range(ϕ_min, ϕ_max, no_ϕ)
    no_τ = length(τ_range)
    no_ϕ = length(ϕ_range)
    ϕ_array = zeros(no_ϕ, no_τ)
    τ_array = zeros(no_ϕ, no_τ)
    map_type_1a = zeros(no_ϕ, no_τ)
    map_type_1b = zeros(no_ϕ, no_τ)
    map_type_2a = zeros(no_ϕ, no_τ)
    map_type_2b = zeros(no_ϕ, no_τ)
    map_type_3 = zeros(no_ϕ, no_τ)

    line_1 = zeros(no_ϕ, no_τ)
    line_2 = zeros(no_ϕ, no_τ)

    map_type_0_ϕ = []
    map_type_0_τ = []

    for τ in τ_range
        points = compute_map_type_0(τ, ϕ_max+1)
        for ϕ in points
            push!(map_type_0_ϕ, ϕ)
            push!(map_type_0_τ, τ)
        end
    end

    map_type_0_ϕ = Float64.(map_type_0_ϕ)
    map_type_0_τ = Float64.(map_type_0_τ)


    for τ_i in eachindex(τ_range)
        
        τ = τ_range[τ_i]
        
        for ϕ_i in eachindex(ϕ_range)

            ϕ = ϕ_range[ϕ_i]

            ϕ_array[ϕ_i, τ_i] = ϕ
            τ_array[ϕ_i, τ_i] = τ
        
            ϕ_inv = 1/ϕ
        
            K = construct_K(τ, ϕ)
            J = construct_J(τ, ϕ)

            C = (τ/(K[1] + 1 - τ))


            if ϕ >= C
                if mid.(K[1]*ϕ) ∉ J
                    map_type_1a[ϕ_i, τ_i] =true
                end
            elseif ϕ<C && mid.(K[1]*ϕ) ∉ J
                map_type_1b[ϕ_i, τ_i] = precise_less_than_or_equal(ϕ, J[1]/K[1]) && precise_more_than_or_equal(ϕ, J[end]/K[2])        
                map_type_2a[ϕ_i, τ_i] = (ϕ > J[1]/K[1]) && precise_more_than_or_equal(ϕ, J[end]/K[2]) 
                map_type_2b[ϕ_i, τ_i] = precise_less_than_or_equal(ϕ, J[1]/K[1]) && (ϕ < J[end]/K[2])   
                map_type_3[ϕ_i, τ_i] = J[1]/K[1] < ϕ < J[end]/K[2]
                line_1[ϕ_i, τ_i] = precise_equal(ϕ, J[1]/K[1])
                line_2[ϕ_i, τ_i] = precise_less_than_or_equal(ϕ, J[1]/K[1])

            end
        end
    end

    map_type_arrays = Dict(
                        "map_type_0_τ" => map_type_0_τ,
                        "map_type_0_ϕ" => map_type_0_ϕ,
                        "map_type_1a" => map_type_1a,
                        "map_type_1b" => map_type_1b,
                        "map_type_2a" => map_type_2a,
                        "map_type_2b" => map_type_2b,
                        "map_type_3" => map_type_3,
                        "ϕ_array" => ϕ_array,
                        "τ_array" => τ_array,
                        "map_type_line_1" => line_1,
                        "map_type_line_2" => line_2
    )
    return map_type_arrays
end

  

function compute_g₁(τ, ϕ)
    """
    This functions checks if the g₄ equality is satisfied 
    for a given τ, ϕ.
    """
    j₁ = compute_j₁_float(τ, ϕ)
    jₘ = compute_jₘ_float(τ, ϕ)
    ϕ_k₁ = ϕ*compute_k₂_float(τ, ϕ)
    M = jₘ - j₁
    
    count = 0
    count += ϕ_k₁ - j₁ == 0.0
    #for c in range(0, )
    #    count += ϕ_k₁ - jₘ - c == 0.0
    #end

    return count

end

function compute_g₁(τ, ϕ)
    """
    This functions checks if the g₄ equality is satisfied 
    for a given τ, ϕ.
    """
    j₁ = compute_j₁_float(τ, ϕ)
    ϕ_k₁ = ϕ*compute_k₁_float(τ, ϕ)
    
    count = ϕ_k₁ - j₁ == 0.0

    return count

end

function compute_g₂(τ, ϕ)
    """
    """
    j₁ = compute_j₁_float(τ, ϕ)
    jₘ = compute_jₘ_float(τ, ϕ)
    ϕ_k₁ = ϕ*compute_k₁_float(τ, ϕ)
    M = jₘ - j₁
    
    count = 0
    for c in range(0, M)
        count += ϕ_k₁ - j₁ - c == 0.0
    end

    return count

end

function compute_g₃(τ, ϕ)
    """
    This functions checks if the g₄ equality is satisfied 
    for a given τ, ϕ.
    """
    jₘ = compute_jₘ_float(τ, ϕ)
    ϕ_k₂ = ϕ*compute_k₂_float(τ, ϕ)

    count = ϕ_k₂ - jₘ == 0.0

    return count

end

function compute_g_roots(τ_range, ϕ_range, g_i)
    """

    """
    ϕ_points = []
    τ_points = []

    τ_min = minimum(τ_range)
    τ_max = maximum(τ_range)
    ϕ_min = minimum(ϕ_range)
    ϕ_max = maximum(ϕ_range)

    min_jₘ = compute_j₁_float(τ_min, ϕ_min)
    max_jₘ = compute_jₘ_float(τ_max, ϕ_max)
    min_k₁ = compute_k₁_float(τ_min, ϕ_max)
    max_k₁ = compute_k₁_float(τ_max, ϕ_min)

    ϕ_range = []
    for num=min_jₘ:max_jₘ, den=min_k₁:max_k₁
        if num/den > 1 && num/den <= ϕ_max
            push!(ϕ_range, num/den)
        end
    end

    ϕ_range = Float64.(ϕ_range)

    if g_i=="g_1"
        for ϕ in ϕ_range
            for τ in τ_range[findall(x->x==1, compute_g₁.(τ_range, ϕ))]
                push!(ϕ_points, ϕ)
                push!(τ_points, τ)
            end
        end
    elseif g_i=="g_2"
        for ϕ in ϕ_range
            for τ in τ_range[findall(x->x==1, compute_g₂.(τ_range, ϕ))]
                push!(ϕ_points, ϕ)
                push!(τ_points, τ)
            end
        end
    else
        for ϕ in ϕ_range
            for τ in τ_range[findall(x->x==1, compute_g₃.(τ_range, ϕ))]
                push!(ϕ_points, ϕ)
                push!(τ_points, τ)
            end
        end
    end

    return Float64.(τ_points), Float64.(ϕ_points)
        
end 


function compute_k₁(τ, ϕ)

    inner = τ*(1//1 + 1//ϕ)
    return ceil(inner) + isinteger(inner)

end

function compute_k₂(τ, ϕ)
    inner = τ*(1//1 + 1//ϕ) + 1//ϕ
    return ceil(inner)
end

function compute_j₁(τ, ϕ)

    inner = τ*(1//1 + ϕ)
    return ceil(inner) 

end

function compute_jₘ(τ, ϕ)

    inner = τ*(1//1+ϕ) + ϕ
    return ceil(inner) 

end



function compute_k₁_float(τ, ϕ)

    inner = τ*(1.0 + 1/ϕ)
    return ceil(inner) + isinteger(inner)

end

function compute_k₂_float(τ, ϕ)
    return compute_k₁_float(τ, ϕ) + 1.0
end

function compute_j₁_float(τ, ϕ)

    inner = τ*(1.0 + ϕ)
    return ceil(inner) + isinteger(inner)

end

function compute_jₘ_float(τ, ϕ)

    inner = τ*(1.0+ϕ) + ϕ
    return ceil(inner) + isinteger(inner)

end


function initialise_plot()

    fig = CairoMakie.Figure(resolution = (1400, 1100), font = "Avenir Light")
    fontsize_theme = CairoMakie.Theme(fontsize = 30)
    CairoMakie.set_theme!(fontsize_theme)
    ax = CairoMakie.Axis(
        fig[1, 1],
    )
    return fig, ax
end

function compute_t_roots(func, interval)

    g_roots = roots(func, interval, Bisection)

    t_points = []
    if !isempty(g_roots)
        for root in g_roots
            push!(t_points, mid(root.interval))
        end
    end

    return t_points

end

#### TODO - Modifying this - Can we extend the function slightly so that we are sure we have a root?


function compute_roots(func, interval, τ, τ_points, ϕ_points, type)

    if type=="Newton"
        g_roots = roots(func, interval, Newton)
    else
        g_roots = roots(func, interval, Bisection)
    end
    if !isempty(g_roots)
        for root in g_roots
            push!(τ_points, τ)
            push!(ϕ_points, mid(root.interval))
        end
    end

end




function plot_bifurcation_curve(τ_points, ϕ_points, label)

    fig = CairoMakie.Figure()
    fontsize_theme = CairoMakie.Theme(fontsize = 30)
    CairoMakie.set_theme!(fontsize_theme)
    ax = CairoMakie.Axis(
        fig[1, 1],
        xlabel="τ",
        ylabel="ϕ",
        xlabelsize=40,
        ylabelsize=40,
    )


    CairoMakie.scatter!(  
        τ_points,
        ϕ_points,
        color=:black,
        markersize=2,
        transparency=false,
    )

    save(label, fig; pt_per_unit=2)
end


function plot_map_types(axis, map_data, ϕ_array, τ_array, colour)

    indices = findall(vcat(map_data...).==1)

    τ_points = vcat(τ_array...)[indices]
    ϕ_points = vcat(ϕ_array...)[indices]

    map_plot_object = CairoMakie.scatter!(
        axis, 
        τ_points,
        ϕ_points,
        color=colour,
        markersize=1.75,
        transparency=false,
    )
    return map_plot_object
end


function plot_full_bifurcation_diagram(data, parameters, overlay, τ_max, ϕ_max)

    τ_array = parameters[1]
    ϕ_array = parameters[2]

    indices = findall(vcat(data[1]...).==1)
    τ_points = vcat(τ_array...)[indices]
    ϕ_points = vcat(ϕ_array...)[indices]

    colours = cgrad(:batlow, 80, categorical = true, scale= :exponential)
    colours_2 = cgrad(:roma, 30, categorical = true, scale= :exponential)

    fig = CairoMakie.Figure(resolution = (1400, 1100), font = "Avenir Light")
    fontsize_theme = CairoMakie.Theme(fontsize = 30)
    CairoMakie.set_theme!(fontsize_theme)
    ax = CairoMakie.Axis(
        fig[1, 1],
        xlabel="τ",
        ylabel="ϕ",
        xlabelsize=50,
        xticklabelsize=50,
        ylabelsize=50,
        yticklabelsize=50
    )


    
    map_2 = CairoMakie.scatter!(
                                ax,
                                τ_points,
                                ϕ_points,
                                color=colours[80],
                                markersize=1.75,
                                transparency=false,
    )


    map_3 = plot_map_types(ax, data[2], ϕ_array, τ_array, colours[20])
    map_4 = plot_map_types(ax, data[3], ϕ_array, τ_array, colours[40])
    map_5 = plot_map_types(ax, data[4], ϕ_array, τ_array, colours[60])
    map_6 = plot_map_types(ax, data[5], ϕ_array, τ_array, colours[1])


    if overlay

        indices = findall(vcat(data[6]...).==1)
        τ_points = vcat(τ_array...)[indices]
        ϕ_points = vcat(ϕ_array...)[indices]


        sample_indices = sample(1:length(indices), 50000)

        CairoMakie.scatter!(
            ax, 
            τ_points[sample_indices],
            ϕ_points[sample_indices],
            color=:black,
            markersize=1,
            transparency=true,
        )

    end



    colors = [:black, colours[80], colours[20], colours[40], colours[60], colours[1], :red]

    labels = [  
                "M(τ, ϕ) = 1",
                "M(τ, ϕ) = 2",
                "M(τ, ϕ) = 3",
                "M(τ, ϕ) = 4",
                "M(τ, ϕ) = 5",
                "M(τ, ϕ) = 6",
                "M(τ, ϕ) = 7"
    ]
    group_color = [CairoMakie.PolyElement(color = color, strokecolor = :transparent)
        for color in colors]
    
    legend = CairoMakie.Legend(
                                fig[1,2],
                                group_color,
                                labels
    )

    ϕ_range = LinRange(1, ϕ_max, 10000)
    τ_range = LinRange(1, τ_max, 10000)

    τ_points_z, ϕ_points_z = compute_zeros(τ_range, ϕ_range)

    identity_map_indices(τ, ϕ) = compute_k₁(τ, ϕ)*ϕ >  1 + τ*(ϕ+1)

    g₂_id_indices = identity_map_indices.(τ_points_z["g₂"], ϕ_points_z["g₂"])
    
    CairoMakie.scatter!(
        ax,
        τ_points_z["g₁"],
        ϕ_points_z["g₁"],
        color=:red,
        markersize=2,
        transparency=false,
    )
    
    CairoMakie.scatter!(
        ax,
        τ_points_z["g₂"][g₂_id_indices],
        ϕ_points_z["g₂"][g₂_id_indices],
        color=:black,
        markersize=2,
        transparency=false,
    )
    
    CairoMakie.scatter!(
        ax,
        τ_points_z["g₃"],
        ϕ_points_z["g₃"],
        color=:red,
        markersize=2,
        transparency=false,
    )

    return fig, ax
end



function compute_full_bifurcation_diagram(ϕ_range, τ_range, no_sample_points, condition)

    map_2_data = zeros(no_sample_points, no_sample_points)
    map_3_data = zeros(no_sample_points, no_sample_points)
    map_4_data = zeros(no_sample_points, no_sample_points)
    map_5_data = zeros(no_sample_points, no_sample_points)
    map_6_data = zeros(no_sample_points, no_sample_points)
    map_g_data = zeros(no_sample_points, no_sample_points)
    ϕ_array = zeros(no_sample_points, no_sample_points) 
    τ_array = zeros(no_sample_points, no_sample_points) 

    for τ_i in eachindex(τ_range)
        for ϕ_i in eachindex(ϕ_range)
            
            τ = τ_range[τ_i]
            ϕ = ϕ_range[ϕ_i]
            
            M_1 = compute_jₘ(τ, ϕ) - compute_j₁(τ, ϕ) 
            map_2_data[τ_i, ϕ_i] = compute_k₁(τ, ϕ)*ϕ > 1 + τ*(ϕ+1)
            map_3_data[τ_i, ϕ_i] = (compute_k₁(τ, ϕ)*ϕ < compute_j₁(τ, ϕ)) && (compute_k₂(τ, ϕ)*ϕ > compute_jₘ(τ, ϕ)) && (compute_k₁(τ, ϕ)*ϕ < 1 + τ*(ϕ+1) )
            map_4_data[τ_i, ϕ_i] = (compute_k₁(τ, ϕ)*ϕ > compute_j₁(τ, ϕ)) && (compute_k₂(τ, ϕ)*ϕ > compute_jₘ(τ, ϕ)) && (compute_k₁(τ, ϕ)*ϕ < 1 + τ*(ϕ+1) )
            map_5_data[τ_i, ϕ_i] = (compute_k₁(τ, ϕ)*ϕ < compute_j₁(τ, ϕ)) && (compute_k₂(τ, ϕ)*ϕ < compute_jₘ(τ, ϕ)) && (compute_k₁(τ, ϕ)*ϕ < 1 + τ*(ϕ+1) )
            map_6_data[τ_i, ϕ_i] = (compute_k₁(τ, ϕ)*ϕ > compute_j₁(τ, ϕ)) && (compute_k₂(τ, ϕ)*ϕ < compute_jₘ(τ, ϕ)) && (compute_k₁(τ, ϕ)*ϕ < 1 + τ*(ϕ+1) )

            if condition == "greater than"
                map_g_data[τ_i, ϕ_i] = ϕ > M_1
            else
                map_g_data[τ_i, ϕ_i] = ϕ < M_1
            end

            ϕ_array[τ_i, ϕ_i] = ϕ
            τ_array[τ_i, ϕ_i] = τ

        end
    end

    data = [map_2_data, map_3_data, map_4_data, map_5_data, map_6_data, map_g_data]
    parameters = [τ_array, ϕ_array]

    return data, parameters
end



function compute_empirical_transition_matrix(partition_list)
    """
    This function computes the empirical transition matrix from the list of visited partitions.
    """
    state = partition_list
    
    state_lag = circshift(partition_list, 1)
    
    df = DataFrame(state=state, state_lag=state_lag)

    df = DataFrames.transform(DataFrames.groupby(df, [:state_lag, :state]), nrow => :transition_count)

    df = DataFrames.transform(DataFrames.groupby(df, :state_lag), nrow => :state_count) 

    df[:,:prob] = df[:,:transition_count]./df[:,:state_count]

    trans = unique(df[:,[:state_lag, :state, :prob]])

    no_states = length(keys(counter(partition_list)))
    state_counts = zeros(no_states)
    state_count = df[!, "state_count"]
    
    Q = zeros(no_states, no_states)

    for row_index in 1:size(trans)[1]
        row = trans[row_index, :]
        Q[row["state_lag"], row["state"]] = row["prob"]
    end

    for i in range(1, no_states)
        
        if isempty(findall(df[!, "state_lag"].==i))
            state_counts[i] = 0
        else
            state_counts[i] = state_count[findall(df[!, "state_lag"].==i)[1]]
        end

    end

    return Q, state_counts
end


function compute_bisections(τ_range, ϕ_range)


    k₁_jump_τ_points = Float64[]
    k₁_jump_ϕ_points = Float64[]
    k₁_jump_τ_points = Float64[]
    k₁_jump_ϕ_points = Float64[]
    j₁_jump_τ_points = Float64[]
    j₁_jump_ϕ_points = Float64[]
    jₘ_jump_τ_points = Float64[]
    jₘ_jump_ϕ_points = Float64[]

    interval = Interval(minimum(ϕ_range), maximum(ϕ_range))

    for τ in τ_range
    
        k₁_jump(ϕ) = 1 - (τ*(1/ϕ+1) -floor(τ*(1/ϕ+1)))
        k₂_jump(ϕ) = 1 - (τ*(1/ϕ+1) + 1/ϕ -floor(τ*(1/ϕ+1) + 1/ϕ))

        j₁_jump(ϕ) = 1 - (τ*(ϕ+1) -floor(τ*(ϕ+1)))
        jₘ_jump(ϕ) = 1 - (τ*(ϕ+1)+ϕ - floor(τ*(ϕ+1) + ϕ))
      
        compute_roots(k₁_jump, interval, τ, k₁_jump_τ_points, k₁_jump_ϕ_points, "Bisection")
        compute_roots(k₂_jump, interval, τ, k₂_jump_τ_points, k₂_jump_ϕ_points, "Bisection")

        compute_roots(j₁_jump , interval, τ, j₁_jump_τ_points, j₁_jump_ϕ_points, "Bisection")
        compute_roots(jₘ_jump , interval, τ, jₘ_jump_τ_points, jₘ_jump_ϕ_points, "Bisection")

    end


    τ_points = Dict(
                    "k₁_jump"=> k₁_jump_τ_points,
                    "k₂_jump"=> k₂_jump_τ_points,
                    "j₁_jump"=> j₁_jump_τ_points,
                    "jₘ_jump"=> jₘ_jump_τ_points,
    )
    ϕ_points = Dict(
                    "k₁_jump"=> k₁_jump_ϕ_points,
                    "k₂_jump"=> k₂_jump_τ_points,
                    "j₁_jump"=> j₁_jump_ϕ_points,
                    "jₘ_jump"=> jₘ_jump_ϕ_points,
    )
    return τ_points, ϕ_points
end

function plot_G_jumps(τ_points, ϕ_points)

    fig = CairoMakie.Figure(resolution = (1400, 1100), font = "Avenir Light")
    fontsize_theme = CairoMakie.Theme(fontsize = 60)
    CairoMakie.set_theme!(fontsize_theme)
    ax = CairoMakie.Axis(
        fig[1, 1],
        xlabel="τ",
        ylabel="ϕ",
        xlabelsize=50,
        xticklabelsize=50,
        ylabelsize=50,
        yticklabelsize=50
    )

    CairoMakie.scatter!(
        ax,
        τ_points["j₁_jump"],
        ϕ_points["j₁_jump"],
        color=:red,
        markersize=3,
        transparency=false,
    )
    CairoMakie.scatter!(
        ax,
        τ_points["jₘ_jump"],
        ϕ_points["jₘ_jump"],
        color=:blue,
        markersize=3,
        transparency=false,
    )

    CairoMakie.scatter!(
        ax,
        τ_points["k₁_jump"],
        ϕ_points["k₁_jump"],
        color=:black,
        markersize=3,
        transparency=false,
    )

    colors = [:black, :red, :blue]
    
    labels = [  
                L"G_{k_1}",
                L"G_{j_1}",
                L"G_{j_m}",
    ]

    group_color = [CairoMakie.PolyElement(color = color, strokecolor = :transparent)
        for color in colors]

    legend = CairoMakie.Legend(
                                fig[1,2],
                                group_color,
                                labels
    )

    return fig, ax
end




function plot_G_jumps(τ_points, ϕ_points)

    fig = CairoMakie.Figure(resolution = (1400, 1100), font = "Avenir Light")
    fontsize_theme = CairoMakie.Theme(fontsize = 60)
    CairoMakie.set_theme!(fontsize_theme)
    ax = CairoMakie.Axis(
        fig[1, 1],
        xlabel="τ",
        ylabel="ϕ",
        xlabelsize=50,
        xticklabelsize=50,
        ylabelsize=50,
        yticklabelsize=50
    )

    CairoMakie.scatter!(
        ax,
        τ_points["j₁_jump"],
        ϕ_points["j₁_jump"],
        color=:red,
        markersize=3,
        transparency=false,
    )
    CairoMakie.scatter!(
        ax,
        τ_points["jₘ_jump"],
        ϕ_points["jₘ_jump"],
        color=:blue,
        markersize=3,
        transparency=false,
    )

    CairoMakie.scatter!(
        ax,
        τ_points["k₁_jump"],
        ϕ_points["k₁_jump"],
        color=:black,
        markersize=3,
        transparency=false,
    )

    colors = [:black, :red, :blue]
    
    labels = [  
                L"G_{k_1}",
                L"G_{j_1}",
                L"G_{j_m}",
    ]

    group_color = [CairoMakie.PolyElement(color = color, strokecolor = :transparent)
        for color in colors]

    legend = CairoMakie.Legend(
                                fig[1,2],
                                group_color,
                                labels
    )

    return fig, ax
end




function add_bifurcation_lines(fig, ax, τ_min, τ_max, ϕ_min, ϕ_max)

    
    τ_min = floor(τ_min)
    τ_max = ceil(τ_max)
    ϕ_min = floor(ϕ_min)
    ϕ_max = ceil(ϕ_max)
    τ_points_b, ϕ_points_b = compute_bisections(
                                                LinRange(Float64(τ_min), Float64(τ_max), 10000),
                                                LinRange(Float64(ϕ_min), Float64(ϕ_max), 10000)
    )
    τ_points_z_0, ϕ_points_z_0, τ_points_z_0_d, ϕ_points_z_0_d = compute_zeros(τ_min,τ_max, 
                                                                                ϕ_min, ϕ_max)
        
    
    τ_points_z = vcat(τ_points_z_0, τ_points_z_0_d) 
    ϕ_points_z = vcat(ϕ_points_z_0, ϕ_points_z_0_d) 

   
    τ_indices = τ_min .< τ_points_z .< τ_max
    ϕ_indices = ϕ_min .< ϕ_points_z .< ϕ_max

    if !isempty(τ_indices)
        CairoMakie.scatter!(
                    ax_1,
                    τ_points_z[τ_indices .+ ϕ_indices .== 2],
                    ϕ_points_z[τ_indices .+ ϕ_indices .== 2],
                    color=:black,
                    markersize=4,
                    transparency=false,
        )
    end


    for key in keys(τ_points_b) 
        if !isempty(τ_points_b[key])

            τ_points = τ_points_b[key]
            ϕ_points = ϕ_points_b[key]
                    
            indices = compute_k₁_float.(τ_points, ϕ_points).*ϕ_points .<= 1 .+ τ_points.*(ϕ_points .+ 1)
                    
            τ_indices = τ_min .< τ_points .< τ_max
            ϕ_indices = ϕ_min .< ϕ_points .< ϕ_max

            
            if !isempty(indices)
                CairoMakie.scatter!(
                    ax,
                    τ_points[τ_indices .+ ϕ_indices .+ indices .== 3],
                    ϕ_points[τ_indices .+ ϕ_indices .+ indices .== 3],
                    color=:black,
                    markersize=4,
                    transparency=false,
                )
            end
        end
    end
    return fig, ax
end


function compute_k₁_float(τ, ϕ)

    inner = τ*(1.0+1.0/ϕ)
    return ceil(inner) + isinteger(inner)

end

function compute_bisections(τ_range, ϕ_range)


    k₁_jump_τ_points = Float64[]
    k₁_jump_ϕ_points = Float64[]
    k₂_jump_τ_points = Float64[]
    k₂_jump_ϕ_points = Float64[]
   
    j₁_jump_τ_points = Float64[]
    j₁_jump_ϕ_points = Float64[]
    jₘ_jump_τ_points = Float64[]
    jₘ_jump_ϕ_points = Float64[]

    interval = Interval(minimum(ϕ_range), maximum(ϕ_range))

    for τ in τ_range
    
        k₁_jump(ϕ) = 1 - (τ*(1/ϕ+1) -floor(τ*(1/ϕ+1)))
        k₂_jump(ϕ) = 1 - (τ*(1/ϕ+1) + 1/ϕ -floor(τ*(1/ϕ+1) + 1/ϕ))

        j₁_jump(ϕ) = 1 - (τ*(ϕ+1) -floor(τ*(ϕ+1)))
        jₘ_jump(ϕ) = 1 - (τ*(ϕ+1)+ϕ - floor(τ*(ϕ+1) + ϕ))
      
        compute_roots(k₁_jump, interval, τ, k₁_jump_τ_points, k₁_jump_ϕ_points, "Bisection")
        compute_roots(k₂_jump, interval, τ, k₂_jump_τ_points, k₂_jump_ϕ_points, "Bisection")

        compute_roots(j₁_jump , interval, τ, j₁_jump_τ_points, j₁_jump_ϕ_points, "Bisection")
        compute_roots(jₘ_jump , interval, τ, jₘ_jump_τ_points, jₘ_jump_ϕ_points, "Bisection")

    end


    τ_points = Dict(
                    "k₁_jump"=> k₁_jump_τ_points,
                    "k₂_jump"=> k₂_jump_τ_points,
                    "j₁_jump"=> j₁_jump_τ_points,
                    "jₘ_jump"=> jₘ_jump_τ_points,
    )
    ϕ_points = Dict(
                    "k₁_jump"=> k₁_jump_ϕ_points,
                    "k₂_jump"=> k₂_jump_ϕ_points,
                    "j₁_jump"=> j₁_jump_ϕ_points,
                    "jₘ_jump"=> jₘ_jump_ϕ_points,
    )
    return τ_points, ϕ_points
end


function compute_zeros(τ_min, τ_max, ϕ_min, ϕ_max)

    τ_range = LinRange(1, τ_max, 10000)

    k₁_min = Int(compute_k₁(τ_min, ϕ_min))
    k₂_min = Int(compute_k₂(τ_min, ϕ_min))

    k₁_max = Int(compute_k₁(τ_max, ϕ_min))
    k₂_max = Int(compute_k₂(τ_max, ϕ_min))

    τ_points_type_0 = []
    ϕ_points_type_0 = []

    τ_points_type_0_d = []
    ϕ_points_type_0_d = []

    for k in k₁_min:k₁_max

        ϕ_values = LinRange(ϕ_min, ϕ_max, Int((ϕ_max - ϕ_min)*k) + 1)

        for ϕ in ϕ_values
            
            indices_0 = findall(compute_k₁.(τ_range, ϕ) .== compute_k₂.(τ_range, ϕ) .&&  k-1 .<= compute_k₁.(τ_range, ϕ) .< k)
            indices_0_d = findall(compute_k₁.(τ_range, ϕ) .!= compute_k₂.(τ_range, ϕ) .&& k - 1 .<= compute_k₁.(τ_range, ϕ) .< k)

            if ϕ - floor(ϕ) > 1e-10
                push!(τ_points_type_0, τ_range[indices_0])
                push!(ϕ_points_type_0, repeat([ϕ], length(indices_0)))
                push!(τ_points_type_0_d, τ_range[indices_0_d])
                push!(ϕ_points_type_0_d, repeat([ϕ], length(indices_0_d)))
            end
        end
    end

    for k in k₂_min:k₂_max 

        ϕ_values = LinRange(ϕ_min, ϕ_max, Int((ϕ_max - ϕ_min)*(k+1) + 1))

        for ϕ in ϕ_values
            indices_0_d = findall(compute_k₁.(τ_range, ϕ) .!= compute_k₂.(τ_range, ϕ) .&& k+1 .<= compute_k₂.(τ_range, ϕ) .< k+2)

            if ϕ - floor(ϕ) > 1e-10
                push!(τ_points_type_0_d, τ_range[indices_0_d])
                push!(ϕ_points_type_0_d, repeat([ϕ], length(indices_0_d)))
            end
        end
    end

    for ϕ in Int(ϕ_min):Int(ϕ_max)
        push!(τ_points_type_0, τ_range)
        push!(ϕ_points_type_0, repeat([ϕ], length(τ_range)))
    end

    τ_points_type_0 = vcat(τ_points_type_0...)
    ϕ_points_type_0 = vcat(ϕ_points_type_0...)
    τ_points_type_0_d = vcat(τ_points_type_0_d...)
    ϕ_points_type_0_d = vcat(ϕ_points_type_0_d...)


    return τ_points_type_0, ϕ_points_type_0, τ_points_type_0_d, ϕ_points_type_0_d
end


function compute_roots(func, interval, τ, τ_points, ϕ_points, type)

    if type=="Newton"
        g_roots = roots(func, interval, Newton)
    else
        g_roots = roots(func, interval, Bisection)
    end
    if !isempty(g_roots)
        for root in g_roots
            push!(τ_points, τ)
            push!(ϕ_points, mid(root.interval))
        end
    end

end