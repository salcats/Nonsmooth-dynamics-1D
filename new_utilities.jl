using IntervalArithmetic
using IntervalUnionArithmetic
using StatsBase
using Combinatorics
using DataStructures
using CairoMakie
using DataFrames

const map_y_tick_label = Dict(
                                2=> ["c₁"],
                                3=> ["c₂"],
                                4=> ["c₃", "c₂"],
                                5=> ["c₂", "c₄"],
                                6=> ["c₃", "c₂", "c₄"]
)

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



function convert_Zₖ_sequence_to_string(Zₖ_sequence)

    Zₖ_string = ""
    for k in Zₖ_sequence
        Zₖ_string = Zₖ_string*"Z"*subscriptnumber(k)*"-"
    end

    return String(collect(Zₖ_string)[1:end-1])
end


function generate_rationals(min_val::Int64, max_val::Int64, max_den::Int64)
    """
    This functions generates all the rational numbers between 
    the minimum and maximum value with bounded denomiantor.
    """
    numbers =  Vector{Rational{Int64}}()

    for den in 1:max_den
        for num in den*min_val:den*max_val
            push!(numbers, num//den)
        end
    end
    return sort!(unique(numbers))
end


function generate_rationals_farey(min_val::Int64, max_val::Int64, max_den::Int64, n::Int64)
    """
    This functions computes the farey numbers between min val and max val
    with denominator smaller than max den.
    """
    farey_numbers = sort!(1 ./ farey(n))
    farey_numbers = farey_numbers[min_val .< farey_numbers .< max_val]
    farey_numbers = farey_numbers[denominator.(farey_numbers) .<= max_den]

    return farey_numbers

end

function find_n_branch_maps(numbers::Vector{Rational{Int64}}, n::Int64)
    """
    This functions finds all maps computed from all pairs in numbers 
    that have n branches.
    """
    τ_values = Vector{Rational{Int64}}()
    ϕ_values = Vector{Rational{Int64}}()
    
    for i in numbers
        for j in numbers
            map_parameters = compute_map(i, j)
            if map_parameters["no_branches"] == n
                push!(τ_values, i)
                push!(ϕ_values, j)
            end
        end
    end
    return τ_values, ϕ_values
end



function compute_k₁(τ::Rational{Int64}, ϕ::Rational{Int64})
    """
    Compute k₁.
    """
    inner = τ*(1//1 + 1//ϕ)
    return Int(ceil(inner) + isinteger(inner))

end

function compute_k₂(τ::Rational{Int64}, ϕ::Rational{Int64})
    """
    Compute k₂.
    """
    inner = τ*(1//1 + 1//ϕ) + 1//ϕ
    return Int(ceil(inner))
end

function compute_j₁(τ::Rational{Int64}, ϕ::Rational{Int64})
    """
    Compute j₁.
    """
    inner = τ*(1//1 + ϕ)
        return Int(ceil(inner))
end

function compute_jₘ(τ::Rational{Int64}, ϕ::Rational{Int64})
    """
    Compute jₘ.
    """
    inner = τ*(1//1+ϕ) + ϕ
    return Int(ceil(inner) )

end


function compute_map_topology(
                                τ::Rational{Int64},
                                ϕ::Rational{Int64},
                                k₁::Int64, 
                                k₂::Int64, 
                                j₁::Int64, 
                                jₘ::Int64, 
    )
    """
    Compute the map topology.
    """
    condition_1 = k₁ == k₂
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



function compute_discontinuity_points(
        τ::Rational{Int64},
        ϕ::Rational{Int64},
        k₁::Int64, 
        k₂::Int64, 
        j₁::Int64, 
        jₘ::Int64, 
    )
    """
    Compute the discontinuity points.
    """
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



function compute_1_map(
        τ::Rational{Int64},
        ϕ::Rational{Int64},
        k₁::Int64, 
        k₂::Int64, 
        j₁::Int64, 
        jₘ::Int64, 
        c_points_dict::Dict
    )

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

function compute_2_map(τ::Rational{Int64},
        ϕ::Rational{Int64},
        k₁::Int64, 
        k₂::Int64, 
        j₁::Int64, 
        jₘ::Int64, 
        c_points_dict::Dict
    )

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

function compute_3_map(
        τ::Rational{Int64},
        ϕ::Rational{Int64},
        k₁::Int64, 
        k₂::Int64, 
        j₁::Int64, 
        jₘ::Int64, 
        c_points_dict::Dict
    )
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

function compute_4_map(
                τ::Rational{Int64},
                ϕ::Rational{Int64},
                k₁::Int64, 
                k₂::Int64, 
                j₁::Int64, 
                jₘ::Int64, 
                c_points_dict::Dict
    )

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

function compute_5_map(
        τ::Rational{Int64},
        ϕ::Rational{Int64},
        k₁::Int64, 
        k₂::Int64, 
        j₁::Int64, 
        jₘ::Int64, 
        c_points_dict::Dict
    )
    
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

function compute_6_map(
        τ::Rational{Int64},
        ϕ::Rational{Int64},
        k₁::Int64, 
        k₂::Int64, 
        j₁::Int64, 
        jₘ::Int64, 
        c_points_dict::Dict
    )
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



function compute_map(τ, ϕ)
    """ 
    This functions computes the map for a given τ and ϕ.
    """
    # Compute the kᵢ and jᵢ functions.
    k₁ = compute_k₁(τ, ϕ)
    k₂ = compute_k₂(τ, ϕ)
    j₁ = compute_j₁(τ, ϕ)
    jₘ = compute_jₘ(τ, ϕ)


    # Compute the map topology.
    map_topology = compute_map_topology(τ, ϕ, k₁, k₂, j₁, jₘ)
    
    # Compute the discontinuity points.
    c_points_dict = compute_discontinuity_points(τ, ϕ, k₁, k₂, j₁, jₘ)

    # Compute the map.
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


function run_map!(x, S, map_parameters)
    """
    Run the map.
    """
    x₀ = map_parameters["x₀"]
    vector_fields = map_parameters["vector_fields"]
    part_start_points = map_parameters["partitions"][:, 1]
    no_samples = map_parameters["no_samples"]
    no_attractor_samples = map_parameters["no_attractor_samples"]

    xₙ = x₀
    xₙ₊₁ = x₀
    
    attractor_starting_point = no_samples - no_attractor_samples

    for i in range(2, attractor_starting_point)
        
        Sᵢ = searchsortedlast(part_start_points, xₙ)
        
        xₙ₊₁ = xₙ + vector_fields[Sᵢ]
 
        xₙ = xₙ₊₁
    end
       
    x[1] = xₙ

    for i in range(2, no_attractor_samples)
        
        Sᵢ = searchsortedlast(part_start_points, x[i-1])
        x[i] =  x[i-1] + vector_fields[Sᵢ]
        S[i] = Sᵢ
    end
end


function run_map_floating_point(map_parameters)
    """
    Run the map.
    """
    x₀ = Float64.(map_parameters["x₀"])
    vector_fields =  Float64.(map_parameters["vector_fields"])
    part_start_points =  Float64.(map_parameters["partitions"][:, 1])
    no_samples = map_parameters["no_samples"]
    no_attractor_samples = map_parameters["no_attractor_samples"]

    x = Vector{Float64}(undef, no_attractor_samples)
    xₙ = x₀
    xₙ₊₁ = x₀
    
    attractor_starting_point = no_samples - no_attractor_samples

    for i in range(2, attractor_starting_point)
        
        Sᵢ = searchsortedlast(part_start_points, xₙ)
        
        xₙ₊₁ = xₙ + vector_fields[Sᵢ]
 
        xₙ = xₙ₊₁
    end
       
    x[1] = xₙ

    for i in range(2, no_attractor_samples)
        
        Sᵢ = searchsortedlast(part_start_points, x[i-1])
        x[i] =  x[i-1] + vector_fields[Sᵢ]
    end


    S = map(x-> searchsortedlast(part_start_points, x), x)
   
    return x, S

end


function compute_Z₀_intervals(map_parameters)
    """
    Compute the Z₀ intervals.
    """
    vector_fields = map_parameters["vector_fields"]
    parts = map_parameters["partitions"][:, 1]
    no_branches = map_parameters["no_branches"]
    c_points = map_parameters["c_points"]
    
    Z_intervals, _ = compute_Z_intervals(map_parameters)
    
    # Initialise relevant data structures.
    Aₙ = Array{Interval{Rational{Int64}}}(undef, 0)
    Aₙ₊₁ = Array{Interval{Rational{Int64}}}(undef, 0)
    
    Z₀_intervals = Z_intervals["Z₀"]
    Z₁_intervals = Z_intervals["Z₁"]
    Z₂_intervals = Z_intervals["Z₂"]
    
    
    Z₀ = popfirst!(Z₀_intervals)
    
    # Check if any cᵢ lie in Z₀ and then sort the ones that do.
    
    cᵢ_members =  sort(c_points[findall(map(cᵢ -> cᵢ ∈ Z₀, c_points))])
    
    if !isempty(cᵢ_members)
        sorted_cᵢ_members = [inf(Z₀)]
        append!(sorted_cᵢ_members, cᵢ_members)
        append!(sorted_cᵢ_members, sup(Z₀))
    
        # Push split intervals back into Z₀.  
        for i in range(1, length(sorted_cᵢ_members)-1)
            push!(Z₀_intervals, Interval(sorted_cᵢ_members[i], sorted_cᵢ_members[i+1]))
        end
    else
        push!(Z₀_intervals, Z₀)
    end
    
    
    for Z₀ ∈ Z₀_intervals
        # Compute forward image and store it.
        push!(Aₙ, compute_interval_image(Z₀, parts, vector_fields))
    end
    
    
    ################# Step 1 #################
    while !isempty(Aₙ) 
        current_no_intervals = length(Aₙ)
        for i in range(1, current_no_intervals)
        
            A = popfirst!(Aₙ)
        
            # Check if any cᵢ lie in A and then sort the ones that do.
            sorted_cᵢ_members = Vector{Rational{Int64}}()
            push!(sorted_cᵢ_members, inf(A))
            append!(sorted_cᵢ_members, sort!(c_points[findall(map(cᵢ -> cᵢ ∈ A, c_points))]))
            push!(sorted_cᵢ_members, sup(A))
        
            if !isempty(sorted_cᵢ_members)
                # Push new intervals to Aₙ₊₁  
                for i in range(1, length(sorted_cᵢ_members)-1)
                    push!(Aₙ, Interval(sorted_cᵢ_members[i], sorted_cᵢ_members[i+1]))            
                end
            end
        end
        
        ################# Step 2 #################
        current_no_intervals = length(Aₙ)
        
        for i in range(1, current_no_intervals)
            A = popfirst!(Aₙ)
        
            # Check for intersections with other intervals.
            Z₀_intersections = filter!(x -> x != ∅, Z₀_intervals .∩ A)
            Z₀_intersections = filter!(x -> diam(x) > 0, Z₀_intersections)

            Z₁_intersections = filter!(x -> x != ∅, Z₁_intervals .∩ A)
            Z₁_intersections = filter!(x -> diam(x) > 0, Z₁_intersections)
        
            Z₂_intersections = filter!(x -> x != ∅, Z₂_intervals .∩ A)
            Z₂_intersections = filter!(x -> diam(x) > 0, Z₂_intersections)
        
            for intersection ∈ vcat(Z₀_intersections, Z₁_intersections, Z₂_intersections)
                push!(Aₙ, intersection)        
            end
        end
        
        ################# Step 3 #################
        current_no_intervals = length(Aₙ)
        
        for i in range(1, current_no_intervals)
            A = popfirst!(Aₙ)
        
            Z₂_inclusions = findall(A .⊆ Z₂_intervals)
        
            if !isempty(Z₂_inclusions)
                # Extract the Z₂ interval.
                Z₂ = Z₂_intervals[Z₂_inclusions[1]]
        
                # Remove it from the list.
                deleteat!(Z₂_intervals, Z₂_inclusions)
        
                # Add new intervals.
                push!(Z₂_intervals, Interval(inf(Z₂), inf(A))) 
                push!(Z₁_intervals, A)
                push!(Z₂_intervals, Interval(sup(A), sup(Z₂))) 
            else 
                push!(Aₙ, A)
            end
        end
        
        filter!(x-> diam(x) > 0.0, Z₂_intervals)
        
        ################# Step 4 #################
        current_no_intervals = length(Aₙ)
        
        for i in range(1, current_no_intervals)
        
            A = popfirst!(Aₙ)
        
            Z₁_inclusions = findall(A .⊆ Z₁_intervals)
        
            if !isempty(Z₁_inclusions)
        
                # Extract the Z₁ interval.
                Z₁ = Z₁_intervals[Z₁_inclusions][1]
        
                # Remove it from the list.
                deleteat!(Z₁_intervals, Z₁_inclusions)
        
                # Add new intervals.
                push!(Z₁_intervals, Interval(inf(Z₁), inf(A))) 
                push!(Z₀_intervals, A)
                push!(Z₁_intervals, Interval(sup(A), sup(Z₁))) 
        
                push!(Aₙ, compute_interval_image(A, parts, vector_fields))
            end
        end
        
        filter!(x-> diam(x) > 0.0, Z₁_intervals)

    end
    return Z₀_intervals
end



function compute_Z_intervals(map_parameters)
    """
    Compute all Z intervals.
    """
    vector_fields = map_parameters["vector_fields"]
    parts = map_parameters["partitions"][:, 1]

    partition_intervals = map(x -> Interval(x[1], x[2]), eachrow(map_parameters["partitions"]))
    partition_images = map(x -> compute_interval_image(x, parts, vector_fields), partition_intervals)
    partition_images = map(x -> intervalUnion(inf(x), sup(x)), partition_images)
    full_interval = intervalUnion(0//1, 1//1)

    Z₀_intervals = Array{Interval{Rational{Int64}}}(undef, 0)
    Z₁_intervals = Array{Interval{Rational{Int64}}}(undef, 0)
    Z₂_intervals = Array{Interval{Rational{Int64}}}(undef, 0)

    for i in (full_interval \ compute_sequential_union(partition_images)).v
        if diam(i) > 0
            push!(Z₀_intervals, i)
        end
    end

    for (i, inter_1) in enumerate(partition_images)
        for (j, inter_2) in enumerate(partition_images)

            intersection = inter_1 ∩ inter_2
            
            if  i != j && length(intersection.v) > 0
                push!(Z₂_intervals, intersection[1])
            end
        end
    end

    Z₂_intervals = unique(Z₂_intervals)

    if !isempty(Z₀_intervals) && !isempty(Z₂_intervals)
        Z₁_union = full_interval \ (compute_sequential_union(Z₀_intervals) ∪ compute_sequential_union(Z₂_intervals))
        for Z₁ in Z₁_union.v
            push!(Z₁_intervals, Z₁)
        end
    elseif !isempty(Z₀_intervals) && isempty(Z₂_intervals)
        Z₁_union = full_interval \ (compute_sequential_union(Z₀_intervals))
        for Z₁ in Z₁_union.v
            push!(Z₁_intervals, Z₁)
        end
    elseif isempty(Z₀_intervals) && !isempty(Z₂_intervals)
        Z₁_union = full_interval \ (compute_sequential_union(Z₂_intervals))
        for Z₁ in Z₁_union.v
            push!(Z₁_intervals, Z₁)
        end
    elseif !isempty(map_parameters["c_points"])
        push!(Z₁_intervals, Interval(0//1, map_parameters["c_points"][1]))
        push!(Z₁_intervals, Interval(map_parameters["c_points"][1], 1//1))
    else
        push!(Z₁_intervals, full_interval.v[1])
    end



    
    Z_intervals = Dict(
        "Z₀" => Z₀_intervals,
        "Z₁" => Z₁_intervals,
        "Z₂" => Z₂_intervals,

    )

    no_intervals = sum(map(x -> length(Z_intervals[x]), collect(keys(Z_intervals))))
    Z_sups = []
    Z_labels = []

    for inter in Z_intervals["Z₀"]
        push!(Z_sups, sup(inter))
        push!(Z_labels, 0)
    end

    for inter in Z_intervals["Z₁"]
        push!(Z_sups, sup(inter))
        push!(Z_labels, 1)
    end


    for inter in Z_intervals["Z₂"]
        push!(Z_sups, sup(inter))
        push!(Z_labels, 2)    
    end


    return Z_intervals, convert_Zₖ_sequence_to_string(Z_labels[sortperm(Z_sups)])
end


function convert_Zₖ_sequence_to_string(Zₖ_sequence)

    Zₖ_string = ""
    for k in Zₖ_sequence
        Zₖ_string = Zₖ_string*"Z"*subscriptnumber(k)*"-"
    end

    return String(collect(Zₖ_string)[1:end-1])
end

function compute_Z₀_gaps(Z₀_intervals, map_partitions)
    
    Z₀_gaps_total = Vector{Interval{Rational{Int64}}}()

    for map_partition in map_partitions
        Z₀_intersections = compute_Z₀_intersections(Z₀_intervals, map_partition)

        if !isempty(Z₀_intersections)
            sort!(Z₀_intersections)
            
            Z₀_gaps = map(i -> Interval(sup(Z₀_intersections[i-1]), inf(Z₀_intersections[i])), 2:length(Z₀_intersections))
            pushfirst!(Z₀_gaps, Interval(inf(map_partition), inf(Z₀_intersections[1])))
            push!(Z₀_gaps, Interval(sup(Z₀_intersections[end]), sup(map_partition)))


        else
            Z₀_gaps = map_partition
        end

        for gap in Z₀_gaps
            push!(Z₀_gaps_total, gap)
        end

    end

    #### CLEAN UP
    """
    Z₀_gaps_total = unique(vcat(Z₀_gaps_total...))
    Z₀_gaps = Vector{Interval{Rational{Int64}}}()

    while length(Z₀_gaps_total) > 1

        interval_1 = popfirst!(Z₀_gaps_total)
        interval_2 = popfirst!(Z₀_gaps_total)

        if interval_1 ∩ interval_2 == ∅
            push!(Z₀_gaps, interval_1)
            pushfirst!(Z₀_gaps_total, interval_2)
        else
            pushfirst!(Z₀_gaps_total, (interval_1 ∪ interval_2)[1])
        end
    end
    push!(Z₀_gaps, Z₀_gaps_total[1])
    """
    Z₀_gaps = unique(vcat(Z₀_gaps_total...))
    return Z₀_gaps
end

function compute_upper_bound(Z₀_intervals, map_partition, g_cd)
    
    Z₀_intersections = compute_Z₀_intersections(Z₀_intervals, map_partition)

    if !isempty(Z₀_intersections)
        sort!(Z₀_intersections)
        
        Z₀_gaps = map(i -> inf(Z₀_intersections[i]) - sup(Z₀_intersections[i-1]), 2:length(Z₀_intersections))
        pushfirst!(Z₀_gaps, inf(Z₀_intersections[1]) - inf(map_partition))
        push!(Z₀_gaps, sup(map_partition) - sup(Z₀_intersections[end]) )
        
        return ceil.(Z₀_gaps./g_cd) 

    else 
        return ceil.(diam(map_partition)./g_cd) 
    end
end

function compute_upper_bounds(Z₀_union, map_partitions, g_cd)
    return map(x-> compute_upper_bound(Z₀_union, x, g_cd), map_partitions)
end


function plot_intervals(fig, ax, intervals, colour, alpha)
    
    for interval in intervals
        band!(ax, LinRange(inf(interval), sup(interval), 10), 0, 1, color=(colour, alpha))
    end
    
    return fig, ax
end



function compute_Z₀_intersections(Z₀_intervals, partition)

    intersections = filter!(x -> x != ∅, map(x -> x ∩ partition, Z₀_intervals))
    filter!(x -> diam(x) > 0, intersections)
    return intersections
end

function convert_to_interval(vector)

    return Interval(vector[1], vector[2])
end

function compute_interval_image(interval, parts, vector_fields)
    s = searchsortedlast(parts, mid(interval))
    return Interval(interval + vector_fields[s])
end


function linear_diophantine(a::Int, b::Int, c::Int)
    # Compute the gcd of a and b using the Euclidean algorithm
    d, x, y = gcdx(a, b)
    
    if c % d != 0
        return [0, 0]
    else
        # Compute the solution using Bezout's identity
        x *= c ÷ d
        y *= c ÷ d
        
        # Return the solution
        return [x, y]
    end
end


function check_solution(x, y, x_ub, y_ub, a, b, g, m)

    x_sol = x+m*(b/g) 
    y_sol = y-m*(a/g)

    solution = false

    if 0.0 <= x_sol <= x_ub && 0.0 <= y_sol <= y_ub
        solution = true
    end

    return solution, x_sol, y_sol
end



function compute_Z₀_gaps_count(Z₀_intervals, map_partitions)
    
    Z₀_gaps_total = Vector{Interval{Rational{Int64}}}()

    for map_partition in map_partitions
        Z₀_intersections = compute_Z₀_intersections(Z₀_intervals, map_partition)

        if !isempty(Z₀_intersections)
            sort!(Z₀_intersections)
            
            Z₀_gaps = map(i -> Interval(sup(Z₀_intersections[i-1]), inf(Z₀_intersections[i])), 2:length(Z₀_intersections))
            pushfirst!(Z₀_gaps, Interval(inf(map_partition), inf(Z₀_intersections[1])))
            push!(Z₀_gaps, Interval(sup(Z₀_intersections[end]), sup(map_partition)))


        else
            Z₀_gaps = map_partition
        end

        for gap in Z₀_gaps
            push!(Z₀_gaps_total, gap)
        end

    end

    Z₀_gaps_total = unique(vcat(Z₀_gaps_total...))
    Z₀_gaps = Vector{Interval{Rational{Int64}}}()

    while length(Z₀_gaps_total) > 1

        interval_1 = popfirst!(Z₀_gaps_total)
        interval_2 = popfirst!(Z₀_gaps_total)

        if interval_1 ∩ interval_2 == ∅
            push!(Z₀_gaps, interval_1)
            pushfirst!(Z₀_gaps_total, interval_2)
        else
            pushfirst!(Z₀_gaps_total, (interval_1 ∪ interval_2)[1])
        end
    end
    push!(Z₀_gaps, Z₀_gaps_total[1])
   
    return Z₀_gaps
end


function compute_actual_upper_bound(map_parameters, upper_bounds)

    v = map_parameters["vector_fields"]
    v = v.*lcm(denominator.(v))
    v = numerator.(v)
    i_star = findmin(abs.(upper_bounds))[2]
    i_other = filter(x -> x != i_star, 1:3)

    v[i_star] = v[i_star]*-1

    v_1 = v[i_other[1]]
    v_2 = v[i_other[2]]
    v_3 = v[i_star]
    g = gcd(v_1, v_2)
    ub_1 = upper_bounds[i_other[1]]
    ub_2 = upper_bounds[i_other[2]]

    solutions = []

    for c_mult in 0:upper_bounds[i_star]

        base_solution = linear_diophantine(v_1, v_2, c_mult*v_3)

        base_solution_1_bound = abs((ub_1 - base_solution[1])*(g/v_2))
        base_solution_2_bound = abs((ub_2 - base_solution[2])*(g/v_1))

        m_bound = ceil(maximum([base_solution_1_bound, base_solution_2_bound]))
        
        possible_solutions = check_solution.(
                                                base_solution[1], 
                                                base_solution[2],
                                                ub_1, 
                                                ub_2, 
                                                v_1, 
                                                v_2, 
                                                g,
                                                -m_bound:m_bound
        )

        filter!(x -> x[1], possible_solutions)

        for sol in possible_solutions
            solution = zeros(1, 4)
            solution[i_other[1]] = sol[2]
            solution[i_other[2]] = sol[3]
            solution[i_star] = c_mult
            solution[4] = sum(solution[1:3])
            push!(solutions, solution)
        end
    end

    return findmax(x -> x[4], solutions)[1]
end


function diophantine_equation(coeffs::Vector{Int}, d::Int, upper_bounds::Vector{Int})
    n = length(coeffs)
    solutions = []
    for x in Iterators.product((0:upper_bounds[i] for i = 1:n)...)
        if dot(x, coeffs) == d && all(x .>= 0) && sum(x) > 0
            push!(solutions, x)
        end
    end
    return solutions
end


function compute_all_solutions(map_parameters, upper_bounds)
    v = map_parameters["vector_fields"]
    v = numerator.(v.*lcm(denominator.(v)))
    solutions = diophantine_equation(v, 0, upper_bounds)
    return solutions
end
function compute_actual_bound(map_parameters, upper_bounds)
    v = map_parameters["vector_fields"]
    v = numerator.(v.*lcm(denominator.(v)))
    solutions = diophantine_equation(v, 0, upper_bounds)
    return findmax(x-> sum(x), solutions)[1]
end


function simulate_period_upper_bound(map_parameters, no_samples)

    temp_periods = zeros(Int64, no_samples)
    x = Vector{Rational{Int64}}(undef, map_parameters["no_attractor_samples"])
    S = Vector{Int64}(undef,  map_parameters["no_attractor_samples"])

    @inbounds for k in 1:no_samples
        # Run the map.
        run_map!(x, S, map_parameters)

        temp_periods[k] = compute_dynamics(x, S, map_parameters)["period"]
    end
    return maximum(temp_periods)
end

function compute_period_upper_bound(map_parameters)
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

function compute_dynamics(x, S, map_parameters)
    """
    This function determines the asymptotic behaviour of the map and 
    characterises whether the attractor is cyclic etc.
    """
    # Determine the cycle count.
    period = length(keys(counter(x)))

    # Determine the symbolic itinerary.
    branch_symbols = determine_branch_symbols(map_parameters["vector_fields"])

    Z_intervals, Zₖ_sequence = compute_Z_intervals(map_parameters)

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



function cobweb_plot!(ax, x, no_points, colour, line_style)

    for i in range(1, no_points)
        x_p = x[end-(no_points+3) + i]
        y_p = x[end-(no_points+3) + i+1]
        f_y_p = x[end-(no_points+3) + i+2]

        CairoMakie.lines!(ax, [x_p, y_p], [y_p, y_p],  color=colour, alpha=0.01, linestyle=line_style)
        CairoMakie.lines!(ax, [y_p, y_p], [y_p, f_y_p], color=colour, alpha=0.01, linestyle=line_style)
      
        
    end
    
end



function add_bifurcation_lines(fig, ax, τ_min, τ_max, ϕ_min, ϕ_max, τ_range, ϕ_range)

    τ_range = LinRange(1, 2, 10000)
    ϕ_range = LinRange(1, 2, 10000)
    
    τ_points_b, ϕ_points_b = compute_bisections(τ_range, ϕ_range)
    τ_points_z, ϕ_points_z = compute_zeros(τ_range, ϕ_range)
        
    
    for key in keys(τ_points_z) 
        if !isempty(τ_points_z[key])

            τ_points = τ_points_z[key]
            ϕ_points = ϕ_points_z[key]
                    
            τ_indices = τ_min .< τ_points .< τ_max
            ϕ_indices = ϕ_min .< ϕ_points .< ϕ_max

            if !isempty(τ_indices)
                CairoMakie.scatter!(
                    ax,
                    τ_points[τ_indices .+ ϕ_indices .== 2],
                    ϕ_points[τ_indices .+ ϕ_indices .== 2],
                    color=:black,
                    markersize=2,
                    transparency=false,
                )
            end
        end
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





function compute_periods_batch!(τ, periods, numbers)

    @inbounds for (index, ϕ) in enumerate(numbers)

        map_parameters = compute_map(τ, ϕ)   
        map_parameters["no_samples"] = Int(1e5)
        map_parameters["no_attractor_samples"] = Int(1e4)    
        map_parameters["x₀"] = rationalize(rand(1)[1])

        periods[index] = compute_period_upper_bound(map_parameters)
    end
end


function simulate_periods_batch!(τ, periods, numbers)


    @Threads.threads for i in range(1, length(numbers))


        map_parameters = compute_map(τ, numbers[i])   
        map_parameters["no_samples"] = Int(1e5)
        map_parameters["no_attractor_samples"] = Int(1e4)
        map_parameters["x₀"] = rationalize(rand(1)[1])

        periods[i] = simulate_period_upper_bound(map_parameters, 25)

    end
end
