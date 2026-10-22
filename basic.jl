include("new_utilities.jl")
using StatsBase
using LinearAlgebra
using CairoMakie
using IntervalRootFinding

min_val = 1
max_val = 5
max_den = 200

numbers = generate_rationals(min_val, max_val, max_den)
τ_values, ϕ_values = find_n_branch_maps(numbers[1:1000], 3)

index = sample(collect(1:length(τ_values)))

τ = 10//1
ϕ = ϕ_values[index]

map_parameters = compute_map(τ, ϕ)   
map_parameters["no_samples"] = Int(1e5)
map_parameters["no_attractor_samples"] = Int(1e4)
map_parameters["x₀"] = rationalize(rand(1)[1])

# Run the map.
@time x, S = run_map(map_parameters)

#dynamics_data = compute_dynamics(x, S, map_parameters)
#indices =  sample(collect(1:length(τ_values)), 10)
#compute_period_upper_bound(map_parameters)
#simulate_period_upper_bound(map_parameters, 100)



fig, ax = plot_map(map_parameters)

cobweb_plot!(ax, x, 100, :green)
cobweb_plot!(ax, Dᵣ_points, 100, :red)

fig



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
