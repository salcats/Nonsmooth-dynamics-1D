include("new_utilities.jl")


A(x::Float64, constant::Float64, coeff::Float64) = constant + coeff*abs(x)^2
B(x::Float64, constant::Float64, coeff::Float64) = constant + coeff*abs(x)^2



min_val = 1
max_val = 5
max_den = 200

numbers = generate_rationals(min_val, max_val, max_den)
τ_values, ϕ_values = find_n_branch_maps(numbers[1:1000], 3)


τ = 5//1



index = sample(collect(1:length(τ_values)))

ϕ = ϕ_values[index]

map_parameters = compute_map(τ, ϕ)   
map_parameters["no_samples"] = Int(1e5)
map_parameters["no_attractor_samples"] = Int(1e4)
map_parameters["x₀"] = sample(τ_values)



simulation_parameters = Dict(
                                "A_constant" => 1/numerator(ϕ),
                                "B_constant" => 1/denominator(ϕ),
                                "no_samples" => 1000000,
                                "τ" => Int(τ),
                                "initial_condition" => map_parameters["x₀"] ,
                                "A" => A,
                                "B" => B
)

A_coeff = rand(1)[1]
B_coeff = rand(1)[1]
perturbation_magnitude = 1e-8
#(1/100000)

simulation_parameters["A_coeff"] = perturbation_magnitude*A_coeff
simulation_parameters["B_coeff"] = perturbation_magnitude*B_coeff

fig, ax = plot_difference(simulation_parameters, map_parameters, perturbation_magnitude)
fig

i = 1

save("plots/"*string(i)*".png")



function plot_difference(simulation_parameters, map_parameters, perturbation_magnitude)

    x_nonconstant, S_nonconstant = run_nonconstant_sim(simulation_parameters)
    x_nonconstant = x_nonconstant[findall(x-> x .> 0.0, diff(S_nonconstant))]
    x_nonconstant = x_nonconstant./simulation_parameters["A_constant"] .- Int(τ)

    simulation_parameters["A_coeff"] = 0.0
    simulation_parameters["B_coeff"] = 0.0

    x_constant, S_constant = run_nonconstant_sim(simulation_parameters)
    x_constant = x_constant[findall(x-> x .> 0.0, diff(S_constant))]
    x_constant = x_constant./simulation_parameters["A_constant"] .- Int(τ)

    fig, ax = plot_map(map_parameters)

    cobweb_plot!(ax, x_constant, 100, :green, :solid)
    cobweb_plot!(ax, x_nonconstant, 100, :red, :dot)

    ax.title = "x^2 coeff: "*string(perturbation_magnitude)
    return fig, ax
end










function run_nonconstant_sim(simulation_parameters)
    
    no_samples = simulation_parameters["no_samples"]
    τ = simulation_parameters["τ"]
    initial_condition = simulation_parameters["initial_condition"]
    A_coeff = simulation_parameters["A_coeff"]
    B_coeff = simulation_parameters["B_coeff"]
    A_constant = simulation_parameters["A_constant"]
    B_constant = simulation_parameters["B_constant"]
    A = simulation_parameters["A"]
    B = simulation_parameters["B"]

    x = zeros(no_samples)
    vec_mode = zeros(no_samples)

    x[1:τ] .= initial_condition
    vec_mode[1:τ] .= 1*(x[1]<0) + 2*(x[1]>0)

    run_inner_loop!(x, vec_mode, τ, no_samples, A_coeff, B_coeff, A_constant, B_constant, A , B)
    
    return x, Int.(vec_mode)

end

function run_inner_loop!(
                            x::Vector{Float64}, 
                            vec_mode::Vector{Float64}, 
                            τ::Int64, 
                            no_samples::Int64, 
                            A_coeff::Float64, 
                            B_coeff::Float64, 
                            A_constant::Float64,
                            B_constant::Float64, 
                            A::Function,
                            B::Function,
    )

    @inbounds for i in range(τ, no_samples-1)
        x[i+1] = x[i] + (x[i-τ] < 0)*A(x[i-τ], A_constant, A_coeff) - (x[i-τ]>0)*B(x[i-τ], B_constant, B_coeff)
        vec_mode[i+1] = 1*(x[i-τ] < 0 ) + 2*(x[i-τ] > 0)
    end

    return x, vec_mode

end






function plot_nonconstant(x, vec_mode, no_points_cobweb, time_series_segment_length)


    pos_turning_point_indices = findall(x-> x > 0.0, diff(vec_mode))
    neg_turning_point_indices = findall(x-> x < 0.0, diff(vec_mode))


    fig = CairoMakie.Figure()
    fontsize_theme = CairoMakie.Theme(fontsize = 20)
    CairoMakie.set_theme!(fontsize_theme)

    ax_1 = CairoMakie.Axis(
        fig[1, 1],
        title="Left TPR"
    )

    ax_2 = CairoMakie.Axis(
        fig[1, 2],
        title="Right TPR",
    )


    ax_3 = CairoMakie.Axis(
        fig[2, 1],
        title="Time-series segment",
    )


    ax_4 = CairoMakie.Axis(
        fig[2, 2],
        title="Cobweb",
        xlabel="xₙ",
        ylabel="xₙ₊₁",
    )

    CairoMakie.scatter!(
        ax_1,
        range(1, length(x[neg_turning_point_indices][1000:end])),
        x[neg_turning_point_indices][1000:end],
        markersize=1
    )

    CairoMakie.scatter!(
        ax_2,
        range(1, length(x[pos_turning_point_indices][1000:end])),
        x[pos_turning_point_indices][1000:end],
        markersize=1
    )    #
    CairoMakie.lines!(ax_3, x[end-time_series_segment_length:end], markersize=0.1)

    diagonal_points = CairoMakie.LinRange(minimum(x[end-no_points_cobweb:end]), maximum(x[end-no_points_cobweb:end]), 1000)
    CairoMakie.lines!(ax_4, diagonal_points, diagonal_points, color=:gray, alpha=0.2)

    colour=:red
    for i in range(1, no_points_cobweb)
        x_p = x[end-(no_points_cobweb+3) + i]
        y_p = x[end-(no_points_cobweb+3) + i+1]
        f_y_p = x[end-(no_points_cobweb+3) + i+2]

        CairoMakie.lines!(ax_4, [x_p, y_p], [y_p, y_p],  color=colour, alpha=0.01)
        CairoMakie.lines!(ax_4, [y_p, y_p], [y_p, f_y_p], color=colour, alpha=0.01)
    end

end

