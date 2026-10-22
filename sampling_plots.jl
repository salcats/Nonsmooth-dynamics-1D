include("new_utilities.jl")
using StatsBase
using LinearAlgebra
using CairoMakie
using IntervalRootFinding


τ = 3//1

min_val = 1
max_val = 5
max_den = 50
numbers_1 = generate_rationals(min_val, max_val, max_den)

min_val = 1
max_val = 5
max_den = 100
numbers_2 = generate_rationals(min_val, max_val, max_den)

numbers_3 = sort!(sample(numbers_2, 1000))
numbers_4 = rationalize.(LinRange(1, 5, 100))
numbers_5 = rationalize.(LinRange(1, 5, length(numbers_2)))



periods_1 = zeros(length(numbers_1))
periods_2 = zeros(length(numbers_2))
periods_3 = zeros(length(numbers_3))
periods_4 = zeros(length(numbers_4))
periods_5 = zeros(length(numbers_5))


@time simulate_periods_batch!(τ, periods_5, numbers_5)





#### Will take about 60 minutes or less to do 10000 simulations so we can just run that later on 
#### at work

compute_periods_batch!(τ, periods_1, numbers_1)
compute_periods_batch!(τ, periods_2, numbers_2)
compute_periods_batch!(τ, periods_3, numbers_3)
compute_periods_batch!(τ, periods_4, numbers_4)


#@time compute_periods_batch!(τ, periods_4, numbers_5)


fig = CairoMakie.Figure(resolution = (2000, 1000), font = "Avenir Light")
fontsize_theme = CairoMakie.Theme(fontsize = 50)
CairoMakie.set_theme!(fontsize_theme)
ax_11 = CairoMakie.Axis(
    fig[1, 1],
    xlabel="ϕ",
    ylabel="Period",
    xlabelsize=50,
    xticklabelsize=50,
    ylabelsize=50,
    yticklabelsize=50,
    yticks=([0, 50, 100], ["0", "50", "100"])
)

ax_21 = CairoMakie.Axis(
    fig[2, 1],
    xlabel="ϕ",
    ylabel="Period",
    xlabelsize=50,
    xticklabelsize=50,
    ylabelsize=50,
    yticks=([0, 1525, 3050], ["0", "1525", "3050"])
)

CairoMakie.scatter!(ax_11, numbers_4, periods_4, markersize=5)
CairoMakie.scatter!(ax_21, numbers_5, periods_5, markersize=3)
fig

save("plots/new_sample_4.png", fig)