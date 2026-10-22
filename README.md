# Nonsmooth dynamical systems — one-dimensional rational case

Julia code for studying one-dimensional switched (piecewise-constant) dynamical
systems with rational parameters, from my PhD research in applied dynamical
systems at the University of Bristol.


## What's here

- `utilities.jl`, `new_utilities.jl` — the core: rational number generation,
  branch maps and their topology, discontinuity points, interval-arithmetic
  bounds, period detection, and Markov chain analysis of the mode dynamics.
- `basic.jl` — generates rational vector fields and finds their branch maps.
- `single_example.jl` — a worked single case with its transition matrix.
- `sampling_plots.jl` — figures over sampled parameter ranges.
- `nonconstant_vec_fields.jl` — extension to non-constant fields.

Note: `utilities.jl` and `new_utilities.jl` are successive versions of the core
routines, with some overlap. `new_utilities.jl` is the later rewrite, adding
interval-arithmetic bounds and period detection. Both are kept because different
scripts depend on each.

## Running it

```julia
using Pkg; Pkg.activate("."); Pkg.instantiate()
include("basic.jl")
```

## Built with

IntervalArithmetic and IntervalRootFinding for rigorous numerics, JuMP with
HiGHS for the optimisation subproblems, ModelingToolkit and Polynomials for
symbolic work, DiscreteMarkovChains and DynamicalSystems for the analysis, and
CairoMakie for figures.

## Related papers

> Catsis, S., Hall, C.L., Jeffrey, M.R. *The hidden sensitivity of non-smooth
> dynamics.* Physica D 463: 134165, 2024.

> Catsis, S., Hall, C.L., Jeffrey, M.R. *Perturbation and determinacy of
> nonsmooth systems.* Physica D 477: 134679, 2025.
