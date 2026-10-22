# Switched dynamical systems — one-dimensional rational case

Julia code for studying one-dimensional switched (piecewise-constant) dynamical
systems with rational parameters, from my PhD research in applied dynamical
systems at the University of Bristol.

The underlying question is whether an idealised switching model remains a
reliable guide to a real system that cannot switch instantaneously.

## The idea

For a vector field with rational entries, choosing a lattice with spacing
`1 // lcm(denominators)` means one integration step maps a whole grid cell
exactly onto another grid cell. The continuous dynamics then collapse onto a
finite directed graph on cell indices, where every node has out-degree one.

That graph decomposes into cycles with trees of transients feeding into them:
the cycles are the periodic orbits — found exactly rather than approximated by
simulation — and the transients hanging off each cycle are its basin of
attraction.

Exactness is the point, which is why the code works in exact rational arithmetic
throughout rather than floating point, and uses interval arithmetic where
rigorous bounds are needed.

## What's here

- `utilities.jl`, `new_utilities.jl` — the core: rational lattice construction,
  branch maps, graph construction and analysis, Markov chain analysis of the
  resulting mode dynamics.
- `basic.jl` — generates rational vector fields and finds their branch maps.
- `single_example.jl` — a worked single case with its transition matrix.
- `sampling_plots.jl` — figures over sampled parameter ranges.
- `nonconstant_vec_fields.jl` — extension to non-constant fields.

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
