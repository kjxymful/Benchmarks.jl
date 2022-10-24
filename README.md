# Benchmarks

Create the most commonly used Dynamical Systems Benchmark systems

## Getting started
Clone the directory, go to directory, open the package manager in the Julia REPL using `"]"`, and add the package:
```
(@v1.7) pkg> add .
```

Or enter the package manager in the Juila REPL, and add the package with:
```
(@v1.7) pkg> add github_link
```

## Generate Data
from directory do
```
julia --project generate_benchmarks.jl
```
with arguments specified in settings/benchmark_defaults.json, i.e. you need a folder settings with benchmark_defaults.json in the parent directory

or in commandline with
```
julia --project generate_benchmarks.jl --kwarg values
```
choose the benchmark system with --name 

If you want to create your own system, use the Dynamical Systems in either src/ds_models.jl, src/ns_models.jl with your parameters, and generate a series with generate_trajectories in src/utils.jl

## Implemented Systems
### Standard Benchmarks
- standard_bursting: Bursting Neuron as used in all Durstewitz Lab papers
- standard_lorenz: Lorenz, with ρ = 28
- lorenz_limit_cycle: A limit cycle of the Lorenz System with ρ=24
- bursting_limit_cycle: The Bursting neuron limit cycle at gₙₘ₀ₐ=10

### Regimes
- bursting_neuron_regimes

### Non-Stationary Benchmarks
- ExplodingLorenz : Starts with a limit cycle and ends in the 
well known chaotic attractor; ρ=22->28
- ShiftingLorenz : Starts with the chaotic attractor and shifts it "forward"; ρ=28->22
- ShrinkingLorenz : Starts with the chaotic attractor and shrinks and shifts it a bit; ρ=28->23, σ=10->5, β=8/3->0.5
- PaperLorenzBigChange : The ns system used in Patel et al. 2022 with a quick parameter change
- PaperLorenzSmallChange : The ns system used in Patel et al. 2022 with a slow parameter change

### Trial Benchmark Systems
- trial_lorenz : a lorenz with parameter shifting from 22->28 across trials
- split_lorenz : a specified ns lorenz split into num_trials

see docs for more information on what kwargs are availabe