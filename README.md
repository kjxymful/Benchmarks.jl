# Benchmarks

Create the most commonly used Dynamical Systems Benchmark systems

## Getting started
Clone the directory, go to directory, open the package manager in the Julia REPL using `"]"`, and add the package:
```
(@v1.8) pkg> add .
```

Or enter the package manager in the Juila REPL, and add the package with:
```
(@v1.8) pkg> add github_link
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

#### Usage in scripts
```
using Benchmarks

generate_benchmarks(parse_commandline())
```
that's it, adjust settings in settings/benchmark_defaults.json

If you want to create your own system follow the steps in template.jl

*Little Warning, process noise increases computation time drastically*

## Implemented Systems
### Standard Benchmarks
- standard_bursting: Bursting Neuron as used in all Durstewitz Lab papers
- standard_lorenz: Lorenz, with ρ = 28
- lorenz_limit_cycle: A limit cycle of the Lorenz System with ρ=24
- bursting_limit_cycle: The Bursting neuron limit cycle at gₙₘ₀ₐ=10

### Regimes
- bursting_neuron_regimes

### Non-Stationary Benchmarks
#### Lorenz
- ShrinkingLorenz : Starts with the chaotic attractor and shrinks plus shifts it a bit; ρ=28->23, σ=10->5, β=8/3->0.5 (transient Time should not be bigger than 50, as it leads to weird regimes)

- PaperLorenzBigChange : The ns system used in Patel et al. 2022 with a quick parameter change

- PaperLorenzSmallChange : The ns system used in Patel et al. 2022 with a slow parameter change

- ExplodingLorenz : Starts with a limit cycle and ends in the 
well known chaotic attractor; ρ=22->28 (not a very stable system, i.e. sensitive to initial conditions, makes it harder to train on)

- ShiftingLorenz : Starts with the chaotic attractor and shifts it "forward"; ρ=28->22 (limit cycle does not appear in time series, -> not great for snapshot comparison)


#### Bursting Neuron
The timescale is much longer than for the Lorenz, and thus needs a lot more time points

- RampUpBN : Starts with the a cycle, and adds bursting loops; g=2->4

- StopBurstBN : Starts with the complicated cycle at g=9.25 and ends in a limit cycle at g=10.15

### Trial Benchmark Systems
- trial_lorenz : a lorenz with parameter shifting from 22->28 across trials
- split_lorenz : a specified ns lorenz split into num_trials

see docs for more information on what kwargs are availabe
