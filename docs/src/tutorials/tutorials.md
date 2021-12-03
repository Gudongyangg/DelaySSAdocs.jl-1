# Tutorials

This tutorial is an introduction to using DelaySSA to define chemical reaction models, define the problem that we call `DelayJumpProblem` and finally solve and visualize the result. To demonstrate this functionality, we will consider a specific case of a chemical system as follows
```math
S+I\xrightarrow{\rho}E+I\\
I\stackrel{r}{\rightarrow}R
```
<!-- $$
S+I\xrightarrow{\rho}E+I\\
I\stackrel{r}{\rightarrow}R
$$ -->
and $S+I\xrightarrow{\rho} R$ will trigger $E\Rightarrow I$ after $\tau$ time.

## Model Initialisasion

[Catalyst.jl](https://github.com/SciML/Catalyst.jl) provides a comprehensive interface to modelling chemical reaction networks in Julia and can be used to construct models fully-compatible with DelaySSA. For more details on how to do so we recommend reading [Catalyst's tutorial](https://catalyst.sciml.ai/stable/tutorials/using_catalyst/). This way, the model can be defined as:

```julia
rn = @reaction_network begin
    ρ, S+I --> E+I
    r, I --> R
end ρ r
```

## Define `DelayJumpProblem`

We have two routes to define our `DelayJumpProblem`, one way is based on `DiscreteProblem` and `JumpSystem`  , the other is based on `DiscreteProblem` and `JumpSet`

### `JumpSystem + DiscreteProblem `

We can define `DelayJumpProblem` through ` JumpSystem`  and `DiscreteProblem` , then we first define `Jumpsystem`

```julia
jumpsys = convert(JumpSystem, rn, combinatoric_ratelaws=false)
```

The use of function `convert` is to convert the `rn` to a `JumpSystem`

Then we turn to `DiscreteProblem`. We should first set

```julia
u0 = [999,1,0,0]
de_chan0 = [[]]
tf = 400.
tspan = (0,tf)
ps = [1e-4, 1e-2]
τ = 20.
```

So we can define the `DiscreteProblem` by

```julia
dprob = DiscreteProblem(jumpsys,u0,tspan,ps)
```

where `DiscreteProblem` inputs `jumpsys`, and the initial condition of molecular numbers `u0` , the timespan `tspan` and the rate of two reactions `ps`

Since there is a delay reaction, so we must define 

```julia
delay_trigger_affect! = function (de_chan, rng)
    append!(de_chan[1], τ)
end
delay_trigger = Dict(1=>delay_trigger_affect!)
delay_complete = Dict(1=>[2=>1, 3=>-1])
delay_interrupt = Dict()
delayjumpset = DelayJumpSet(delay_trigger, delay_complete, delay_interrupt)
```

- `delay_trigger`  
  - Keys: Indices of reactions defined in `jumpset` that can trigger the delay reaction. Here we have the  reaction $E\stackrel{r}\rightarrow I$ that will trigger the $E$ to degrade and $I$ to appear after time $\tau$.
  
  - Values: A update function that determines how to update the delay channel. In this example, once the delay reaction is trigged, the delay channel will be added a delay time $\tau$
  
- `delay_interrupt`
  - There are no delay interrupt reactions in this example so we set `delay_interrupt = Dict()`
- ```delay_complete:``` 
  - Keys: Indices of delay channel.
  - Values: A vector of `Pair`s, mapping species id to net change of stoichiometric coefficient.

We can see more details in [Defining a `DelayJumpSet`(bursty)](https://palmtree2013.github.io/DelaySSAToolkit.jl/dev/tutorials/bursty/#Defining-a-DelayJumpSet) and [Defining a `DelayJumpSet`(birth-death example)](https://palmtree2013.github.io/DelaySSAToolkit.jl/dev/tutorials/delay_degradation/#Defining-a-DelayJumpSet) or [Defining a `DelayJumpSet`(multi-next-delay example)](https://palmtree2013.github.io/DelaySSAToolkit.jl/dev/tutorials/delay_multidegradation/#Defining-a-DelayJumpSet)

At last, we can define the `DelayJumpProblem` by 

```julia
jprob = DelayJumpProblem(jumpsys, dprob, DelayRejection(), delayjumpset, de_chan0, save_positions=(true,true))
```

where `DelayJumpProblem` inputs `jumpsys`,`DiscreteProblem`, `DelayJumpSet`, the algorithm we choose and the initial condition of the delay channel `de_chan0`.

So we can solve the problem and visualize the result

```julia
sol = solve(jprob, SSAStepper(), seed = 1234)
using Plots
fig = plot(sol, label = ["S" "I" "E" "R"], linewidth = 3, legend = :top, ylabel = "# of individuals", xlabel = "Time", fmt=:svg)
```

### `JumpSet + DiscreteProblem `

We first define `Jumpset` before that we should first define the parameters and the mass-action jump (see [Defining a Mass Action Jump](https://diffeq.sciml.ai/stable/types/jump_types/#Defining-a-Mass-Action-Jump) for details).

```julia 
ρ, r = [1e-4, 1e-2]
rate1 = [ρ, r]
reactant_stoich = [[1=>1,2=>1],[2=>1]]
net_stoich = [[1=>-1,3=>1],[2=>-1,4=>1]]
mass_jump = MassActionJump(rate1, reactant_stoich, net_stoich; scale_rates =false)
jumpset = JumpSet((),(),nothing,[mass_jump])
```

We can also see more details of the parameters in [this example](https://palmtree2013.github.io/DelaySSAToolkit.jl/dev/tutorials/delay_multidegradation/tutorials/@ref Model_definition).

Then we initialise the problem by setting

```julia
u0 = [999,1,0,0]
de_chan0 = [[]]
tf = 400.
tspan = (0,tf)
τ = 20.
```

So we can define the `DiscreteProblem`

```Julia
dprob = DiscreteProblem(u0, tspan)
```

The same as we did before, we must define the  `DelayJumpSet`

```julia
delay_trigger_affect! = function (de_chan, rng)
    append!(de_chan[1], τ)
end
delay_trigger = Dict(1=>delay_trigger_affect!)
delay_complete = Dict(1=>[2=>1, 3=>-1])
delay_interrupt = Dict()
delayjumpset = DelayJumpSet(delay_trigger, delay_complete, delay_interrupt)
```

So we can define the problem

```julia 
djprob = DelayJumpProblem(dprob, DelayRejection(), jumpset, delayjumpset, de_chan0, save_positions=(true,true))
```

where `DelayJumpProblem` inputs `DiscreteProblem`, `JumpSet`,`DelayJumpSet`, the algorithm we choose  and the initial condition of the delay channel `de_chan0`.

At last, we can solve the problem and visualize it.

```julia
sol = solve(djprob, SSAStepper(), seed = 1234)
using Plots
fig = plot(sol, label = ["S" "I" "E" "R"], linewidth = 3, legend = :top, ylabel = "# of individuals", xlabel = "Time", fmt=:svg)
```
![seir](../assets/seir.svg)
