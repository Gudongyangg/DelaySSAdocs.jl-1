# [Tutorials](@id seir_model)

This tutorial is aims to explain how to use DelaySSAToolkit to define chemical reaction models, solve the problem and visualize the results. To demonstrate the functionalities, we will consider a specific case from epideimc modelling as follows
```math
S+I\xrightarrow{\rho}E+I\\
I\stackrel{r}{\rightarrow}R
```
and $S+I\xrightarrow{\rho} E+I$ will trigger $E\Rightarrow I$ after $\tau$ time, where $S$, $I$ and $R$ are the susceptible, infected and removed populations. $E$ represents the exposed population. It means, wtih rate $\rho$, a susceptible contacted by an infected will become an individual that is exposed to the disease and then it takes certain amount of time delay $\tau$ to become an infected inidividual.

# Model
What differs from the Markov process that can be modelled via SSA is the introduction of **delay reactions**. To show how we incorporate the delay reactions into the Markovian system, we first need to define the Markovian part and then its non-Markovian part. These two parts mainly form a `DelayJumpProblem`. Here we show two routes to define our delay system, one way is based on `JumpSystem`, `DiscreteProblem` and `DelayJumpSet`, the other is based on `JumpSet`, `DiscreteProblem` and `DelayJumpSet`.

## First route: `JumpSystem + DiscreteProblem + DelayJumpSet`
### [Markovian part](@id Markovian_part)
[Catalyst.jl](https://github.com/SciML/Catalyst.jl) provides a comprehensive interface to modelling reaction networks in Julia and can be used to construct models fully-compatible with DelaySSAToolkit. For more details on how to construct a reaction network, we recommend reading [Catalyst's tutorial](https://catalyst.sciml.ai/stable/tutorials/using_catalyst/). In our example, the Markovian part (**model without delays**) of the model can be defined as:
```julia
rn = @reaction_network begin
    ρ, S+I --> E+I
    r, I --> R
end ρ r
```
We can easily obtain a `Jumpsystem` from the reaction network `rn` that is previously defined using Catalyst.

```julia
jumpsys = convert(JumpSystem, rn, combinatoric_ratelaws=false)
```
where `combinatoric_ratelaws` is an optional parameter that specifies whether the rate constants correspond to stochastic rate constants in the sense used by Gillespie, and hence need to be rescaled. The default, `combinatoric_ratelaws=true`, corresponds to rescaling the passed in rate constants. The default behavior is to assume rate constants correspond to stochastic rate constants in the sense used by Gillespie [1]. This means that for a reaction such as $2A \overset{k}{\rightarrow} B$, the jump rate function constructed by `MassActionJump` would be `k*A*(A-1)/2!`. For a trimolecular reaction like $3A \overset{k}{\rightarrow} B$ the rate function would be `k*A*(A-1)*(A-2)/3!`. To *avoid* having the reaction rates rescaled (by `1/2` and `1/6` for these two examples), one can pass the `JumpSystem` constructor the optional named parameter `combinatoric_ratelaws=false` see [reaction rate laws used in simulations](https://catalyst.sciml.ai/stable/tutorials/using_catalyst/#Reaction-rate-laws-used-in-simulations) for details.

With the initial conditions, we can then define `DiscreteProblem`
```julia
u0 = [999,1,0,0] # S, I, E, R
de_chan0 = [[]]
tf = 400.
tspan = (0,tf)
ps = [1e-4, 1e-2] # parameters for ρ, r
τ = 20.
dprob = DiscreteProblem(jumpsys,u0,tspan,ps)
```
where `DiscreteProblem` inputs `jumpsys`, the initial condition of reactants `u0`, the simulation timespan `tspan` and the reaction rates `ps`.

### Non-Markovian part
The non-Markovian part consists of three elements:
- delay trigger reactions: those reactions in the [Markovian part](@ref Markovian_part) that trigger the change of the state of the delay channels or/and the state of the reactants upon initiation.
- delay interrupt reactions: those reactions in the [Markovian part](@ref Markovian_part) that change the state of the delay channels or/and the state of the reactants in the middle of on-going delay reactions.
- delay complete reactions: those reactions that are initiated by delay trigger reactions and change the state of the delay channels or/and the state of the reactants upon completion.
  
With these three definitions in mind and based on this particular example, we define the `DelayJumpSet` by
```julia
delay_trigger_affect! = function (integrator, rng)
    append!(integrator.de_chan[1], τ)
end
delay_trigger = Dict(1=>delay_trigger_affect!)
delay_interrupt = Dict()
delay_complete = Dict(1=>[2=>1, 3=>-1])
delayjumpset = DelayJumpSet(delay_trigger, delay_complete, delay_interrupt)
```
- `delay_trigger::Dict`  A dictionary that contains
  - Keys: Indices of reactions defined in the [Markovian part](@ref Markovian_part) that can trigger the delay reaction see [Remark](@ref Remark) on the reaction indices. Here we have the first reaction $S+I\Rightarrow E+ I$ that will trigger the transfer from $E$ to $I$ after time $\tau$, hence the key here is `1`.
  
  - Values: An update function that determines how to update the delay channel. In this example, once the delay reaction is triggered, the first delay channel will be added a delay time $\tau$. The update function has two inputs: 1. `integrator`: which stores the current state of the reactants (`integrator.u`) and the delay channels (`integrator.de_chan`); 2. `rng`: the random seed for a given stochastic simulation.
  
- `delay_interrupt::Dict`
  - There are no delay interrupt reactions in this example so we set `delay_interrupt = Dict()`.
- ```delay_complete::Dict``` 
  - Keys: Indices of delay channels. Here we only have one delay channel.
  - Values: A vector of `Pair`s, mapping species index to net change of stoichiometric coefficient. Here once the delay reaction is completed, the third species *E* is transfered to *I*, thus the net change of the second species is $1$ and the third $-1$. For the details about the order of the species see [Remark](@id Remark).

We refer to [Defining a `DelayJumpSet` (bursty model)](bursty.md/#Defining-a-DelayJumpSet) and [Defining a `DelayJumpSet`(birth-death model)](delay_degradation.md/#Defining-a-DelayJumpSet) for more details.

Now we can define the `DelayJumpProblem` by 
```julia
djprob = DelayJumpProblem(jumpsys, dprob, DelayRejection(), delayjumpset, de_chan0, save_positions=(true,true))
```
where `DelayJumpProblem` inputs `jumpsys`,`DiscreteProblem`, `DelayJumpSet`, the algorithm we choose and the initial condition of the delay channel `de_chan0`. The optional augment `save_positions` is a Boolean tuple for whether to save before and after the event.
Tnen one can use 
```julia
sol = solve(djprob, SSAStepper())
```
to solve the problem.

## [Second route: `JumpSet + DiscreteProblem + DelayJumpSet`](@id second_route)

Now we explain how to define the `DelayJumpProblem` in another way. To that aim, we should first define the parameters and the mass-action jump (see [Defining a Mass Action Jump](https://diffeq.sciml.ai/stable/types/jump_types/#Defining-a-Mass-Action-Jump) for details) and construct a `Jumpset`.
### Markovian part
```julia 
ρ, r = [1e-4, 1e-2]
rates = [ρ, r]
reactant_stoich = [[1=>1,2=>1],[2=>1]]
net_stoich = [[1=>-1,3=>1],[2=>-1,4=>1]]
mass_jump = MassActionJump(rates, reactant_stoich, net_stoich; scale_rates =false)
jumpset = JumpSet((),(),nothing,mass_jump)
```
We briefly explain the notations here:
- `rates` is a vector of rates of reactions.
- `reactant_stoich` is a vector whose `k`th entry is the reactant stoichiometry of the `k`th reaction. The reactant stoichiometry for an individual reaction is assumed to be represented as a vector of `Pair`s, mapping species id to stoichiometric coefficient.
- `net_stoich`  is assumed to have the same type as `reactant_stoich`; a vector whose `k`th entry is the net stoichiometry of the `k`th reaction. The net stoichiometry for an individual reaction is again represented as a vector of `Pair`s, mapping species id to the net change in the species when the reaction occurs.
- `scale_rates` is an optional parameter that specifies whether the rate constants correspond to stochastic rate constants in the sense used by Gillespie, and hence need to be rescaled. *The default, `scale_rates=true`, corresponds to rescaling the passed in rate constants.* When using `MassActionJump` the default behavior is to assume rate constants correspond to stochastic rate constants in the sense used by Gillespie (J. Comp. Phys., 1976, 22 (4)). This means that for a reaction such as $2A \overset{k}{\rightarrow} B$, the jump rate function constructed by `MassActionJump` would be `k*A*(A-1)/2!`. For a trimolecular reaction like $3A \overset{k}{\rightarrow} B$ the rate function would be `k*A*(A-1)*(A-2)/3!`. To *avoid* having the reaction rates rescaled (by `1/2` and `1/6` for these two examples), one can pass the `MassActionJump` constructor the optional named parameter `scale_rates=false`
  
The `JumpSet` consists of four inputs, namely variable jumps, constant rate jumps, regular jumps and mass-action jumps. As far as discrete stochastic simulation is concerned, we only focus on constant rate jumps and mass-action jumps which are the second and fourth slots of `JumpSet` (see [different jump types](https://diffeq.sciml.ai/stable/types/jump_types/)).
Here we only have two mass-action jumps that are wrapped in `mass_jump`.
Then we initialise the problem by setting
```julia
u0 = [999,1,0,0]
de_chan0 = [[]]
tf = 400.
tspan = (0,tf)
τ = 20.
```
As before, we can define the `DiscreteProblem`
```Julia
dprob = DiscreteProblem(u0, tspan)
```
### Non-Markovian part
In the same way, we can define the  `DelayJumpSet` by
```julia
delay_trigger_affect! = function (de_chan, rng)
    append!(de_chan[1], τ)
end
delay_trigger = Dict(1=>delay_trigger_affect!)
delay_complete = Dict(1=>[2=>1, 3=>-1])
delay_interrupt = Dict()
delayjumpset = DelayJumpSet(delay_trigger, delay_complete, delay_interrupt)
```
Now we can define the problem
```julia 
djprob = DelayJumpProblem(dprob, DelayRejection(), jumpset, delayjumpset, de_chan0, save_positions=(true,true)).
```
At last, we can solve the problem and visualize it
```julia
sol = solve(djprob, SSAStepper())
```
![seir](../assets/seir.svg)



### [Remark](@id Remark) 
!!! warning
    To check the order of the species in a reaction network `rn::ReactionSystem` (or a `jumpsys::JumpSystem`), one can call `species(rn)` (or `states(jumpsys)` respectively).
    Converting a `ReactionSystem` defined by Catalyst into a `JumpSystem` might change the order of the reactions that is in your orignal reaction network. Internally, all MassActionJumps are ordered before ConstantRateJumps (with the latter internally ordered in the same order they were passed in). The same principle applies for the construction of a `JumpSet`.  
