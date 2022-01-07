"""
$(TYPEDEF)

A delay jump set that consists of five inputs at most, namely `delay_trigger`, `delay_interrupt`, `delay_complete`, `delay_trigger_set` and `delay_interrupt_set`.

# Fields
$(FIELDS)

# Notes
- `delay_trigger::Dict{Int,Any}`:    those reactions in the Markovian part that trigger the change of the state of the delay channels or/and the state of the reactants upon initiation. 
    - Keys: Indices of reactions defined in the Markovian part that can trigger the delay reactions; 
    - Values: Update functions(or `Pair` type) that decide how to update the delay channel or the state of the reactants.

- `delay_interrupt::Dict{Int,Any}`:  those reactions in the Markovian part that change the state of the delay channels or/and the state of the reactants in the middle of on-going delay reactions. 
    - Keys: Indices of reactions defined in the Markovian part that can interrupt the delay reactions; 
    - Values: Update functions(or `Pair` type) that decide how to update the delay channel or the state of the reactants.

- `delay_complete::Dict{Int,Any}`:     those reactions that are initiated by delay trigger reactions and change the state of the delay channels or/and the state of the reactants upon completion. 
    - Keys: Indices of the delay channel; 
    - Values: Update functions(or `Pair` type) that decide how to update the delay channel or the state of the reactants.

- `delay_trigger_set::Vector{Int}`:  indices of reactions that can trigger the delay reaction.

- `delay_interrupt_set::Vector{Int}`:  indices of reactions that can interrupt the delay reactions.

We take this [model](https://palmtree2013.github.io/DelaySSAdocs.jl/dev/tutorials/delay_degradation/) for example.
```julia
delay_trigger_affect! = function (integrator, rng)
  append!(integrator.de_chan[1], τ)
end

# the 3rd reaction will trigger the delay reaction
delay_trigger = Dict(3=>delay_trigger_affect!)

# the 1st delay reaction will cause the 2nd species of molecule to degrade
delay_complete = Dict(1=>[2=>-1]) 
delay_affect! = function (integrator, rng)
   i = rand(rng, 1:length(integrator.de_chan[1]))
   deleteat!(integrator.de_chan[1],i)
end

# the 4th reaction will interrupt the delay reactions
delay_interrupt = Dict(4=>delay_affect!) 
delaysets = DelayJumpSet(delay_trigger,delay_complete,delay_interrupt)
```


"""
mutable struct DelayJumpSet
    """those reactions in the Markovian part that trigger the change of the state of the delay channels or/and the state of the reactants upon initiation."""
    delay_trigger::Dict{Int,Any}
    """those reactions in the Markovian part that change the state of the delay channels or/and the state of the reactants in the middle of on-going delay reactions."""
    delay_complete::Dict{Int,Any}
    """those reactions that are initiated by delay trigger reactions and change the state of the delay channels or/and the state of the reactants upon completion."""
    delay_interrupt::Dict{Int,Any}
    """keys of `delay_trigger`."""
    delay_trigger_set::Vector{Int}
    """keys of `delay_interrupt`."""
    delay_interrupt_set::Vector{Int}    
end

DelayJumpSet(delay_trigger,delay_complete,delay_interrupt) = DelayJumpSet(delay_trigger,delay_complete,delay_interrupt, collect(keys(delay_trigger)), collect(keys(delay_interrupt)))




"""
    function DelayJumpProblem(prob::DiscreteProblem, aggregator::AbstractDelayAggregatorAlgorithm, jumps::JumpSet, delayjumpsets::DelayJumpSet, de_chan0)
# Fields
- `prob::DiscreteProblem`

    A problem defined by some initial values in Markovian part.
- `aggregator::AbstractDelayAggregatorAlgorithm`

    The algorithm chose to solve the DelaySSA problem.
- `jumps::JumpSet`

    A jump set containing the information of Markovian part.
- `delayjumpsets::DelayJumpSet`

    A jump set containing the information of Non-Markovian part.
- `de_chan0::Vector{Vector{T}}` 

    Initial condition of the delay channel and `T::Float64`.
```julia
djprob = DelayJumpProblem(dprob, DelayRejection(), jumpset, delayjumpset, de_chan0, save_positions=(true,true)).
```
"""
function DelayJumpProblem(prob, aggregator, jumps, delayjumpsets::DelayJumpSet, de_chan0;
                     save_positions = (true,true),
                     rng = Xorshifts.Xoroshiro128Star(rand(UInt64)), scale_rates = true, useiszero = true, spatial_system=nothing, hopping_constants=nothing, kwargs...)

  # initialize the MassActionJump rate constants with the user parameters
  if using_params(jumps.massaction_jump) 
    rates = jumps.massaction_jump.param_mapper(prob.p)
    maj = MassActionJump(rates, jumps.massaction_jump.reactant_stoch, jumps.massaction_jump.net_stoch, 
                         jumps.massaction_jump.param_mapper; scale_rates=scale_rates, useiszero=useiszero, 
                         nocopy=true)
  else
    maj = jumps.massaction_jump
  end

  ## Spatial jumps handling
  if spatial_system !== nothing && hopping_constants !== nothing && !is_spatial(aggregator) # check if need to flatten
    prob, maj = flatten(maj, prob, spatial_system, hopping_constants; kwargs...)
  end
  ## Constant Rate Handling
  t,end_time,u = prob.tspan[1],prob.tspan[2],prob.u0
  if (typeof(jumps.constant_jumps) <: Tuple{}) && (maj === nothing) && !is_spatial(aggregator) # check if there are no jumps
    disc = nothing
    constant_jump_callback = CallbackSet()
  else
    disc = aggregate(aggregator,u,prob.p,t,end_time,jumps.constant_jumps,maj,save_positions,rng; spatial_system = spatial_system, hopping_constants = hopping_constants, kwargs...)
    constant_jump_callback = DiscreteCallback(disc)
  end

  iip = isinplace_jump(prob, jumps.regular_jump)

  ## Variable Rate Handling
  if typeof(jumps.variable_jumps) <: Tuple{}
    new_prob = prob
    variable_jump_callback = CallbackSet()
  else
    new_prob = extend_problem(prob,jumps)
    variable_jump_callback = build_variable_callback(CallbackSet(),0,jumps.variable_jumps...)
  end
  callbacks = CallbackSet(constant_jump_callback,variable_jump_callback)

  DelayJumpProblem{iip,typeof(new_prob),typeof(aggregator),typeof(callbacks),
              typeof(disc),typeof(jumps.variable_jumps),
              typeof(jumps.regular_jump),typeof(maj),typeof(delayjumpsets),typeof(de_chan0)}(
                        new_prob,aggregator,disc,
                        callbacks,
                        jumps.variable_jumps,
                        jumps.regular_jump, maj, delayjumpsets, de_chan0)
end

"""
    function DelayJumpProblem(js::JumpSystem, prob, aggregator, delayjumpset, de_chan0; kwargs...)
# Fields
- `js::JumpSystem`    

    A jump system containing the information of Markovian part.
- `prob::DiscreteProblem`
- `aggregator::AbstractDelayAggregatorAlgorithm`
- `delayjumpsets::DelayJumpSet`
- `de_chan0::Vector{Vector{T}}` 
```julia
djprob = DelayJumpProblem(jumpsys, dprob, DelayRejection(), delayjumpset, de_chan0, save_positions=(true,true))
```
"""
function DelayJumpProblem(js, prob, aggregator, delayjumpset, de_chan0; kwargs...)
    statetoid = Dict(value(state) => i for (i,state) in enumerate(states(js)))
    eqs       = equations(js)
    invttype  = prob.tspan[1] === nothing ? Float64 : typeof(1 / prob.tspan[2])

    # handling parameter substition and empty param vecs
    p = (prob.p isa DiffEqBase.NullParameters || prob.p === nothing) ? Num[] : prob.p

    majpmapper = JumpSysMajParamMapper(js, p; jseqs=eqs, rateconsttype=invttype)
    majs = isempty(eqs.x[1]) ? nothing : assemble_maj(eqs.x[1], statetoid, majpmapper)
    crjs = ConstantRateJump[assemble_crj(js, j, statetoid) for j in eqs.x[2]]
    vrjs = VariableRateJump[assemble_vrj(js, j, statetoid) for j in eqs.x[3]]
    ((prob isa DiscreteProblem) && !isempty(vrjs)) && error("Use continuous problems such as an ODEProblem or a SDEProblem with VariableRateJumps")
    jset = JumpSet(Tuple(vrjs), Tuple(crjs), nothing, majs)

    if needs_vartojumps_map(aggregator) || needs_depgraph(aggregator)
        jdeps = asgraph(js)
        vdeps = variable_dependencies(js)
        vtoj = jdeps.badjlist
        jtov = vdeps.badjlist
        jtoj = needs_depgraph(aggregator) ? eqeq_dependencies(jdeps, vdeps).fadjlist : nothing
    else
        vtoj = nothing; jtov = nothing; jtoj = nothing
    end
    DelayJumpProblem(prob, aggregator, jset, delayjumpset, de_chan0; dep_graph=jtoj, vartojumps_map=vtoj, jumptovars_map=jtov, scale_rates=false, nocopy=true, kwargs...)
end

 
