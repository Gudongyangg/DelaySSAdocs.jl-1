# An example of stochastic delay 

## Model
According to [[1]](https://www.nature.com/articles/s41467-021-22919-1), the model assumes that the gene can switch between active $G$ and inactive $G^*$ states, transcribes nascent mRNA (denoted by  $N$) while in the active state which subsequently is removed a time $\tau$ later. Here the delay $\tau$ can be a random variable. We can define the model
```math
G^*\xrightarrow{\sigma_{\text{on}}} G\\
G\xrightarrow{\sigma_{\text{off}}}G^*\\
G\xrightarrow{\rho}G+N
```
and $G\xrightarrow{\rho}G+N$ will trigger $N\Rightarrow \emptyset$ after $\tau$ time.
We set $\sigma_{\text{on}}=0.0282$, $\sigma_{\text{off}}=0.609$ and $\rho=2.11$.

### Markovian part

We first define the parameters and the mass-action jump (see [Defining a Mass Action Jump](https://diffeq.sciml.ai/stable/types/jump_types/#Defining-a-Mass-Action-Jump) for details).

```julia
rate1 = [0.0282, 0.609, 2.11]
reactant_stoich = [[2=>1],[1=>1],[1=>1]]
net_stoich = [[1=>1,2=>-1],[1=>-1,2=>1],[3=>1]]
mass_jump = MassActionJump(rate1, reactant_stoich, net_stoich; scale_rates =false)
jumpset = JumpSet((),(),nothing,mass_jump)
```

Then we set the Initial value and define a `DiscreteProblem`.

```julia
u0 = [1,0,0]
de_chan0 = [[]]
tf = 2000.
tspan = (0,tf)
dprob = DiscreteProblem(u0, tspan)
```
### Non-Markovian part
Unlike other examples, the elongation time $\tau$ is a random variable sampled from two different lognormal distributions. We assume $\tau\sim \text{LogNormal}(0,2)+120$ and $\tau\sim \text{LogNormal}(1,\sqrt{2})+120$. Here we take  $\tau\sim \text{LogNormal}(1,\sqrt{2})+120$ for example.

```julia
delay_trigger_affect! = function (integrator, rng)
    τ=rand(LogNormal(1,sqrt(2)),1)+120
    append!(integrator.de_chan[1], τ)
end
delay_trigger = Dict(3=>delay_trigger_affect!)
delay_complete = Dict(1=>[3=>-1]) 
delay_interrupt = Dict() 
delayjumpset = DelayJumpSet(delay_trigger,delay_complete,delay_interrupt)
```
So we can define the problem
```julia
djprob = DelayJumpProblem(dprob, DelayRejection(), jumpset, delayjumpset, de_chan0, save_positions=(false,false))
```

## Visualisation

Then we simulate $10^5$ trajectories and calculate the probability distribution.
```julia
using DiffEqJump
ensprob = EnsembleProblem(djprob)
@time ens = solve(ensprob, SSAStepper(), EnsembleThreads(), trajectories=10^5)
```

## Reference

[1] Qingchao Jiang, Xiaoming Fu, Shifu Yan, Runlai Li, Wenli Du, Zhixing Cao, Feng Qian, Ramon Grima, "Neural network aided approximation and parameter inference of non-Markovian models of gene expression". Nature communications, (2021) 12(1), 1-12. https://doi.org/10.1038/s41467-021-22919-1