# Main API
```@meta
CurrentModule = DelaySSAdocs
```

```@docs
DelayJumpSet
DelayJumpProblem
```

## More about defining a `DelayJumpSet`
Further notes on `delay_trigger`, `delay_interrupt`, `delay_complete`.
### `delay_trigger` 
   delay_trigger defines a `Dict` type:

- Keys: Indices of reactions defined in the Markovian part that can trigger the delay reaction;
- Values: A update function or a `Pair` type that determines how to update the delay channel.
 

### `delay_interrupt`
   delay_interrupt defines a `Dict` type:
   
- Keys: Indices of reactions defined in Markovian part that can cause the change in the delay channels;
- Values: A update function that determines how to update the delay channel

### `delay_complete` 
   delay_interrupt defines a `Dict` type:
   
- Keys: Indices of delay channels;
- Values: A vector of `Pair`s, mapping species index to net change of stoichiometric coefficient.

## Types and Algorithms
```@docs
AbstractDelayAggregatorAlgorithm
DelayDirect
DelayRejection
DelayMNRM
DelayDirectCR
```