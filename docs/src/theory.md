# Theory

## Exact Stochastic Simulation Algorithm (SSA) Without Delays

Consider a system consisting of $N \geq 1$ chemical species, $\{X_1,\ldots, X_N\}$, undergoing $M \geq 1$ chemical reactions through reaction channels $\{R_1,\ldots,R_M\}$, each of which is equipped with a propensity function (or intensity function in the mathematics literature), $a_k(X)$. The dynamic state of this chemical system can be described by the state vector $X(t) =[X_1(t),\ldots,X_N(t)]^T$, where $X_n(t),n = 1,\ldots,N,$ is the number of $X_n$ molecules at time $t$, and $[·]^T$ denotes the transpose of the vector in the bracket.

Following Gillespie [1], the dynamics of reaction $R_k$ defined by a state-change vector $\nu_k = [\nu_{1k} ,\ldots,\nu_{Nk}]^T$, where $\nu_{nk}$ gives the changes in the $X_n$ molecular population produced by one $R_k$ reaction, and a propensity function $a_k(t)$ together with the fundamental premise of stochastic chemical kinetics:

```math
\begin{equation}
\begin{aligned}
a_k(t)\Delta t = &\text{ the probability, given } X(t)=\mathbf{x}, \\
&\text{ that one reaction }R_k \text{ will occur in the}\\
&\text{ next infinitesimal time interval }[t,t+\Delta t].
\end{aligned}
\end{equation}
```

Defining the probability rate constant $c_k$ as the probability that a randomly selected combination of $R_k$ reactant molecules reacts in a unit time period, we can calculate  $a_k(t)$ from $c_k$ and the molecular numbers of $R_k$ reactants at time $t$ using the method given by Gillespie.

For a chemical system in a given state $X(t)=\mathbf{x}$ at time $t$, assuming that all reactions occur instantly, Gillespie’s exact SSA answers the following two questions: (i)  when will the next reaction occur?  (ii)  which reaction will occur? Specifically, Gillespie’s exact SSA simulates the following event in each step:

```math
\begin{equation}
\begin{aligned}
\text{E: } & \text{no reaction occurs in the time interval }[t,t+\Delta],\\
& \text{and a reaction }R_\mu \ \text{occurs in the infinitesimal}\\
& \text{time interval }[t+\Delta,t+\Delta+\Delta t].
\end{aligned}
\end{equation}
```

Based upon the fundamental premise Eq. ([1](#mjx-eqn-1)), Gillespie showed that that $\Delta$ and $\mu$ are two independent random variables and have the following probability density functions ($\text{PDFs}$), respectively:

```math
\begin{equation}
f_\Delta(\Delta)=a_0(t) \exp(-a_0(t)\Delta),\ \ \ \  \Delta > 0,
\end{equation}
```

and

```math
\begin{equation}
f_\mu(\mu)={{a_\mu(t)} \over {a_0(t)}},\ \ \ \ \  \mu = 1,\ldots,M,
\end{equation}
```

where $a_0(t)=\begin{matrix} \sum_{k=1}^M a_k(t) \end{matrix}$. According to the $\text{PDF}$ Eq. ([4](#mjx-eqn-4)), a realization of $\mu$ can be generated from a standard uniform random variable $r_2$, by taking $\mu$ to be the integer for which $\sum_{k=1}^{\mu-1} a_k(t)  < r_2 a_0(t) \leq \sum_{k=1}^\mu a_k(t)$;based on the $\text{PDF}$ Eq. ([3](#mjx-eqn-3)), a realization of $\Delta$ can be generated from another standard uniform random variable $r_1$ as $\Delta=−\ln(r_1)/a_0(t)$. Therefore, Gillespie’s exact SSA generates a realization of $\mu$ and $\Delta$ in each step of simulation, and then updates the time and system state as $t\leftarrow t+\Delta$ and  $\mathbf{x} \leftarrow \mathbf{x}+ \mathbf{\nu_\mu}$, respectively.

## Exact SSA For Coupled Chemical Reaction With Delays

### Delay Direct method

As in the derivation of Gillespie’s exact SSA, we first need to find the probability of event Eq. ([2](#mjx-eqn-2)), that is defined as $P(\Delta,\mu)d\Delta$, where $P(\Delta,\mu)$ is the joint PDF of $\Delta$ and $\mu$. Suppose that there are $d$ ongoing reactions at time $t$, which will finish at $t+T_1,\ldots,t+T_{d}$, respectively. Without loss of generality, we assume that $T_1 \leq T_2 \leq \ldots \leq T_d$. Unlike in the reaction system without delays where the propensity functions remain unchanged in the time interval $[t,t+\Delta]$, the propensity functions here change at $t+T_i,i=1,\ldots,d$, due to delayed reactions. We need to take into account such changes in propensity functions when deriving  $P(\Delta,\mu)$.

As in the derivation of Gillespie’s exact SSA, $P(\Delta,\mu)d\Delta$ can be found from the fundamental premise Eq. ([1](#mjx-eqn-1)) as

```math
\begin{equation}
P(\Delta,\mu)d\Delta=P_0(\Delta) a_\mu(t +\Delta)d\Delta,
\end{equation}
```

where $P_0(\Delta)$ is the probability that no reaction will occur in the time interval $[t,t+\Delta]$, while $a_\mu(t+\Delta)d\Delta$ is the probability that a reaction $R_\mu$ occurs in $[t+\Delta,t+\Delta+d\Delta]$. Defining $T_0=0$ and $T_{d+1}=\infty$, we can find $P_0(\Delta)$ for $\Delta$ that lies in different time intervals $[T_i,T_{i+1}),i=0,\ldots,d$. If $\Delta \in [T_i,T_i+1)$, we define the event $E_k$ as the event that no reaction occurs in the time interval $[t+T_k,t+T_{k+1}),k=0,\ldots,k=i−1$, respectively,and the event  $E_i$  as the event that no reaction occurs in the time interval $[t+T_i,t+\Delta)$. Then, we can express $P_0(\Delta)$ as

```math
\begin{equation}
P_0(\Delta)=P(E_0,\ldots,E_i)=P(E_0) \prod_{k=1}^i P(E_k丨E_0,\ldots,E_{k-1}).
\end{equation}
```

From the derivation of Gillespie’s exact SSA,we know that  
```math
P(E_0) = \exp (−a_0(t)T_1)\\
P_0(E_k丨E_0,\ldots,E_{k-1}) = \exp(-a_0(t+T_k)) × (T_{k+1}−T_k),k=0,\ldots,i−1,\\
P(E_i丨E_0,\ldots,E_{i-1}) = \exp(-a_0(t+T_i)(\Delta-T_i)).
```
Notice that propensity functions change at $t+T_k$ after a delayed reaction finishes, and we use $a_0(t+T_k)$ to represent the new $a_0$. The probability $P_0(\Delta)$ is then given by

```math
\begin{equation}
\begin{aligned}
& P_0(\Delta) = \exp \bigg (-\sum_{k=0}^{i-1} a_0(t+T_k)(T_{k+1}-T_k)-a_0(t+T_i)(\Delta-T_i) \bigg ), \\
& \Delta \in [T_i,T_{i+1}), i = 0,\ldots,d,
\end{aligned}
\end{equation}
```

where we assume that the first term of the exponent is equal to zero when $i = 0$. Since $P_0(\Delta)$ does not depend on individual propensity functions, as shown in Eq. ([7](#mjx-eqn-7)), it is seen from Eq. ([5](#mjx-eqn-5)) that $\Delta$ and $\mu$ are independent random variables. Combining Eq. ([5](#mjx-eqn-5)) and Eq. ([7](#mjx-eqn-7)) and noticing that $a_\mu(t+\Delta)=a_\mu(t+T_i)$ for $\Delta \in [T_i,T_i+1)$, we obtain the $\text{PDF}$ of $\Delta$ and $\mu$ as follows:

```math
\begin{equation}
\begin{aligned}
& f_\Delta(\Delta) = a_0(t+T_i) \exp \bigg (-\sum_{k=0}^{i-1} a_0(t+T_k)(T_{k+1}-T_k) - a_0(t+T_i)(\Delta-T_i) \bigg ), \\
& \Delta \in [T_i,T_{i+1}), i = 0,\ldots,d,
\end{aligned}
\end{equation}
```

and

```math
f_\mu(\mu)={ {a_\mu(t+T_i)} \over {a_0(t+T_i)} },\ \ \  \mu = 1,\ldots,M,\ \ \  \Delta \in [T_i,T_{i+1}),
```

It is not difficult to verify that $\int_{0}^{\infty} f_\Delta(\Delta)\, d\Delta = 1$. In simulation, $\mu$ can be generated, from a standard uniform random variable $u_1$, by taking $\mu$ to be the integer for which $\begin{matrix} \sum_{k=1}^{\mu-1} a_k(t+T_i) \end{matrix} < u_1 a_0(t+T_i) ≤ \begin{matrix} \sum_{k=1}^\mu a_k(t+T_i) \end{matrix}$, after $\Delta$ is generated to be in the time interval $[T_i,T_{i+1})$.We next derive the method of generating  $\Delta$ according to its $\text{PDF}$ in Eq. ([8](#mjx-eqn-8)).

The cumulative distribution function of $\Delta$can be found from Eq. ([8](#mjx-eqn-8)) as

```math
\begin{aligned}
& F_\Delta(\Delta)=1 - \exp  \bigg (-\sum_{k=0}^{i-1} a_0(t+T_k)(T_{k+1}-T_k)-a_0(t+T_i)(\Delta-T_i) \bigg ), \\
& \Delta \in [T_i,T_{i+1}), i = 0,\ldots,d,
\end{aligned}
```

Then, we can generate $\Delta$ from a standard uniform random variable $u_2$, by taking $\Delta=F_\Delta^{−1}(u_2)$, where $F_\Delta^{−1}(\cdot)$ represents the inverse of $F_\Delta(\Delta)$. More specifically, we can obtain $\Delta$ as follows:

Find $T_i$ such that  $F_\Delta(T_i) ≤ u_2 ≤ F_\Delta(T_{i+1})$, then calculate  $\Delta$ from

```math
\begin{aligned}
& \Delta = T_i + {{-\ln (1-u_2)-\begin{matrix} \sum_{k=0}^{i-1} a_0(t+T_k)(T_{k+1}-T_k) \end{matrix} } \over {a_0(t+T_i)}} \\
& \Delta \in [T_i,T_{i+1}).
\end{aligned}
```

Since we need $T_1,\ldots,T_d$ to generate $\Delta$ and $\mu$, we define an array of data structures, named *Tstruct*, whose $i$th $(i=1,\ldots,d)$ cell stores $T_i$ and the index, $\mu_i$, of the reaction that $T_i$ is associated with. The reaction index $\mu_i$ is needed during the generation of $\Delta$, when we update the propensity functions affected by the reaction that is delayed but finishes at $t+T_i$. During simulation, we need to generate $\Delta$ and $\mu$, maintain *Tstruct*, and then update the state vector $X(t)$.


### Delay Rejection method
Now let us see whether the rejection method can correctly simulate the event ([2](#mjx-eqn-2)). The rejection algorithm essentially generates $\Delta$ in the event ([2](#mjx-eqn-2)) using a rejection method in an iterative fashion: in the *i*-th iteration, it generates a $\Delta_i$ iaccording to an exponential $\text{PDF}$ with parameter $a_0(t+T_{i−1})$, where we have denoted the $\Delta$ generated in the *i*-th iteration as $\Delta_i$. If $\Delta_i < T_i - T_{i−1}$, then we have $\Delta = \sum_{k=0}^{i-1} T_k+\Delta_i$ and the algorithm continues simulation to generate $\mu$; otherwise, it rejects $\Delta_i$, updates the state vector $X(t+T_i)$, calculates $a_k(t+T_i),k=1,\ldots,M$, and goes to the next iteration. If $\Delta$ is determined in the *(i+1)*-th iteration, where *i*
is a non-negative integer, then we have $\Delta \in [T_i,T_{i+1})$ and *i* delayed reactions finished in the time interval $[t,t+\Delta)$. 

From the iterative procedure of generating $\Delta$ described
above, we can find $P_0(\Delta)$ that the rejection method algorithm produces. Specifically, if $\Delta \in [T_i,T_{i+1})$, we have $P(E_0)=P(\Delta_1 > T_1)$, $P(E_k丨E_0,\ldots,E_{k-1}) = P(\Delta_{k+1} > T_{k+1} - T_k), k=1,\ldots,i−1$, because $\Delta_k,k=1,\ldots,i$, are rejected. Since $\Delta_{k+1}$ is accepted, at least one reaction will occur in the time interval $[t+T_i,t+\Delta)$, if $\Delta_{i+1} < \Delta −T_i$. Thus, $P(E_i丨E_0,\ldots,E_{i-1}) = 1−P(\Delta_{i+1} < \Delta - T_i) = P(\Delta_{i+1} > \Delta - T_i)$. Therefore, for the rejection method, $P_0(\Delta)$ in Eq. ([6](#mjx-eqn-6)) can be written as
```math
\begin{equation}
P_0(\Delta) = P(\Delta_{i+1} > \Delta - T_i) \prod_{k=1}^i P(\Delta_k > T_k - T_{k-1}).
\end{equation}
```
The random variables  $\Delta_k,k=1,\ldots,i+1$, follow an exponential distribution with parameter $a_0(t+T_{k−1})$, and thus we have
```math
\begin{equation}
\begin{aligned}
& P(\Delta_k > T_k - T_{k-1}) = \exp(-a_0(t+T_{k-1})(T_k - T_{k-1})) \\
& k= 1,\ldots,i,
\end{aligned}
\end{equation}
```
and
```math
\begin{equation}
P(\Delta_{i+1} > \Delta −T_i) = \exp(-a_0(t+T_i)(\Delta-T_i)).
\end{equation}
```
Substituting Eqs. ([10](#mjx-eqn-10)) and ([11](#mjx-eqn-11)) into Eq. ([9](#mjx-eqn-9)), we find that $P_0(\Delta)$ in Eq. ([9](#mjx-eqn-9)) is exactly the same as $P_0(\Delta)$ in Eq. ([7](#mjx-eqn-7)) that is derived directly from the event ([2](#mjx-eqn-2)) and the fundamental premise ([1](#mjx-eqn-1)). Since our algorithm generates $\Delta$ and $\mu$ according to $\text{PDF}$s of $\Delta$ and $\mu$ derived from $P_0(\Delta)$ in Eq. ([7](#mjx-eqn-7)), the rejection method is equivalent to our direct method and also is an exact SSA for chemical reaction systems with delays.

We next analyze the complexity of the rejection method algorithm and the direct method algorithm. As we have seen, the difference between two algorithms lies in the generation of $\Delta$. Suppose that both algorithms generate a  $\Delta \in [T_i,T_{i+1})$, where *i* is a non-negative integer. Both algorithms update the state vector **x** at $t+T_k$, $k=1,\ldots,i$, and calculate $a_k(t+T_k)$ and $a_0(t+T_k)$, $k=1,\ldots,i$. The direct method algorithm also calculates $a_k(t+T_{i+1})$ and $a_0(t+T_{i+1})$, but $a_k(t+T_{i+1})$ and $a_0(t+T_{i+1})$ can be reused in generating next $\Delta$. Therefore, two algorithms require the same computation on calculating propensity functions and updating the state vector. The direct method algorithm needs to evaluate the exponential function and calculateat $a_t$ *i+1*  times, while the rejection method does not need such operations. Also, the direct method needs slightly more computation on calculating $\Delta$ from a uniform random variable than the rejection method algorithm. To generate a $\Delta$, the direct method algorithm generates exactly one uniform random variable regardless of the value of *i*, while the rejection method algorithm generates *i*+1 uniform random variable. 

## Reference

[1] Daniel T. Gillespie, "Exact stochastic simulation of coupled chemical reactions", The Journal of Physical Chemistry 1977 81 (25), 2340-2361.
[https://doi.org/10.1021/j100540a008](https://doi.org/10.1021/j100540a008).