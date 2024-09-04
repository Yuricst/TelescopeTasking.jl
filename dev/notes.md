# Misc notes


## 2024.09.04

### Formulation

Suppose we have $n$ observation corresponding to $m$ targets. Note that necessarily $n\geq m$ (and likely $n > m$ since most targets have multiple passes over the ground station over a single night). 
We consider a target to be observed if it has been observed at least over $E$ arcs (i.e. $E$ exposures). 
We consider the decision variables:

- $Y \in \mathbb{B}^n$ : $Y_i = \{0,1\}$ dictates whether to observe along arc $i$
- $\theta \in \mathbb{B}^m$ : $\theta_k = \{0,1\}$ dictates whether target $k$ is observed

We also precompute the following coefficient matrices:

- $A \in \mathbb{B}^{n \times m}$ : $A_{ik} = \{0,1\}$ dictates whether target $k$ is observed by observation arc $i$; note that since $n \geq m$, this matrix is tall (or square if $n = m$).
- $T \in \mathbb{B}^{n \times n}$ : $T_ij = \{0,1\}$ dictates whether slewing from arc $i$ to arc $j$ is possible (from a hardware point of view, taking into account slewing etc.)

Then, the telescope observation scheduling problem (TOSP) is given by

$$
\begin{aligned}
\min_{Y,\theta} \quad& \sum_{k=1}^m \theta_k
\\
\text{such that} \quad
&\sum_{i=1}^n A_{ik} Y_i \geq E \theta_k \quad \forall k=1,\ldots,m
\\
&Y_i + Y_j \leq 1 + T_{ij} \quad \forall i = 1,\ldots,n-1, \quad j = i+1,\ldots,n
\end{aligned}
$$

### Code organization

- Make ordered vector `ObservationArc[]`, where each `ObservationArc` is a struct containing:
    - Initial, max elevation, and final observation time
    - Initial, max elevation, and final observation azimuth & elevation
    - Target index `k` and name (TLE name)