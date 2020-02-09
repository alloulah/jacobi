# jacobi

Years back, I worked on inter-carrier interference (ICI) equalisation for terrestrial broadcast channels. In terrestrial broadcast OFDM, the number of subcarriers is typically much larger than say WiFi e.g. 1024 to 32768. Within the ICI problem formulation, we typically need to invert a large channel matrix in real-time on an embedded baseband transceiver. LU decomposition can be a bit expensive in terms of compute, but luckily approximations such as the Jacobi iterative method come handy.

This repo has a skeleton code for solving **iteratively** a problem of the form

```math
\textbf{y} = \textbf{H} \textbf{x}
```
whose solution is
```math
\textbf{x} = \textbf{H}^{-1} \textbf{y}
```

using Jacobi iteration as described by [Molisch et al.](https://doi.org/10.1109/TVT.2007.897628)
