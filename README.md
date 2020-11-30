### Pseudo Wigner-Ville Distribution for discrete and finite real sampled signals

###### Alejandro R. Urzúa (2020)

_________

This repository contains the research on the application of Pseudo Wigner-Ville Distribution, hereafter **PWVD**, for the analysis of real sampled signals who are discrete ($N$-time points) and finite ($N<\infin$, bounded).

This research started as of 2012 when I was working on my [bachelor's thesis](https://doi.org/10.13140/RG.2.1.3969.5122) implementing PWVD for the analysis of vibrational machine signals. Since then, the algorithms and analysis go stacked. Now, on lated 2020, I want to take again the subject and work it out.

(Possible) Roadmap:

1. Translate the `Fortran` and `Mathematica` codes to `Julia` 

2. Optimize for $N$ values greater that $8192$. (Who is the maximum data points achieved with some reasonable performance with the `Fortran` code. `Mathematica` can manage up to $4096$ data points, but in a ridiculous long time.)

3. Implement some other discrete filters and smooth-windows for the final representation. (Currently, _Hamming_ filter and _Gaussian_ smoothing window are implemented.)

4. Extend the analysis for datasets with $N$ unrestricted to $2^n$. (Condition imposed by the discrete and finite Fast Fourier Transform.)

5. _Some other unreal goals goes here_

~~Theory goes here~~
