# EXAFS-data-reduction
Compiled and tested with Turbo C++
This package was developed to reduce EXAFS data, initially on Turbo C/C++, so DOS-like GUI was provided. 

It includes two parts: 

The fist part commits the pre-edge subtraction, post-edge fitting, E-to-k conversion, 
k-space weighting, windowing, fast Fourier Transform, and inverse Fourier Transform for amplitude and phase function extraction.

The second part is designed to extract local coordiation parameters (including coordination number, coordination distance, disorders,
energy shift) from curve fitting based on simulated annealling approach.
