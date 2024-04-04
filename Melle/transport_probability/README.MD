# Solving the Fokker Planck equation
The Fokker-Planck operator is given by

$$ \mathcal{L}^* \rho = -\partial_x(\rho f) + \frac{\sigma^2}{2}\partial_{xx} \rho

Where in our application $f$ is given by the Stommel 2-box model

$$ \dot{x} = f(x) = \lambda - |1-x|x $$

The discrete Fokker-Plack operator is now given by

$$L = -D_1 \text{diag}(f) + \frac{\sigma^2}{2}D_2$$

Where $D_1$ and $D_2$ are the differences matrices approximating the first and second derivative respectively. If we take period boundary conditions the matrices are given by

$$ D_1 = \begin{pmatrix}
	 &  &  &  &  \\
	 &  &  &  &  \\
	 &  &  &  &  \\
	 &  &  &  &  \\
	 &  &  &  &  \\
\end{pmatrix}
$$
