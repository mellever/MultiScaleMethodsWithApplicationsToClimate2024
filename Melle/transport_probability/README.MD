# Solving the Fokker Planck equation
The Fokker-Planck operator is given by

$$ \mathcal{L}^* \rho = -\partial_x(\rho f) + \frac{\sigma^2}{2}\partial_{xx} \rho

Where in our application $f$ is given by the Stommel 2-box model

$$ \dot{x} = f(x) = \lambda - |1-x|x $$

The discrete Fokker-Plack operator is now given by

$$L = -D_1 \text{diag}(f) + \frac{\sigma^2}{2}D_2$$

Where $D_1$ and $D_2$ are the differences matrices approximating the first and second derivative respectively. We choose to approximate the derivatives by central differences. If we take period boundary conditions the matrices are given by

$$ 
D_1 = \frac{1}{2\Delta x} \begin{pmatrix}
0 & 1 & 0 & \cdots & 0 & 0 & -1 \\
-1 & 0 & 1 & \cdots & 0 & 0 & 0 \\
0 & -1 & 0 & \cdots & 0 & 0 & 0 \\
\vdots & \vdots & \vdots & \ddots & \vdots & \vdots & \vdots \\
0 & 0 & 0 & \cdots & 0 & 1 & 0 \\
0 & 0 & 0 & \cdots & -1 & 0 & 1 \\
1 & 0 & 0 & \cdots & 0 & -1 & 0 \\
\end{pmatrix}
$$

$$ D_2 =  \frac{1}{\Delta x^2} \begin{pmatrix}
-2 & 1 & 0 & \cdots & 0 & 0 & 1 \\
1 & -2 & 1 & \cdots & 0 & 0 & 0 \\
0 & 1 & -2 & \cdots & 0 & 0 & 0 \\
\vdots & \vdots & \vdots & \ddots & \vdots & \vdots & \vdots \\
0 & 0 & 0 & \cdots & -2 & 1 & 0 \\
0 & 0 & 0 & \cdots & 1 & -2 & 1 \\
1 & 0 & 0 & \cdots & 0 & 1 & -2 \\
\end{pmatrix} $$

And diag(f) is simply the evaluation of $f$ on the gridpoints and putting these of the diagonal. Making use of this, we can find the probability density via

$$\rho(t+\Delta t) = \exp{\Delta t L} \rho(t)$$

Note that we have to choose the domain large enough such that the probability density is small at the boundaries.



