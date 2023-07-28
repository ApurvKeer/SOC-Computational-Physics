# Summer Of Code: Hands-on Computational Physics

Details:
Apurv Keer | 210260007 | BTech 3rd Year Engineering Physics

## Numerical Methods in Computational Physics

This involved implementation of some [numerical techniques](/Numerical%20Methods/), like the Euler, Jacobi, Gauss-Seidel, Newton-Raphson, and Secant methods. This is illustrated in a [jupyter notebook](/Numerical%20Methods/Apurv_Numerical_Method.ipynb). There are also Fortran codes for these methods.

## Visualizing Poisson Equation Solutions

This was the first Mini-Project. The project involved solving the [Poisson equation](/Poisson%20Equation/) using the finite difference method. There is a [jupyter notebook](/Poisson%20Equation/Apurv_poisson_solutions.ipynb) and a [Fortran code](/Poisson%20Equation/poisson.f90) for this. The Fortran code produces a text file with pixel position and value.

## Visualizing Hydrogen atom wavefunctions

This was the second Mini-Project. The [project](/Hydrogen%20Atom/) involved numerically solving for eigenstates and eigenvalues of the Schrödinger equation for hydrogen atom. Similar content as mini-project one.

## Visualizing Time Evolution systems

This was the third Mini-Project. The project involves solving [time evolution](/Time%20Evolutions/) systems of equations using Runge-Kutta Algorithms.

## Raytracing in Curved Spacetime (Kerr Blackhole)

This is an open-ended project. The project involves visualizing a Kerr Blackhole by solving null geodesics around the spacetime. I have written the code in C++ and Fortran. The code prints an image in ppm P3 format. To save it, run

```console
$ ./a.out > image.ppm
```

This was the final result.

![Final Render](/final_render_1.png)

##BlackRay

I wrote a better raytracer. This one ray traces a Schwarzschild Black Hole. This was based on the previous work [Starless]{https://rantonels.github.io/starless/}, by Rantonels The geodesics equation was written in the form of a Binet equation given by,

$$ \frac{d^2}{dx^2} \vec{r} = -\frac{3}{2} h^2 \frac{\hat{r}}{r^5} $$

I have uv wrapped two images, one on the accretion disk and the other on the celestial sphere. This was the result.














