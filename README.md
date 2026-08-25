# Julia notebooks

[Pluto](https://plutojl.org) notebooks on geometry and computation. Each one is re-run from a
clean environment on every push, so the output on the site is the output of the code in this
repository.

## Notebooks

### [The Jacobian conjecture in dimension three](https://jitendravjh.in/julia-notebooks/jacobian-conjecture.html)

Keller asked in 1939 whether a polynomial map with an everywhere-invertible derivative must
itself be invertible. The conjecture was refuted for three variables in July 2026, while the
two-variable case he posed remains open.

The notebook verifies the counterexample in exact rational arithmetic, then accounts for it:
multiplication of binary forms is generically three to one, and a tangent sweep shows the same
multiplicity in two dimensions, where it can be drawn.

### [The largest empty circle](https://jitendravjh.in/julia-notebooks/largest-empty-circle.html)

Which point is farthest from every site in a finite set? The answer is centred at a Voronoi
vertex or on the convex hull, which turns an optimisation over the plane into a search over a
few hundred candidates.

Applied to Australian towns it lands in the Great Victoria Desert, 576 km from the nearest, with
three towns tied at that distance because the circle rests on all three.

## Built with

[Meshes.jl](https://github.com/JuliaGeometry/Meshes.jl),
[CoordRefSystems.jl](https://github.com/JuliaEarth/CoordRefSystems.jl),
[GeoArtifacts.jl](https://github.com/JuliaEarth/GeoArtifacts.jl),
[DynamicPolynomials.jl](https://github.com/JuliaAlgebra/DynamicPolynomials.jl) and
[Makie.jl](https://docs.makie.org).

## Running locally

```julia
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. -e 'using Pluto; Pluto.run(notebook="notebooks/jacobian-conjecture.jl")'
```
