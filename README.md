# Julia notebooks

[Pluto](https://plutojl.org) notebooks, mostly about geometry.

## Notebooks

### [A conjecture that was false](https://jitendravjh.in/julia-notebooks/jacobian-conjecture.html)

Keller asked in 1939 whether a polynomial map with an everywhere-invertible derivative must be
invertible. It held for eighty-seven years and broke in July 2026.

The notebook verifies the counterexample exactly: its Jacobian determinant is the constant −2,
and three different points share an image. Then it shows the geometry behind it, where a grid
pushed through a tangent sweep folds along the curve and covers the inside twice.

### [The middle of nowhere](https://jitendravjh.in/julia-notebooks/middle-of-nowhere.html)

The point in Australia farthest from any town. Searching the continent is unnecessary, because
the answer always sits on a Voronoi vertex, and there are only a few hundred of those.

It lands in the Great Victoria Desert, 576 km from the nearest town, with three towns tied at
that distance because the circle is resting on all three.

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
