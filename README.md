# Julia notebooks

[Pluto](https://plutojl.org) notebooks on geometry and computation. Each one is re-run from a
clean environment on every push, so the output on the site is the output of the code here.

## Notebooks

### [The Kepler conjecture](https://jitendravjh.in/julia-notebooks/kepler-conjecture.html)

No arrangement of equal spheres fills more than π/√18 of space. Kepler asserted it in 1611
without argument, Hales proved it in 1998, and the Flyspeck project finished a machine-checked
proof in 2014.

The notebook computes the densities, builds the rhombic dodecahedron that is the Voronoi cell of
the packing, and shows why the obvious local argument cannot work: a single sphere can sit in a
cell filling 0.7547 of its volume, above Kepler's 0.7405.

### [The Jacobian conjecture in dimension three](https://jitendravjh.in/julia-notebooks/jacobian-conjecture.html)

Keller asked in 1939 whether a polynomial map with an everywhere-invertible derivative must
itself be invertible. It was refuted for three variables in July 2026, while the two-variable
case he posed remains open.

The counterexample is verified in exact rational arithmetic, then explained: multiplication of
binary forms is generically three to one, and running a loop through the space of cubics
permutes the three branches, so no global inverse can exist.

### [The largest empty circle](https://jitendravjh.in/julia-notebooks/largest-empty-circle.html)

Which point is farthest from every site in a finite set? The answer sits at a Voronoi vertex or
on the convex hull, which turns an optimisation over the plane into a search over a few hundred
candidates.

Applied to Australian towns it lands in the Great Victoria Desert, 576 km from the nearest, with
three towns tied because the circle rests on all three.

## Built with

[Meshes.jl](https://github.com/JuliaGeometry/Meshes.jl),
[CoordRefSystems.jl](https://github.com/JuliaEarth/CoordRefSystems.jl),
[GeoArtifacts.jl](https://github.com/JuliaEarth/GeoArtifacts.jl),
[DynamicPolynomials.jl](https://github.com/JuliaAlgebra/DynamicPolynomials.jl) and
[Makie.jl](https://docs.makie.org).

## Running locally

```julia
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. -e 'using Pluto; Pluto.run(notebook="notebooks/kepler-conjecture.jl")'
```
