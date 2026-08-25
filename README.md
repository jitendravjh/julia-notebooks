# Geo notebooks

[Pluto](https://plutojl.org) notebooks built on the [JuliaEarth](https://github.com/JuliaEarth)
stack — real geospatial data, real projections, and pictures that explain something.

## Notebooks

### [The shape of the world](https://jitendravjh.in/geo-notebooks/shape-of-the-world.html)

Every flat map of a round planet is a lie, and Gauss proved you cannot avoid it. This notebook
morphs the world continuously between projections, then uses Tissot's indicatrix — small circles
drawn on the globe and projected along with the land — to show exactly what each projection
gives up.

Mercator keeps every circle a circle and lets them grow without limit. Gall–Peters keeps every
circle the same area and squashes them into ellipses. Robinson keeps neither, on purpose.

### [A round-trip atlas](https://jitendravjh.in/geo-notebooks/roundtrip-atlas.html)

Project a point, convert it back, and measure how far it moved. Doing that for every projection
at every point on the globe turns a table of numbers into a map of where each implementation
struggles — and the failures have structure. Robinson's error falls in horizontal bands sitting
on the five-degree spacing of its lookup table. Transverse Mercator's falls in two discs at the
singularity ninety degrees from its central meridian.

## Built with

[CoordRefSystems.jl](https://github.com/JuliaEarth/CoordRefSystems.jl) for projections,
[Meshes.jl](https://github.com/JuliaGeometry/Meshes.jl) for geometry,
[GeoArtifacts.jl](https://github.com/JuliaEarth/GeoArtifacts.jl) for Natural Earth data,
[Makie.jl](https://docs.makie.org) for drawing.

## Running locally

```julia
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. -e 'using Pluto; Pluto.run(notebook="notebooks/shape-of-the-world.jl")'
```
