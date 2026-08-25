# Projection round-trip atlas

Every map projection in [CoordRefSystems.jl](https://github.com/JuliaEarth/CoordRefSystems.jl)
makes a quiet promise: project a latitude and longitude, convert it back, and you should land
where you started.

This is a [Pluto](https://plutojl.org) notebook that measures how well that promise holds, for
every projection, at every point on the globe, and draws the answer as a map.

**Read it here: https://jitendravjh.github.io/projection-atlas**

## What it shows

Round-trip displacement is measured as **great-circle distance in metres** between the original
point and the recovered one. That matters: comparing latitude and longitude componentwise makes
every projection look broken at the poles, where longitude is undefined and any value is the
same physical point.

Measured over a 1° global grid against CoordRefSystems v0.19.24:

| tier | projections | worst displacement |
|---|---|---|
| exact | PlateCarree, WebMercator, Mercator, Sinusoidal, WinkelTripel | 3–6 nm |
| pole-limited | OrthoNorth, EqualEarth, LambertAzimuthal, GallPeters, Albers | 0.10–0.31 m |
| **outliers** | **Robinson**, **TransverseMercator** | **39 m**, **217 m** |

The two outliers are the two projections with open accuracy issues, and the atlas shows *where*
each one fails:

- **Robinson** is defined by a table at every 5° of latitude rather than a formula. Its error
  appears as horizontal bands sitting on the tabulation knots, which is the `Float32`
  interpolation-coefficient discontinuity described in
  [issue #55](https://github.com/JuliaEarth/CoordRefSystems.jl/issues/55).
- **TransverseMercator** shows two bright discs on the equator, 90° east and west of the
  central meridian, where the projection is singular. That is
  [issue #40](https://github.com/JuliaEarth/CoordRefSystems.jl/issues/40).

## Running it

```julia
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. -e 'using Pluto; Pluto.run(notebook="notebook/atlas.jl")'
```

The environment is pinned, so the numbers above are reproducible. Re-running against a newer
release of CoordRefSystems.jl shows what has changed.
