### A Pluto.jl notebook ###
# v0.20.4

#> [frontmatter]
#> title = "A round-trip atlas for map projections"
#> description = "Projecting every point on the globe and converting it back, to measure where each projection in CoordRefSystems.jl loses accuracy."
#> tags = ["julia", "geospatial", "cartography"]

using Markdown
using InteractiveUtils

# ╔═╡ 30ed3460-d80b-41ca-877e-1f2a44f42215
md"""
# A round-trip atlas for map projections

Every map projection in [CoordRefSystems.jl](https://github.com/JuliaEarth/CoordRefSystems.jl)
makes a quiet promise: take a latitude and longitude, project it, convert it back, and you
should land where you started.

This notebook measures how well that promise holds, for every projection, at every point on
the globe. The answer is not uniform. Most projections are exact to the nanometre almost
everywhere. A few are off by tens or hundreds of metres in places, and **where** they go wrong
turns out to say something precise about how they are implemented.
"""

# ╔═╡ 8f5e4671-47d4-467d-8795-49b2a2d2cdb4
begin
    using CoordRefSystems
    using CairoMakie
    using Unitful
    using Unitful: °
end

# ╔═╡ 189825a6-0243-423b-91de-5a7f3454da10
md"""
## Measuring the error honestly

The obvious metric is misleading. Comparing latitude and longitude componentwise makes every
projection look catastrophic at the poles, because longitude is undefined there: a point at the
north pole can come back with any longitude at all and still be the same place.

So the atlas measures **great-circle displacement** instead, in metres: how far the recovered
point actually landed from where it started. Longitude ambiguity at a pole then contributes
nothing, which is correct, and what remains is real error.
"""

# ╔═╡ 69379198-97e3-4396-b15f-801d0c2e9d3d
const RADIUS = 6371008.8  # mean Earth radius, metres

# ╔═╡ 06cf0570-b012-4888-8b3b-eae6f4b955bc
function gcdist(lat₁, lon₁, lat₂, lon₂)
    φ₁, φ₂ = deg2rad(lat₁), deg2rad(lat₂)
    Δφ = φ₂ - φ₁
    Δλ = deg2rad(lon₂ - lon₁)
    a = sin(Δφ / 2)^2 + cos(φ₁) * cos(φ₂) * sin(Δλ / 2)^2
    2RADIUS * atan(sqrt(a), sqrt(max(1 - a, 0)))
end

# ╔═╡ b311d58a-e285-4da7-bdc7-ef0d6ca58950
md"""
Two outcomes deserve to be distinguished from a merely inaccurate answer: a point outside the
projection domain, where the question does not apply, and a conversion that throws.
"""

# ╔═╡ 2b68888b-5690-4196-bc0a-ace96dbe7047
function rterror(CRS, lat, lon)
    ll = LatLon(lat, lon)
    try
        indomain(CRS, ll) || return NaN   # outside the domain
        ll′ = convert(LatLon, convert(CRS, ll))
        gcdist(lat, lon, ustrip(ll′.lat), ustrip(ll′.lon))
    catch
        Inf                               # the conversion threw
    end
end

# ╔═╡ 3cae9ae7-ee02-4882-82f7-69def41d19cc
function atlas(CRS; nlat=181, nlon=361)
    lats = range(-90, 90, length=nlat)
    lons = range(-180, 180, length=nlon)
    lons, lats, [rterror(CRS, lat, lon) for lon in lons, lat in lats]
end

# ╔═╡ 0b6a1ff7-cdbb-46e9-a71a-7b5118197a78
function atlasplot!(ax, CRS; nlat=181, nlon=361, crange=(-10, 3))
    lons, lats, E = atlas(CRS; nlat, nlon)
    L = map(E) do e
        isnan(e) && return NaN
        isinf(e) && return crange[2]
        e <= 0 ? crange[1] : clamp(log10(e), crange...)
    end
    heatmap!(ax, lons, lats, L, colormap=:inferno, colorrange=crange, nan_color=:gray85)
end

# ╔═╡ 6eb7fb0a-b472-4f54-9f9c-8c573797fc06
function atlasgrid(entries; nlat=121, nlon=241, crange=(-10, 3), ncol=2)
    nrow = cld(length(entries), ncol)
    fig = Figure(size=(1000, 250nrow))
    hm = nothing
    for (k, (name, CRS)) in enumerate(entries)
        i, j = fldmod1(k, ncol)
        ax = Axis(fig[i, j], title=name, xticks=-180:90:180, yticks=-90:45:90)
        hm = atlasplot!(ax, CRS; nlat, nlon, crange)
    end
    Colorbar(fig[1:nrow, ncol+1], hm, label="log₁₀ round-trip displacement (m)")
    fig
end

# ╔═╡ 920ca5fd-648a-4e66-b116-bfb6f634d616
md"""
## The atlas

Dark is good: black means the round trip is exact to within a nanometre. Bright is bad.
Grey is outside the projection domain.
"""

# ╔═╡ 54456cb4-90f5-4c07-969e-dc116dbeccfc
projections = [
    "Mercator"           => Mercator,
    "WebMercator"        => WebMercator,
    "PlateCarree"        => PlateCarree,
    "Robinson"           => Robinson,
    "WinkelTripel"       => WinkelTripel,
    "Sinusoidal"         => Sinusoidal,
    "EqualEarth"         => EqualEarth,
    "GallPeters"         => GallPeters,
    "OrthoNorth"         => OrthoNorth,
    "LambertAzimuthal"   => LambertAzimuthal{15.0°},
    "Albers"             => Albers{23.0°,29.5°,45.0°},
    "TransverseMercator" => TransverseMercator{0.9996,0.0°},
]

# ╔═╡ 30528655-4570-4b45-a489-4b1e5b529420
atlasgrid(projections)

# ╔═╡ 0983dfe9-671a-48d9-a3f5-98200bbc233d
md"""
## Robinson: you can see the lookup table

Robinson is not defined by a formula. It is defined by a **table** of values at every five
degrees of latitude, interpolated in between. The atlas shows this directly: the error is
small everywhere except on a set of horizontal lines, and those lines sit on the tabulation
knots.

The cause was diagnosed in [issue #55](https://github.com/JuliaEarth/CoordRefSystems.jl/issues/55):
the reference implementation stores its interpolation coefficients as `Float32`, so in
`Float64` arithmetic the interpolant jumps slightly as you cross from one interval to the next.
"""

# ╔═╡ 1812d260-085d-41d6-9589-96d050ee2529
let
    _, lats, E = atlas(Robinson, nlat=361, nlon=73)
    worst = [maximum(filter(isfinite, E[:, j])) for j in eachindex(lats)]
    p = sortperm(worst, rev=true)[1:10]
    [(lat = round(lats[j], digits=2), worst_m = round(worst[j], sigdigits=4)) for j in p]
end

# ╔═╡ e7843334-a454-41f0-9290-d081a7cca1b5
md"""
Every one of the worst latitudes sits within half a degree of a multiple of five. That is the
table talking.

## TransverseMercator: a singularity you can point at

The two bright discs at the equator, ninety degrees east and west of the central meridian, are
where the transverse Mercator projection is singular. That region is the subject of
[issue #40](https://github.com/JuliaEarth/CoordRefSystems.jl/issues/40).
"""

# ╔═╡ 64b8fb10-6b20-4190-a1e2-f364d4c25d6b
atlasgrid(["TransverseMercator" => TransverseMercator{0.9996,0.0°}], ncol=1, nlat=241, nlon=481)

# ╔═╡ 9edd7ea1-ea55-435e-8031-b4b7c65c2f58
md"""
## Summary

Worst and mean round-trip displacement over the globe, in metres.
"""

# ╔═╡ 7fa883ea-58e6-49bc-abe6-dc874e869e27
let
    rows = map(projections) do (name, CRS)
        _, _, E = atlas(CRS, nlat=91, nlon=181)
        fin = filter(isfinite, E)
        (crs = name,
         worst_m = round(maximum(fin), sigdigits=3),
         mean_m = round(sum(fin) / length(fin), sigdigits=3),
         outside = count(isnan, E),
         threw = count(isinf, E))
    end
    sort(rows, by = r -> -r.worst_m)
end

# ╔═╡ 0af78af4-32ab-40d5-8aa7-6452a69f8c8b
md"""
## What this is for

Round-trip accuracy is easy to assume and rarely checked. Rendering it as a map makes two
things visible at once: which projections have a problem, and where the problem lives. A
horizontal band means a lookup table. A disc means a singularity. A single bright line at a
pole means an arcsine running out of precision.

The notebook is reproducible: it pins its dependencies and can be re-run against any release
of CoordRefSystems.jl to see what has changed.
"""

# ╔═╡ Cell order:
# ╟─30ed3460-d80b-41ca-877e-1f2a44f42215
# ╠═8f5e4671-47d4-467d-8795-49b2a2d2cdb4
# ╟─189825a6-0243-423b-91de-5a7f3454da10
# ╠═69379198-97e3-4396-b15f-801d0c2e9d3d
# ╠═06cf0570-b012-4888-8b3b-eae6f4b955bc
# ╟─b311d58a-e285-4da7-bdc7-ef0d6ca58950
# ╠═2b68888b-5690-4196-bc0a-ace96dbe7047
# ╠═3cae9ae7-ee02-4882-82f7-69def41d19cc
# ╠═0b6a1ff7-cdbb-46e9-a71a-7b5118197a78
# ╠═6eb7fb0a-b472-4f54-9f9c-8c573797fc06
# ╟─920ca5fd-648a-4e66-b116-bfb6f634d616
# ╠═54456cb4-90f5-4c07-969e-dc116dbeccfc
# ╠═30528655-4570-4b45-a489-4b1e5b529420
# ╟─0983dfe9-671a-48d9-a3f5-98200bbc233d
# ╠═1812d260-085d-41d6-9589-96d050ee2529
# ╟─e7843334-a454-41f0-9290-d081a7cca1b5
# ╠═64b8fb10-6b20-4190-a1e2-f364d4c25d6b
# ╟─9edd7ea1-ea55-435e-8031-b4b7c65c2f58
# ╠═7fa883ea-58e6-49bc-abe6-dc874e869e27
# ╟─0af78af4-32ab-40d5-8aa7-6452a69f8c8b
