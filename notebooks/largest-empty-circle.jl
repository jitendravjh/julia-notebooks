### A Pluto.jl notebook ###
# v0.20.4

#> [frontmatter]
#> title = "The largest empty circle"
#> description = "The point in Australia farthest from any town, reduced from an optimisation over the plane to a search over Voronoi vertices."
#> tags = ["julia", "geometry", "geospatial"]

using Markdown
using InteractiveUtils

# ╔═╡ 0dda0a73-ef84-426b-a1d4-c2737c568dab
md"""
# The largest empty circle

Given a finite set of sites in the plane, which point is farthest from all of them? Stated that
way it is an optimisation over a continuum, but it is not one. The largest circle containing no
site, with centre constrained to the convex hull of the sites, is centred either at a vertex of
their Voronoi diagram or on the hull boundary [1, ch. 7].

The reason is short. Grow a circle about a candidate centre until it meets a site. If it touches
one site, the centre can move directly away from it and the circle grows; if it touches two, the
centre can move along their bisector; only when it touches three is it stuck. Three equidistant
nearest sites is precisely the condition defining a Voronoi vertex.

That reduces the search from a continuum to $O(n)$ candidates. Below it is applied to towns in
Australia.
"""

# ╔═╡ f63d16d8-1eeb-4761-a587-2f2be59a20cc
begin
    using GeoArtifacts, Meshes, CoordRefSystems
    using CairoMakie, DataFrames, Unitful
    using Unitful: °
    using GeometryBasics: Polygon
end

# ╔═╡ 416619f4-058e-42db-b4d0-1ce9dc6ab5f7
md"""
## Sites

Natural Earth's populated places [2], restricted to Australia. This is a selected gazetteer
rather than a complete one, so throughout, *town* means a place Natural Earth lists.
"""

# ╔═╡ ac03cdf0-a561-46fa-8fc5-1ab97636c6e4
places = NaturalEarth.populatedplaces()

# ╔═╡ 11c17343-eb56-48ee-9fa2-6b69504fc054
towns = let sel = findall(x -> !ismissing(x) && x == "Australia", places.ADM0NAME)
    DataFrame(
        name = String.(places.NAME[sel]),
        lat  = [ustrip(coords(g).lat) for g in places.geometry[sel]],
        lon  = [ustrip(coords(g).lon) for g in places.geometry[sel]],
    )
end

# ╔═╡ 8ad0ddc1-a88d-48d7-8843-e752195d4486
md"""
Voronoi tessellation is a planar construction, so the sites are projected first. Lambert
azimuthal equal-area, centred on the continent, keeps the distortion across Australia small.
"""

# ╔═╡ 1bbb2c4d-42b2-4e2f-ad3c-94db0547fd1d
CRS = CoordRefSystems.shift(LambertAzimuthal{-25.0°, WGS84Latest}, lonₒ = 134.0°)

# ╔═╡ 7411f4d6-45cd-429e-ac64-a6722cb9584c
project(lat, lon) = let p = convert(CRS, LatLon(lat, lon))
    (ustrip(p.x), ustrip(p.y))
end

# ╔═╡ e3a91009-90e3-40b3-939b-28e59f7b8403
unproject(x, y) = let g = convert(LatLon, CRS(x * u"m", y * u"m"))
    (ustrip(g.lat), ustrip(g.lon))
end

# ╔═╡ 29eac646-7ac5-40e2-9796-16e3c88a0aa4
md"""
## The diagram

Each cell is the set of points closer to its own site than to any other.
"""

# ╔═╡ 04d92bba-70a8-4dc0-95d1-5b15235bbf27
sites = PointSet([Meshes.Point(project(r.lat, r.lon)...) for r in eachrow(towns)])

# ╔═╡ 6a61508e-fcd9-4c9c-8059-be4e26aac578
cells = tesselate(sites, VoronoiTesselation())

# ╔═╡ 9e25fb7f-fdcf-4708-b04c-c045999165e3
coast = let ct = NaturalEarth.countries()
    g = ct.geometry[findfirst(==("Australia"), ct.NAME)]
    mainland(lat, lon) = -45 ≤ lat ≤ -9 && 110 ≤ lon ≤ 155
    rs = [[(ustrip(coords(p).lat), ustrip(coords(p).lon)) for p in vertices(r)]
          for poly in parent(g) for r in Meshes.rings(poly)]
    rs = filter(r -> length(r) ≥ 4 && all(q -> mainland(q...), r), rs)
    [[project(a, b) for (a, b) in r] for r in rs]
end

# ╔═╡ a5284318-5553-4b25-9d02-e9990d13ece2
let
    fig = Figure(size = (1000, 860), backgroundcolor = PAPER)
    ax = plainaxis(fig)
    drawcoast!(ax)
    viz!(ax, cells, showsegments = true, segmentcolor = CELL, color = RGBAf(0, 0, 0, 0))
    drawtowns!(ax)
    frame!(ax)
    fig
end

# ╔═╡ ac472fbe-8253-4566-a523-918453e49d05
md"""
Cell size is a direct reading of how sparse the sites are: small along the settled coast, large
across the interior.

## Candidate centres

The Voronoi vertices, restricted to those lying on land, since a centre in the Southern Ocean
answers a different question.
"""

# ╔═╡ 12b9a8f4-1685-4022-8685-1296f827ed8d
corners = [(ustrip(coords(p).x), ustrip(coords(p).y)) for p in vertices(cells)]

# ╔═╡ 31685a46-1350-4f74-bbc6-3ff4a4387220
onland = let polys = [PolyArea(Ring(Meshes.Point.(r)...)) for r in coast]
    (x, y) -> any(pa -> Meshes.Point(x, y) ∈ pa, polys)
end

# ╔═╡ 0ac0b54a-3475-4096-89c1-fc239f6beb80
md"""
Distances are then measured on the sphere rather than in the projection, so the reported figure
is a real distance and not a projected one.
"""

# ╔═╡ 0b1af604-d396-44d1-bea2-973927d9c605
function nearesttown(lat, lon)
    d, i = findmin(greatcircle(lat, lon, r.lat, r.lon) for r in eachrow(towns))
    (distance = d, town = towns.name[i])
end

# ╔═╡ e5984119-30bc-4b19-b492-ac544ded3f39
remote = let candidates = [c for c in corners if onland(c...)]
    scored = map(candidates) do (x, y)
        lat, lon = unproject(x, y)
        (lat = lat, lon = lon, x = x, y = y, nearesttown(lat, lon)...)
    end
    argmax(s -> s.distance, scored)
end

# ╔═╡ 7ed5bbcc-fa70-41ac-82ea-27ab1d842945
DataFrame(
    candidates_on_land = count(c -> onland(c...), corners),
    latitude  = round(remote.lat, digits = 4),
    longitude = round(remote.lon, digits = 4),
    km_to_nearest = round(remote.distance / 1000, digits = 1),
)

# ╔═╡ f45cef38-a2f7-408b-9f6f-76617da7190e
md"""
The argument above predicts three sites tied at the minimum distance. They are.
"""

# ╔═╡ 7d4333b3-56a9-4e47-84ca-5d859fdc51f5
DataFrame(town = towns.name,
          km = round.([greatcircle(remote.lat, remote.lon, r.lat, r.lon) / 1000
                       for r in eachrow(towns)], digits = 1)) |>
    df -> first(sort(df, :km), 5)

# ╔═╡ 9704adf4-7117-4692-96fa-ebb398abc462
md"""
Three towns within about two kilometres of one another, then a clear gap to the fourth. The
circle is resting on all three.

They are not exactly equal because the diagram was built in the projection while the distances
were measured on the sphere. Two kilometres in five hundred and seventy-six is the magnitude of
that inconsistency, and it is the price of using a planar tessellation on a curved surface. A
spherical Voronoi construction would remove it.
"""

# ╔═╡ cb7fcb4a-95d6-40dc-a2c0-52689fb2e04b
let
    fig = Figure(size = (1000, 860), backgroundcolor = PAPER)
    ax = plainaxis(fig)
    drawcoast!(ax)
    viz!(ax, cells, showsegments = true, segmentcolor = CELLFAINT, color = RGBAf(0, 0, 0, 0))
    drawtowns!(ax)

    ring = smallcircle(remote.lat, remote.lon, remote.distance)
    lines!(ax, [Point2f(project(a, b)...) for (a, b) in ring], color = CIRCLE, linewidth = 2.5)
    scatter!(ax, [Point2f(remote.x, remote.y)], color = CIRCLE, markersize = 13,
             strokecolor = :black, strokewidth = 1)
    frame!(ax)
    fig
end

# ╔═╡ bcb7ca7e-f19e-4b2e-ae0e-e624defe1830
md"""
The centre sits where three cell boundaries meet, and the circle encloses no site.

The construction is not specific to towns. Any question of the form *what is farthest from all
of these* over a fixed finite set — clinics, transmitters, fuel stops, sensors — has the same
finite candidate set, and the Voronoi diagram produces it.
"""

# ╔═╡ d3a08ed1-6b07-47cf-b022-e98bdd822239
md"""
---
## Appendix
"""

# ╔═╡ 1d7657ec-d57d-4b0a-817b-c5d34edd7a85
begin
    const PAPER     = RGBf(0.99, 0.98, 0.96)
    const LANDFILL  = RGBf(0.91, 0.88, 0.80)
    const EDGE      = RGBf(0.45, 0.42, 0.38)
    const CELL      = RGBAf(0.25, 0.45, 0.72, 0.45)
    const CELLFAINT = RGBAf(0.25, 0.45, 0.72, 0.22)
    const TOWN      = RGBf(0.20, 0.25, 0.32)
    const CIRCLE    = RGBf(0.85, 0.30, 0.15)
end

# ╔═╡ 7ceb27a3-cfb4-48e4-a2ad-84d7ca98fde8
const EARTH = 6371008.8

# ╔═╡ e45718b8-237f-4101-b6cb-d3f1e5018e2c
"great-circle distance in metres on a sphere of mean Earth radius"
function greatcircle(lat₁, lon₁, lat₂, lon₂)
    φ₁, φ₂ = deg2rad(lat₁), deg2rad(lat₂)
    Δφ, Δλ = φ₂ - φ₁, deg2rad(lon₂ - lon₁)
    h = sin(Δφ / 2)^2 + cos(φ₁) * cos(φ₂) * sin(Δλ / 2)^2
    2EARTH * atan(sqrt(h), sqrt(max(1 - h, 0)))
end

# ╔═╡ 7662f964-70b3-4084-b77e-6110c6bde60f
"the locus of points at a fixed great-circle distance from a centre"
function smallcircle(lat₀, lon₀, metres; n = 240)
    ρ = metres / EARTH
    φ₀, λ₀ = deg2rad(lat₀), deg2rad(lon₀)
    map(range(0, 2π, length = n)) do α
        φ = asin(clamp(sin(φ₀) * cos(ρ) + cos(φ₀) * sin(ρ) * cos(α), -1, 1))
        λ = λ₀ + atan(sin(α) * sin(ρ) * cos(φ₀), cos(ρ) - sin(φ₀) * sin(φ))
        (rad2deg(φ), rad2deg(λ))
    end
end

# ╔═╡ 437d252c-c7ea-4a9b-97fa-9ec0fff8e5f3
function plainaxis(fig)
    ax = Axis(fig[1, 1], aspect = DataAspect(), backgroundcolor = PAPER)
    hidedecorations!(ax)
    hidespines!(ax)
    ax
end

# ╔═╡ 140a95d6-66fa-42cb-a38d-497092e450a6
drawcoast!(ax) = poly!(ax, [Polygon(Point2f.(r)) for r in coast],
                       color = LANDFILL, strokecolor = EDGE, strokewidth = 0.6)

# ╔═╡ a639cdf5-a050-4da5-b094-6aa4d2fce493
drawtowns!(ax) = scatter!(ax, [Point2f(project(r.lat, r.lon)...) for r in eachrow(towns)],
                          color = TOWN, markersize = 4)

# ╔═╡ ee9cdbfb-8ee5-4b0d-88d6-86ab78fbea10
function frame!(ax; pad = 2.0e5)
    xs = [p[1] for r in coast for p in r]
    ys = [p[2] for r in coast for p in r]
    limits!(ax, minimum(xs) - pad, maximum(xs) + pad, minimum(ys) - pad, maximum(ys) + pad)
end

# ╔═╡ 4e6008ce-f257-48c4-a0cd-a73445b3c406
md"""
---
## References and data

[1] M. de Berg, O. Cheong, M. van Kreveld and M. Overmars, *Computational Geometry: Algorithms
and Applications*, 3rd ed., Springer, 2008. Chapter 7 treats Voronoi diagrams; the largest
empty circle is a standard application.

[2] [Natural Earth](https://www.naturalearthdata.com), populated places and country boundaries,
public domain, loaded through
[GeoArtifacts.jl](https://github.com/JuliaEarth/GeoArtifacts.jl).

Tessellation from [Meshes.jl](https://github.com/JuliaGeometry/Meshes.jl), projection from
[CoordRefSystems.jl](https://github.com/JuliaEarth/CoordRefSystems.jl), figures with
[Makie.jl](https://docs.makie.org).
"""

# ╔═╡ Cell order:
# ╟─0dda0a73-ef84-426b-a1d4-c2737c568dab
# ╠═f63d16d8-1eeb-4761-a587-2f2be59a20cc
# ╟─416619f4-058e-42db-b4d0-1ce9dc6ab5f7
# ╠═ac03cdf0-a561-46fa-8fc5-1ab97636c6e4
# ╠═11c17343-eb56-48ee-9fa2-6b69504fc054
# ╟─8ad0ddc1-a88d-48d7-8843-e752195d4486
# ╠═1bbb2c4d-42b2-4e2f-ad3c-94db0547fd1d
# ╠═7411f4d6-45cd-429e-ac64-a6722cb9584c
# ╠═e3a91009-90e3-40b3-939b-28e59f7b8403
# ╟─29eac646-7ac5-40e2-9796-16e3c88a0aa4
# ╠═04d92bba-70a8-4dc0-95d1-5b15235bbf27
# ╠═6a61508e-fcd9-4c9c-8059-be4e26aac578
# ╠═9e25fb7f-fdcf-4708-b04c-c045999165e3
# ╠═a5284318-5553-4b25-9d02-e9990d13ece2
# ╟─ac472fbe-8253-4566-a523-918453e49d05
# ╠═12b9a8f4-1685-4022-8685-1296f827ed8d
# ╠═31685a46-1350-4f74-bbc6-3ff4a4387220
# ╟─0ac0b54a-3475-4096-89c1-fc239f6beb80
# ╠═0b1af604-d396-44d1-bea2-973927d9c605
# ╠═e5984119-30bc-4b19-b492-ac544ded3f39
# ╠═7ed5bbcc-fa70-41ac-82ea-27ab1d842945
# ╟─f45cef38-a2f7-408b-9f6f-76617da7190e
# ╠═7d4333b3-56a9-4e47-84ca-5d859fdc51f5
# ╟─9704adf4-7117-4692-96fa-ebb398abc462
# ╠═cb7fcb4a-95d6-40dc-a2c0-52689fb2e04b
# ╟─bcb7ca7e-f19e-4b2e-ae0e-e624defe1830
# ╟─d3a08ed1-6b07-47cf-b022-e98bdd822239
# ╠═1d7657ec-d57d-4b0a-817b-c5d34edd7a85
# ╠═7ceb27a3-cfb4-48e4-a2ad-84d7ca98fde8
# ╠═e45718b8-237f-4101-b6cb-d3f1e5018e2c
# ╠═7662f964-70b3-4084-b77e-6110c6bde60f
# ╠═437d252c-c7ea-4a9b-97fa-9ec0fff8e5f3
# ╠═140a95d6-66fa-42cb-a38d-497092e450a6
# ╠═a639cdf5-a050-4da5-b094-6aa4d2fce493
# ╠═ee9cdbfb-8ee5-4b0d-88d6-86ab78fbea10
# ╟─4e6008ce-f257-48c4-a0cd-a73445b3c406
