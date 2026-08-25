### A Pluto.jl notebook ###
# v0.20.4

#> [frontmatter]
#> title = "The middle of nowhere"
#> description = "The point in Australia farthest from any town, found with a Voronoi diagram instead of a search."
#> tags = ["julia", "geometry", "geospatial"]

using Markdown
using InteractiveUtils

# ╔═╡ 34588ddc-a4ea-4cc9-8f9c-80182628e993
md"""
# The middle of nowhere

Somewhere in Australia there is a point on land that is farther from any town than any other
point on land. Finding it sounds like it needs a search over the whole continent. It does not.

The answer is always a corner of a Voronoi diagram.
"""

# ╔═╡ afb5b070-8aa3-443b-adb1-95b8ae304d15
begin
    using GeoArtifacts, Meshes, CoordRefSystems
    using CairoMakie, DataFrames, Unitful
    using Unitful: °
    using GeometryBasics: Polygon
end

# ╔═╡ bb9b3781-8927-48d8-8031-49456f33816c
md"""
## The towns

Natural Earth's populated places, restricted to Australia.
"""

# ╔═╡ 4b09519c-c8be-454b-b3c3-76a0b634d6d8
places = NaturalEarth.populatedplaces()

# ╔═╡ 2271370a-d1aa-4916-af94-56f2609378aa
towns = let sel = findall(x -> !ismissing(x) && x == "Australia", places.ADM0NAME)
    DataFrame(
        name = String.(places.NAME[sel]),
        lat  = [ustrip(coords(g).lat) for g in places.geometry[sel]],
        lon  = [ustrip(coords(g).lon) for g in places.geometry[sel]],
    )
end

# ╔═╡ 926a9351-b555-4fc3-8fa4-d59f74a08dba
md"""
Voronoi diagrams are planar, so the towns need projecting first. Lambert azimuthal equal-area,
centred on the continent.
"""

# ╔═╡ 6c534f8b-7628-4500-b595-ec0b7e1221fb
CRS = CoordRefSystems.shift(LambertAzimuthal{-25.0°, WGS84Latest}, lonₒ = 134.0°)

# ╔═╡ fe76bf59-30ff-467b-89dd-68e190c0ff85
project(lat, lon) = let p = convert(CRS, LatLon(lat, lon))
    (ustrip(p.x), ustrip(p.y))
end

# ╔═╡ 53ad0950-2fe3-4f57-885a-565b61bb07bd
unproject(x, y) = let g = convert(LatLon, CRS(x * u"m", y * u"m"))
    (ustrip(g.lat), ustrip(g.lon))
end

# ╔═╡ 7ebd9281-05a0-4be2-a69e-a9e264cf9bfe
md"""
## The diagram

Each cell is the set of points closer to that town than to any other.
"""

# ╔═╡ bb394935-3e3e-4e5b-a5fd-d54c602e0367
sites = PointSet([Meshes.Point(project(r.lat, r.lon)...) for r in eachrow(towns)])

# ╔═╡ f4c28be3-e031-4210-a082-71af2790e7ff
cells = tesselate(sites, VoronoiTesselation())

# ╔═╡ 509c69ff-04f0-49bf-8226-b7e8699d8b48
coast = let ct = NaturalEarth.countries()
    g = ct.geometry[findfirst(==("Australia"), ct.NAME)]
    mainland(lat, lon) = -45 ≤ lat ≤ -9 && 110 ≤ lon ≤ 155
    rs = [[(ustrip(coords(p).lat), ustrip(coords(p).lon)) for p in vertices(r)]
          for poly in parent(g) for r in Meshes.rings(poly)]
    rs = filter(r -> length(r) >= 4 && all(q -> mainland(q...), r), rs)
    [[project(a, b) for (a, b) in r] for r in rs]
end

# ╔═╡ 84e8b09d-b46d-4d22-b6bb-e7dabab5044e
let
    fig = Figure(size = (1000, 860), backgroundcolor = PAPER)
    ax = plainaxis(fig)
    drawcoast!(ax)
    viz!(ax, cells, showsegments = true, segmentcolor = CELL, color = RGBAf(0, 0, 0, 0))
    scatter!(ax, [Point2f(project(r.lat, r.lon)...) for r in eachrow(towns)],
             color = TOWN, markersize = 4)
    frameaustralia!(ax)
    fig
end

# ╔═╡ 292c904e-bfcb-415d-a3a4-1750ee33d72f
md"""
## Where the empty circle is largest

Grow a circle from any point until it touches a town. The largest such circle over the whole
continent has its centre at a Voronoi vertex, because only there can it touch three towns at
once, and touching fewer means you can still push it further.

So the search is over a few hundred corners, not the whole map.
"""

# ╔═╡ 0db4c88a-7fbe-48aa-b734-017873a58796
corners = [(ustrip(coords(p).x), ustrip(coords(p).y)) for p in vertices(cells)]

# ╔═╡ 50f12a18-d845-4099-963b-fd1138b6c11f
md"""
Keeping only the corners that fall on land, and measuring the real distance on the ellipsoid
rather than the distance in the projection.
"""

# ╔═╡ 3322bf44-b416-448b-8edc-7c3f7e089f74
onland = let polys = [PolyArea(Ring(Meshes.Point.(r)...)) for r in coast]
    (x, y) -> any(pa -> Meshes.Point(x, y) ∈ pa, polys)
end

# ╔═╡ 0f6d5849-0333-44c9-9d9f-975f36eec645
function nearesttown(lat, lon)
    d, i = findmin(greatcircle(lat, lon, r.lat, r.lon) for r in eachrow(towns))
    (distance = d, town = towns.name[i])
end

# ╔═╡ 65473e0f-4324-4339-9080-24c40ceafcb0
remote = let cands = [c for c in corners if onland(c...)]
    scored = map(cands) do (x, y)
        lat, lon = unproject(x, y)
        (lat = lat, lon = lon, x = x, y = y, nearesttown(lat, lon)...)
    end
    argmax(s -> s.distance, scored)
end

# ╔═╡ a6efdd0b-ca74-4172-9a92-b52876579b06
md"""
The point, and how far it is from anywhere.
"""

# ╔═╡ a98d0f36-0e0e-4d52-9e17-635e4959a7b0
DataFrame(
    latitude  = round(remote.lat, digits = 4),
    longitude = round(remote.lon, digits = 4),
    km_to_nearest_town = round(remote.distance / 1000, digits = 1),
    nearest = remote.town,
)

# ╔═╡ e8f592fd-101a-4aa1-9486-01bea3d9e9d3
md"""
The theorem says the circle should touch three towns. It does.
"""

# ╔═╡ b80ef110-570a-4b89-b26c-0753ad0076c2
DataFrame(
    town = towns.name,
    km   = round.([greatcircle(remote.lat, remote.lon, r.lat, r.lon) / 1000
                   for r in eachrow(towns)], digits = 1),
) |> df -> sort(df, :km)[1:5, :]

# ╔═╡ 62e183b8-8044-4e62-b815-9e5e5d1a16af
md"""
Three towns within two kilometres of each other in distance, and then a gap. The circle is
resting on all three.

They do not come out exactly equal because the diagram was built in a projection while the
distances were measured on the ellipsoid. Two kilometres in five hundred and seventy-six is the
size of that disagreement.
"""

# ╔═╡ 71e35ee7-701c-481e-8e09-996521533fef
let
    fig = Figure(size = (1000, 860), backgroundcolor = PAPER)
    ax = plainaxis(fig)
    drawcoast!(ax)
    viz!(ax, cells, showsegments = true, segmentcolor = CELLFAINT, color = RGBAf(0, 0, 0, 0))
    scatter!(ax, [Point2f(project(r.lat, r.lon)...) for r in eachrow(towns)],
             color = TOWN, markersize = 4)

    ring = smallcircle(remote.lat, remote.lon, remote.distance)
    lines!(ax, [Point2f(project(a, b)...) for (a, b) in ring], color = CIRCLE, linewidth = 2.5)
    scatter!(ax, [Point2f(remote.x, remote.y)], color = CIRCLE, markersize = 13,
             strokecolor = :black, strokewidth = 1)
    frameaustralia!(ax)
    fig
end

# ╔═╡ 708fc37d-eb24-4c46-be32-1c52c42fd8a1
md"""
Nothing inside the circle. Five hundred and seventy-six kilometres in every direction before
you reach a town, in the Great Victoria Desert.

The same three lines of code work for any country in the dataset, and for anything else you can
put on a map: hospitals, transmitters, fuel stops. The question "what is farthest from all of
these" has a finite answer set, and Voronoi hands it to you.
"""

# ╔═╡ b665ef59-3265-4799-b6d0-c4761cc1b3dd
md"""
---
## Appendix
"""

# ╔═╡ 0c3db225-f3b2-41e6-91d6-a7fc74e1f27f
begin
    const PAPER     = RGBf(0.99, 0.98, 0.96)
    const LAND      = RGBf(0.91, 0.88, 0.80)
    const EDGE      = RGBf(0.45, 0.42, 0.38)
    const CELL      = RGBAf(0.25, 0.45, 0.72, 0.45)
    const CELLFAINT = RGBAf(0.25, 0.45, 0.72, 0.20)
    const TOWN      = RGBf(0.20, 0.25, 0.32)
    const CIRCLE    = RGBf(0.85, 0.30, 0.15)
end

# ╔═╡ acc26733-429c-4bd4-b275-c98507b0f33e
const EARTH = 6371008.8

# ╔═╡ 342e69e9-3429-4df9-88ac-15b1f2c163a4
"great-circle distance in metres"
function greatcircle(lat₁, lon₁, lat₂, lon₂)
    φ₁, φ₂ = deg2rad(lat₁), deg2rad(lat₂)
    Δφ, Δλ = φ₂ - φ₁, deg2rad(lon₂ - lon₁)
    h = sin(Δφ / 2)^2 + cos(φ₁) * cos(φ₂) * sin(Δλ / 2)^2
    2EARTH * atan(sqrt(h), sqrt(max(1 - h, 0)))
end

# ╔═╡ 4eec2963-c612-4380-b607-9335ec445a5a
"the set of points at a fixed great-circle distance from a centre"
function smallcircle(lat₀, lon₀, metres; n = 240)
    ρ = metres / EARTH
    φ₀, λ₀ = deg2rad(lat₀), deg2rad(lon₀)
    map(range(0, 2π, length = n)) do α
        φ = asin(clamp(sin(φ₀) * cos(ρ) + cos(φ₀) * sin(ρ) * cos(α), -1, 1))
        λ = λ₀ + atan(sin(α) * sin(ρ) * cos(φ₀), cos(ρ) - sin(φ₀) * sin(φ))
        (rad2deg(φ), rad2deg(λ))
    end
end

# ╔═╡ 54bbd347-220d-4e83-9b10-3ca87ccc5a17
function plainaxis(fig)
    ax = Axis(fig[1, 1], aspect = DataAspect(), backgroundcolor = PAPER)
    hidedecorations!(ax)
    hidespines!(ax)
    ax
end

# ╔═╡ a0fbce1d-a849-4e74-aa3d-c980ec35ee0d
drawcoast!(ax) = poly!(ax, [Polygon(Point2f.(r)) for r in coast],
                       color = LAND, strokecolor = EDGE, strokewidth = 0.6)

# ╔═╡ 0d9646a8-d2e8-49f9-8178-470d0419c144
function frameaustralia!(ax)
    xs = [p[1] for r in coast for p in r]
    ys = [p[2] for r in coast for p in r]
    pad = 2.0e5
    limits!(ax, minimum(xs) - pad, maximum(xs) + pad, minimum(ys) - pad, maximum(ys) + pad)
end

# ╔═╡ f9e413b1-1541-4377-b291-5c9ace538f0f
md"""
---
## References and data

Populated places and country boundaries from [Natural Earth](https://www.naturalearthdata.com)
(public domain), loaded with
[GeoArtifacts.jl](https://github.com/JuliaEarth/GeoArtifacts.jl). Natural Earth's populated
places are a selected set, not a complete gazetteer, so "town" here means "town Natural Earth
lists".

The result that the largest empty circle is centred on a Voronoi vertex, on a site, or on the
convex hull boundary is standard computational geometry; see M. de Berg, O. Cheong, M. van
Kreveld and M. Overmars, *Computational Geometry: Algorithms and Applications*, 3rd ed.,
Springer, 2008, chapter 7.

Voronoi tessellation from [Meshes.jl](https://github.com/JuliaGeometry/Meshes.jl), projection
from [CoordRefSystems.jl](https://github.com/JuliaEarth/CoordRefSystems.jl), drawing with
[Makie.jl](https://docs.makie.org).
"""

# ╔═╡ Cell order:
# ╟─34588ddc-a4ea-4cc9-8f9c-80182628e993
# ╠═afb5b070-8aa3-443b-adb1-95b8ae304d15
# ╟─bb9b3781-8927-48d8-8031-49456f33816c
# ╠═4b09519c-c8be-454b-b3c3-76a0b634d6d8
# ╠═2271370a-d1aa-4916-af94-56f2609378aa
# ╟─926a9351-b555-4fc3-8fa4-d59f74a08dba
# ╠═6c534f8b-7628-4500-b595-ec0b7e1221fb
# ╠═fe76bf59-30ff-467b-89dd-68e190c0ff85
# ╠═53ad0950-2fe3-4f57-885a-565b61bb07bd
# ╟─7ebd9281-05a0-4be2-a69e-a9e264cf9bfe
# ╠═bb394935-3e3e-4e5b-a5fd-d54c602e0367
# ╠═f4c28be3-e031-4210-a082-71af2790e7ff
# ╠═509c69ff-04f0-49bf-8226-b7e8699d8b48
# ╠═84e8b09d-b46d-4d22-b6bb-e7dabab5044e
# ╟─292c904e-bfcb-415d-a3a4-1750ee33d72f
# ╠═0db4c88a-7fbe-48aa-b734-017873a58796
# ╟─50f12a18-d845-4099-963b-fd1138b6c11f
# ╠═3322bf44-b416-448b-8edc-7c3f7e089f74
# ╠═0f6d5849-0333-44c9-9d9f-975f36eec645
# ╠═65473e0f-4324-4339-9080-24c40ceafcb0
# ╟─a6efdd0b-ca74-4172-9a92-b52876579b06
# ╠═a98d0f36-0e0e-4d52-9e17-635e4959a7b0
# ╟─e8f592fd-101a-4aa1-9486-01bea3d9e9d3
# ╠═b80ef110-570a-4b89-b26c-0753ad0076c2
# ╟─62e183b8-8044-4e62-b815-9e5e5d1a16af
# ╠═71e35ee7-701c-481e-8e09-996521533fef
# ╟─708fc37d-eb24-4c46-be32-1c52c42fd8a1
# ╟─b665ef59-3265-4799-b6d0-c4761cc1b3dd
# ╠═0c3db225-f3b2-41e6-91d6-a7fc74e1f27f
# ╠═acc26733-429c-4bd4-b275-c98507b0f33e
# ╠═342e69e9-3429-4df9-88ac-15b1f2c163a4
# ╠═4eec2963-c612-4380-b607-9335ec445a5a
# ╠═54bbd347-220d-4e83-9b10-3ca87ccc5a17
# ╠═a0fbce1d-a849-4e74-aa3d-c980ec35ee0d
# ╠═0d9646a8-d2e8-49f9-8178-470d0419c144
# ╟─f9e413b1-1541-4377-b291-5c9ace538f0f
