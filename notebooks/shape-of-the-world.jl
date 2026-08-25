### A Pluto.jl notebook ###
# v0.20.4

#> [frontmatter]
#> title = "The shape of the world"
#> description = "Morphing real coastlines between map projections, and using Tissot's indicatrix to show exactly what each one gives up."
#> tags = ["julia", "geospatial", "cartography"]

using Markdown
using InteractiveUtils

# ╔═╡ b3d68a38-bb27-40e2-a529-5b4304a320bb
md"""
# The shape of the world

Every flat map of a round planet is a lie. The only question is which lie you prefer.

You cannot flatten a sphere without tearing or stretching it. That is Gauss's *Theorema
Egregium*, and it is not a limitation of our cleverness but a fact about curved surfaces. So
every projection has to give something up. Some preserve angles. Some preserve areas. None
preserve both, and none ever will.

This notebook makes the tradeoff visible, using real coastlines and the projections in
[CoordRefSystems.jl](https://github.com/JuliaEarth/CoordRefSystems.jl).
"""

# ╔═╡ 3da10fc9-0e94-4426-a8ca-fae3ff10d583
begin
    using GeoArtifacts, Meshes, CoordRefSystems
    using CairoMakie, Unitful, Base64
    using Unitful: °
    using GeometryBasics: Polygon
end

# ╔═╡ 00226b63-87b7-472d-a7a5-73f0f4df95bc
md"""
## The data

Natural Earth country boundaries, straight from `GeoArtifacts`: every country as a multipolygon
on the `🌐` manifold, in geodetic latitude and longitude. We keep the rings as plain `(lat, lon)`
tuples so we can push them through any projection we like.
"""

# ╔═╡ efb0cafc-e855-456f-bc01-b895b994cf35
land = filter(r -> length(r) >= 8, decimate(worldrings(), 6));

# ╔═╡ 23b3f9e1-ca54-49d9-af2d-5568e7b335b3
(rings = length(land), points = sum(length, land))

# ╔═╡ ef070602-d668-42ba-8125-8e5d7bcdf74b
begin
    grid    = graticule()
    frame   = [frameline()]
    circles = tissot()
end;

# ╔═╡ 4f6b3239-e7a9-4be7-aa1e-4032a11f87b6
md"""
## Watching the world change its mind

Every point gets projected twice, once with each projection, and the two results are blended.
The world morphs continuously from one way of being flat to another.

Watch the graticule rather than the coastlines. In Plate Carrée the meridians are parallel
straight lines. By Sinusoidal they have collapsed into curves meeting at the poles. The land
has no choice but to follow.
"""

# ╔═╡ 6ada4791-c3bd-4339-8047-a8879b2f11af
morphprojections = [
    Proj("Plate Carrée",  PlateCarree),
    Proj("Robinson",      Robinson),
    Proj("Winkel Tripel", WinkelTripel),
    Proj("Equal Earth",   EqualEarth),
    Proj("Sinusoidal",    Sinusoidal),
    Proj("Gall–Peters",   GallPeters),
]

# ╔═╡ 895f4784-cc99-4d04-93b5-448581375280
video(morphmovie(tempname() * ".mp4", morphprojections, land, grid, frame))

# ╔═╡ 07e85657-41b2-47f4-a5e1-e8e0c013d917
md"""
## Tissot's indicatrix

The morph shows *that* shapes change. It does not show *how much*, or in which direction. For
that there is a device from 1859.

Draw small circles of equal size on the globe, then project them along with everything else.
Whatever happens to a circle is exactly what the projection does to anything small at that
spot. A circle that stays circular means angles survive. A circle that keeps its area means
areas survive.
"""

# ╔═╡ 6db80fbd-3c1d-4f95-bed3-4cbd56b1f099
function tissotmap(P; title = P.name)
    fig = Figure(size = (1150, 660), backgroundcolor = OCEAN)
    ax = mapaxis(fig)
    ax.title = title
    ax.titlecolor = :white
    ax.titlesize = 20
    drawmap!(ax, project(P, land), project(P, grid), project(P, frame))
    drawtissot!(ax, project(P, circles))
    fig
end

# ╔═╡ b6cb65d5-452c-4e29-a530-24d1528be4bd
md"""
### Mercator preserves angles

Every circle stays a perfect circle, at every latitude. That is what *conformal* means, and it
is why you navigate with Mercator: a constant compass bearing is a straight line on the page.

The price is size. The circles grow without limit toward the poles. Greenland is not the size
of Africa. Mercator only makes it look that way.
"""

# ╔═╡ 26247572-b2d9-472c-9431-43fd3dbae91c
tissotmap(Proj("Mercator", Mercator, 82.0))

# ╔═╡ 870cf06b-7ff5-4270-a37d-8b72015c50b3
md"""
### Gall–Peters preserves areas

Now the circles are ellipses, so angles are wrecked. But look at the *area* of each one: it is
constant. Tall and narrow near the equator, short and wide near the poles, every one covering
the same amount of ink.

Same planet, same data, same code. The entire difference is in what was given up.
"""

# ╔═╡ f67bca02-ba46-4624-b0d6-fcb8ee3cfd98
tissotmap(Proj("Gall–Peters", GallPeters))

# ╔═╡ cbb59117-a527-409a-8002-d6cb2dfa849d
md"""
### Robinson preserves nothing, deliberately

Robinson is a compromise. The circles are neither perfectly round nor perfectly equal, because
it was designed by eye to *look* right rather than to satisfy a theorem. For decades it was
what National Geographic put on the wall.

Sometimes the honest answer to an impossible constraint is to miss both targets by a little.
"""

# ╔═╡ d6844ed7-6b74-47c1-aa4c-8a96952f5021
tissotmap(Proj("Robinson", Robinson))

# ╔═╡ e0b7639b-ded3-42f4-bd8d-4650105c4760
md"""
## The library already knows

None of this is folklore. CoordRefSystems.jl carries the properties as queryable traits, so you
can ask a projection what it promises before you trust it with your data.
"""

# ╔═╡ 7798f437-6d69-4805-a5be-94bbee12c175
let
    ps = [("Mercator", Mercator), ("Plate Carrée", PlateCarree), ("Robinson", Robinson),
          ("Winkel Tripel", WinkelTripel), ("Sinusoidal", Sinusoidal),
          ("Equal Earth", EqualEarth), ("Gall–Peters", GallPeters)]
    mark(b) = b ? "**yes**" : "no"
    rows = map(ps) do (name, T)
        "| $name | $(mark(CoordRefSystems.isconformal(T))) " *
        "| $(mark(CoordRefSystems.isequalarea(T))) " *
        "| $(mark(CoordRefSystems.isequidistant(T))) |"
    end
    Markdown.parse(join(["| projection | conformal | equal-area | equidistant |",
                         "|---|:--:|:--:|:--:|", rows...], "\n"))
end

# ╔═╡ ce8ecfb6-e1c4-42d5-96cd-661bc7af4841
md"""
No row says yes twice in the first two columns, and no row ever can. That is Gauss's theorem
showing up as a column of *no*.

Which one you should use depends entirely on the question you are asking. Navigating? Mercator.
Comparing how much land two countries cover? Anything equal-area. Putting a world map on a wall
and wanting it to look like the world? Robinson or Winkel Tripel, and accept that every number
you read off it is slightly wrong.
"""

# ╔═╡ 4646a8ed-68ec-4e3d-bd8c-1193d057cc85
md"""
---
## Appendix: the machinery

Everything above runs on the handful of functions below. Nothing here is specific to any one
projection, which is rather the point: `CoordRefSystems.convert` does the hard part, and the
rest is turning tuples into pictures.
"""

# ╔═╡ cc43fc0b-70ef-4568-b743-8ab2af2006db
md"""
### Reading the world
"""

# ╔═╡ 39a8f29c-e96b-49da-b0a6-c15e661aa3ad
"every country ring as a plain vector of (lat, lon)"
function worldrings(; minpts = 12)
    gt = NaturalEarth.countries()
    out = Vector{Vector{Tuple{Float64,Float64}}}()
    for geom in gt.geometry, poly in parent(geom), ring in Meshes.rings(poly)
        pts = [(ustrip(coords(p).lat), ustrip(coords(p).lon)) for p in vertices(ring)]
        length(pts) >= minpts && push!(out, pts)
    end
    out
end

# ╔═╡ 21b4a30e-8555-432e-a1cf-9fe0b34247f5
decimate(rs, stride) = [length(r) <= 20 ? r : r[1:stride:end] for r in rs]

# ╔═╡ 0079048a-cdc3-467e-9bd9-d5525b3de382
"meridians and parallels, as (lat, lon) polylines"
function graticule(; dlon = 30, dlat = 30, n = 181)
    g = Vector{Vector{Tuple{Float64,Float64}}}()
    for lon in -180:dlon:180
        push!(g, [(lat, float(lon)) for lat in range(-89.5, 89.5, length = n)])
    end
    for lat in -60:dlat:60
        push!(g, [(float(lat), lon) for lon in range(-180, 180, length = n)])
    end
    g
end

# ╔═╡ 36db0674-73aa-4bc0-9e00-97f591c954d9
"the outer boundary of the map: the antimeridian both ways, plus the poles"
function frameline(; n = 361)
    [[(lat, -179.99) for lat in range(-89.5, 89.5, length = n)];
     [(89.5, lon) for lon in range(-179.99, 179.99, length = n)];
     [(lat, 179.99) for lat in range(89.5, -89.5, length = n)];
     [(-89.5, lon) for lon in range(179.99, -179.99, length = n)]]
end

# ╔═╡ f3caa581-116e-4275-8d6c-3561915ebcba
md"""
### Projecting

`project` is the only place that touches CoordRefSystems. Results are normalised to a common
scale afterwards, so two projections with wildly different units can be blended and drawn in
the same axis.
"""

# ╔═╡ 1448c824-e90f-4e5d-9f55-25fc4a7e022d
begin
    "a projection, plus the latitude beyond which it blows up"
    struct Proj
        name::String
        crs::Any
        latclip::Float64
    end
    Proj(name, crs) = Proj(name, crs, 89.5)
end

# ╔═╡ 07ea1ca3-c1a7-4e7b-9c0d-a4449eaa64a6
function project(P::Proj, rs)
    out = map(rs) do pts
        map(pts) do (lat, lon)
            p = convert(P.crs, LatLon(clamp(lat, -P.latclip, P.latclip), lon))
            (ustrip(p.x), ustrip(p.y))
        end
    end
    s = maximum(abs(v) for r in out for xy in r for v in xy)
    [[(x / s, y / s) for (x, y) in r] for r in out]
end

# ╔═╡ cf08360c-58cb-4461-ba18-94d006815508
blend(A, B, t) = [[((1 - t) * a[1] + t * b[1], (1 - t) * a[2] + t * b[2])
                   for (a, b) in zip(ra, rb)] for (ra, rb) in zip(A, B)]

# ╔═╡ 84688259-9341-4506-bbd2-230134f917d1
md"""
### Drawing

Two thousand coastline rings is enough that one `poly!` per ring is far too slow to animate.
Batching them into a single `poly!` with a vector of polygons, and the whole graticule into one
`lines!` with `NaN` separators, takes a frame from thirty seconds to under a tenth of a second.
"""

# ╔═╡ ef944087-9b31-42e2-965d-6cae8d127eab
begin
    const LAND   = RGBf(0.93, 0.90, 0.83)
    const EDGE   = RGBf(0.35, 0.33, 0.30)
    const OCEAN  = RGBf(0.05, 0.09, 0.13)
    const GRID   = RGBf(0.22, 0.32, 0.40)
    const TISSOT = RGBf(0.95, 0.45, 0.25)
end

# ╔═╡ 8fb80132-9f6f-47e4-8c96-2342c5533c73
"one lines! call for the graticule, one poly! call for the land"
function drawmap!(ax, land, grid, frame)
    pts = Point2f[]
    for r in grid
        append!(pts, Point2f.(r)); push!(pts, Point2f(NaN, NaN))
    end
    lines!(ax, pts, color = GRID, linewidth = 0.5)

    fpts = Point2f[]
    for r in frame
        append!(fpts, Point2f.(r)); push!(fpts, Point2f(NaN, NaN))
    end
    lines!(ax, fpts, color = GRID, linewidth = 1.4)

    poly!(ax, [Polygon(Point2f.(r)) for r in land if length(r) >= 3],
          color = LAND, strokecolor = EDGE, strokewidth = 0.25)
end

# ╔═╡ 47c8329c-2277-41e5-a377-9016ddc8c3b6
function mapaxis(fig)
    ax = Axis(fig[1, 1], backgroundcolor = OCEAN, aspect = DataAspect())
    hidedecorations!(ax)
    hidespines!(ax)
    ax
end

# ╔═╡ fc89d040-0f32-4649-b228-63808a6de745
drawtissot!(ax, circles) =
    poly!(ax, [Polygon(Point2f.(r)) for r in circles],
          color = (TISSOT, 0.55), strokecolor = TISSOT, strokewidth = 0.8)

# ╔═╡ cbf07ada-c616-4bb6-b9cc-dca51efde261
md"""
### Animating

`ease` is smoothstep, so each transition starts and ends at rest instead of snapping. The frame
list is built up front as `(from, to, t)` triples, which keeps the recording loop trivial.
"""

# ╔═╡ 6fd83ed7-1a9d-4dd6-a406-8efb2838d3d7
ease(t) = t * t * (3 - 2t)

# ╔═╡ 1ef673b2-e962-426f-ad5a-279999f36d8e
function morphmovie(path, projs, land, grid, frame;
                    nframes = 36, hold = 6, size = (1200, 700), framerate = 24)
    L = [project(P, land)  for P in projs]
    G = [project(P, grid)  for P in projs]
    F = [project(P, frame) for P in projs]

    fig = Figure(; size, backgroundcolor = OCEAN)
    ax = mapaxis(fig)
    ax.titlecolor = :white
    ax.titlesize = 24

    steps = Tuple{Int,Int,Float64}[]
    for k in eachindex(projs)
        j = mod1(k + 1, length(projs))
        append!(steps, [(k, k, 0.0) for _ in 1:hold])
        append!(steps, [(k, j, ease(t)) for t in range(0, 1, length = nframes)])
    end

    record(fig, path, steps; framerate) do (i, j, t)
        empty!(ax)
        ax.title = t < 0.5 ? projs[i].name : projs[j].name
        limits!(ax, -1.05, 1.05, -1.05, 1.05)
        drawmap!(ax, blend(L[i], L[j], t), blend(G[i], G[j], t), blend(F[i], F[j], t))
    end
    path
end

# ╔═╡ 7cc30cee-01b5-4fa2-a708-a8785491ef3d
md"""
### Distortion, and embedding the result
"""

# ╔═╡ baf317c0-7f32-4f14-b1e1-1a6e3e4fe2f1
"small circles of equal angular radius, laid out on a grid over the sphere"
function tissot(; dlat = 30, dlon = 40, radius = 5.0, n = 48)
    ρ = deg2rad(radius)
    circles = Vector{Vector{Tuple{Float64,Float64}}}()
    for lat0 in -60:dlat:60, lon0 in -160:dlon:160
        φ₀, λ₀ = deg2rad(lat0), deg2rad(lon0)
        pts = map(range(0, 2π, length = n)) do α
            φ = asin(clamp(sin(φ₀) * cos(ρ) + cos(φ₀) * sin(ρ) * cos(α), -1, 1))
            λ = λ₀ + atan(sin(α) * sin(ρ) * cos(φ₀), cos(ρ) - sin(φ₀) * sin(φ))
            (rad2deg(φ), rad2deg(λ))
        end
        push!(circles, pts)
    end
    circles
end

# ╔═╡ b97ed1f7-eb82-4ed4-8350-69a27434302c
"embed an mp4 so the exported page stays self-contained"
function video(path; width = "100%")
    data = base64encode(read(path))
    HTML("""<video width="$width" autoplay loop muted playsinline
              style="border-radius:6px; display:block">
              <source src="data:video/mp4;base64,$data" type="video/mp4">
            </video>""")
end

# ╔═╡ a06df991-9bfe-4f1b-8d98-2122d40419e1
md"""
---

*Country boundaries from [Natural Earth](https://www.naturalearthdata.com) via
[GeoArtifacts.jl](https://github.com/JuliaEarth/GeoArtifacts.jl). Geometry from
[Meshes.jl](https://github.com/JuliaGeometry/Meshes.jl). Projections from
[CoordRefSystems.jl](https://github.com/JuliaEarth/CoordRefSystems.jl). Drawing with
[Makie.jl](https://docs.makie.org).*
"""

# ╔═╡ Cell order:
# ╟─b3d68a38-bb27-40e2-a529-5b4304a320bb
# ╠═3da10fc9-0e94-4426-a8ca-fae3ff10d583
# ╟─00226b63-87b7-472d-a7a5-73f0f4df95bc
# ╠═efb0cafc-e855-456f-bc01-b895b994cf35
# ╠═23b3f9e1-ca54-49d9-af2d-5568e7b335b3
# ╠═ef070602-d668-42ba-8125-8e5d7bcdf74b
# ╟─4f6b3239-e7a9-4be7-aa1e-4032a11f87b6
# ╠═6ada4791-c3bd-4339-8047-a8879b2f11af
# ╠═895f4784-cc99-4d04-93b5-448581375280
# ╟─07e85657-41b2-47f4-a5e1-e8e0c013d917
# ╠═6db80fbd-3c1d-4f95-bed3-4cbd56b1f099
# ╟─b6cb65d5-452c-4e29-a530-24d1528be4bd
# ╠═26247572-b2d9-472c-9431-43fd3dbae91c
# ╟─870cf06b-7ff5-4270-a37d-8b72015c50b3
# ╠═f67bca02-ba46-4624-b0d6-fcb8ee3cfd98
# ╟─cbb59117-a527-409a-8002-d6cb2dfa849d
# ╠═d6844ed7-6b74-47c1-aa4c-8a96952f5021
# ╟─e0b7639b-ded3-42f4-bd8d-4650105c4760
# ╠═7798f437-6d69-4805-a5be-94bbee12c175
# ╟─ce8ecfb6-e1c4-42d5-96cd-661bc7af4841
# ╟─4646a8ed-68ec-4e3d-bd8c-1193d057cc85
# ╟─cc43fc0b-70ef-4568-b743-8ab2af2006db
# ╠═39a8f29c-e96b-49da-b0a6-c15e661aa3ad
# ╠═21b4a30e-8555-432e-a1cf-9fe0b34247f5
# ╠═0079048a-cdc3-467e-9bd9-d5525b3de382
# ╠═36db0674-73aa-4bc0-9e00-97f591c954d9
# ╟─f3caa581-116e-4275-8d6c-3561915ebcba
# ╠═1448c824-e90f-4e5d-9f55-25fc4a7e022d
# ╠═07ea1ca3-c1a7-4e7b-9c0d-a4449eaa64a6
# ╠═cf08360c-58cb-4461-ba18-94d006815508
# ╟─84688259-9341-4506-bbd2-230134f917d1
# ╠═ef944087-9b31-42e2-965d-6cae8d127eab
# ╠═8fb80132-9f6f-47e4-8c96-2342c5533c73
# ╠═47c8329c-2277-41e5-a377-9016ddc8c3b6
# ╠═fc89d040-0f32-4649-b228-63808a6de745
# ╟─cbf07ada-c616-4bb6-b9cc-dca51efde261
# ╠═6fd83ed7-1a9d-4dd6-a406-8efb2838d3d7
# ╠═1ef673b2-e962-426f-ad5a-279999f36d8e
# ╟─7cc30cee-01b5-4fa2-a708-a8785491ef3d
# ╠═baf317c0-7f32-4f14-b1e1-1a6e3e4fe2f1
# ╠═b97ed1f7-eb82-4ed4-8350-69a27434302c
# ╟─a06df991-9bfe-4f1b-8d98-2122d40419e1
