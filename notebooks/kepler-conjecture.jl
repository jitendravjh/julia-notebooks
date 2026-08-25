### A Pluto.jl notebook ###
# v0.20.4

#> [frontmatter]
#> title = "The Kepler conjecture"
#> description = "Computing the quantities behind the densest sphere packing, and showing why a statement that sounds obvious resisted proof for 387 years."
#> tags = ["julia", "mathematics", "geometry"]

using Markdown
using InteractiveUtils

# ╔═╡ 71ae9d75-c1b5-472b-a5cf-5354ed91fe7d
md"""
# The Kepler conjecture

In 1611 Kepler wrote a short essay on why snowflakes have six corners. In passing he asserted
that no arrangement of equal spheres fills space more densely than the way greengrocers stack
oranges, which is $\pi/\sqrt{18} \approx 0.7405$ of the available volume [1].

He gave no argument. Neither did anyone else for the next three hundred and eighty-seven years.
Hales announced a proof in 1998, it appeared in the *Annals* in 2005 with the referees unable
to certify it completely, and the Flyspeck project finished a formal machine-checked proof in
2014 [2, 3, 4].

This notebook computes the quantities the conjecture is about, and shows why a statement that
sounds obvious was so hard to establish.
"""

# ╔═╡ 2b15d80b-6e4d-4f45-a169-8390174f6779
begin
    using Meshes, CairoMakie, DataFrames, LinearAlgebra, Unitful
    import GeometryBasics
end

# ╔═╡ 1db29877-6af6-41b4-9c1b-8dd358f5b715
begin
    const PAPER = RGBf(0.99, 0.98, 0.96)
    const CA    = RGBf(0.85, 0.42, 0.20)
    const CB    = RGBf(0.26, 0.45, 0.68)
    const CC    = RGBf(0.32, 0.58, 0.42)
    const CD    = RGBf(0.62, 0.36, 0.60)
    const INK   = RGBf(0.20, 0.22, 0.25)

    function flat(fig, pos = fig[1, 1]; kw...)
        ax = Axis(pos; aspect = DataAspect(), backgroundcolor = PAPER, titlesize = 15, kw...)
        hidedecorations!(ax); hidespines!(ax)
        ax
    end

    function solid(fig, pos = fig[1, 1]; az = 0.35π, el = 0.20π, kw...)
        ax = Axis3(pos; backgroundcolor = PAPER, aspect = :data, titlesize = 15,
                   azimuth = az, elevation = el, perspectiveness = 0.3, kw...)
        hidedecorations!(ax); hidespines!(ax)
        ax
    end

    disc(c, r; n = 64) = [Point2f(c[1] + r*cos(t), c[2] + r*sin(t)) for t in range(0, 2π, length = n)]
end

# ╔═╡ 8f00326e-4094-417c-98e3-0eee9e71bacd
md"""
## The stack

The arrangement Kepler had in mind is the one a square pyramid of cannonballs already has. Each
layer sits in the hollows of the one below.
"""

# ╔═╡ 3dc048c2-2ef8-496a-9445-63a0e7a86299
"a square-based pyramid of unit-diameter spheres, returned with each sphere's layer index"
function cannonballs(n)
    centres = Point3f[]; layer = Int[]
    for k in 0:n-1, i in 0:n-k-1, j in 0:n-k-1
        push!(centres, Point3f(i + k/2, j + k/2, k / sqrt(2)))
        push!(layer, k)
    end
    centres, layer
end

# ╔═╡ 10f6730d-769f-4e22-ae1c-6581b9d704de
stack, stacklayer = cannonballs(4)

# ╔═╡ 70d58803-2d9c-4247-b4d4-1db9eaa387ed
let
    fig = Figure(size = (880, 720), backgroundcolor = PAPER)
    ax = solid(fig; az = 0.62π, el = 0.18π)
    for (p, k) in zip(stack, stacklayer)
        mesh!(ax, GeometryBasics.Sphere(p, 0.5f0), color = (CA, CB, CC, CD)[k+1])
    end
    fig
end

# ╔═╡ a81f3344-6de1-4c5b-9012-7efab3c62905
md"""
Nothing so far is a claim about density; it is a claim that these spheres do not overlap. The
smallest distance between two centres should be exactly one.
"""

# ╔═╡ 22cb83df-e3cd-481b-8f47-2d1d31c46672
closest(P) = minimum(norm(P[i] - P[j]) for i in 1:length(P) for j in i+1:length(P))

# ╔═╡ 67d8f37f-5758-48ad-b3e6-5db1c81e2c0e
closest(stack)

# ╔═╡ a95a9ec7-431e-4a7f-b18e-74170997d7df
md"""
## Two dimensions first

The planar version of the question is easier and was settled long before: among packings of
equal circles, the hexagonal arrangement is densest, at $\pi/\sqrt{12} \approx 0.9069$. Thue
claimed it in 1892 and Fejes Tóth gave a complete proof in 1940 [5].
"""

# ╔═╡ 83885bd9-c204-4352-a642-92439cfc521d
hexsites   = [(i + j/2, j*sqrt(3)/2) for i in -6:6, j in -6:6][:]

# ╔═╡ dce5da8e-594c-4f1d-850a-57e9642bfe15
sqsites    = [(float(i), float(j)) for i in -6:6, j in -6:6][:]

# ╔═╡ 812b6aa9-420f-4c21-9ede-1da6c344b6a1
let
    fig = Figure(size = (1000, 520), backgroundcolor = PAPER)
    for (k, (sites, name, col)) in enumerate(((hexsites, "hexagonal", CB), (sqsites, "square", CA)))
        ax = flat(fig, fig[1, k], title = name)
        for c in sites
            poly!(ax, disc(c, 0.5), color = (col, 0.75), strokecolor = INK, strokewidth = 0.7)
        end
        limits!(ax, -3.2, 3.2, -3.2, 3.2)
    end
    fig
end

# ╔═╡ 9049ca55-9bc2-448e-8f96-08f7e6a91e62
md"""
The density of a periodic packing is the fraction of a single Voronoi cell that its own circle
occupies, so the cells are worth looking at directly.
"""

# ╔═╡ 5298da1f-635f-4600-b14c-31c0a0f24652
hexcells = tesselate(PointSet([Meshes.Point(p...) for p in hexsites]), VoronoiTesselation())

# ╔═╡ e0d88893-f038-4b18-ab8d-a3d83ba2f23b
let
    fig = Figure(size = (760, 700), backgroundcolor = PAPER)
    ax = flat(fig)
    viz!(ax, hexcells, showsegments = true, segmentcolor = RGBAf(0.2, 0.22, 0.25, 0.55),
         color = RGBAf(0, 0, 0, 0))
    for c in hexsites
        poly!(ax, disc(c, 0.5), color = (CB, 0.7), strokecolor = INK, strokewidth = 0.7)
    end
    limits!(ax, -3.0, 3.0, -3.0, 3.0)
    fig
end

# ╔═╡ 09505137-b286-416a-97a8-46f8f298782b
md"""
Taking only the cells well inside the cloud, so that none of them is truncated by the boundary:
"""

# ╔═╡ b2e13e13-8bed-434a-a0ef-e365024a1053
interior = let keep = findall(p -> norm(collect(p)) < 3.0, hexsites)
    DataFrame(
        sides = [nvertices(element(hexcells, k)) for k in keep],
        cellarea = [ustrip(Meshes.area(element(hexcells, k))) for k in keep],
    )
end

# ╔═╡ 13a4dea6-081f-44aa-8b57-0b3f3ec572ef
DataFrame(
    cells          = nrow(interior),
    distinct_sides = join(unique(interior.sides), ", "),
    cell_area      = round(only(unique(round.(interior.cellarea, digits = 9))), digits = 6),
    sqrt3_over_2   = round(sqrt(3)/2, digits = 6),
    density        = round(π * 0.25 / (sqrt(3)/2), digits = 6),
    pi_over_sqrt12 = round(π / sqrt(12), digits = 6),
)

# ╔═╡ e42a4a1a-b286-4409-bac0-27efd9342e04
md"""
Every interior cell is a hexagon of the same area, and the circle fills $\pi/\sqrt{12}$ of it.

## Three dimensions

For a lattice packing the same calculation is immediate: put a sphere at each lattice point with
radius half the shortest nonzero lattice vector, and divide its volume by the volume of a
fundamental cell.
"""

# ╔═╡ 6f3e6a40-d020-4921-877a-4afde4b47dfb
"length of the shortest nonzero vector of the lattice spanned by the columns of B"
function shortest(B, dim)
    best = Inf
    for c in Iterators.product(ntuple(_ -> -3:3, dim)...)
        all(iszero, c) && continue
        n = norm(B * collect(Float64, c))
        n > 1e-9 && n < best && (best = n)
    end
    best
end

# ╔═╡ abfc83a6-cbc8-419e-9618-f33c819df861
function latticedensity(B, dim)
    r = shortest(B, dim) / 2
    ball = dim == 2 ? π * r^2 : (4/3) * π * r^3
    ball / abs(det(B))
end

# ╔═╡ 447a2fa0-4862-4482-853c-d2af682abc96
lattices = DataFrame(
    lattice = ["hexagonal (2D)", "square (2D)", "face-centred cubic",
               "body-centred cubic", "simple cubic"],
    dimension = [2, 2, 3, 3, 3],
    density = round.([
        latticedensity([1.0 0.5; 0.0 sqrt(3)/2], 2),
        latticedensity([1.0 0.0; 0.0 1.0], 2),
        latticedensity([0.0 0.5 0.5; 0.5 0.0 0.5; 0.5 0.5 0.0], 3),
        latticedensity([1.0 0.0 0.5; 0.0 1.0 0.5; 0.0 0.0 0.5], 3),
        latticedensity(Matrix{Float64}(I, 3, 3), 3),
    ], digits = 6),
    closed_form = ["π/√12", "π/4", "π/√18", "π√3/8", "π/6"],
)

# ╔═╡ 7611a54b-0d6a-49f0-82a8-3a6c094b7c32
md"""
Face-centred cubic is the best of these, but it is not the only packing achieving
$\pi/\sqrt{18}$. Stacking close-packed layers in the order $ABAB$ instead of $ABCABC$ gives
hexagonal close packing, which is not a lattice at all yet has exactly the same density.
"""

# ╔═╡ 5dbefa38-337f-4116-a110-179863bc46b0
"close-packed layers; offsets chosen by the stacking word, e.g. \"ABC\" or \"ABA\""
function stacking(word; rings = 2)
    base = [(i + j/2, j*sqrt(3)/2) for i in -rings:rings, j in -rings:rings][:]
    base = filter(p -> norm(collect(p)) < rings + 0.4, base)
    off = Dict('A' => (0.0, 0.0), 'B' => (0.5, sqrt(3)/6), 'C' => (0.0, sqrt(3)/3))
    centres = Point3f[]; layer = Int[]
    for (k, ch) in enumerate(word), (x, y) in base
        push!(centres, Point3f(x + off[ch][1], y + off[ch][2], (k-1) * sqrt(2/3)))
        push!(layer, k)
    end
    centres, layer
end

# ╔═╡ c0ee3494-fe07-4ba5-9921-7ad12bc8046f
DataFrame(
    stacking = ["ABC (face-centred cubic)", "ABA (hexagonal close packed)"],
    closest_centres = [round(closest(first(stacking("ABC"))), digits = 9),
                       round(closest(first(stacking("ABA"))), digits = 9)],
)

# ╔═╡ 45d6d3ff-68db-429f-ae8c-a236b9f1d704
md"""
Both are genuine packings: no two centres are closer than one diameter. The difference is only
where the third layer goes.
"""

# ╔═╡ c81ed778-eb15-4ee1-8d7c-6e5857f7cb1f
let
    fig = Figure(size = (1060, 560), backgroundcolor = PAPER)
    for (k, (word, name)) in enumerate((("ABC", "ABC — face-centred cubic"),
                                        ("ABA", "ABA — hexagonal close packed")))
        ax = solid(fig, fig[1, k]; az = 0.20π, el = 0.10π, title = name)
        centres, layer = stacking(word)
        for (p, l) in zip(centres, layer)
            mesh!(ax, GeometryBasics.Sphere(p, 0.5f0), color = (CA, CB, CC)[l])
        end
    end
    fig
end

# ╔═╡ c8485810-04fd-4688-acd3-485bd198aff9
md"""
## The cell around one sphere

In the face-centred cubic packing the Voronoi cell of a sphere is a rhombic dodecahedron: twelve
rhombic faces, one perpendicular to each of the twelve nearest neighbours.
"""

# ╔═╡ 8b9207ee-5d78-40e7-ab4b-3adf8d39dae1
rhombic = vcat([Float64[x, y, z] for x in (-1, 1), y in (-1, 1), z in (-1, 1)][:],
               [Float64[2,0,0], Float64[-2,0,0], Float64[0,2,0],
                Float64[0,-2,0], Float64[0,0,2], Float64[0,0,-2]])

# ╔═╡ 8d81e9ad-696b-43ef-89b0-7989b08f3da0
facenormals = [Float64[a, b, 0] for a in (-1, 1), b in (-1, 1)][:]

# ╔═╡ 822ff338-453e-4939-be8f-01c5744b0fa8
allnormals = vcat(facenormals,
                  [Float64[a, 0, b] for a in (-1, 1), b in (-1, 1)][:],
                  [Float64[0, a, b] for a in (-1, 1), b in (-1, 1)][:])

# ╔═╡ c9ad0bff-38d0-4968-9eac-f1166560f1a2
rhombicedges = let d = minimum(norm(rhombic[i] - rhombic[j])
                               for i in 1:length(rhombic) for j in i+1:length(rhombic))
    [(i, j) for i in 1:length(rhombic) for j in i+1:length(rhombic)
            if abs(norm(rhombic[i] - rhombic[j]) - d) < 1e-9]
end

# ╔═╡ e7d7c2b4-480d-48dc-ba17-2ebf8b11cd90
DataFrame(
    vertices = length(rhombic),
    faces = length(allnormals),
    edges = length(rhombicedges),
    every_vertex_on_a_face = all(p -> any(n -> isapprox(dot(n, p), 2; atol = 1e-9), allnormals), rhombic),
    none_outside = all(p -> all(n -> dot(n, p) ≤ 2 + 1e-9, allnormals), rhombic),
)

# ╔═╡ cfc27217-7055-4e40-9181-466372a4134e
md"""
Its inradius is $\sqrt2$ and its volume is 16, so a sphere inscribed in it fills
"""

# ╔═╡ d29979c5-b811-4d9d-928e-ede577701ac4
let ρ = 2 / sqrt(2), V = 12 * (1/3) * (0.5 * 2 * 2sqrt(2)) * ρ
    DataFrame(inradius = ρ, volume = V,
              density = (4/3) * π * ρ^3 / V, pi_over_sqrt18 = π / sqrt(18))
end

# ╔═╡ 2d08cc2f-85d2-45b9-a494-4722eefd8404
md"""
which is Kepler's number again, arrived at from a single cell rather than from a lattice. The
rotating outline below is that cell, with its inscribed sphere; under orthographic projection
the sphere stays a fixed circle however the cell turns.
"""

# ╔═╡ b6cdbc71-ce38-4ca6-b205-09dddc56293c
"""
An orthographic wireframe spun about a tilted axis, emitted as an animated SVG. Edge endpoints
are precomputed per frame and interpolated by SMIL, and the first frame is written as static
attributes so the figure still reads where SMIL is unavailable.
"""
function spin(V, E, r; nframes = 72, dur = 16, scale = 74, pad = 20)
    tilt = 0.42
    Rx = [1 0 0; 0 cos(tilt) -sin(tilt); 0 sin(tilt) cos(tilt)]
    frames = map(range(0, 2π, length = nframes + 1)) do θ
        Rz = [cos(θ) -sin(θ) 0; sin(θ) cos(θ) 0; 0 0 1]
        [Rx * Rz * v for v in V]
    end
    W = round(Int, 2 * maximum(norm.(V)) * scale + 2pad); cx = W / 2
    px(f, i) = round(cx + frames[f][i][1] * scale, digits = 1)
    py(f, i) = round(cx - frames[f][i][2] * scale, digits = 1)
    anim(a, f) = string("<animate attributeName=\"", a, "\" dur=\"", dur,
                        "s\" repeatCount=\"indefinite\" values=\"",
                        join((f(k) for k in 1:length(frames)), ";"), "\"/>")
    segs = map(E) do (i, j)
        string("<line x1=\"", px(1,i), "\" y1=\"", py(1,i), "\" x2=\"", px(1,j),
               "\" y2=\"", py(1,j), "\" stroke=\"#2f4f6f\" stroke-width=\"1.5\" ",
               "stroke-linecap=\"round\">",
               anim("x1", k -> px(k,i)), anim("y1", k -> py(k,i)),
               anim("x2", k -> px(k,j)), anim("y2", k -> py(k,j)), "</line>")
    end
    HTML(string("<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 ", W, " ", W,
        "\" width=\"100%\" style=\"max-width:", W, "px;display:block;margin:0 auto\">",
        "<rect width=\"", W, "\" height=\"", W, "\" fill=\"#fcfaf5\"/>",
        "<circle cx=\"", cx, "\" cy=\"", cx, "\" r=\"", round(r * scale, digits = 1),
        "\" fill=\"#d9822b\" fill-opacity=\"0.5\" stroke=\"#b35f10\" stroke-width=\"1.5\"/>",
        join(segs), "</svg>"))
end

# ╔═╡ 7c98a30b-0b91-49f4-8f08-6b9570c359ac
spin(rhombic, rhombicedges, sqrt(2))

# ╔═╡ 76f050d5-446e-4fe8-ae3a-6bdbf2b3a7f3
md"""
## Twelve neighbours, and Newton's thirteenth

In the face-centred cubic packing every sphere touches twelve others, their centres at the
vertices of a cuboctahedron. In 1694 Newton and Gregory disagreed over whether thirteen was
possible; Newton said no, Gregory thought yes, and the question was not settled until Schütte
and van der Waerden in 1953 [6].

Newton was right, but Gregory's intuition was not unreasonable, because twelve spheres can be
arranged much more evenly than the cuboctahedron manages.
"""

# ╔═╡ 71d4e112-8902-4ba2-b389-35b274f4db88
cuboctahedron = vcat(
    [normalize(Float64[x, y, 0]) for x in (-1, 1), y in (-1, 1)][:],
    [normalize(Float64[x, 0, z]) for x in (-1, 1), z in (-1, 1)][:],
    [normalize(Float64[0, y, z]) for y in (-1, 1), z in (-1, 1)][:])

# ╔═╡ 05d036fd-3881-4629-aece-ecb07f734321
icosahedron = let φ = (1 + sqrt(5)) / 2, out = Vector{Vector{Float64}}()
    for (s, t) in Iterators.product((-1, 1), (-1, 1))
        push!(out, normalize(Float64[0, s, t*φ]))
        push!(out, normalize(Float64[s, t*φ, 0]))
        push!(out, normalize(Float64[t*φ, 0, s]))
    end
    out
end

# ╔═╡ 74cc1e10-b09e-4e8d-a3a8-fcf40d9b7ecb
shortedges(P) = let d = minimum(norm(P[i] - P[j]) for i in 1:length(P) for j in i+1:length(P))
    [(i, j) for i in 1:length(P) for j in i+1:length(P) if abs(norm(P[i] - P[j]) - d) < 1e-6]
end

# ╔═╡ ac390f9c-e608-421c-bab4-a0d56d9e26f0
let
    fig = Figure(size = (1040, 520), backgroundcolor = PAPER)
    for (k, (P, name, col)) in enumerate(((cuboctahedron, "cuboctahedral — as in the packing", CB),
                                          (icosahedron, "icosahedral — more evenly spread", CC)))
        ax = solid(fig, fig[1, k]; az = 0.30π, el = 0.16π, title = name)
        for (i, j) in shortedges(P)
            lines!(ax, [Point3f(P[i]...), Point3f(P[j]...)], color = (INK, 0.35), linewidth = 1.4)
        end
        scatter!(ax, [Point3f(p...) for p in P], color = col, markersize = 17,
                 strokecolor = :black, strokewidth = 1)
    end
    fig
end

# ╔═╡ 464b25db-0e26-4e1b-8a4a-12035767c751
minangle(P) = minimum(acosd(clamp(dot(P[i], P[j]), -1, 1))
                      for i in 1:length(P) for j in i+1:length(P))

# ╔═╡ c9aa5a60-c739-48f3-839b-fc6e85362187
DataFrame(
    arrangement = ["cuboctahedral", "icosahedral"],
    points = [length(cuboctahedron), length(icosahedron)],
    least_angular_separation = round.([minangle(cuboctahedron), minangle(icosahedron)], digits = 4),
    note = ["neighbours touch each other", "arccos(1/√5), neighbours do not touch"],
)

# ╔═╡ f36f3ef3-8760-4dae-8671-2116c2b28167
md"""
Three and a half degrees of slack, spread over twelve spheres, is enough to make a thirteenth
look plausible. It is not enough to fit one.

## Why it took three hundred and eighty-seven years

The natural strategy is local: show that the Voronoi cell around any sphere in any packing has
volume at least that of the rhombic dodecahedron, and the global bound follows immediately. The
strategy fails, and it fails for a specific reason.

A regular dodecahedron with inradius equal to the sphere's radius is *smaller* than the rhombic
dodecahedron.
"""

# ╔═╡ a2511064-7fd0-41df-a14c-564e23ef1f1f
dodecahedron_local_density = let
    a = 1 / (0.5 * sqrt((25 + 11*sqrt(5)) / 10))   # edge length giving inradius 1
    V = (15 + 7*sqrt(5)) / 4 * a^3
    (4/3) * π / V
end

# ╔═╡ 77d18078-5de1-4411-87cf-931b157edbb6
DataFrame(
    cell = ["regular dodecahedron", "rhombic dodecahedron"],
    local_density = round.([dodecahedron_local_density, π / sqrt(18)], digits = 6),
    comment = ["not extendable to a packing of space", "the face-centred cubic cell"],
)

# ╔═╡ 1ec73653-ba41-47d7-b66f-12cfdcd78ddd
md"""
A single sphere can therefore sit in a cell filling $0.7547$ of its volume, comfortably above
Kepler's $0.7405$. No local argument can succeed, because locally the bound is simply false. The
dodecahedral configuration cannot be repeated throughout space — regular dodecahedra do not
tile — but ruling out every partial imitation of it is what makes the problem hard.

That statement, that the regular dodecahedron is the extreme case, is itself a theorem: the
dodecahedral conjecture, proved by Hales and McLaughlin in 1998 and published in 2010 [7].

## The proof

Hales reduced the conjecture to a finite optimisation over local configurations, scored by a
function chosen so that the face-centred cubic value is optimal, and then bounded the score of
every configuration by computer: interval arithmetic for the nonlinear parts and linear
programming for the rest, running to some three gigabytes of output [2]. The *Annals* referees
reported after four years that they were ninety-nine per cent certain of correctness but could
not fully verify the computations, and the paper was published with that caveat [3].

Hales responded by formalising the whole argument. The Flyspeck project completed a proof
checked end to end by HOL Light and Isabelle in August 2014, with the formal-proof paper
appearing in 2017 [4]. The conjecture is now a theorem in the ordinary sense, and one of the
reasons anyone outside discrete geometry has heard of proof assistants.
"""

# ╔═╡ 5f1a2ecb-6960-4925-862b-807922d92b3f
md"""
---
## References

[1] J. Kepler, *Strena seu de nive sexangula*, Frankfurt, 1611. Translated as *The Six-Cornered
Snowflake*, Oxford University Press, 1966. The assertion is made in passing, without argument.

[2] T. C. Hales, *A proof of the Kepler conjecture*,
[Annals of Mathematics **162** (2005), 1065–1185](https://annals.math.princeton.edu/2005/162-3/p01).

[3] T. C. Hales et al., *A revision of the proof of the Kepler conjecture*, Discrete &
Computational Geometry **44** (2010), 1–34. [arXiv:0902.0350](https://arxiv.org/abs/0902.0350).

[4] T. C. Hales et al., *A formal proof of the Kepler conjecture*, Forum of Mathematics, Pi
**5** (2017), e2. [arXiv:1501.02155](https://arxiv.org/abs/1501.02155). The Flyspeck project,
completed 10 August 2014.

[5] L. Fejes Tóth, *Über die dichteste Kugellagerung*, Mathematische Zeitschrift **48** (1943),
676–684; the planar result dates to A. Thue, 1892 and 1910.

[6] K. Schütte and B. L. van der Waerden, *Das Problem der dreizehn Kugeln*, Mathematische
Annalen **125** (1953), 325–334. Settles the Newton–Gregory question at twelve.

[7] T. C. Hales and S. McLaughlin, *The dodecahedral conjecture*, Journal of the American
Mathematical Society **23** (2010), 299–344. [arXiv:math/9811079](https://arxiv.org/abs/math/9811079).

Voronoi tessellation from [Meshes.jl](https://github.com/JuliaGeometry/Meshes.jl), figures with
[Makie.jl](https://docs.makie.org).
"""

# ╔═╡ Cell order:
# ╟─71ae9d75-c1b5-472b-a5cf-5354ed91fe7d
# ╠═2b15d80b-6e4d-4f45-a169-8390174f6779
# ╟─1db29877-6af6-41b4-9c1b-8dd358f5b715
# ╟─8f00326e-4094-417c-98e3-0eee9e71bacd
# ╠═3dc048c2-2ef8-496a-9445-63a0e7a86299
# ╠═10f6730d-769f-4e22-ae1c-6581b9d704de
# ╠═70d58803-2d9c-4247-b4d4-1db9eaa387ed
# ╟─a81f3344-6de1-4c5b-9012-7efab3c62905
# ╠═22cb83df-e3cd-481b-8f47-2d1d31c46672
# ╠═67d8f37f-5758-48ad-b3e6-5db1c81e2c0e
# ╟─a95a9ec7-431e-4a7f-b18e-74170997d7df
# ╠═83885bd9-c204-4352-a642-92439cfc521d
# ╠═dce5da8e-594c-4f1d-850a-57e9642bfe15
# ╠═812b6aa9-420f-4c21-9ede-1da6c344b6a1
# ╟─9049ca55-9bc2-448e-8f96-08f7e6a91e62
# ╠═5298da1f-635f-4600-b14c-31c0a0f24652
# ╠═e0d88893-f038-4b18-ab8d-a3d83ba2f23b
# ╟─09505137-b286-416a-97a8-46f8f298782b
# ╠═b2e13e13-8bed-434a-a0ef-e365024a1053
# ╠═13a4dea6-081f-44aa-8b57-0b3f3ec572ef
# ╟─e42a4a1a-b286-4409-bac0-27efd9342e04
# ╠═6f3e6a40-d020-4921-877a-4afde4b47dfb
# ╠═abfc83a6-cbc8-419e-9618-f33c819df861
# ╠═447a2fa0-4862-4482-853c-d2af682abc96
# ╟─7611a54b-0d6a-49f0-82a8-3a6c094b7c32
# ╠═5dbefa38-337f-4116-a110-179863bc46b0
# ╠═c0ee3494-fe07-4ba5-9921-7ad12bc8046f
# ╟─45d6d3ff-68db-429f-ae8c-a236b9f1d704
# ╠═c81ed778-eb15-4ee1-8d7c-6e5857f7cb1f
# ╟─c8485810-04fd-4688-acd3-485bd198aff9
# ╠═8b9207ee-5d78-40e7-ab4b-3adf8d39dae1
# ╠═8d81e9ad-696b-43ef-89b0-7989b08f3da0
# ╠═822ff338-453e-4939-be8f-01c5744b0fa8
# ╠═c9ad0bff-38d0-4968-9eac-f1166560f1a2
# ╠═e7d7c2b4-480d-48dc-ba17-2ebf8b11cd90
# ╟─cfc27217-7055-4e40-9181-466372a4134e
# ╠═d29979c5-b811-4d9d-928e-ede577701ac4
# ╟─2d08cc2f-85d2-45b9-a494-4722eefd8404
# ╟─b6cdbc71-ce38-4ca6-b205-09dddc56293c
# ╠═7c98a30b-0b91-49f4-8f08-6b9570c359ac
# ╟─76f050d5-446e-4fe8-ae3a-6bdbf2b3a7f3
# ╠═71d4e112-8902-4ba2-b389-35b274f4db88
# ╠═05d036fd-3881-4629-aece-ecb07f734321
# ╟─74cc1e10-b09e-4e8d-a3a8-fcf40d9b7ecb
# ╠═ac390f9c-e608-421c-bab4-a0d56d9e26f0
# ╠═464b25db-0e26-4e1b-8a4a-12035767c751
# ╠═c9aa5a60-c739-48f3-839b-fc6e85362187
# ╟─f36f3ef3-8760-4dae-8671-2116c2b28167
# ╠═a2511064-7fd0-41df-a14c-564e23ef1f1f
# ╠═77d18078-5de1-4411-87cf-931b157edbb6
# ╟─1ec73653-ba41-47d7-b66f-12cfdcd78ddd
# ╟─5f1a2ecb-6960-4925-862b-807922d92b3f
