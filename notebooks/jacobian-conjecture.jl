### A Pluto.jl notebook ###
# v0.20.4

#> [frontmatter]
#> title = "The Jacobian conjecture in dimension three"
#> description = "Verifying Alpöge's 2026 counterexample in exact arithmetic, and accounting for how a polynomial map can be étale everywhere and still fail to be injective."
#> tags = ["julia", "mathematics", "geometry"]

using Markdown
using InteractiveUtils

# ╔═╡ 01e4e41a-2bca-4bfc-93aa-1b4656309416
md"""
# The Jacobian conjecture in dimension three

Let $F : \mathbb{C}^n \to \mathbb{C}^n$ be a polynomial map whose Jacobian determinant
$\det DF$ is a nonzero constant. By the inverse function theorem $F$ is then locally invertible
at every point. The Jacobian conjecture asserts that it is globally invertible, with a
polynomial inverse. Keller posed it for $n = 2$ in 1939 [1].

The case $n = 1$ is elementary. The conjecture was refuted for $n = 3$ in July 2026 [2], and
the case $n = 2$ that Keller actually asked about remains open [5].

What follows verifies the counterexample in exact arithmetic, then accounts for it: how a
polynomial map can be locally injective everywhere and still identify distinct points.
"""

# ╔═╡ 23e1f022-4d85-49d8-b3d8-dbdae4099f31
begin
    using DynamicPolynomials, DataFrames
    using Meshes, CairoMakie, Unitful
    import CoordRefSystems: Cartesian
end

# ╔═╡ cfbfb450-2f31-4654-8830-b6a9cdb5af3b
begin
    const PAPER = RGBf(0.99, 0.98, 0.96)
    const CURVE = RGBf(0.80, 0.15, 0.25)
    const LCOL  = RGBf(0.85, 0.35, 0.13)
    const QCOL  = RGBf(0.20, 0.35, 0.60)
    const TCOL  = RGBf(0.20, 0.50, 0.34)
    const WASH  = RGBf(0.82, 0.87, 0.93)

    function plainaxis(fig, pos = fig[1, 1]; kw...)
        ax = Axis(pos; aspect = DataAspect(), backgroundcolor = PAPER, titlesize = 15, kw...)
        hidedecorations!(ax)
        hidespines!(ax)
        ax
    end
end

# ╔═╡ 907985e9-6c0a-41f2-afbe-9e3913730444
md"""
## The map

Alpöge's counterexample has degree seven in three variables. The form below is the one given in
Tao's account [3].
"""

# ╔═╡ 09ae5760-dbc4-411b-b4a7-c896afc766b2
@polyvar z₁ z₂ z₃ ζ

# ╔═╡ 623f5e3b-acfa-40e2-824c-019490a4bf98
F = let u = 1 + z₁ * z₂
    [ u^3 * z₃ + z₂^2 * u * (4 + 3z₁ * z₂),
      z₂ + 3z₁ * u^2 * z₃ + 3z₁ * z₂^2 * (4 + 3z₁ * z₂),
      2z₁ - 3z₁^2 * z₂ - z₁^3 * z₃ ]
end

# ╔═╡ 003cd4cd-799e-466f-954a-78b2523d404b
md"""
The entries of $DF$ have degree at most six, so a priori $\det DF$ could have degree as high as
eighteen. It does not.
"""

# ╔═╡ 59887c43-e632-47fa-aff0-daaccda272ff
DF = [differentiate(F[i], v) for i in 1:3, v in (z₁, z₂, z₃)]

# ╔═╡ fba1c3f1-f6b1-43e3-9a7c-8d448cd9d0a7
detDF = DF[1,1] * (DF[2,2]*DF[3,3] - DF[2,3]*DF[3,2]) -
        DF[1,2] * (DF[2,1]*DF[3,3] - DF[2,3]*DF[3,1]) +
        DF[1,3] * (DF[2,1]*DF[3,2] - DF[2,2]*DF[3,1])

# ╔═╡ 1f524210-d3a5-474e-bebc-f155967ab847
md"""
Every non-constant coefficient cancels. The map is therefore étale: its derivative is invertible
at every point of $\mathbb{C}^3$.

## Failure of injectivity

The fibre over $(-1/4, 0, 0)$ contains at least three points.
"""

# ╔═╡ 3f10066b-5035-4424-93eb-87c4c2b79dd4
fibre = let pts = [(0//1, 0//1, -1//4), (1//1, -3//2, 13//2), (-1//1, 3//2, 13//2)]
    DataFrame(
        z₁ = first.(pts),
        z₂ = getindex.(pts, 2),
        z₃ = last.(pts),
        image = [Tuple(f(z₁ => p[1], z₂ => p[2], z₃ => p[3]) for f in F) for p in pts],
    )
end

# ╔═╡ 68bab77f-f7fa-4765-8426-37f0631ad6fa
md"""
The coordinates are rationals and the arithmetic is exact, so this settles the matter rather
than providing numerical evidence for it. Adjoining identity coordinates carries the
counterexample into every dimension above three [4].
"""

# ╔═╡ 59ace16f-27dd-4300-9bb9-366f656fb2e6
DataFrame(
    n          = ["1", "2", "≥ 3"],
    conjecture = ["true", "open", "false"],
    reference  = ["elementary", "Keller 1939 [1]", "Alpöge 2026 [2]"],
)

# ╔═╡ 096786f4-2b68-4dc3-8020-77057f60159a
md"""
## Where the extra preimages come from

Tao accounts for the construction in terms of multiplication of binary forms [3]. Writing
$\mathrm{Sym}^k$ for the homogeneous polynomials of degree $k$ in two variables, consider

$$F(L, Q) = LQ, \qquad \mathrm{Sym}^1 \times \mathrm{Sym}^2 \to \mathrm{Sym}^3 .$$

A binary cubic with distinct roots factors as a linear form times a quadratic form in exactly
three ways, one for each choice of the root assigned to the linear factor.
"""

# ╔═╡ 8efa5cfa-02af-4de0-9b17-1737c5756f13
roots₃ = [0.8 + 0.3im, -1.1 + 0.9im, 0.2 - 1.4im]

# ╔═╡ 67c5717a-346c-46c1-9dc1-9b9b1d1705f4
let
    fig = Figure(size = (1040, 380), backgroundcolor = PAPER)
    for k in 1:3
        ax = plainaxis(fig, fig[1, k], title = "L takes root $k")
        lines!(ax, [-1.8, 1.6], [0, 0], color = (:black, 0.18), linewidth = 1)
        lines!(ax, [0, 0], [-2.0, 1.6], color = (:black, 0.18), linewidth = 1)
        j = setdiff(1:3, k)
        scatter!(ax, [Point2f(reim(roots₃[i])...) for i in j], color = QCOL,
                 markersize = 15, strokecolor = :black, strokewidth = 1)
        scatter!(ax, [Point2f(reim(roots₃[k])...)], color = LCOL,
                 markersize = 23, strokecolor = :black, strokewidth = 1)
        limits!(ax, -1.8, 1.6, -2.0, 1.6)
    end
    fig
end

# ╔═╡ 62b2b319-70bd-4bbd-9ce5-b1d5b8fc32cc
md"""
The three roots in the complex plane, with the one assigned to the linear factor in orange and
the two left to the quadratic in blue. Same cubic, three preimages.

The pair $(L, Q)$ is determined only up to $(\lambda L, \lambda^{-1} Q)$, which leaves the
product unchanged. For $L = az + bw$ and $Q = cz^2 + dzw + ew^2$ the resultant

$$\mathrm{Res}(L, Q) = a^2 e - abd + b^2 c = a^2 c\,(\alpha - \beta_1)(\alpha - \beta_2)$$

scales by $\lambda$ under that action, so imposing $\mathrm{Res} = 1$ selects one representative
per orbit.
"""

# ╔═╡ 569d5117-8695-4a97-ac9e-7ebb8d458e6e
lead = 2.0 + 0.5im

# ╔═╡ a8725f4c-4b10-4e4b-b27a-d52ac6651690
"a binary cubic, dehomogenised, from its leading coefficient and roots"
cubic(f, α) = f * prod(ζ - r for r in α)

# ╔═╡ fec8ec57-332c-4977-862e-83ba3b4a82ae
"the factorisation assigning root k to the linear form, normalised to Res(L,Q) = 1"
function split(f, α, k)
    j = setdiff(1:3, k)
    Δ = (α[k] - α[j[1]]) * (α[k] - α[j[2]])
    (a = 1 / (f * Δ), c = f^2 * Δ, root_of_L = α[k], roots_of_Q = α[j])
end

# ╔═╡ 23e40b54-7fd2-48ae-bfa7-cd1b2bfca44a
"largest surviving coefficient, so an exact cancellation reads as zero"
gap(p) = (c = coefficients(p); isempty(c) ? 0.0 : maximum(abs.(c)))

# ╔═╡ 95639604-b2b9-41dc-9589-9f3fa83a3f8a
splits = let target = cubic(lead, roots₃)
    DataFrame(map(1:3) do k
        s = split(lead, roots₃, k)
        L = s.a * (ζ - s.root_of_L)
        Q = s.c * prod(ζ - r for r in s.roots_of_Q)
        (assigned_root = k,
         leading_coefficient_of_L = round(s.a, digits = 4),
         resultant = round(s.a^2 * s.c * (s.root_of_L - s.roots_of_Q[1]) *
                                         (s.root_of_L - s.roots_of_Q[2]), digits = 12),
         residual = gap(L * Q - target))
    end)
end

# ╔═╡ beb269fd-edc2-481e-ab62-ee643f8ecbb1
md"""
Three distinct normalised pairs, one product.

Local injectivity follows from the same resultant. It vanishes precisely when $L$ and $Q$ share
a root, so $\mathrm{Res} = 1$ keeps the root of $L$ away from the two roots of $Q$. Perturbing
the product moves all three roots slightly and cannot reassign them, which recovers $(L, Q)$
from $LQ$ on a neighbourhood. The map is therefore étale and generically three to one.

## Why the local inverses do not patch together

That the three branches never assemble into one global inverse is visible directly. Take the
family $C_\theta(z) = z^3 - e^{i\theta}$ and let $\theta$ run from $0$ to $2\pi$. The cubic
returns to itself, but its roots do not.
"""

# ╔═╡ 37e4d3c4-eea9-4b27-9346-586c8c687dc0
rootpaths = let θs = range(0, 2π, length = 400)
    [[exp(im * (θ + 2π * k) / 3) for θ in θs] for k in 0:2]
end

# ╔═╡ 70ef98cd-f9b6-4575-91c6-821772d88e6d
let
    fig = Figure(size = (720, 660), backgroundcolor = PAPER)
    ax = plainaxis(fig)
    lines!(ax, [-1.35, 1.35], [0, 0], color = (:black, 0.18), linewidth = 1)
    lines!(ax, [0, 0], [-1.35, 1.35], color = (:black, 0.18), linewidth = 1)
    for (path, col) in zip(rootpaths, (LCOL, QCOL, TCOL))
        lines!(ax, [Point2f(reim(z)...) for z in path], color = col, linewidth = 3.5)
        scatter!(ax, [Point2f(reim(path[1])...)], color = col, markersize = 18,
                 strokecolor = :black, strokewidth = 1)
        scatter!(ax, [Point2f(reim(path[end])...)], color = col, markersize = 13,
                 marker = :rect, strokecolor = :black, strokewidth = 1)
    end
    limits!(ax, -1.35, 1.35, -1.35, 1.35)
    fig
end

# ╔═╡ 797f7070-15e1-49d4-aec3-1cdb9e0dd6be
md"""
Circles mark $\theta = 0$, squares mark $\theta = 2\pi$. Each root has travelled a third of a
turn and landed on its neighbour's starting position, so the loop permutes the fibre cyclically.
"""

# ╔═╡ 279db3a3-9fe5-4d03-8c80-f58d10fd3ce6
monodromy = let start = first.(rootpaths), finish = last.(rootpaths)
    [argmin(abs.(start .- f)) for f in finish]
end

# ╔═╡ 55df08eb-c49e-4808-bd3d-93f0bfe8f19b
md"""
Following a branch continuously around a loop returns a different branch than it started on.
There is therefore no continuous choice, over the whole space of cubics, of which root belongs
to the linear factor. Local inverses exist everywhere and no global one does, which is precisely
the shape of the counterexample.

Turning this into an actual counterexample requires cutting down to a three-dimensional slice
isomorphic to affine space by polynomial changes of variable. That step is where the difficulty
lies, and Tao's post works it through [3].

## A two-dimensional model

Speyer reads the counterexample as a sweep of the tangent lines of a plane curve, where
projective duality forces a generic point to lie on several tangents [4]. The mechanism is
already visible for a parabola, whose tangents envelope it.
"""

# ╔═╡ 2753796b-f581-4d9d-9c38-6b2cecd93a75
let
    fig = Figure(size = (860, 560), backgroundcolor = PAPER)
    ax = plainaxis(fig)
    for t in range(-2.6, 2.6, length = 45)
        xs = [-3.4, 3.4]
        lines!(ax, xs, 2t .* xs .- t^2, color = (QCOL, 0.38), linewidth = 1)
    end
    xs = range(-3.4, 3.4, length = 300)
    lines!(ax, xs, xs .^ 2, color = CURVE, linewidth = 3)
    limits!(ax, -3.4, 3.4, -3.0, 7.4)
    fig
end

# ╔═╡ 04887f86-cd2a-46f0-9cf2-f3b320524c76
md"""
For a curve $\gamma$, the tangent sweep is

$$(t, s) \mapsto \gamma(t) + s\,\gamma'(t) ,$$

which is the map sending a point of the curve and a distance along its tangent to the resulting
point of the plane. Pushing a Cartesian grid through it makes the covering visible.
"""

# ╔═╡ e059d8f3-3556-4945-9885-8a3094c89099
sweep(c) = let t = ustrip(c.x), s = ustrip(c.y)
    Cartesian(t + s, t^2 + 2t * s)
end

# ╔═╡ 19ddbab1-4085-4841-ab2b-338d56e817ea
sheet = CartesianGrid((-2.0, -1.2), (2.0, 1.2), dims = (60, 40)) |> Morphological(sweep)

# ╔═╡ 9200db66-f46c-431c-86cd-8f0242fc5a20
let
    fig = Figure(size = (900, 640), backgroundcolor = PAPER)
    ax = plainaxis(fig)
    viz!(ax, sheet, showsegments = true,
         segmentcolor = RGBAf(0.23, 0.43, 0.65, 0.45), color = RGBAf(0.91, 0.57, 0.18, 0.14))
    xs = range(-3.3, 3.3, length = 300)
    lines!(ax, xs, xs .^ 2, color = CURVE, linewidth = 3)
    limits!(ax, -3.5, 3.5, -1.6, 9.2)
    fig
end

# ╔═╡ c8998e62-2b46-4925-8384-a93941b2edae
md"""
The sheet folds along the curve, which is exactly the envelope above.

A point $(X, Y)$ is in the image when $t^2 - 2Xt + Y = 0$ has a real root, that is when
$X^2 \ge Y$, and it then has one preimage per root. The convex side of the parabola is not
covered at all; everything below it is covered twice.
"""

# ╔═╡ 0564bfd4-0665-48f4-95d1-cc043514d200
let
    X = range(-3.5, 3.5, length = 400)
    Y = range(-1.6, 9.2, length = 400)
    fig = Figure(size = (600, 780), backgroundcolor = PAPER)
    ax = plainaxis(fig)
    heatmap!(ax, X, Y, [x^2 > y ? 1.0 : 0.0 for x in X, y in Y],
             colormap = [PAPER, WASH])
    xs = range(-3.5, 3.5, length = 300)
    lines!(ax, xs, xs .^ 2, color = CURVE, linewidth = 3)
    text!(ax, 0, 6.0, text = "no preimages", align = (:center, :center), fontsize = 17)
    text!(ax, 0, -0.9, text = "two preimages", align = (:center, :center), fontsize = 17)
    limits!(ax, -3.5, 3.5, -1.6, 9.2)
    fig
end

# ╔═╡ 1d5be1f3-bf80-4c39-a5f9-42aa38a7b088
md"""
This model map is not itself a counterexample, and the reason is worth stating.
"""

# ╔═╡ ebf6b9a3-e4ad-49b9-a8d6-62dd7b7fc860
"the Jacobian determinant of the tangent sweep, in its own parameters"
function sweepjacobian()
    @polyvar t s
    S = [t + s, t^2 + 2t * s]
    J = [differentiate(S[i], v) for i in 1:2, v in (t, s)]
    J[1,1] * J[2,2] - J[1,2] * J[2,1]
end

# ╔═╡ 5d56ea19-2474-45ba-97ce-997214d23693
sweepjacobian()

# ╔═╡ 9183d589-eca5-4b94-82fd-6caffa986595
md"""
The determinant vanishes on $s = 0$, which is the curve itself, so the sweep is ramified there
rather than étale. Alpöge's map achieves the multiplicity without the ramification; injectivity
fails only through points escaping to infinity [4].

The two preimages of a point below the curve are the two tangents through it.
"""

# ╔═╡ 821c8d94-1f25-45fd-bf8d-5b1c0cb8ecdf
orbit(θ) = (0.8cos(θ), 0.45sin(θ) - 0.7)

# ╔═╡ 8640d6e0-dd54-462d-ab3b-5e6fbfb18a04
"the two points of tangency on y = x² for the tangents through (X, Y)"
tangencies(X, Y) = let r = sqrt(X^2 - Y); (X - r, X + r) end

# ╔═╡ 6cebb7dd-9fe9-4799-87d8-372fbfa56168
"""
The tangent construction as an animated SVG: frames are precomputed and interpolated by SMIL,
so the figure is a few kilobytes of vector graphics and falls back to its first frame where SMIL
is unavailable.
"""
function tangentsvg(; nframes = 72, dur = 9, xlim = (-2.4, 2.4), ylim = (-1.5, 4.4), scale = 108)
    x₀, x₁ = xlim; y₀, y₁ = ylim
    W = round(Int, (x₁ - x₀) * scale); H = round(Int, (y₁ - y₀) * scale)
    sx(x) = round((x - x₀) * scale, digits = 1)
    sy(y) = round((y₁ - y) * scale, digits = 1)

    θs = range(0, 2π, length = nframes + 1)
    P  = [orbit(θ) for θ in θs]
    T  = [tangencies(p...) for p in P]

    anim(attr, f) = """<animate attributeName="$attr" dur="$(dur)s" repeatCount="indefinite" """ *
                    """values="$(join((f(i) for i in eachindex(θs)), ";"))"/>"""
    curve = join(("$(sx(x)),$(sy(x^2))" for x in range(x₀, x₁, length = 260)), " ")

    line(k) = """<line x1="$(sx(P[1][1]))" y1="$(sy(P[1][2]))" x2="$(sx(T[1][k]))"
                  y2="$(sy(T[1][k]^2))" stroke="#3a6ea5" stroke-width="1.6">
        $(anim("x1", i -> sx(P[i][1])))$(anim("y1", i -> sy(P[i][2])))
        $(anim("x2", i -> sx(T[i][k])))$(anim("y2", i -> sy(T[i][k]^2)))</line>"""
    dot(k) = """<circle cx="$(sx(T[1][k]))" cy="$(sy(T[1][k]^2))" r="5"
                 fill="#e8912f" stroke="#1a1a1a" stroke-width="1">
        $(anim("cx", i -> sx(T[i][k])))$(anim("cy", i -> sy(T[i][k]^2)))</circle>"""

    HTML("""<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 $W $H" width="100%"
              style="max-width:$(W)px;display:block;margin:0 auto">
      <rect width="$W" height="$H" fill="#fcfaf5"/>
      <polyline points="$curve" fill="none" stroke="#cc2637" stroke-width="2.6"/>
      $(line(1))$(line(2))$(dot(1))$(dot(2))
      <circle cx="$(sx(P[1][1]))" cy="$(sy(P[1][2]))" r="4.5" fill="#1a1a1a">
        $(anim("cx", i -> sx(P[i][1])))$(anim("cy", i -> sy(P[i][2])))</circle>
    </svg>""")
end

# ╔═╡ 988df9c6-fdc9-4737-b575-d358c6fd2be4
tangentsvg()

# ╔═╡ c6bc0664-0cd3-40f4-8bc1-ce34f5daa324
md"""
The tangency points move continuously and never coincide, so either preimage can be followed
locally. There are nonetheless always two, and as with the cubic there is no way to choose one
of them consistently everywhere.

## What remains open

The case $n = 2$.
"""

# ╔═╡ 8a687849-e7ff-414e-ab32-2389914460e8
md"""
---
## References

[1] O.-H. Keller, *Ganze Cremona-Transformationen*, Monatshefte für Mathematik und Physik **47**
(1939), 299–306. [doi:10.1007/BF01695502](https://doi.org/10.1007/BF01695502). The conjecture
is posed here, for two variables.

[2] L. Alpöge, counterexample in dimension three, announced 19–20 July 2026. The map used above
is this one. Alpöge credited the AI model Fable with producing it; see [3] and [4].

[3] T. Tao, [*A digestion of the Jacobian conjecture counterexample*](https://terrytao.wordpress.com/2026/07/21/a-digestion-of-the-jacobian-conjecture-counterexample/),
21 July 2026. Source of the explicit map, of $\det DF = -2$, of the three points in the fibre,
and of the multiplication-of-forms account together with the resultant normalisation.

[4] S. Gao, [*Counterexamples to the Jacobian conjecture in dimensions greater than two*](https://arxiv.org/abs/2608.00222),
arXiv:2608.00222 [math.AG], 31 July 2026. Source for attribution and dates (Alpöge 19 July,
Gallagher 20 July, Speyer 23 July), for the tangent-sweep reading, and for the description of
these maps as étale coverings that fail to be proper.

[5] nLab, [*Jacobian conjecture*](https://ncatlab.org/nlab/show/Jacobian+conjecture). Status by
dimension, and the observation that the open case is the one Keller stated.

Polynomial arithmetic with
[DynamicPolynomials.jl](https://github.com/JuliaAlgebra/DynamicPolynomials.jl), the swept grid
with [Meshes.jl](https://github.com/JuliaGeometry/Meshes.jl), figures with
[Makie.jl](https://docs.makie.org).
"""

# ╔═╡ Cell order:
# ╟─01e4e41a-2bca-4bfc-93aa-1b4656309416
# ╠═23e1f022-4d85-49d8-b3d8-dbdae4099f31
# ╟─cfbfb450-2f31-4654-8830-b6a9cdb5af3b
# ╟─907985e9-6c0a-41f2-afbe-9e3913730444
# ╠═09ae5760-dbc4-411b-b4a7-c896afc766b2
# ╠═623f5e3b-acfa-40e2-824c-019490a4bf98
# ╟─003cd4cd-799e-466f-954a-78b2523d404b
# ╠═59887c43-e632-47fa-aff0-daaccda272ff
# ╠═fba1c3f1-f6b1-43e3-9a7c-8d448cd9d0a7
# ╟─1f524210-d3a5-474e-bebc-f155967ab847
# ╠═3f10066b-5035-4424-93eb-87c4c2b79dd4
# ╟─68bab77f-f7fa-4765-8426-37f0631ad6fa
# ╠═59ace16f-27dd-4300-9bb9-366f656fb2e6
# ╟─096786f4-2b68-4dc3-8020-77057f60159a
# ╠═8efa5cfa-02af-4de0-9b17-1737c5756f13
# ╠═67c5717a-346c-46c1-9dc1-9b9b1d1705f4
# ╟─62b2b319-70bd-4bbd-9ce5-b1d5b8fc32cc
# ╠═569d5117-8695-4a97-ac9e-7ebb8d458e6e
# ╠═a8725f4c-4b10-4e4b-b27a-d52ac6651690
# ╠═fec8ec57-332c-4977-862e-83ba3b4a82ae
# ╟─23e40b54-7fd2-48ae-bfa7-cd1b2bfca44a
# ╠═95639604-b2b9-41dc-9589-9f3fa83a3f8a
# ╟─beb269fd-edc2-481e-ab62-ee643f8ecbb1
# ╠═37e4d3c4-eea9-4b27-9346-586c8c687dc0
# ╠═70ef98cd-f9b6-4575-91c6-821772d88e6d
# ╟─797f7070-15e1-49d4-aec3-1cdb9e0dd6be
# ╠═279db3a3-9fe5-4d03-8c80-f58d10fd3ce6
# ╟─55df08eb-c49e-4808-bd3d-93f0bfe8f19b
# ╠═2753796b-f581-4d9d-9c38-6b2cecd93a75
# ╟─04887f86-cd2a-46f0-9cf2-f3b320524c76
# ╠═e059d8f3-3556-4945-9885-8a3094c89099
# ╠═19ddbab1-4085-4841-ab2b-338d56e817ea
# ╠═9200db66-f46c-431c-86cd-8f0242fc5a20
# ╟─c8998e62-2b46-4925-8384-a93941b2edae
# ╠═0564bfd4-0665-48f4-95d1-cc043514d200
# ╟─1d5be1f3-bf80-4c39-a5f9-42aa38a7b088
# ╠═ebf6b9a3-e4ad-49b9-a8d6-62dd7b7fc860
# ╠═5d56ea19-2474-45ba-97ce-997214d23693
# ╟─9183d589-eca5-4b94-82fd-6caffa986595
# ╠═821c8d94-1f25-45fd-bf8d-5b1c0cb8ecdf
# ╠═8640d6e0-dd54-462d-ab3b-5e6fbfb18a04
# ╟─6cebb7dd-9fe9-4799-87d8-372fbfa56168
# ╠═988df9c6-fdc9-4737-b575-d358c6fd2be4
# ╟─c6bc0664-0cd3-40f4-8bc1-ce34f5daa324
# ╟─8a687849-e7ff-414e-ab32-2389914460e8
