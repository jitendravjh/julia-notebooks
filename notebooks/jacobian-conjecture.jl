### A Pluto.jl notebook ###
# v0.20.4

#> [frontmatter]
#> title = "The Jacobian conjecture in dimension three"
#> description = "Verifying Alpöge's 2026 counterexample in exact arithmetic, and accounting for how a polynomial map can be étale everywhere and still fail to be injective."
#> tags = ["julia", "mathematics", "geometry"]

using Markdown
using InteractiveUtils

# ╔═╡ 413d2b53-1bb0-4bec-9fb1-2faa6dcabc31
md"""
# The Jacobian conjecture in dimension three

Let $F : \mathbb{C}^n \to \mathbb{C}^n$ be a polynomial map whose Jacobian determinant
$\det DF$ is a nonzero constant. By the inverse function theorem $F$ is then locally invertible
at every point. The Jacobian conjecture asserts that it is globally invertible, with a
polynomial inverse. Keller posed it for $n = 2$ in 1939 [1].

The case $n = 1$ is elementary. The conjecture was refuted for $n = 3$ in July 2026 [2], and
the case $n = 2$ that Keller actually asked about remains open [5].

What follows verifies the counterexample in exact arithmetic and then accounts for it: how a
polynomial map can be locally injective everywhere and still identify distinct points.
"""

# ╔═╡ 959e82aa-04a2-495b-ac99-4c9b4f8ba050
begin
    using DynamicPolynomials, DataFrames
    using Meshes, CairoMakie, Unitful
    import CoordRefSystems: Cartesian
end

# ╔═╡ 12f3e5e2-b82b-4950-ae37-34a56a70e6a1
md"""
## The map

Alpöge's counterexample has degree seven in three variables. The form below is the one given in
Tao's account of it [3].
"""

# ╔═╡ 786e44ab-1fdb-414e-9a32-f129c4447c07
@polyvar z₁ z₂ z₃ ζ

# ╔═╡ b042d7a5-d30c-4511-881d-61dc6cfd909e
F = let u = 1 + z₁ * z₂
    [ u^3 * z₃ + z₂^2 * u * (4 + 3z₁ * z₂),
      z₂ + 3z₁ * u^2 * z₃ + 3z₁ * z₂^2 * (4 + 3z₁ * z₂),
      2z₁ - 3z₁^2 * z₂ - z₁^3 * z₃ ]
end

# ╔═╡ 32175586-b51a-47a1-b08a-7f64755dbf8d
md"""
The entries of $DF$ have degree at most six, so a priori $\det DF$ could have degree as high as
eighteen. It does not.
"""

# ╔═╡ f9db59f0-298e-49d4-9a6f-a7c8592a3bd0
DF = [differentiate(F[i], v) for i in 1:3, v in (z₁, z₂, z₃)]

# ╔═╡ 4b78a35f-f624-4d4e-8056-de737ae7e74e
detDF = DF[1,1] * (DF[2,2]*DF[3,3] - DF[2,3]*DF[3,2]) -
        DF[1,2] * (DF[2,1]*DF[3,3] - DF[2,3]*DF[3,1]) +
        DF[1,3] * (DF[2,1]*DF[3,2] - DF[2,2]*DF[3,1])

# ╔═╡ 030fe359-319c-47e6-b75d-6b02a9dbeac6
md"""
Every non-constant coefficient cancels, leaving $\det DF = -2$. The map is therefore étale: its
derivative is invertible at every point of $\mathbb{C}^3$.

## Failure of injectivity

The fibre over $(-1/4, 0, 0)$ contains at least three points.
"""

# ╔═╡ 2c48a5c2-7b09-4772-b1e5-b9538cee65e1
fibre = let pts = [(0//1, 0//1, -1//4), (1//1, -3//2, 13//2), (-1//1, 3//2, 13//2)]
    DataFrame(
        z₁ = first.(pts),
        z₂ = getindex.(pts, 2),
        z₃ = last.(pts),
        image = [Tuple(f(z₁ => p[1], z₂ => p[2], z₃ => p[3]) for f in F) for p in pts],
    )
end

# ╔═╡ c29a4626-b9b9-4a5c-b9c5-8952371e919e
md"""
The coordinates are rationals and the arithmetic is exact, so this settles the matter rather
than providing numerical evidence for it. Adjoining identity coordinates carries the
counterexample into every dimension above three [4].
"""

# ╔═╡ 6aa4e3a6-a29b-4b8c-a622-6178212a996f
DataFrame(
    n         = ["1", "2", "≥ 3"],
    conjecture = ["true", "open", "false"],
    reference  = ["elementary", "Keller 1939 [1]", "Alpöge 2026 [2]"],
)

# ╔═╡ 9c87827a-fe62-4a6d-b2a3-74cfb76b6584
md"""
## Where the extra preimages come from

Tao accounts for the construction in terms of multiplication of binary forms [3]. Writing
$\mathrm{Sym}^k$ for the homogeneous polynomials of degree $k$ in two variables, consider

$$F(L, Q) = LQ, \qquad \mathrm{Sym}^1 \times \mathrm{Sym}^2 \to \mathrm{Sym}^3 .$$

A binary cubic with distinct roots factors as a linear form times a quadratic form in exactly
three ways, one for each choice of the root assigned to the linear factor. The generic fibre of
multiplication is therefore a triple.

The pair $(L, Q)$ is determined only up to $(\lambda L, \lambda^{-1} Q)$, which leaves the
product unchanged. For $L = az + bw$ and $Q = cz^2 + dzw + ew^2$ the resultant

$$\mathrm{Res}(L, Q) = a^2 e - abd + b^2 c = a^2 c\,(\alpha - \beta_1)(\alpha - \beta_2)$$

scales by $\lambda$ under that action, so imposing $\mathrm{Res} = 1$ selects one representative
per orbit.
"""

# ╔═╡ 81aa3f66-bd2a-471b-a395-bab0e86ff608
lead = 2.0 + 0.5im

# ╔═╡ 736460ad-7b5f-4286-9ad5-7887f458d2e9
roots₃ = [0.8 + 0.3im, -1.1 + 0.9im, 0.2 - 1.4im]

# ╔═╡ 23fadc95-64db-4776-9878-51e53a22b071
"the factorisation assigning root k to the linear form, normalised to Res(L,Q) = 1"
function split(f, α, k)
    j = setdiff(1:3, k)
    Δ = (α[k] - α[j[1]]) * (α[k] - α[j[2]])
    (a = 1 / (f * Δ), c = f^2 * Δ, root_of_L = α[k], roots_of_Q = α[j])
end

# ╔═╡ f71c5308-877f-417e-b167-0594bddc77be
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

# ╔═╡ 48948931-f6a8-495a-a434-de3b00bfe6a3
md"""
Three distinct normalised pairs, one product.

Local injectivity follows from the same resultant. It vanishes precisely when $L$ and $Q$ share
a root, so $\mathrm{Res} = 1$ keeps the root of $L$ away from the two roots of $Q$. Perturbing
the product moves all three roots slightly and cannot reassign them, which recovers $(L, Q)$
from $LQ$ on a neighbourhood. The multiplication map is thus étale and generically three to one.

Turning this into a counterexample requires cutting down to a three-dimensional slice that is
isomorphic to affine space by polynomial changes of variable. That step is where the difficulty
lies, and Tao's post works it through [3].

## A two-dimensional model

Speyer reads the counterexample as a sweep of the tangent lines of a plane curve, where
projective duality forces a generic point to lie on several tangents [4]. The mechanism is
already visible for a parabola. For a curve $\gamma$, the tangent sweep is

$$(t, s) \mapsto \gamma(t) + s\,\gamma'(t) .$$
"""

# ╔═╡ 5d0b7ff5-26d9-4101-bf8f-a733e2f25498
sweep(c) = let t = ustrip(c.x), s = ustrip(c.y)
    Cartesian(t + s, t^2 + 2t * s)
end

# ╔═╡ a4bd1e04-d844-45d8-b3af-898d65bbc7ca
sheet = CartesianGrid((-2.0, -1.2), (2.0, 1.2), dims = (60, 40)) |> Morphological(sweep)

# ╔═╡ 4a19763f-b2f5-4765-9ed6-ea27726dc9f4
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

# ╔═╡ ee18a3be-0750-418d-a85b-2c154f1ab2e0
md"""
The image folds along the curve. A point $(X, Y)$ lies in it exactly when $t^2 - 2Xt + Y = 0$
has a real root, and it then has one preimage per root, so the interior of the parabola is
covered twice.

This model map is not itself a counterexample, and the reason is worth stating.
"""

# ╔═╡ 5376b54c-a364-44f1-bfaf-2c9f45ff6b92
"the Jacobian determinant of the tangent sweep, as a polynomial in its own parameters"
function sweepjacobian()
    @polyvar t s
    S = [t + s, t^2 + 2t * s]
    J = [differentiate(S[i], v) for i in 1:2, v in (t, s)]
    J[1,1] * J[2,2] - J[1,2] * J[2,1]
end

# ╔═╡ 62ba6957-5d07-4279-9bbc-b06f1a5244ba
sweepjacobian()

# ╔═╡ 970c0bad-81e5-48dc-9093-1acc60693e98
md"""
The determinant vanishes on $s = 0$, which is the curve itself, so the sweep is ramified there
rather than étale. Alpöge's map achieves the multiplicity without the ramification; injectivity
fails only through points escaping to infinity [4].

The two preimages of a point are the two tangent lines through it.
"""

# ╔═╡ 18d323ec-a6bc-4de9-a582-cbc9dd4c4532
tangentsvg()

# ╔═╡ c427d0e8-f2eb-4e36-b35f-351ec4c9ebc5
md"""
The tangency points move continuously and never coincide, so either preimage can be followed
locally. There are nonetheless always two.

## What remains open

The case $n = 2$.
"""

# ╔═╡ d5b336bf-2da7-4840-8820-a6517694fdc9
md"""
---
## Appendix
"""

# ╔═╡ ca67603b-1aca-4c4f-8dc2-62abf056973f
begin
    const PAPER = RGBf(0.99, 0.98, 0.96)
    const CURVE = RGBf(0.80, 0.15, 0.25)
end

# ╔═╡ a05406b4-5a11-4a5d-be8b-29d5088dda42
function plainaxis(fig)
    ax = Axis(fig[1, 1], aspect = DataAspect(), backgroundcolor = PAPER)
    hidedecorations!(ax)
    hidespines!(ax)
    ax
end

# ╔═╡ a2a90b50-6986-40bb-9a28-e8162f44fa2f
"a binary cubic, dehomogenised, given by its leading coefficient and roots"
cubic(f, α) = f * prod(ζ - r for r in α)

# ╔═╡ 3065aba4-b4cd-4e33-96f9-5830a6061584
"largest surviving coefficient, so an exact cancellation reads as zero"
gap(p) = (c = coefficients(p); isempty(c) ? 0.0 : maximum(abs.(c)))

# ╔═╡ 726ccc33-6fb9-4816-b7f3-20a4bff57d95
orbit(θ) = (0.8cos(θ), 0.45sin(θ) - 0.7)

# ╔═╡ dc14bef5-971e-4bd1-b59d-6d4424e2f478
"the two points of tangency on y = x² for the tangents through (X, Y)"
tangencies(X, Y) = let r = sqrt(X^2 - Y); (X - r, X + r) end

# ╔═╡ 028b87ae-ced7-48da-9585-5c2a6d15049f
"""
The tangent construction as an animated SVG. Frames are precomputed and interpolated by SMIL,
so the figure is a few kilobytes of vector graphics and degrades to its first frame where SMIL
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

# ╔═╡ 895d7e7f-f865-4b9c-85e5-c3f7b59da9e5
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
# ╟─413d2b53-1bb0-4bec-9fb1-2faa6dcabc31
# ╠═959e82aa-04a2-495b-ac99-4c9b4f8ba050
# ╟─12f3e5e2-b82b-4950-ae37-34a56a70e6a1
# ╠═786e44ab-1fdb-414e-9a32-f129c4447c07
# ╠═b042d7a5-d30c-4511-881d-61dc6cfd909e
# ╟─32175586-b51a-47a1-b08a-7f64755dbf8d
# ╠═f9db59f0-298e-49d4-9a6f-a7c8592a3bd0
# ╠═4b78a35f-f624-4d4e-8056-de737ae7e74e
# ╟─030fe359-319c-47e6-b75d-6b02a9dbeac6
# ╠═2c48a5c2-7b09-4772-b1e5-b9538cee65e1
# ╟─c29a4626-b9b9-4a5c-b9c5-8952371e919e
# ╠═6aa4e3a6-a29b-4b8c-a622-6178212a996f
# ╟─9c87827a-fe62-4a6d-b2a3-74cfb76b6584
# ╠═81aa3f66-bd2a-471b-a395-bab0e86ff608
# ╠═736460ad-7b5f-4286-9ad5-7887f458d2e9
# ╠═23fadc95-64db-4776-9878-51e53a22b071
# ╠═f71c5308-877f-417e-b167-0594bddc77be
# ╟─48948931-f6a8-495a-a434-de3b00bfe6a3
# ╠═5d0b7ff5-26d9-4101-bf8f-a733e2f25498
# ╠═a4bd1e04-d844-45d8-b3af-898d65bbc7ca
# ╠═4a19763f-b2f5-4765-9ed6-ea27726dc9f4
# ╟─ee18a3be-0750-418d-a85b-2c154f1ab2e0
# ╠═5376b54c-a364-44f1-bfaf-2c9f45ff6b92
# ╠═62ba6957-5d07-4279-9bbc-b06f1a5244ba
# ╟─970c0bad-81e5-48dc-9093-1acc60693e98
# ╠═18d323ec-a6bc-4de9-a582-cbc9dd4c4532
# ╟─c427d0e8-f2eb-4e36-b35f-351ec4c9ebc5
# ╟─d5b336bf-2da7-4840-8820-a6517694fdc9
# ╠═ca67603b-1aca-4c4f-8dc2-62abf056973f
# ╠═a05406b4-5a11-4a5d-be8b-29d5088dda42
# ╠═a2a90b50-6986-40bb-9a28-e8162f44fa2f
# ╠═3065aba4-b4cd-4e33-96f9-5830a6061584
# ╠═726ccc33-6fb9-4816-b7f3-20a4bff57d95
# ╠═dc14bef5-971e-4bd1-b59d-6d4424e2f478
# ╠═028b87ae-ced7-48da-9585-5c2a6d15049f
# ╟─895d7e7f-f865-4b9c-85e5-c3f7b59da9e5
