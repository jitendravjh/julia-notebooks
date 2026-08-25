### A Pluto.jl notebook ###
# v0.20.4

#> [frontmatter]
#> title = "A conjecture that was false"
#> description = "The Jacobian conjecture stood for 87 years and broke in July 2026. Verifying the counterexample exactly, and seeing the geometry that makes it work."
#> tags = ["julia", "mathematics", "geometry"]

using Markdown
using InteractiveUtils

# ╔═╡ 0b70acc1-5803-4156-a3a8-7b9360df1a16
md"""
# A conjecture that was false

In 1939 Ott-Heinrich Keller asked a question about polynomial maps of the plane. If a
polynomial map has a derivative that is invertible at every point, must the map itself be
invertible, with a polynomial inverse?

Local invertibility everywhere, in exchange for global invertibility. It sounds too good, and
for eighty-seven years nobody could prove it or break it.

In July 2026 it broke.
"""

# ╔═╡ 14a8566d-b0bd-40d1-901c-f6864eeb1f20
begin
    using DynamicPolynomials, DataFrames
    using Meshes, CairoMakie, Unitful, Base64
    import CoordRefSystems: Cartesian
end

# ╔═╡ 8527040f-9f26-450e-b4c3-b3036d3aa7b6
md"""
## The map

Three variables, degree seven.
"""

# ╔═╡ 52e5bf0c-8520-43c5-bcd1-661f98239884
@polyvar z₁ z₂ z₃ ζ

# ╔═╡ e890d0e8-1266-449e-93a4-5844dca94b41
F = let u = 1 + z₁ * z₂
    [ u^3 * z₃ + z₂^2 * u * (4 + 3z₁ * z₂),
      z₂ + 3z₁ * u^2 * z₃ + 3z₁ * z₂^2 * (4 + 3z₁ * z₂),
      2z₁ - 3z₁^2 * z₂ - z₁^3 * z₃ ]
end

# ╔═╡ ea4ddbcc-11b3-40bd-9e47-97aef9b8882c
md"""
Its Jacobian matrix is a three by three array of polynomials.
"""

# ╔═╡ 56feac3f-6387-43ed-b7af-92a782aa09cf
DF = [differentiate(F[i], v) for i in 1:3, v in (z₁, z₂, z₃)]

# ╔═╡ 3960672f-2ca1-4711-9918-17e550eba1d8
detDF = DF[1,1] * (DF[2,2]*DF[3,3] - DF[2,3]*DF[3,2]) -
        DF[1,2] * (DF[2,1]*DF[3,3] - DF[2,3]*DF[3,1]) +
        DF[1,3] * (DF[2,1]*DF[3,2] - DF[2,2]*DF[3,1])

# ╔═╡ 991377f7-efb5-4c16-bbd1-4e53d2f5710d
md"""
Every variable cancels. A degree seven map in three variables has a Jacobian that could reach
degree eighteen, and all of it vanishes except the constant.

So `F` is locally invertible at every single point of $\mathbb{C}^3$. The conjecture says it
must then be globally invertible.
"""

# ╔═╡ 2b83e97a-d559-403f-b617-813beea16b32
md"""
## It is not

Three different points, one image.
"""

# ╔═╡ d9594023-f015-44eb-b7d9-92a45310f1cf
collisions = let
    pts = [(0//1, 0//1, -1//4), (1//1, -3//2, 13//2), (-1//1, 3//2, 13//2)]
    DataFrame(
        z₁ = [p[1] for p in pts],
        z₂ = [p[2] for p in pts],
        z₃ = [p[3] for p in pts],
        image = [Tuple(f(z₁ => p[1], z₂ => p[2], z₃ => p[3]) for f in F) for p in pts],
    )
end

# ╔═╡ 589805f6-f86b-4107-a976-8fc78819feef
md"""
Exact rational arithmetic, no rounding. The conjecture is false.
"""

# ╔═╡ b5776039-3736-4ae3-b943-039dddf9deb0
DataFrame(
    dimension = ["n = 1", "n = 2", "n ≥ 3"],
    status    = ["true", "open", "false"],
    note      = ["elementary",
                 "the case Keller actually asked about",
                 "Alpöge, July 2026"],
)

# ╔═╡ 27784306-d92b-4cdb-b62f-a8f9137daf5e
md"""
Keller posed the problem for two variables. That case is still open.

## Why such a map can exist

Take a linear form and a quadratic form in two variables and multiply them:

$$F(L, Q) = L \cdot Q$$

A cubic form has three roots. Splitting it back into linear times quadratic means choosing
which root goes to the linear factor. Three choices, three preimages.
"""

# ╔═╡ 1dc2c5fb-82dc-49e2-ae05-8a62af592b34
roots₃ = [0.8 + 0.3im, -1.1 + 0.9im, 0.2 - 1.4im]

# ╔═╡ c3dbbf49-ed39-4462-a396-4d3b1c2f5bd0
lead = 2.0 + 0.5im

# ╔═╡ cf1392ac-c1bf-4f93-a96d-b699d2e2ad9f
md"""
The pair $(L, Q)$ is only defined up to $(\lambda L, \lambda^{-1} Q)$, which leaves the product
alone. Fixing the resultant to 1 spends that freedom and pins each split to a single point.
"""

# ╔═╡ 375d9b0c-b67f-47bb-823b-0f7d01383c55
function split(f, α, k)
    j = setdiff(1:3, k)
    Δ = (α[k] - α[j[1]]) * (α[k] - α[j[2]])
    (a = 1 / (f * Δ), c = f^2 * Δ, root_of_L = α[k], roots_of_Q = α[j])
end

# ╔═╡ 600ef344-275e-4a40-90c7-9c2d4f0102bb
cubic(f, α) = f * prod(ζ - r for r in α)

# ╔═╡ 30c1003f-e84a-41b4-9740-c2a16204fde1
splits = let target = cubic(lead, roots₃)
    rows = map(1:3) do k
        s = split(lead, roots₃, k)
        L = s.a * (ζ - s.root_of_L)
        Q = s.c * prod(ζ - r for r in s.roots_of_Q)
        (which_root = k,
         L_coefficient = round(s.a, digits = 4),
         resultant = round(s.a^2 * s.c * (s.root_of_L - s.roots_of_Q[1]) *
                                         (s.root_of_L - s.roots_of_Q[2]), digits = 12),
         error_vs_target = gap(L * Q - target))
    end
    DataFrame(rows)
end

# ╔═╡ 1a44946c-9d2a-40a2-beff-81d6f12ae046
md"""
Three genuinely different pairs, all with resultant 1, all multiplying back to the same cubic
to machine precision. That is the non-injectivity.

Local injectivity comes from the same resultant. It vanishes exactly when $L$ and $Q$ share a
root, so forcing it to 1 keeps the linear factor's root away from the other two. Perturb the
cubic slightly and the three roots move slightly; you can still tell which one belongs to $L$.
Locally you can always undo the multiplication. Globally there are three ways to do it.

## The same thing, in a picture

Speyer's reading of the counterexample is that it sweeps the tangent lines of a plane curve.
That is easier to see one dimension down. Take a parabola and push a grid through

$$(t, s) \mapsto \gamma(t) + s\,\gamma'(t)$$
"""

# ╔═╡ f4b8a114-696c-4c0f-af15-a6717ccd260b
γ(t) = (t, t^2)

# ╔═╡ 559bb2ad-a645-4b60-9745-e66cd5c6d997
sweep(c) = let t = ustrip(c.x), s = ustrip(c.y)
    Cartesian(t + s, t^2 + 2t * s)
end

# ╔═╡ a499db77-72fb-4921-addc-e4f5afdb3f78
sheet = CartesianGrid((-2.0, -1.2), (2.0, 1.2), dims = (60, 40)) |> Morphological(sweep)

# ╔═╡ 86ed2006-d92a-49ef-b155-dd84f365a7ed
let
    fig = Figure(size = (900, 640), backgroundcolor = PAPER)
    ax = plainaxis(fig)
    viz!(ax, sheet, showsegments = true,
         segmentcolor = RGBAf(0.25, 0.45, 0.72, 0.45), color = RGBAf(0.95, 0.62, 0.2, 0.15))
    xs = range(-3.3, 3.3, length = 300)
    lines!(ax, xs, xs .^ 2, color = CURVE, linewidth = 3)
    limits!(ax, -3.5, 3.5, -1.6, 9.2)
    fig
end

# ╔═╡ e811f863-11f4-4893-a552-d5601d89c64f
md"""
The sheet folds along the curve itself. Inside the fold every point is covered twice, and the
grid lines crossing there belong to two different parts of the sheet.

Those two coverings are the two tangent lines you can draw to a parabola from a point below it.
"""

# ╔═╡ b19ba15e-df5f-416d-832d-0ed728af743a
tangencies(X, Y) = let r = sqrt(X^2 - Y); (X - r, X + r) end

# ╔═╡ 14794518-aca4-4d3c-a1bf-21729c22dd6f
orbit(θ) = (1.0cos(θ), 0.55sin(θ) - 0.85)

# ╔═╡ 2f4d17ed-b4f1-457d-b6cc-8becd0a35bda
video(tangentmovie("tangents.mp4"))

# ╔═╡ 1452d9d3-dfa3-4479-a4a1-8a5bc134d970
md"""
The two tangency points move continuously and never collide. Follow either one and the map
looks invertible. Step back and there are always two.

Arranging that in three dimensions, while keeping the Jacobian constant, is what took
eighty-seven years.

## What is still open

Two variables. Keller's original question.
"""

# ╔═╡ 26761d1c-8f90-4950-b98a-f12a28b2c84f
md"""
---
## Appendix
"""

# ╔═╡ 540413fa-7add-40df-af75-fb0ccb1c5320
begin
    const PAPER = RGBf(0.99, 0.98, 0.96)
    const CURVE = RGBf(0.80, 0.15, 0.25)
    const LINE  = RGBf(0.25, 0.45, 0.72)
    const MARK  = RGBf(0.95, 0.55, 0.15)
end

# ╔═╡ 95dbece4-8d82-4915-8158-fc2e776004b4
"largest coefficient left over, so an exact match reads as zero"
gap(p) = (c = coefficients(p); isempty(c) ? 0.0 : maximum(abs.(c)))

# ╔═╡ 996d175c-97ab-434f-9949-fa86ccc15eaf
function plainaxis(fig)
    ax = Axis(fig[1, 1], aspect = DataAspect(), backgroundcolor = PAPER)
    hidedecorations!(ax)
    hidespines!(ax)
    ax
end

# ╔═╡ 3136fa17-9958-4d71-814f-8cd9515fc0a9
function tangentmovie(path; nframes = 120, framerate = 24)
    fig = Figure(size = (900, 640), backgroundcolor = PAPER)
    ax = plainaxis(fig)
    limits!(ax, -2.8, 2.8, -1.8, 6.4)

    xs = range(-2.6, 2.6, length = 300)
    lines!(ax, xs, xs .^ 2, color = CURVE, linewidth = 3)

    θ = Observable(0.0)
    touch = lift(θ) do a
        X, Y = orbit(a); t₁, t₂ = tangencies(X, Y)
        [Point2f(t₁, t₁^2), Point2f(t₂, t₂^2)]
    end
    here = lift(a -> [Point2f(orbit(a)...)], θ)

    linesegments!(ax, lift((p, q) -> [q[1], p[1], q[1], p[2]], touch, here),
                  color = LINE, linewidth = 2)
    scatter!(ax, touch, color = MARK, markersize = 15, strokecolor = :black, strokewidth = 1)
    scatter!(ax, here, color = :black, markersize = 12)

    record(fig, path, range(0, 2π, length = nframes); framerate) do a
        θ[] = a
    end
    path
end

# ╔═╡ 088649c4-fcbf-4ced-a3c6-a3473a79ed5e
function video(path; width = "100%")
    data = base64encode(read(path))
    HTML("""<video width="$width" autoplay loop muted playsinline
              style="border-radius:6px; display:block">
              <source src="data:video/mp4;base64,$data" type="video/mp4">
            </video>""")
end

# ╔═╡ 3cdb960a-1362-4f90-9e6d-d3721c842a62
md"""
---
## References

1. O.-H. Keller, *Ganze Cremona-Transformationen*, Monatshefte für Mathematik und Physik **47**
   (1939), 299–306. [doi:10.1007/BF01695502](https://doi.org/10.1007/BF01695502).
   The conjecture is posed here, for two variables.

2. L. Alpöge, counterexample in dimension three, announced 19–20 July 2026. The explicit map
   used above is this one. Alpöge credited the AI model Fable with producing it; see
   references 3 and 4.

3. T. Tao, [*A digestion of the Jacobian conjecture counterexample*](https://terrytao.wordpress.com/2026/07/21/a-digestion-of-the-jacobian-conjecture-counterexample/),
   21 July 2026. Source of the map, of $\det DF = -2$, of the three colliding points, and of the
   multiplication-of-forms explanation and the resultant normalisation used here.

4. S. Gao, [*Counterexamples to the Jacobian conjecture in dimensions greater than two*](https://arxiv.org/abs/2608.00222),
   arXiv:2608.00222 [math.AG], 31 July 2026. Source for the attribution and dates
   (Alpöge 19 July, Gallagher 20 July, Speyer 23 July) and for the tangent-sweep reading.

5. nLab, [*Jacobian conjecture*](https://ncatlab.org/nlab/show/Jacobian+conjecture).
   Source for the status by dimension, and for the fact that $n = 2$ remains open.

Computation with [DynamicPolynomials.jl](https://github.com/JuliaAlgebra/DynamicPolynomials.jl)
and [Polynomials.jl](https://github.com/JuliaMath/Polynomials.jl), geometry with
[Meshes.jl](https://github.com/JuliaGeometry/Meshes.jl), drawing with
[Makie.jl](https://docs.makie.org).
"""

# ╔═╡ Cell order:
# ╟─0b70acc1-5803-4156-a3a8-7b9360df1a16
# ╠═14a8566d-b0bd-40d1-901c-f6864eeb1f20
# ╟─8527040f-9f26-450e-b4c3-b3036d3aa7b6
# ╠═52e5bf0c-8520-43c5-bcd1-661f98239884
# ╠═e890d0e8-1266-449e-93a4-5844dca94b41
# ╟─ea4ddbcc-11b3-40bd-9e47-97aef9b8882c
# ╠═56feac3f-6387-43ed-b7af-92a782aa09cf
# ╠═3960672f-2ca1-4711-9918-17e550eba1d8
# ╟─991377f7-efb5-4c16-bbd1-4e53d2f5710d
# ╟─2b83e97a-d559-403f-b617-813beea16b32
# ╠═d9594023-f015-44eb-b7d9-92a45310f1cf
# ╟─589805f6-f86b-4107-a976-8fc78819feef
# ╠═b5776039-3736-4ae3-b943-039dddf9deb0
# ╟─27784306-d92b-4cdb-b62f-a8f9137daf5e
# ╠═1dc2c5fb-82dc-49e2-ae05-8a62af592b34
# ╠═c3dbbf49-ed39-4462-a396-4d3b1c2f5bd0
# ╟─cf1392ac-c1bf-4f93-a96d-b699d2e2ad9f
# ╠═375d9b0c-b67f-47bb-823b-0f7d01383c55
# ╠═600ef344-275e-4a40-90c7-9c2d4f0102bb
# ╠═30c1003f-e84a-41b4-9740-c2a16204fde1
# ╟─1a44946c-9d2a-40a2-beff-81d6f12ae046
# ╠═f4b8a114-696c-4c0f-af15-a6717ccd260b
# ╠═559bb2ad-a645-4b60-9745-e66cd5c6d997
# ╠═a499db77-72fb-4921-addc-e4f5afdb3f78
# ╠═86ed2006-d92a-49ef-b155-dd84f365a7ed
# ╟─e811f863-11f4-4893-a552-d5601d89c64f
# ╠═b19ba15e-df5f-416d-832d-0ed728af743a
# ╠═14794518-aca4-4d3c-a1bf-21729c22dd6f
# ╠═2f4d17ed-b4f1-457d-b6cc-8becd0a35bda
# ╟─1452d9d3-dfa3-4479-a4a1-8a5bc134d970
# ╟─26761d1c-8f90-4950-b98a-f12a28b2c84f
# ╠═540413fa-7add-40df-af75-fb0ccb1c5320
# ╠═95dbece4-8d82-4915-8158-fc2e776004b4
# ╠═996d175c-97ab-434f-9949-fa86ccc15eaf
# ╠═3136fa17-9958-4d71-814f-8cd9515fc0a9
# ╠═088649c4-fcbf-4ced-a3c6-a3473a79ed5e
# ╟─3cdb960a-1362-4f90-9e6d-d3721c842a62
