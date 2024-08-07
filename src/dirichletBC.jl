function dirichlet_hom_E(f, σ, z)
    ω = 2 * pi * f
    μ₀ = pi * 4e-7
    kk = sqrt(-im * ω * μ₀ * σ)
    E = z >= 0 ? - im * kk / σ * exp.(-im * kk * z) : -im * kk ./ σ * (1.0 .- im
 * kk * z)

    return E
end

function dirichlet_hom_H(f, σ, z)
    ω = 2 * pi * f
    μ₀ = pi * 4e-7
    kk = sqrt(-im * ω * μ₀ * σ)
    H = z >= 0 ? exp.(-im * kk * z) : ComplexF64(1.0)

    return H
end

