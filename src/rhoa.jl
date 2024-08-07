function get_rhoa_phase(
    u::Tuple{SingleFieldFEFunction, SingleFieldFEFunction}, 
    f::Float64, obs::Vector{Float64}, ρ_field::SingleFieldFEFunction)
    uh_e, uh_h = u
    μ₀ = pi * 4e-7
    rhoa = zeros(length(obs), 2)
    phi = zeros(length(obs), 2)

    for (k, obs) in enumerate(obs)
            iω = im * 2 * pi * f
            ex = (uh_e)(Point(obs, 0.0))
            ey = (ρ_field)(Point(obs, 0.01)) * Gridap.Fields.gradient(uh_h)(Point(obs, 0.01))[2]
            hx = (uh_h)(Point(obs, 0.0))
            hy = -Gridap.Fields.gradient(uh_e)(Point(obs, 0.0))[2] / iω / μ₀
            hz = +Gridap.Fields.gradient(uh_e)(Point(obs, 0.0))[1] / iω / μ₀

            Z_e = ex ./ hy
            Z_h = ey ./ hx
            T = hz ./ hy

            rhoa_e = abs.(Z_e).^2 / (2 * pi * f * μ₀)
            rhoa_h = abs.(Z_h).^2 / (2 * pi * f * μ₀)
            phi_e = atan.(imag(Z_e) ./ real(Z_e)) * 180 / pi
            phi_h = atan.(imag(Z_h) ./ real(Z_h)) * 180 / pi

            rhoa[k, 1] = rhoa_e
            rhoa[k, 2] = rhoa_h
            phi[k, 1] = phi_e
            phi[k, 2] = phi_h
    end

    return rhoa, phi
end
