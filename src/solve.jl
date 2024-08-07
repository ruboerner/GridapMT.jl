function solve_eh(freq, spaces, maps, σ₀)
	Vₑ, Vₕ, dΩ = spaces
	σ_field, ρ_field = maps
	iω = im * 2 * pi * freq
	μ₀ = pi * 4e-7
	
	g(x) = -dirichlet_hom_E(freq, σ₀, x[2])
	h(x) = dirichlet_hom_H(freq, σ₀, x[2])
	
	Uₑ = TrialFESpace(Vₑ, [g, g])
	Uₕ = TrialFESpace(Vₕ, [h, h, h])

	# bilinear forms
	ae(u, v) = ∫( ∇(v) ⋅ ∇(u) / μ₀ + iω * (σ_field * v) * u ) * dΩ
	ah(u, v) = ∫( ρ_field * ∇(v) ⋅ ∇(u) + iω * μ₀ * v * u ) * dΩ

    # linear form
	b(v) = 0
	
	# linear FE operators
    op_e = AffineFEOperator(ae, b, Uₑ, Vₑ)
    op_h = AffineFEOperator(ah, b, Uₕ, Vₕ)

    # solve for E- and H-pol
    uh_e = solve(op_e)
    uh_h = solve(op_h)
		
	return uh_e, uh_h
end

function solve_e(freq, spaces, maps, σ₀)
    Vₑ, dΩ = spaces
	σ_field = maps
	iω = im * 2 * pi * freq
	μ₀ = pi * 4e-7
	
	g(x) = -dirichlet_hom_E(freq, σ₀, x[2])
	
	Uₑ = TrialFESpace(Vₑ, [g, g])

	# bilinear form
	ae(u, v) = ∫( ∇(v) ⋅ ∇(u) / μ₀ + iω * (σ_field * v) * u ) * dΩ

    # linear form
	b(v) = 0
	
	# linear FE operator
    op_e = AffineFEOperator(ae, b, Uₑ, Vₑ)

    # solve
    uh_e = solve(op_e)
		
	return uh_e
end

function solve_h(freq, spaces, maps, σ₀)
	Vₕ, dΩ = spaces
	ρ_field = maps
	iω = im * 2 * pi * freq
	μ₀ = pi * 4e-7
	
	h(x) = dirichlet_hom_H(freq, σ₀, x[2])
	
	Uₕ = TrialFESpace(Vₕ, [h, h, h])

	# bilinear form
	ah(u, v) = ∫( ρ_field * ∇(v) ⋅ ∇(u) + iω * μ₀ * v * u ) * dΩ

    # linear form
	b(v) = 0
	
	# linear FE operator
    op_h = AffineFEOperator(ah, b, Uₕ, Vₕ)

    # solve
    uh_h = solve(op_h)
		
	return uh_h
end
