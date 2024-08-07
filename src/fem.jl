function FE_setup(;model=model, order=2)
	
	reffe = ReferenceFE(lagrangian, Float64, order)
	
	Vₑ = TestFESpace(model, reffe,
    	vector_type=Vector{ComplexF64},
    	conformity=:H1, 
    	dirichlet_tags=["Dirichlet", "Dirichlet_Nodes"])

	Vₕ = TestFESpace(model, reffe,
    	vector_type=Vector{ComplexF64},
    	conformity=:H1, 
    	dirichlet_tags=["Dirichlet", "Dirichlet_Nodes", "Air"])

	Ω = Triangulation(model)
	dΩ = Measure(Ω, 2 * order)

	reffe_0 = ReferenceFE(lagrangian, Float64, 0)
	V₀ = TestFESpace(model, reffe_0)

	return FEProblem(order, Vₑ, Vₕ, V₀, Ω, dΩ )
end

