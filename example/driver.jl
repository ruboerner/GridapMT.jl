using GridapMT
# domain
L = 100000.0
H = 100000.0
# prism
xc = 0.0
dx = 2000.0
zc = 2000.0
dz = 2000.0

# conductivities
const σ₀ = 1e-2 # background
const σ₁ = 1e-1 # prism
# mean freq to get reasonable skin depth estimate
freq = 1.0
skindepth = 503 / sqrt(freq * σ₀)
lc = 0.7 * skindepth # length scale

PrismGenerator(L, H, xc, dx, zc, dz, lc, "prism_4211.msh")

model = loadModel("prism_4211.msh")

Vₑ, Vₕ, Ω, dΩ, V₀ = FE_setup(model=model, order=2)
f = exp10.(range(-3, 1, length=5))
obs = [-5000.0:1000.0:5000.0;]
σ₀ = 1e-2

mtsurvey = MT(model, f, obs, σ₀)

grid = createTensorGrid(domain=(0, 1, 0, 1), partition=(12, 12), value=1.0)

(σ_field, ρ_field) = interpolate_grid_to_fe_space(grid, σ₀, V₀)

(uh_e, uh_h) = solve_eh(f[1], (Vₑ, Vₕ, dΩ), (σ_field, ρ_field), σ₀)

get_rhoa_phase((uh_e, uh_h), f[1], obs, ρ_field)