mutable struct TensorGrid
    grid
    xi
    yi
end

function createTensorGrid(;domain=(0, 1, 0, 1), partition=(2, 2), value=1.0)
	xi = collect(LinRange(domain[1], domain[2], 1 + partition[1]))
	yi = collect(LinRange(domain[3], domain[4], 1 + partition[2]))
	
	nx = length(xi)
	ny = length(yi)
	grid = value * ones(Float64, ny - 1, nx - 1)
	
	return TensorGrid(grid, xi, yi)
end

function createGrid(;domain=(0, 1, 0, 1), partition=(2, 2), value=1.0)
	xi = collect(LinRange(domain[1], domain[2], 1 + partition[1]))
	yi = collect(LinRange(domain[3], domain[4], 1 + partition[2]))
	
	nx = length(xi)
	ny = length(yi)
	grid = value * ones(Float64, ny - 1, nx - 1)
	
	return (grid, xi, yi)
end

function checkerboard_index(
    x::Float64, y::Float64, 
    x_grid_lines::Vector{Float64}, y_grid_lines::Vector{Float64})

	# Find the column index by finding the largest x_grid_lines
	# value less than or equal to x
	col = findlast(g -> g <= x, x_grid_lines)

	# Find the row index by finding the largest y_grid_lines 
	# value less than or equal to y
	row = findlast(g -> g <= y, y_grid_lines)

	# If col or row is nothing, it means the coordinate is outside the grid
	if col === nothing || row === nothing
		return (nothing, nothing)
	end
	
	return (row, col)
end

function interpolate_grid_to_fe_space(
    grid::Matrix{Float64}, 
    xi::Vector{Float64}, yi::Vector{Float64}, σ₀::Float64, U::FESpace)
	function σ_interp(x, grid, xi, yi)
		if (x[1] >= xi[1]) & (x[2] >= yi[1]) & (x[1] <= xi[end]) & (x[2] <= yi[end])
			row, col = checkerboard_index(x[1], x[2], xi, yi)
			sigma = grid[row, col]
		else
			sigma = x[2] < 0.0 ? 1e-14 : σ₀
		end
		return sigma
	end
	σ_field = Gridap.interpolate(x -> σ_interp(x, grid, xi, yi), U)
	ρ_field = Gridap.interpolate(x -> 1.0 / σ_interp(x, grid, xi, yi), U)
	
	return (σ_field, ρ_field)
end

function interpolate_grid_to_fe_space(grid::TensorGrid, σ₀::Float64, U::FESpace)
	function σ_interp(x, grid)
		if (x[1] >= grid.xi[1]) & (x[2] >= grid.yi[1]) & (x[1] <= grid.xi[end]) & (x[2] <= grid.yi[end])
			row, col = checkerboard_index(x[1], x[2], grid.xi, grid.yi)
			sigma = grid.grid[row, col]
		else
			sigma = x[2] < 0.0 ? 1e-14 : σ₀
		end
		return sigma
	end
	σ_field = Gridap.interpolate(x -> σ_interp(x, grid), U)
	ρ_field = Gridap.interpolate(x -> 1.0 / σ_interp(x, grid), U)
	
	return (σ_field, ρ_field)
end
