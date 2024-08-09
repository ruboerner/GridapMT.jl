#
function getPT(fem::FEProblem)
	nodes = fem.Ω.grid.node_coordinates
	nP = length(nodes)
	P = [[nodes[i][1]; nodes[i][2]] for i in 1:nP]
	P = hcat(P...)
	T = hcat(get_cell_node_ids(fem.Ω)...)
	return P, T
end


