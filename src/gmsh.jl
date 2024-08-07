function PrismGenerator(L, H, x, lx, z, lz, lc, filename)

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.clear()
    gmsh.model.add("geometry")

    # Add points
	# Air and Earth halfspaces
    gmsh.model.geo.addPoint(-L / 2, -H / 2, 0, lc, 1)
    gmsh.model.geo.addPoint(L / 2, -H / 2, 0, lc, 2)
    gmsh.model.geo.addPoint(L / 2, 0, 0, lc, 3)
    gmsh.model.geo.addPoint(L / 2, H / 2, 0, lc, 4)
    gmsh.model.geo.addPoint(-L / 2, H / 2, 0, lc, 5)
    gmsh.model.geo.addPoint(-L / 2, 0, 0, lc, 6)

	# Prism
    gmsh.model.geo.addPoint(x - lx / 2, z - lz / 2, 0, lc, 7)
    gmsh.model.geo.addPoint(x - lx / 2, z + lz / 2, 0, lc, 8)
    gmsh.model.geo.addPoint(x + lx / 2, z + lz / 2, 0, lc, 9)
    gmsh.model.geo.addPoint(x + lx / 2, z - lz / 2, 0, lc, 10)

    gmsh.model.geo.synchronize()

    # Add lines
	# -/- -> +/-
    gmsh.model.geo.addLine(1, 2, 1)
	# +/- -> +/0
    gmsh.model.geo.addLine(2, 3, 2)
	# +/0 -> +/+
    gmsh.model.geo.addLine(3, 4, 3)
	# +/+ -> -/+ 
    gmsh.model.geo.addLine(4, 5, 4)
	# -/+ -> -/0
    gmsh.model.geo.addLine(5, 6, 5)
	# -/0 -> -/-
    gmsh.model.geo.addLine(6, 1, 6)
	# -/0 -> +/0 Air-Earth interface
    gmsh.model.geo.addLine(6, 3, 7)

	# Prism:
    gmsh.model.geo.addLine(7, 8, 8)
    gmsh.model.geo.addLine(8, 9, 9)
    gmsh.model.geo.addLine(9, 10, 10)
    gmsh.model.geo.addLine(10, 7, 11)

	# Add closed polygons; sign indicates orientation
	#
    gmsh.model.geo.addCurveLoop([1, 2, -7, 6], 1) # air
    gmsh.model.geo.addCurveLoop([7, 3, 4, 5], 2) # background
    gmsh.model.geo.addCurveLoop([8, 9, 10, 11], 3) # prism

    gmsh.model.geo.synchronize()

    gmsh.model.geo.addPlaneSurface([1], 1) # air
    gmsh.model.geo.addPlaneSurface([2, 3], 2) # background
    gmsh.model.geo.addPlaneSurface([3], 3) # prism

    gmsh.model.geo.synchronize()


	# (dim=1, [lines], number)
    gmsh.model.addPhysicalGroup(1, [1, 2, 3, 4, 5, 6], 1)
    gmsh.model.setPhysicalName(1, 1, "Dirichlet")

	# (dim=0, [points], number)
    gmsh.model.addPhysicalGroup(0, [1, 2, 3, 4, 5, 6], 2)
    gmsh.model.setPhysicalName(0, 2, "Dirichlet_Nodes")

	# (dim=2, [planesurface], number)
    gmsh.model.addPhysicalGroup(2, [2], 3)
    gmsh.model.setPhysicalName(2, 3, "Earth")

    gmsh.model.addPhysicalGroup(2, [1], 4)
    gmsh.model.setPhysicalName(2, 4, "Air")

    gmsh.model.addPhysicalGroup(2, [3], 5)
    gmsh.model.setPhysicalName(2, 5, "Prism")

    gmsh.model.geo.synchronize()

	# Local refinement: Place cylinder axis at (0,0)
	# Given: Radius, Vin, Vout
    gmsh.model.mesh.field.add("Cylinder", 1)	
    gmsh.model.mesh.field.setNumber(1, "Radius", 2 * lc)
    gmsh.model.mesh.field.setNumber(1, "VIn", lc / 16)
    gmsh.model.mesh.field.setNumber(1, "VOut", lc)
    gmsh.model.mesh.field.setNumber(1, "XAxis", L / 10)
    gmsh.model.mesh.field.setNumber(1, "XCenter", 0)
    gmsh.model.mesh.field.setNumber(1, "YCenter", 0)
    gmsh.model.mesh.field.setNumber(1, "ZAxis", 0)
    gmsh.model.mesh.field.setNumber(1, "YAxis",  0)
    gmsh.model.mesh.field.setAsBackgroundMesh(1)
	
    gmsh.model.mesh.generate(2)

    # ... and save it to disk
    gmsh.write(filename)
    gmsh.finalize()
end


function loadModel(filename::String)
    model = GmshDiscreteModel(filename)
    return model
end