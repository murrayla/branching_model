"""
Author: Liam Murray, murrayla@student.unimelb.edu.au
Descrption: test openCMISS-iron implementation via hexahedral element 
                cylinder inflation.
Input: runtime_files/
                    cylinder_hexa.ele
                    cylinder_hexa.nodes
        Files contain element and node data for hexahedral cylinder.
Output: vtk_files/cylinder_hexa.vtk
        vtk files of deformation under inflation
"""

import numpy as np
# import meshio
# from opencmiss.iron import iron

DIM = 3
N_N_EL = 27
QUAD_ORDER = 4
X, Y, Z = (1, 2, 3)

def coordinate_setup():
    coord_sys_n = 1
    coord_sys = iron.CoordinateSystem()
    coord_sys.CreateStart(coord_sys_n)
    coord_sys.DimensionSet(DIM)
    coord_sys.CreateFinish()
    return coord_sys

def node_setup(test_name):
    n_idx = []
    n_xyz = []
    with open("opt/test_files/runtime_files/" + test_name + ".nodes", 'r') as n_file:
        for i, line in enumerate(n_file):
            if i == 0: continue
            line = line.strip().split('\t')
            n_idx.append(int(line[0]))
            n_xyz.append(line[1:])

    return np.array(n_xyz).astype(float), n_idx, i

def element_setup(test_name):
    e_idx = []
    e_map = []
    with open("opt/test_files/runtime_files/" + test_name + ".ele", 'r') as e_file:
        for i, line in enumerate(e_file):
            if i == 0: continue
            line = line.strip().split('\t')
            e_idx.append(int(line[2]))
            e_map.append(line[3:])

    return np.array(e_map).astype(int), e_idx, i

def basis_setup():
    xi_n = DIM
    basis_n = 1
    basis = iron.Basis()
    basis.CreateStart(basis_n)
    basis.type = iron.BasisTypes.LAGRANGE_HERMITE_TP
    basis.numberOfXi = xi_n
    basis.interpolationXi = [iron.BasisInterpolationSpecifications.CUBIC_HERMITE] * numberOfXi
    basis.quadratureOrder = QUAD_ORDER
    basis.CreateFinish()
    return basis

def region_setup(coord_sys):
    region_n = 1
    region = iron.Region()
    region.CreateStart(region_n, iron.WorldRegion)
    region.CoordinateSystemSet(coord_sys)
    region.LabelSet("Region")
    region.CreateFinish()
    return region

def mesh_setup(region, n_n, e_n, basis, e_idx, e_np_map):
    nodes = iron.Nodes()
    nodes.CreateStart(region, n_n)
    nodes.CreateFinish()
    #
    mesh = iron.Mesh()
    mesh_n = 1
    mesh.CreateStart(mesh_n, region, DIM)
    mesh.NumberOfElementsSet(e_n)
    mesh.NumberOfComponentsSet(1)
    #
    # Initialise Elements
    mesh_e = iron.MeshElements()
    mesh_e.CreateStart(mesh, 1, basis)

    for i in range(e_n):
        element = e_idx[i][0]
        nodesList = list(
        map(int,[e_np_map[i][:]]))
        mesh_e.NodesSet(int(element), nodesList)

    mesh_e.CreateFinish()
    mesh.CreateFinish()

    return mesh_e, mesh

def decomposition_setup(mesh):
    decomp_n = 1
    decomp = iron.Decomposition()
    decomp.CreateStart(decomp_n, mesh)
    decomp.CreateFinish()

    return decomp

def geometric_setup(region, decomp, n_n, n_idx, n_xyz):
    geo_field_n = 1
    geo_field = iron.Field()
    geo_field.CreateStart(geo_field_n, region) #notice that the geometric field is associated with region in this function call.
    geo_field.LabelSet('Geometry')
    geo_field.MeshDecompositionSet(decomp)
    geo_field.CreateFinish()

    computationalNodeNumber = iron.ComputationalNodeNumberGet()
    for idx in range(n_n):
        n_id = n_n[idx][0]
        n_x, n_y, n_z = (n_xyz[idx][0], n_xyz[idx][1], n_xyz[idx][2])
        geo_field.ParameterSetUpdateNodeDP(
            iron.FieldVariableTypes.U, 
            iron.FieldParameterSetTypes.VALUES,
            1, 
            1, 
            n_id, 
            1, 
            n_x
        )
        geo_field.ParameterSetUpdateNodeDP(
            iron.FieldVariableTypes.U, 
            iron.FieldParameterSetTypes.VALUES,
            1, 
            1, 
            n_id, 
            2, 
            n_y
        )
        geo_field.ParameterSetUpdateNodeDP(
            iron.FieldVariableTypes.U, 
            iron.FieldParameterSetTypes.VALUES,
            1, 
            1, 
            n_id, 
            3, 
            n_z
        )

    geo_field.ParameterSetUpdateStart(
        iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES
    )
    geo_field.ParameterSetUpdateFinish(
        iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES
    )

def main(test_name):
    coord_sys = coordinate_setup()
    n_np_xyz, n_idx, n_n = node_setup(test_name)
    e_np_map, e_idx, e_n = element_setup(test_name)
    basis = basis_setup()
    region = region_setup(coord_sys)
    mesh_e, mesh = mesh_setup(region, n_n, e_n, basis, e_idx, e_np_map)
    decomp = decomposition_setup(mesh)
    geo_field = geometric_setup(region, decomp)
    print("In Progress")

if __name__ == '__main__':
    test_name = "cylinder_hexa"
    main(test_name)