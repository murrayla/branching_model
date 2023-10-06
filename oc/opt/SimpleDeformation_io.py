## ----
## AUTHOR: Liam Murray
## DESCRIPTION: OpenCMISS-Iron script which takes in node and element files from gmsh and outputs unaxial deformation
## INPUTS: gmsh.msh file
## OUTPUTS: .vtk paraview file and deformation date 
## ----

## Setup
import numpy as np
import meshio
from opencmiss.iron import iron
from cvtmsh2numpy import msh2numpy

## Inputs
file_name = "SOTD.msh"
msh2numpy(file_name)
# File name
file_name_open = 'gmsh2iron'
#Inputing element file
element_file = open(file_name_open + '.ele', 'r')
node_file = open(file_name_open + '.node', 'r')

## Constants setup
X, Y, Z = (1, 2, 3)             # Constants to indicate coordinate direction
number_of_xi = 3                # 
number_of_dimensions = 3        #
interpolation_type = 8          # 7,8,9 = Simplex Elements
number_of_gauss_per_xi = 3      # Integration points along basis
use_pressure_basis = False      # 
number_of_load_increments = 10  # 
quadrature_order = 5            #
component_number = 1            #
xi_interpolation = [            
    iron.BasisInterpolationSpecifications.QUADRATIC_SIMPLEX,
    iron.BasisInterpolationSpecifications.QUADRATIC_SIMPLEX,
    iron.BasisInterpolationSpecifications.QUADRATIC_SIMPLEX]
number_of_guass_xi = [6, 6, 6]  #
variable_num = 2                #

## Coordinate Field
coordinate_system_user_number = 1
coordinate_system = iron.CoordinateSystem()
coordinate_system.CreateStart(coordinate_system_user_number)
coordinate_system.DimensionSet(number_of_dimensions)
coordinate_system.CreateFinish()

## Take files, convert data to lists, and ensure no loose nodes
element_list = list()
element_nodes = list()
# Loop through each element and store node data
for i, line in enumerate(element_file):
    element_list.append(line)
    if element_list[i].split()[0] == "second-order-tetrahedron":
        element_nodes.extend(element_list[i].split()[3:13])
# Find used nodes     
elem_nodes_used = set(element_nodes)
# Loop through each node remove each node that is not present 
node_list = list()
for i, line in enumerate(node_file):
    node_list.append(line)
    if (node_list[-1][:].split()[0] in elem_nodes_used) or (i == 0):
        continue
    else: 
        node_list.pop()

## Breakdown Node Data
num_nodes = len(elem_nodes_used)
num_coords= number_of_dimensions
# Creating variables to store node number
node_ID = [[0 for m in range(1)] for n in range(num_nodes)]
# Creating array to store coordinates
node_coords = [[0 for m in range(num_coords)] for n in range(num_nodes)]
# Store IDs
nodeIDs = list()
# Loop through to store data
for i in range(num_nodes):
    node_ID[i][0], node_coords[i][0], node_coords[i][1], node_coords[i][2] = node_list[i+1].split()  #node number, x, y, z, boundary marker
    # Converting from string to appropriate type
    node_ID[i][0] = int(node_ID[i][0])
    node_coords[i][0] = float(node_coords[i][0])
    node_coords[i][1] = float(node_coords[i][1])
    node_coords[i][2] = float(node_coords[i][2])
    nodeIDs.append(node_ID[i][0])

# Closing the file
node_file.close()

## Breakdown Element Data
# Set Element Types
types = {1: "line", 2: "triangle", 3: "quadrangle", 4: "tetrahedron", 
         5: "hexahedron", 6: "prism", 7: "pyramid", 8: "second-order-line", 
         9: "second-order-triangle", 10: "second-order-quadrangle", 
         11: "second-order-tetrahedron", 15: "point"}
# types = {11: "second-order-tetrahedron"}
element_type = list()
# Loop through file and store data
for i, line in enumerate(element_list):
    if line.split()[0] in types.values():
        element_type.append(line.split()[0])

# Store data
element_type_set = set(element_type)
element_type_count = list()
element_type_count = map(lambda x: element_type.count(x), element_type_set)
elements_breakdown = dict(zip(element_type_set, list(element_type_count)))

# element counts
num_elements = elements_breakdown["second-order-tetrahedron"]
# node counts
num_nodes_per_ele = 10
# Creating variable to store the element map
ele_node_map = [[0 for x in range(num_nodes_per_ele)] for y in range(num_elements)]
# Element indexes
ele_IDs = [[0 for m in range(1)] for n in range(num_elements)]
# Reading element data from elemfile
ele_IDs_count = 0

# Loop through lines and store maping of nodes to elements
# Here the ordering of the nodes being input is intentional
# We are converting fromg GMSH ordering which can be found here: https://gmsh.info/doc/texinfo/gmsh.html#Node-ordering
#   to OpenCMISS-iron numbering which is based off Zienkiwicz Winding
for line in element_list:
    if line.split()[0] == "second-order-tetrahedron":
        _, _, _, \
            ele_node_map[ele_IDs_count][0], ele_node_map[ele_IDs_count][1], ele_node_map[ele_IDs_count][2], \
                ele_node_map[ele_IDs_count][3],  ele_node_map[ele_IDs_count][4], ele_node_map[ele_IDs_count][7], \
                    ele_node_map[ele_IDs_count][5],  ele_node_map[ele_IDs_count][6], ele_node_map[ele_IDs_count][8], \
                        ele_node_map[ele_IDs_count][9] = line.split()
        ele_IDs_count += 1
        ele_IDs[ele_IDs_count-1] = ele_IDs_count

# for line in element_list:
#     if line.split()[0] == "second-order-tetrahedron":
#         _, _, _, \
#             ele_node_map[ele_IDs_count][2], ele_node_map[ele_IDs_count][0], ele_node_map[ele_IDs_count][3], \
#                 ele_node_map[ele_IDs_count][1],  ele_node_map[ele_IDs_count][5], ele_node_map[ele_IDs_count][6], \
#                     ele_node_map[ele_IDs_count][9],  ele_node_map[ele_IDs_count][7], ele_node_map[ele_IDs_count][8], \
#                         ele_node_map[ele_IDs_count][4] = line.split()
#         ele_IDs_count += 1
#         ele_IDs[ele_IDs_count-1] = ele_IDs_count

# # UNCOMMENT FOR UNMODIFIED NODE NUMBERING
# for line in element_list:
#     if line.split()[0] == "second-order-tetrahedron":
#         _, _, _, \
#             ele_node_map[ele_IDs_count][0], ele_node_map[ele_IDs_count][1], ele_node_map[ele_IDs_count][2], \
#                 ele_node_map[ele_IDs_count][3],  ele_node_map[ele_IDs_count][4], ele_node_map[ele_IDs_count][5], \
#                     ele_node_map[ele_IDs_count][6],  ele_node_map[ele_IDs_count][7], ele_node_map[ele_IDs_count][8], \
#                         ele_node_map[ele_IDs_count][9] = line.split()
#         ele_IDs_count += 1
#         ele_IDs[ele_IDs_count-1] = ele_IDs_count

# Closing the file
element_file.close()

## Set up basis
basis_user_number = 1
pressure_basis_user_number = 2
# Define basis
basis = iron.Basis()
basis.CreateStart(basis_user_number)
basis.TypeSet(iron.BasisTypes.SIMPLEX)
basis.NumberOfXiSet(number_of_xi)
basis.InterpolationXiSet(
    [iron.BasisInterpolationSpecifications.QUADRATIC_SIMPLEX]*number_of_xi)
basis.QuadratureOrderSet(quadrature_order) # CHECK
basis.CreateFinish()

## Define region
region_user_number = 1
region = iron.Region()
region.CreateStart(region_user_number, iron.WorldRegion)
region.CoordinateSystemSet(coordinate_system)
region.LabelSet("Region")
region.CreateFinish()

## Define nodes
nodes = iron.Nodes()
nodes.CreateStart(region, num_nodes)
nodes.CreateFinish()

print(num_elements, num_coords, num_nodes)
## Initialise Mesh
mesh = iron.Mesh()
mesh_user_number = 1
mesh.CreateStart(mesh_user_number, region, num_coords)
mesh.NumberOfElementsSet(num_elements)
mesh.NumberOfComponentsSet(component_number) # CHECK
## Initialise Elements
meshElements = iron.MeshElements()
meshElements.CreateStart(mesh, 1, basis)

# Loop through each element and store the associted nodes
for i in range(num_elements):
    element = ele_IDs[i]
    mapped_nodes = list(
    map(int,[ele_node_map[i][0], ele_node_map[i][1], ele_node_map[i][2], 
                ele_node_map[i][3], ele_node_map[i][4], ele_node_map[i][5],
                ele_node_map[i][6], ele_node_map[i][7], ele_node_map[i][8],
                ele_node_map[i][9]]))
    meshElements.NodesSet(int(element), mapped_nodes)
 
meshElements.CreateFinish()
mesh.CreateFinish()

## Define decomposition 
number_of_computational_nodes = iron.ComputationalNumberOfNodesGet()
# Perform mesh decomposition.
decomposition_user_number = 1
decomposition = iron.Decomposition()
decomposition.CreateStart(decomposition_user_number, mesh)
decomposition.TypeSet(iron.DecompositionTypes.CALCULATED)
decomposition.NumberOfDomainsSet(number_of_computational_nodes)
decomposition.CreateFinish()

## Define Geometry
geometric_field_user_number = 1
geometric_field = iron.Field()
geometric_field.CreateStart(geometric_field_user_number, region) #notice that the geometric field is associated with region in this function call.
geometric_field.MeshDecompositionSet(decomposition)
geometric_field.TypeSet(iron.FieldTypes.GEOMETRIC)
geometric_field.VariableLabelSet(iron.FieldVariableTypes.U,"Geometry")
geometric_field.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,1) # CHECK
geometric_field.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,1) # CHECK
geometric_field.ComponentMeshComponentSet(iron.FieldVariableTypes.U,3,1) # CHECK
geometric_field.CreateFinish()
# Contribute mesh to geometric field
computationalNodeNumber = iron.ComputationalNodeNumberGet()
for idx in range(num_nodes):
    curr_node_id = node_ID[idx][0]
    node_x = node_coords[idx][0]
    node_y = node_coords[idx][1]
    node_z = node_coords[idx][2]
    geometric_field.ParameterSetUpdateNodeDP(
        iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES,
        1, 1, curr_node_id, 1, node_x)
    geometric_field.ParameterSetUpdateNodeDP(
        iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES,
        1, 1, curr_node_id, 2, node_y)
    geometric_field.ParameterSetUpdateNodeDP(
        iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES,
        1, 1, curr_node_id, 3, node_z)
    
# Update the geometric field.
geometric_field.ParameterSetUpdateStart(
        iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
geometric_field.ParameterSetUpdateFinish(
        iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)

## Define dependent field 
dependent_field_user_number = 2
dependent_field = iron.Field()
dependent_field.CreateStart(dependent_field_user_number, region)
dependent_field.MeshDecompositionSet(decomposition)
dependent_field.TypeSet(iron.FieldTypes.GEOMETRIC_GENERAL) # CHECK
dependent_field.GeometricFieldSet(geometric_field)
dependent_field.DependentTypeSet(iron.FieldDependentTypes.DEPENDENT) 
dependent_field.VariableLabelSet(iron.FieldVariableTypes.U, "Dependent")
dependent_field.NumberOfVariablesSet(variable_num) # CHECK
dependent_field.NumberOfComponentsSet(iron.FieldVariableTypes.U, 4) # CHECK
dependent_field.NumberOfComponentsSet(iron.FieldVariableTypes.DELUDELN, 4) # CHECK
# tend to pressure # CHECK all of the below
if use_pressure_basis:
    dependent_field.ComponentInterpolationSet(
        iron.FieldVariableTypes.U, 4, iron.FieldInterpolationTypes.GAUSS_POINT_BASED) # CHECK
    dependent_field.ComponentMeshComponentSet(
        iron.FieldVariableTypes.U, 4, 1)
    # DELUDELN variable (forces)
    dependent_field.ComponentInterpolationSet(
        iron.FieldVariableTypes.DELUDELN, 4,
        iron.FieldInterpolationTypes.NODE_BASED)
    dependent_field.ComponentMeshComponentSet(
        iron.FieldVariableTypes.DELUDELN, 4, 1)
else:
    # Set the hydrostatic pressure to be constant within each element.
    dependent_field.ComponentInterpolationSet(
        iron.FieldVariableTypes.U, 4,
        iron.FieldInterpolationTypes.ELEMENT_BASED)
    dependent_field.ComponentInterpolationSet(
        iron.FieldVariableTypes.DELUDELN, 4,
        iron.FieldInterpolationTypes.ELEMENT_BASED)
dependent_field.CreateFinish()

## Initialise dependent field from undeformed geometry
iron.Field.ParametersToFieldParametersComponentCopy(
    geometric_field, iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1,
    dependent_field, iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1) # CHECK
iron.Field.ParametersToFieldParametersComponentCopy(
    geometric_field, iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 2,
    dependent_field, iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 2) # CHECK
iron.Field.ParametersToFieldParametersComponentCopy(
    geometric_field, iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 3, 
    dependent_field, iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 3) # CHECK
iron.Field.ComponentValuesInitialiseDP(
    dependent_field, iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 4, 0.0) # CHECK

## Define material field # CHECK
material_field_user_number = 3
material_field = iron.Field()
material_field.CreateStart(material_field_user_number, region)
material_field.TypeSet(iron.FieldTypes.MATERIAL)
material_field.MeshDecompositionSet(decomposition)
material_field.GeometricFieldSet(geometric_field)
material_field.VariableLabelSet(iron.FieldVariableTypes.U, "Material")
# Set the number of components for the Mooney Rivlin constitutive equation (2).
material_field.NumberOfComponentsSet(iron.FieldVariableTypes.U, 2) # CHECK
# Loop through components
for component in [1, 2]:
    material_field.ComponentInterpolationSet(
      iron.FieldVariableTypes.U, component,
      iron.FieldInterpolationTypes.ELEMENT_BASED)
    if interpolation_type == 4:
        # Set arc length scaling for cubic-Hermite elements.
        material_field.FieldScalingTypeSet(iron.FieldScalingTypes.ARITHMETIC_MEAN)
material_field.CreateFinish()
# Set Mooney-Rivlin constants c10 and c01 respectively.
material_field.ComponentValuesInitialiseDP(
    iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, 2.0)
material_field.ComponentValuesInitialiseDP(
    iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 2, 1.0)

## Define equation set
equations_set_field_user_number = 4
equations_set_field = iron.Field()
equations_set = iron.EquationsSet()
equations_set_user_number = 1
# Finite elasticity equation specification.
equations_set_specification = [iron.ProblemClasses.ELASTICITY,
    iron.ProblemTypes.FINITE_ELASTICITY,
    iron.EquationsSetSubtypes.MOONEY_RIVLIN]
equations_set.CreateStart(
    equations_set_user_number, region, geometric_field, equations_set_specification,
    equations_set_field_user_number, equations_set_field)
# Add the dependent field that we created earlier.
equations_set.DependentCreateStart(dependent_field_user_number, dependent_field)
equations_set.DependentCreateFinish()
# Add the material field that we created earlier.
equations_set.MaterialsCreateStart(material_field_user_number, material_field)
equations_set.MaterialsCreateFinish()
equations_set.CreateFinish()

## Define equations.
equations = iron.Equations()
equations_set.EquationsCreateStart(equations)
equations.SparsityTypeSet(iron.EquationsSparsityTypes.SPARSE)
equations.OutputTypeSet(iron.EquationsOutputTypes.NONE)
equations_set.EquationsCreateFinish()

## Define problem
problem_user_number = 1
problem = iron.Problem()
problem_specification = (
  [iron.ProblemClasses.ELASTICITY,
   iron.ProblemTypes.FINITE_ELASTICITY,
   iron.ProblemSubtypes.NONE])
problem.CreateStart(problem_user_number, problem_specification)
problem.CreateFinish()

## Create the problem control loop.
problem.ControlLoopCreateStart()
control_loop = iron.ControlLoop()
problem.ControlLoopGet([iron.ControlLoopIdentifiers.NODE], control_loop)
control_loop.MaximumIterationsSet(5)
problem.ControlLoopCreateFinish()

nonlinear_solver = iron.Solver()
linear_solver = iron.Solver()
problem.SolversCreateStart()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE], 1, nonlinear_solver)
nonlinear_solver.OutputTypeSet(iron.SolverOutputTypes.PROGRESS)
nonlinear_solver.NewtonJacobianCalculationTypeSet(
  iron.JacobianCalculationTypes.EQUATIONS)
nonlinear_solver.NewtonLinearSolverGet(linear_solver)
linear_solver.LinearTypeSet(iron.LinearSolverTypes.DIRECT)
problem.SolversCreateFinish()

solver = iron.Solver()
solver_equations = iron.SolverEquations()
problem.SolverEquationsCreateStart()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE], 1, solver)
solver.SolverEquationsGet(solver_equations)
solver_equations.SparsityTypeSet(iron.SolverEquationsSparsityTypes.SPARSE)
_ = solver_equations.EquationsSetAdd(equations_set)
problem.SolverEquationsCreateFinish()

## Define boundary conditions
boundary_conditions = iron.BoundaryConditions()
solver_equations.BoundaryConditionsCreateStart(boundary_conditions)

## MINI MYO
# Find all nodes with z = 0 or z = 1
for i, coords in enumerate(node_coords):
    boundary_conditions.AddNode(dependent_field, iron.FieldVariableTypes.U,\
                                    1,1,i+1,Y,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundary_conditions.AddNode(dependent_field, iron.FieldVariableTypes.U,\
                                    1,1,i+1,X,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundary_conditions.AddNode(dependent_field, iron.FieldVariableTypes.U,\
                                    1,1,i+1,Z,iron.BoundaryConditionsTypes.FIXED,0.0)
    # if coords[1] == 0:
    #     boundary_conditions.AddNode(dependent_field, iron.FieldVariableTypes.U,\
    #                                 1,1,i+1,Y,iron.BoundaryConditionsTypes.FIXED,0.0)
    #     if (coords[0] == 0) or (coords[0] == 0.5):
    #         boundary_conditions.AddNode(dependent_field, iron.FieldVariableTypes.U,\
    #                                   1,1,i+1,X,iron.BoundaryConditionsTypes.FIXED,0.0)
    #     if (coords[2] == 0) or (coords[2] == 0.3):
    #         boundary_conditions.AddNode(dependent_field, iron.FieldVariableTypes.U,\
    #                                  1,1,i+1,Z,iron.BoundaryConditionsTypes.FIXED,0.0)
    # elif coords[1] == 2:
    #     boundary_conditions.AddNode(dependent_field, iron.FieldVariableTypes.U,\
    #                                 1,1,i+1,Y,iron.BoundaryConditionsTypes.FIXED,0.0)
    #     if (coords[0] == -0.1) or (coords[0] == 0.15):
    #         boundary_conditions.AddNode(dependent_field, iron.FieldVariableTypes.U,\
    #                                 1,1,i+1,X,iron.BoundaryConditionsTypes.FIXED,0.0)
    #     if (coords[2] == 0) or (coords[2] == 0.3):
    #         boundary_conditions.AddNode(dependent_field, iron.FieldVariableTypes.U,\
    #                                 1,1,i+1,Z,iron.BoundaryConditionsTypes.FIXED,0.0)
    # elif coords[1] == 2.5:
    #     boundary_conditions.AddNode(dependent_field, iron.FieldVariableTypes.U,\
    #                                 1,1,i+1,Y,iron.BoundaryConditionsTypes.FIXED,0.0)
    #     if (coords[0] == 0.3) or (coords[0] == 0.75):
    #         boundary_conditions.AddNode(dependent_field, iron.FieldVariableTypes.U,\
    #                                 1,1,i+1,X,iron.BoundaryConditionsTypes.FIXED,0.0)
    #     if (coords[2] == 0) or (coords[2] == 0.3):
    #         boundary_conditions.AddNode(dependent_field, iron.FieldVariableTypes.U,\
    #                                 1,1,i+1,Z,iron.BoundaryConditionsTypes.FIXED,0.0)
    
        
solver_equations.BoundaryConditionsCreateFinish()
try:
    problem.Solve()
except:
    print("Failed")

deformed_field_user_number = 5
deformed_field = iron.Field()
deformed_field.CreateStart(deformed_field_user_number, region)
deformed_field.MeshDecompositionSet(decomposition)
deformed_field.TypeSet(iron.FieldTypes.GEOMETRIC)
deformed_field.VariableLabelSet(iron.FieldVariableTypes.U, "DeformedGeometry")
for component in [1, 2, 3]:
    deformed_field.ComponentMeshComponentSet(
            iron.FieldVariableTypes.U, component, 1)
deformed_field.CreateFinish()

pressure_field_user_number = 6
pressure_field = iron.Field()
pressure_field.CreateStart(pressure_field_user_number, region)
pressure_field.MeshDecompositionSet(decomposition)
pressure_field.VariableLabelSet(iron.FieldVariableTypes.U, "Pressure")
pressure_field.ComponentMeshComponentSet(iron.FieldVariableTypes.U, 1, 1)
pressure_field.ComponentInterpolationSet(
  iron.FieldVariableTypes.U, 1, iron.FieldInterpolationTypes.ELEMENT_BASED)
pressure_field.NumberOfComponentsSet(iron.FieldVariableTypes.U, 1)
pressure_field.CreateFinish()

# Copy deformed geometry into deformed field.
for component in [1, 2, 3]:
    dependent_field.ParametersToFieldParametersComponentCopy(
        iron.FieldVariableTypes.U,
        iron.FieldParameterSetTypes.VALUES, component,
        deformed_field, iron.FieldVariableTypes.U,
        iron.FieldParameterSetTypes.VALUES, component)

# Copy the hydrostatic pressure solutions from the dependent field into the
# pressure field.
dependent_field.ParametersToFieldParametersComponentCopy(
  iron.FieldVariableTypes.U,
  iron.FieldParameterSetTypes.VALUES, 4,
  pressure_field, iron.FieldVariableTypes.U,
  iron.FieldParameterSetTypes.VALUES, 1)


meshNodes = iron.MeshNodes()
mesh.NodesGet(1,meshNodes)
num_nodes= meshNodes.NumberOfNodesGet()
nodes_list = [[0 for _ in range(0,num_coords+2)] for _ in range(0,num_nodes)]

node_ID_hold        = [0 for _ in range(num_nodes)]
nodes_before_def    = [[0 for _ in range(num_coords)] for node in range(num_nodes)]
nodes_after_def     = [[0 for _ in range(num_coords)] for node in range(num_nodes)]

#create list object containing nodes and another one for elements
for i in range(0,num_nodes):
    node_ID_hold[i] = i+1
    nodes_before_def[i][0]    = geometric_field.ParameterSetGetNodeDP(
      iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,node_ID_hold[i],1)
    nodes_before_def[i][1]    = geometric_field.ParameterSetGetNodeDP(
      iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,node_ID_hold[i],2)
    nodes_before_def[i][2]    = geometric_field.ParameterSetGetNodeDP(
      iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,node_ID_hold[i],3)
    nodes_after_def[i][0]     = dependent_field.ParameterSetGetNodeDP(
      iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,node_ID_hold[i],1)
    nodes_after_def[i][1]     = dependent_field.ParameterSetGetNodeDP(
      iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,node_ID_hold[i],2)
    nodes_after_def[i][2]     = dependent_field.ParameterSetGetNodeDP(
      iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,1,node_ID_hold[i],3)
    nodes_list[i] = [node_ID_hold[i], nodes_before_def[i][0], nodes_before_def[i][1],
                     nodes_before_def[i][2], nodes_after_def[i][0], 
                     nodes_after_def[i][1], nodes_after_def[i][2]]

node_file = open('output_mesh.node', 'w')
node_file.writelines([str(line) + "\n" for line in nodes_list])
# Closing the file
node_file.close()

num_elem = mesh.NumberOfElementsGet()
elements_list = [[0 for i in range(9)] for elem in range(num_elem)]
for j in range(0,num_elem):
    elemnodes = meshElements.NodesGet(1+j,10)
    elements_list[j] = [
      elemnodes[0], elemnodes[1], elemnodes[2], 
      elemnodes[3], elemnodes[4], elemnodes[7],
      elemnodes[5], elemnodes[6], elemnodes[9], 
      elemnodes[8],
      ]

# for j in range(0,num_elem):
#     elemnodes = meshElements.NodesGet(1+j,10)
#     elements_list[j] = [
#       elemnodes[0], elemnodes[1], elemnodes[2], 
#       elemnodes[3], elemnodes[4], elemnodes[5],
#       elemnodes[6], elemnodes[7], elemnodes[9], 
#       elemnodes[8],
#       ]

elem_file=open('output_mesh.ele','w')
elem_file.writelines([str(line) + "\n" for line in elements_list])
elem_file.close()

before_def = np.array(nodes_before_def)
np.save('output_mesh_before_def.npy',before_def)
after_def = np.array(nodes_after_def)
np.save('output_mesh_after_def.npy',after_def)

# Obtain the nodeal values and their coordinates
points = before_def
point_data = {"deformed": after_def}

cells_array = np.array(elements_list)[:,:] - 1
cells = [("tetra10", cells_array)] + [("tetra10", cells_array)]

meshio.write_points_cells("output_mesh.vtk",points, cells, point_data)

meshio.write_points_cells("output_mesh_iron.vtk",points,cells, point_data)

problem.Destroy()
coordinate_system.Destroy()
region.Destroy()
basis.Destroy()
iron.Finalise()


## CHECK
def orient_posneg(nodes, node_coords):
    #x_A, y_A, z_A, x_B, y_B, z_B, x_C, y_C, z_C, x_D, y_D, z_D, x_E, y_E, z_E, x_F, y_F, z_F, x_G, y_G, z_G, x_H, y_H, z_H, x_I, y_I, z_I, x_J, y_J, z_J):
    # Define vectors AB, AC, AD, AE, AF, BC, BD, BE, BF, CD, CE, CF, DE, DF, EF
    AB = [x_B - x_A, y_B - y_A, z_B - z_A]
    AC = [x_C - x_A, y_C - y_A, z_C - z_A]
    AD = [x_D - x_A, y_D - y_A, z_D - z_A]
    AE = [x_E - x_A, y_E - y_A, z_E - z_A]
    AF = [x_F - x_A, y_F - y_A, z_F - z_A]
    BC = [x_C - x_B, y_C - y_B, z_C - z_B]
    BD = [x_D - x_B, y_D - y_B, z_D - z_B]
    BE = [x_E - x_B, y_E - y_B, z_E - z_B]
    BF = [x_F - x_B, y_F - y_B, z_F - z_B]
    CD = [x_D - x_C, y_D - y_C, z_D - z_C]
    CE = [x_E - x_C, y_E - y_C, z_E - z_C]
    CF = [x_F - x_C, y_F - y_C, z_F - z_C]
    DE = [x_E - x_D, y_E - y_D, z_E - z_D]
    DF = [x_F - x_D, y_F - y_D, z_F - z_D]
    EF = [x_F - x_E, y_F - y_E, z_F - z_E]
    
    # Calculate six determinants of submatrices formed by the vectors
    det_1 = np.linalg.det([AB, AC, AD])
    det_2 = np.linalg.det([AC, AE, AF])
    det_3 = np.linalg.det([AD, AE, AF])
    det_4 = np.linalg.det([BD, BE, BF])
    det_5 = np.linalg.det([CD, CE, CF])
    det_6 = np.linalg.det([DE, DF, EF])
    
    vectors = {"AB": 0, "AC": 0, "AD": 0, "AE": 0, "AF": 0, 
               "BC": 0, "BD": 0, "BE": 0, "BF": 0, "CD": 0, 
               "CE": 0, "CF": 0, "DE": 0, "DF": 0, "EF": 0}
    for i, vector_name in enumerate(vectors.keys()):
        vectors[vector_name] = [node_coords[nodes[i+1]-1][0] - node_coords[nodes[i]-1][0], 
            node_coords[nodes[i+1]-1][1] - node_coords[nodes[i]-1][1], 
            node_coords[nodes[i+1]-1][2] - node_coords[nodes[i]-1][2]]
    
    # Calculate six determinants of submatrices formed by the vectors
    det_1 = np.linalg.det([vectors["AB"], vectors["AC"], vectors["AD"]])
    det_2 = np.linalg.det([vectors["AC"], vectors["AE"], vectors["AF"]])
    det_3 = np.linalg.det([vectors["AD"], vectors["AE"], vectors["AF"]])
    det_4 = np.linalg.det([vectors["BD"], vectors["BE"], vectors["BF"]])
    det_5 = np.linalg.det([vectors["CD"], vectors["CE"], vectors["CF"]])
    det_6 = np.linalg.det([vectors["DE"], vectors["DF"], vectors["EF"]])
    
    # Determine orientation based on signs of the determinants
    if (det_1 > 0 and det_6 > 0) or (det_1 < 0 and det_6 < 0):
        return "Positive"
    elif (det_2 > 0 and det_4 > 0 and det_5 > 0) or (det_2 < 0 and det_4 < 0 and det_5 < 0):
        return "Positive"
    elif (det_3 > 0 and det_4 < 0 and det_5 < 0) or (det_3 < 0 and det_4 > 0 and det_5 > 0):
        return "Positive"
    else:
        return "Negative"
