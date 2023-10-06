#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy

# Intialise OpenCMISS-Iron.
from opencmiss.iron import iron


# In[2]:


# Set constants
X, Y, Z = (1, 2, 3)

# Set model number to solve (these specify different loading conditions).
model = 1

# Specify the number of local element directions.
number_of_xi = 3

number_of_dimensions = 3
interpolation_type = 8
#number_of_gauss_per_xi = 3 # Gauss points along each local element coordinate direction (xi).
use_pressure_basis = False
number_of_load_increments = 12


# In[3]:


# Create coordinate system.
coordinate_system_user_number = 1
coordinate_system = iron.CoordinateSystem()
coordinate_system.CreateStart(coordinate_system_user_number)
coordinate_system.DimensionSet(3)
coordinate_system.CreateFinish()


# In[4]:


# File name
file_name = 'gmsh2iron'

#Inputing element file
elem_file = open(file_name + '.ele', 'r')

ele_file_list = list()
element_nodes = list()
for i, line in enumerate(elem_file):
    ele_file_list.append(line)
    if ele_file_list[i].split()[0] == "second-order-tetrahedron":
        element_nodes.append(ele_file_list[i].split()[3])
        element_nodes.append(ele_file_list[i].split()[4])
        element_nodes.append(ele_file_list[i].split()[5])
        element_nodes.append(ele_file_list[i].split()[6])
        element_nodes.append(ele_file_list[i].split()[7])
        element_nodes.append(ele_file_list[i].split()[8])
        element_nodes.append(ele_file_list[i].split()[9])
        element_nodes.append(ele_file_list[i].split()[10])
        element_nodes.append(ele_file_list[i].split()[11])
        element_nodes.append(ele_file_list[i].split()[12])
    else:
        continue
        
elem_nodes_used = set(element_nodes)
# Reading node file
node_file = open(file_name + '.node', 'r')

node_file_list = list()
for i, line in enumerate(node_file):
    node_file_list.append(line)
    if (node_file_list[-1][:].split()[0] in elem_nodes_used) or (i == 0):
        continue
    else: 
        node_file_list.pop()


# In[5]:


# As formatted to gmsh documentation
# numEntityBlocks(size_t) numNodes(size_t) minNodeTag(size_t) maxNodeTag(size_t)

# Reading the mesh details from the first line of the node file
numEntityBlocks, numNodes, minNodeTag, maxNodeTag = node_file_list[0].split()
num_nodes = len(elem_nodes_used)
num_coords= 3
num_attributes=0
boundary_markers=0
# Creating variables to store node number & boundary marker
NodeNums = [[0 for m in range(2)] for n in range(num_nodes)]
# Creating array to store x and y coordinates
NodeCoords = [[0 for m in range(num_coords)] for n in range(num_nodes)]

nodeIDs = list()
# Reading details from node fileb
for i in range(num_nodes):
    NodeNums[i][0], NodeCoords[i][0], NodeCoords[i][1], NodeCoords[i][2] = node_file_list[i+1].split()  #node number, x, y, z, boundary marker
    # Converting from string to appropriate type
    NodeNums[i][0] = int(NodeNums[i][0])
    NodeCoords[i][0] = float(NodeCoords[i][0])
    NodeCoords[i][1] = float(NodeCoords[i][1])
    NodeCoords[i][2] = float(NodeCoords[i][2])
    nodeIDs.append(NodeNums[i][0])

# Closing the file
node_file.close()


# In[6]:


# gmsh .msh file format:
# numEntityBlocks(size_t) numElements(size_t) minElementTag(size_t) maxElementTag(size_t)

#Reading the values of the first line
numEntityBlocks, numElements, minElementTag, maxElementTag = ele_file_list[0].split()

# Set Element Types
# types = {1: "line", 2: "triangle", 3: "quadrangle", 4: "tetrahedron", 
#          5: "hexahedron", 6: "prism", 7: "pyramid", 8: "second-order-line", 
#          9: "second-order-triangle", 10: "second-order-quadrangle", 
#          11: "second-order-tetrahedron", 15: "point"}
types = {4: "tetrahedron", 11: "second-order-tetrahedron"}

element_type = list()
count = 0
for i, line in enumerate(ele_file_list):
    if line.split()[0] in types.values():
        element_type.append(line.split()[0])

# element_count.append(count)
element_type_count = list()
element_type_set = set(element_type)

element_type_count = map(lambda x: element_type.count(x), element_type_set)

elements_breakdown = dict(zip(element_type_set, list(element_type_count)))

ele_attributes = 0

# Condition on the different amount of mesh elements (3 dimensional)
if len(elements_breakdown.values()) >= 3:
    # element counts
    num_ele_tetra = elements_breakdown["tetrahedron"]
    num_ele_quadr = elements_breakdown["quadrangle"]
    num_ele_pyram = elements_breakdown["pyramid"]
    # node counts
    num_nodes_tetra = 4
    num_nodes_quadr = 4
    num_nodes_pyram = 5
    # Creating variable to store the element map
    ele_map_tetra = [[0 for x in range(num_nodes_tetra + ele_attributes)] for y in range(num_ele_tetra)]
    ele_map_quadr = [[0 for x in range(num_nodes_quadr + ele_attributes)] for y in range(num_ele_quadr)]
    ele_map_pyram = [[0 for x in range(num_nodes_pyram + ele_attributes)] for y in range(num_ele_pyram)]
    # Element indexes
    ele_idx_tetra = [[0 for m in range(1)] for n in range(num_ele_tetra)]
    ele_idx_quadr = [[0 for m in range(1)] for n in range(num_ele_quadr)]
    ele_idx_pyram = [[0 for m in range(1)] for n in range(num_ele_pyram)]
    # Reading element data from elemfile
    ele_idx_count_tetra = 0
    ele_idx_count_quadr = 0
    ele_idx_count_pyram = 0
    check = 0
    for line in ele_file_list:
        if line.split()[0] == "tetrahedron":
            _, _, _, ele_map_tetra[ele_idx_count_tetra][0], ele_map_tetra[ele_idx_count_tetra][1], ele_map_tetra[ele_idx_count_tetra][2], ele_map_tetra[ele_idx_count_tetra][3] = line.split()
            ele_idx_count_tetra += 1
            ele_idx_tetra[ele_idx_count_tetra-1][0] = ele_idx_count_tetra
        elif line.split()[0] == "quadrangle":
            _, _, _, ele_map_quadr[ele_idx_count_quadr][0], ele_map_quadr[ele_idx_count_quadr][1], ele_map_quadr[ele_idx_count_quadr][2], ele_map_quadr[ele_idx_count_quadr][3] = line.split()
            ele_idx_count_quadr += 1
            ele_idx_quadr[ele_idx_count_quadr-1][0] = ele_idx_count_quadr
        elif line.split()[0] == "pyramid":
            _, _, _, ele_map_pyram[ele_idx_count_pyram][0], ele_map_pyram[ele_idx_count_pyram][1], ele_map_pyram[ele_idx_count_pyram][2], ele_map_pyram[ele_idx_count_pyram][3], ele_map_pyram[ele_idx_count_pyram][4] = line.split()
            ele_idx_count_pyram += 1
            ele_idx_pyram[ele_idx_count_pyram-1][0] = ele_idx_count_pyram
    
elif len(elements_breakdown.values()) == 2:
    # element counts
    num_ele_tetra = elements_breakdown["tetrahedron"]
    num_ele_quadr = elements_breakdown["pyramid"]
    # node counts
    num_nodes_tetra = 4
    num_nodes_quadr = 4
    # Creating variable to store the element map
    ele_map_tetra = [[0 for x in range(num_nodes_tetra + ele_attributes)] for y in range(num_ele_tetra)]
    ele_map_quadr = [[0 for x in range(num_nodes_quadr + ele_attributes)] for y in range(num_ele_quadr)]
    # Element indexes
    ele_idx_tetra = [[0 for m in range(1)] for n in range(num_ele_tetra)]
    ele_idx_quadr = [[0 for m in range(1)] for n in range(num_ele_quadr)]
    # Reading element data from elemfile
    ele_idx_count_tetra = 0
    ele_idx_count_quadr = 0
    check = 0
    for line in ele_file_list:
        if line.split()[0] == "tetrahedron":
            _, _, _, ele_map_tetra[ele_idx_count_tetra][0], ele_map_tetra[ele_idx_count_tetra][1], ele_map_tetra[ele_idx_count_tetra][2], ele_map_tetra[ele_idx_count_tetra][3] = line.split()
            ele_idx_count_tetra += 1
            ele_idx_tetra[ele_idx_count_tetra-1][0] = ele_idx_count_tetra
        elif line.split()[0] == "quadrangle":
            _, _, _, ele_map_quadr[ele_idx_count_quadr][0], ele_map_quadr[ele_idx_count_quadr][1], ele_map_quadr[ele_idx_count_quadr][2], ele_map_quadr[ele_idx_count_quadr][3] = line.split()
            ele_idx_count_quadr += 1
            ele_idx_quadr[ele_idx_count_quadr-1][0] = ele_idx_count_quadr
            
if len(elements_breakdown.values()) == 1:
    num_ele_2Otetra = 0
    num_ele_tetra = 0
    print(elements_breakdown)
    # quadrangle
    if "second-order-tetrahedron" in elements_breakdown.keys():
        # element counts
        num_ele_2Otetra = elements_breakdown["second-order-tetrahedron"]
        # node counts
        num_nodes_2Otetra = 10
        # Creating variable to store the element map
        ele_map_2Otetra = [[0 for x in range(num_nodes_2Otetra + ele_attributes)] for y in range(num_ele_2Otetra)]
        # Element indexes
        ele_idx_2Otetra = [[0 for m in range(1)] for n in range(num_ele_2Otetra)]
        # Reading element data from elemfile
        ele_idx_count = 0
        check = 0
        for line in ele_file_list:
            if line.split()[0] == "second-order-tetrahedron":
                _, _, _, \
                    ele_map_2Otetra[ele_idx_count][0], ele_map_2Otetra[ele_idx_count][1], ele_map_2Otetra[ele_idx_count][2], \
                        ele_map_2Otetra[ele_idx_count][3],  ele_map_2Otetra[ele_idx_count][4], ele_map_2Otetra[ele_idx_count][5], \
                            ele_map_2Otetra[ele_idx_count][6],  ele_map_2Otetra[ele_idx_count][7], ele_map_2Otetra[ele_idx_count][8], \
                                ele_map_2Otetra[ele_idx_count][9] = line.split()
                ele_idx_count += 1
                ele_idx_2Otetra[ele_idx_count-1][0] = ele_idx_count

    # tetrahedral
    elif "tetrahedron" in elements_breakdown.keys():
        # element counts
        num_ele_tetra = elements_breakdown["tetrahedron"]
        # node counts
        num_nodes_tetra = 4
        # Creating variable to store the element map
        ele_map_tetra = [[0 for x in range(num_nodes_tetra + ele_attributes)] for y in range(num_ele_tetra)]
        # Element indexes
        ele_idx_tetra = [[0 for m in range(1)] for n in range(num_ele_tetra)]
        # Reading element data from elemfile
        ele_idx_count = 0
        check = 0
        for line in ele_file_list:
            if line.split()[0] == "tetrahedron":
                _, _, _, ele_map_tetra[ele_idx_count][0], ele_map_tetra[ele_idx_count][1], ele_map_tetra[ele_idx_count][2], ele_map_tetra[ele_idx_count][3] = line.split()
                ele_idx_count += 1
                ele_idx_tetra[ele_idx_count-1][0] = ele_idx_count

# # converting 2d list into 1d
# num_tetra_nodes = len(set(chain.from_iterable(ele_map_tetra)))

# Closing the file
elem_file.close()


# In[7]:

basis_user_number = 1
pressure_basis_user_number = 2

if num_ele_2Otetra > 0:
    xi_interpolation = [
        iron.BasisInterpolationSpecifications.QUADRATIC_SIMPLEX,
        iron.BasisInterpolationSpecifications.QUADRATIC_SIMPLEX,
        iron.BasisInterpolationSpecifications.QUADRATIC_SIMPLEX]
    number_of_guass_xi = [3, 3, 3]

    # Define geometric basis.
    basis = iron.Basis()
    basis.CreateStart(basis_user_number)
    if interpolation_type in (1,2,3,4):
        basis.TypeSet(iron.BasisTypes.LAGRANGE_HERMITE_TP)
    elif interpolation_type in (7,8,9):
        basis.TypeSet(iron.BasisTypes.SIMPLEX)
    basis.NumberOfXiSet(number_of_xi)
    basis.InterpolationXiSet(
        [iron.BasisInterpolationSpecifications.QUADRATIC_SIMPLEX]*number_of_xi)
    # if number_of_gauss_per_xi>0:
    #     basis.QuadratureNumberOfGaussXiSet( [number_of_gauss_per_xi]*number_of_xi)
    basis.QuadratureOrderSet(5)

elif num_ele_tetra > 0:
    xi_interpolation = [
            iron.BasisInterpolationSpecifications.LINEAR_SIMPLEX,
            iron.BasisInterpolationSpecifications.LINEAR_SIMPLEX,
            iron.BasisInterpolationSpecifications.LINEAR_SIMPLEX]
    number_of_guass_xi = [3, 3, 3]

    # Define geometric basis.
    basis = iron.Basis()
    basis.CreateStart(basis_user_number)
    basis.TypeSet(iron.BasisTypes.SIMPLEX)
    basis.NumberOfXiSet(number_of_xi)
    basis.InterpolationXiSet(xi_interpolation)
    # basis.QuadratureNumberOfGaussXiSet(number_of_guass_xi)

basis.CreateFinish()


# In[8]:


# Create region
region_user_number = 1
region = iron.Region()
region.CreateStart(region_user_number, iron.WorldRegion)
region.CoordinateSystemSet(coordinate_system)
region.LabelSet("Region")
region.CreateFinish()


# In[9]:

# Initialise Nodes
nodes = iron.Nodes()
nodes.CreateStart(region, num_nodes)
nodes.CreateFinish()


# In[10]:


# Initialise Mesh
mesh = iron.Mesh()
mesh_user_number=1
mesh.CreateStart(mesh_user_number, region, num_coords)
mesh.NumberOfElementsSet(num_ele_2Otetra)
mesh.NumberOfComponentsSet(1)


# In[11]:


# Initialise Elements
meshElements = iron.MeshElements()
meshElements.CreateStart(mesh, 1, basis)


# In[12]:


# ChatGPT reordering function
def get_quadratic_simplex_node_numbers(labels, coordinates):
    # Sort nodes by x-coordinate, breaking ties with y- and z-coordinates
    nodes = [(label, coordinates[label-1]) for label in labels]
    nodes.sort(key=lambda n: (n[1][0], n[1][1], n[1][2]))
    # Find first node as vertex with smallest x-coordinate
    vertex_label, vertex_coords = nodes[0]
    vertex_index = labels.index(vertex_label) + 1
    
    # Find remaining nodes in lexicographic order
    remaining_nodes = [(label, coords) for label, coords in nodes if label != vertex_label]
    edge_nodes = sorted(remaining_nodes, key=lambda n: (n[1][1], n[1][2], n[1][0]))
    edge_node_indices = [labels.index(label) + 1 for label, coords in edge_nodes]
    # Construct final list of node numbers
    node_numbers = [vertex_index] + edge_node_indices
    return node_numbers



# In[12]:
if num_ele_2Otetra > 0:

    for i in range(num_ele_2Otetra):
        element = ele_idx_2Otetra[i][0]
        
#         try:
        nodesList = list(
        map(int,[ele_map_2Otetra[i][0], ele_map_2Otetra[i][1], ele_map_2Otetra[i][2], 
                 ele_map_2Otetra[i][3], ele_map_2Otetra[i][4], ele_map_2Otetra[i][5],
                 ele_map_2Otetra[i][6], ele_map_2Otetra[i][7], ele_map_2Otetra[i][8],
                 ele_map_2Otetra[i][9]]))
        # reorder the nodes to convert from .msh to iron standard
        new_order_indexes = get_quadratic_simplex_node_numbers(nodesList, NodeCoords)
        iron_order = [nodesList[i-1] for i in new_order_indexes]
        meshElements.NodesSet(int(element), iron_order)
#         except:
#             print("error at {}".format(i+638))
#             print("with nodesList {}".format(nodesList))
        
elif num_ele_tetra > 0:

    for i in range(num_ele_tetra):
        element = ele_idx_tetra[i][0]
        
        try: 
            nodesList = list(
            map(int,[ele_map_tetra[i][0], ele_map_tetra[i][1], ele_map_tetra[i][2], ele_map_tetra[i][3]]))
            meshElements.NodesSet(int(element), nodesList)   
        except:
            print("error at {}".format(i+638))
            print("with nodesList {}".format(nodesList))
            print((nodesList[0] in nodeIDs), (nodesList[1] in nodeIDs),(nodesList[2] in nodeIDs),(nodesList[3] in nodeIDs),)
        
meshElements.CreateFinish()


# In[13]:


# Finalise Mesh
mesh.CreateFinish()


# In[14]:


# Get the number of computational nodes.
number_of_computational_nodes = iron.ComputationalNumberOfNodesGet()

# Perform mesh decomposition.
decomposition_user_number = 1
decomposition = iron.Decomposition()
decomposition.CreateStart(decomposition_user_number, mesh)

decomposition.TypeSet(iron.DecompositionTypes.CALCULATED)
decomposition.NumberOfDomainsSet(number_of_computational_nodes)
decomposition.CreateFinish()


# In[15]:


# Create a field for the geometry.
geometric_field_user_number = 1
geometric_field = iron.Field()
geometric_field.CreateStart(geometric_field_user_number, region) #notice that the geometric field is associated with region in this function call.
geometric_field.MeshDecompositionSet(decomposition)
geometric_field.TypeSet(iron.FieldTypes.GEOMETRIC)
geometric_field.VariableLabelSet(iron.FieldVariableTypes.U,"Geometry")
geometric_field.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,1)
geometric_field.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,1)
geometric_field.ComponentMeshComponentSet(iron.FieldVariableTypes.U,3,1)
geometric_field.CreateFinish()


# In[16]:


# Set geometric field values from customized mesh.
computationalNodeNumber = iron.ComputationalNodeNumberGet()
for node_idx in range(num_nodes):
    try:
        node_id = NodeNums[node_idx][0]
        node_x = NodeCoords[node_idx][0]
        node_y = NodeCoords[node_idx][1]
        node_z = NodeCoords[node_idx][2]
        geometric_field.ParameterSetUpdateNodeDP(
          iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES,
          1, 1, node_id, 1, node_x)
        geometric_field.ParameterSetUpdateNodeDP(
          iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES,
          1, 1, node_id, 2, node_y)
        geometric_field.ParameterSetUpdateNodeDP(
          iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES,
          1, 1, node_id, 3, node_z)
    except:
        print("error at node_idx {}".format(node_idx))
        print("with node coordinates {}".format([node_x, node_y, node_z]))


# In[17]:


# Update the geometric field.
geometric_field.ParameterSetUpdateStart(
  iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)
geometric_field.ParameterSetUpdateFinish(
  iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES)


# ## /////////////////

# In[18]:


dependent_field_user_number = 2

dependent_field = iron.Field()
dependent_field.CreateStart(dependent_field_user_number, region)
dependent_field.MeshDecompositionSet(decomposition)
dependent_field.TypeSet(iron.FieldTypes.GEOMETRIC_GENERAL)
dependent_field.GeometricFieldSet(geometric_field)
dependent_field.DependentTypeSet(iron.FieldDependentTypes.DEPENDENT)
dependent_field.VariableLabelSet(iron.FieldVariableTypes.U, "Dependent")
dependent_field.NumberOfVariablesSet(2)

# Set the number of componets for the U variable (position) and the DELUDELN
# (forces).
dependent_field.NumberOfComponentsSet(iron.FieldVariableTypes.U, 4)
dependent_field.NumberOfComponentsSet(iron.FieldVariableTypes.DELUDELN, 4)

if use_pressure_basis:
    # Set the hydrostatic pressure to be nodally based and use the second mesh component.
    # U variable (position)
    dependent_field.ComponentInterpolationSet(
        iron.FieldVariableTypes.U, 4, iron.FieldInterpolationTypes.GAUSS_POINT_BASED)
    dependent_field.ComponentMeshComponentSet(
        iron.FieldVariableTypes.U, 4, 1)

    # DELUDELN variable (forces)
    dependent_field.ComponentInterpolationSet(
        iron.FieldVariableTypes.DELUDELN, 4,
        iron.FieldInterpolationTypes.NODE_BASED)
    dependent_field.ComponentMeshComponentSet(
        iron.FieldVariableTypes.DELUDELN, 4, 1)

    if interpolation_type == 4:
        # Set arc length scaling for cubic-Hermite elements.
        dependent_field.FieldScalingTypeSet(
            iron.FieldScalingTypes.ARITHMETIC_MEAN)
else:
    # Set the hydrostatic pressure to be constant within each element.
    dependent_field.ComponentInterpolationSet(
        iron.FieldVariableTypes.U, 4,
        iron.FieldInterpolationTypes.ELEMENT_BASED)
    dependent_field.ComponentInterpolationSet(
        iron.FieldVariableTypes.DELUDELN, 4,
        iron.FieldInterpolationTypes.ELEMENT_BASED)
dependent_field.CreateFinish()


# In[19]:


# Initialise dependent field from undeformed geometry and displacement bcs and set hydrostatic pressure.
iron.Field.ParametersToFieldParametersComponentCopy(
    geometric_field, iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1,
    dependent_field, iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1)
iron.Field.ParametersToFieldParametersComponentCopy(
    geometric_field, iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 2,
    dependent_field, iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 2)
iron.Field.ParametersToFieldParametersComponentCopy(
    geometric_field, iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 3,
    dependent_field, iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 3)
iron.Field.ComponentValuesInitialiseDP(
    dependent_field, iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 4, 0.0)


# In[20]:


# Create the material field.
material_field_user_number = 3

material_field = iron.Field()
material_field.CreateStart(material_field_user_number, region)
material_field.TypeSet(iron.FieldTypes.MATERIAL)
material_field.MeshDecompositionSet(decomposition)
material_field.GeometricFieldSet(geometric_field)
material_field.VariableLabelSet(iron.FieldVariableTypes.U, "Material")

# Set the number of components for the Mooney Rivlin constitutive equation (2).
material_field.NumberOfComponentsSet(iron.FieldVariableTypes.U, 2)

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
    iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, 5.0)
material_field.ComponentValuesInitialiseDP(
    iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 2, 3.0)


# In[21]:


# Equation set field.
equations_set_field_user_number = 4
equations_set_field = iron.Field()
equations_set = iron.EquationsSet()


# In[22]:


equations_set_user_number = 1

# Finite elasticity equation specification.
equations_set_specification = [iron.ProblemClasses.ELASTICITY,
    iron.ProblemTypes.FINITE_ELASTICITY,
    iron.EquationsSetSubtypes.MOONEY_RIVLIN]

# Add the geometric field and equations set field that we created earlier (note
# that while we defined the geometric field above, we only initialised an empty
# field for the equations_set_field. When an empty field is provided to the 
# equations_set, it will automatically populate it with default values).
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


# In[23]:


# Create equations.
equations = iron.Equations()
equations_set.EquationsCreateStart(equations)
equations.SparsityTypeSet(iron.EquationsSparsityTypes.SPARSE)
equations.OutputTypeSet(iron.EquationsOutputTypes.NONE)
equations_set.EquationsCreateFinish()


# In[24]:


# Define the problem.
problem_user_number = 1

problem = iron.Problem()
problem_specification = (
  [iron.ProblemClasses.ELASTICITY,
   iron.ProblemTypes.FINITE_ELASTICITY,
   iron.ProblemSubtypes.NONE])
problem.CreateStart(problem_user_number, problem_specification)
problem.CreateFinish()


# In[25]:


# Create the problem control loop.
problem.ControlLoopCreateStart()
control_loop = iron.ControlLoop()
problem.ControlLoopGet([iron.ControlLoopIdentifiers.NODE], control_loop)
control_loop.MaximumIterationsSet(number_of_load_increments)
problem.ControlLoopCreateFinish()


# In[26]:


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


# In[27]:


solver = iron.Solver()
solver_equations = iron.SolverEquations()
problem.SolverEquationsCreateStart()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE], 1, solver)
solver.SolverEquationsGet(solver_equations)
solver_equations.SparsityTypeSet(iron.SolverEquationsSparsityTypes.SPARSE)
_ = solver_equations.EquationsSetAdd(equations_set)
problem.SolverEquationsCreateFinish()


# In[28]:

# Prescribe boundary conditions (absolute nodal parameters).
boundary_conditions = iron.BoundaryConditions()
solver_equations.BoundaryConditionsCreateStart(boundary_conditions)


# Find all nodes with z = 0 or z = 1
node_idx_z0 = list()
node_idx_z1 = list()
for i, coords in enumerate(NodeCoords):
    if coords[2] == 0:
        node_idx_z0.append(i+1)
        boundary_conditions.AddNode(dependent_field, iron.FieldVariableTypes.U,\
                                    1,1,i+1,Z,iron.BoundaryConditionsTypes.FIXED,0.0)
        if (coords[1] == 0) or (coords[1] == 0.5):
            boundary_conditions.AddNode(dependent_field, iron.FieldVariableTypes.U,\
                                    1,1,i+1,Y,iron.BoundaryConditionsTypes.FIXED,0.0)
        if (coords[0] == 0) or (coords[0] == 0.5):
            boundary_conditions.AddNode(dependent_field, iron.FieldVariableTypes.U,\
                                    1,1,i+1,X,iron.BoundaryConditionsTypes.FIXED,0.0)
    elif coords[2] == 1:
        node_idx_z1.append(i+1)
        boundary_conditions.AddNode(dependent_field, iron.FieldVariableTypes.U,\
                                    1,1,i+1,Z,iron.BoundaryConditionsTypes.FIXED,0.1)
        if (coords[1] == 0) or (coords[1] == 0.5):
            boundary_conditions.AddNode(dependent_field, iron.FieldVariableTypes.U,\
                                    1,1,i+1,Y,iron.BoundaryConditionsTypes.FIXED,0.0)
        if (coords[0] == 0) or (coords[0] == 0.5):
            boundary_conditions.AddNode(dependent_field, iron.FieldVariableTypes.U,\
                                    1,1,i+1,X,iron.BoundaryConditionsTypes.FIXED,0.0)
    if (coords[0] == 0.125) and (coords[1] == 0.125):
        boundary_conditions.AddNode(dependent_field, iron.FieldVariableTypes.U,\
                                    1,1,i+1,X,iron.BoundaryConditionsTypes.FIXED,0.0)
        boundary_conditions.AddNode(dependent_field, iron.FieldVariableTypes.U,\
                                    1,1,i+1,Y,iron.BoundaryConditionsTypes.FIXED,0.0)
        


# if model == 1:
#     boundary_conditions.AddNode(dependent_field, iron.FieldVariableTypes.U,1,1,1,X,iron.BoundaryConditionsTypes.FIXED,0.0)
#     # boundary_conditions.AddNode(dependent_field, iron.FieldVariableTypes.U,1,1,3,X,iron.BoundaryConditionsTypes.FIXED,0.0)
#     # # boundary_conditions.AddNode(dependent_field, iron.FieldVariableTypes.U,1,1,5,X,iron.BoundaryConditionsTypes.FIXED,0.0)
#     # boundary_conditions.AddNode(dependent_field, iron.FieldVariableTypes.U,1,1,7,X,iron.BoundaryConditionsTypes.FIXED,0.0)
#     # boundary_conditions.AddNode(dependent_field, iron.FieldVariableTypes.U,1,1,2,X,iron.BoundaryConditionsTypes.FIXED,0.5)
#     # boundary_conditions.AddNode(dependent_field, iron.FieldVariableTypes.U,1,1,4,X,iron.BoundaryConditionsTypes.FIXED,0.5)
#     # boundary_conditions.AddNode(dependent_field, iron.FieldVariableTypes.U,1,1,6,X,iron.BoundaryConditionsTypes.FIXED,0.5)
#     # boundary_conditions.AddNode(dependent_field, iron.FieldVariableTypes.U,1,1,8,X,iron.BoundaryConditionsTypes.FIXED,0.5)

#     boundary_conditions.AddNode(dependent_field, iron.FieldVariableTypes.U,1,1,1,Y,iron.BoundaryConditionsTypes.FIXED,0.0)
#     # boundary_conditions.AddNode(dependent_field, iron.FieldVariableTypes.U,1,1,2,Y,iron.BoundaryConditionsTypes.FIXED,0.0)
#     # boundary_conditions.AddNode(dependent_field, iron.FieldVariableTypes.U,1,1,5,Y,iron.BoundaryConditionsTypes.FIXED,0.0)
#     # boundary_conditions.AddNode(dependent_field, iron.FieldVariableTypes.U,1,1,6,Y,iron.BoundaryConditionsTypes.FIXED,0.0)

#     boundary_conditions.AddNode(dependent_field, iron.FieldVariableTypes.U,1,1,1,Z,iron.BoundaryConditionsTypes.FIXED,0.0)
#     # boundary_conditions.AddNode(dependent_field, iron.FieldVariableTypes.U,1,1,2,Z,iron.BoundaryConditionsTypes.FIXED,0.0)
#     # boundary_conditions.AddNode(dependent_field, iron.FieldVariableTypes.U,1,1,3,Z,iron.BoundaryConditionsTypes.FIXED,0.0)
#     # boundary_conditions.AddNode(dependent_field, iron.FieldVariableTypes.U,1,1,4,Z,iron.BoundaryConditionsTypes.FIXED,0.0)

# elif model == 2:
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,1,X,iron.BoundaryConditionsTypes.FIXED,0.0)
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,3,X,iron.BoundaryConditionsTypes.FIXED,0.0)
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,5,X,iron.BoundaryConditionsTypes.FIXED,0.0)
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,7,X,iron.BoundaryConditionsTypes.FIXED,0.0)
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,2,X,iron.BoundaryConditionsTypes.FIXED,0.25)
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,4,X,iron.BoundaryConditionsTypes.FIXED,0.25)
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,6,X,iron.BoundaryConditionsTypes.FIXED,0.25)
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,8,X,iron.BoundaryConditionsTypes.FIXED,0.25)

#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,1,Y,iron.BoundaryConditionsTypes.FIXED,0.0)
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,2,Y,iron.BoundaryConditionsTypes.FIXED,0.0)
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,5,Y,iron.BoundaryConditionsTypes.FIXED,0.0)
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,6,Y,iron.BoundaryConditionsTypes.FIXED,0.0)
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,3,Y,iron.BoundaryConditionsTypes.FIXED,0.25)
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,4,Y,iron.BoundaryConditionsTypes.FIXED,0.25)
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,7,Y,iron.BoundaryConditionsTypes.FIXED,0.25)
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,8,Y,iron.BoundaryConditionsTypes.FIXED,0.25)

#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,1,Z,iron.BoundaryConditionsTypes.FIXED,0.0)
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,2,Z,iron.BoundaryConditionsTypes.FIXED,0.0)
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,3,Z,iron.BoundaryConditionsTypes.FIXED,0.0)
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,4,Z,iron.BoundaryConditionsTypes.FIXED,0.0)

# elif model == 3:
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,1,X,iron.BoundaryConditionsTypes.FIXED,0.0)
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,3,X,iron.BoundaryConditionsTypes.FIXED,0.0)
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,5,X,iron.BoundaryConditionsTypes.FIXED,0.5)
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,7,X,iron.BoundaryConditionsTypes.FIXED,0.5)
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,2,X,iron.BoundaryConditionsTypes.FIXED,0.0)
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,4,X,iron.BoundaryConditionsTypes.FIXED,0.0)
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,6,X,iron.BoundaryConditionsTypes.FIXED,0.5)
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,8,X,iron.BoundaryConditionsTypes.FIXED,0.5)

#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,1,Y,iron.BoundaryConditionsTypes.FIXED,0.0)
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,2,Y,iron.BoundaryConditionsTypes.FIXED,0.0)
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,5,Y,iron.BoundaryConditionsTypes.FIXED,0.0)
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,6,Y,iron.BoundaryConditionsTypes.FIXED,0.0)

#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,1,Z,iron.BoundaryConditionsTypes.FIXED,0.0)
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,2,Z,iron.BoundaryConditionsTypes.FIXED,0.0)
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,3,Z,iron.BoundaryConditionsTypes.FIXED,0.0)
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,4,Z,iron.BoundaryConditionsTypes.FIXED,0.0)
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,5,Z,iron.BoundaryConditionsTypes.FIXED,0.0)
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,6,Z,iron.BoundaryConditionsTypes.FIXED,0.0)
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,7,Z,iron.BoundaryConditionsTypes.FIXED,0.0)
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,8,Z,iron.BoundaryConditionsTypes.FIXED,0.0)

# elif model == 4:
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,1,X,iron.BoundaryConditionsTypes.FIXED,0.0)
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,3,X,iron.BoundaryConditionsTypes.FIXED,0.0)
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,6,X,iron.BoundaryConditionsTypes.FIXED,0.5)
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,8,X,iron.BoundaryConditionsTypes.FIXED,0.5)

#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,1,Y,iron.BoundaryConditionsTypes.FIXED,0.0)
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,2,Y,iron.BoundaryConditionsTypes.FIXED,0.0)
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,5,Y,iron.BoundaryConditionsTypes.FIXED,0.0)
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,6,Y,iron.BoundaryConditionsTypes.FIXED,0.0)

#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,1,Z,iron.BoundaryConditionsTypes.FIXED,0.0)
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,3,Z,iron.BoundaryConditionsTypes.FIXED,0.0)
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,6,Z,iron.BoundaryConditionsTypes.FIXED,0.5)
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,8,Z,iron.BoundaryConditionsTypes.FIXED,0.5)

# elif model == 5:
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,1,X,iron.BoundaryConditionsTypes.FIXED,0.0)
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,3,X,iron.BoundaryConditionsTypes.FIXED,0.0)
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,5,X,iron.BoundaryConditionsTypes.FIXED,0.5)
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,7,X,iron.BoundaryConditionsTypes.FIXED,0.5)
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,2,X,iron.BoundaryConditionsTypes.FIXED,0.25)
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,4,X,iron.BoundaryConditionsTypes.FIXED,0.25)
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,6,X,iron.BoundaryConditionsTypes.FIXED,0.75)
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,8,X,iron.BoundaryConditionsTypes.FIXED,0.75)

#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,1,Y,iron.BoundaryConditionsTypes.FIXED,0.0)
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,2,Y,iron.BoundaryConditionsTypes.FIXED,0.0)
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,5,Y,iron.BoundaryConditionsTypes.FIXED,0.0)
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,6,Y,iron.BoundaryConditionsTypes.FIXED,0.0)

#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,1,Z,iron.BoundaryConditionsTypes.FIXED,0.0)
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,2,Z,iron.BoundaryConditionsTypes.FIXED,0.0)
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,3,Z,iron.BoundaryConditionsTypes.FIXED,0.0)
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,4,Z,iron.BoundaryConditionsTypes.FIXED,0.0)
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,5,Z,iron.BoundaryConditionsTypes.FIXED,0.0)
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,6,Z,iron.BoundaryConditionsTypes.FIXED,0.0)
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,7,Z,iron.BoundaryConditionsTypes.FIXED,0.0)
#     boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,8,Z,iron.BoundaryConditionsTypes.FIXED,0.0)


# In[29]:


solver_equations.BoundaryConditionsCreateFinish()


# In[30]:


# Solve the problem.
problem.Solve()


# In[ ]:


problem.Destroy()
coordinate_system.Destroy()
region.Destroy()
basis.Destroy()
iron.Finalise()


# In[ ]:




