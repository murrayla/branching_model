

import numpy as np

# Initialise OpenCMISS-Iron.
from opencmiss.iron import iron

# Set constants
X, Y, Z = (1, 2, 3)

# Specify the number of local element directions.
dimensions = 3
number_of_xi = dimensions

# Specify the number of elements along each element direction.
# Note that the number of elements in each direction were set to create a mesh with uniform square elements.
number_global_x_elements = 5
number_global_y_elements = 5
number_global_z_elements = 5

# Set dimensions of the cube (in mm).
width = 40
length = 40
height = 90

interpolation_type = 1 # Linear Lagrange interpolation.
number_of_gauss_per_xi = 3 # Gauss points along each local element coordinate direction (xi).
use_pressure_basis = True
number_of_load_increments = 40

# Version and derivative numbers for this specific problem as constants.
VersionNumber = 1
DerivativeNumber = 1

# Create a 3D rectangular cartesian coordinate system.
coordinate_system_user_number = 1
coordinate_system = iron.CoordinateSystem()
coordinate_system.CreateStart(coordinate_system_user_number)
coordinate_system.DimensionSet(dimensions)
coordinate_system.CreateFinish()

basis_user_number = 1
pressure_basis_user_number = 2

# Define basis parameters.
if dimensions == 1:
    xi_interpolation = [iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE]
    number_of_guass_xi = [3]
elif dimensions == 2:
    xi_interpolation = [
        iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE,
        iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE]
    number_of_guass_xi = [3, 3]
elif dimensions == 3:
    xi_interpolation = [
        iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE,
        iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE,
        iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE]
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
    [iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*number_of_xi)
if number_of_gauss_per_xi>0:
    basis.QuadratureNumberOfGaussXiSet( [number_of_gauss_per_xi]*number_of_xi)
basis.CreateFinish()

if use_pressure_basis:
    # Define hydrostatic pressure basis.
    pressure_basis = iron.Basis()
    pressure_basis.CreateStart(pressure_basis_user_number)
    if interpolation_type in (1,2,3,4):
        pressure_basis.TypeSet(iron.BasisTypes.LAGRANGE_HERMITE_TP)
    elif interpolation_type in (7,8,9):
        pressure_basis.TypeSet(iron.BasisTypes.SIMPLEX)
    pressure_basis.NumberOfXiSet(number_of_xi)
    pressure_basis.InterpolationXiSet(
        [iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE]*number_of_xi)
    if number_of_gauss_per_xi > 0:
        pressure_basis.QuadratureNumberOfGaussXiSet(
          [number_of_gauss_per_xi]*number_of_xi)
    pressure_basis.CreateFinish()

# Create a region and assign the coordinate system to the region.
region_user_number = 1
region = iron.Region()
region.CreateStart(region_user_number,iron.WorldRegion)
region.LabelSet("Region")
region.CoordinateSystemSet(coordinate_system)
region.CreateFinish()

# Start the creation of a generated mesh in the region.
generated_mesh_user_number = 1
generated_mesh = iron.GeneratedMesh()
generated_mesh.CreateStart(generated_mesh_user_number, region)
generated_mesh.TypeSet(iron.GeneratedMeshTypes.REGULAR)
if use_pressure_basis:
    generated_mesh.BasisSet([basis, pressure_basis])
else:
    generated_mesh.BasisSet([basis])
    generated_mesh.ExtentSet([width, length, height])
    generated_mesh.NumberOfElementsSet(
        [number_global_x_elements, number_global_y_elements, number_global_z_elements])

# Finish the creation of a generated mesh in the region.
mesh_user_number = 1
mesh = iron.Mesh()
generated_mesh.CreateFinish(mesh_user_number,mesh)

# Get the number of computational nodes.
number_of_computational_nodes = iron.ComputationalNumberOfNodesGet()

# Create a decomposition for the mesh.
decomposition_user_number = 1
decomposition = iron.Decomposition()
decomposition.CreateStart(decomposition_user_number,mesh)
decomposition.TypeSet(iron.DecompositionTypes.CALCULATED)
decomposition.NumberOfDomainsSet(number_of_computational_nodes)
decomposition.CreateFinish()

# Create a field for the geometry.
geometric_field_user_number = 1

geometric_field = iron.Field()
geometric_field.CreateStart(geometric_field_user_number,region)
geometric_field.MeshDecompositionSet(decomposition)
geometric_field.TypeSet(iron.FieldTypes.GEOMETRIC)
geometric_field.VariableLabelSet(iron.FieldVariableTypes.U,"Geometry")
geometric_field.ComponentMeshComponentSet(iron.FieldVariableTypes.U,1,1)
geometric_field.ComponentMeshComponentSet(iron.FieldVariableTypes.U,2,1)
geometric_field.ComponentMeshComponentSet(iron.FieldVariableTypes.U,3,1)
if interpolation_type == 4:
    # Set arc length scaling for cubic-Hermite elements.
    geometric_field.ScalingTypeSet(iron.FieldScalingTypes.ARITHMETIC_MEAN)
geometric_field.CreateFinish()

# Update the geometric field parameters from generated mesh.
generated_mesh.GeometricParametersCalculate(geometric_field)

dependent_field_user_number = 2

dependent_field = iron.Field()
dependent_field.CreateStart(dependent_field_user_number, region)
dependent_field.MeshDecompositionSet(decomposition)
dependent_field.TypeSet(iron.FieldTypes.GEOMETRIC_GENERAL)
dependent_field.GeometricFieldSet(geometric_field)
dependent_field.DependentTypeSet(iron.FieldDependentTypes.DEPENDENT)
dependent_field.VariableLabelSet(iron.FieldVariableTypes.U, "Dependent")
dependent_field.NumberOfVariablesSet(2)

# Set the number of components for the U variable (position) and the DELUDELN
# (forces).
dependent_field.NumberOfComponentsSet(iron.FieldVariableTypes.U, 4)
dependent_field.NumberOfComponentsSet(iron.FieldVariableTypes.DELUDELN, 4)

if use_pressure_basis:
    # Set the hydrostatic pressure to be nodally based and use the second mesh component.
    # U variable (position)
    dependent_field.ComponentInterpolationSet(
        iron.FieldVariableTypes.U, 4, iron.FieldInterpolationTypes.NODE_BASED)
    dependent_field.ComponentMeshComponentSet(
        iron.FieldVariableTypes.U, 4, 2)

    # DELUDELN variable (forces)
    dependent_field.ComponentInterpolationSet(
        iron.FieldVariableTypes.DELUDELN, 4,
        iron.FieldInterpolationTypes.NODE_BASED)
    dependent_field.ComponentMeshComponentSet(
        iron.FieldVariableTypes.DELUDELN, 4, 2)

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

# Density parameters
#density=0.0009 #in g mm^-3
density = 0.00102 #9345

# Create the material field.
material_field_user_number = 3
material_field = iron.Field()
material_field.CreateStart(material_field_user_number, region)
material_field.TypeSet(iron.FieldTypes.MATERIAL)
material_field.MeshDecompositionSet(decomposition)
material_field.GeometricFieldSet(geometric_field)
material_field.NumberOfVariablesSet(2)
material_field.VariableTypesSet([iron.FieldVariableTypes.U, iron.FieldVariableTypes.V])
material_field.VariableLabelSet(iron.FieldVariableTypes.U, "Material")
material_field.VariableLabelSet(iron.FieldVariableTypes.V, "Density")

# Set the number of components for the Mooney Rivlin constitutive equation (2).
material_field.NumberOfComponentsSet(iron.FieldVariableTypes.U, 2)

# Set a new 1-component variable for density.
material_field.NumberOfComponentsSet(iron.FieldVariableTypes.V, 1)

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
    iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, 2.5)
material_field.ComponentValuesInitialiseDP(
    iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 2, 0.0)
material_field.ComponentValuesInitialiseDP(
    iron.FieldVariableTypes.V,iron.FieldParameterSetTypes.VALUES, 1, density)

# Equation set field.
equations_set_field_user_number = 4
equations_set_field = iron.Field()
equations_set = iron.EquationsSet()


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

# Set up source field for gravity.
sourceFieldUserNumber = 5

# Gravity parameters.
angle = np.deg2rad(60)
gravity= [0.0, -9.81*np.sin(angle), -9.81*np.cos(angle)] #in ms^-2.
print("Gravity:", gravity)

#Create the source field with the gravity vector
sourceField = iron.Field()
equations_set.SourceCreateStart(sourceFieldUserNumber,sourceField)
if interpolation_type == 4:
    sourceField.fieldScalingType = iron.FieldScalingTypes.ARITHMETIC_MEAN
else:
    sourceField.fieldScalingType = iron.FieldScalingTypes.UNIT
equations_set.SourceCreateFinish()

#Set the gravity vector component values
sourceField.ComponentValuesInitialiseDP(
    iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,1,gravity[0])
sourceField.ComponentValuesInitialiseDP(
    iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,2,gravity[1])
sourceField.ComponentValuesInitialiseDP(
    iron.FieldVariableTypes.U,iron.FieldParameterSetTypes.VALUES,3,gravity[2])

# Finish creating equations set.
equations_set.CreateFinish() # once this is called, cant add anything else to equations_set.

# Create equations.
equations = iron.Equations()
equations_set.EquationsCreateStart(equations)
equations.SparsityTypeSet(iron.EquationsSparsityTypes.SPARSE)
equations.OutputTypeSet(iron.EquationsOutputTypes.NONE)
equations_set.EquationsCreateFinish()

# Define the problem.
problem_user_number = 1

problem = iron.Problem()
problem_specification = (
  [iron.ProblemClasses.ELASTICITY,
   iron.ProblemTypes.FINITE_ELASTICITY,
   iron.ProblemSubtypes.NONE])
problem.CreateStart(problem_user_number, problem_specification)
problem.CreateFinish()

# Create the problem control loop.
problem.ControlLoopCreateStart()
control_loop = iron.ControlLoop()
problem.ControlLoopGet([iron.ControlLoopIdentifiers.NODE], control_loop)
control_loop.MaximumIterationsSet(number_of_load_increments)
problem.ControlLoopCreateFinish()

# Number of iterations the solver runs for.
iteration_num = 1000

# Solving process.
nonlinear_solver = iron.Solver()
linear_solver = iron.Solver()
problem.SolversCreateStart()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE], 1, nonlinear_solver)
nonlinear_solver.OutputTypeSet(iron.SolverOutputTypes.NONE)
nonlinear_solver.NewtonJacobianCalculationTypeSet(
  iron.JacobianCalculationTypes.EQUATIONS)
nonlinear_solver.NewtonLinearSolverGet(linear_solver)
nonlinear_solver.NewtonMaximumIterationsSet(iteration_num)
linear_solver.LinearTypeSet(iron.LinearSolverTypes.DIRECT)
problem.SolversCreateFinish()

solver = iron.Solver()
solver_equations = iron.SolverEquations()
problem.SolverEquationsCreateStart()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE], 1, solver)
solver.SolverEquationsGet(solver_equations)
solver_equations.SparsityTypeSet(iron.SolverEquationsSparsityTypes.SPARSE)
_ = solver_equations.EquationsSetAdd(equations_set)
#solver_equations.Finalise()
problem.SolverEquationsCreateFinish()
#problem.Finalise()

nodes = iron.MeshNodes()
mesh.NodesGet(1, nodes)
num_nodes = nodes.NumberOfNodesGet()
node_nums = (np.arange(num_nodes)+1).tolist()

elements = iron.MeshElements()
mesh.ElementsGet(1, elements)
num_elements = mesh.NumberOfElementsGet()
element_nums = (np.arange(num_elements)+1).tolist()

boundary_conditions = iron.BoundaryConditions()
solver_equations.BoundaryConditionsCreateStart(boundary_conditions)

# Applying zero-displacement boundary condition on the nodes of the bottom face.
for node_id, node in enumerate(range(1,num_nodes+1)):
    value = dependent_field.ParameterSetGetNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, VersionNumber, DerivativeNumber, int(node), 3)
    if np.allclose(value, 0.0, atol = 1e-06):
        for component_idx, component in enumerate(range(X,Z+1)):
            boundary_conditions.AddNode(dependent_field, iron.FieldVariableTypes.U, VersionNumber, DerivativeNumber, int(node), component, iron.BoundaryConditionsTypes.FIXED, 0.0)
            fixed_nodes = int(node)

# Number of fixed nodes.
print(fixed_nodes)

# Initialise nodes and reference coordinate lists.
node_value = []
reference = []

# Reference configuration coordinates for all nodes.
for node_id, node in enumerate(range(1,num_nodes+1)):
    for component_idx, component in enumerate(range(X,Z+1)):
        coordinate = dependent_field.ParameterSetGetNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, VersionNumber, DerivativeNumber, int(node), component)
        if component == 1:
            x_value = coordinate
        elif component == 2:
            y_value = coordinate
        elif component == 3:
            z_value = coordinate
        node_value.append(node)
        reference.append(coordinate)

    #print(x_value, y_value, z_value)

    # Let's pick a corner node at the free end, with nodal coordinates (0,0,height).
    if (x_value, y_value, z_value) == (0, 0, float(height)):
        node_want = int(node)

# Convert list to array.
node_value = np.array(node_value)
reference = np.array(reference)

solver_equations.BoundaryConditionsCreateFinish()

# Solve the problem.
try:
    problem.Solve()
    print("Problem solved.")

except:
    print("Too few iterations.")

# Initialise nodes and deformed coordinates lists.
deformed = []

for node_id, node in enumerate(range(1,num_nodes+1)):
    for component_idx, component in enumerate(range(X,Z+1)):
        coordinate = dependent_field.ParameterSetGetNodeDP(iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, VersionNumber, DerivativeNumber, int(node), component)
        deformed.append(coordinate)

# Convert list to array.
deformed = np.array(deformed)

# Nodal displacement info.
displacement = deformed - reference
node_displacement = displacement[3*(node_want - 1):3*node_want]

# Outputting the relevant info.
print("Number of elements:", num_elements)
print("Node:", node_want)
print('Displacements (mm):', abs(node_displacement))

DoFs = 3*(num_nodes - fixed_nodes) + num_elements
print("Degrees of Freedom:", DoFs)

deformed_field_user_number = 6
deformed_field = iron.Field()
deformed_field.CreateStart(deformed_field_user_number, region)
deformed_field.MeshDecompositionSet(decomposition)
deformed_field.TypeSet(iron.FieldTypes.GEOMETRIC)
deformed_field.VariableLabelSet(iron.FieldVariableTypes.U, "DeformedGeometry")
for component in [1, 2, 3]:
    deformed_field.ComponentMeshComponentSet(
            iron.FieldVariableTypes.U, component, 1)
if interpolation_type == 4:
    # Set arc length scaling for cubic-Hermite elements.
    deformed_field.ScalingTypeSet(iron.FieldScalingTypes.ARITHMETIC_MEAN)
deformed_field.CreateFinish()

pressure_field_user_number = 7
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

# Export results to exnode and exelem format.
fields = iron.Fields()
fields.CreateRegion(region)
fields.NodesExport("cube", "FORTRAN")
fields.ElementsExport("cube", "FORTRAN")
fields.Finalise()

problem.Destroy()
coordinate_system.Destroy()
region.Destroy()
basis.Destroy()
iron.Finalise()

