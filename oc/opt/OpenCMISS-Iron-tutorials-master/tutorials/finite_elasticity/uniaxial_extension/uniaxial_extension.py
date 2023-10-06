#!/usr/bin/env python
# coding: utf-8

# # Finite elasticity - uniaxial extension
# 
# ## Introduction
# 
# This tutorial demonstrates how to setup and solve a nonlinear continuum mechanics problem using OpenCMISS-Iron in python. For the purpose of this tutorial, we will be solving a 3D finite elasticity problem. This tutorial has been set up to solve the deformation of an isotropic, unit cube under a range of loading conditions.
# 
# See the [OpenCMISS-Iron tutorial documenation page](https://opencmiss-iron-tutorials.readthedocs.io/en/latest/tutorials.html#how-to-run-tutorials) for instructions on how to run this tutorial.
# 
# ## Learning outcomes
# At the end of this tutorial you will:
# 
# - Know the steps involved in setting up and solving a nonlinear finite elasticity simulation with OpenCMISS-Iron.
# 
# - Know how different mechanical loads can be applied to a geometry.
# 
# - Know how information about the deformation can be extracted for analysis.
# 
# ## Problem summary
# 
# ### The finite elasticity stress equilibrium equation
# In this example we are solving the stress equilibrium equation for nonlinear finite elasticity. The stress equilibrium equation represents the principle of linear momentum as follows:
# 
# $$\displaystyle \nabla \sigma + \rho \mathbf{b}-\rho \mathbf {a}=0. \qquad \text{in} \qquad  \Omega $$
# 
# where $\sigma$ is the Cauchy stress tensor, $\mathbf{b}$ is the body force vector and $\mathbf{a}$ is the vector representing the acceleration due to any unbalanced forces.
# 
# The boundary equations for the stress equilibrium equation are partitioned into the Dirichlet boundary conditions representing a fixed displacement $u_d$ over the boundary $\Gamma_d$, and the Neumann boundary conditions representing the traction forces $\mathbf{t}$ applied on the boundary $\Gamma_t$ along the normal $\mathbf{n}$ as follows:
# 
# $$
# \begin{aligned}
# \displaystyle u &= u_d \quad &\text{on} \quad \Gamma_d \\
# \displaystyle \sigma^T \mathbf{n} &= \mathbf{t} \quad &\text{on} \quad \Gamma_t 
# \end{aligned}
# $$
# 
# In this example we will solve the stress equilibrium equation for a incompressible material, which we shall perform using **the incompressibility constraint**:
# 
# $$\displaystyle (I_3-1) = 0$$
# 
# We will also use a Mooney Rivlin constitutive equation to describe the behaviour of the material.
# 
# ### Solution variables
# The problem we are about to solve includes 4 dependent variables. The first three variables represent each of the 3D coordinates $(x,y,z)$ of the mesh nodes in the deformed state. The incompressibility constraint is satisfied using lagrange multipliers, which are represented as a scalar variable, $p$. This variable is often referred to as the **hydrostatic pressure**.
# 
# ### Constitutive relationship
# Another important set of equations that need to be included when solving a mechanics problem are the stress-strain relationships or **constitutive relationships**. There are a set of constitutive equations that have already been included within the OpenCMISS-Iron library as [EquationSet subtypes](http://cmiss.bioeng.auckland.ac.nz/OpenCMISS/doc/user/group___o_p_e_n_c_m_i_s_s___equations_set_subtypes.html). We will demonstrate how these constants are used to incorporate the constitutive equation in the simulation  when setting up the equation set.
# 
# ### Loading conditions
# We will see in this tutorial how five different types of mechanical loads can be applied as follows:  
# 
# * Model 1 (Uniaxial extension of a unit cube)
# * Model 2 (Equibiaxial extension of a unit cube)
# * Model 3 (Simple shear of a unit cube)
# * Model 4 (Shear of a unit cube)
# * Model 5 (Extension and shear of a unit cube)
# 

# ## Setup
# 
# ### Loading the OpenCMISS-Iron library
# In order to use OpenCMISS-Iron we have to first import the opencmiss.iron module from the OpenCMISS-Iron package.

# In[1]:


import numpy

# Intialise OpenCMISS-Iron.
from opencmiss.iron import iron


# ### Set up parameters
# We first specify a set of variables that will be used throughout the tutorial. This is a common practice among experienced users to set these variables up at the start because we can then change these easily and re-run the rest of the code as required.

# In[2]:


# Set constants
X, Y, Z = (1, 2, 3)

# Set model number to solve (these specify different loading conditions).
model = 1

# Specify the number of local element directions.
number_of_xi = 3

# Specify the number of elements along each element direction.
number_global_x_elements = 1
number_global_y_elements = 1
number_global_z_elements = 1

# Set dimensions of the cube.
width = 1.0
length = 1.0
height = 1.0

number_of_dimensions = 3
interpolation_type = 1
number_of_gauss_per_xi = 2 # Gauss points along each local element coordinate direction (xi).
use_pressure_basis = True
number_of_load_increments = 1


# - `interpolation_type` is an integer variable that can be changed to choose one of the nine basis interpolation types defined in OpenCMISS-Iron [Basis Interpolation Specifications Constants](http://opencmiss.org/documentation/apidoc/iron/latest/python/classiron_1_1_basis_interpolation_specifications.html)
# 
# - `number_of_gauss_per_xi` is the number of Gauss points used along a given local element direction for integrating the equilibrium equations. Note that there is a minimum number of gauss points required for a given interpolation scheme (e.g. linear Lagrange interpolation requires at least 2 Gauss points along each local element direction). Using less than the minimum number will result in the solver not converging.
# 
# - `use_pressure_basis` is a boolean variable that we can set to true to set up node-based interpolation scheme for the incompressibility constraint or set to false if we want the pressure basis to be only varying between elements. This choice is typically determined by the basis function used to interpolate the geometric variables. The pressure basis is typically one order lower than the geometric and displacement field interpolation schemes. If you choose linear lagrange for displacements, then the lower order is a constant, element based interoplation. If you choose quadratic lagrange for the geometry then a linear interpolation can chosen used.  
# 
# - `number_of_load_increments` sets the number of load steps to be taken to solve the nonlinear mechanics problem. This concept is briefly explained below.
# 
#   The weak form finite element matrix equations for the stress equilibrium equations above have been implemented in the OpenCMISS-Iron libraries. Numerical treatment of the equations will show that the equations are nonlinear and, specifically, the determination of the deformed configuration coordinates turn out to be like a root finding exercise. We need to find the deformed coordinates $\mathbf{x}$ such that $f(x)=0$. These equations are therefore solved using a nonlinear iterative solver, **Newton Raphson technique** to be exact. Nonlinear solvers require an initial guess at the root of the equation and if we are far away from the actual solution then the solvers can diverge and give spurious $x$ values. Therefore, in most simulations it is best to split an entire mechanical load into lots of smaller loads. This is what `number_of_load_increments` allows us to do. You will see it come to use in the control loop section of the tutorial below.

# The next section describes how we can interact with the OpenCMISS-Iron library through it's object-oriented API to create and solve the mechanics problem.
# 
# ## Step by step guide
# 
# ### 1. Creating a coordinate system
# 
# First we construct a coordinate system that will be used to describe the geometry in our problem. The 3D geometry will exist in a 3D space, so we need a 3D coordinate system.

# In[3]:


# Create a 3D rectangular cartesian coordinate system.
coordinate_system_user_number = 1
coordinate_system = iron.CoordinateSystem()
coordinate_system.CreateStart(coordinate_system_user_number)
coordinate_system.DimensionSet(3)
coordinate_system.CreateFinish()


# ### 2. Creating basis functions
# 
# The finite element description of our fields requires a basis function to interpolate field values over elements. In OpenCMISS-Iron, we start by initalising a `basis` object on which can specific an interpolation scheme. In the following section, you can choose from a number of basis function interpolation types. These are defined in the OpenCMISS-Iron [Basis Interpolation Specifications Constants](http://opencmiss.org/documentation/apidoc/iron/latest/python/classiron_1_1_basis_interpolation_specifications.html). The `interpolation_type` variable that was set at the top of the tutorial will be used to determine which basis interpolation will be used in this simulation. Note that you will need to ensure that an appropriate `number_of_gauss_per_xi` have been set if you change the basis interpolation.
# 
# We have also initialised a second `pressure_basis` object for the hydrostatic pressure. In mechanics theory, it is generally well understood that the interpolation scheme for the hydrostatic pressure basis should be one order lower than the geometric basis. `use_pressure_basis`, which is set at the top of the tutorial, determines whether the incompressibility constraint equation will be interpolated using a nodally based interpolation scheme or interpolated as a constant within elements. If the pressure basis is not defined then the simulation is set up with element based interpolation. 
# 
# Note that if you want to describe nearly incompressible or compressible materials then you need to choose a different equation subtype that describes a constitutive equation for compressible/nearly incompressible materials. 

# In[4]:


basis_user_number = 1
pressure_basis_user_number = 2

xi_interpolation = [
        iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE,
        iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE,
        iron.BasisInterpolationSpecifications.LINEAR_LAGRANGE]
number_of_guass_xi = [3, 3, 3]

# Define geometric basis.
basis = iron.Basis()
basis.CreateStart(basis_user_number)
basis.TypeSet(iron.BasisTypes.LAGRANGE_HERMITE_TP)
basis.NumberOfXiSet(number_of_xi)
basis.InterpolationXiSet(xi_interpolation)
basis.QuadratureNumberOfGaussXiSet(number_of_guass_xi)
basis.CreateFinish()

if use_pressure_basis:
    # Define hydrostatic pressure basis.
    pressure_basis = iron.Basis()
    pressure_basis.CreateStart(pressure_basis_user_number)
    pressure_basis.TypeSet(iron.BasisTypes.LAGRANGE_HERMITE_TP)
    pressure_basis.NumberOfXiSet(number_of_xi)
    pressure_basis.InterpolationXiSet(xi_interpolation)
    pressure_basis.QuadratureNumberOfGaussXiSet(number_of_guass_xi)
    pressure_basis.CreateFinish()


# ### 3. Creating a region
# 
# Next we create a region that our fields will be defined on, and tell it to use the 3D coordinate system we created previously.

# In[5]:


# Create a region and assign the coordinate system to the region.
region_user_number = 1
region = iron.Region()
region.CreateStart(region_user_number,iron.WorldRegion)
region.LabelSet("Region")
region.CoordinateSystemSet(coordinate_system)
region.CreateFinish()


# ### 4. Setting up a simple cuboid mesh
# 
# In this example we will use the `iron.GeneratedMesh()` class of OpenCMISS to automatically create a 3D geometric mesh on which to solve the mechanics problem. We will create a regular mesh of size `width` (defined along x), `height` (defined along y) and `depth` (defined along z) and divide the mesh into `number_global_x_elements` in the X direction, `number_global_y_elements` in the Y direction and  `number_global_z_elements` in the Z direction. We will then tell it to use the basis we created previously.

# In[6]:


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
        [number_global_x_elements,
         number_global_y_elements,
         number_global_z_elements])
# Finish the creation of a generated mesh in the region.
mesh_user_number = 1
mesh = iron.Mesh()
generated_mesh.CreateFinish(mesh_user_number,mesh)


# ### 5. Decomposing the mesh
# Once the mesh has been created we can decompose it into a number of domains in order to allow for parallelism. We choose the options to let OpenCMISS-Iron to calculate the best way to break up the mesh. We also set the number of domains to be equal to the number of computational nodes this example is running on. In this example, we will only be using a single domain. Look for our parallelisation example for a demonstratoin of how to execute simulations using parallel processing techniques.

# In[7]:


# Get the number of computational nodes.
number_of_computational_nodes = iron.ComputationalNumberOfNodesGet()

# Create a decomposition for the mesh.
decomposition_user_number = 1
decomposition = iron.Decomposition()
decomposition.CreateStart(decomposition_user_number,mesh)
decomposition.TypeSet(iron.DecompositionTypes.CALCULATED)
decomposition.NumberOfDomainsSet(number_of_computational_nodes)
decomposition.CreateFinish()


# ### 6. Creating a geometric field
# Now that the mesh has been decomposed we are in a position to create fields. The first field we need to create is the geometric field. Once we have finished creating the field, we can change the field degrees of freedom (DOFs) to give us our geometry. Since the mesh was constructed using the OpenCMISS-Iron `GenerateMesh` class, we can use it's `GeometricParametersCalculate` method to automatically calculate and populate the geometric field parameters of the regular mesh.

# In[8]:


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
    geometric_field.FieldScalingTypeSet(iron.FieldScalingTypes.ARITHMETIC_MEAN)
geometric_field.CreateFinish()

# Update the geometric field parameters from generated mesh.
generated_mesh.GeometricParametersCalculate(geometric_field)


# ### Visualising the geometry
# We now visualise the geometry using pythreejs. Note that this visualisation currently only supports elements with linear Lagrange interpolation. This includes the node numbers for all elements.

# In[9]:


# ### 7. Creating fields
# 
# #### Dependent field
# 
# When solving the mechanics equations set, we require somewhere to store the deformed geometry (our solution). In OpenCMISS-Iron, we store the solutions to equations sets in a dependent field that contains our dependent variables. Note that the dependent field has been pre-defined in the OpenCMISS-Iron library to contain four components when solving `ProblemTypes.FINITE_ELASTICITY` with a `EquationsSetSubTypes.MOONEY_RIVLIN` subtype. The first three components store the deformed coordinates and the fourth stores the hydrostatic pressure.
# 
# Remember that `use_pressure_basis` can be set to true or false to switch between element based or nodally based interpolation for the incompressibility constraint.
# 
# One can either make the hydrostatic pressure constant within an element by calling the `dependent_field.ComponentInterpolationSet` method with the `iron.FieldInterpolationTypes.ELEMENT_BASED` option. Alternatively, we can make the hydrostatic pressure variable across different nodes by calling the `dependent_field.ComponentInterpolationSet` method with the `iron.FieldInterpolationTypes.NODE_BASED` option. In this tutorial, we have chosen the, hydrostatic pressure to be nodally interpolated if `use_pressure_basis` is true, or use element based interpolation if `use_pressure_basis` is false.

# In[10]:


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


# This dependent field needs to be initialised before the simulation is run. To this end, we copy the values of the coordinates from the geometric field into the dependent field in the below code snippet. The hydrostatic pressure field is set to 0.0.

# In[11]:


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


# #### Material field
# We now set up a new field called the material field, which will store the constitutive equation parameters of the Mooney Rivlin equation. This field can be set to have the same values throughout the mesh to represent a homogeneous material as shown below. If you want to describe heterogeneous materials, you can set the values of these parameters differently across the mesh (e.g. either using an nodally varying variation, or a element constant variation). This is not shown here but we'll give you a little example of this at the end of this tutorial. Below, the ```ComponentValuesInitialiseDP``` function sets all the nodal values to be the same: 1.0 for `c10` and 0.2 for `c01`.

# In[12]:


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
    iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, 1.0)
material_field.ComponentValuesInitialiseDP(
    iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 2, 0.2)


# #### Equation set field

# In[13]:


# Equation set field.
equations_set_field_user_number = 4
equations_set_field = iron.Field()
equations_set = iron.EquationsSet()


# ### 8. Defining the finite elasticity equations set
# We now specify that we want to solve a finite elasticity equation, and identify the specific material constitutive equation that we wish to use to describe the mechanical behaviour of the cube.
# 
# The key constants used to define the equation set are:
# - `ProblemClasses.ELASTICITY` defines that the equation set is of the elasticity class.
# - `ProblemTypes.FINITE_ELASTICITY` defines that the finite elasticity equations set will be used.
# - `EquationsSetSubTypes.MOONEY_RIVLIN` defines that the Mooney Rivlin constitutive equation that should be used from the range of constitutive equations that are implemented within the OpenCMISS-Iron library. You can find more information on these by browsing the OpenCMISS-Iron [Equations Set Subtypes Constants](http://opencmiss.org/documentation/apidoc/iron/latest/python/classiron_1_1_equations_set_subtypes.html). Future tutorials will demonstrate how you can dynamically specify constitutive relations using CellML.

# In[14]:


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


# Once the equations set is defined, we create the equations that use our fields to construct equations matrices and vectors.

# In[15]:


# Create equations.
equations = iron.Equations()
equations_set.EquationsCreateStart(equations)
equations.SparsityTypeSet(iron.EquationsSparsityTypes.SPARSE)
equations.OutputTypeSet(iron.EquationsOutputTypes.NONE)
equations_set.EquationsCreateFinish()


# ### 9. Defining the problem
# 
# Now that we have defined the equations we can now create our problem to be solved by OpenCMISS-Iron. We create a standard finite elasticity problem, which is a member of the elasticity field problem class. The problem control loop uses the default load increment loop and hence does not require a subtype.

# In[16]:


# Define the problem.
problem_user_number = 1

problem = iron.Problem()
problem_specification = (
  [iron.ProblemClasses.ELASTICITY,
   iron.ProblemTypes.FINITE_ELASTICITY,
   iron.ProblemSubtypes.NONE])
problem.CreateStart(problem_user_number, problem_specification)
problem.CreateFinish()


# ### 10. Defining control loops
# 
# The problem type defines a control loop structure that is used when solving the problem. The OpenCMISS-Iron control loop is a "supervisor" for the computational process. We may have multiple control loops with nested sub loops, and control loops can have different types, for example load incremented loops or time loops for dynamic problems. These control loops have been defined in the OpenCMISS-Iron library for the finite elasticity type of equations as a load increment loop. If we wanted to access the control loop and modify it we would use the `problem.ControlLoopGet` method before finishing the creation of the control loops. In the below code snippet we get the control loop to set the number of load increments to be used to solve the problem using the variable `number_of_load_increments`.

# In[17]:


# Create the problem control loop.
problem.ControlLoopCreateStart()
control_loop = iron.ControlLoop()
problem.ControlLoopGet([iron.ControlLoopIdentifiers.NODE], control_loop)
control_loop.MaximumIterationsSet(number_of_load_increments)
problem.ControlLoopCreateFinish()


# ### 11. Defining solvers
# 
# After defining the problem structure we can create the solvers that will be run to actually solve our problem. As the finite elasticity equations are nonlinear, we require a nonlinear solver. Nonlinear solvers typically involve a linearisation step, and therefore a linear solver is also required. In OpenCMISS-Iron, we can start the creation of solvers by calling the `problem.SolversCreateStart()` method, whose properties can be specified, before they are finalised using a call to the `problem.SolversCreateFinish()` method. Once finalised, only solver parameter (.e.g tolerances) can be changed, however, fundamental properties (e.g. which solver library to use) cannot be changed. If an additional solver is required, the existing solver can be destroyed and recreated, or another solver can be constructed.

# In[18]:


nonlinear_solver = iron.Solver()
linear_solver = iron.Solver()
problem.SolversCreateStart()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE], 1, nonlinear_solver)
nonlinear_solver.OutputTypeSet(iron.SolverOutputTypes.NONE)
nonlinear_solver.NewtonJacobianCalculationTypeSet(
  iron.JacobianCalculationTypes.EQUATIONS)
nonlinear_solver.NewtonLinearSolverGet(linear_solver)
linear_solver.LinearTypeSet(iron.LinearSolverTypes.DIRECT)
problem.SolversCreateFinish()


# ### 12. Defining solver equations
# After defining our solver we can create the equations for the solver to solve. This is achived by adding our equations sets to an OpenCMISS-Iron `solver_equations` object. In this example, we have just one equations set to add, however, for coupled problems, we may have multiple equations sets in the solver equations.

# In[19]:


solver = iron.Solver()
solver_equations = iron.SolverEquations()
problem.SolverEquationsCreateStart()
problem.SolverGet([iron.ControlLoopIdentifiers.NODE], 1, solver)
solver.SolverEquationsGet(solver_equations)
solver_equations.SparsityTypeSet(iron.SolverEquationsSparsityTypes.SPARSE)
_ = solver_equations.EquationsSetAdd(equations_set)
problem.SolverEquationsCreateFinish()


# ### 13. Defining the boundary conditions
# 
# The final step in configuring the problem is to define the boundary conditions for the simulations. Here, as stated at the top of the tutorial, we have set up five different boundary conditions settings to represent four independent loading conditions on the cube geometry.
# 
# - Model 1 (Uniaxial extension of a unit cube)
# - Model 2 (Equibiaxial extension of a unit cube)
# - Model 3 (Simple shear of a unit cube)
# - Model 4 (Shear of a unit cube)
# - Model 5 (Extension and shear of a unit cube)
# 
# The variable `model` set at the top of the tutorial program can be used to switch between these deformations. Each line of code below sets a dirichlet boundary conditions that prescribe a nodal coordinate value. The constant `iron.BoundaryConditionsTypes.FIXED` indicates that the value needs to be fixed to a certain value specified in the final argument of the `boundary_conditions.AddNode` method.

# In[20]:


# Prescribe boundary conditions (absolute nodal parameters).
boundary_conditions = iron.BoundaryConditions()
solver_equations.BoundaryConditionsCreateStart(boundary_conditions)
if model == 1:
    boundary_conditions.AddNode(dependent_field, iron.FieldVariableTypes.U,1,1,1,X,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundary_conditions.AddNode(dependent_field, iron.FieldVariableTypes.U,1,1,3,X,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundary_conditions.AddNode(dependent_field, iron.FieldVariableTypes.U,1,1,5,X,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundary_conditions.AddNode(dependent_field, iron.FieldVariableTypes.U,1,1,7,X,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundary_conditions.AddNode(dependent_field, iron.FieldVariableTypes.U,1,1,2,X,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundary_conditions.AddNode(dependent_field, iron.FieldVariableTypes.U,1,1,4,X,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundary_conditions.AddNode(dependent_field, iron.FieldVariableTypes.U,1,1,6,X,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundary_conditions.AddNode(dependent_field, iron.FieldVariableTypes.U,1,1,8,X,iron.BoundaryConditionsTypes.FIXED,0.0)

    boundary_conditions.AddNode(dependent_field, iron.FieldVariableTypes.U,1,1,1,Y,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundary_conditions.AddNode(dependent_field, iron.FieldVariableTypes.U,1,1,2,Y,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundary_conditions.AddNode(dependent_field, iron.FieldVariableTypes.U,1,1,5,Y,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundary_conditions.AddNode(dependent_field, iron.FieldVariableTypes.U,1,1,6,Y,iron.BoundaryConditionsTypes.FIXED,0.0)

    boundary_conditions.AddNode(dependent_field, iron.FieldVariableTypes.U,1,1,1,Z,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundary_conditions.AddNode(dependent_field, iron.FieldVariableTypes.U,1,1,2,Z,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundary_conditions.AddNode(dependent_field, iron.FieldVariableTypes.U,1,1,3,Z,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundary_conditions.AddNode(dependent_field, iron.FieldVariableTypes.U,1,1,4,Z,iron.BoundaryConditionsTypes.FIXED,0.0)

elif model == 2:
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,1,X,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,3,X,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,5,X,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,7,X,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,2,X,iron.BoundaryConditionsTypes.FIXED,0.25)
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,4,X,iron.BoundaryConditionsTypes.FIXED,0.25)
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,6,X,iron.BoundaryConditionsTypes.FIXED,0.25)
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,8,X,iron.BoundaryConditionsTypes.FIXED,0.25)

    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,1,Y,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,2,Y,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,5,Y,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,6,Y,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,3,Y,iron.BoundaryConditionsTypes.FIXED,0.25)
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,4,Y,iron.BoundaryConditionsTypes.FIXED,0.25)
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,7,Y,iron.BoundaryConditionsTypes.FIXED,0.25)
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,8,Y,iron.BoundaryConditionsTypes.FIXED,0.25)

    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,1,Z,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,2,Z,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,3,Z,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,4,Z,iron.BoundaryConditionsTypes.FIXED,0.0)

elif model == 3:
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,1,X,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,3,X,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,5,X,iron.BoundaryConditionsTypes.FIXED,0.5)
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,7,X,iron.BoundaryConditionsTypes.FIXED,0.5)
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,2,X,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,4,X,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,6,X,iron.BoundaryConditionsTypes.FIXED,0.5)
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,8,X,iron.BoundaryConditionsTypes.FIXED,0.5)

    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,1,Y,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,2,Y,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,5,Y,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,6,Y,iron.BoundaryConditionsTypes.FIXED,0.0)

    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,1,Z,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,2,Z,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,3,Z,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,4,Z,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,5,Z,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,6,Z,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,7,Z,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,8,Z,iron.BoundaryConditionsTypes.FIXED,0.0)

elif model == 4:
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,1,X,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,3,X,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,6,X,iron.BoundaryConditionsTypes.FIXED,0.5)
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,8,X,iron.BoundaryConditionsTypes.FIXED,0.5)

    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,1,Y,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,2,Y,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,5,Y,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,6,Y,iron.BoundaryConditionsTypes.FIXED,0.0)

    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,1,Z,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,3,Z,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,6,Z,iron.BoundaryConditionsTypes.FIXED,0.5)
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,8,Z,iron.BoundaryConditionsTypes.FIXED,0.5)

elif model == 5:
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,1,X,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,3,X,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,5,X,iron.BoundaryConditionsTypes.FIXED,0.5)
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,7,X,iron.BoundaryConditionsTypes.FIXED,0.5)
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,2,X,iron.BoundaryConditionsTypes.FIXED,0.25)
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,4,X,iron.BoundaryConditionsTypes.FIXED,0.25)
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,6,X,iron.BoundaryConditionsTypes.FIXED,0.75)
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,8,X,iron.BoundaryConditionsTypes.FIXED,0.75)

    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,1,Y,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,2,Y,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,5,Y,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,6,Y,iron.BoundaryConditionsTypes.FIXED,0.0)

    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,1,Z,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,2,Z,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,3,Z,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,4,Z,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,5,Z,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,6,Z,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,7,Z,iron.BoundaryConditionsTypes.FIXED,0.0)
    boundary_conditions.AddNode(dependent_field,iron.FieldVariableTypes.U,1,1,8,Z,iron.BoundaryConditionsTypes.FIXED,0.0)


# We then construct the solver matrices and vectors by making a call to the `solver_equations.BoundaryConditionsCreateFinish()` method.

# In[21]:


solver_equations.BoundaryConditionsCreateFinish()


# ### 14. Solving the problem
# After our problem solver equations have been fully defined, we are now ready to solve our problem. When we call the Solve method of the problem it will loop over the control loops and control loop solvers to solve our problem.

# In[22]:


# Solve the problem.
problem.Solve()


# ## Visualising results
# We can now visualise the resulting deformation as animation using pythreejs.

# In[23]:
# ## Exporting results
# Before we export the results in Cmgui format, we will first create a new `deformed_field` and `pressure_field` to separately hold the solution for the deformed geometry and hydrostatic pressure for visualisation in Cmgui (this enables simplified access to these fields in Cmgui visualisation scripts).

# In[24]:


deformed_field_user_number = 5
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

# Export results to exnode and exelem format.
fields = iron.Fields()
fields.CreateRegion(region)
fields.NodesExport("cube", "FORTRAN")
fields.ElementsExport("cube", "FORTRAN")
fields.Finalise()


# ## Evaluating mechanics tensor fields

# In[25]:


results = {}
elementNumber = 1
xiPosition = [0.5, 0.5, 0.5]
F = equations_set.TensorInterpolateXi(
    iron.EquationsSetDerivedTensorTypes.DEFORMATION_GRADIENT,
    elementNumber, xiPosition,(3,3))
results['Deformation Gradient Tensor'] = F
print("Deformation Gradient Tensor")
print(F)


# In[26]:


C = equations_set.TensorInterpolateXi(
    iron.EquationsSetDerivedTensorTypes.R_CAUCHY_GREEN_DEFORMATION,
    elementNumber, xiPosition,(3,3))
results['Right Cauchy-Green Deformation Tensor'] = C
print("Right Cauchy-Green Deformation Tensor")
print(C)


# In[27]:


E = equations_set.TensorInterpolateXi(
    iron.EquationsSetDerivedTensorTypes.GREEN_LAGRANGE_STRAIN,
    elementNumber, xiPosition,(3,3))
results['Green-Lagrange Strain Tensor'] = E
print("Green-Lagrange Strain Tensor")
print(E)


# In[28]:


I1=numpy.trace(C)
I2=0.5*(numpy.trace(C)**2.-numpy.tensordot(C,C))
I3=numpy.linalg.det(C)
results['Invariants'] = [I1, I2, I3]
print("Invariants")
print("I1={0}, I2={1}, I3={2}".format(I1,I2,I3))


# In[29]:


TC = equations_set.TensorInterpolateXi(
    iron.EquationsSetDerivedTensorTypes.CAUCHY_STRESS,
    elementNumber, xiPosition,(3,3))
results['Cauchy Stress Tensor'] = TC
print("Cauchy Stress Tensor")
print(TC)


# In[30]:


# Calculate the Second Piola-Kirchhoff Stress Tensor from
# TG=J*F^(-1)*TC*F^(-T), where T denotes the
# matrix transpose, and assumes J=1.
TG = numpy.dot(numpy.linalg.inv(F),numpy.dot(
        TC,numpy.linalg.inv(numpy.matrix.transpose(F))))
results['Second Piola-Kirchhoff Stress Tensor'] = TG
print("Second Piola-Kirchhoff Stress Tensor")
print(TG)


# In[31]:


p = dependent_field.ParameterSetGetElement(
    iron.FieldVariableTypes.U,
    iron.FieldParameterSetTypes.VALUES,elementNumber,4)
results['Hydrostatic pressure'] = p
print("Hydrostatic pressure")
print(p)


# ## Finalising session

# In[32]:


problem.Destroy()
coordinate_system.Destroy()
region.Destroy()
basis.Destroy()
iron.Finalise()


# ## Additional exercises

# ### Changing boundary conditions
# 
# Change the value of **model** and re-run the notebook to see the effect of altered boundary conditions.
# 
