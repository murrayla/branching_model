from ipywidgets import embed
import pythreejs as pjs
from IPython.display import display
import numpy as np
import matplotlib.cm as cm
import matplotlib
import mesh_tools

from pythreejs.materials.SpriteMaterial_autogen import SpriteMaterial
from pythreejs.textures.TextTexture_autogen import TextTexture
from pythreejs.objects.Sprite_autogen import Sprite

def visualise(
        mesh, geometric_field, number_of_dimensions, xi_interpolation,
        dependent_field=None, variable=None, mechanics_animation=False,
        colour_map_dependent_component_number=None,
        cmap='gist_rainbow', resolution=1, node_labels=False):

    if number_of_dimensions != 3:
        print(
            'Warning: Only visualisation of 3D meshes is currently supported.')
        return

    if xi_interpolation != [1, 1, 1]:
        print(
            'Warning: Only visualisation of 3D elements with linear Lagrange \
            interpolation along all coordinate directions is currently \
            supported.')
        return

    view_width = 600
    view_height = 600

    debug = False
    if debug:
        vertices = [
            [0, 0, 0],
            [0, 0, 1],
            [0, 1, 0],
            [0, 1, 1],
            [1, 0, 0],
            [1, 0, 1],
            [1, 1, 0],
            [1, 1, 1]
        ]

        faces = [
            [0, 1, 3],
            [0, 3, 2],
            [0, 2, 4],
            [2, 6, 4],
            [0, 4, 1],
            [1, 4, 5],
            [2, 3, 6],
            [3, 7, 6],
            [1, 5, 3],
            [3, 5, 7],
            [4, 6, 5],
            [5, 6, 7]
        ]

        vertexcolors = ['#000000', '#0000ff', '#00ff00', '#ff0000',
                        '#00ffff', '#ff00ff', '#ffff00', '#ffffff']
    else:
        # Get mesh topology information.
        num_nodes = mesh_tools.num_nodes_get(mesh, mesh_component=1)
        node_nums = list(range(1, num_nodes + 1))
        num_elements, element_nums = mesh_tools.num_element_get(
            mesh, mesh_component=1)

        # Convert geometric field to a morphic mesh and export to json
        mesh = mesh_tools.OpenCMISS_to_morphic(
            mesh, geometric_field, element_nums, node_nums, dimension=3,
            interpolation='linear')
        vertices, faces, _, xi_element_nums, xis = get_faces(
            mesh, res=resolution, exterior_only=True, include_xi=True)

        vertices = vertices.tolist()
        faces = faces.tolist()

    centroid = np.mean(vertices, axis=0)
    max_positions = np.max(vertices, axis=0)
    min_positions = np.min(vertices, axis=0)
    range_positions = max_positions - min_positions

    if (dependent_field is not None) and (colour_map_dependent_component_number is not None):

        solution = np.zeros(xis.shape[0])
        for idx, (xi, xi_element_num) in enumerate(zip(xis, xi_element_nums)):
            solution[idx] = mesh_tools.interpolate_opencmiss_field_xi(
                dependent_field, xi, element_ids=[xi_element_num], dimension=3,
                deriv=1)[colour_map_dependent_component_number-1]

        minima = min(solution)
        maxima = max(solution)

        import matplotlib
        norm = matplotlib.colors.Normalize(vmin=minima, vmax=maxima, clip=True)
        mapper = cm.ScalarMappable(
            norm=norm, cmap=cm.get_cmap(name=cmap))

        vertex_colors = np.zeros((len(vertices), 3), dtype='float32')
        for idx, v in enumerate(solution):
            vertex_colors[idx, :] = mapper.to_rgba(v, alpha=None)[:3]
        # else:
        #     raise ValueError('Visualisation not supported.')
    else:
        vertex_colors = np.tile(
            np.array([0.5, 0.5, 0.5], dtype='float32'),
            (len(vertices), 1))

    geometry = pjs.BufferGeometry(attributes=dict(
        position=pjs.BufferAttribute(vertices, normalized=False),
        index=pjs.BufferAttribute(
            np.array(faces).astype(dtype='uint16').ravel(),
            normalized=False),
        color=pjs.BufferAttribute(vertex_colors),
    ))

    if mechanics_animation:
        deformed_vertices = np.zeros((xis.shape[0], 3), dtype='float32')
        for idx, (xi, xi_element_num) in enumerate(
                zip(xis, xi_element_nums)):
            deformed_vertices[idx, :] = \
            mesh_tools.interpolate_opencmiss_field_xi(
                dependent_field, xi, element_ids=[xi_element_num],
                dimension=3,
                deriv=1)[0][:3]
        geometry.morphAttributes = {'position': [
            pjs.BufferAttribute(deformed_vertices),
        ]}

        geometry.exec_three_obj_method('computeFaceNormals')
        geometry.exec_three_obj_method('computeVertexNormals')

        surf1 = pjs.Mesh(
            geometry,
            pjs.MeshPhongMaterial(
                color='#ff3333',
                shininess=150,
                morphTargets=True,
                side='FrontSide'), name='A')
        surf2 = pjs.Mesh(
            geometry,
            pjs.MeshPhongMaterial(
                color='#ff3333',
                shininess=150,
                morphTargets=True,
                side='BackSide'), name='B')
        surf = pjs.Group(children=[surf1, surf2])

        # camera = pjs.PerspectiveCamera(
        #     fov=20, position=[range_positions[0] * 10,
        #                       range_positions[1] * 10,
        #                       range_positions[2] * 10],
        #     width=view_width,
        #     height=view_height, near=1,
        #     far=max(range_positions) * 10)

        camera = pjs.PerspectiveCamera(
            position=[range_positions[0] * 3,
                      range_positions[1] * 3,
                      range_positions[2] * 3],
            aspect=view_width / view_height)
        camera.up = [0, 0, 1]
        camera.lookAt(centroid.tolist())

        scene3 = pjs.Scene(children=[surf1, surf2, camera,
                                     pjs.DirectionalLight(position=[3, 5, 1],
                                                          intensity=0.6),
                                     pjs.AmbientLight(intensity=0.5)])
        axes = pjs.AxesHelper(size=range_positions[0] * 2)
        scene3.add(axes)

        A_track = pjs.NumberKeyframeTrack(
            name='scene/A.morphTargetInfluences[0]', times=[0, 3],
            values=[0, 1])
        B_track = pjs.NumberKeyframeTrack(
            name='scene/B.morphTargetInfluences[0]', times=[0, 3],
            values=[0, 1])
        pill_clip = pjs.AnimationClip(tracks=[A_track, B_track])
        pill_action = pjs.AnimationAction(
            pjs.AnimationMixer(scene3), pill_clip, scene3)

        renderer3 = pjs.Renderer(
            camera=camera, scene=scene3,
            controls=[pjs.OrbitControls(controlling=camera)],
            width=view_width, height=view_height)


        display(renderer3, pill_action)

    else:
        geometry.exec_three_obj_method('computeFaceNormals')
        geometry.exec_three_obj_method('computeVertexNormals')

        surf1 = pjs.Mesh(geometry=geometry,
                     material=pjs.MeshLambertMaterial(
                         vertexColors='VertexColors',
                         side='FrontSide'))  # Center the cube.
        surf2 = pjs.Mesh(geometry=geometry,
                     material=pjs.MeshLambertMaterial(
                         vertexColors='VertexColors',
                         side='BackSide'))  # Center the cube.
        surf = pjs.Group(children=[surf1, surf2])

        camera = pjs.PerspectiveCamera(position=[range_positions[0] * 3,
                              range_positions[1] * 3,
                              range_positions[2] * 3],
                                        aspect=view_width / view_height)

        camera.up = [0, 0, 1]
        camera.lookAt(centroid.tolist())

        # if perspective:
        #     camera.mode = 'perspective'
        # else:
        #     camera.mode = 'orthographic'

        lights = [
            pjs.DirectionalLight(
                position=[range_positions[0] * 16,
                          range_positions[1] * 12,
                          range_positions[2] * 17], intensity=0.5),
            pjs.AmbientLight(intensity=0.8),
        ]
        orbit = pjs.OrbitControls(
            controlling=camera, screenSpacePanning=True, target=centroid.tolist())

        scene = pjs.Scene()
        axes = pjs.AxesHelper(size=max(range_positions) * 2)
        scene.add(axes)
        scene.add(surf1)
        scene.add(surf2)
        scene.add(lights)

        if node_labels:
            # Add text labels for each mesh node.
            v, ids = mesh.get_node_ids(group='_default')
            for idx, v in enumerate(v):
                text = make_text(str(ids[idx]), position=(v[0], v[1], v[2]))
                scene.add(text)

        # Add text for axes labels.
        x_axis_label = make_text(
            'x', position=(max(range_positions) * 2, 0, 0))
        y_axis_label = make_text(
            'y', position=(0, max(range_positions) * 2, 0))
        z_axis_label = make_text(
            'z', position=(0, 0, max(range_positions) * 2))
        scene.add(x_axis_label)
        scene.add(y_axis_label)
        scene.add(z_axis_label)

        renderer = pjs.Renderer(
            scene=scene, camera=camera, controls=[orbit],
            width=view_width, height=view_height)
        camera.zoom = 1
        display(renderer)

    return vertices, faces

def make_text(text, position=(0, 0, 0), colour='white', size=0.05):
    """
    Return a text object at the specified location with a given height
    """
    sm = SpriteMaterial(map=TextTexture(string=text, color=colour, size=100, squareTexture=False))
    return Sprite(material=sm, position=position, scaleToTexture=True, scale=[size, size*1.5, size])

def export_html(renderer):
    embed.embed_minimal_html('export.html', views=renderer, title='Renderer')

def get_faces(
        mesh, res=8, exterior_only=True, include_xi=False, elements=None):
    mesh.generate()

    if elements == None:
        Faces = mesh.faces
    else:
        Faces = []
        for face in mesh.faces:
            for element_face in face.element_faces:
                if element_face[0] in elements:
                    Faces.append(face)
                    break

    if exterior_only:
        Faces = [face for face in Faces if len(face.element_faces) == 1]

    basis = mesh.elements[mesh.get_element_ids()[0]].basis

    if Faces[0].shape == 'quad':
        face_metadata = {
            'face_xi_idxs': {
                0: [0, 1],
                1: [0, 1],
                2: [0, 2],
                3: [0, 2],
                4: [1, 2],
                5: [1, 2]
            },
            'face_normal_xi_idx': {
                0: 2,
                1: 2,
                2: 1,
                3: 1,
                4: 0,
                5: 0
            },
            'face_normal_xi_value': {
                0: 0,
                1: 1,
                2: 0,
                3: 1,
                4: 0,
                5: 1
            }
        }

    XiT, TT = xi_grid(shape='tri', res=res)
    XiQ, TQ = xi_grid(shape='quad', res=res)
    NPT, NTT = XiT.shape[0], TT.shape[0]
    NPQ, NTQ = XiQ.shape[0], TQ.shape[0]

    XiQ0 = np.zeros(NPQ)
    XiQ1 = np.ones(NPQ)

    NP, NT = 0, 0
    for face in Faces:
        if face.shape == 'tri':
            NP += NPT
            NT += NTT
        elif face.shape == 'quad':
            NP += NPQ
            NT += NTQ

    X = np.zeros((NP, 3))  #######TODO#####face.nodes[0].num_fields))
    T = np.zeros((NT, 3), dtype='uint32')
    if include_xi:
        Xi = np.zeros((NP, 2))
        Xe = np.zeros(NP, dtype='uint32')
        Xi3D = np.zeros((NP, 3))
    npp, nt = 0, 0
    for face in Faces:
        if face.shape == 'tri':
            X[npp:npp + NPT, :] = mesh._core.evaluate(face.cid, XiT)
            if include_xi:
                Xi[npp:npp + NPT, :] = XiT
            T[nt:nt + NTT, :] = TT + npp
            npp += NPT
            nt += NTT
        elif face.shape == 'quad':
            elem = mesh.elements[face.element_faces[0][0]]
            face_index = face.element_faces[0][1]
            if face_index == 0:
                X[npp:npp + NPQ, :] = mesh._core.evaluate(
                    elem.cid,
                    np.array([XiQ[:, 0],
                              XiQ[:, 1],
                              XiQ0]).T)
            elif face_index == 1:
                X[npp:npp + NPQ, :] = mesh._core.evaluate(
                    elem.cid,
                    np.array([XiQ[:, 0],
                              XiQ[:, 1],
                              XiQ1]).T)
            elif face_index == 2:
                X[npp:npp + NPQ, :] = mesh._core.evaluate(
                    elem.cid,
                    np.array(
                        [XiQ[:, 0], XiQ0,
                         XiQ[:, 1]]).T)
            elif face_index == 3:
                X[npp:npp + NPQ, :] = mesh._core.evaluate(
                    elem.cid,
                    np.array(
                        [XiQ[:, 0], XiQ1,
                         XiQ[:, 1]]).T)
            elif face_index == 4:
                X[npp:npp + NPQ, :] = mesh._core.evaluate(
                    elem.cid,
                    np.array(
                        [XiQ0, XiQ[:, 0],
                         XiQ[:, 1]]).T)
            elif face_index == 5:
                X[npp:npp + NPQ, :] = mesh._core.evaluate(
                    elem.cid,
                    np.array(
                        [XiQ1, XiQ[:, 0],
                         XiQ[:, 1]]).T)

            T[nt:nt + NTQ, :] = TQ + npp
            if include_xi:
                Xi[npp:npp + NPQ, :] = XiQ
                Xe[npp:npp + NPQ] = elem.id
                Xi3D[npp:npp + NPQ, face_metadata['face_xi_idxs'][face_index]] = XiQ
                Xi3D[npp:npp + NPQ, face_metadata['face_normal_xi_idx'][face_index]] = face_metadata['face_normal_xi_value'][face_index]
            npp += NPQ
            nt += NTQ
    if include_xi:
        return X, T, Xi, Xe, Xi3D
    return X, T


def xi_grid(shape='quad', res=[8, 8], units='div', method='fit'):
    if units == 'div':
        if isinstance(res, int):
            divs = [res, res]
        else:
            divs = res
    elif units == 'xi':
        raise TypeError('Unimplemented units')

    nx = divs[0] + 1
    dx = 0.5 / nx
    if method == 'fit':
        xi = np.linspace(0, 1, divs[0] + 1)
    elif method == 'center':
        xi = np.linspace(dx, 1 - dx, divs[0] + 1)
    else:
        xi = np.linspace(0, 1, divs[0] + 1)

    if shape == 'quad':
        NPQ = int(nx * nx)
        NTQ = int(2 * (divs[0] * divs[0]))

        xi1, xi2 = np.meshgrid(xi, xi)
        xi1 = xi1.reshape([xi1.size])
        xi2 = xi2.reshape([xi2.size])
        XiQ = np.array([xi1, xi2]).T
        TQ = np.zeros((NTQ, 3), dtype='uint32')
        npp = 0
        for row in range(divs[0]):
            for col in range(divs[0]):
                NPPR = row * nx
                TQ[npp, :] = [NPPR + col, NPPR + col + 1, NPPR + col + nx]
                npp += 1
                TQ[npp, :] = [NPPR + col + 1, NPPR + col + nx + 1,
                             NPPR + col + nx]
                npp += 1

        return XiQ, TQ

    elif shape == 'tri':
        NPT = int(0.5 * nx * (nx - 1) + nx)
        NTT = int(divs[0] * divs[0])

        XiT = np.zeros([NPT, 2])
        TT = np.zeros((NTT, 3), dtype='uint32')
        NodesPerLine = range(divs[0], 0, -1)
        npp = 0
        for row in range(nx):
            for col in range(nx - row):
                XiT[npp, 0] = xi[col]
                XiT[npp, 1] = xi[row]
                npp += 1

        npp = 0
        ns = 0
        for row in range(divs[0]):
            for col in range(divs[0] - row):
                TT[npp, :] = [ns, ns + 1, ns + nx - row]
                npp += 1
                if col != divs[0] - row - 1:
                    TT[npp, :] = [ns + 1, ns + nx - row + 1, ns + nx - row]
                    npp += 1
                ns += 1
            ns += 1

        return XiT, TT


if __name__ == '__main__':
    # Intialise OpenCMISS-Iron.
    from opencmiss.iron import iron

    # Create coordinate system.
    coordinate_system_user_number = 1
    coordinate_system = iron.CoordinateSystem()
    coordinate_system.CreateStart(coordinate_system_user_number)
    coordinate_system.DimensionSet(3)
    coordinate_system.CreateFinish()

    # Create region.
    region_user_number = 1
    region = iron.Region()
    region.CreateStart(region_user_number, iron.WorldRegion)
    region.CoordinateSystemSet(coordinate_system)
    region.CreateFinish()

    # Create basis functions.
    basis_user_number = 1
    basis = iron.Basis()
    basis.CreateStart(basis_user_number)
    basis.CreateFinish()

    # Define mesh parameters.
    number_global_x_elements = 1
    number_global_y_elements = 1
    number_global_z_elements = 1
    height = 1.0
    width = 1.0
    length = 1.0

    # Create mesh.
    generated_mesh_user_number = 1
    generated_mesh = iron.GeneratedMesh()
    generated_mesh.CreateStart(generated_mesh_user_number, region)
    generated_mesh.TypeSet(iron.GeneratedMeshTypes.REGULAR)
    generated_mesh.BasisSet([basis])
    generated_mesh.ExtentSet([width, height, length])
    generated_mesh.NumberOfElementsSet(
        [number_global_x_elements,
         number_global_y_elements,
         number_global_z_elements])
    mesh = iron.Mesh()
    mesh_user_number = 1
    generated_mesh.CreateFinish(mesh_user_number, mesh)

    # Perform mesh decomposition.
    decomposition_user_number = 1
    decomposition = iron.Decomposition()
    decomposition.CreateStart(decomposition_user_number, mesh)
    decomposition.CreateFinish()

    # Create geometric field.
    geometric_field_user_number = 1
    geometric_field = iron.Field()
    geometric_field.CreateStart(geometric_field_user_number, region)
    geometric_field.MeshDecompositionSet(decomposition)
    geometric_field.CreateFinish()

    # Set geometric field values from the generated mesh.
    generated_mesh.GeometricParametersCalculate(geometric_field)

    import sys

    sys.path.insert(1, '../../tools/')
    import threejs_visualiser

    threejs_visualiser.visualise(
        mesh, geometric_field, variable=iron.FieldVariableTypes.U)

    equations_set_user_number = 1
    equations_set_field_user_number = 2
    equations_set_field = iron.Field()
    equations_set = iron.EquationsSet()
    equations_set_specification = [
        iron.EquationsSetClasses.CLASSICAL_FIELD,
        iron.EquationsSetTypes.LAPLACE_EQUATION,
        iron.EquationsSetSubtypes.STANDARD_LAPLACE]
    equations_set.CreateStart(
        equations_set_user_number, region, geometric_field,
        equations_set_specification, equations_set_field_user_number,
        equations_set_field)
    equations_set.CreateFinish()
    # DOC-END equation set

    # DOC-START dependent field
    # Create dependent field.
    dependent_field_user_number = 3
    dependent_field = iron.Field()
    equations_set.DependentCreateStart(
        dependent_field_user_number, dependent_field)
    equations_set.DependentCreateFinish()

    # DOC-START initialise dependent field
    # Initialise dependent field.
    dependent_field.ComponentValuesInitialiseDP(
        iron.FieldVariableTypes.U, iron.FieldParameterSetTypes.VALUES, 1, 0.5)
    # DOC-END initialise dependent field

    # DOC-START equations
    # Create equations.
    equations = iron.Equations()
    equations_set.EquationsCreateStart(equations)
    equations_set.EquationsCreateFinish()
    # DOC-END equations

    # DOC-START problem
    # Create problem.
    problem_user_number = 1
    problem = iron.Problem()
    problem_specification = [
        iron.ProblemClasses.CLASSICAL_FIELD,
        iron.ProblemTypes.LAPLACE_EQUATION,
        iron.ProblemSubtypes.STANDARD_LAPLACE]
    problem.CreateStart(problem_user_number, problem_specification)
    problem.CreateFinish()
    # DOC-END problem

    # DOC-START control loops
    # Create control loops.
    problem.ControlLoopCreateStart()
    problem.ControlLoopCreateFinish()
    # DOC-END control loops

    # DOC-START problem solver
    # Create problem solver.
    solver = iron.Solver()
    problem.SolversCreateStart()
    problem.SolverGet([iron.ControlLoopIdentifiers.NODE], 1, solver)
    solver.OutputTypeSet(iron.SolverOutputTypes.SOLVER)
    problem.SolversCreateFinish()
    # DOC-END problem solver

    # DOC-START solver equations
    # Create solver equations and add equations set to solver equations.
    solver = iron.Solver()
    solver_equations = iron.SolverEquations()
    problem.SolverEquationsCreateStart()
    problem.SolverGet([iron.ControlLoopIdentifiers.NODE], 1, solver)
    solver.SolverEquationsGet(solver_equations)
    solver_equations.EquationsSetAdd(equations_set)
    problem.SolverEquationsCreateFinish()
    # DOC-END solver equations

    # DOC-START boundary condition nodes
    # Identify first and last node number.
    firstNodeNumber = 1
    nodes = iron.Nodes()
    region.NodesGet(nodes)
    lastNodeNumber = nodes.NumberOfNodesGet()
    # DOC-END boundary condition nodes

    # DOC-START boundary conditions
    # Create boundary conditions and set first and last nodes to 0.0 and 1.0
    boundary_conditions = iron.BoundaryConditions()
    solver_equations.BoundaryConditionsCreateStart(boundary_conditions)
    boundary_conditions.SetNode(
        dependent_field, iron.FieldVariableTypes.U, 1, 1, firstNodeNumber,
        1, iron.BoundaryConditionsTypes.FIXED, 0.0)
    boundary_conditions.SetNode(
        dependent_field, iron.FieldVariableTypes.U, 1, 1, lastNodeNumber,
        1, iron.BoundaryConditionsTypes.FIXED, 1.0)
    solver_equations.BoundaryConditionsCreateFinish()
    # DOC-END boundary conditions

    # DOC-START solve
    problem.Solve()
    # DOC-END solve

    threejs_visualiser.visualise(
        mesh, geometric_field, dependent_field=dependent_field,
        variable=iron.FieldVariableTypes.U, component=1, perspective=True,
        animate=True)


