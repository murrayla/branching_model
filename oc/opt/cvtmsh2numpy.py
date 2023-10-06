# Author: Liam Murray
# Date: 29/03/2023
# File Name: msh2numpy.py
# Description: Script to convert '.msh' output form gmsh.app to numpy arrays to be input into OpenCMISS-iron

# BEGIN
# Read in the file 
def msh2numpy(file_name):
    msh_file    = open(file_name, 'r')

    # Find breaks which are indicated by '$'
    # Also retain information requires for Entities, Nodes and Elements
    # Set up containers and checks
    break_container = dict()
    entities_list   = list()
    ent_check       = 0
    nodes_list      = list()
    nod_check       = 0
    elements_list   = list()
    ele_check       = 0
    # Iterate
    for i, line in enumerate(msh_file):
        # Store break information
        if line[0][0] == '$':
            break_container[line[1:-1]] = i
        # Store information
        if ent_check:
            entities_list.append(line[:-1])
        if nod_check:
            nodes_list.append(line[:-1])
        if ele_check:
            elements_list.append(line[:-1])
        # Checks on if we are in the right section
        if line[1:-1] == 'Entities':
            ent_check = 1
            continue
        elif line[1:-1] == 'EndEntities':
            ent_check = 0
            continue
        elif line[1:-1] == 'Nodes':
            nod_check = 1
            continue
        elif line[1:-1] == 'EndNodes':
            nod_check = 0
            continue
        elif line[1:-1] == 'Elements':
            ele_check = 1 
            continue
        elif line[1:-1] == 'EndElements':
            ele_check = 0
            continue

    # Remove the last value of each list to compensate for the title added
    entities_list.pop()
    nodes_list.pop()
    elements_list.pop()

    # Augement Data Structure

    #   # ENTITIES #s
    # Break down the entities list into components
    entity_param_values = entities_list[0].split(" ")
    entity_params = {"numPoints": int(entity_param_values[0]), "numCurves": int(entity_param_values[1]), 
                    "numSurfaces": int(entity_param_values[2]), "numVolumes": int(entity_param_values[3])}
    

    #   # NODES #

    # Create a file to write to: gmsh2iron.node
    node_file = open('gmsh2iron.node', 'w')
    # Break down the node list into components
    node_param_values = nodes_list[0].split(" ")
    node_params = {"numEntityBlocks": int(node_param_values[0]), "numNodes": int(node_param_values[1]), 
                    "minNodeTag": int(node_param_values[2]), "maxNodeTag": int(node_param_values[3])}
    
    node_file.write(nodes_list[0] + "\n")
    # Loop through nodes and blocks
    node_positions = dict()
    count = 1
    for i, block in enumerate(nodes_list):
        if count:
            count -= 1
            continue
        # print(block.split(" ")[3])
        for j in range(int(block.split(" ")[3])):
            node_positions[int(nodes_list[i+j+1])] = (float(nodes_list[int(block.split(" ")[3])+i+j+1].split(" ")[0]), 
                                                    float(nodes_list[int(block.split(" ")[3])+i+j+1].split(" ")[1]), 
                                                    float(nodes_list[int(block.split(" ")[3])+i+j+1].split(" ")[2]))
            node_file.write(nodes_list[i+j+1]+"\t")
            node_file.write(nodes_list[int(block.split(" ")[3])+i+j+1].split(" ")[0] + "\t" +
                            nodes_list[int(block.split(" ")[3])+i+j+1].split(" ")[1] + "\t" +
                            nodes_list[int(block.split(" ")[3])+i+j+1].split(" ")[2] + "\n")
            count +=2

    # print(node_positions)

    #   # ELEMENTS #

    # Create a file to store the element information: gmsh2iron.ele
    element_file = open('gmsh2iron.ele', 'w')
    # Break down the elements list into components
    element_param_values = elements_list[0].split(" ")
    element_params = {"numEntityBlocks": int(element_param_values[0]), "numElements": int(element_param_values[1]), 
                    "minElementTag": int(element_param_values[2]), "maxElementTag": int(element_param_values[3])}
    

    element_file.write(elements_list[0] + "\n")
    # Set Element Types
    types = {1: "line", 2: "triangle", 3: "quadrangle", 4: "tetrahedron", 
            5: "hexahedron", 6: "prism", 7: "pyramid", 8: "second-order-line", 
            9: "second-order-triangle", 10: "second-order-quadrangle", 
            11: "second-order-tetrahedron", 15: "point"}

    # Loop through nodes and blocks
    element_info    = dict()
    element_group   = dict()
    count           = 1
    curr_shape      = "avocado"
    for i, block in enumerate(elements_list):
        if count:
            count -= 1
            continue
        for j in range(int(block.split(" ")[3])):
            element_info[(int(block.split(" ")[0]), int(block.split(" ")[1]))] = types[int(block.split(" ")[2])]
            element_group[(int(block.split(" ")[0]), int(block.split(" ")[1]), int(elements_list[i+j+1].split()[0]))] = [float(x) for x in elements_list[i+j+1].split()][1::]
            count +=1
            # if types[int(block.split(" ")[2])] != curr_shape:
            #     element_file.write(types[int(block.split(" ")[2])] + "\n")
            #     curr_shape = types[int(block.split(" ")[2])]
            element_file.write(types[int(block.split(" ")[2])] + "\t")
            element_file.write(block.split(" ")[1] + "\t")
            for value in elements_list[i+j+1].split():
                element_file.write(value + "\t")
            element_file.write("\n")


