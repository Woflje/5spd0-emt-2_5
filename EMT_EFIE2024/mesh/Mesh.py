from mesh.Object import Object
from mesh.Triangle import Triangle
from mesh.Edge import Edge
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
import os
import re

def calculate_triangle_area(vertices):
    a = np.linalg.norm(vertices[0] - vertices[1])
    b = np.linalg.norm(vertices[1] - vertices[2])
    c = np.linalg.norm(vertices[2] - vertices[0])
    s = (a + b + c) / 2
    return np.sqrt(s * (s - a) * (s - b) * (s - c))

def restructure_mesh_rwg_data(r_vect, dunavant_values):
    N = len(r_vect)
    # Flatten r_vect to get all vertices and deduplicate to get unique vertices
    all_vertices = r_vect.reshape(-1, 3)
    unique_vertices, inverse_indices = np.unique(all_vertices, axis=0, return_inverse=True)
    
    # Map vertices in r_vect to unique vertices indices
    mapped_indices = inverse_indices.reshape(r_vect.shape[0], r_vect.shape[1])

    # Initialize structures
    rwgs = np.zeros((N, 6), dtype=int)  # Mapping indices
    areas = []  # To store areas of each triangle

    triangle_tracker = {}  # Tracks whether a triangle has been processed

    for i, rwg in enumerate(mapped_indices):
        rwgs[i, :4] = rwg[:4]  # Fill in rwgs for vertices
        
        for j in range(2, 4):  # Loop over the unique vertices for each triangle
            triangle_key = tuple(sorted([rwg[0], rwg[1], rwg[j]]))  # Create a unique key for the triangle

            if triangle_key not in triangle_tracker:
                vertices = unique_vertices[list(triangle_key)]
                area = calculate_triangle_area(vertices)
                areas.append(area)

                # Track this triangle's index
                triangle_tracker[triangle_key] = len(areas) - 1

            # Fill in rwgs for area and positions indices
            rwgs[i, 4 + (j - 2)] = triangle_tracker[triangle_key]  # Triangle index

    areas = np.array(areas)

    dunavant_positions = np.zeros((N,2,dunavant_values.shape[0],3))
    for n in range(N):
        for t in range(2):
            dunavant_positions[n,t] = np.dot(dunavant_values,unique_vertices[rwgs[n,[0,1,t+2]]])
    
    return unique_vertices, rwgs, areas, dunavant_positions

def check_triangle_pair_singularity(rwgs):
    # Step 1: Expand rwgs to triangle vertices
    # rwgs structure: [common_vertex_0, common_vertex_1, non_common_vertex_1st_triangle, non_common_vertex_2nd_triangle]
    N = rwgs.shape[0]
    triangle_pairs = np.ones((N, N, 4, 3), dtype=int)  # Shape: (N, N, 4 pairs, 3 vertices per triangle)
    
    # Extract vertices indices for each triangle in every RWG pair
    for n in range(N):
        for i in range(n,N):
            triangle_pairs[n, i, 0, :] = [rwgs[n, 0], rwgs[n, 1], rwgs[n, 2]]  # Triangle 1 of RWG n
            triangle_pairs[n, i, 1, :] = [rwgs[n, 0], rwgs[n, 1], rwgs[n, 3]]  # Triangle 2 of RWG n
            triangle_pairs[n, i, 2, :] = [rwgs[i, 0], rwgs[i, 1], rwgs[i, 2]]  # Triangle 1 of RWG i
            triangle_pairs[n, i, 3, :] = [rwgs[i, 0], rwgs[i, 1], rwgs[i, 3]]  # Triangle 2 of RWG i

    # Step 2: Vectorized Singularity Check
    singularities_map = np.zeros((N, N, 4), dtype=bool)  # Shape: (N, N, 4 pairs)
    
    for t1 in range(2):
        for t2 in range(2, 4):
            t1_vertices = triangle_pairs[:, :, t1, :].reshape(N, N, 1, 3)
            t2_vertices = triangle_pairs[:, :, t2, :].reshape(N, N, 3, 1)
            # Check if any vertex is shared between t1 and t2 triangles of each RWG pair
            singularities_map[:, :, 2*t1+t2-2] = np.any(np.any(t1_vertices == t2_vertices, axis=-1), axis=-1)

    return singularities_map

class Mesh:
    """Class holding all triangles from a mesh file"""

    def __init__(self, path, verbose=False, properties=None):

        if verbose:
            print("Opening mesh file...")

        if not re.match(r"(.*\.ply$)|(.*\.obj$)|(.*\.dat$)", path, re.IGNORECASE):
            raise NotImplementedError("Currently only .dat .obj or .ply files are supported")

        filepath, extension = os.path.splitext(path)
        extension = extension.lower()

        if extension == ".obj":
            raise NotImplementedError(".obj still not supported")

        # Open mesh file
        if not isinstance(path, str):
            raise TypeError("Mesh path should be a string")
        file = open(path, 'r')  # can raise FileNotFoundError

        if verbose:
            print("Parsing mesh file header...")

        lines_raw = file.readlines()

        # Just to stop the IDE from complaining
        header_length = None
        spacing = None
        splitter = None
        n_vertices = None
        n_faces = None
        n_objects = None
        normal = None
        offset = None
        current_object = None
        lines = []
        vertices = []
        normals = []
        triangles = []
        edge_objects = []
        self.objects = []
        self.array1_reshape = []
        self.plot_limits = []
        edges = []

        if extension == ".ply":
            # Removing all comments from the .ply file
            for line in lines_raw:
                if not line[0:8] == "comment ":
                    lines.append(line[0:len(line) - 1])

            # Asserting a correct .ply file header
            header_length = 12
            spacing = 0
            splitter = " "
            offset = 0

            assert lines[0] == "ply", "Magic number does not match, current magic number: " + lines[0]
            assert lines[1] == "format ascii 1.0", "Format not ASCII 1.0, received: " + lines[1]
            assert lines[2][
                   0:14] == "element vertex", "Element vertex expected in line 3 of the mesh file, but received: " + \
                                              lines[2]
            assert int(lines[2].split(splitter)[
                           2]), "Amount of vertex is not a number, or incorrectly defined in line 3 of mesh file: " + \
                                line[2]
            assert lines[
                       3] == "property float x", "property float x expected in line 4 of the mesh file, but received: " + \
                                                 lines[3]
            assert lines[
                       4] == "property float y", "property float x expected in line 5 of the mesh file, but received: " + \
                                                 lines[4]
            assert lines[
                       5] == "property float z", "property float x expected in line 6 of the mesh file, but received: " + \
                                                 lines[5]
            assert lines[
                       6] == "property float nx", "property float x expected in line 7 of the mesh file, but received: " + \
                                                  lines[6]
            assert lines[
                       7] == "property float ny", "property float x expected in line 8 of the mesh file, but received: " + \
                                                  lines[7]
            assert lines[
                       8] == "property float nz", "property float x expected in line 9 of the mesh file, but received: " + \
                                                  lines[8]
            assert lines[9][
                   0:12] == "element face", "Element vertex expected in line 10 of the mesh file, but received: " + \
                                            lines[9]
            assert int(lines[9].split(splitter)[
                           2]), "Amount of vertex is not a number, or incorrectly defined in line 10  of mesh file: " + \
                                line[9]
            assert lines[
                       10] == "property list uchar uint vertex_indices", "property list uchar uint vertex_indices expected, but received: " + \
                                                                         lines[10]
            assert lines[11] == "end_header", "header end expected, but received: " + lines[11]

            n_vertices = int(lines[2].split(splitter)[2])
            n_faces = int(lines[9].split(splitter)[2])
            n_objects = 1

        elif extension == ".dat":
            # getting all relevant header information from .dat file
            for line in lines_raw:
                lines.append(line.rstrip("\r\n"))

            header_length = 14
            spacing = 5
            splitter = None
            offset = 1

            # asserting correct header file format
            spacer = r"-*"
            error = ".dat file header is incorrect"

            for i in [0, 2, 10, 12]:
                assert re.match(spacer, lines[i]), error
            assert lines[1] == "Dimensions", error + ": " + lines[1] + " is not \"Dimensions\""
            assert lines[3] == "Number of nodes", error + ": " + lines[3] + " is not \"Number of nodes\""
            assert lines[5] == "Number of elements", error + ": " + lines[5] + " is not \"Number of elements\""
            assert lines[7] == "Number of objects", error + ": " + lines[7] + " is not \"Number of objects\""
            assert lines[9] == "", error + ": line 10 is expected to be blank, but read: " + lines[9]
            assert lines[13].split(splitter) == ["Nodenr", "x-coordinate", "y-coordinate", "z-coordinate"], error

            n_vertices = int(lines[4].split(splitter)[0])
            n_faces = int(lines[6].split(splitter)[0])
            n_objects = int(lines[8].split(splitter)[0])

            for i in [1, 3]:
                assert re.match(spacer, lines[header_length + n_vertices + i]), error

            assert lines[header_length + n_vertices] == "", error + ": line 10 is expected to be blank, but read: " + \
                                                            lines[header_length + n_vertices]
            assert lines[header_length + n_vertices + 2] == "Elements", error + ": " + lines[
                header_length + n_vertices + 2] + " is not \"Elements\""
            assert lines[header_length + n_vertices + 4].split(splitter) == ["Elementnr", "node1", "node2", "node3",
                                                                             "Objectnr"], error

        assert len(lines) == header_length + n_vertices + spacing + n_faces, "length of file is not as expected"
        assert properties is None or len(
            properties) == n_objects, "list of properties does not match the amount of objects"

        for i in range(n_objects):
            self.objects.append(Object())

        if verbose:
            print("Creating all vertices...")

        # scan all lines containing vertices and store them, for .ply: also store associated normals
        for line in lines[header_length:header_length + n_vertices]:
            line = line.split(splitter)
            line = [float(number) for number in line]

            if extension == ".ply":
                vertices.append(line[0:3])
                normals.append(line[3:6])
            elif extension == ".dat":
                vertices.append(line[1:4])

        if verbose:
            print("Creating all edges and faces...")

        # scan all lines containing the face information
        for line in lines[header_length + n_vertices + spacing:header_length + n_vertices + spacing + n_faces]:
            line = line.split(splitter)
            line = [int(number) for number in line]

            # .ply files can also contain other types of polygons, assert that only 3 vertices are used in faces
            if extension == ".ply":
                assert line[0] == 3, "Mesh contains a non-triangle face"

            current_edges = []
            current_edges_obj = []
            current_new_edges = []

            # check all different combinations of edges between 2 vertices in the triangle
            for i in [[1, 2], [2, 3], [3, 1]]:
                vertex1 = vertices[line[i[0]] - offset]
                vertex2 = vertices[line[i[1]] - offset]
                current_edges = np.append(current_edges, np.array(vertex1))

                # check if the edge already exist (or the reverse edge), if not, make a new one
                if [vertex1, vertex2] in edges:
                    current_edges_obj.append(edge_objects[edges.index([vertex1, vertex2])])
                elif [vertex2, vertex1] in edges:
                    current_edges_obj.append(edge_objects[edges.index([vertex2, vertex1])])
                else:
                    edge = Edge(np.array(vertex1), np.array(vertex2))
                    edges.append([vertex1, vertex2])
                    edge_objects.append(edge)
                    current_edges_obj.append(edge)
                    current_new_edges.append(edge)

            # calculate or look up the normal of the triangle
            if extension == ".ply":
                assert normals[line[1]] == normals[line[2]] and normals[line[2]] == normals[
                    line[3]], "Normals do not match"
                normal = np.array(normals[line[1]])
            elif extension == ".dat":
                # .dat file have right handed way of assuming the normal vector
                edge1 = np.array(vertices[line[2] - offset]) - np.array(vertices[line[1] - offset])
                edge2 = np.array(vertices[line[3] - offset]) - np.array(vertices[line[1] - offset])
                normal = np.cross(edge1, edge2)
                normal = normal / np.linalg.norm(normal)

            # make a new triangle object
            if extension == ".ply":
                current_object = 0
            elif extension == ".dat":
                current_object = line[4] - offset

            current_edges = np.reshape(current_edges, (3, 3))
            current_triangle = Triangle(current_edges, normal)
            self.objects[current_object].add_triangle(current_triangle)

            for edge in current_new_edges:
                self.objects[current_object].add_edge(edge)

            # associate all triangle edges with the new triangle
            for edge in current_edges_obj:
                edge.add_triangle(current_triangle)

        if verbose:
            print("Linking neighboring triangles...")

        # associate all triangles with the correct neighbor
        for edge in edge_objects:
            assert len(edge.triangles) == 2, "There exist edges in the mesh that do not belong to 2 triangles"
            triangle1 = edge.triangles[0]
            triangle2 = edge.triangles[1]
            triangle1.add_neighbor(triangle2)
            triangle2.add_neighbor(triangle1)

        # just a check: make sure that all triangles have 3 neighbors
        for triangle in triangles:
            assert len(triangle.neighbors) == 3, "There exist triangles that do not have 3 neighbors"

        self.n_edges = len(edge_objects)

    def __str__(self):
        return "Mesh structure containing " + str(len(self.objects)) + " objects."

    def prepare_plot_data(self, normals=False):
        # Initialize lists to accumulate all coordinates for determining plot limits
        all_x, all_y, all_z = [], [], []
        
        self.plot_data = []
        edge_lengths = []  # Collect all edge lengths to determine padding based on the largest edge length

        for obj in self.objects:
            vertices_data, mean_normals = [], []
            
            for triangle in obj.triangles:
                vertices = np.array([point for point in triangle.points])
                all_x.extend(vertices[:, 0])
                all_y.extend(vertices[:, 1])
                all_z.extend(vertices[:, 2])
                
                if normals:
                    mean_x, mean_y, mean_z = np.mean(vertices[:, 0]), np.mean(vertices[:, 1]), np.mean(vertices[:, 2])
                    normal = triangle.normal
                    edge_length = np.sqrt(np.sum(np.square(obj.edges[0].vertex1 - obj.edges[0].vertex2)))
                    edge_lengths.append(edge_length)
                    mean_normals.append((mean_x, mean_y, mean_z, *normal, edge_length * 0.5))

                vertices_data.append({'vertices': vertices, 'normal': triangle.normal if normals else None})

            self.plot_data.append({'vertices_data': vertices_data, 'mean_normals': mean_normals if normals else None})

        if edge_lengths:
            max_edge_length = max(edge_lengths)
        else:
            max_edge_length = 1  # Default to 1 if no edges (unlikely, but safe fallback)

        # Determine the plot limits based on accumulated coordinates
        boundary_padding = max_edge_length * 0.1
        self.plot_limits = {
            'x': [min(all_x) - boundary_padding, max(all_x) + boundary_padding],
            'y': [min(all_y) - boundary_padding, max(all_y) + boundary_padding],
            'z': [min(all_z) - boundary_padding, max(all_z) + boundary_padding]
        }

        # Prepare self.array1_reshape based on all vertices
        if all_x and all_y and all_z:  # Ensure there is data to process
            all_vertices = np.array([all_x, all_y, all_z]).T
            self.array1_reshape = all_vertices.reshape(-1, 9)  # Assuming this matches the original intent
        
    def getmatrix(self):

        # TRIANGLE AREAS

        # Line 441 - 446: make empty matrices which can be filled with values later on
        s = (len(self.array1_reshape), 3)
        vector1 = np.zeros(s)
        vector2 = np.zeros(s)
        cross = np.zeros(s)
        magnitude = np.zeros(len(self.array1_reshape))
        area_triangle = np.zeros(len(self.array1_reshape))

        # Make two 3xM matrices with M = amount of triangles where for vector 1: column 1 = x_2 - x_1,
        # column 2 = y_2 - y_1, column 3 = z_2 - z_1. And for vector2: column 1 = x_3 - x_1, column 2 = y_3 - y_1,
        # column 3 = z_3 - z_1. This is needed for further calculation of the triangle area.
        for l in range(len(self.array1_reshape)):
            vector1[l, 0] = self.array1_reshape[l, 3] - self.array1_reshape[l, 0]
            vector1[l, 1] = self.array1_reshape[l, 4] - self.array1_reshape[l, 1]
            vector1[l, 2] = self.array1_reshape[l, 5] - self.array1_reshape[l, 2]
            vector2[l, 0] = self.array1_reshape[l, 6] - self.array1_reshape[l, 0]
            vector2[l, 1] = self.array1_reshape[l, 7] - self.array1_reshape[l, 1]
            vector2[l, 2] = self.array1_reshape[l, 8] - self.array1_reshape[l, 2]

            # Take the cross product of the vectors which were just constructed.
            cross[l] = np.cross(vector1[l], vector2[l])
            # Take the magnitude of the cross product.
            magnitude[l] = np.linalg.norm(cross[l])
            # The half of this magnitude is equal to the area triangle.
            area_triangle[l] = 0.5 * magnitude[l]

        # INNER EDGES

        # Line 468 - 475: Defining empty variables, arrays and matrices for the inner edge calculation.
        count = np.zeros(len(self.array1_reshape))
        count2 = 0
        s2 = (len(self.array1_reshape) * 3 + 1, 6)
        s3 = (1, 6)
        inner_edges = np.zeros(s2)
        no_inner_edges = np.zeros(s2)
        temp = np.zeros(s3)

        # Inner edge for-loops.
        for k in range(len(self.array1_reshape)):  # Amount triangles = len(self.array1_reshape).
            for m in range(len(self.array1_reshape)):  # Update the temporary variables to be zero for each iteration.
                count2 = 0  # Variable used to count if two coordinates of a triangle are similar.
                temp = np.zeros(s3) # Variable to store equal coordinate values.
                if k != m and count[k] < 3 and count[m] < 3:  # Check if both triangles which will be compared don't
                    # already have 3 found inner edges.
                    if np.array_equal(self.array1_reshape[k, 0:3], self.array1_reshape[m, 0:3]): # Check if the first
                        # coordinates of the triangle are equal.
                        count2 = count2 + 1 # Update the count, 1 more similar coordinate found!

                        if np.array_equal(temp[0, 0:3], np.array([0, 0, 0])) and count2 == 1:   # Check if there's already a coordinate
                            # stored in the temporary matrix.
                            temp[0, 0:3] = self.array1_reshape[k, 0:3]  # If the first 3 values are still empty, store
                            # the found coordinates here.
                        else:
                            temp[0, 3:6] = self.array1_reshape[k, 0:3]  # Otherwise, store them in the other 3 columns.

                    # This process is repeated for the other coordinates of the triangle until all coordinates are
                    # compared.
                    if np.array_equal(self.array1_reshape[k, 0:3], self.array1_reshape[m, 3:6]):
                        count2 = count2 + 1

                        if np.array_equal(temp[0, 0:3], np.array([0, 0, 0])) and count2 == 1:
                            temp[0, 0:3] = self.array1_reshape[k, 0:3]
                        else:
                            temp[0, 3:6] = self.array1_reshape[k, 0:3]

                    if np.array_equal(self.array1_reshape[k, 0:3], self.array1_reshape[m, 6:9]):
                        count2 = count2 + 1

                        if np.array_equal(temp[0, 0:3], np.array([0, 0, 0])) and count2 == 1:
                            temp[0, 0:3] = self.array1_reshape[k, 0:3]
                        else:
                            temp[0, 3:6] = self.array1_reshape[k, 0:3]

                    if np.array_equal(self.array1_reshape[k, 3:6], self.array1_reshape[m, 0:3]):
                        count2 = count2 + 1

                        if np.array_equal(temp[0, 0:3], np.array([0, 0, 0])) and count2 == 1:
                            temp[0, 0:3] = self.array1_reshape[k, 3:6]
                        else:
                            temp[0, 3:6] = self.array1_reshape[k, 3:6]

                    if np.array_equal(self.array1_reshape[k, 3:6], self.array1_reshape[m, 3:6]):
                        count2 = count2 + 1

                        if np.array_equal(temp[0, 0:3], np.array([0, 0, 0])) and count2 == 1:
                            temp[0, 0:3] = self.array1_reshape[k, 3:6]
                        else:
                            temp[0, 3:6] = self.array1_reshape[k, 3:6]

                    if np.array_equal(self.array1_reshape[k, 3:6], self.array1_reshape[m, 6:9]):
                        count2 = count2 + 1

                        if np.array_equal(temp[0, 0:3], np.array([0, 0, 0])) and count2 == 1:
                            temp[0, 0:3] = self.array1_reshape[k, 3:6]
                        else:
                            temp[0, 3:6] = self.array1_reshape[k, 3:6]

                    if np.array_equal(self.array1_reshape[k, 6:9], self.array1_reshape[m, 0:3]):
                        count2 = count2 + 1

                        if np.array_equal(temp[0, 0:3], np.array([0, 0, 0])) and count2 == 1:
                            temp[0, 0:3] = self.array1_reshape[k, 6:9]
                        else:
                            temp[0, 3:6] = self.array1_reshape[k, 6:9]

                    if np.array_equal(self.array1_reshape[k, 6:9], self.array1_reshape[m, 3:6]):
                        count2 = count2 + 1

                        if np.array_equal(temp[0, 0:3], np.array([0, 0, 0])) and count2 == 1:
                            temp[0, 0:3] = self.array1_reshape[k, 6:9]
                        else:
                            temp[0, 3:6] = self.array1_reshape[k, 6:9]

                    if np.array_equal(self.array1_reshape[k, 6:9], self.array1_reshape[m, 6:9]):
                        count2 = count2 + 1

                        if np.array_equal(temp[0, 0:3], np.array([0, 0, 0])) and count2 == 1:
                            temp[0, 0:3] = self.array1_reshape[k, 6:9]
                        else:
                            temp[0, 3:6] = self.array1_reshape[k, 6:9]

                    # Now all the coordinates of the 2 triangles are compared. If 2 similar coordinates are found, we
                    # update the count variable for both triangles. This is the variable which counts how many inner
                    # edges were found. If 3 inner edges are found for a triangle, all inner edges are found (since a
                    # triangle only has 3 edges in total).
                    if count2 == 2:
                        count[k] = count[k] + 1
                        count[m] = count[m] + 1

                        # Since the amount of inner edges N is larger than the amount of triangles M (N = 1.5*M), we
                        # must be careful with placing values in the matrix. Here we update the parameter z to keep
                        # track of where in the matrix the inner edges must be stored, such that no previous
                        # found inner edges are overwritten.
                        if count[k] == 1:
                            z = k * 3
                        elif count[k] == 2:
                            z = k * 3 + 1
                        elif count[k] == 3:
                            z = k * 3 + 2

                        # Place the temp array in a row of the inner_edges matrix. Note that we use z here.
                        inner_edges[z] = temp

                        # These if statements are for updating the no_inner_edge matrix. This is a matrix which contains
                        # the coordinates of the 2 triangles which have an inner edge together, but are not part of the
                        # inner edge. So the remaining coordinates of the triangles. It is a bit hard to explain in
                        # words, so if this idea is not clear to you, we've made a figure of how the coordinates are
                        # stored and placed it in the report!
                        # We first have if statements for triangle k, and then we have the same if statements for
                        # triangle m. These are the if statements for triangle k:
                        if not np.array_equal(self.array1_reshape[k, 0:3],
                                              inner_edges[z, 0:3]) and not np.array_equal(
                            self.array1_reshape[k, 0:3], inner_edges[z, 3:6]):  # If the first coordinate of the
                            # triangle k is not equal to the coordinates stored in the inner_edge matrix-row, then this
                            # is the missing coordinate of triangle k!
                            no_inner_edges[z, 0:3] = self.array1_reshape[k, 0:3]    # Put the missing coordinates in
                            # the no_inner_edges matrix.
                        # If one of the coordinates was the same, then we continue with the following if statements
                        # which work the same as the first if statement, but now compare the other coordinates of the
                        # triangle k.
                        elif not np.array_equal(self.array1_reshape[k, 3:6],
                                                inner_edges[z, 0:3]) and not np.array_equal(
                            self.array1_reshape[k, 3:6], inner_edges[z, 3:6]):
                            no_inner_edges[z, 0:3] = self.array1_reshape[k, 3:6]
                        elif not np.array_equal(self.array1_reshape[k, 6:9],
                                                inner_edges[z, 0:3]) and not np.array_equal(
                            self.array1_reshape[k, 6:9], inner_edges[z, 3:6]):
                            no_inner_edges[z, 0:3] = self.array1_reshape[k, 6:9]

                        # These if statements are for triangle m. They will be stored in the last three columns of
                        # the no_inner_edges matrix instead of the first three columns where the coordinate of triangle
                        # k is stored. Furthermore, they work exactly the same as the earlier if statements for
                        # triangle k.
                        if not np.array_equal(self.array1_reshape[m, 0:3],
                                              inner_edges[z, 0:3]) and not np.array_equal(
                            self.array1_reshape[m, 0:3], inner_edges[z, 3:6]):
                            no_inner_edges[z, 3:6] = self.array1_reshape[m, 0:3]
                        elif not np.array_equal(self.array1_reshape[m, 3:6],
                                                inner_edges[z, 0:3]) and not np.array_equal(
                            self.array1_reshape[m, 3:6], inner_edges[z, 3:6]):
                            no_inner_edges[z, 3:6] = self.array1_reshape[m, 3:6]
                        elif not np.array_equal(self.array1_reshape[m, 6:9],
                                                inner_edges[z, 0:3]) and not np.array_equal(
                            self.array1_reshape[m, 6:9], inner_edges[z, 3:6]):
                            no_inner_edges[z, 3:6] = self.array1_reshape[m, 6:9]

                    # If the count2 is not equal to 2. Then NO inner edge is found. Update the temporary variables to be
                    # zero again, such that the for loop can be done correctly again for another triangle combination.
                    else:
                        count2 = 0
                        temp = np.zeros(s3)

        # Sometimes rows with only zero's are constructed (because of the z's). These rows are deleted by using the
        # following comments.
        inner_edges = inner_edges[~np.all(inner_edges == 0, axis=1)]
        no_inner_edges = no_inner_edges[~np.all(no_inner_edges == 0, axis=1)]

        # LENGTH OF INNER EDGES

        # Make a zero array to store the length of the inner edges in.
        s4 = (len(inner_edges))
        inner_edges_length = np.zeros(s4)

        # In this for loop the length of each inner edge is calculated using the standard equation for this.
        for q in range(len(inner_edges)):
            inner_edges_length[q] = np.sqrt(
                (inner_edges[q, 0] - inner_edges[q, 3]) ** 2 + (inner_edges[q, 1] - inner_edges[q, 4]) ** 2 + (
                        inner_edges[q, 2] - inner_edges[q, 5]) ** 2)

        return inner_edges_length, inner_edges, no_inner_edges, area_triangle

    def plot_current(self, currents_on_edge, e_vertex, dir,polarization): # Function to plot the current/Efield values over the edges

        # Create a colour map to assign colours to the current value
        mappie = cm.get_cmap('plasma')

        # Use the absolute current and normalize it
        currents_on_edge = np.absolute(currents_on_edge)
        currents_on_edge = (currents_on_edge - np.min(currents_on_edge)) / (np.max(currents_on_edge))# - np.min(currents_on_edge))

        # Create a figure to plot all edges and directions of the incoming field in
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        # Plot each edge with its colour based on the colour map
        for i in range(len(e_vertex)):
            ax.plot([e_vertex[i,0], e_vertex[i,3]], [e_vertex[i,1], e_vertex[i,4]], [e_vertex[i,2], e_vertex[i,5]], c=mappie(float(currents_on_edge[i])), linewidth=2, linestyle='solid')

        # Calculate the k-vector and polarization of the incoming field to also plot this for a better understanding of the situation

        # Determine wave vector using input angles phi (dir[0]) and theta (dir[1])
        kx = np.sin(dir[0]) * np.cos(dir[1])
        ky = np.sin(dir[0]) * np.sin(dir[1])
        kz = np.cos(dir[0])

        # Find values for a vector with a slightly different phi value
        kx_phi = np.sin(dir[0] + 0.01) * np.cos(dir[1])
        ky_phi = np.sin(dir[0] + 0.01) * np.sin(dir[1])
        kz_phi = np.cos(dir[0] + 0.01)
        # Cross product of the wave vector and the changed phi vector to find a perpendicular vector
        Vect1 = np.cross(np.array([kx_phi, ky_phi, kz_phi]), np.array([kx, ky, kz]))

        # Find values for a vector with a slightly different theta value
        kx_th = np.sin(dir[0]) * np.cos(dir[1] + 0.01)
        ky_th = np.sin(dir[0]) * np.sin(dir[1] + 0.01)
        kz_th = np.cos(dir[0])
        # Cross product of the wave vector and the changed theta vector to find a perpendicular vector
        Vect2 = np.cross(np.array([kx_th, ky_th, kz_th]), np.array([kx, ky, kz]))

        # Normalize the perpendicular vectors
        pol_vect1 = np.multiply((1 / np.linalg.norm(Vect1)), Vect1)

        pol_vect2 = np.multiply((1 / np.linalg.norm(Vect2)), Vect2)

        # Allow user to choose polarization by using input angle between the two normalized polarization vectors
        polarization = np.sin(polarization) * pol_vect2 + np.cos(polarization) * pol_vect1
        polarization = np.multiply((1 / np.linalg.norm(polarization)), polarization)

        # Plot the k-vector in green, from the origin in the direction of the wave, and the E field polarization in red.
        ax.quiver(0, 0, 0, -kx, -ky, -kz, length=1, color='green') #plot -k because of the -j convention in the exponent
        ax.quiver(0,0,0,polarization[0],polarization[1],polarization[2], length=0.3, color='red')

        # Plot the colormap
        a = matplotlib.colors.Normalize(vmin=min(currents_on_edge), vmax=max(currents_on_edge))
        mapp = plt.cm.ScalarMappable(cmap=mappie, norm=a)
        mapp.set_array([])

        

        # Set the boundaries based on earlier calculated limits (see mesh.plot())
        ax.set_xlim(*self.plot_limits['x'])
        ax.set_ylim(*self.plot_limits['y'])
        ax.set_zlim(*self.plot_limits['z'])
        ax.set_xlabel('X axis')
        ax.set_ylabel('Y axis')
        ax.set_zlabel('Z axis')

        fig.colorbar(mapp, ax=ax, ticks=np.linspace(min(currents_on_edge[0]), max(currents_on_edge[0]), 20))
        plt.show()