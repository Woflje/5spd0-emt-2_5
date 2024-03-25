import numpy as np

def restructure_mesh_rwg_data(r_vect, dunavant_values):
    # Flatten r_vect to get all vertices and deduplicate to get unique vertices
    all_vertices = r_vect.reshape(-1, 3)
    unique_vertices, inverse_indices = np.unique(all_vertices, axis=0, return_inverse=True)
    
    # Map vertices in r_vect to unique vertices indices
    mapped_indices = inverse_indices.reshape(r_vect.shape[0], r_vect.shape[1])

    # Initialize structures
    rwgs = np.zeros((len(r_vect), 6), dtype=int)  # Mapping indices
    areas = []  # To store areas of each triangle
    positions = []  # To store positions (P, 3) for each triangle

    triangle_tracker = {}  # Tracks whether a triangle has been processed

    for i, rwg in enumerate(mapped_indices):
        rwgs[i, :4] = rwg[:4]  # Fill in rwgs for vertices
        
        for j in range(2, 4):  # Loop over the unique vertices for each triangle
            triangle_key = tuple(sorted([rwg[0], rwg[1], rwg[j]]))  # Create a unique key for the triangle

            if triangle_key not in triangle_tracker:
                vertices = unique_vertices[list(triangle_key)]
                area = calculate_triangle_area(vertices)
                areas.append(area)
                # Directly use the corresponding dunavant positions data, assuming it matches the processing order
                positions.append(np.dot(dunavant_values, vertices))

                # Track this triangle's index
                triangle_tracker[triangle_key] = len(areas) - 1

            # Fill in rwgs for area and positions indices
            rwgs[i, 4 + (j - 2)] = triangle_tracker[triangle_key]  # Triangle index

    areas = np.array(areas)
    dunavant_positions = np.array(positions)  # Convert to a structured NumPy array

    return unique_vertices, rwgs, areas, dunavant_positions

def restructure_mesh_rwg_data2(r_vect, dunavant_values):
    # Flatten r_vect to get all vertices and deduplicate to get unique vertices
    all_vertices = r_vect.reshape(-1, 3)
    unique_vertices, inverse_indices = np.unique(all_vertices, axis=0, return_inverse=True)
    
    # Map vertices in r_vect to unique vertices indices
    mapped_indices = inverse_indices.reshape(r_vect.shape[0], r_vect.shape[1])

    # Initialize structures
    rwgs = np.zeros((len(r_vect), 6), dtype=int)  # Mapping indices

    triangle_tracker = {}  # Tracks whether a triangle has been processed

    triangle_id = 0

    for i, rwg in enumerate(mapped_indices):
        rwgs[i, :4] = rwg[:4]  # Fill in rwgs for vertices
        
        for j in range(2, 4):  # Loop over the unique vertices for each triangle
            triangle_key = tuple(sorted([rwg[0], rwg[1], rwg[j]]))  # Create a unique key for the triangle

            if triangle_key not in triangle_tracker:

                # Track this triangle's index
                triangle_tracker[triangle_key] = triangle_id
                triangle_id += 1

            # Fill in rwgs for area and positions indices
            rwgs[i, 4 + (j - 2)] = triangle_tracker[triangle_key]  # Triangle index

    triangle_tracker = set()

    areas = []
    positions = []

    for n in range(0,len(rwgs)):
        for i in range(4,6):
            triangle_id = rwgs[n,i]
            if triangle_id not in triangle_tracker:
                vertices = unique_vertices[rwgs[n,[0,1,i-2]]]
                areas.append(calculate_triangle_area(vertices))
                positions.append(np.dot(dunavant_values, vertices))
                triangle_tracker.add(triangle_id)

    areas = np.array(areas)
    dunavant_positions = np.array(positions)  # Convert to a structured NumPy array

    return unique_vertices, rwgs, areas, dunavant_positions

def calculate_triangle_area(vertices):
    a = np.linalg.norm(vertices[0] - vertices[1])
    b = np.linalg.norm(vertices[1] - vertices[2])
    c = np.linalg.norm(vertices[2] - vertices[0])
    s = (a + b + c) / 2
    return np.sqrt(s * (s - a) * (s - b) * (s - c))

def check_triangle_pair_singularity(rwgs):
    N = rwgs.shape[0]
    # Initialize an array to hold singularity results for each triangle pair
    # Shape (N, N, 4) where the last dimension holds boolean flags for each triangle pair comparison
    singularities = np.ones((N, N, 4), dtype=bool)
    
    # Generate all triangle vertex indices from rwgs
    # Triangle 1 indices: common edge + non-common vertex of the first triangle
    # Triangle 2 indices: common edge + non-common vertex of the second triangle
    triangle_indices = np.zeros((N, 2, 3), dtype=int)
    triangle_indices[:, 0, :2] = rwgs[:, :2]  # Common edge for the first triangle
    triangle_indices[:, 1, :2] = rwgs[:, :2]  # Common edge for the second triangle
    triangle_indices[:, 0, 2] = rwgs[:, 2]  # Non-common vertex for the first triangle
    triangle_indices[:, 1, 2] = rwgs[:, 3]  # Non-common vertex for the second triangle
    
    # Check for shared vertices between each pair of triangles in rwgs
    for n in range(N):
        for i in range(N):
            if i!=n:
                for t1 in range(2):
                    for t2 in range(2):
                        # Extract vertex indices for the triangles being compared
                        vertices_n = triangle_indices[n, t1]
                        vertices_i = triangle_indices[i, t2]
                        # Check if there's any shared vertex
                        shared_vertex = np.intersect1d(vertices_n, vertices_i).size > 0
                        singularities[n, i, 2*t1 + t2] = shared_vertex
    
    return singularities