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

def calculate_triangle_area(vertices):
    a = np.linalg.norm(vertices[0] - vertices[1])
    b = np.linalg.norm(vertices[1] - vertices[2])
    c = np.linalg.norm(vertices[2] - vertices[0])
    s = (a + b + c) / 2
    return np.sqrt(s * (s - a) * (s - b) * (s - c))

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



