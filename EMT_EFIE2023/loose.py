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