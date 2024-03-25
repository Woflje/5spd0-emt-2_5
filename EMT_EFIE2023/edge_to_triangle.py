
import numpy as np
from collections import defaultdict

from Mesh import Mesh
from parameters import parameters as param

mesh = Mesh(param["input_mesh_path"], True)
mesh.prepare_plot_data()  # Always run .plot() because the edge sorting is done here
# mesh.plot_objects()

# Load all edges in r_vect with a xyz for each of the 3 vertices of a Triangle
#  r_vect has N elements for each inner edge with 4 vertices: 0 and 1 for the inner edge and 2 and 3 for the two other vertices that make the 2 triangles

[length, e_vertex, other_vertex, area] = mesh.getmatrix()
N = len(e_vertex)


# Assuming 'data' is your original (N, 4, 3) array
r_vect = np.array([e_vertex[:, :3], e_vertex[:, 3:], other_vertex[:, :3], other_vertex[:, 3:]]).transpose(1, 0, 2)


# Step 1: Extract Unique Vertices
vertices, unique_indices = np.unique(r_vect.reshape(-1, 3), axis=0, return_inverse=True)

# Step 2: Correctly Create Triangles Array - Corrected
# Each row in data represents two triangles sharing an edge, so we correctly map these to unique vertex indices.
# Since unique_indices is a flat array mapping each vertex in 'data' to a unique vertex in 'vertices',
# we reshape it back to match 'data's shape to easily identify the triangles.
triangles_indices_flat = unique_indices.reshape(-1, 4)
# No need to double the triangles this time; just map directly.
triangles = np.vstack([triangles_indices_flat[:, [0, 1, 2]], triangles_indices_flat[:, [0, 2, 3]]])

# Step 3: Create Edge to Triangle Mapping
edges_to_triangles = defaultdict(list)

for i, tri in enumerate(triangles):
    # Sorting ensures the edge is always represented in the same order regardless of triangle orientation
    edges = [tuple(sorted((tri[j], tri[(j + 1) % 3]))) for j in range(3)]
    for edge in edges:
        edges_to_triangles[edge].append(i)

# Optionally, convert defaultdict to a regular dict for edges that only appear once
edges_to_triangles = {k: (v if len(v) > 1 else v + [-1]) for k, v in edges_to_triangles.items()}

print("Vertices shape:", vertices.shape)
print("Triangles shape:", triangles.shape)
print("Number of unique edges:", len(edges_to_triangles))

# Example: Get the triangles for an edge
edge = list(edges_to_triangles.keys())[0]  # Example: first edge in the dictionary
print("Triangles sharing the edge", edge, ":", edges_to_triangles[edge])

# To access the coordinates of the vertices of the triangles for a given edge:
triangles_for_edge = edges_to_triangles[edge]
for tri_idx in triangles_for_edge:
    if tri_idx != -1:  # Check if there's an actual triangle (for boundary edges)
        print("Triangle", tri_idx, "vertices:", vertices[triangles[tri_idx]])