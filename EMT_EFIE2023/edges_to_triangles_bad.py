import numpy as np

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

def process_r_vect(r_vect):
    # Step 1: Flatten r_vect to get all vertices, find unique vertices
    all_vertices = r_vect.reshape(-1, 3)
    unique_vertices, inverse_indices = np.unique(all_vertices, axis=0, return_inverse=True)
    
    # Step 2: Re-map r_vect vertices to unique vertices indices
    mapped_indices = inverse_indices.reshape(-1, 4)  # Remap each RWG vertices to unique vertices indices
    
    # Step 3: Construct triangles and RWGs arrays
    # Initialize triangles list to keep track of unique triangles as (index1, index2, index3)
    triangles_dict = {}
    triangles = []  # This will store the triangle vertex indices
    rwgs = np.zeros((len(r_vect), 2), dtype=int)  # This will map each RWG to two triangles
    
    for i, rwg in enumerate(mapped_indices):
        # Define each triangle by its vertices indices, including shared edge for uniqueness
        tri1 = tuple(sorted([rwg[0], rwg[1], rwg[2]]))
        tri2 = tuple(sorted([rwg[0], rwg[1], rwg[3]]))
        
        # Add triangles to the dictionary if not already present, tracking their index
        for tri in [tri1, tri2]:
            if tri not in triangles_dict:
                triangles_dict[tri] = len(triangles)
                triangles.append(list(tri))
                
        # Map RWG to its triangles
        rwgs[i] = [triangles_dict[tri1], triangles_dict[tri2]]
    
    triangles = np.array(triangles)
    
    return unique_vertices, triangles, rwgs


unique_vertices, triangles, rwgs = process_r_vect(r_vect)

# Verifying the shapes and contents
print("Unique vertices shape:", unique_vertices.shape)
print("Triangles shape:", triangles.shape)
print("RWGs shape:", rwgs.shape)
print("RWGs:\n", rwgs)

print(unique_vertices[triangles[rwgs[0]]])
print(r_vect[0])