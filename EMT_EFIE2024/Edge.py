########################################################################################
#
# This code is part of the Eindhoven ElectroMagnetics Solver
#
#   Function name: Edge class
#
#   Description: Creates edge object and stores relevant triangle information
#
#   Input: numpy vector of first and second vertex
#
#   Output: none
#
#   Documentation: see documentation
#
#   Author name(s): Sjoerd Aker, Leroy Driessen
#
#   Date: 8-3-2020
#
# The above authors declare that the following code is free of plagiarism.
#
# Maintainer/Supervisor: Roeland J. Dilz (r.dilz@tue.nl)
#
# This code will be published under GPL v 3.0 https://www.gnu.org/licenses/gpl-3.0.en.html
#
#######################################################################################


class Edge:
    """Class holding edge information"""

    def __init__(self, vertex1, vertex2):
        self.vertex1 = vertex1
        self.vertex2 = vertex2
        self.triangles = []

    def __str__(self):
        return str(self.vertex1)+"->"+str(self.vertex2)

    def vector(self):
        return self.vertex2 - self.vertex1

    def add_triangle(self, triangle):
        if len(self.triangles) >= 2:
            raise ValueError("Edge already has two associated triangles")
        self.triangles.append(triangle)

    def get_associated_triangles(self):
        return self.triangles

    def get_opposite_triangle(self, triangle):
        if triangle not in self.triangles:
            raise ValueError("Triangle does not contain this edge")
        if triangle == self.triangles[0]:
            return self.triangles[1]
        else:
            return self.triangles[0]
