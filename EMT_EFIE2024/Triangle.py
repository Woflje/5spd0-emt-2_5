########################################################################################
#
# This code is part of the Eindhoven ElectroMagnetics Solver
#
#   Function name: Triangle class
#
#   Description: Stores triangles and neighbors
#
#   Input: Points of triangle, normal of triangle
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


class Triangle:
    """Class holding all information about a triangle"""

    def __init__(self, points, normal):
        self.points = points
        self.normal = normal
        self.neighbors = []

    def __str__(self):
        return "Triangle with vertices "+str(self.points[0])+", "+str(self.points[1])+", and "+str(self.points[2])

    def add_neighbor(self, neighbor):
        if len(self.neighbors) >= 3:
            raise ValueError("Triangle already has three neighbors")
        self.neighbors.append(neighbor)
