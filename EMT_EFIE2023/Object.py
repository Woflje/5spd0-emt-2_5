########################################################################################
#
# This code is part of the Eindhoven ElectroMagnetics Solver
#
#   Function name: Object class
#
#   Description: Creates Object object and can output RWGs
#
#   Input: (object properties)
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

import numpy as np


class Object:
    """Object containing edge information, properties can be added in the future"""

    def __init__(self, properties=None):
        self.edges = []
        self.triangles = []
        self.properties = properties
        self.rwg = []

    def __str__(self):
        return "Object containing "+str(len(self.edges))+" edges."

    def add_edge(self, edge):
        self.edges.append(edge)

    def add_triangle(self, triangle):
        self.triangles.append(triangle)

    def get_rwgs(self):
        if len(self.rwg) == 0:
            for edge in self.edges:
                self.rwg = np.append(self.rwg, np.append(edge.vertex1, edge.vertex2))
                self.rwg = np.append(self.rwg, np.append(edge.triangles[0].points, edge.triangles[1].points))
                self.rwg = np.append(self.rwg, np.append(edge.triangles[0].normal, edge.triangles[1].normal))

            self.rwg = np.reshape(self.rwg, (-1, 10, 3))
        return self.rwg
