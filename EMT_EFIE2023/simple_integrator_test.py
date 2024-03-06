import numpy as np
class integrator_test():
    def initialize(self):
        print("Initialized test integrator")
        print("This integrator was not well tested!!!!")
        print("Results might be wrong \n")

    def int_test(self,fun,ding1):
        midpoint = (ding1[0] + ding1[1] + ding1[2])/3.
        funval = fun(midpoint)
        direction = midpoint - ding1[2]
        lvec =ding1[1] - ding1[0]
        lleng = np.sqrt(lvec[0]**2 + lvec[1] **2 + lvec[[2]]**2)
        #direction /= 2*np.sqrt(direction[0]**2+direction[1]**2+ direction[2]**2)
        #print(direction)
        #print(funval)
        #print(lleng)
        fundir = direction[0]*funval[0] + direction[1]*funval[1] + direction[2]*funval[2]
        side1 = ding1[1] - ding1[0]
        side2 = ding1[2] - ding1[0]
        #area=.5*np.sqrt((side1[0]*side2[1] - side1[1]* side2[0])**2 + \
        #        (side1[1]*side2[2] - side1[2]* side2[1])**2 +
        #        (side1[2]*side2[0] - side1[0]* side2[2])**2 )

        return fundir * lleng /(2)


    def int_triangle(self,fun,ding1):
        midpoint = (ding1[0] + ding1[1] + ding1[2])/3.
        funval = fun(midpoint)
        side1 = ding1[1] - ding1[0]
        side2 = ding1[2] - ding1[0]
        area=.5*np.sqrt((side1[0]*side2[1] - side1[1]* side2[0])**2 + \
                (side1[1]*side2[2] - side1[2]* side2[1])**2 +
                (side1[2]*side2[0] - side1[0]* side2[2])**2 )
        return funval*area

