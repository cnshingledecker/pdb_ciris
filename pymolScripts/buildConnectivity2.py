from pymol import cmd
import numpy as np
        

def scantoview(newview, steps, imgstart):
    imgval = imgstart
    oldview = cmd.get_view()
    delta = np.subtract(newview, oldview)
    delta = delta / steps
    for i in range(1,steps):
        currentview = np.add(oldview, delta * i)
        cmd.set_view(currentview)
        cmd.png("img{:04}".format(imgval))
        imgval += 1
    return imgval
    
class Gradient():
    def __init__(self, e1, rgb1, e2, rgb2):
        self.slope = np.subtract(rgb2, rgb1) / (e2 - e1)
        self.interc = np.subtract(rgb2, e2 * self.slope)
        
    def color(self, energy):
        return list(self.slope * energy + self.interc)
        
    
class AtomNode():
    def __init__(self, startIt, energyColorGradient):
        #Needed declarations
        self.it = startIt - 1
        self.colNum = 0
        self.pastAD = []
        self.currentAD = []
        
        self.advanceAtom(energyColorGradient)
               
    def advanceAtom(self, energyColorGradient):
        self.it += 1
        #Adding next atom and advancing
        myspace = {"pastAD": [], "currentAD": []}
        cmd.select("past", "id " + str(self.it))
        cmd.iterate("past", "pastAD += [resi, q, b]", space=myspace)
        cmd.iterate_state(-1, "past", "pastAD += [x, y, z]", space=myspace)
        cmd.select("current", "id " + str(self.it + 1))
        cmd.iterate("current", "currentAD += [resi, q, b]", space=myspace)
        cmd.iterate_state(-1, "current", "currentAD += [x, y, z]", space=myspace)
        
        self.pastAD = myspace["pastAD"]
        self.currentAD = myspace["currentAD"]
        
        #Showing past atom
        cmd.show("spheres", "past")
        cmd.set_color("color" + str(self.colNum), energyColorGradient.color(self.pastAD[2]))
        cmd.color("color" + str(self.colNum), "past")
        self.colNum += 1
        
    def addNext(self, energyColorGradient):
        returnVal = None
        if self.currentAD == []:
            return "END"
        elif self.pastAD[0] != self.currentAD[0]:
            self.advanceAtom(energyColorGradient)
            return "newProt"
        elif self.pastAD[1] !=  self.currentAD[1]:
            self.advanceAtom(energyColorGradient)
            return "newBranch"
        else:
            x1, y1, z1 = self.pastAD[3], self.pastAD[4], self.pastAD[5]
            rgb1 = energyColorGradient.color(self.pastAD[2])             
            x2, y2, z2 = self.currentAD[3], self.currentAD[4], self.currentAD[5]
            rgb2 = energyColorGradient.color(self.currentAD[2])
            #Making connector
            radius = .5
            holdview = cmd.get_view()
            cmd.load_cgo([9.0, x1, y1, z1, x2, y2, z2, radius, rgb1[0], rgb1[1], rgb1[2], rgb2[0], rgb2[1], rgb2[2]], "c" + str(self.it))    
            cmd.set_view(holdview)
            self.advanceAtom(energyColorGradient)
            return None
        
        
     
        
def main():
    #Load in structure
    cmd.load("..\\testfiles\\trackplot.pdb")
    cmd.hide("everything")
    cmd.bg_color("white")
    
    #Sets up energy color gradient (this can be changed at any time)
    energyColor = Gradient(86.56, (1, 0 , 0), 0, (.70, .70, .70))
    
    cmd.set_view((\
        0.957694948,    0.028183326,   -0.286400944,\
        -0.138612419,    0.917322457,   -0.373235583,\
        0.252203047,    0.397143424,    0.882423639,\
        -0.000240028,   -0.004241563, -948.997741699,\
        -249.833526611, -132.563354492, 2972.489990234,\
        -40.742408752, 1938.937011719,  -20.000000000 ))
    
    atomInfo = AtomNode(1, energyColor)
    breaker = None
    while(breaker == None):
        breaker = atomInfo.addNext(energyColor)
    print(breaker)
    breaker = None
    while(breaker == None):
        breaker = atomInfo.addNext(energyColor)
    print(breaker)
    atomInfo.advanceAtom(energyColor)
    breaker = None
    while(breaker == None):
        breaker = atomInfo.addNext(energyColor)
    print(breaker)
        
        
    
    '''
    
    
    imgval = scantoview((\
        0.301146746,    0.683516979,   -0.664916396,\
        0.748767972,   -0.601271391,   -0.278967589,\
        -0.590473652,   -0.413857192,   -0.692866027,\
        -0.000024408,    0.000082886, -268.568328857,\
        111.845207214,  103.542381287,  915.251525879,\
        141.994079590,  395.165588379,  -20.000000000 ), 20, imgval)
    '''
    
    

    
    
if __name__ == "__main__":
    main()