from pymol import cmd
import numpy as np
        
#Linearly shifts view to a new point. Give it the new view matrix, the Pngprinter, the number of imgs to print and the number the img starts at.
def scantoview(newview, steps, pngprinter):
    oldview = cmd.get_view()
    delta = np.subtract(newview, oldview)
    delta = delta / steps
    for i in range(1, steps):
        currentview = np.add(oldview, delta * i)
        cmd.set_view(currentview)
        pngprinter.out()
    
#Stores a gradient function for use in coloring
class Gradient():
    def __init__(self, e1, rgb1, e2, rgb2):
        self.slope = np.subtract(rgb2, rgb1) / (e2 - e1)
        self.interc = np.subtract(rgb2, e2 * self.slope)
      
    #Given an energy value, gives a RGB value
    def color(self, energy):
        return list(self.slope * energy + self.interc)

#Class for storing image print settings
class Pngprinter():
    def __init__(self):
        self.filebase = "img"
        self.width = "6in"
        self.height = "6in"
        self.dpi = 150
        self.ray = 1
        self.quiet = 1
        self.imgval = 0
        
    #Used to print an image
    def out(self):
        cmd.png("{}{:05d}".format(self.filebase, self.imgval), self.width, self.height, self.dpi, self.ray, self.quiet)
        self.imgval += 1
  
#Device for traversing electron path  
class AtomNode():
    def __init__(self, energyColorGradient):
        #Needed declarations
        self.it = 0
        self.colNum = 0
        self.conRad = 0.5
        self.pastAD = []
        self.currentAD = []
        self.ECG = energyColorGradient
        self.advanceAtom()
        
    def setdata(self, selectionName, pdbID):
        myspace = {"atomdata": []}
        cmd.select(selectionName, "id " + str(pdbID))
        cmd.iterate(selectionName, "atomdata += [resi, q, b,]", space=myspace)
        cmd.iterate_state(-1, selectionName, "atomdata += [x, y, z]", space=myspace)
        return myspace["atomdata"]
       
    #stores the new atom data and shows the new atom
    def advanceAtom(self):
        self.it += 1
        #Adding next atom and advancing
        self.pastAD = self.setdata("past", self.it)
        self.currentAD = self.setdata("current", self.it + 1)
        
        #Showing past atom
        cmd.show("spheres", "past")
        cmd.set_color("color" + str(self.colNum), self.ECG.color(self.pastAD[2]))
        cmd.color("color" + str(self.colNum), "past")
        self.colNum += 1
        
    #Makes connector
    def connect(self):
        print(self.it)
        x1, y1, z1 = self.pastAD[3], self.pastAD[4], self.pastAD[5]
        rgb1 = self.ECG.color(self.pastAD[2])             
        x2, y2, z2 = self.currentAD[3], self.currentAD[4], self.currentAD[5]
        rgb2 = self.ECG.color(self.currentAD[2])
        holdview = cmd.get_view()
        cmd.load_cgo([9.0, x1, y1, z1, x2, y2, z2, self.conRad, rgb1[0], rgb1[1], rgb1[2], rgb2[0], rgb2[1], rgb2[2]], "c" + str(self.it))    
        cmd.set_view(holdview)
        
    #method which adds the next atom and connects it 
    def addNext(self):
        if self.currentAD == []:
            return "END"
        elif self.pastAD[0] != self.currentAD[0]:            
            self.advanceAtom()
            return "newProt"
        elif int(round(self.pastAD[1]*100)) != int(round(self.currentAD[1]*100)):
            myspace = {"near": []}            
            cmd.select("poss","current around 13")
            cmd.iterate("poss", "near.append([ID, resi, q, b])", space=myspace)
            near = myspace["near"]
            near.sort(key=lambda x: x[2])
            self.pastAD = self.setdata("past", near[0][0])
            self.connect()            
            self.advanceAtom()
            return "newBranch"
        else:
            self.connect()
            self.advanceAtom()
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
    
    atomInfo = AtomNode(energyColor)
    breaker = None
    while(breaker == None):
        breaker = atomInfo.addNext()
    print(breaker)
    breaker = None
    while(breaker == None):
        breaker = atomInfo.addNext()
    print(breaker)
    atomInfo.advanceAtom()
    breaker = None
    while(breaker == None):
        breaker = atomInfo.addNext()
    print(breaker)
    
    

    
    
if __name__ == "__main__":
    main()