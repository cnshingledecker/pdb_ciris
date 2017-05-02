from pymol import cmd
import numpy as np
        
#Stores a gradient function for use in coloring
class Gradient():
    def __init__(self, e1, rgb1, e2, rgb2):
        self.eMax = e1
        self.slope = np.subtract(rgb2, rgb1) / (e2 - e1)
        self.interc = np.subtract(rgb2, e2 * self.slope)
      
    #Given an energy value, gives a RGB value
    def color(self, energy):
        if energy > self.eMax:
            return [0, 0, 0]
        else:
            return list(self.slope * energy + self.interc)

#Class for storing image print settings
class Pngprinter():
    def __init__(self, actuallyPrint = "Do it"):
        self.filebase = "img"
        self.width = "6in"
        self.height = "6in"
        self.dpi = 150
        self.ray = 1
        self.quiet = 1
        self.imgval = 0
        self.actuallyPrint = actuallyPrint
        
    #Used to print an image
    def out(self):
        if self.actuallyPrint == "Do it":
            cmd.png("{}{:05d}".format(self.filebase, self.imgval), self.width, self.height, self.dpi, self.ray, self.quiet)
        elif self.actuallyPrint == "Quick":
            cmd.save("{}{:05d}.png".format(self.filebase, self.imgval))
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
    
    #Creates a selection of one atom at pdbID and outputs a List of important data from the selection
    def setdata(self, selectionName, pdbID):
        myspace = {"atomdata": []}
        cmd.select(selectionName, "id " + str(pdbID))
        cmd.iterate(selectionName, "atomdata += [ID, resi, q, b,]", space=myspace)
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
        cmd.set_color("color" + str(self.colNum), self.ECG.color(self.pastAD[3]))
        cmd.color("color" + str(self.colNum), "past")
        self.colNum += 1
        
    #Makes connector
    def connect(self):
        x1, y1, z1 = self.pastAD[4], self.pastAD[5], self.pastAD[6]
        rgb1 = self.ECG.color(self.pastAD[3])
        x2, y2, z2 = self.currentAD[4], self.currentAD[5], self.currentAD[6]
        rgb2 = self.ECG.color(self.currentAD[3])
        holdview = cmd.get_view()
        cmd.load_cgo([9.0, x1, y1, z1, x2, y2, z2, self.conRad, rgb1[0], rgb1[1], rgb1[2], rgb2[0], rgb2[1], rgb2[2]], "c" + str(self.it))    
        cmd.set_view(holdview)        
        
    #method which adds the next atom and connects it 
    def addNext(self):
        if self.currentAD == []:
            return "END"
        elif self.pastAD[1] != self.currentAD[1]:
            #Finding last proton
            myspace = {"near": []} 
            cmd.select("poss", "resn PRO")
            cmd.iterate("poss", "near.append([ID, resi, q, b])", space=myspace)
            near = myspace["near"]
            near.sort(key=lambda x: x[0])
            for i in range(0, len(near)):
                if near[i][0] == self.currentAD[0]:
                    i -= 1
                    break
            self.pastAD = self.setdata("past", near[i][0])
            #Advancing atom
            self.connect()
            self.advanceAtom()
            return "newProt"
        elif int(round(self.pastAD[2]*100)) != int(round(self.currentAD[2]*100)):
            #Finding highest energy "close" atom
            myspace = {"near": []}
            cmd.select("poss","current around 13")
            cmd.iterate("poss", "near.append([ID, resi, q, b])", space=myspace)
            near = myspace["near"]
            near.sort(key=lambda x: x[2])
            self.pastAD = self.setdata("past", near[0][0])
            #Advancing atom
            self.connect()
            self.advanceAtom()
            return "newBranch"
        else:
            self.connect()
            self.advanceAtom()
            return None
            
    #Adds the next whole branch
    def addNextBranch(self, pngprinter):
        breaker = None
        while(breaker == None):
            breaker = self.addNext()
            pngprinter.out()
        return breaker
        
    def addNextAtoms(self, pngprinter, start, stop):
        for i in range(start, stop):
            self.addNext()
            pngprinter.out()

#Linearly shifts view to a new point. Give it the new view matrix, the Pngprinter, the number of imgs to print and the number the img starts at.
def scantoview(newview, steps, pngprinter):
    oldview = cmd.get_view()
    delta = np.subtract(newview, oldview)
    delta = delta / steps
    for i in range(0, steps + 1):
        currentview = np.add(oldview, delta * i)
        cmd.set_view(currentview)
        pngprinter.out()
        
def main():
    #Load in structure and set stage
    cmd.load("..\\testfiles\\trackplotClean.pdb")
    cmd.hide("everything")
    cmd.bg_color("white")
    cmd.clip("slab", 1000)
    
    #Creates important objects: atom node and pngprinter
    atomInfo = AtomNode(Gradient(35, (1, 0 , 0), 0, (.70, .70, .70)))
    pngprinter = Pngprinter("Quick")
    
    #Setting up for quick prints
    #pngprinter.ray = 0
    #pngprinter.dpi = 100
    
    #Creates the movie
    cmd.set_view((\
        0.659889638,    0.127799481,    0.740413308,\
        -0.385417193,    0.903480828,    0.187554732,\
        -0.644981921,   -0.409131795,    0.645453334,\
        -0.006196454,    0.004652441, -788.275756836,\
        63.789550781,  -30.806030273, 2916.111083984,\
        -32534.138671875, 34107.878906250,  -20.000000000 ))
        
    atomInfo.addNextAtoms(pngprinter, 1, 117)
    
    scantoview((\
        0.659820437,    0.127603829,    0.740508616,\
        -0.385310411,    0.903510094,    0.187633067,\
        -0.645116448,   -0.409128159,    0.645321131,\
        -0.006378829,    0.004803181, -1291.816406250,\
        -1.017456055,    1.457458496, 2915.541503906,\
        -32029.906250000, 34612.128906250,  -20.000000000 ), 5, pngprinter)
    
    atomInfo.addNextAtoms(pngprinter, 117, 754)
    
    scantoview((\
        0.659625113,    0.135221750,    0.739315271,\
        -0.385568261,    0.905246913,    0.178425491,\
        -0.645161927,   -0.402749211,    0.649261534,\
        -0.008741498,    0.046658099, -1614.154052734,\
        188.652221680, -151.970001221, 2749.275390625,\
        -31419.833984375, 34638.425781250,  -20.000000000 ), 5, pngprinter)
    
    #cmd.hide("everything")
    atomInfo.addNextAtoms(pngprinter, 754, 1178)
    
    
    #print(atomInfo.addNextBranch(pngprinter))
    #pngprinter.actuallyPrint = "Do it"
    #pngprinter.out()    
    print(atomInfo.pastAD)
    print(atomInfo.currentAD)
    
    
    
        
    
if __name__ == "__main__":
    main()