from pymol import cmd
import numpy as np
        

def scantoview(newview, steps, imgstart):
    imgval = imgstart
    oldview = cmd.get_view()
    delta = np.subtract(newview, oldview)
    delta = delta / steps
    for i in range(1,steps):
        print(imgval)
        currentview = np.add(oldview, delta * i)
        cmd.set_view(currentview)
        cmd.png("img{:04}".format(imgval))
        imgval += 1
    return imgval
    
class branch():
    def __init__(self, startIt):
        #Needed declarations
        self.it = startIt - 1
        self.startE = -1
        self.rm = -1
        self.rb = -1
        self.gm = -1
        self.gb = -1
        self.bm = -1
        self.bb = -1
        self.pastAD = []
        self.currentAD = []
        
        self.advanceAtom()
               
    def advanceAtom(self):
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
        cmd.color("gray70", "past")
        
    def addNext(self):
        print(self.it)
        if self.currentAD == []:
            print("END")
            pass
        elif self.pastAD[0] !=  self.currentAD[0]:
            print("new prot")
            pass
        elif self.pastAD[1] !=  self.currentAD[1]:
            print("new Branch")
            pass
        else:
            x1, y1, z1 = self.pastAD[3], self.pastAD[4], self.pastAD[5]
            r1, g1, b1 = 70, 70, 70
            x2, y2, z2 = self.currentAD[3], self.currentAD[4], self.currentAD[5]
            r2, g2, b2 = 70, 70, 70
            #Making connector
            radius = .4
            holdview = cmd.get_view()
            cmd.load_cgo([9.0, x1, y1, z1, x2, y2, z2, radius, r1, g1, b1, r2, g2, b2], "c" + str(self.it))    
            cmd.set_view(holdview)
            
        self.advanceAtom()
        
        
        '''
        while(True):
            
            
            #Getting starting energy and colors
            if self.startE < 0:
                self.startE = myspace["pastAD"][2]
                self.rm = (.7 - 1) / (0 - self.startE)
                self.rb = .7
                self.gm = (.7 - 0) / (0 - self.startE)
                self.gb = .7
                self.bm = (.7 - 0) / (0 - self.startE)
                self.bb = .7
            
            #Coloring past atom and storing data.
            past = myspace["pastAD"]            
            x1, y1, z1 = past[3], past[4], past[5]
            r1, g1, b1 = past[2] * self.rm + self.rb, past[2] * self.gm + self.gb, past[2] * self.bm + self.bb
            cmd.set_color("chold"+str(self.it), [r1, g1, b1])
            
            #printing img
            cmd.png("img{:04d}".format(self.imgval))
            self.imgval += 1
            
            #Setting up current atom
            if myspace["currentAD"] == []:
                break
            elif myspace["pastAD"][0] !=  myspace["currentAD"][0]:
                break
            elif myspace["pastAD"][1] !=  myspace["currentAD"][1]:
                break
        
            current = myspace["currentAD"]
            x2, y2, z2 = current[3], current[4], current[5]
            r2, g2, b2 = current[2] * self.rm + self.rb, current[2] * self.gm + self.gb, current[2] * self.bm + self.bb
        
            #Making connector
            radius = .4
            holdview = cmd.get_view()
            cmd.load_cgo([9.0, x1, y1, z1, x2, y2, z2, radius, r1, g1, b1, r2, g2, b2], "c" + str(self.it))    
            cmd.set_view(holdview)
            self.it += 1
        '''
        
def main():
    cmd.load("..\\testfiles\\trackplot.pdb")
    cmd.hide("everything")
    
    cmd.set_view((\
        0.989454925,    0.119257055,    0.082188576,\
        -0.078681566,    0.919020355,   -0.386279583,\
        -0.121599913,    0.375739068,    0.918713212,\
        -0.000024408,    0.000082886, -1595.923706055,\
        111.845207214,  103.542381287,  915.251525879,\
        1469.349487305, 1722.521118164,  -20.000000000 ))
    
    branch1 = branch(1)
    for i in range(1,30):
        branch1.addNext()
        cmd.png("img{:04d}".format(i))
    
    '''
    it, imgval = 1, 1
    #Making the movie
    cmd.set_view((\
        0.989454925,    0.119257055,    0.082188576,\
        -0.078681566,    0.919020355,   -0.386279583,\
        -0.121599913,    0.375739068,    0.918713212,\
        -0.000024408,    0.000082886, -1595.923706055,\
        111.845207214,  103.542381287,  915.251525879,\
        1469.349487305, 1722.521118164,  -20.000000000 ))
    
    holdbranch = branch(it, imgval)
    it = holdbranch.it + 1
    imgval = holdbranch.imgval
    holdbranch = branch(it, imgval)
    it = holdbranch.it + 1
    imgval = holdbranch.imgval
    holdbranch = branch(it, imgval)
    it = holdbranch.it + 1
    imgval = holdbranch.imgval
    
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