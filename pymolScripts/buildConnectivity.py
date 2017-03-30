from pymol import cmd

class progroup():
    def __init__(self, startIt):
        self.it = startIt
        

class branch():
    def __init__(self, startIt):
        #Needed declarations
        self.it = startIt
        self.startE = -1
        self.rm = -1
        self.rb = -1
        self.gm = -1
        self.gb = -1
        self.bm = -1
        self.bb = -1
        
        while(True):
            
        
            #Setting up objects in pymol and collecting atom data
            myspace = {"pastAD": [], "currentAD": []}
            cmd.select("past", "id " + str(self.it))
            cmd.iterate("past", "pastAD += [resi, q, b]", space=myspace)
            cmd.iterate_state(-1, "past", "pastAD += [x, y, z]", space=myspace)
            cmd.select("current", "id " + str(self.it + 1))
            cmd.iterate("current", "currentAD += [resi, q, b]", space=myspace)
            cmd.iterate_state(-1, "current", "currentAD += [x, y, z]", space=myspace)
            cmd.show("spheres", "past or current")
            
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
            
            #Setting up current atom
            if myspace["currentAD"] == []:
                print("breakEnd")
                break
            elif myspace["pastAD"][0] !=  myspace["currentAD"][0]:
                print("breakPro")
                break
            elif myspace["pastAD"][1] !=  myspace["currentAD"][1]:
                print("breakEle")
                break
        
            current = myspace["currentAD"]
            x2, y2, z2 = current[3], current[4], current[5]
            r2, g2, b2 = current[2] * self.rm + self.rb, current[2] * self.gm + self.gb, current[2] * self.bm + self.bb
        
            #Making connector
            radius = .4
            cmd.load_cgo([9.0, x1, y1, z1, x2, y2, z2, radius, r1, g1, b1, r2, g2, b2], "c" + str(self.it))    
            self.it += 1
            
    def buildConnectivity(self, past, current):
        
        
        #Fixing periodic boundaries
        translate = False
        
        '''
        if x2 - x1 > 4:
            if x2 > 0:
                x2 += 1000
                cmd.translate([1000, 0, 0], 'current')
            elif x2 < 0:
                x2 += -1000
                cmd.translate([-1000, 0, 0], 'current')
        if y2 - y1 > 4:
            if y2 > 0:
                y2 += 1000
                cmd.translate([0, 1000, 0], 'current')
            if y2 < 0:
                y2 += -1000
                cmd.translate([0, -1000, 0], 'current')
        '''
        
        
    
    
    

def main():
    cmd.load("..\\testfiles\\trackplot.pdb")
    cmd.hide("everything")
    
    it = 1
    while it < 36367:
        holdbranch = branch(it)
        it = holdbranch.it
        it += 1
    
    
if __name__ == "__main__":
    main()