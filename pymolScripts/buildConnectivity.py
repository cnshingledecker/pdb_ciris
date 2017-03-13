from pymol import cmd

def buildConnectivity(past, current, it):
    x1, y1, z1 = past[3], past[4], past[5]
    r1, g1, b1 = past[2]/84.43, 0, 0
    x2, y2, z2 = current[3], current[4], current[5]
    r2, g2, b2 = current[2]/84.43, 0, 0
    radius = .4
    cmd.load_cgo([9.0, x1, y1, z1, x2, y2, z2, radius, r1, g1, b1, r2, g2, b2], "c" + str(it))
    

def main():
    cmd.load("..\\testfiles\\trackplotShort.pdb")
    cmd.hide("everything")
    
    it = 4
    while(True):
        myspace = {"pastAD": [], "currentAD": []}
        cmd.select("past", "id " + str(it - 1))
        cmd.iterate("past", "pastAD += [resi, q, b]", space=myspace)
        cmd.iterate_state(-1, "past", "pastAD += [x, y, z]", space=myspace)
        print(myspace["pastAD"][3])
        cmd.select("current", "id " + str(it))
        cmd.iterate("current", "currentAD += [resi, q, b]", space=myspace)
        cmd.iterate_state(-1, "current", "currentAD += [x, y, z]", space=myspace)
        cmd.show("spheres", "past or current")
        if myspace["pastAD"][0] !=  myspace["currentAD"][0]:
            break
        elif myspace["pastAD"][1] !=  myspace["currentAD"][1]:
            break
        elif it == 622:
            break
        buildConnectivity(myspace["pastAD"], myspace["currentAD"], it)
        it += 1
    
if __name__ == "__main__":
    main()