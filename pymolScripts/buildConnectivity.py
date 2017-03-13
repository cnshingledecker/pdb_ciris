from pymol import cmd

def buildConnectivity():
    check = True
    cmd.select("current", "id 1")
    it = 1
    myspace = {"atomdata": []}
    while(check):
        it += 1
        cmd.select("current","id " + str(it))
        cmd.show("spheres", "current")
        cmd.iterate("current", "atomdata.append([ID, rank, index, elem])", space = myspace)
        cmd.iterate_state(-1, "current", "atomdata.append([x, y, z])", space = myspace)
        print(myspace["atomdata"])
        myspace["atomdata"] = []
        if it >= 10:
            check = False

def main():
    cmd.load("..\\testfiles\\trackplotShort.pdb")
    cmd.hide("everything")
    buildConnectivity()
    
if __name__ == "__main__":
    main()