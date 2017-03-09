import argparse

#Stores a wsv entry and can convert to an ATOM entry
class wsvEntry(object):
    def __init__(self, entry):
        parts = entry.split()
        self.ptype = parts[0]
        self.energy = float(parts[1])
        self.z = float(parts[2])
        self.y = float(parts[3])
        self.x = float(parts[4])
        self.pro_coll_num = int(parts[5])
        self.generation = int(parts[6])
    
    def __str__(self):
        return "{} {} {:8.3f} {:8.3f} {:8.3f} {} {}".format(self.ptype, self.energy, self.z, self.y, self.x, self.pro_coll_num, self.generation)
    
    #Outputs a ATOM entry based on a wsv entry. 
    def toAtom(self):
        outAtom = atom()
        outAtom.name = self.ptype
        outAtom.resName = self.ptype
        outAtom.resSeq = self.pro_coll_num      
        outAtom.x = self.x
        outAtom.y = self.y
        outAtom.z = self.z + 3000
        outAtom.occupancy = self.generation
        outAtom.tempFactor = self.energy/10
        if self.ptype == "ELE":
            outAtom.element = "O"
            outAtom.tempFactor = self.energy
        else:
            outAtom.element = "H"
            outAtom.tempFactor = self.energy/100
        return outAtom      

#Stores a wsv file and has methods to convert to a pdb file.
class wsvfile(object):
    def __init__(self, filename):
        self.entries = []
        with open(filename, "r") as fle:
            for line in fle:
                self.entries.append(wsvEntry(line))

    def __str__(self):
        outstring = ""
        for entry in self.entries:
            outstring += str(entry) + "\n"
        return outstring[:-1]
    
    #Outputs a pdbfile based on the wsv file
    def toPdb(self):
        outpdb = pdbfile()
        serialVal = 1
        for wsventry in self.entries:
            outpdb.entries.append(wsventry.toAtom())
            outpdb.entries[-1].serial = serialVal
            serialVal += 1
        return outpdb
        
#Stores all entries of an Atom entry and formats output
class atom(object):
    def __init__(self):
        self.serial = 0
        self.name = ""
        self.altLoc = ""
        self.resName = ""
        self.chainID = "X"
        self.resSeq = 0
        self.iCode = ""
        self.x = 0
        self.y = 0
        self.z = 0
        self.occupancy = 0
        self.tempFactor = 0
        self.element = ""
        self.charge = ""

    def __str__(self):
        return "ATOM  {:>5} {:<3}{:1}{:<3} {:1}{:>4}{:1}   {:>8.3f}{:>8.3f}{:>8.3f}{:>6.2f}{:>6.2f}          {:>2}{:>2}".format(self.serial, self.name, self.altLoc, self.resName, self.chainID, self.resSeq, self.iCode, self.x, self.y, self.z, self.occupancy, self.tempFactor, self.element, self.charge)

#Stores all entries of a pdbfile and formats output
class pdbfile(object):
    def __init__(self):
        self.entries = []
        
    def __str__(self):
        outstring = ""
        for entry in self.entries:
            outstring += str(entry) + "\n"
        return outstring[:-1]
        
def main():
    #Argument parsing
    parser = argparse.ArgumentParser(description='Convert wsv file to pdb file')
    parser.add_argument("wsvfile", help='wsv file location')
    parser.add_argument("-o", "--outfile", help='output pdb file location')
    args = parser.parse_args()
    
    #Control flow
    outpdb = wsvfile(args.wsvfile).toPdb()
    if args.outfile:
        with open(args.outfile, "w") as outgoing:
            outgoing.write(str(outpdb))
    else:
        print(outpdb)
    
if __name__ == "__main__":
    main()