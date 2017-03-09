import argparse

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

class wsvFile(object):
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
        
def main():
    parser = argparse.ArgumentParser(description='Convert wsv file to pdb file')
    parser.add_argument("wsvfile")
    args = parser.parse_args()
    print(wsvFile(args.wsvfile))
    
if __name__ == "__main__":
    main()