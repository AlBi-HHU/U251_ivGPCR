import sys

def main():
    out = open(sys.argv[2],"w")
    genenames=[]
    with open(sys.argv[1]) as map:
        for line in map:
            linearr = line.split("\t")
            genenames.append(linearr[0])
    for gene in genenames:
        out.write("%s\n" % gene)
    out.close()

if __name__ == '__main__':
    main()
