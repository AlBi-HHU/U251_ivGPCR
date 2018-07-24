import sys
import re
import csv
import itertools

# Input: All "go_and_kegg_annotation" files you want to merge to one List with FDRs vs Scores.
# Output: .csv file FDRs with all their Scores.

def main():
    #check how many inputfiles and set name
    filecounter = len(sys.argv)-1
    outputfilename = re.search("([US28UL3]+)_vs_[a-z]+_([A-Z]+)",sys.argv[1])
    outputfilename= "GO_"+outputfilename.group(1)+outputfilename.group(2)+".csv"
    #open all files one by one
    scorelist=[]
    datalist=[]
    for i in range (1,filecounter+1):
        csv_columnname = re.search("[US28UL3]+_vs_[a-z]+_[A-Z]+_([0-9.e-]+)",sys.argv[i])
        data = open(sys.argv[i])
        #read in GO line
        for line in data:
            GO_entry_score = re.search(r"GO:[0-9]+\s+[\w, -/]+\s([\d//.e-]+)",line)
            if GO_entry_score:
                scorelist.append(csv_columnname.group(1)+","+GO_entry_score.group(1))
                #datalist.append(csv_columnname.group(1))
    #datalist= list(itertools.zip_longest(*datalist))
    print(scorelist)
    file = open(outputfilename, 'w')
    file.write("FDR,Score\n")
    for i in range(0,len(scorelist)):
        file.write(scorelist[i]+"\n")

if __name__ == '__main__':
    main()
