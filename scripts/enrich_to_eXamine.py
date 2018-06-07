import sys
import re

def main():
    data = open(sys.argv[1])
    print("ID\tP-value\tLabel")
    for line in data:
        GO_entry = re.search(r"(GO:[0-9]+)\s+\S+\s+[a-z]\s+([\w, -/]+)\s+[\d//.]+\s+[\d//.]+\s+([\d//.]+)",line)
        if re.search("GO:[0-9]+",line):
            print(GO_entry.group(1)+"\t",end="")
            print(GO_entry.group(3)+"\t",end="")
            print(GO_entry.group(2)+"\t")
            
if __name__ == '__main__':
    main()
