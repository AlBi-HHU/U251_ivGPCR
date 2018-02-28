# python mappingTopGO.py GO_biomart.txt 3-group-CCS.txt all_go_uniq_ancestors.txt > 3-group-CCS_gsymb2go.map
import sys

f_root = open(sys.argv[3])
go_root = {}
for line in f_root:
  s = line.rstrip("\n").split("\t")
  go_root[s[0]] = set(s[2:-1])

go = {}
genes = set()
f = open(sys.argv[2])
for line in f:
  s = line.split()
  genes.add(s[0])

f = open(sys.argv[1])
for line in f:
  s = line.rstrip("\n").split("\t")
  if s[0] not in genes:
    continue

  if s[0] not in go:
    go[s[0]] = set()
  if s[3] != "":
    go[s[0]].add(s[1])
    go[s[0]] |= go_root[s[1]]

for key in go:
  print(key + "\t" + ", ".join(go[key]))
