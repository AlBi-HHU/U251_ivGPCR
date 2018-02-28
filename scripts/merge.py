#!/usr/bin/python
import sys
import subprocess

if len(sys.argv) not in [5]:
  #sys.stderr.write("Usage: " + sys.argv[0] + " <INPUT_FILE_PVAL> <INPUT_FILE_SCORE> <INPUT_FILE_MODULE> <INPUT_FILE_MGI> <INPUT_FILE_KEGG_MAPPING> <INPUT_GO_FILE> <INPUT_GO_ANCESTOR>\n")
  sys.stderr.write("Usage: " + sys.argv[0] + " <INPUT_FILE_PVAL> <INPUT_FILE_MODULE> <INPUT_GO_FILE> <INPUT_GO_ANCESTOR>\n")
  sys.exit(1)

# read module file
module = set()

f_module = open(sys.argv[2])
for line in f_module:
  if line.startswith("#"): continue
  s = line.split()
  if s[1] != "NaN": module.add(s[0])

# # read mgi file
# mgi = {}
#
# f_mgi = open(sys.argv[4])
# for line in f_mgi:
#   s = line.rstrip("\n").split("\t")
#   mgi[s[0]] = [s[1], s[2]]
#
# # read kegg mapping file
# kegg = {}
#
# f_kegg = open(sys.argv[5])
# for line in f_kegg:
#   s = line.split()
#   if s[1] in kegg:
#     kegg[s[1]] += " "
#   else:
#     kegg[s[1]] = ""
#
#   kegg[s[1]] += s[0]

# read go file
go = {}

f_go = open(sys.argv[3])
for line in f_go:
  s = line.rstrip("\n").split("\t")
  if len(s) > 3:
      geneid = s[0]
      domain = s[3]
      term = s[1]
  if term == "": continue

  if geneid not in go:
    go[geneid] = [set(), set(), set()]
  if domain == "biological_process":
    go[geneid][0].add(term)
  elif domain == "molecular_function":
    go[geneid][1].add(term)
  elif domain == "cellular_component":
    go[geneid][2].add(term)

# add all the way up the root
f_root = open(sys.argv[4])
go_root = {}
for line in f_root:
  s = line.rstrip("\n").split("\t")
  go_root[s[0]] = set(s[2:-1])

for ensg in go:
  bp = set()
  for go_term in go[ensg][0]:
    bp |= go_root[go_term]
  go[ensg][0] |= bp

  mf = set()
  for go_term in go[ensg][1]:
    mf |= go_root[go_term]
  go[ensg][1] |= mf

  cc = set()
  for go_term in go[ensg][2]:
    cc |= go_root[go_term]
  go[ensg][2] |= cc

# # extract score from input_file_score
output = {}
f_score = open(sys.argv[1])
for line in f_score:
  if line.startswith("#"):
    continue
  s = line.split()
  if s[0] in module:
    output[s[0]] = ["True", s[1]]
  else:
    output[s[0]] = ["False", s[1]]

# extract logFC from input_file_pval
f_pval = open(sys.argv[1])
for line in f_pval:
  if line.startswith("#"):
    continue
  s = line.split()
  if len(s) > 2: output[s[0]] += [s[2]]
  else: output[s[0]] += " "

# # add mgi_symbol and mgi_id and kegg_id to output
# for key in output:
#   output[key] += mgi[key]
#   output[key] += [kegg[key]]

# kegg_pathway = {}
# f_pathway = open(sys.argv[8])
# for line in f_pathway:
#   s = line.rstrip("\n").split("\t")
#   kegg_pathway[s[0]] = s[1:]

# print output
print("ID\tModule\tP-value\tFC\tGO Process\tGO Function\tGO Component\tLabel")
for key in output:
  if key in go:
    go_bp = "|".join(go[key][0])
    go_mf = "|".join(go[key][1])
    go_cc = "|".join(go[key][2])
  else:
    go_bp = ""
    go_mf = ""
    go_cc = ""
  # if key in kegg_pathway:
  #   p = " ".join(kegg_pathway[key])
  # else:
  #   p = ""

  print("\t".join([key] + output[key] + [go_bp, go_mf, go_cc] + [key]))
