import sys
from util import runCommand
from rna_type import *
import os
import json

rfam_path = "Rfam.seed"
if not os.path.exists(rfam_path):
    seed_ftp = "ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.seed.gz"
    code = run_command("wget " + seed_ftp)
    if code != 0:
        print("Could not download Rfam family informations.")
        quit()

    code = run_command("gzip -d Rfam.seed.gz")
    if code != 0:
        print("Could not extract Rfam family informations from gzip")
        quit()

rfam2rnacentral_path = "rfam.tsv"
if not os.path.exists(rfam2rnacentral_path):
    rfam2rnacentral_ftp = ("ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/"
                        +"current_release/id_mapping/database_mappings/rfam.tsv")
    code = run_command("wget " + rfam2rnacentral_ftp)
    if code != 0:
        print("Could not download rnacentral rfam mappings.")
        quit()

output_dir = sys.argv[1]

lines_writen = set()
family_type = dict()

rnacentral_rfam2type = {}
print("Making rnacentral2rfam from " + rfam2rnacentral_path)
with open(rfam2rnacentral_path,'r') as input_stream:
    with open(output_dir+"/rnacentral2rfam.tsv", 'w') as output_stream:
        for line in input_stream.readlines():
            cells = line.rstrip("\n").split()
            rnacentral = cells[0]
            rfam = cells[2]
            tp = cells[4]
            rnacentral_rfam2type[rfam] = tp
            new_line = rnacentral[3:]+"\t"+rfam[2:]+"\n"
            if not new_line in lines_writen:
                output_stream.write(rnacentral[3:]+"\t"+rfam[2:]+"\n")
                lines_writen.add(new_line)

print("Making rfam2type from " + rfam_path)
rfam_types = []
rfam2type = {}
with open(rfam_path, 'r', encoding='Windows-1252') as input_stream:
    current_rfam = ""
    ACs = []
    TPs = []
    looking_for_type = False
    for raw_line in input_stream.readlines():
        line = raw_line.strip()
        if looking_for_type:
            if "#=GF TP" in line:
                TPs.append(line.replace(" ", "").split("#=GFTP")[-1])
                looking_for_type = False
        elif "#=GF AC" in line:
            ACs.append(line.replace(" ", "").split("#=GFAC")[-1])
            looking_for_type = True
    for i in range(len(ACs)):
        rfam_id = ACs[i]
        full_type = TPs[i].rstrip(";").lstrip(";").split(";")
        rfam2type[rfam_id] = full_type
        rfam_types.append((rfam_id, ";".join(full_type)))

with open(output_dir+"/rfam2type.tsv",'w') as output_stream:
    for rfam_id, full_type in rfam_types:
        output_stream.write(rfam_id+"\t"+full_type+"\n")

print("Creating ncRNA Type Tree")

nodes = {}
nodes["scRNA"] = make_node([], parent="Gene")
for rfam_id, type_vec in rfam2type.items():
    for i in range(len(type_vec)):
        type_name = type_vec[i]
        if not type_name in nodes:
            nodes[type_name] = make_node([], parent=None)
        if i > 0:
            nodes[type_name]["parent"] = type_vec[i-1]
        if i < len(type_vec)-1:
            nodes[type_name]["children"].add(type_vec[i+1])

print("Showing type tree:")
for node_name, data in nodes.items():
    if data["parent"] == None:
        print_tree(node_name, nodes)

parent_sets = {}

for rfam_id, type_name in rnacentral_rfam2type.items():
    if type_name in aliases:
        type_name = aliases[type_name]
    if not type_name in nodes:
        if not type_name in parent_sets:
            parent_sets[type_name] = set()
        if rfam_id in rfam2type:
            real_rfam_type = rfam2type[rfam_id][-1]
            #parent_type = nodes[real_rfam_type]["parent"]
            parent_sets[type_name].add((real_rfam_type))

def choose_parent(parent_set):
    parent_list = list(parent_set)
    parent_list.sort(key = lambda x: node_height(x, nodes))
    return parent_list[0]

best_parent = {type_name: choose_parent(parents) 
                for type_name, parents in parent_sets.items()}

best_parent["other"] = None

for new_type in parent_sets.keys():
    print("New type: " + new_type 
        + ". Possible parents:\n\t"+str(parent_sets[new_type])
        +".\nBest parent: " + str(best_parent[new_type]))

for new_type, parent in best_parent.items():
    nodes[new_type] = make_node([], parent=parent)
    if parent != None:
        nodes[parent]["children"].add(new_type)

print("\nShowing type tree (expanded):")
for node_name, data in nodes.items():
    if data["parent"] == None:
        print_tree(node_name, nodes)

print("Saving tree")
with open(output_dir+"/ncRNA_type-tree.json", 'w') as fp:
    for node_name in nodes.keys():
        nodes[node_name]["children"] = list(nodes[node_name]["children"])
    json.dump(nodes, fp, indent=4)