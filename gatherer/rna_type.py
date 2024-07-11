import gzip
import os
import json

rfam2type = {}
rfam2desc = {}
rnacentral2rfam = {}
rnacentral2go = {}
type_tree = {}

aliases = {"antisense_RNA": "antisense",
            "precursor_RNA": "miRNA"}

def make_node(children, parent = None):
    node = {"children": set(children),
            "parent": parent}
    return node

def print_tree(start_node, node_dict, height = 0):
    print("\t"*height + start_node)
    node = node_dict[start_node]
    for children in node["children"]:
        print_tree(children, node_dict, height + 1)

def node_height(node_name, node_dict):
    node = node_dict[node_name]
    if node["parent"] != None:
        return node_height(node["parent"], node_dict) + 1
    else:
        return 0

def load_rna_types():
    global_data = os.path.dirname(os.path.realpath(__file__)) + "/../data"
    with open(global_data+"/ncRNA_type-tree.json", 'r') as in_stream:
        json_str = in_stream.read()
        loaded_tree = json.loads(json_str)
        for key, value in loaded_tree.items():
            type_tree[key] = value

def get_type_list(rna_type):
    parent_type = type_tree[rna_type]["parent"]
    if parent_type != None:
        return get_type_list(parent_type) + [parent_type]
    else:
        return []

def get_full_type(rna_type):
    if rna_type in aliases:
        rna_type = aliases[rna_type]
    if len(type_tree.keys()) == 0:
        load_rna_types()

    if rna_type in type_tree:
        return ";".join(get_type_list(rna_type) + [rna_type])
    else:
        return ";".join([rna_type])

def load_rnacentral2rfam(load_to = {}):
    global_data = os.path.dirname(os.path.realpath(__file__)) + "/../data"
    with open(global_data+"/rnacentral2rfam.tsv",'r') as input_stream:
        for line in input_stream.readlines():
            cells = line.rstrip("\n").split()
            rnacentral = "URS"+cells[0]
            rfam = "RF"+cells[1]
            load_to[rnacentral] = rfam
    
    with gzip.open(global_data+"/rnacentral2rfam2.tsv.gz",'rt') as input_stream:
        for line in input_stream.readlines():
            cells = line.rstrip("\n").split()
            rnacentral = "URS"+cells[0]
            rfam = "RF"+cells[1]
            load_to[rnacentral] = rfam
    
    return load_to

def load_rnacentral2go():
    global_data = os.path.dirname(os.path.realpath(__file__)) + "/../data"
    with gzip.open(global_data+"/rnacentral2go.tsv.gz",'rt') as input_stream:
        for line in input_stream.readlines():
            cells = line.rstrip("\n").split()
            rnacentral = "URS"+(cells[0].split('_')[0])
            goid = "GO:"+cells[1]
            if not rnacentral in rnacentral2go:
                rnacentral2go[rnacentral] = set()
            rnacentral2go[rnacentral].add(goid)

def get_rfam_from_rnacentral(id_, print_response=False):
    if len(rnacentral2rfam.keys()) == 0:
        load_rnacentral2rfam(load_to=rnacentral2rfam)

    if id_ in rnacentral2rfam:
        return rnacentral2rfam[id_]
    elif id_.split("_")[0] in rnacentral2rfam:
        return rnacentral2rfam[id_.split("_"[0])]
    else:
        return None
    
def get_goid_from_rnacentral(id_, print_response=False):
    if len(rnacentral2go.keys()) == 0:
        load_rnacentral2go()

    if id_ in rnacentral2go:
        return list(rnacentral2go[id_])
    elif id_.split("_")[0] in rnacentral2go:
        return list(rnacentral2go[id_.split("_"[0])])
    else:
        return []

def load_rfam2type():
    print("Loading rfam2type")
    global_data = os.path.dirname(os.path.realpath(__file__)) + "/../data"
    with open(global_data+"/rfam2type.tsv",'r') as input_stream:
        for line in input_stream.readlines():
            cells = line.rstrip("\n").split()
            rfam = cells[0]
            tp = cells[1]
            rfam2type[rfam] = tp

def load_rfam2desc():
    print("Loading rfam2desc")
    global_data = os.path.dirname(os.path.realpath(__file__)) + "/../data"
    with gzip.open(global_data+"/rfam2desc.tsv.gz",'rt') as input_stream:
        for line in input_stream.readlines():
            cells = line.rstrip("\n").split()
            rfam = cells[0]
            desc = cells[1]
            rfam2desc[rfam] = desc

def get_rna_type(id_):
    if len(rfam2type.keys()) == 0:
        load_rfam2type()

    if id_ in aliases:
        id_ = aliases[id_]
    if id_ in rfam2type:
        return rfam2type[id_]
    else:
        return None
    
def get_rna_desc(id_):
    if len(rfam2desc.keys()) == 0:
        load_rfam2desc()

    if id_ in aliases:
        id_ = aliases[id_]
    if id_ in rfam2desc:
        return rfam2desc[id_]
    else:
        return None

if __name__ == '__main__':
    '''Testing the use of files from ../data'''
    print("Loading rnacentral2rfam")
    load_rnacentral2rfam()
    print("Loading rfam2type")
    load_rfam2type()
    print("Loading RNA Types")
    load_rna_types()

    print("CD-box full type: " + get_full_type("CD-box"))
    print("Intron full type: " + get_full_type("Intron"))
    print("lncRNA full type: " + get_full_type("lncRNA"))