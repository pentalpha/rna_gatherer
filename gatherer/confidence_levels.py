from gatherer.util import get_subdirs, file_name
import os
import json
def load_confidence(intervals_file):
    #print("Reading " + intervals_file)
    with open(intervals_file, 'r') as stream:
        lines = [line.rstrip("\n").split("\t") for line in stream.readlines()]
        th_lines = lines[2:]
        confidences = [{} for i in range(len(th_lines[0])-1)]
        for cells in th_lines:
            metric = cells[0]
            for i in range(1,len(cells)):
                if cells[i] == "None":
                    confidences[i-1][metric] = None
                else:
                    confidences[i-1][metric] = float(cells[i])
        #for i in range(len(confidences)):
        #    print("Confidence "+ str(i)+": " + str(confidences[i]))
        return confidences
    return []

def get_available_species():
    confs_dir = (os.path.dirname(os.path.realpath(__file__)) 
                    + "/../data/confidence_levels/")
    species_dirs = get_subdirs(confs_dir)
    return [file_name(f) for f in species_dirs]

def get_species_dir(species):
    confs_dir = (os.path.dirname(os.path.realpath(__file__)) 
                    + "/../data/confidence_levels/")
    return confs_dir+"/"+species

def load_confidence_levels(species):
    species_dir = get_species_dir(species)
    confs = {"MF": load_confidence(species_dir+"/confidence_intervals-MF.tsv"),
            "BP": load_confidence(species_dir+"/confidence_intervals-BP.tsv"),
            "CC": load_confidence(species_dir+"/confidence_intervals-CC.tsv")}
    return confs

'''def load_metrics(table_path):
    metric_combinations = {}
    with open(table_path,'r') as stream:
        for raw_line in stream.readlines():
            cells = raw_line.rstrip("\n").split("\t")
            metric_combinations[cells[0]] = {
                "biological_process": cells[1].split("-"),
                "molecular_function": cells[2].split("-"),
                "cellular_component": cells[3].split("-")
            }
    return metric_combinations'''

def get_best_metrics(species):
    short_names = ["MF","BP","CC"]
    full_names = ["molecular_function", 
        "biological_process", 
        "cellular_component"]
    species_dir = get_species_dir(species)
    config_files = [species_dir+"/"+ont+"-config.tsv"
                    for ont in short_names]
    metric_combinations = {}
    for i in range(len(full_names)):
        ont_name = full_names[i]
        with open(config_files[i], 'r') as stream:
            for line in stream.readlines():
                cells = line.rstrip("\n").split("\t")
                if cells[0] != "confidenceLevel":
                    conf_level = cells[0]
                    if not conf_level in metric_combinations:
                        metric_combinations[conf_level] = {}
                    metrics_list = cells[1].split("-")
                    metrics_list.sort()
                    metric_combinations[conf_level][ont_name] = metrics_list
    return metric_combinations

def get_preset_confs(species):
    optimal_path = get_species_dir(species) + "/optimal.json"
    optimal_str = open(optimal_path, 'r').read()
    optimal_obj = json.loads(optimal_str)
    return optimal_obj


def geometric_filter(value, th):
    if value <= th:
        return value
    else:
        return None

def normal_filter(value, th):
    if value >= th:
        return value
    else:
        return None

def geometric_pass(value, th):
    return value <= th

def normal_pass(value, th):
    return value >= th

def compare_to_th(value, th, metric):
    if metric != "SOB" and metric != "FSH":
        return normal_pass(value, th)
    else:
        return geometric_pass(value, th)
