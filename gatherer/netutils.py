import hashlib
import requests
import time
import json
from gatherer.util import *
from gatherer.bioinfo import *
from concurrent.futures import ThreadPoolExecutor, as_completed
from gatherer.rna_type import *
import sys
from tqdm import tqdm
import obonet
from urllib import request

def download_to(url, final_path):
    print("Downloading", url, "to", final_path)
    try:
        request.urlretrieve(url, final_path)
    except Exception as err:
        print(err)
        return False
    print("Download with success")
    return True

def get_md5(sequence):
    """
    Calculate md5 for an RNA sequence
    """
    # RNAcentral stores DNA md5 hashes
    sequence = sequence.replace('U','T')
    # get the md5 digest
    m = hashlib.md5()
    m.update(sequence.encode('utf-8'))
    return m.hexdigest()

def get_rnacentral_id_by(argument, value, rna_central_api):
    """
    Parse json output and return the RNAcentral id.
    """
    url = rna_central_api
    r = requests.get(url, params = {argument: value})
    data = r.json()
    if data['count'] > 0:
        #print("Got response: \n\t" + str(data['results'][0]))
        #print(value + " matchs an ID in RNACentral: " + str(data['results'][0]['rnacentral_id']))
        return data['results'][0]['rnacentral_id']
    else:
        #print(value + "Does not match an real ID in RNACentral")
        return None

def get_rnacentral_json(value, rna_central_api):
    url = rna_central_api+'/'+value
    r = requests.get(url, params = {"format": "json"})
    #print(str(r.json()))
    data = r.json()
    if 'rnacentral_id' in data:
        return data
    else:
        return None

def confirm_rnacentral_id(value, rna_central_api):
    """
    Parse json output and return the RNAcentral id.
    """
    url = rna_central_api+'/'+value.split("_")[0]
    r = requests.get(url, params = {"format": "json"})
    #print(str(r.json()))
    data = r.json()
    if 'rnacentral_id' in data:
        #print(value + " matchs an real ID in RNACentral")
        return data['rnacentral_id']
    else:
        print(value + " does not match an real ID in RNACentral")
        return None

def retrieve_rnacentral_id(seq_id, seq, rna_central_api):
    #print("Trying to retrieve " + seq_id)
    result = None
    try:
        if "URS" in seq_id:
            result = confirm_rnacentral_id(seq_id, rna_central_api)
        if result == None:
            result = get_rnacentral_id_by("external_id",seq_id, rna_central_api)
        if result == None and seq != "":
            result = get_rnacentral_id_by("md5",get_md5(seq), rna_central_api)
    except Exception:
        print("Could not retrieve information for " + seq_id)
    return result, seq_id, seq

def retrieve_quickgo_annotations(chunk, api_url, taxon_id):
    result = []
    gene_lines = []
    gene_ids = ",".join(chunk)
    for i in range(len(chunk)):
        if "URS" in chunk[i]:
            if not "_"+taxon_id in chunk[i]:
                chunk[i] = chunk[i]+"_"+taxon_id
    gene_ids = ",".join(chunk)
    for aspect in ["biological_process","molecular_function","cellular_component"]:
        requestURL = (api_url+"?selectedFields=geneProductId&selectedFields=goId&geneProductId="
                        +gene_ids+"&taxonId="+taxon_id+"&aspect="+aspect)
        tries = 0
        response_lines = []
        max_tries = 5
        while tries < max_tries:
            if tries > 0:
                #print("Trying again, waiting " + str(tries*tries))
                time.sleep(tries*tries)
            #print("Requesting:\n\t"+requestURL)
            try:
                response = requests.get(requestURL, headers={"Accept":"text/tsv"})
                if response.ok:
                    #print("Got okay response")
                    tries = max_tries
                    text = response.text
                    lines = text.split("\n")[1:]
                    response_lines = [line.split("\t")+[aspect] for line in lines]
                    tries = max_tries
                else:
                    #print("Response not okay.")
                    #print(response.text)
                    try:
                        json_response = json.loads(response.text)
                        msgs = json_response["messages"]
                        #print("Analyzing error msgs")
                        for msg in msgs:
                            if "The 'Gene Product ID' parameter contains in" in msg:
                                invalid_ids = msg.split(": ")[-1].split(", ")
                                #print("Invalid ids: "+ str(invalid_ids))
                                #invalid_ids.append(invalid_id)
                                for id_ in invalid_ids:
                                    #print("Replacing " + id_)
                                    requestURL = requestURL.replace(id_,"")
                                    requestURL = requestURL.replace(",,",",")
                                #print("New request URL:\n\t" + requestURL)
                        tries += 1
                    except exception:
                        tries += 1
            except:
                tries += 1
        gene_lines += response_lines
    added_line = False
    for line in gene_lines:
        valid = True
        valid = valid and len(line) > 1
        if valid:
            valid = valid and line[0] != ""
        if valid:
            result += [line[1:]]
            added_line = True
    '''if added_line:
        #print(str(len(result)) + " annotations retrieved")
    else:
        print("No annotations retrieved for " + gene_ids)'''
    return result

def get_gene_info(gene_name, confs, sequence):
    new_id, old_id, seq = retrieve_rnacentral_id(gene_name, sequence,
                                    confs["rna_central_api"])
    retrieval_id = new_id if new_id != None else gene_name
    description = None
    rna_type = None
    if "URS" in retrieval_id:
        short_id = retrieval_id.split("_")[0]
        rnacentral_json = get_rnacentral_json(short_id, confs["rna_central_api"])
        if "description" in rnacentral_json:
            description = rnacentral_json["description"]
        if "rna_type" in rnacentral_json:
            rna_type = rnacentral_json["rna_type"]
    
    return gene_name, (new_id, description, rna_type)

def parallel_rnacentral_requester(to_retrieve, seqs_dict, confs, tries = 3):
    info_by_id = {}
    while tries > 0:
        print("Trying to retrieve " + str(len(to_retrieve)) + " ids.")
        failed = 0
        to_retrieves = [list(chunk) for chunk in chunks(list(to_retrieve), 100)]
        for chunk in tqdm(to_retrieves):
            processes = []
            with ThreadPoolExecutor(max_workers=80) as executor:
                for seq_id in to_retrieve:
                    if seq_id in chunk:
                        if seq_id in seqs_dict:
                            processes.append(
                                executor.submit(get_gene_info, seq_id, confs, seqs_dict[seq_id]))
                        else:
                            processes.append(
                                executor.submit(get_gene_info, seq_id, confs, ""))

                for task in as_completed(processes):
                    if task.result() != None:
                        seq_id, info = task.result()
                        #print(result)
                        to_retrieve.remove(seq_id)
                        info_by_id[seq_id] = info
                    else:
                        failed += 1
        print("\t" + str(failed) + " failed.")
        tries -= 1
        if failed == 0:
            tries = 0
    return info_by_id

def update_with_info(annotation_path, output_path, confs, 
    sep_id_by_dot = True, seqs_dict = None):
    annotation = pd.read_csv(annotation_path, sep="\t", header=None, names=["seqname", 
                "source", "feature", "start", "end", "score", "strand", "frame", "attribute"])
    raw_ids = [get_gff_attributes(attr_str)["ID"] 
    for attr_str in annotation["attribute"].tolist()]
    to_retrieve = [".".join(_id.split(".")[:-1]) for _id in raw_ids] if sep_id_by_dot else raw_ids

    if seqs_dict == None:
        seqs_dict = {}

    info_by_id = parallel_rnacentral_requester(to_retrieve, seqs_dict, confs)

    retrieval_stats = {"total": 0, "rfams_attributed": 0,
            "descriptions_attributed": 0,
            "types_attributed": 0,
            "IDs_attributed": 0}
    
    def update_attribute(row, retrieval_stats):
        attr_str = row["attribute"]
        attributes = None
        try:
            attributes = get_gff_attributes(attr_str)
        except:
            print("Row could not have attributes parsed:\n\t"+str(row))
        if attributes != None:
            short_id = ".".join(attributes["ID"].split(".")[:-1])
            rfam_id = get_rfam_from_rnacentral(short_id.split("_")[0])
            if not "rfam" in attributes and rfam_id != None:
                if len(rfam_id) == 7:
                    attributes["rfam"] = rfam_id
                    retrieval_stats["rfams_attributed"] += 1
            if short_id in info_by_id:
                new_name, description, rna_type = info_by_id[short_id]
                if not "description" in attributes and description != None:
                    attributes["description"] = description
                    retrieval_stats["descriptions_attributed"] += 1
                if not "type" in attributes and rna_type != None:
                    attributes["type"] = rna_type
                    retrieval_stats["types_attributed"] += 1
                if new_name != None:
                    attributes["ID"] = attributes["ID"].replace(short_id, new_name)
                    retrieval_stats["IDs_attributed"] += 1
        retrieval_stats["total"] += 1
        return get_gff_attributes_str(attributes)
    annotation["attribute"] = annotation.apply(lambda row: update_attribute(row,
                                                                    retrieval_stats),
                                            axis=1)
    annotation.to_csv(output_path, sep="\t", index=False, header=False)

    return retrieval_stats

def retrieve_func_annotation(annotation_path, output, confs, taxon_id):
    annotation = pd.read_csv(annotation_path, sep="\t", header=None, 
                names=["seqname", "source", "feature", "start", "end", 
                        "score", "strand", "frame", "attribute"])

    to_retrieve = get_ids_from_annotation(annotation)
    id_chunks = chunks(list(to_retrieve), 75)
    results = set()
    annotations_retrieved = False

    for chunk in tqdm(list(id_chunks)):
        annotations = retrieve_quickgo_annotations(chunk, confs["quickgo_api"], taxon_id)
        for id_, go, aspect in annotations:
            results.add((id_, go, aspect))
            annotations_retrieved = True

    if annotations_retrieved:
        with open(output, 'w') as stream:
            for id_, go, aspect in results:
                stream.write(id_+"\t"+go+"\t"+aspect+"\n")

def get_gene_annotation(gene_name, base_url, taxid):
    if base_url[-1] != '/':
        base_url += "/"
    #url = base_url+gene_name.split("_")[0]+"/go-annotations/"+gene_name.split("_")[1]
    url = base_url+gene_name+"/go-annotations/"+str(taxid)
    r = requests.get(url, params = {"format": "json"})
    if r.status_code == 200:
        data = r.json()
        associations = []
        for result in data:
            associations.append((result["rna_id"], result["go_term_id"]))
        return gene_name, associations, r.status_code
    else:
        return gene_name, None, r.status_code

def go_from_rnacentral(id_list, api_url, taxid, tries = 3):
    associations = set()
    new_ids = list(set([x.split("_")[0].split(".")[0] for x in id_list]))
    while tries > 0:
        print("Trying to retrieve " + str(len(id_list)) + " ids.")
        failed = 0
        to_retrieves = [list(chunk) for chunk in chunks(list(new_ids), 100)]
        too_many = False
        for chunk in tqdm(to_retrieves):
            processes = []
            with ThreadPoolExecutor(max_workers=25) as executor:
                for seq_id in new_ids:
                    if seq_id in chunk:
                        processes.append(
                            executor.submit(get_gene_annotation, seq_id, api_url, taxid))

                for task in as_completed(processes):
                    gene_name, new_associations, response = task.result()
                    if new_associations != None:
                        if len(new_associations) > 0:
                            associations.update(new_associations)
                            #print(str(len(new_associations)) + " annotations for " + gene_id)
                        else:
                            print("No annotations for " + gene_name)
                        new_ids.remove(gene_name)
                    else:
                        print(gene_name + " failed with " + str(response))
                        if response == 429:
                            too_many = True
                        failed += 1
        
        print("\t" + str(failed) + " failed.")
        if too_many:
            time.sleep(3)
        else:
            tries -= 1
        if failed == 0:
            tries = 0
    return associations

def get_term_ontology(go_obo):
    learned = {}
    with open(go_obo, 'r') as stream:
        ids_to_assign = []
        namespace = ""
        for line in stream:
            if "[Term]" in line:
                if namespace != "":
                    for id_ in ids_to_assign:
                        learned[id_] = namespace
                        namespace = ""
                ids_to_assign = []
            elif "namespace: " in line:
                namespace = line.replace("namespace: ", "").rstrip("\n")
            elif "id: " in line:
                new_id = line.replace("id: ", "").rstrip("\n")
                ids_to_assign.append(new_id)
            elif "alt_id: " in line:
                new_id = line.replace("alt_id: ", "").rstrip("\n")
                ids_to_assign.append(new_id)
        if namespace != "":
            for id_ in ids_to_assign:
                learned[id_] = namespace
                namespace = ""
    return learned

def retrieve_func_annotation_rnacentral(annotation_path, output, confs, taxon_id):
    annotation = pd.read_csv(annotation_path, sep="\t", header=None, 
                names=["seqname", "source", "feature", "start", "end", 
                        "score", "strand", "frame", "attribute"])
    
    to_retrieve = get_ids_from_annotation(annotation)
    term_ontologies = get_term_ontology(confs["go_obo"])
    results = go_from_rnacentral(to_retrieve, confs['rna_central_api'], taxon_id)
    final_results = []
    for id_, go in results:
        if go in term_ontologies:
            final_results.append((id_,go,term_ontologies[go]))
        else:
            print(go, "from", id_, "is an unknown term.")
    '''id_chunks = chunks(list(to_retrieve), 75)
    results = set()
    annotations_retrieved = False

    for chunk in tqdm(list(id_chunks)):
        annotations = retrieve_quickgo_annotations(chunk, confs["rna_central_api"], taxon_id)
        for id_, go, aspect in annotations:
            results.add((id_, go, aspect))
            annotations_retrieved = True

    if annotations_retrieved:'''
    if len(final_results) > 0:
        with open(output, 'w') as stream:
            for id_, go, aspect in final_results:
                stream.write(id_+"\t"+go+"\t"+aspect+"\n")