import gzip
import hashlib
from bs4 import BeautifulSoup
import requests
import time
import json
import ftplib
from gatherer.final_steps import read_rfam2go
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
    try:
        data = json.loads(r.text)
        if data['count'] > 0:
            #print("Got response: \n\t" + str(data['results'][0]))
            #print(value + " matchs an ID in RNACentral: " + str(data['results'][0]['rnacentral_id']))
            return data['results'][0]['rnacentral_id']
        else:
            #print(value + "Does not match an real ID in RNACentral")
            return None
    except json.decoder.JSONDecodeError as err:
        print('Could not retrieve json from', url, r, r.text, file=sys.stderr)
        print(value, file=sys.stderr)
        print(err, file=sys.stderr)
        return None

def get_rnacentral_json(value, rna_central_api):
    url = rna_central_api+'/'+value
    r = requests.get(url, params = {"format": "json"})
    #print(str(r.json()))
    try:
        data = json.loads(r.text)
        if 'rnacentral_id' in data:
            return data
        else:
            return None
    except json.decoder.JSONDecodeError as err:
        print('Error parsing', value, url, r, r.text, file=sys.stderr)
        print(value, file=sys.stderr)
        print(err, file=sys.stderr)
        raise(err)

def confirm_rnacentral_id(value, rna_central_api):
    """
    Parse json output and return the RNAcentral id.
    """
    url = rna_central_api+'/'+value.split("_")[0]
    r = requests.get(url, params = {"format": "json"})
    #print(str(r.json()))
    try:
        data = json.loads(r.text)
        if 'rnacentral_id' in data:
            #print(value + " matchs an real ID in RNACentral")
            return data['rnacentral_id']
        else:
            print(value + " does not match an real ID in RNACentral", file=sys.stderr)
            return None
    except json.decoder.JSONDecodeError as err:
        print('Could not retrieve json from', url, r, r.text, file=sys.stderr)
        print(value, file=sys.stderr)
        print(err, file=sys.stderr)
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
    except Exception as err:
        print("Could not retrieve information for " + seq_id, file=sys.stderr)
        print(err, file=sys.stderr)
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
                    except Exception:
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

def load_rnacentral_details(rnacentral_ids):
    short_ids = [x.split('_')[0] for x in rnacentral_ids]
    global_data = os.path.dirname(os.path.dirname(os.path.realpath(__file__))) + "/data"
    input_path = global_data + '/rnacentral_details.tsv.gz'

    tps = {}
    descs = {}
    for rawline in gzip.open(input_path, 'rt'):
        cells = rawline.rstrip('\n').split('\t')
        if len(cells) == 3:
            short_id = cells[0].split('_')[0]
            if short_id in short_ids:
                tps[short_id] = cells[1]
                descs[short_id] = cells[2]
    
    return tps, descs

def load_all_rnacentral_details():
    global_data = os.path.dirname(os.path.dirname(os.path.realpath(__file__))) + "/data"
    input_path = global_data + '/rnacentral_details.tsv.gz'

    tps = {}
    descs = {}
    for rawline in gzip.open(input_path, 'rt'):
        cells = rawline.rstrip('\n').split('\t')
        if len(cells) == 3:
            short_id = cells[0].split('_')[0]
            tps[short_id] = cells[1]
            descs[short_id] = cells[2]
    
    return tps, descs

def update_with_info(annotation_path, output_path, 
        sep_id_by_dot = True):
    print('Reading', annotation_path, file=sys.stderr)
    annotation = pd.read_csv(annotation_path, sep="\t", header=None, names=["seqname", 
                "source", "feature", "start", "end", "score", "strand", "frame", "attribute"])
    
    info_fields1 = annotation["attribute"].tolist()

    print('Loading RNACentral tables', file=sys.stderr)
    all_rnacentral_rfam = load_rnacentral2rfam()
    all_rnacentral_rna_type, all_rnacentral_rna_desc = load_all_rnacentral_details()

    print('Getting info from tables', file=sys.stderr)
    info_fields2 = []
    retrieval_stats = {"total": 0, "rfams_attributed": 0,
        "descriptions_attributed": 0,
        "types_attributed": 0,
        "IDs_attributed": 0}
    for attr_str in tqdm(info_fields1):
        attributes = get_gff_attributes(attr_str)
        raw_id = attributes["ID"]
        first_part = ".".join(raw_id.split(".")[:-1]) if sep_id_by_dot else raw_id
        rna_id = first_part.split('_')[0]
        
        if not "rfam" in attributes and rna_id in all_rnacentral_rfam:
            attributes["rfam"] = all_rnacentral_rfam[rna_id]
            retrieval_stats["rfams_attributed"] += 1
        if not "description" in attributes and rna_id in all_rnacentral_rna_desc:
            attributes["description"] = all_rnacentral_rna_desc[rna_id]
            retrieval_stats["descriptions_attributed"] += 1
        if not "type" in attributes and rna_id in all_rnacentral_rna_type:
            attributes["type"] = all_rnacentral_rna_type[rna_id]
            retrieval_stats["types_attributed"] += 1

        retrieval_stats["total"] += 1
        attr_str2 = get_gff_attributes_str(attributes)
        info_fields2.append(attr_str2)
    
    print('Saving to', output_path)
    annotation["attribute"] = info_fields2
    annotation.to_csv(output_path, sep="\t", index=False, header=False)

    return retrieval_stats

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

def retrieve_func_annotation(annotation_path, output, confs):
    annotation = pd.read_csv(annotation_path, sep="\t", header=None, 
                names=["seqname", "source", "feature", "start", "end", 
                        "score", "strand", "frame", "attribute"])
    
    info_fields1 = annotation["attribute"].tolist()

    print('Loading RNACentral tables', file=sys.stderr)
    global_data = os.path.dirname(os.path.realpath(__file__)) + "/../data"
    rnacentral2go = {}
    with gzip.open(global_data+"/rnacentral2go.tsv.gz",'rt') as input_stream:
        for line in input_stream.readlines():
            cells = line.rstrip("\n").split()
            rnacentral = "URS"+(cells[0].split('_')[0])
            goid = "GO:"+cells[1]
            if not rnacentral in rnacentral2go:
                rnacentral2go[rnacentral] = set()
            rnacentral2go[rnacentral].add(goid)

    print("Loading rfam2go associations")
    rfam2go = read_rfam2go(global_data + "/rfam2go")

    print('Getting info from tables', file=sys.stderr)
    results = set()
    retrieval_stats = {"total_with_func": 0, "by_rfam": 0,
        "by_rnacentral": 0, "no_func": 0}
    
    for attr_str in tqdm(info_fields1):
        attributes = get_gff_attributes(attr_str)
        raw_id = attributes["ID"]

        first_part = ".".join(raw_id.split(".")[:-1])
        rna_id = first_part.split('_')[0]
        results_before = retrieval_stats['by_rnacentral'] + retrieval_stats['by_rfam']
        if rna_id in rnacentral2go:
            for goid in rnacentral2go[rna_id]:
                new_item = (raw_id, goid)
                results.add(new_item)
            retrieval_stats['by_rnacentral'] += 1
        
        if 'rfam' in attributes:
            rfam_id = attributes['rfam']
            if rfam_id in rfam2go:
                added = False
                for goid in rfam2go[rfam_id]:
                    new_item = (raw_id, goid)
                    if not new_item in results:
                        added = True
                    results.add(new_item)
                if added:
                    retrieval_stats['by_rfam'] += 1
        
        results_after = retrieval_stats['by_rnacentral'] + retrieval_stats['by_rfam']
        if results_after > results_before:
            retrieval_stats['total_with_func'] += 1
        else:
            retrieval_stats['no_func'] += 1
    
    print(retrieval_stats)

    term_ontologies = get_term_ontology(confs["go_obo"])

    if len(results) > 0:
        results = list(results)
        results.sort()
        with open(output, 'w') as stream:
            for id_, go in results:
                ont = term_ontologies[go] if go in term_ontologies else ''
                stream.write(id_+"\t"+go+"\t"+ont+"\n")
    
    return retrieval_stats

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


def listFD(url, ext=''):
    page = requests.get(url).text
    soup = BeautifulSoup(page, 'html.parser')
    return [url + '/' + node.get('href') for node in soup.find_all('a') if node.get('href').endswith(ext)]

def download_rnacentral_details(rnacentral_details_output):
    data_jsons_urls = listFD('https://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/json/', 'json')
    #data_jsons_urls = [data_jsons_urls[0], data_jsons_urls[1]]
    filenames = ['./'+os.path.basename(d) for d in data_jsons_urls]
    dumps = ['data/'+os.path.basename(d).replace('.json', '.tsv') for d in filenames]
    for i in tqdm(range(len(data_jsons_urls))):
        d = data_jsons_urls[i]
        dump = dumps[i]
        filename = filenames[i]
        if not os.path.exists(dump):
            output = open(dump, 'w')
            runCommand('wget --quiet '+d)
            for entry in json.load(open(filename, 'r')):
                if 'rna_type' in entry and 'description' in entry and 'rnacentral_id' in entry:
                    tp = entry['rna_type']
                    desc = entry['description']
                    id = entry['rnacentral_id']
                    output.write(id+'\t'+tp+'\t'+desc+'\n')
            output.close()
        #runCommand('rm '+filename)
    output.close()

    main_output = gzip.open(rnacentral_details_output, 'wt')
    for i in range(len(data_jsons_urls)):
        d = data_jsons_urls[i]
        dump = dumps[i]
        filename = filenames[i]
        for rawline in open(dump, 'r'):
            main_output.write(rawline)
        runCommand('rm ' + filename + ' ' + dump)
    main_output.close()
    quit()

def download_rnacentral2rfam2(rnacentral2rfam2_path):
    download_raw_path = os.path.dirname(rnacentral2rfam2_path) + '/rfam_annotations.tsv.gz'
    if not os.path.exists(download_raw_path):
        runCommand("cd "+os.path.dirname(rnacentral2rfam2_path)
            +"&&  wget --quiet https://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/rfam/rfam_annotations.tsv.gz")
    
    annots2 = set()
    rfam2desc = {}
    print('Reading', download_raw_path)
    for rawline in gzip.open(download_raw_path, 'rt'):
        cells = rawline.rstrip("\n").split('\t')
        rnacentral = cells[0]
        rfam = cells[1]
        desc = cells[8]
        rfam2desc[rfam] = desc
        annots2.add(rnacentral[3:]+"\t"+rfam[2:]+"\n")
    
    print('Sorting annotations read')
    annots2 = list(annots2)
    annots2.sort()

    print('Saving in compressed format')
    with gzip.open(rnacentral2rfam2_path, 'wt') as output:
        for rawline in annots2:
            output.write(rawline)

    print('Saving rfam descriptions')
    with gzip.open(os.path.dirname(rnacentral2rfam2_path) + '/rfam2desc.tsv.gz', 'wt') as output:
        for rfam, desc in rfam2desc.items():
            output.write(rfam + '\t' + desc + '\n')

def download_rnacentral2go(rnacentral2go_path):
    datadir = os.path.dirname(rnacentral2go_path)
    download_raw_path =  datadir + '/rnacentral_rfam_annotations.tsv.gz'
    if not os.path.exists(download_raw_path):
        runCommand("cd "+os.path.dirname(datadir)
            +"&&  wget --quiet https://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/go_annotations/rnacentral_rfam_annotations.tsv.gz")
    
    annots = set()
    print('Reading', download_raw_path)
    for rawline in gzip.open(download_raw_path, 'rt'):
        cells = rawline.rstrip("\n").split('\t')
        if len(cells) >= 2:
            rnacentral = cells[0]
            goid = cells[1]
            
            annots.add(rnacentral[3:]+"\t"+goid[3:]+"\n")
    
    print('Sorting go annotations read')
    annots = list(annots)
    annots.sort()

    print('Saving in compressed format at', rnacentral2go_path)
    with gzip.open(rnacentral2go_path, 'wt') as output:
        for rawline in annots:
            output.write(rawline)

if __name__ == "__main__":
    download_rnacentral_details('data/rnacentral_details.tsv.gz')
