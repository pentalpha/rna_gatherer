#Module: config.py
#Rename this file as "config.py" and customize it
import multiprocessing
import os
import json
import psutil
from gatherer.util import runCommand
from gatherer.netutils import download_to
from pathlib import Path

project_dir = os.path.dirname(os.path.realpath(__file__))
user_data_dir = str(Path.home())+"/.rna-gatherer_downloads"
config_json = project_dir + "/config.json"
if not os.path.exists(config_json):
    print("config.json file not found, creating a new one.")
    runCommand("cp " + project_dir+"/config.dummy.json " + config_json)

config_str = open(config_json, 'r').read()
config_str = config_str.replace("PROJECT_DIR", project_dir)
configs = json.loads(config_str)
if not "threads" in configs:
    configs["threads"] = str(max(2, multiprocessing.cpu_count()-1))
if not "max_mem" in configs:
    configs["max_mem"] = str(int((psutil.virtual_memory().total/1024/1024)*0.8))

print("Using " + configs["threads"] + " threads and " 
    + configs["max_mem"] + " MB of RAM.")

def missing_from_config(possibly_missing):
    missing = []
    for x in possibly_missing:
        if not x in configs:
            missing.append(x)
        else:
            if configs[x] == "":
                missing.append(x)
            elif x == "rna_dbs":
                if len(configs[x]) == 0:
                    missing.append(x)
    return missing

def download_go_obo():
    download_location = user_data_dir + "/" + "go.obo"
    success = download_to(configs['go_obo_url'], download_location)
    return success, download_location

def download_rfam_cm():
    download_location = user_data_dir + "/" + "Rfam.cm.gz"
    success = download_to(configs['rfam_cm_url'], download_location)
    clanin_path = configs['rfam_cm_url'].replace(".cm.gz",".clanin")
    clanin_location = user_data_dir + "/" + "Rfam.clanin"
    if success:
        success = download_to(clanin_path, clanin_location)
    rfam_location = user_data_dir + "/" + "Rfam.cm"
    if success:
        #press cm scan file
        runCommand("gzip -d " + download_location)
        rfam_location = download_location.replace(".gz","")
        runCommand("cmpress " + rfam_location)

    return success, rfam_location

def download_missing_files(missing):
    locations = {}
    failed = []
    if not os.path.exists(user_data_dir):
        os.mkdir(user_data_dir)
    for f in missing:
        success = False
        location = ''
        if 'obo' in f:
            success, location = download_go_obo()
        elif 'cm' in f:
            success, location = download_rfam_cm()
        if success:
            locations[f] = location
    
    for key in locations.keys():
        missing.remove(key)

    return missing, locations

def rename_config(locations):
    for key, location in locations.items():
        modified = ""
        with open(config_json, 'r') as in_stream:
            for line in in_stream:
                if "\""+key+"\"" in line or "'"+key+"'" in line:
                    updated_field = "\""+key+"\": \"" + location + "\",\n"
                    print("Renamed", line.rstrip("\n"), "to", updated_field.rstrip("\n"))
                    modified += updated_field
                else:
                    modified += line
        open(config_json,'w').write(modified)

def require_files(possibly_missing, mandatory = True):
    missing = missing_from_config(possibly_missing)
    to_download = []
    if len(missing) > 0 and not(mandatory):
        print("The following optional files are missing: " + str(missing))
        print("If you consider their data necessary, insert their locations in the " 
                + "configuration file:\n\t" + config_json)
    elif len(missing) > 0 and mandatory:
        print("The following mandatory files are missing: " + str(missing))
        print("Download them or manually insert their paths in the configuration file?")
        print("\tY (download) | N (define custom paths): ")
        downloaded = False
        answer = input('>')
        if answer in ['y','Y','Yes','yes']:
            still_missing, locations = download_missing_files(missing)
            if len(locations.keys()) > 0:
                rename_config(locations)
            if len(still_missing) == 0:
                downloaded = True
            else:
                print(str(still_missing) + " could not be downloaded.")
        if not downloaded:
            print("Insert their locations in the configuration file:\n\t" + config_json)
            quit()
            