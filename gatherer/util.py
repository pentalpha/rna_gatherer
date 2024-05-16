import subprocess
import os
from sys import getsizeof
from tqdm import tqdm

def runCommand(cmd, print_cmd=True):
        if print_cmd:
                print("\t> " + cmd)
        process = subprocess.call(cmd, shell=True)
        return process

def sliceString(string, max_len):
        slices = []
        for x in range(0,int(len(string)/max_len)+1):
                slices.append(string[x*max_len:(x+1)*max_len])
        return slices

def write_file(content, path):
        with open(path, 'w') as stream:
                stream.write(content)

def getFilesWith(directory, name_part, ending=False):
        files = []
        # r=root, d=directories, f = files
        for r, d, f in os.walk(directory):
                for file in f:
                        success = (name_part in file and not ending) or (file.endswith(name_part) and ending)
                        if success:
                                files.append(os.path.join(r, file))
        return files

def rm_last_part(path, c):
        parts = path.split(c)
        if len(parts) > 1:
                return c.join(parts[:-1])
        else:
                return parts[0]

def replace_last(original, to_repl, repl):
        k = original.rfind(to_repl)
        new_string = original[:k] + repl + original[k+len(to_repl):]
        return new_string

def file_name(path):
        no_file_extension = rm_last_part(os.path.basename(path), ".")
        return no_file_extension

def get_subdirs(dir):
        subs = os.listdir(dir)
        subdirs = []
        for sub in subs:
                path = dir + "/" + sub
                if os.path.isdir(path):
                        subdirs.append(path)
        return subdirs

def chunks(lst, n):
        """Yield successive n-sized chunks from lst."""
        for i in range(0, len(lst), n):
                yield lst[i:i + n]

def get_iterator(my_list, show = False):
        if show:
                return tqdm(range(len(my_list)))
        else:
                return range(len(my_list))

def get_cache(usage=1.0):
        default = 4*1024*1024

        lscpu_out = subprocess.check_output("lscpu | grep 'cache'", shell=True).decode("utf-8").split("\n")
        #print(str(lscpu_out))
        kilo_byte_notations = ["KB","K","k","kb","KiB","kib"]
        mega_byte_notations = ["MB","M","m","mb","MiB","mib"]
        if len(lscpu_out) >= 2:
                raw_nums = []
                for line in lscpu_out:
                        if len(line) > 0:
                                last = line.split()[-1]
                                if last in kilo_byte_notations or last in mega_byte_notations:
                                        last = line.split()[-2]+line.split()[-1]
                                raw_nums.append(last)
                #print(str(raw_nums))
                nums = []
                
                if len(raw_nums) > 0:
                        for raw_num in raw_nums:
                                for notation in kilo_byte_notations:
                                        if raw_num.endswith(notation):
                                                try:
                                                        nums.append(int(raw_num.replace(notation, ""))*1024)
                                                except ValueError:
                                                        print("Could not convert '" + raw_num + "' to bytes.")
                                for notation in mega_byte_notations:
                                        if raw_num.endswith(notation):
                                                try:
                                                        nums.append(int(raw_num.replace(notation, ""))*1024*1024)
                                                except ValueError:
                                                        print("Could not convert '" + raw_num + "' to bytes.")
                        greater = 0
                        for num in nums:
                                if num > greater:
                                        greater = num
                        if greater > 0:
                                default = greater
                                #print("Available cache memory is " + str(default))
        #quit()
        #print(lscpu_out)
        return int((default)*usage)

# input - df: a Dataframe, chunkSize: the chunk size
# output - a list of DataFrame
# purpose - splits the DataFrame into smaller of max size chunkSize (last is smaller)
def splitDataFrameIntoSmaller(df, chunkSize = 10000): 
        listOfDf = list()
        numberChunks = len(df) // chunkSize + 1
        for i in range(numberChunks):
                listOfDf.append(df[i*chunkSize:(i+1)*chunkSize])
        if len(listOfDf[-1]) == 0:
                return listOfDf[:-1]
        else:
                return listOfDf

def split_df_to_max_mem(df, available_size, max_lines=2000, min_lines=1000):
        original_ids = set(df[df.columns[0]].tolist())
        print("\tAvailable KB: " + str(available_size/1024))
        lines_per_part = len(df)
        counts_size = getsizeof(df)
        print("\tCounts table size: " + str(counts_size/1024) + "KB")
        dfs = [df]
        percent_per_part = float(available_size) / float(counts_size)
        print("\tpercent per part: " + str(percent_per_part))
        lines_per_part = int(len(df)*percent_per_part)
        if lines_per_part > max_lines:
                lines_per_part = max_lines
        elif lines_per_part < min_lines:
                lines_per_part = min_lines
        print("\tlines per part: " + str(lines_per_part))
        dfs = splitDataFrameIntoSmaller(df, chunkSize=lines_per_part)
        print("\tSplitted into " + str(len(dfs)))
        id_found = set()
        for data in dfs:
                id_found.update(data[data.columns[0]].tolist())
        print("\tIDs expected/found:", len(original_ids), "/", len(id_found))
        return dfs

def read_to_list(file_path):
        if os.path.exists(file_path):
                with open(file_path,'r') as stream:
                        return [line.rstrip("\n") for line in stream.readlines()]
        else:
                return None

def delete_if_empty(file_path, min_cells = 1, sep = "\t"):
        if os.path.exists(file_path):
                empty = True
                with open(file_path,'r') as stream:
                        line = stream.readline()
                        while len(line) > 0:
                                cells = line.split(sep)
                                if len(cells) >= min_cells:
                                        empty = False
                                        break
                                line = stream.readline()
                if empty:
                        #print("Removing " + file_path + " because it's empty.")
                        os.remove(file_path)

def join_files_in_one(inputs, output):
        wrote = False
        with open(output, 'w') as out_stream:
                for input_file in inputs:
                        if os.path.exists(input_file):
                                wrote = True
                                with open(input_file, 'r') as in_stream:
                                        for line in in_stream:
                                                out_stream.write(line)
        if not wrote and os.path.exists(output):
                os.remove(output)
                return False
        return wrote

