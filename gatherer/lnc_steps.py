from gatherer.bioinfo import *
from gatherer.util import *
import pandas as pd
import os
from tqdm import tqdm

def filter_small_sequences(args, confs, tmpDir, stepDir):
    if "transcriptome" in args:
        seqs = readSeqsFromFasta(args["transcriptome"])
        longSeqs = list()

        smallestLength = 9999999
        filtered = list()
        for entry in seqs:
            if len(entry[1]) < smallestLength:
                smallestLength = len(entry[1])
            if len(entry[1]) >= 200:
                longSeqs.append(entry)
            else:
                filtered.append(entry)
        print(str(len(filtered)) + " transcripts were too short to be a lncRNA. Smallest: " + str(smallestLength))
        print(str(len(longSeqs)) + " long enough to be lncRNA.")
        print("Writing them down")
        output_fasta = tmpDir + "/long_transcripts.fasta"
        writeFastaSeqs(longSeqs, output_fasta)
    return True

def filter_long_orfs(args, confs, tmpDir, stepDir):
    #tmpDir = tempDir[1]
    longOrfs = confs["long_orfs"]
    long_transcripts = stepDir["filter_small_sequences"] + "/long_transcripts.fasta"
    if(os.path.exists(long_transcripts)):
        print("Finding ORFs with TransDecoder.LongORFs")
        longOrfsCmd = longOrfs + " -t " + long_transcripts
        code = runCommand("cd " + tmpDir + " && " + longOrfsCmd)
        if code != 0:
            return False
        orfsFile = tmpDir + "/long_transcripts.fasta.transdecoder_dir/longest_orfs.pep"
        orfs = readSeqsFromFasta(orfsFile)
        longOrfs = set()
        for entry in orfs:
            header = entry[0]
            seqLen = len(entry[1])
            name = header.rstrip("\n").lstrip(">").split()[0].split("|")[0].split(".p")[0]
            if seqLen > int(args["orf_length"]):
                longOrfs.add(name)
        print("Some longOrfs names: " + str(list(longOrfs)[:6]))
        #filtering out transcripts with long orfs
        seqs = readSeqsFromFasta(long_transcripts)
        noLongOrf = list()
        filtered = list()
        filtered, noLongOrf = filterSeqs(seqs, longOrfs)

        print(str(len(filtered)) + " transcripts filtered out")
        writeFastaSeqs(noLongOrf, tmpDir + "/no_orfs.fasta")
        return True
    else:
        return True

def test_coding_potential(args, confs, tmpDir, stepDir):
    long_transcripts_path = stepDir["filter_long_orfs"] + "/no_orfs.fasta"
    if os.path.exists(long_transcripts_path):
        output = tmpDir + "/samba.tsv"
        cmd = " ".join([confs["rnasamba"], "classify", output, 
                long_transcripts_path, confs["rnasamba_model"]])
        code = runCommand(cmd)
        if code != 0:
            return False
    return True

def parse_coding_potential(args, confs, tmpDir, stepDir):
    long_transcripts_path = stepDir["filter_long_orfs"] + "/no_orfs.fasta"
    output = stepDir["test_coding_potential"] + "/samba.tsv"
    if os.path.exists(output):
        non_coding_ids = set()
        with open(output, 'r') as stream:
            lines = [raw_line.rstrip("\n").split("\t") for raw_line in stream.readlines()]
            for cells in lines:
                if not cells[0].startswith("sequence_name"):
                    if cells[-1] == "noncoding":
                        ID = header_to_id(cells[0])
                        non_coding_ids.add(ID)
        seqs = readSeqsFromFasta(long_transcripts_path)
        noncoding, coding = filterSeqs(seqs, non_coding_ids)
        print(str(len(coding)) + " transcripts with coding potential filtered out")
        writeFastaSeqs(noncoding, tmpDir + "/no_coding_potential.fasta")
    return True

def nr_alignment(args, confs, tmpDir, stepDir):
    query_fasta = stepDir["parse_coding_potential"] + "/no_coding_potential.fasta"
    if os.path.exists(query_fasta):
        if "non_redundant" in confs:
            db = confs["non_redundant"]
            if os.path.exists(db):
                print("Starting diamond to blastx transcripts against NR proteins DB")
                diamond_cmd = [confs["diamond"],"blastx","--db",db,'-q',query_fasta,
                        "--out",tmpDir + "/blast.tsv","--outfmt","6",
                        "--threads",str(confs['threads']),"--evalue","0.0001"]
                code = runCommand("cd " + tmpDir + " && " + " ".join(diamond_cmd))
                return code == 0
            else:
                print("No NR database location defined.")
        else:
            print("No NR database location defined.")
    else:
        print("No coding potential evaluation, not doing lncRNA detection.")
    return True

def read_nr_alignment(args, confs, tmpDir, stepDir):
    blastFile = stepDir["nr_alignment"] + "/blast.tsv"
    if os.path.exists(blastFile):
        query_fasta = stepDir["parse_coding_potential"] + "/no_coding_potential.fasta"
        print("Filtering according to blastx results")

        seqs = readSeqsFromFasta(query_fasta)
        blast_df = pd.read_csv(blastFile,sep='\t',header=None)
        ids = set(blast_df[0].tolist())
        invalidSeqs, validSeqs = filterSeqs(seqs, ids)

        print(str(len(invalidSeqs)) + " transcripts with protein match filtered")
        print(str(len(validSeqs)) + " transcripts remaining")
        writeFastaSeqs(validSeqs, tmpDir + "/not_protein.fasta")
    else:
        print("Skiping NR alignment parsing")
    return True

def lnc_alignment_minimap(args, confs, tmpDir, stepDir):
    lnc_fasta = stepDir["read_nr_alignment"] + "/not_protein.fasta"
    if not os.path.exists(lnc_fasta):
        lnc_fasta = stepDir["parse_coding_potential"] + "/no_coding_potential.fasta"
    if os.path.exists(lnc_fasta):
        cmd = " ".join([confs["minimap2"],
            "-x splice:hq -uf -t", str(confs["threads"]),
            args["genome_index"], lnc_fasta, 
            ">", tmpDir + "/genome_mapping.paf"])
        code = runCommand(cmd)
        return code == 0 
    else:
        return True

def lnc_alignment_parsing(args, confs, tmpDir, stepDir):
    genome_alignment = stepDir["lnc_alignment_minimap"] + "/genome_mapping.paf"
    lnc_fasta = stepDir["read_nr_alignment"] + "/not_protein.fasta"
    if not os.path.exists(lnc_fasta):
        lnc_fasta = stepDir["parse_coding_potential"] + "/no_coding_potential.fasta"
    if os.path.exists(lnc_fasta) and os.path.exists(genome_alignment):
        output_gff = tmpDir + "/lncRNA_annotation.gff"
        annotated_fasta = tmpDir + "/lncRNA_in_genome.fasta"
        mapped_fasta = minimap_annotation(genome_alignment, output_gff,
                        annotated_fasta, source="rnasamba", mol_type="lncRNA",
                        query_file=lnc_fasta)
        return True
    else:
        print("Skiping lnc mapping")
    return True