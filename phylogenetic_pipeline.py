import re
import os
import pandas as pd
import subprocess
import sys
import multiprocessing
from ete3 import Tree, TreeStyle
from Bio import Entrez, SeqIO, AlignIO
from urllib import error
from typing import List
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.Phylo import write
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceCalculator


def mk_dirs(dir_name: str, a_dir: str):
    """
    Creates directories essential for pipeline to work properly.
    :param dir_name: str, name of directory to create
    :param a_dir: str, name of the directory where the analysis is currently being carried out.
    :return:
    """
    if not os.path.exists(f"{a_dir}/{dir_name}"):
        os.mkdir(f"{a_dir}/{dir_name}")
        print(f"Directory '{dir_name}' status: created")
    else:
        print(f"Directory '{dir_name}' status: already exists")


def read_list_of_ids(file: str) -> List[str]:
    """
    Reads given file and creates List with unique IDs.
    :param file: str, name of file.
    :raise: Exception when can't find a file.
    :return: List with IDs.
    """
    try:
        data = ''
        with open(file, 'r') as f:
            for line in f.readlines():
                line = re.sub('[!@/#$|:]', ' ', line).strip()
                # Substitutes above signs with space ensuring your
                # file with IDs will be always read correctly.
                data = data + ' ' + line
            data = data.split()
            data = list(dict.fromkeys(data))
        return data
    except FileNotFoundError:
        raise Exception("!!!ERROR!!!\n"
                        "Are you sure, you passed right filename?"
                        "If it is not in the same directory as this script, you must specify the path.")


def download_ncbi(ncbi_id: str, a_dir: str, c_dir: str):
    """
    Fetching files using NCBI IDs in Entrez module and saving them in the indicated directory. In case of http error
    or after confirming that entry is obsolete, appends currently handled ID to the appropriate list.
    :param ncbi_id: str, ID confirmed to be from the NCBI database.
    :param c_dir: str, data directory, where the downloaded files are stored
    :param a_dir: str, name of the directory where the analysis is currently being carried out.
    """
    try:
        if os.path.exists(f'{a_dir}/{c_dir}/{ncbi_id}.fasta') is False:
            handle = Entrez.efetch(db="nuccore", id=ncbi_id, rettype="fasta_cds_aa", retmode="text")
            record = handle.read()

            with open(f'{a_dir}/{c_dir}/{ncbi_id}.fasta', 'w+') as f:
                f.write(record)
            if os.path.exists(f'{a_dir}/{c_dir}/{ncbi_id}.fasta') is False:
                print(f'Error in file {ncbi_id}')

    except error.HTTPError or TimeoutError or error.URLError or TypeError:
        print(f'Error in file {ncbi_id}')
    except FileNotFoundError:
        print(f'Error in file {ncbi_id}')


def clustering(df, a_dir: str, cut: int, indexes: list):
    """
    Gather genomic data from NCBI, creates 3 dictionaries that are used to retrieve information about clusters content
    and then places it into separate files for each cluster. Cluster files contain organism names and aa gene sequence.
    One per each organism.
    :param df: DataFrame, read mcl result file with clusters
    :param a_dir: str, name of the directory where the analysis is currently being carried out.
    :param cut: int, clusters with less genes than cutoff are dropped
    :param indexes: list of NCBI ids
    :return:
    """
    names = {}
    seqs = {}
    name_gene = {}
    count = 0

    print("Gathering NCBI genomic data...")
    counter = 1
    for index in indexes:
        handle = Entrez.efetch(db="nucleotide", id=index, rettype="gb", retmode="text")
        record = SeqIO.read(handle, 'genbank')
        record = re.sub(' ', '_', record.annotations['organism'])

        try:
            if index not in names.keys():
                names[index] = record
        except KeyError:
            names[index] = f"Unknown_{counter}"
            counter += 1

    for record in SeqIO.parse(f'{a_dir}/{family_name}.fasta', "fasta"):
        if record.id not in seqs.keys():
            seqs[record.id] = record.seq

    for index in names.keys():
        for prot in seqs.keys():
            x = re.search('([A-Z]){2,3}\w{7}\.\d{1}', prot)
            if index == x.group(0):
                if f"{names[index]}" in name_gene:
                    name_gene[f"{names[index]}"].append(prot)
                else:
                    name_gene[f"{names[index]}"] = [prot]

    print("Gathering clusters...")
    for index, series in df.iterrows():
        if str(series[cut]) != "nan" and str(series[cut]) != "None":
            count += 1
            unique = {}
            if os.path.exists(f'{a_dir}/clusters/cluster_{count}.fasta') is False:
                with open(f'{a_dir}/clusters/cluster_{count}.fasta', 'a+') as cluster:
                    for s in series:
                        if s in seqs.keys():
                            for name in name_gene.keys():
                                if s in name_gene[name]:
                                    if name not in unique.keys():
                                        unique[f">{name}"] = seqs[s]
                    for key in unique.keys():
                        cluster.write(f"{key}\n{unique[key]}\n\n")


def clustalo(clusters: int, a_dir: str):
    """
    Preforms MSA of every cluster file using ClustalOmega
    :param clusters: number of clusters in the current analysis.
    :param a_dir: str, name of the directory where the analysis is currently being carried out.
    """
    print("Performing MSA with ClustalOmega...")
    for cluster in range(1, clusters + 1):
        if os.path.exists(f"{a_dir}/MSA/MSA_cluster_{cluster}.fasta") is False:
            in_file = f"{a_dir}/clusters/cluster_{cluster}.fasta"
            out_file = f"{a_dir}/MSA/MSA_cluster_{cluster}.fasta"
            clustalo_cline = ClustalOmegaCommandline(infile=in_file, outfile=out_file, auto=True, verbose=False)
            clustalo_cline()


def construct_nj_trees(clusters: int, a_dir: str):
    """
    Constructs NJ trees with given MSAs and saves them into newick-format file.
    :param clusters: Number of clusters in current analysis.
    :param a_dir: str, name of the directory where the analysis is currently being carried out.
    """
    print("Constructing NJ trees...")
    if os.path.exists(f"{a_dir}/supertree"):
        with open(f"{a_dir}/supertree/nj_source_trees_{family_name}.newick", 'a') as nj_file:
            for cluster in range(1, clusters + 1):
                msa = AlignIO.read(open(f"{a_dir}/MSA/MSA_cluster_{cluster}.fasta"), 'fasta')
                calculator = DistanceCalculator('blosum62')
                dm = calculator.get_distance(msa)
                constructor = DistanceTreeConstructor()
                tree = constructor.nj(dm)
                with open(f"{a_dir}/supertree/temp_{family_name}.newick", 'w') as nj_temp:
                    write(tree, nj_temp, format='newick')
                with open(f"{a_dir}/supertree/temp_{family_name}.newick", 'r') as temp:
                    newick = temp.read()
                    newick = re.sub('-([0-9])\.+([0-9])\d+', '0.0', newick)
                nj_file.write(newick)
    else:
        raise Exception("!!!ERROR!!!\n"
                        "I don't see 'supertree' directory in your analysis directory.\n"
                        "Please run mode 'p' or create suitable directory yourself.")


def bash_process(command: str):
    """
    Runs some of the 3rd party programs that have to be run from command-line.
    :param command: Bash command with arguments to run.
    """
    process = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
    stdout = process.communicate()
    print(f"{command.split()[0]} done.")
    return stdout


def draw_supertree(method: str, scale: int, style: str = "r"):
    """
    Drawing supertree using ete3 module and saving it to a .png file.
    :param method: str, method used to build supertree
    :param scale: int, scale branch lengths
    :param style: str, "c" for circular tree, "r" for regular
    """
    with open(f'{method}_supertree.txt', 'r') as super_file:
        text = super_file.readlines()
        st = re.sub('-([0-9])\.+([0-9])\d+', '0.0', text[0])
    t = Tree(f"{st.split(';')[0]};")
    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.show_branch_length = True
    ts.scale = scale
    ts.mode = style
    t.render(f"{method}_supertree.png", tree_style=ts)


if __name__ == '__main__':
    try:
        analysis_dir = sys.argv[1]
        filename = sys.argv[2]
        family_name = sys.argv[3]
        cutoff = sys.argv[4]
        supertree = sys.argv[5]
        mode = sys.argv[6]
        mode = [char for char in mode]

        ids = read_list_of_ids(filename)
        clear_ids = list(dict.fromkeys(ids))
        Entrez.email = ''
        data_dir = "data"

        if 'p' in mode:
            if not os.path.exists(analysis_dir):
                os.mkdir(analysis_dir)
                print(f"Directory '{analysis_dir}' status: created")
            else:
                print(f"Directory '{analysis_dir}' status: already exists")

            dirs = ['data', 'clusters', 'MSA', 'supertree']
            for single_dir in dirs:
                mk_dirs(single_dir, analysis_dir)

        if 'd' in mode:
            if os.path.exists(f"{analysis_dir}/data"):
                print("Downloading...")
                for idx in clear_ids:
                    if re.match(r'([A-Z]){2,3}\w+', idx) or re.match(r'([A-Z]){2,3}\w+\.\w+', idx):
                        download_ncbi(idx, analysis_dir, data_dir)
                print("Downloading done.")
            else:
                raise Exception("!!!ERROR!!!\n"
                                "I don't see 'data' directory in your analysis directory.\n"
                                "Please run mode 'p' or create suitable directory yourself.")

        if 'b' in mode:
            if os.path.exists(f"{analysis_dir}/data"):
                if len(os.listdir(f"{analysis_dir}/data")) != 0:
                    os.system(f"cat {analysis_dir}/data/* > {analysis_dir}/{family_name}.fasta")
                    print("Cat done.")
                    bdb_command = f"makeblastdb -in {analysis_dir}/{family_name}.fasta -dbtype prot"
                    bash_process(bdb_command)
                    print("BLASTing...")
                    blast_command = f"blastp -query {analysis_dir}/{family_name}.fasta " \
                                    f"-db {analysis_dir}/{family_name}.fasta " \
                                    f"-out {analysis_dir}/blast_{family_name}.fasta -outfmt 6"
                    bash_process(blast_command)
                else:
                    raise Exception("!!!ERROR!!!\n"
                                    "It looks like your 'data' directory is empty.\n"
                                    "You can't blast if there is nothing to blast."
                                    "Please place your fastas in 'data' directory or use mode 'd'.")
            else:
                raise Exception("!!!ERROR!!!\n"
                                "I don't see 'data' directory in your analysis directory.\n"
                                "Please run mode 'p' or create suitable directory yourself.")

        if 'c' in mode:
            if os.path.exists(f"{analysis_dir}/clusters"):
                if os.path.exists(f"{analysis_dir}/blast_{family_name}.fasta"):
                    os.system(f"cut -f 1,2,11 {analysis_dir}/blast_{family_name}.fasta "
                              f"> {analysis_dir}/{family_name}.abc")
                    print("Cut done.")
                    print("Clustering...")
                    mcl_command = f"mcl {analysis_dir}/{family_name}.abc --abc  " \
                                  f"-o {analysis_dir}/mcl_{family_name}.txt " \
                                  f"-te {str(multiprocessing.cpu_count())} -V all"
                    bash_process(mcl_command)
                    mcl = pd.read_csv(f"{analysis_dir}/mcl_{family_name}.txt", header=None, sep="\t", engine='python')
                    clustering(mcl, analysis_dir, int(cutoff), clear_ids)
                    cluster_count = len([name for name in os.listdir(f"{analysis_dir}/clusters") if
                                         os.path.isfile(os.path.join(f"{analysis_dir}/clusters", name))])
                    print("Clustering done.")
                else:
                    raise Exception("!!!ERROR!!!\n"
                                    "It looks like there is no blast file to form clusters.\n"
                                    "Please place blast file in your analysis directory or use mode 'b'.")
            else:
                raise Exception("!!!ERROR!!!\n"
                                "I don't see 'clusters' directory in your analysis directory.\n"
                                "Please run mode 'p' or create suitable directory yourself.")

        if 'm' in mode:
            if os.path.exists(f"{analysis_dir}/clusters"):
                if len(os.listdir(f"{analysis_dir}/clusters")) != 0:
                    if os.path.exists(f"{analysis_dir}/MSA"):
                        cluster_count = len([name for name in os.listdir(f"{analysis_dir}/clusters") if
                                             os.path.isfile(os.path.join(f"{analysis_dir}/clusters", name))])
                        clustalo(cluster_count, analysis_dir)
                        print("MSA performed.")
                    else:
                        raise Exception("!!!ERROR!!!\n"
                                        "I don't see 'MSA' directory in your analysis directory.\n"
                                        "Please run mode 'p' or create suitable directory yourself.")
                else:
                    raise Exception("!!!ERROR!!!\n"
                                    "It looks like your 'clusters' directory is empty.\n"
                                    "You can't perform MSA if there are no sequences to align."
                                    "Please place right files in 'clusters' directory or use mode 'c'.")
            else:
                raise Exception("!!!ERROR!!!\n"
                                "I don't see 'clusters' directory in your analysis directory.\n"
                                "Please run mode 'p' or create suitable directory yourself.")

        if 't' in mode:
            if os.path.exists(f"{analysis_dir}/MSA"):
                if len(os.listdir(f"{analysis_dir}/MSA")) != 0:
                    cluster_count = len([name for name in os.listdir(f"{analysis_dir}/MSA") if
                                         os.path.isfile(os.path.join(f"{analysis_dir}/MSA", name))])
                    construct_nj_trees(cluster_count, analysis_dir)
                    os.system(f"rm {analysis_dir}/supertree/temp_{family_name}.newick")
                    print("NJ trees constructed.")
                else:
                    raise Exception("!!!ERROR!!!\n"
                                    "It looks like your 'MSA' directory is empty.\n"
                                    "You can't build trees if you don't have files with MSA"
                                    "Please place right files in 'MSA' directory or use mode 'm'.")
            else:
                raise Exception("!!!ERROR!!!\n"
                                "I don't see 'MSA' directory in your analysis directory.\n"
                                "Please run mode 'p' or create suitable directory yourself.")

        if 's' in mode:
            if os.path.exists(f"{analysis_dir}/supertree"):
                if os.path.exists(f"{analysis_dir}/supertree/nj_source_trees_{family_name}.newick"):
                    os.chdir(f"{analysis_dir}/supertree")
                    with open(f"clann_commands_{supertree}.txt", 'w+') as f:
                        f.write(f"exe nj_source_trees_{family_name}.newick\n"
                                f"{supertree} savetrees={supertree}_supertree.txt\nquit")
                    print("Building supertree...")
                    clann_command = f"clann --c clann_commands_{supertree}.txt"
                    bash_process(clann_command)
                    draw_supertree(supertree, 500)
                    print(f"{supertree} supertree constructed.")
                else:
                    raise Exception("!!!ERROR!!!\n"
                                    "It looks like there is no newick file to build a supertree.\n"
                                    "Please place newick file in the supertree directory or use mode 't'.")
            else:
                raise Exception("!!!ERROR!!!\n"
                                "I don't see 'supertree' directory in your analysis directory.\n"
                                "Please run mode 'p' or create suitable directory yourself.")
    except IndexError:
        raise Exception("!!!ERROR!!!\n"
                        "Something wrong with arguments.\n"
                        "Check README and pass argument correctly.")
