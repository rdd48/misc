import sys
import math
import subprocess
import pickle

if len(sys.argv) < 2:
    print('USAGE: python efi_cacao.py INPUT_FASTA')
    exit()
fasta_file = sys.argv[1]

def process_fasta(fasta_file):
    with open(fasta_file) as f:

        fasta_dict = {}

        lines = f.readlines()
        for i in range(len(lines)):
            if lines[i][0] == '>':
                name = lines[i][1:]
                seq = ''
                new_index = i + 1
                while lines[new_index][0] != '>':
                    seq += lines[new_index]
                    new_index += 1
                    if new_index == len(lines):
                        break
                    
                if ' ' in name:
                    name_only = name.split()[0]
                    if '[' in name:
                        species = name.split('[')[-1].replace(']', '').strip()
                    else:
                        species = ''

                    if len(name.split()) > 1:
                        descript = ' '.join(name.split()[1:])
                        descript = descript.replace('[' + species + ']', '').strip()
                    else:
                        descript = ''

                    fasta_dict[name_only] = (seq.replace('\n', ''), species, descript)

                else:
                    fasta_dict[name.strip()] = (seq.replace('\n', ''), '', '')
    
    return fasta_dict

def fasta_to_dicts(input_fasta):

    # keep info for xgmml
    nodes, edges = {}, {}

    fasta_dict = process_fasta(input_fasta)

    out = subprocess.check_output("blastp -query " + input_fasta + " -subject " + input_fasta + " -matrix BLOSUM62 -outfmt '6 qseqid qlen sseqid slen bitscore length pident'", shell=True)

    for line in out.decode('utf-8').split('\n'):
        if line:
            n1, l1, n2, l2, bitscore, align_len, pident = line.split('\t')

            # nodes dict structure = name: (length,  sequence, species, descript)
            if n1 not in nodes:
                seq, species, descript = fasta_dict[n1]
                nodes[n1] = (l1, seq, species, descript)
            if n2 not in nodes:
                seq, species, descript = fasta_dict[n2]
                nodes[n2] = (l2, seq, species, descript)

            if n1 != n2:

                # https://github.com/EnzymeFunctionInitiative/EFITools/blob/master/lib/EFI/Util/AlignmentScore.pm
                # is formula from ^ actually?:
                align_score = -(math.log(float(l1) * float(l2)) / math.log(10)) + (float(bitscore) * math.log(2) / math.log(10)) # log must be ln
                
                # edges dict structure = (name 1, name 2): (percent id, align_score, align_len)
                if align_score > 0.:
                    if (n1, n2) and (n2, n1) not in edges:
                        edges[(n1, n2)] = (pident, align_score, align_len)
                    elif (n1, n2) in edges:
                        if align_score > edges[(n1, n2)][1]:
                            edges[(n1, n2)] = (pident, align_score, align_len)
                    elif (n2, n1) in edges:
                        if align_score > edges[(n2, n1)][1]:
                            edges[(n2, n1)] = (pident, align_score, align_len)
        
    return nodes, edges

def dataset_analysis(nodes, edges):

    # keep info for descriptive graphs
    num_seqs, score_by_len, score_by_id, score_by_edge = {}, {}, {}, {}

    for v in nodes.values():
        seqlen = int(v[0])
        if seqlen not in num_seqs:
            num_seqs[seqlen] = 1
        else:
            num_seqs[seqlen] += 1
    
    for v in edges.values():
        pident, align_score, align_len = v

        pident = float(pident)
        align_score = int(align_score)
        align_len = int(align_len)

        align_score = round(align_score)

        if align_score not in score_by_len:
            score_by_len[align_score] = [align_len]
        else:
            v_copy = score_by_len[align_score][:]
            v_copy.append(align_len)
            score_by_len[align_score] = v_copy

        if align_score not in score_by_id:
            score_by_id[align_score] = [pident]
        else:
            v_copy = score_by_id[align_score][:]
            v_copy.append(pident)
            score_by_id[align_score] = v_copy

        if align_score not in score_by_edge:
            score_by_edge[align_score] = 1
        else:
            score_by_edge[align_score] += 1
    
    return num_seqs, score_by_len, score_by_id, score_by_edge

def pickle_dict(d, fname):
    with open(fname, 'wb') as fout:
        pickle.dump(d, fout)

def write_all_dicts(network_name, nodes, edges, num_seqs, score_by_len, score_by_id, score_by_edge):

    # # this would also work, but would create multiple files. 
    # pickle_dict(nodes, f'{network_name}_nodes.pickle')
    # pickle_dict(edges, f'{network_name}_edges.pickle')

    dict_of_dicts = {
        'nodes': nodes,
        'edges': edges,
        'num_seqs': num_seqs,
        'score_by_len': score_by_len,
        'score_by_id': score_by_id,
        'score_by_edge': score_by_edge
    }

    pickle_dict(dict_of_dicts, network_name + '_data.pickle')

if __name__ == '__main__':
    nodes, edges = fasta_to_dicts(fasta_file)
    network_name = fasta_file.split('/')[-1].split('.fasta')[0].strip()
    num_seqs, score_by_len, score_by_id, score_by_edge = dataset_analysis(nodes, edges)

    write_all_dicts(network_name, nodes, edges, num_seqs, score_by_len, score_by_id, score_by_edge)
