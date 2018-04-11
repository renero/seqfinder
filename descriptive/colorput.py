import sys
import csv
from termcolor import colored, cprint

def read_patterns(sequence_file):
    sequences = []
    with open(sequence_file, 'r') as seqs_file:
        seqs_reader = csv.reader(seqs_file)
        for line in seqs_reader:
            sequences.append(line)
    return sequences


symbols = ['??','__','BB','MA','MM','AL','OD','CP','GM','DD','GU','GA','GB','Gt','Gp','GQ','GH','Gc','GD','Gd',
    'Ge','GF','Gf','Gg','GI','Gi','GL','Ga','Gu','GR','Gr','Gm','GT','Gw','Re','Rg','Rp','Rl','Rq','RW','RS','RC',
    'Rt','RF','Rs','RA','Ra','RB','RG','Rc','RT','RM','Rw','Rv','Ri','FO','Fc','FR','FA','Fx','FC','FN','FB','FD',
    'FI','FP','Fp','SI','Sn','SE','Si','SN','AA','AD','AP','NN','NP','Nn','Np']

fg = ['grey','red','green','yellow','blue','magenta','cyan','white']
bg = ['on_grey','on_red','on_green','on_yellow','on_blue','on_magenta','on_cyan','on_white']

if len(sys.argv) is not 2:
    print("usage: ", sys.argv[0], "<patterns_file>")
    exit(1)

lines = read_patterns(sys.argv[1])
for line in lines:
    pat = line[1]
    for pos in range(0,len(pat),2):
        gene = '{:s}'.format(pat[pos:pos+2])
        col_idx = symbols.index(gene)
        col_pos = col_idx % len(fg)
        cprint(gene, fg[col_pos], end='')
    print()

