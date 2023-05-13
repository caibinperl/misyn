import os
import argparse
import utils


parser = argparse.ArgumentParser()
parser.add_argument('meta', help='metadata file')
parser.add_argument('output', help='output directory')
parser.add_argument('-f', '--flank', help='flank flanking gene number',
                    type=int, default=50)
parser.add_argument('-b', '--blast_e', help='e-value for BLAST', type=float,
                    default=0.01)
parser.add_argument('-t', '--threads', help='threads', type=int, default=1)
args = parser.parse_args()

metadata_file = args.meta
output_dir = args.output
flank_gene_num = args.flank
blast_e = args.blast_e

chr_len_dir = os.path.join(output_dir, 'chr_len')
pep_file = os.path.join(output_dir, 'pep.fasta')
gene_list_file = os.path.join(output_dir, 'gene_list.tsv')
homology_file = os.path.join(output_dir, 'homology.tsv')

if not os.path.isdir(output_dir):
    os.mkdir(output_dir)

if not os.path.isdir(chr_len_dir):
    os.mkdir(chr_len_dir)

list_fout = open(gene_list_file, 'wt')
pep_fout = open(pep_file, 'wt')

lines = [line.strip() for line in open(metadata_file, "rt")]
for i in range(1, len(lines)):
    line = lines[i]
    if line.startswith('#'):
        continue
    ss = line.split('\t')
    if len(ss) != 4:
        continue

    spe, fam_gene_id_file, gff_file, genome_file = ss
    flank_genes = utils.get_flank_genes(spe, fam_gene_id_file, gff_file, 
        flank_gene_num)
    for gene in flank_genes:
        fam_gene_id, gene_id, spe, chr_id, start, end, strand = gene
        list_fout.write(f'{fam_gene_id}\t{gene_id}\t{spe}\t{chr_id}\t{start}\t{end}\t{strand}\n')

    prot_records = utils.ex_pep(flank_genes, gff_file, genome_file)

    for prot_record in prot_records:
        pep_fout.write('>' + prot_record.id + '\n')
        pep_fout.write(str(prot_record.seq) + '\n')

    chr2len = utils.get_chr_len(genome_file)
    chr_len_fout = open(os.path.join(chr_len_dir, f'{spe}_len.tsv'), 'wt')
    for chr_id in sorted(chr2len.keys()):
        chr_len = chr2len[chr_id]
        chr_len_fout.write(f'{chr_id}\t{chr_len}\n')
    chr_len_fout.close()

pep_fout.close()

blast_out_file = os.path.join(output_dir, 'blast_output.tsv')
cmd = f'makeblastdb -in {pep_file} -dbtype prot'
os.system(cmd)
cmd = f'blastp -query {pep_file} -db {pep_file} -out {blast_out_file} -outfmt 6 -evalue {blast_e}'
os.system(cmd)

homology_fout = open(homology_file, 'wt')

for line in open(blast_out_file, 'rt'):
    line = line.strip()
    ss = line.split('\t')
    if len(ss) != 12:
        continue
    if ss[0] == ss[1]:
        continue
    homology_fout.write(ss[0] + '\t' + ss[1] + '\t' + ss[10]  + '\n')

homology_fout.close()

cmd = f'./misyn -l {gene_list_file} -s {homology_file} -o {output_dir} -c {chr_len_dir}'
os.system(cmd)

print("Finished!")