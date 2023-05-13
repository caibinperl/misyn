import sys
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class Gff3(object):
    def __init__(self, gff_file):
        self.chr2genes = {}
        self.gene2rnas = {}
        self.gene2long_rna = {}
        self.rna2cdss = {}
        self.rna2utrs = {}
        self.rna2exons = {}
        self.gene_id2gene = {}

        for line in open(gff_file, 'rt'):
            line = line.strip()
            items = line.split('\t')
            if len(items) != 9:
                continue
            items[3] = int(items[3])
            items[4] = int(items[4])
            chr_id = items[0]
            source = items[1]
            feature = items[2]
            start = items[3]
            end = items[4]
            score = items[5]
            strand = items[6]
            phase = items[7]
            attr = items[8]
            if start > end:
                print('position erro for line: ' + line)
                # sys.exit(1)
            if feature == 'gene':
                m = re.search('ID=([^;]+)', attr)
                if m:
                    gene_id = m.group(1)
                    gene = items[0:8]
                    gene.append(gene_id)
                    self.gene_id2gene[gene_id] = gene
                    if chr_id in self.chr2genes:
                        self.chr2genes[chr_id].append(gene)
                    else:
                        self.chr2genes[chr_id] = [gene]
                else:
                    print('Error gene: ' + line)
            if feature == 'mRNA':
                m = re.search('ID=([^;]+).*Parent=([^;]+)', attr)
                if m:
                    rna_id = m.group(1)
                    gene_id = m.group(2)
                    rna = items[0:8]
                    rna.append(rna_id)
                    if gene_id in self.gene2rnas:
                        self.gene2rnas[gene_id].append(rna)
                    else:
                        self.gene2rnas[gene_id] = [rna]

            if feature == 'CDS':
                m = re.search('Parent=([^;]+)', attr)
                if m:
                    rna_id = m.group(1)
                    cds = items[0:8]
                    cds.append(rna_id)
                    if rna_id in self.rna2cdss:
                        self.rna2cdss[rna_id].append(cds)
                    else:
                        self.rna2cdss[rna_id] = [cds]
                    if rna_id in self.rna2exons:
                        self.rna2exons[rna_id].append(cds)
                    else:
                        self.rna2exons[rna_id] = [cds]
                else:
                    print('Error CDS: ' + line)
            if feature.upper().find('UTR') != -1:
                m = re.search('Parent=([^;]+)', attr)
                if m:
                    rna_id = m.group(1)
                    utr = items[0:8]
                    utr[2] = 'UTR'
                    utr.append(rna_id)
                    if rna_id in self.rna2exons:
                        self.rna2exons[rna_id].append(utr)
                    else:
                        self.rna2exons[rna_id] = [utr]
                    if rna_id in self.rna2utrs:
                        self.rna2utrs[rna_id].append(utr)
                    else:
                        self.rna2utrs[rna_id] = [utr]
                else:
                    print(sys.stderr, 'Error UTR: ' + line)
        for chr_id in self.chr2genes:
            genes = self.chr2genes[chr_id]
            genes.sort(key=lambda e: e[3])
        for gene_id in self.gene2rnas:
            rnas = self.gene2rnas[gene_id]
            rnas.sort(key=lambda e: e[3])
        for rna_id in self.rna2cdss:
            cdss = self.rna2cdss[rna_id]
            cdss.sort(key=lambda e: e[3])
        for rna_id in self.rna2exons:
            exons = self.rna2exons[rna_id]
            exons.sort(key=lambda e: e[3])
        for rna_id in self.rna2utrs:
            utrs = self.rna2utrs[rna_id]
            utrs.sort(key=lambda e: e[3])
        self.ex_longest_rna()

    def ex_longest_rna(self):
        for gene_id in self.gene2rnas:
            rnas = self.gene2rnas[gene_id]
            max_len = 0
            long_rna = None
            for rna in rnas:
                rna_id = rna[8]
                cdss = self.rna2cdss[rna_id]
                gene_len = 0
                for cds in cdss:
                    cds_len = cds[4] - cds[3] + 1
                    gene_len += cds_len
                if gene_len > max_len:
                    max_len = gene_len
                    long_rna = rna
            self.gene2long_rna[gene_id] = long_rna


def fine_gene_index(genes, gene_id):
    for i in range(0, len(genes)):
        if genes[i][8] == gene_id:
            return i
    return -1


def get_flank_genes(spe, fam_gene_id_file, gff_file, flank_gene_num):
    fam_gene_ids = [line.strip() for line in open(fam_gene_id_file, 'rt')]
    gff = Gff3(gff_file)
    flank_genes = []
    for fam_gene_id in fam_gene_ids:
        fam_gene = gff.gene_id2gene[fam_gene_id]
        genes = gff.chr2genes[fam_gene[0]]
        fam_gene_i = fine_gene_index(genes, fam_gene_id)
        if fam_gene_i == -1:
            print(f'Can\'t find gene {fam_gene_id}')
            sys.exit(1)
        start = fam_gene_i - flank_gene_num
        if start < 0:
            start = 0
        end = fam_gene_i + flank_gene_num
        if end > len(genes):
            end = len(genes)
        for gene in genes[start: end]:
            flank_genes.append(
                [fam_gene_id, gene[8], spe, gene[0], gene[3], gene[4], gene[6]])

    return flank_genes


def ex_pep(genes, gff_file, genome_file):
    prot_records = []
    gff = Gff3(gff_file)
    id2seq = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))
    for gene in genes:
        gene_id = gene[1]
        chr_id = gene[3]
        long_rna = gff.gene2long_rna[gene_id]
        cdss = gff.rna2cdss[long_rna[8]]
        chr_seq = id2seq[chr_id].seq
        rna_seq = Seq('')
        for cds in cdss:
            start = cds[3] - 1
            end = cds[4]
            cds_seq = chr_seq[start: end]
            rna_seq += cds_seq

        if long_rna[6] == '-':
            rna_seq = rna_seq.reverse_complement()
        prot_seq = rna_seq.translate()
        prot_record = SeqRecord(prot_seq, id=gene_id)
        prot_records.append(prot_record)
    return prot_records


def get_chr_len(genome_file):
    chr2len = {}
    id2seq = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))
    for chr_id in id2seq:
        chr2len[chr_id] = len(id2seq[chr_id].seq)
    return chr2len
