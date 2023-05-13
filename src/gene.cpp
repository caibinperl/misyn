#include "gene.h"
#include <iostream>

Gene::Gene(const string &gene_id, const string &spe_id, const string &chr_id, long start,
           long end, int strand)
{
    gene_id_ = gene_id;
    spe_id_ = spe_id;
    chr_id_ = chr_id;
    start_ = start;
    end_ = end;
    strand_ = strand;
    is_tandem_ = false;
    is_tandem_represent_ = false;
    tandem_represent_ = "";
    is_remapped_ = false;
    is_fam_gene_ = false;
}

Gene::~Gene()
{
}

// Add a new gene to homologyIds

void Gene::add_homo_ids(const string &gene_id)
{
    homo_ids_.insert(gene_id);
}

void Gene::add_homo_gene(Gene *gene)
{
    homo_genes_.push_back(gene);
}

void Gene::add_neibor_fam_ids(const string &gene_id)
{
    neibor_fam_ids_.insert(gene_id);
}

// fills in the list with integers corresponding to the position of its homologous genes in gene list the gene

void Gene::homo_positions(const vector<Gene *> &genes, vector<int> &poss)
{
    for (int i = 0; i < genes.size(); i++)
    {
        if (is_homo_with(genes[i]))
        {
            poss.push_back(i);
        }
    }
}

// remap the this gene onto other gene

void Gene::remap_to(Gene *gene)
{
    is_tandem_ = true;
    is_tandem_represent_ = false;
    tandem_represent_ = gene->gene_id();
    is_remapped_ = true;
    gene->is_tandem_ = true;
    gene->is_tandem_represent_ = true;
    gene->tandem_represent_ = gene->gene_id();
    gene->is_remapped_ = false;
    gene->tandem_ids_.insert(gene_id_);

    set<string>::iterator it = tandem_ids_.begin();
    while (it != tandem_ids_.end())
    {
        gene->tandem_ids_.insert(*it);
        it++;
    }
}

void Gene::set_homos(set<string> &gene_ids)
{
    set<string>::iterator it = gene_ids.begin();
    while (it != gene_ids.end())
    {
        homo_ids_.insert(*it);
    }
}

bool Gene::is_homo_with(Gene *gene)
{
    set<string>::iterator it = homo_ids_.find(gene->gene_id());
    if (it != homo_ids_.end())
    {
        return true;
    }
    else
    {
        return false;
    }
}

void Gene::print()
{
    cerr << "Gene ID: " << gene_id_ << ",";
    cerr << "Chromosome ID: " << chr_id_ << ",";
    cerr << "Start: " << start_ << ",";
    cerr << "End: " << end_ << ",";
    cerr << "Strand: " << strand_ << ",";
    cerr << "Is tandem: " << is_tandem_ << ",";
    cerr << "Is tandem represent: " << is_tandem_represent_ << ",";
    cerr << "Tandem represent:" << tandem_represent_ << ",";
    cerr << "Is remammped: " << is_remapped_ << ",";
    cerr << "Is family gene: " << is_fam_gene_ << ",";
    cerr << "Num of homologous gene ID: " << homo_ids_.size() << ",";
    cerr << "Number of tandem: " << tandem_ids_.size() << ",";
    cerr << "Number of neibor family gene ID: " << neibor_fam_ids_.size() << endl;
}
