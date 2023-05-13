#include <iostream>
#include <cstdlib>
#include "chr.h"

Chr::Chr(const string &fam_gene_id, const vector<Gene *> &chr_genes)
{
    fam_gene_id_ = fam_gene_id;
    chr_genes_ = chr_genes;
    chr_len_ = chr_genes.size();
    fam_gene_index_ = get_index_by_gene_id(fam_gene_id);
}

Chr::~Chr()
{
}

// according to coordince

const string &Chr::get_gene_id_by_index(int i)
{
    if (i >= 0 && i < chr_genes_.size())
    {
        return chr_genes_[i]->gene_id();
    }
    return NULL;
}

Gene *Chr::get_gene_by_id(const string &gene_id)
{
    for (int i = 0; i < chr_genes_.size(); i++)
    {
        Gene *gene = chr_genes_[i];
        if (gene->gene_id() == gene_id)
        {
            return gene;
        }
    }
    return nullptr;
}

int Chr::get_index_by_gene_id(const string &gene_id)
{
    for (int i = 0; i < chr_genes_.size(); i++)
    {
        if (chr_genes_[i]->gene_id() == gene_id)
        {
            return i;
        }
    }
    return -1;
}

Gene *Chr::get_gene_by_index(int i)
{
    if (i >= 0 && i < chr_genes_.size())
    {
        return chr_genes_[i];
    }
    else
    {
        return nullptr;
    }
}

void Chr::index_genes()
{
    for (int i = 0; i < chr_genes_.size(); i++)
    {
        chr_genes_[i]->set_index(i);
    }
}
