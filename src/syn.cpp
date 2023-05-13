#include "syn.h"

Syn::Syn(Cluster *cluster)
{
    cluster_ = cluster;
    genes_x_ = cluster->genes_x();
    genes_y_ = cluster->genes_y();
    homos_ = cluster->homos();
    prob_ = cluster->prob();
    orientation_ = cluster->orientation();
    fam_gene_x_id_ = cluster->fam_gene_x_id();
    fam_gene_y_id_ = cluster->fam_gene_y_id();
    fam_gene_x_index_ = -1;
    fam_gene_y_index_ = -1;

    vector<Homo *>::const_iterator it = homos_.begin();
    while (it != homos_.end())
    {
        Homo *homo = *it;
        int x = homo->x();
        int y = homo->y();
        homo->set_genes(genes_x_[x], genes_y_[y]);
        it++;
    }
}

Syn::~Syn()
{
}
