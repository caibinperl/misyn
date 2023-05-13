#ifndef SYN_H
#define SYN_H

#include "base.h"
#include "gene.h"
#include "homo.h"
#include "cluster.h"

class Syn
{
private:
    Cluster *cluster_;
    vector<Gene *> genes_x_;
    vector<Gene *> genes_y_;
    vector<Homo *> homos_;
    double prob_;
    int orientation_;
    string fam_gene_x_id_, fam_gene_y_id_;
    int fam_gene_x_index_, fam_gene_y_index_;

public:
    Syn(Cluster *cluster);
    ~Syn();

    int get_homo_count()
    {
        return homos_.size();
    }

    Cluster *cluster() const
    {
        return cluster_;
    }

    double prob() const
    {
        return prob_;
    }

    const string &fam_gene_x_id() const
    {
        return fam_gene_x_id_;
    }

    const string &fam_gene_y_id() const
    {
        return fam_gene_y_id_;
    }

    const vector<Homo *> &homos() const
    {
        return homos_;
    }
};

#endif /* SYN_H_ */
