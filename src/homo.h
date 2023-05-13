#ifndef HOMO_H
#define HOMO_H

#include "gene.h"

class Cluster;

class Homo
{
private:
    Gene *gene_x_;
    Gene *gene_y_;
    Cluster *cluster_;
    int x_, y_;
    bool nil_;

public:
    Homo(int x, int y, Cluster *cluster);
    ~Homo();

    void set_genes(Gene *gene_x, Gene *gene_y)
    {
        gene_x_ = gene_x;
        gene_y_ = gene_y;
    }

    Gene *gene_x() const
    {
        return gene_x_;
    }

    Gene *gene_y() const
    {
        return gene_y_;
    }

    int x() const
    {
        return x_;
    }

    int y() const
    {
        return y_;
    }

    void set_nil(bool nil)
    {
        nil_ = nil;
    }

    bool nil() const
    {
        return nil_;
    }
};

#endif /* HOMO_H_ */
