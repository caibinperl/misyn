#ifndef CHR_H
#define CHR_H

#include "base.h"
#include "gene.h"

class Chr
{
private:
    string fam_gene_id_;
    vector<Gene *> chr_genes_;
    int chr_len_;
    int fam_gene_index_;

public:
    Chr(const string &gene_fam_id, const vector<Gene *> &chr_genes);
    ~Chr();

    const string &get_gene_id_by_index(int i);

    Gene *get_gene_by_id(const string &gene_id);
    int get_index_by_gene_id(const string &gene_id);
    Gene *get_gene_by_index(int i);

    int chr_len() const
    {
        return chr_len_;
    }

    int fam_gene_index() const
    {
        return fam_gene_index_;
    }

    const string &fam_gene_id() const
    {
        return fam_gene_id_;
    }

    const vector<Gene *> &chr_genes() const
    {
        return chr_genes_;
    }

    void index_genes();
};

#endif /* CHR_H */
