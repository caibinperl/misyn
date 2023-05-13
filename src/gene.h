#ifndef GENE_H
#define GENE_H

#include "base.h"

class Gene
{
private:
    string gene_id_;
    string spe_id_;
    string chr_id_;
    long start_;
    long end_;
    int strand_; // 1, 0
    bool is_tandem_;
    bool is_tandem_represent_;
    string tandem_represent_;
    bool is_remapped_;
    bool is_fam_gene_;
    set<string> homo_ids_;
    set<string> tandem_ids_;
    vector<Gene *> homo_genes_;
    set<string> neibor_fam_ids_;
    int index_;

public:
    Gene(const string &gene_id, const string &spe_id, const string &chr_id, long start, long end, int strand);
    ~Gene();

    void add_homo_ids(const string &gene_id);
    void add_homo_gene(Gene *gene);
    void add_neibor_fam_ids(const string &gene_id);
    void homo_positions(const vector<Gene *> &genes, vector<int> &poss);

    void invert_orientation()
    {
        strand_ = strand_ == 1 ? 0 : 1;
    }

    bool has_tandems()
    {
        return !tandem_ids_.empty();
    }

    const set<string> &tandem_ids() const
    {
        return tandem_ids_;
    }

    bool has_homos()
    {
        return !homo_ids_.empty();
    }

    void remap_to(Gene *gene);

    void set_homos(set<string> &gene_ids);
    bool is_homo_with(Gene *gene);

    void print();

    const string &gene_id() const
    {
        return gene_id_;
    }

    const string &spe_id() const
    {
        return spe_id_;
    }

    const string &chr_id() const
    {
        return chr_id_;
    }

    long start() const
    {
        return start_;
    }

    long end() const
    {
        return end_;
    }

    void set_is_fam_gene(bool is_fam_gene)
    {
        is_fam_gene_ = is_fam_gene;
    }

    bool is_fam_gene() const
    {
        return is_fam_gene_;
    }

    bool is_remapped() const
    {
        return is_remapped_;
    }

    int strand() const
    {
        return strand_;
    }

    int index() const
    {
        return index_;
    }

    void set_index(int index)
    {
        index_ = index;
    }
};

#endif /* GENE_H */
