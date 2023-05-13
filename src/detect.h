
#ifndef MISYN_H
#define MISYN_H

#include <iostream>
#include <fstream>
#include <map>
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>
#include "base.h"
#include "config.h"
#include "gene.h"
#include "chr.h"
#include "cluster.h"
#include "syn.h"

using boost::bad_lexical_cast;
using boost::lexical_cast;

typedef map<string, vector<Gene *>> ChrGenesMap;

void load_chr2genes(const string &list_file, ChrGenesMap &fam2chr_genes);

void print_chr2genes(const ChrGenesMap &fam2chr_genes);

void load_fam(const ChrGenesMap &fam2chr_genes, vector<string> &fam_ids);

void generate_id2gene(const ChrGenesMap &fam2chr_genes, map<string, Gene *> &id2gene);

void print_id2gene(const map<string, Gene *> &id2gene);

void print_fam(const vector<string> &fam_ids);

void mark_fam_genes(const vector<string> &fam_ids, map<string, Gene *> &id2gene);

void print_fam_genes(const map<string, Gene *> &id2gene);

void load_homo(const string &homo_file, const vector<string> &fam_ids,
               map<string, Gene *> &id2gene);

void remap_tandems(int tandem_gap, ChrGenesMap &fam2chr_genes);

void remap_fam_tandems(int tandem_gap, ChrGenesMap &fam2chr_genes);

void generate_chr(const ChrGenesMap &fam2chr_genes, map<string, Chr *> &fam2chrs);

void print_chr(const map<string, Chr *> &fam2chrs);

void generate_fam(const map<string, Chr *> &fam2chrs, vector<string> &remap_fam_ids);

void build_matrix(Chr *chr_x, Chr *chr_y, Matrix matrixs[2]);

int point_count(const Matrix &matrix);

void print_matrix(const Matrix &matrix);

Cluster *detect_cluster(Chr *chr_x, Chr *chr_y, int orientation,
                        const Matrix &matrix, int gap, double q_value);

/* searches for anchorpoints to be inserted into a basesynteny that are
 located insided the gap size
 */
void enrich_clusters(const Matrix &matrix, int gap, double q_value, Cluster *cluster);

Syn *generate_syn(Cluster *clusters[2]);

void filter_syn(double e_value, int homo_count, Syn **syn);

Syn *detect_syn(Chr *chr_x, Chr *chr_y, int gap, double q_value, double e_value,
                int homo_count);

void detect(const vector<string> &remap_fam_ids, map<string, Chr *> &fam2chr,
            const Config &config, vector<Syn *> &syns);

void sort_syns(vector<Syn *> &syns);

void output_syn(const vector<Syn *> &syns, const string &out_file);

void output_homo(const vector<Syn *> &syns, const string &out_file);

void output_matrix(const vector<Syn *> &syns, const string &out_file);

void output_tandem(const vector<string> &remap_fam_ids,
                   map<string, Gene *> &id2gene,
                   const string &out_file);

void circos_label(const vector<string> &remap_fam_ids,
                  map<string, Gene *> &id2gene,
                  const string &out_file);

void circos_link(const vector<Syn *> &syns, map<string, Gene *> &id2gene,
                 const string &out_file);

void circos_karyotype(const vector<string> &remap_fam_ids,
                      map<string, Gene *> &id2gene,
                      const string &chr_len_dir, const string &out_file);

void delete_cluster(Cluster *cluster);

void delete_syn(Syn *syn);

#endif /* MISYN_H */
