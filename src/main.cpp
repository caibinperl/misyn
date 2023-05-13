#include "detect.h"

int main(int argc, char *argv[])
{
    vector<string> fam_ids;
    ChrGenesMap fam2chr_genes;
    map<string, Gene *> id2gene;
    vector<string> remap_fam_ids;
    map<string, Chr *> fam2chrs;
    vector<Syn *> syns;

    Config config;
    config.parse_arg(argc, argv);
    config.print_config();

    cout << "loading chromosome list!" << endl;
    load_chr2genes(config.list_file(), fam2chr_genes);
    if (DEBUG)
        print_chr2genes(fam2chr_genes);

    cout << "loading family gene id!" << endl;
    load_fam(fam2chr_genes, fam_ids);
    if (DEBUG)
        print_fam(fam_ids);

    cout << "Generating id2gene" << endl;
    generate_id2gene(fam2chr_genes, id2gene);
    if (DEBUG)
        print_id2gene(id2gene);

    cout << "Mark famaliy gene in id2gene" << endl;
    mark_fam_genes(fam_ids, id2gene);
    if (DEBUG)
        print_fam_genes(id2gene);

    cout << "Add homologous genes" << endl;
    load_homo(config.homo_file(), fam_ids, id2gene);
    if (DEBUG)
        print_id2gene(id2gene);

    cout << "Remap tandem genes" << endl;
    remap_tandems(config.tandem_gap(), fam2chr_genes);
    if (DEBUG)
        print_id2gene(id2gene);

    cout << "Remap tandem family genes" << endl;
    remap_fam_tandems(config.tandem_gap(), fam2chr_genes);
    if (DEBUG)
        print_id2gene(id2gene);

    cout << "Generate chromoseom" << endl;
    generate_chr(fam2chr_genes, fam2chrs);
    if (DEBUG)
        print_chr(fam2chrs);

    cout << "Generate family gene id" << endl;
    generate_fam(fam2chrs, remap_fam_ids);
    if (DEBUG)
        print_fam(remap_fam_ids);

    cout << "Detect cluster" << endl;
    detect(remap_fam_ids, fam2chrs, config, syns);


    cout << "Output results" << endl;
    if (!syns.empty())
    {
        sort_syns(syns);

        output_syn(syns, config.out_syn_file());

        output_homo(syns, config.out_homo_file());

        output_matrix(syns, config.out_matrix_file());

        output_tandem(remap_fam_ids, id2gene, config.out_tandem_file());

        circos_label(remap_fam_ids, id2gene, config.circos_label_file());

        circos_link(syns, id2gene, config.circos_link_file());

        circos_karyotype(remap_fam_ids, id2gene, config.chr_len_dir(),
                         config.circos_karyotype_file());
    }

    for (ChrGenesMap::iterator it = fam2chr_genes.begin(); it != fam2chr_genes.end(); it++)
    {
        vector<Gene *> &genes = it->second;
        for (int i = 0; i < genes.size(); i++)
        {
            delete genes[i];
        }
    }

    for (map<string, Chr *>::iterator it = fam2chrs.begin(); it != fam2chrs.end(); it++)
    {
        Chr *chr = it->second;
        delete chr;
    }

    for (int i = 0; i < syns.size(); i++)
    {
        delete_syn(syns[i]);
    }

    return 0;
}
