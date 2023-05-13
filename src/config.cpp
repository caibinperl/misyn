#include <stdlib.h>
#include <vector>
#include <iostream>
#include <cmath>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>

#include "config.h"

using boost::bad_lexical_cast;
using boost::lexical_cast;

Config::Config()
{
    list_file_ = "";
    homo_file_ = "";
    output_dir_ = "";
    gap_ = 15;
    tandem_gap_ = 2;
    q_value_ = 0.8;
    homo_count_ = 3;
    e_value_ = 0.01;

    out_syn_file_ = "";
    out_homo_file_ = "";
    out_tandem_file_ = "";
    out_matrix_file_ = "";
    circos_label_file_ = "";
    circos_link_file_ = "";
    circos_karyotype_file_ = "";
}

Config::~Config(){};

void Config::print_usag_()
{
    const string USAGE = "misyn [option] args\n"
                         "-l Segment list file\n"
                         "-s Homology gene file\n"
                         "-o Output directory\n"
                         "-c Directory of chromosomes' length  for circos visualization\n"
                         "-e Expected Value\n";

    cout << USAGE << endl;
    exit(1);
}

void Config::parse_arg(int argc, char *argv[])
{
    int c;
    while ((c = getopt(argc, argv, "l:s:o:c:e:h")) != -1)
    {
        switch (c)
        {
        case 'l':
            list_file_ = optarg;
            break;
        case 's':
            homo_file_ = optarg;
            break;
        case 'o':
            output_dir_ = optarg;
            break;
        case 'c':
            chr_len_dir_ = optarg;
            break;
        case 'e':
            try
            {
                e_value_ = lexical_cast<double>(optarg);
            }
            catch (bad_lexical_cast &)
            {
                print_usag_();
            }
            break;
        case 'h':
            print_usag_();
            break;
        case '?':
            print_usag_();
            break;
        default:
            print_usag_();
        }
    }
    check_arg_();
    check_file_();
}

void Config::check_file_()
{
    boost::filesystem::path list_path(list_file_);
    boost::filesystem::path homo_path(homo_file_);
    boost::filesystem::path output_path(output_dir_);
    boost::filesystem::path chr_len_path(chr_len_dir_);
    if (!exists(list_path))
    {
        cerr << "Genome list file doesn't exist" << endl;
    }
    if (!exists(homo_path))
    {
        cerr << "Homology file doesn't exist" << endl;
    }

    if (!exists(chr_len_path))
    {
        cerr << "Chromosome length directory doesn't exist" << endl;
    }

    if (!exists(output_path))
    {
        cerr << "Output directory doesn't exist" << endl;
    }
    else
    {
        out_syn_file_ = (output_path / "syn_gene.tsv").string();
        out_homo_file_ = (output_path / "syn_homo.tsv").string();
        out_tandem_file_ = (output_path / "syn_tandem.txt").string();
        out_matrix_file_ = (output_path / "syn_matrix.tsv").string();
        circos_label_file_ = (output_path / "circos_label.txt").string();
        circos_link_file_ = (output_path / "circos_link.txt").string();
        circos_karyotype_file_ = (output_path / "circos_karyotype.txt").string();
    }
}

void Config::check_arg_()
{
    bool check = true;
    if (list_file_.empty())
    {
        cerr << "Genome list file was not set!" << endl;
        check = false;
    }

    if (homo_file_.empty())
    {
        cerr << "Homologous gene file was not set!" << endl;
        check = false;
    }

    if (output_dir_.empty())
    {
        cerr << "Output directory was not set!" << endl;
        check = false;
    }

    if (homo_count_ <= 0)
    {
        cerr << "Homology points size must > 0 !" << endl;
        check = false;
    }

    if (e_value_ <= 0)
    {
        cerr << "Expected Value: must > 0 !" << endl;
        check = false;
    }

    if (gap_ <= 0)
    {
        cerr << "Gap size must > 0 !" << endl;
        check = false;
    }

    if (tandem_gap_ <= 0)
    {
        cerr << "Tandem gap size must > 0 !" << endl;
        check = false;
    }

    if (!check)
    {
        exit(1);
    }
}

void Config::print_config()
{
    cout << "Genome list file: " << list_file_ << endl;
    cout << "Homology gene file: " << homo_file_ << endl;
    cout << "Output directory: " << output_dir_ << endl;
    cout << "Gap size: " << gap_ << endl;
    cout << "Tandem gap: " << tandem_gap_ << endl;
    cout << "Q value: " << q_value_ << endl;
    cout << "Homology pair cout: " << homo_count_ << endl;
    cout << "Probality cutoff: " << e_value_ << endl;
    cout << "Output tandem file: " << out_tandem_file_ << endl;
    cout << "Output synteny file: " << out_syn_file_ << endl;
    cout << "Output homology file: " << out_homo_file_ << endl;
    cout << "Output matrix file: " << out_matrix_file_ << endl;
}
