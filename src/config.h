#ifndef CONFIG_H
#define CONFIG_H

#include <string>

using namespace std;

class Config
{
private:
    string list_file_;
    string homo_file_;
    string chr_len_dir_;
    string output_dir_;
    int gap_;
    int tandem_gap_;
    double q_value_;
    int homo_count_;
    double e_value_;

    string out_syn_file_;
    string out_homo_file_;
    string out_tandem_file_;
    string out_matrix_file_;

    string circos_label_file_;
    string circos_link_file_;
    string circos_karyotype_file_;

    int gap_sizes_[10];

public:
    Config();
    ~Config();

    void parse_arg(int argc, char *argv[]);
    void print_config();

    int gap() const
    {
        return gap_;
    }

    int homo_count() const
    {
        return homo_count_;
    }

    const string &homo_file() const
    {
        return homo_file_;
    }

    const string &list_file() const
    {
        return list_file_;
    }

    void set_list_file(const string &list_file)
    {
        list_file_ = list_file;
    }

    const string &out_syn_file() const
    {
        return out_syn_file_;
    }

    const string &out_homo_file() const
    {
        return out_homo_file_;
    }

    const string &out_matrix_file() const {
        return out_matrix_file_;
    }

    const string &out_tandem_file() const
    {
        return out_tandem_file_;
    }

    const string &circos_label_file() const
    {
        return circos_label_file_;
    }

    const string &circos_link_file() const
    {
        return circos_link_file_;
    }

    const string &circos_karyotype_file() const
    {
        return circos_karyotype_file_;
    }

    const string &output_dir() const
    {
        return output_dir_;
    }

    const string &chr_len_dir() const
    {
        return chr_len_dir_;
    }

    double e_value() const
    {
        return e_value_;
    }

    double q_value() const
    {
        return q_value_;
    }

    int tandem_gap() const
    {
        return tandem_gap_;
    }

    const int *gap_sizes() const
    {
        return gap_sizes_;
    }

private:
    void print_usag_();
    void check_arg_();
    void check_file_();
};

#endif /* CONFIG_H_ */
