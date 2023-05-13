#include "detect.h"

bool cmp_gene(Gene *g1, Gene *g2)
{
    if (g1->start() < g2->start())
        return true;
    else
        return false;
}

void load_chr2genes(const string &list_file, ChrGenesMap &fam2chr_genes)
{
    ifstream fin;
    fin.open(list_file.c_str());
    if (!fin.is_open())
    {
        cerr << "Could not open " << list_file << endl;
        exit(1);
    }
    // AT3G07770.TAIR10	AT3G08710.TAIR10	ath	Chr3	2645167	2646497	+
    string line;
    while (getline(fin, line))
    {
        boost::trim(line);
        vector<string> items;
        boost::split(items, line, boost::is_any_of("\t"));
        if (items.size() != 7)
        {
            cerr << "Error in " << list_file << endl;
            exit(1);
        }
        string fam_gene_id = items[0];
        string gene_id = items[1];
        string spe_id = items[2];
        string chr_id = items[3];
        string start_s = items[4];
        string end_s = items[5];
        string strand_s = items[6];
        if (!fam_gene_id.empty() && !gene_id.empty() && !spe_id.empty() &&
            !chr_id.empty() && !start_s.empty() && !end_s.empty() && !strand_s.empty())
        {
            int strand = 1;
            if (strand_s == "-")
            {
                strand = 0;
            }
            long start, end;
            try
            {
                start = lexical_cast<long>(start_s);
            }
            catch (bad_lexical_cast &)
            {
                cerr << "Wrong format for start in " << line << endl;
                exit(1);
            }
            try
            {
                end = lexical_cast<long>(end_s);
            }
            catch (bad_lexical_cast &)
            {
                cerr << "Wrong format for end in " << line << endl;
                exit(1);
            }
            // fill genes_info
            Gene *gene = new Gene(gene_id, spe_id, chr_id, start, end, strand);
            ChrGenesMap::iterator it = fam2chr_genes.find(fam_gene_id);
            if (it != fam2chr_genes.end())
            {
                it->second.push_back(gene);
            }
            else
            {
                vector<Gene *> genes = vector<Gene *>({gene});
                fam2chr_genes[fam_gene_id] = genes;
            }
        }
    }
    fin.close();
    // Sort genes in genes_info acoording to positon in GFF
    ChrGenesMap::iterator it;
    for (it = fam2chr_genes.begin(); it != fam2chr_genes.end(); it++)
    {
        vector<Gene *> &genes = it->second;
        sort(genes.begin(), genes.end(), cmp_gene);
    }
}

void print_chr2genes(const ChrGenesMap &fam2chr_genes)
{
    ChrGenesMap::const_iterator it;
    for (it = fam2chr_genes.begin(); it != fam2chr_genes.end(); it++)
    {
        cerr << "Chromosome: " << it->first << " --- Gene number: "
             << it->second.size() << endl;
    }
}

void load_fam(const ChrGenesMap &fam2chr_genes, vector<string> &fam_ids)
{
    ChrGenesMap::const_iterator it;
    for (it = fam2chr_genes.begin(); it != fam2chr_genes.end(); it++)
    {
        fam_ids.push_back(it->first);
    }
}

void print_fam(const vector<string> &fam_ids)
{
    cerr << "Family gene id:" << endl;
    vector<string>::const_iterator it;
    for (it = fam_ids.begin(); it != fam_ids.end(); it++)
    {
        cerr << *it << endl;
    }
}

void generate_id2gene(const ChrGenesMap &fam2chr_genes, map<string, Gene *> &id2gene)
{
    ChrGenesMap::const_iterator it1;
    for (it1 = fam2chr_genes.begin(); it1 != fam2chr_genes.end(); it1++)
    {
        const vector<Gene *> &genes = it1->second;
        vector<Gene *>::const_iterator it2;
        for (it2 = genes.begin(); it2 != genes.end(); it2++)
        {
            id2gene[(*it2)->gene_id()] = *it2;
        }
    }
}

void print_id2gene(const map<string, Gene *> &id2gene)
{
    cerr << "Id2gene: " << endl;
    map<string, Gene *>::const_iterator it;
    for (it = id2gene.begin(); it != id2gene.end(); it++)
    {
        it->second->print();
    }
}

void mark_fam_genes(const vector<string> &fam_ids, map<string, Gene *> &id2gene)
{
    vector<string>::const_iterator it1;
    for (it1 = fam_ids.begin(); it1 != fam_ids.end(); it1++)
    {
        map<string, Gene *>::iterator it2 = id2gene.find(*it1);
        if (it2 != id2gene.end())
        {
            it2->second->set_is_fam_gene(true);
        }
    }
}

void print_fam_genes(const map<string, Gene *> &id2gene)
{
    cerr << "Print family gene: " << endl;
    map<string, Gene *>::const_iterator it;
    for (it = id2gene.begin(); it != id2gene.end(); it++)
    {
        Gene *gene = it->second;
        if (gene->is_fam_gene())
            gene->print();
    }
}

void load_homo(const string &homo_file, const vector<string> &fam_ids,
               map<string, Gene *> &id2gene)
{
    ifstream fin;
    fin.open(homo_file.c_str());
    if (!fin.is_open())
    {
        cerr << "Could not open " << homo_file << endl;
        exit(1);
    }
    string line;
    while (getline(fin, line))
    {
        boost::trim(line);
        vector<string> items;
        boost::split(items, line, boost::is_any_of("\t"));
        if (items.size() == 3)
        {
            string gene_id_a = items[0];
            string gene_id_b = items[1];
            if (!gene_id_a.empty() && !gene_id_b.empty() && gene_id_a != gene_id_b)
            {
                vector<string>::const_iterator it_a = find(fam_ids.begin(),
                                                           fam_ids.end(), gene_id_a);
                vector<string>::const_iterator it_b = find(fam_ids.begin(),
                                                           fam_ids.end(), gene_id_b);

                if (it_a != fam_ids.end() || it_b != fam_ids.end())
                    continue;
                map<string, Gene *>::iterator it1 = id2gene.find(gene_id_a);
                map<string, Gene *>::iterator it2 = id2gene.find(gene_id_b);
                if (it1 != id2gene.end() && it2 != id2gene.end())
                {
                    it1->second->add_homo_ids(gene_id_b);
                    it2->second->add_homo_ids(gene_id_a);
                }
            }
        }
    }
    fin.close();
}

void remap_tandems(int tandem_gap, ChrGenesMap &fam2chr_genes)
{
    ChrGenesMap::iterator it;
    for (it = fam2chr_genes.begin(); it != fam2chr_genes.end(); it++)
    {
        vector<Gene *> &genes = it->second;
        for (int i = genes.size() - 1; i > 0; i--)
        {
            Gene *gene_i = genes[i];
            if (gene_i->has_homos())
            {
                int j = i - 1;
                bool found = false;
                while (!found && j >= 0 && j >= i - (tandem_gap + 1))
                {
                    Gene *gene_j = genes[j];
                    if (!gene_i->is_fam_gene() && !gene_j->is_fam_gene() && gene_i->is_homo_with(gene_j))
                    {
                        gene_i->remap_to(gene_j);
                        found = true;
                    }
                    j--;
                }
            }
        }
    }
}

void remap_fam_tandems(int tandem_gap, ChrGenesMap &fam2chr_genes)
{
    ChrGenesMap::iterator it;
    for (it = fam2chr_genes.begin(); it != fam2chr_genes.end(); it++)
    {
        vector<Gene *> &genes = it->second;
        for (int i = genes.size() - 1; i > 0; i--)
        {
            Gene *gene_i = genes[i];
            if (gene_i->is_fam_gene())
            {
                int j = i - 1;
                bool found = false;
                while (!found && j >= 0 && j >= i - (tandem_gap + 1))
                {
                    Gene *gene_j = genes[j];
                    if (gene_j->is_fam_gene())
                    {
                        gene_i->remap_to(gene_j);
                        found = true;
                    }
                    j--;
                }
            }
        }
    }
}

void generate_chr(const ChrGenesMap &fam2chr_genes, map<string, Chr *> &fam2chrs)
{
    ChrGenesMap::const_iterator it1;
    for (it1 = fam2chr_genes.begin(); it1 != fam2chr_genes.end(); it1++)
    {
        string fam_gene_id = it1->first;
        const vector<Gene *> &genes = it1->second;
        vector<Gene *> remap_genes;
        vector<Gene *>::const_iterator it2;
        for (it2 = genes.begin(); it2 != genes.end(); it2++)
        {
            Gene *gene = (*it2);
            if (!gene->is_remapped())
            {
                remap_genes.push_back(gene);
            }
        }
        sort(remap_genes.begin(), remap_genes.end(), cmp_gene);
        fam2chrs[fam_gene_id] = new Chr(fam_gene_id, remap_genes);
    }

    for (map<string, Chr *>::iterator it = fam2chrs.begin();
         it != fam2chrs.end(); it++)
    {
        Chr *chr = it->second;
        chr->index_genes();
    }
}

void print_chr(const map<string, Chr *> &fam2chrs)
{
    map<string, Chr *>::const_iterator it;
    for (it = fam2chrs.begin(); it != fam2chrs.end(); it++)
    {
        Chr *chr = it->second;
        cout << it->first << ": " << chr->chr_len() << endl;
    }
}

void generate_fam(const map<string, Chr *> &fam2chrs, vector<string> &remap_fam_ids)
{
    map<string, Chr *>::const_iterator it1;
    for (it1 = fam2chrs.begin(); it1 != fam2chrs.end(); it1++)
    {
        const vector<Gene *> &genes = it1->second->chr_genes();
        vector<Gene *>::const_iterator it2;
        for (it2 = genes.begin(); it2 != genes.end(); it2++)
        {
            Gene *gene = *it2;
            if (gene->is_fam_gene())
            {
                remap_fam_ids.push_back(gene->gene_id());
            }
        }
    }
}

void build_matrix(Chr *chr_x, Chr *chr_y, Matrix matrixs[2])
{
    //<x, <y1,y2,y3> >
    const vector<Gene *> &genes_x = chr_x->chr_genes();
    const vector<Gene *> &genes_y = chr_y->chr_genes();
    for (int x = 0; x < genes_x.size(); x++)
    {
        Gene *gene_x = genes_x[x];
        if (gene_x->has_homos())
        {
            vector<int> poss;
            gene_x->homo_positions(genes_y, poss);
            for (int i = 0; i < poss.size(); i++)
            {
                int y = poss[i];
                Gene *gene_y = genes_y[y];
                int orientation = 0;
                if (gene_x->strand() == gene_y->strand())
                {
                    orientation = 1;
                }
                Matrix::iterator it = matrixs[orientation].find(x);
                if (it != matrixs[orientation].end())
                {
                    vector<int> &temp = matrixs[orientation][x];
                    temp.push_back(y);
                }
                else
                {
                    vector<int> temp;
                    temp.push_back(y);
                    matrixs[orientation][x] = temp;
                }
            }
        }
    }
}

int point_count(const Matrix &matrix)
{
    int num = 0;
    Matrix::const_iterator it1;
    for (it1 = matrix.begin(); it1 != matrix.end(); it1++)
    {
        const vector<int> &ys = it1->second;
        num += ys.size();
    }
    return num;
}

void print_matrix(const Matrix &matrix)
{
    Matrix::const_iterator it1;
    for (it1 = matrix.begin(); it1 != matrix.end(); it1++)
    {
        int x = it1->first;
        const vector<int> &ys = it1->second;
        vector<int>::const_iterator it2;
        cout << "x:" << x << " -- y: ";
        for (it2 = ys.begin(); it2 != ys.end(); it2++)
        {
            cout << *it2 << ",";
        }
        cout << endl;
    }
}

Cluster *detect_cluster(Chr *chr_x, Chr *chr_y, int orientation,
                        const Matrix &matrix, int gap, double q_value)
{
    int fam_gene_index_x = chr_x->fam_gene_index();
    int fam_gene_index_y = chr_y->fam_gene_index();
    if (fam_gene_index_x == -1)
    {
        cerr << "Failed to find the coord of genes!" << endl;
        exit(1);
    }

    if (fam_gene_index_y == -1)
    {
        cerr << "Failed to find the coord of genes!" << endl;
        exit(1);
    }

    int sign = 1;
    if (orientation == 0)
    {
        sign = -1;
    }

    Cluster *cluster = new Cluster(chr_x, chr_y, orientation, matrix);

    // ref_x, ref_y start point every  search cirle
    int ref_x = fam_gene_index_x;
    int ref_y = fam_gene_index_y;
    cluster->add_homo_point(fam_gene_index_x, fam_gene_index_y);

    while (ref_x >= 0)
    {
        bool found = false;
        int closest_x = 0, closest_y = 0;
        double closest_distance = gap + 1;
        for (int x = ref_x + 1; x <= ref_x + gap; x++)
        {
            Matrix::const_iterator it1 = matrix.find(x);
            if (it1 != matrix.end())
            {
                const vector<int> &ys = it1->second;
                for (int y = ref_y + sign; abs(ref_y - y) <= gap; y += sign)
                {
                    vector<int>::const_iterator it2 = find(ys.begin(), ys.end(), y);
                    if (it2 != ys.end())
                    {
                        double distance = cluster->dpd(ref_x, ref_y, x, y);
                        if (distance < closest_distance)
                        {
                            closest_distance = distance;
                            closest_x = x;
                            closest_y = y;
                            found = true;
                        }
                    }
                }
            }
        }
        if (found)
        {
            cluster->add_homo_point(closest_x, closest_y);
            ref_x = closest_x;
            ref_y = closest_y;
        }
        else
        {
            ref_x = -1;
        }
    }

    ref_x = fam_gene_index_x;
    ref_y = fam_gene_index_y;

    while (ref_x >= 0)
    {
        bool found = false;
        int closest_x = 0, closest_y = 0;
        double closest_distance = gap + 1;
        for (int x = ref_x - 1; x >= ref_x - gap; x--)
        {
            Matrix::const_iterator it1 = matrix.find(x);
            if (it1 != matrix.end())
            {
                const vector<int> &ys = it1->second;
                for (int y = ref_y - sign; abs(ref_y - y) <= gap; y -= sign)
                {
                    vector<int>::const_iterator it2 = find(ys.begin(), ys.end(), y);
                    if (it2 != ys.end())
                    {
                        double distance = cluster->dpd(ref_x, ref_y, x, y);
                        if (distance < closest_distance)
                        {
                            closest_distance = distance;
                            closest_x = x;
                            closest_y = y;
                            found = true;
                        }
                    }
                }
            }
        }

        if (found)
        {
            cluster->add_reverse_homo_point(closest_x, closest_y);
            ref_x = closest_x;
            ref_y = closest_y;
        }
        else
        {
            ref_x = -1;
        }
    }

    if (cluster->get_homo_count() > 2)
    {
        cluster->set_fam_gene_coordinate(fam_gene_index_x, fam_gene_index_y);
        cluster->confidence_interval();
        return cluster;
    }
    else
    {
        delete_cluster(cluster);
        return nullptr;
    }
}

void enrich_clusters(const Matrix &matrix, int gap, double q_value, Cluster *cluster)
{
    if (cluster != nullptr)
    {
        cluster->confidence_interval();
        Matrix::const_iterator it1;
        for (it1 = matrix.begin(); it1 != matrix.end(); it1++)
        {
            int x = it1->first;
            const vector<int> &ys = it1->second;
            for (int i = ys.size() - 1; i >= 0; i--)
            {
                int y = ys[i];
                if (cluster->is_exist_point(x, y))
                {
                    continue;
                }
                double closest_distance = gap + 1;
                double distance = cluster->distance_to_point(x, y);
                if (distance < closest_distance && cluster->r_squared(x, y) >= q_value && cluster->in_interval(x, y))
                {
                    cluster->add_homo_point(x, y);
                    cluster->confidence_interval();
                }
            }
        }
    }
}

Syn *generate_syn(Cluster *clusters[2])
{
    Syn *syn = nullptr;
    if (clusters[0] != nullptr && clusters[1] != nullptr)
    {
        double prob0, prob1;
        prob0 = clusters[0]->prob();
        prob1 = clusters[1]->prob();
        if (prob1 <= prob0)
        {
            syn = new Syn(clusters[1]);
            delete_cluster(clusters[0]);
        }
        else
        {
            syn = new Syn(clusters[0]);
            delete_cluster(clusters[1]);
        }
    }
    else if (clusters[0] != nullptr && clusters[1] == nullptr)
    {
        syn = new Syn(clusters[0]);
    }
    else if (clusters[0] == nullptr && clusters[1] != nullptr)
    {
        syn = new Syn(clusters[1]);
    }
    return syn;
}

void filter_syn(double e_value, int homo_count, Syn **syn)
{
    if (*syn != nullptr)
    {
        double probe = (*syn)->prob();
        if (probe > e_value || (*syn)->get_homo_count() < homo_count)
        {
            delete_syn(*syn);
            *syn = nullptr;
        }
    }
}

Syn *detect_syn(Chr *chr_x, Chr *chr_y, int gap, double q_value, double e_value,
                int homo_count)
{
    Syn *syn = nullptr;
    Cluster *clusters[2];
    for (int i = 0; i < 2; i++)
    {
        clusters[i] = nullptr;
    }

    Matrix matrixs[2];

    build_matrix(chr_x, chr_y, matrixs);

    if (DEBUG)
    {
        cout << "Matrix 0 for chr_x: " << chr_x->fam_gene_id() << " and chr_y: "
             << chr_y->fam_gene_id() << endl;
        print_matrix(matrixs[0]);

        cout << "Matrix 1 for chr_x: " << chr_x->fam_gene_id() << " and chr_y: "
             << chr_y->fam_gene_id() << endl;
        print_matrix(matrixs[1]);
    }

    for (int i = 0; i < 2; i++)
    {
        int n = point_count(matrixs[i]);
        if (n < 3)
        {
            continue;
        }
        Cluster *cluster = detect_cluster(chr_x, chr_y, i, matrixs[i], gap, q_value);
        enrich_clusters(matrixs[i], gap, q_value, cluster);
        int area = chr_x->chr_len() * chr_y->chr_len();
        if (cluster != nullptr)
        {
            cluster->calculate_prob(area, n);
            clusters[i] = cluster;
        }
    }

    syn = generate_syn(clusters);
    filter_syn(e_value, homo_count, &syn);
    return syn;
}

void detect(const vector<string> &remap_fam_ids, map<string, Chr *> &fam2chr,
            const Config &config, vector<Syn *> &syns)
{
    for (int i = 0; i < remap_fam_ids.size() - 1; i++)
    {
        string gene_id_x = remap_fam_ids[i];

        for (int j = i + 1; j < remap_fam_ids.size(); j++)
        {
            string gene_id_y = remap_fam_ids[j];
            map<string, Chr *>::const_iterator it1 = fam2chr.find(gene_id_x);
            if (it1 == fam2chr.end())
            {
                cerr << "No corresponding list for " << gene_id_x << endl;
                exit(1);
            }
            map<string, Chr *>::const_iterator it2 = fam2chr.find(gene_id_y);
            if (it2 == fam2chr.end())
            {
                cerr << "No corresponding list for " << gene_id_y << endl;
                exit(1);
            }
            Chr *chr_x = fam2chr[gene_id_x];
            Chr *chr_y = fam2chr[gene_id_y];

            Syn *syn = detect_syn(chr_x, chr_y, config.gap(), config.q_value(),
                                  config.e_value(), config.homo_count());
            if (syn != nullptr)
            {
                syns.push_back(syn);
            }
        }
    }
}

static bool cmp_syn_(Syn *s1, Syn *s2)
{
    if (s1->prob() < s2->prob())
        return true;
    else
        return false;
}

void sort_syns(vector<Syn *> &syns)
{
    sort(syns.begin(), syns.end(), cmp_syn_);
}

void output_syn(const vector<Syn *> &syns, const string &out_file)
{
    ofstream fout(out_file.c_str());
    if (!fout.is_open())
    {
        cerr << "Cannot open " << out_file << " file for output" << endl;
    }
    fout.precision(2);
    fout << "syn_id\tgene_x_id\tgene_y_id\thomo_count\tprob" << endl;
    for (int i = 0; i < syns.size(); i++)
    {
        Syn *syn = syns[i];
        fout << i + 1 << "\t" << syn->fam_gene_x_id() << "\t" << syn->fam_gene_y_id()
             << "\t" << syn->get_homo_count() << "\t" << syn->prob() << endl;
    }
    fout.close();
}

void output_homo(const vector<Syn *> &syns, const string &out_file)
{
    ofstream fout(out_file.c_str());
    if (!fout.is_open())
    {
        cerr << "Cannot open " << out_file << "file for output" << endl;
    }
    fout << "syn_id\thomo_x_id\thomo_y_id\thomo_x_index\thomo_y_index" << endl;
    fout.precision(2);
    for (int i = 0; i < syns.size(); i++)
    {
        Syn *syn = syns[i];
        const vector<Homo *> &homos = syn->homos();
        for (int j = 0; j < homos.size(); j++)
        {
            Homo *homo = homos[j];
            Gene *gene_x = homo->gene_x();
            Gene *gene_y = homo->gene_y();
            fout << i + 1 << "\t" << gene_x->gene_id() << "\t" << gene_y->gene_id()
                 << "\t" << gene_x->index() << "\t" << gene_y->index() << endl;
        }
    }
    fout.close();
}

void output_matrix(const vector<Syn *> &syns, const string &out_file)
{
    ofstream fout(out_file.c_str());
    if (!fout.is_open())
    {
        cerr << "Cannot open " << out_file << "file for output" << endl;
    }
    fout << "syn_id\tx\ty" << endl;
    fout.precision(2);
    for (int i = 0; i < syns.size(); i++)
    {
        Syn *syn = syns[i];
        Matrix matrix = syn->cluster()->matrix();
        for (auto it = matrix.begin(); it != matrix.end(); it++)
        {
            int x = it->first;
            vector<int> ys = it->second;
            for (auto it2 = ys.begin(); it2 != ys.end(); it2++)
            {
                int y = *it2;
                fout << i + 1 << "\t" << x << "\t" << y << endl;
            }
        }
    }
    fout.close();
}

void output_tandem(const vector<string> &remap_fam_ids,
                   map<string, Gene *> &id2gene,
                   const string &out_file)
{

    ofstream fout(out_file.c_str());
    if (!fout.is_open())
    {
        cerr << "Cannot open " << out_file << "file for output" << endl;
    }
    fout.precision(2);
    for (int i = 0; i < remap_fam_ids.size(); i++)
    {
        string gene_id = remap_fam_ids[i];
        Gene *gene = id2gene[gene_id];
        if (gene->has_tandems())
        {
            fout << gene->gene_id();
            const set<string> &tandem_ids = gene->tandem_ids();
            set<string>::const_iterator it = tandem_ids.begin();
            while (it != tandem_ids.end())
            {
                fout << "\t" << *it;
                it++;
            }
            fout << endl;
        }
    }
    fout.close();
}

void circos_label(const vector<string> &remap_fam_ids,
                  map<string, Gene *> &id2gene,
                  const string &out_file)
{
    ofstream fout(out_file.c_str());
    if (!fout.is_open())
    {
        cerr << "Cannot open " << out_file << "file for output" << endl;
    }
    for (int i = 0; i < remap_fam_ids.size(); i++)
    {
        string gene_id = remap_fam_ids[i];
        Gene *gene = id2gene[gene_id];
        fout << gene->spe_id() << " " << gene->chr_id() << " " << gene->start()
             << " " << gene->end() << " " << gene->gene_id() << endl;
    }
    fout.close();
}

void circos_link(const vector<Syn *> &syns, map<string, Gene *> &id2gene,
                 const string &out_file)
{
    ofstream fout(out_file.c_str());
    if (!fout.is_open())
    {
        cerr << "Cannot open " << out_file << "file for output" << endl;
    }
    for (int i = 0; i < syns.size(); i++)
    {
        Syn *syn = syns[i];
        fout.precision(2);
        Gene *gene_x = id2gene[syn->fam_gene_x_id()];
        Gene *gene_y = id2gene[syn->fam_gene_y_id()];
        fout << gene_x->spe_id() << " " << gene_x->chr_id() << " "
             << gene_x->start() << " " << gene_x->end() << " "
             << gene_y->spe_id() << " " << gene_y->chr_id() << " "
             << gene_y->start() << " " << gene_y->end() << " color=red"
             << endl;
    }
    fout.close();
}

void circos_karyotype(const vector<string> &remap_fam_ids,
                      map<string, Gene *> &id2gene,
                      const string &chr_len_dir, const string &out_file)
{

    map<string, set<string>> spe2chr_ids;
    ofstream fout(out_file.c_str());
    if (!fout.is_open())
    {
        cerr << "Cannot open " << out_file << "file for output" << endl;
    }
    for (int i = 0; i < remap_fam_ids.size(); i++)
    {
        string gene_id = remap_fam_ids[i];
        Gene *gene = id2gene[gene_id];
        string spe_id = gene->spe_id();
        string chr_id = gene->chr_id();
        map<string, set<string>>::const_iterator it = spe2chr_ids.find(spe_id);
        if (it != spe2chr_ids.end())
        {
            spe2chr_ids[spe_id].insert(chr_id);
        }
        else
        {
            spe2chr_ids[spe_id] = set<string>{chr_id};
        }
    }

    int i = 0;

    for (map<string, set<string>>::iterator it = spe2chr_ids.begin();
         it != spe2chr_ids.end(); it++)
    {
        string spe_id = it->first;
        set<string> &chr_ids = it->second;
        i = i % 20 + 1;

        boost::filesystem::path chr_len_path(chr_len_dir);
        string chr_len_file_ = (chr_len_path / (spe_id + "_len.tsv")).string();
        ifstream fin;
        fin.open(chr_len_file_.c_str());
        if (!fin.is_open())
        {
            cerr << "Could not open " << chr_len_file_ << endl;
            exit(1);
        }
        string line;
        while (getline(fin, line))
        {
            boost::trim(line);
            vector<string> ss;
            boost::split(ss, line, boost::is_any_of("\t"));
            if (ss.size() == 2)
            {
                string chr_id = ss[0];
                string chr_len = ss[1];
                set<string>::iterator it = chr_ids.find(chr_id);
                if (it != chr_ids.end())
                {
                    fout << "chr - " << spe_id + chr_id << " " << chr_id << " 0 "
                         << chr_len << " chr" << i << endl;
                }
            }
        }
        fin.close();
        i += 1;
    }
    fout.close();
}

void delete_cluster(Cluster *cluster)
{
    vector<Homo *> homos = cluster->homos();
    for (int j = 0; j < homos.size(); j++)
    {
        delete homos[j];
    }
    delete cluster;
}

void delete_syn(Syn *syn)
{
    Cluster *cluster = syn->cluster();
    delete_cluster(cluster);
    delete syn;
}