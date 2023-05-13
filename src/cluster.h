#ifndef CLUSTER_H
#define CLUSTER_H

#include "base.h"
#include "gene.h"
#include "chr.h"
#include "homo.h"

class Cluster
{
private:
    /*
    -----start_x_-------end_x_-----
    -----start_y_-------end_y_-----
    -----start_y_hat_-------end_y_hat_-----
     */
    int start_x_, end_x_, start_y_, end_y_;
    double start_y_hat_, end_y_hat_;

    double a_, b_;
    double avg_x_, var_x_, mse_;

    double prob_;
    int orientation_;
    string fam_gene_x_id_, fam_gene_y_id_;
    int fam_gene_x_index_, fam_gene_y_index_;

    vector<Gene *> genes_x_;
    vector<Gene *> genes_y_;
    vector<Homo *> homos_;
    Matrix matrix_;

public:
    Cluster(Chr *chr_x, Chr *chr_y, int orientation, const Matrix &matrix);
    ~Cluster();

    void add_homo_point(int x, int y);

    void add_reverse_homo_point(int x, int y);

    void sort_homo_points();

    double dpd(int x1, int y1, int x2, int y2);

    // Returns for a given x coordinate the upper and lower y coordinates of the confidence interval
    void interval_bounds(int x, double ys[2]);

    /*
     *Returns a boolean saying the specified x and y coordinate is
     *located inside the confidence interval of the cluster
     */
    bool in_interval(int x, int y);

    double distance_to_point(int x, int y);

    /*Calculates the squared Pearson value for all x and y coordinates of all homologys,
     *including the specified values.
     */
    double r_squared(int X, int Y);

    bool is_exist_point(int x, int y);

    void set_fam_gene_coordinate(int fam_gene_x_index, int fam_gene_y_index);

    int get_lowest_x();

    int get_highest_x();

    int get_lowest_y();

    int get_highest_y();

    void set_bounds();

    /* regression coefficient
     *avg_x_, a_, b_
     */
    void regression();

    /*
     *when b_ == 0, n>=3, update a_, b_, avg_x_, var_x_, mse_, start_x_, end_x_,
     * start_y_hat_, end_y_hat_
     */
    void confidence_interval();

    /*
     * Calculates the probability to be generated by chance.
     * the area and number of points in the GHM needs to be passed as arguments
     */
    void calculate_prob(int area, int point_count);

    bool overlap_interval(Cluster *cluster);

    int overlap_coordinates(Cluster *cluster);

    int get_homo_count();

    void print();

    int orientation() const
    {
        return orientation_;
    }

    double avg_x() const
    {
        return avg_x_;
    }

    const string &fam_gene_x_id() const
    {
        return fam_gene_x_id_;
    }

    int fam_gene_x_index() const
    {
        return fam_gene_x_index_;
    }

    const string &fam_gene_y_id() const
    {
        return fam_gene_y_id_;
    }

    int fam_gene_y_index() const
    {
        return fam_gene_y_index_;
    }

    const vector<Gene *> &genes_x() const
    {
        return genes_x_;
    }

    const vector<Gene *> &genes_y() const
    {
        return genes_y_;
    }

    const vector<Homo *> &homos() const
    {
        return homos_;
    }

    double a() const
    {
        return a_;
    }

    double prob() const
    {
        return prob_;
    }

    double b() const
    {
        return b_;
    }

    double var_x() const
    {
        return var_x_;
    }

    double mes() const
    {
        return mse_;
    }

    int end_x() const
    {
        return end_x_;
    }

    int end_y() const
    {
        return end_y_;
    }

    int end_y_hat() const
    {
        return end_y_hat_;
    }

    int start_x() const
    {
        return start_x_;
    }

    int start_y() const
    {
        return start_y_;
    }

    int start_y_hat() const
    {
        return start_y_hat_;
    }

    const Matrix &matrix() const
    {
        return matrix_;
    }
};

#endif