#include <algorithm>
#include <cmath>
#include <gsl/gsl_cdf.h>
#include "homo.h"
#include "cluster.h"

Cluster::Cluster(Chr *chr_x, Chr *chr_y, int orientation, const Matrix &matrix)
{
    genes_x_ = chr_x->chr_genes();
    genes_y_ = chr_y->chr_genes();
    orientation_ = orientation;
    fam_gene_x_id_ = chr_x->fam_gene_id();
    fam_gene_y_id_ = chr_y->fam_gene_id();
    fam_gene_x_index_ = chr_x->fam_gene_index();
    fam_gene_y_index_ = chr_y->fam_gene_index();

    start_x_ = 0;
    end_x_ = 0;
    start_y_ = 0;
    end_y_ = 0;
    start_y_hat_ = 0.0;
    end_y_hat_ = 0.0;

    a_ = 0.0;
    b_ = 0.0;
    avg_x_ = 0.0;
    var_x_ = 0.0;
    mse_ = 0.0;
    prob_ = 0.0;

    matrix_ = matrix;
}

Cluster::~Cluster()
{
}

void Cluster::add_homo_point(int x, int y)
{
    homos_.push_back(new Homo(x, y, this));
    b_ = 0;
}

// add into the front position in list

void Cluster::add_reverse_homo_point(int x, int y)
{
    homos_.insert(homos_.begin(), new Homo(x, y, this));
    b_ = 0;
}

static bool cmp_homo_(Homo *h1, Homo *h2)
{
    if (h1->x() < h2->x())
        return true;
    else
        return false;
}

void Cluster::sort_homo_points()
{
    sort(homos_.begin(), homos_.end(), cmp_homo_);
}

double Cluster::dpd(int start_x, int start_y, int end_x, int end_y)
{
    if (abs(start_x - end_x) >= abs(start_y - end_y))
    {
        return (2 * abs(start_x - end_x) - abs(start_y - end_y));
    }
    else
    {
        return (2 * abs(start_y - end_y) - abs(start_x - end_x));
    }
}

bool Cluster::in_interval(int x, int y)
{
    double ys[2];
    ys[0] = -1;
    ys[1] = -1;
    interval_bounds(x, ys);
    if (ys[0] == -1 || ys[1] == -1)
    {
        return false;
    }
    if (y >= ys[0] && y <= ys[1])
    {
        return true;
    }
    else
    {
        return false;
    }
}

void Cluster::interval_bounds(int x, double ys[2])
{
    confidence_interval();
    int n = homos_.size();
    if (n > 2)
    {
        double t = gsl_cdf_tdist_Qinv(0.005, (double)(homos_.size() - 2));
        t = abs(t);
        double sdev_y_est = sqrt(mse_) * sqrt(1.0 / n + pow(x - avg_x_, 2) / var_x_);
        ys[0] = (a_ + b_ * x) - (t * sdev_y_est);
        ys[1] = (a_ + b_ * x) + (t * sdev_y_est);
    }
}

// Calculates the distance from the basecluster to the specified x and y coordinate

double Cluster::distance_to_point(int x, int y)
{
    if (x < start_x_ && y < start_y_hat_)
    {
        return dpd(x, y, start_x_, start_y_hat_);
    }
    else if (x > end_x_ && y > end_y_hat_)
    {
        return dpd(x, y, end_x_, end_y_hat_);
    }
    else
    {
        return 0;
    }
}

double Cluster::r_squared(int X, int Y)
{
    double sum_x = 0, sum_y = 0, sum_xy = 0, sum_x2 = 0, sum_y2 = 0;
    int n = 0;
    for (int i = 0; i < homos_.size(); i++)
    {
        int x = homos_[i]->x();
        int y = homos_[i]->y();
        sum_x += x;
        sum_y += y;
        sum_xy += x * y;
        sum_x2 += x * x;
        sum_y2 += y * y;
        n++;
    }
    // if X and Y are given as arguments, count them too
    if (X != 0 || Y != 0)
    {
        sum_x += X;
        sum_y += Y;
        sum_xy += X * Y;
        sum_x2 += X * X;
        sum_y2 += Y * Y;
        n++;
    }

    double denominator = ((sum_x2 - sum_x * sum_x / n) * (sum_y2 - sum_y * sum_y / n));
    if (denominator == 0)
    {
        return -1;
    }
    else
    {
        double r = (sum_xy - (sum_x * sum_y / n)) / sqrt(denominator);
        return (r * r);
    }
}

/*
 * n!/x!(n-x)! p^x (1-p)^(n-x)
 * x = 1, np(1-p)^(n-1)
 */
void Cluster::calculate_prob(int area, int point_count)
{
    double n;
    double p = (double)point_count / area;
    double prob_i;
    double prob = 1;
    double distance;
    sort_homo_points();
    int fam_gene_index = -1;

    for (int i = 0; i < homos_.size(); i++)
    {
        if (homos_[i]->x() == fam_gene_x_index_)
        {
            fam_gene_index = i;
        }
    }

    if (fam_gene_index == -1)
    {
        cerr << "No gene index was found for " << fam_gene_x_id_ << endl;
        exit(1);
    }

    for (int i = fam_gene_index; i > 0; i--)
    {
        if (!homos_[i - 1]->nil() && !homos_[i]->nil())
        {
            distance = dpd(homos_[i - 1]->x(), homos_[i - 1]->y(), homos_[i]->x(),
                           homos_[i]->y());
            if (distance == 0)
            {
                homos_[i - 1]->set_nil(true);
            }
            else
            {
                n = ceil(distance * distance / 2.0);
                prob_i = (double)(n * p * pow(1 - p, n - 1));
                prob *= prob_i;
            }
        }
    }

    for (int i = fam_gene_index; i < homos_.size() - 1; i++)
    {
        if (!homos_[i + 1]->nil() && !homos_[i]->nil())
        {
            distance = dpd(homos_[i + 1]->x(), homos_[i + 1]->y(), homos_[i]->x(),
                           homos_[i]->y());
            if (distance == 0)
            {
                homos_[i + 1]->set_nil(true);
            }
            else
            {
                n = ceil((double)distance * distance / 2);
                prob_i = (double)(n * p * pow(1 - p, n - 1));
                prob *= prob_i;
            }
        }
    }

    vector<Homo *>::iterator it = homos_.begin();
    while (it != homos_.end())
    {
        if ((*it)->nil())
            it = homos_.erase(it);
        else
            it++;
    }

    prob_ = prob;
}

void Cluster::confidence_interval()
{
    if (b_ == 0.0)
    {
        regression();
        var_x_ = 0;
        mse_ = 0;
        for (int i = 0; i < homos_.size(); i++)
        {
            int x = homos_[i]->x();
            int y = homos_[i]->y();
            double y_hat = a_ + b_ * x;
            mse_ += (y - y_hat) * (y - y_hat);
            var_x_ += (x - avg_x_) * (x - avg_x_);
        }
        mse_ = mse_ / (homos_.size() - 2);
        start_x_ = get_lowest_x();
        end_x_ = get_highest_x();
        start_y_hat_ = a_ + start_x_ * b_;
        end_y_hat_ = a_ + end_x_ * b_;
    }
}

void Cluster::regression()
{
    if (homos_.size() > 2)
    {
        double sum_x = 0, sum_y = 0, sum_xy = 0, sum_x2 = 0;
        for (int i = 0; i < homos_.size(); i++)
        {
            double x = homos_[i]->x();
            double y = homos_[i]->y();
            sum_x += x;
            sum_y += y;
            sum_xy += x * y;
            sum_x2 += x * x;
        }
        avg_x_ = sum_x / homos_.size();
        double sxy = homos_.size() * sum_xy - sum_x * sum_y;
        double sxx = homos_.size() * sum_x2 - sum_x * sum_x;
        b_ = sxy / sxx;
        a_ = sum_y / homos_.size() - b_ * sum_x / homos_.size();
    }
}

bool Cluster::is_exist_point(int x, int y)
{
    vector<Homo *>::iterator it = homos_.begin();
    while (it != homos_.end())
    {
        if ((*it)->x() == x || (*it)->y() == y)
        {
            return true;
        }
        it++;
    }
    return false;
}

// Set

void Cluster::set_fam_gene_coordinate(int fam_gene_x_index, int fam_gene_y_index)
{
    fam_gene_x_index_ = fam_gene_x_index;
    fam_gene_y_index_ = fam_gene_y_index;
}

void Cluster::set_bounds()
{
    start_x_ = get_lowest_x();
    end_x_ = get_highest_x();
    start_y_ = get_lowest_y();
    end_y_ = get_highest_y();
}

int Cluster::get_lowest_x()
{
    if (homos_.size() > 0)
    {
        int lowest_x = homos_[0]->x();
        for (int i = 1; i < homos_.size(); i++)
        {
            int x = homos_[i]->x();
            if (x < lowest_x)
                lowest_x = x;
        }
        return lowest_x;
    }
    else
    {
        return -1;
    }
}

int Cluster::get_lowest_y()
{
    if (homos_.size() > 0)
    {
        int lowest_y = homos_[0]->y();
        for (int i = 1; i < homos_.size(); i++)
        {
            int y = homos_[i]->y();
            if (y < lowest_y)
            {
                lowest_y = y;
            }
        }
        return lowest_y;
    }
    else
    {
        return -1;
    }
}

int Cluster::get_highest_x()
{
    if (homos_.size() > 0)
    {
        int highest_x = homos_[0]->x();
        for (int i = 1; i < homos_.size(); i++)
        {
            int x = homos_[i]->x();
            if (x > highest_x)
            {
                highest_x = x;
            }
        }
        return highest_x;
    }
    else
    {
        return -1;
    }
}

int Cluster::get_highest_y()
{
    if (homos_.size() > 0)
    {
        int highest_y = homos_[0]->y();
        for (int i = 1; i < homos_.size(); i++)
        {
            int y = homos_[i]->y();
            if (y > highest_y)
            {
                highest_y = y;
            }
        }
        return highest_y;
    }
    else
    {
        return -1;
    }
}

int Cluster::get_homo_count()
{
    return homos_.size();
}

void Cluster::print()
{
    cout << "Num of homo: " << homos_.size() << ", Start x: " << start_x_
         << ", End x: " << end_x_ << ", Start y: " << start_y_ << ", End y: " << end_y_ << endl;
}
