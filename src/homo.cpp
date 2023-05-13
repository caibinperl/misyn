#include "homo.h"
#include "cluster.h"

Homo::Homo(int x, int y, Cluster *cluster)
{
    cluster_ = cluster;
    x_ = x;
    y_ = y;
    gene_x_ = nullptr;
    gene_y_ = nullptr;
    nil_ = false;
}

Homo::~Homo()
{
}
