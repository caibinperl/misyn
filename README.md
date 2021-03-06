# MiSyn: A tool for detection of synteny between two genes

MiSyn is designed to detect synteny between two genes.

## Algorithm
### Searching strategy
To study the relationship between two genes (eg. g1 and g2), their neighboring genes are considered. 
For example, the genomic fragment X contains gene g1 and its neighboring genes; the genomic fragment Y contains g2 and its neighboring genes. 
To survey the colinearity between g1 and g2, MiSyn uses the homologous relationship between genes on fragments X and Y.
The algorithm takes the gene pair (gene g1 and g2) as a starting point to search neighboring homologous pairs.
In the process of searching, the software always determines the best adjacent homologous gene pair which remains the best colinearity on the genomic fragments.
To determine the best adjacent gene pair, we adopt the distance (di) function that was introduced in ADhoRe[^1]. 
Consider one homologous pair with coordinate ($x_i$, $y_i$) and another pair with ($x_{i+1}$, $y_{i+1}$). 
The $d_i$ between two adjacent homology gene pairs ($x_i$, $y_i$) and ($x_{i+1}$, $y_{i+1}$) is given by:
$$2max(|x_{i+1} - x_i|, |y_{i+1} - y_i|) - min(|x_{i+1} - x_i|, |y_{i+1} - y_i|)$$

The distance controls the extent of the colinearity for a series of homologous pairs between two fragments.
In this way, the homologous pairs that are located in minimum distance are recorded and then taken as the new point for searching the adjacent homologous pair.
Finally, all the fitted pairs are stored into one cluster that contains gene pairs. 
The cluster represents the conserved region between two genomic fragments that have evolution relationship. 
The searching process is mainly controlled by two parameters: 
the maximum distance, $max(|x_{i+1} - x_i|, |y_{i+1} - y_i|)$, between two adjacent genes in each fragment and the actual value of $d_i$ between two adjacent homologous pairs.

### Statistical validation of conserved genomic fragment
To discard a negative cluster that is likely generated by chance, a statistical assessment was developed. 
Consider two fragments with number of $m$ and $n$ genes, respectively, and the number of homologous gene pairs represented by $l$.
The probability of finding a homologous pair in the two fragments is given by:

$$a = \frac{l}{m \times n}$$

Next, considering the detected cluster containing homologous pairs and a searching range of ${d_i}^2/2$ gene pairs from two fragments, within the context of a binomial distribution, the gene pair is observed by chance with the probability:

$$p_i = \frac{{d_i}^2}{2} a(1-a)^{\frac{{d_i}^2}{2} - 1}$$

Thus, when MicroSyn starts searching from a specific gene pair, it detects homologous pairs in an area of $\sum_{1}^{k}{\frac{d_i^2}{2}}$, where $k$ equals the number of homologous gene pairs in the cluster. Thus, the final probability to find the cluster containing homologous gene pairs by chance is given by:

$$p_c = \prod_{1}^{k} p_i$$

In other words, the value of $p_c$ is determined by both the number of homologous gene pairs in a given detected cluster and the total number of genes in two genomic fragments. 
So, if $p_c$ of a cluster exceeds a threshold value, the cluster is considered to be negative and should be discarded.

## Installation

### Prerequisites

* Python 3.10.8
* build-essential
* libgsl-dev 
* libboost-all-dev
* blast

### compile isomir
Move to c++ source code directory and compile using:

```bash
cd misyn/src
make
```

The compiled binary is called `misyn`

## Usage
We demonstrate the utility of MiSyn by applying it to exploration of HSP90 genes in _Arabidopsis_ and _Populus_.

1. Download test data from Phytozome <https://phytozome-next.jgi.doe.gov/>, including:
    * the genome chromosome sequence for _Arabidopsis_ and _Populus_
    * GFF3 for _Arabidopsis_ and _Populus_

    Store HSP90 gene ids of _Arabidopsis_ and _Populus_ in gene list file.

2. Edit metadata.tsv

3. Run `misyn.py`
    ```bash
    python misyn.py metadata.tsv output
    ```

## References

[^1] Simillion C, Vandepoele K, Saeys Y, Van de Peer Y: Building genomic profiles for uncovering segmental homology in the twilight zone. Genome Res 2004, 14(6):1095-1106.