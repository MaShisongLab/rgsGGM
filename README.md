# rgsGGM

**Gene co-expression network analysis based on the graphical Gaussian model.**

This is a MATLAB algorithm for gene co-expression network analysis based on the graphical Gaussian model (GGM). The algorithm is modified from a previously published method ([Ma *et al*, 2007](https://github.com/MaShisongLab/rgsGGM#References)). It takes an gene expression matrix as input and uses a procedure consisting of ~ 20,000 iterations to calculate partial correlation coefficients (pcors) between gene pairs. In each iteration, 2,000 genes are randomly selected and pcors between these genes are calculated. After all iterations, every gene pair is selected in multiple iterations with multiple pcors calculated, and the pcor with lowest absolution value is chosen as the gene pair's final pcor. All the gene pairs with pcors larger or equal to a selected cutoff value are then used to build a GGM gene co-expression network. Hence, the algorithm is named as <B>rgsGGM</B> (random gene sampling-based graphical Gaussian model). For more details on rgsGGM, please refer to [Wang et al, 2022](https://github.com/MaShisongLab/rgsGGM#References).

## Table of Contents
- [Install](https://github.com/MaShisongLab/rgsGGM#Install)
- [Usage](https://github.com/MaShisongLab/rgsGGM#Usage)
- [References](https://github.com/MaShisongLab/rgsGGM#References)

## Install
This algorithm requires [MATLAB](https://www.mathworks.com/products/matlab.html). Copy the file `rgsGGM.m` to your working folder or the MATLAB scripts folder and start using it.

## Usage

The algorithm takes a gene expression matrix, the number of iterations, and the names of the genes within the matrix as inputs. The expression matrix should samples in rows and genes in columns. The sample numbers should be large and the low-expression genes should be filtered out first. 

<B>ggm = rgsGGM(`expression_matrix`,`number_of_iterations`,`gene_name`)</B>

`expression_matrix` - the gene expression martrix, samples in rows, genes in columns<br/>
`number_of_iterations` - the number of iteration used for pcor calculation, usually 20000<br/>
`gene_name` - the names for the genes in the matrix


Below, we use a human gene expression matrix extracted from the ARCHS4 database ([Lachmann *et al*, 2018](https://github.com/MaShisongLab/rgsGGM#References)) as an example to demonstrate how to conduct GGM gene co-expression network analysis via rgsGGM. A human gene expression file `human_transcript_v7.h5` was downloaded from [the ARCHS4 database](https://maayanlab.cloud/archs4/download.html) and processed into CPM gene expression values. The expression values are then log-transformed via log<sub>2</sub>(CPM + 1). After quality control, 24,959 human genes in 59,097 bulk RNA-seq samples are selected to generate a human gene expression matrix, contained within a file `human_expression_extracted_from_archs4_v7.h5`, which is available via [Figshare](https://figshare.com/s/ec58e5b149c3060e1a6f). 

Download the file `human_expression_extracted_from_archs4_v7.h5` from [Figshare](https://figshare.com/s/ec58e5b149c3060e1a6f) and place it into your working folder. Conduct GGM gene co-expression network analysis using the following commands within MATLAB.
```matlab
% read in the gene expression matrix
expression_matrix = h5read('human_expression_extracted_from_archs4_v7.h5','/expression');

% read in the gene names
gene_name = h5read ('human_expression_extracted_from_archs4_v7.h5','/gene_name');

% conduct GGM gene co-expression network analysis. This step might take a long time.
ggm = rgsGGM( expression_matrix, 20000, gene_name );

% inspect the results
ggm
ggm.SigEdges(1:10,:)

```

Here is a glimpse into the results:
```shell
ggm

ggm = 

  rgsGGM with properties:

             gene_num: 24959
            gene_name: {24959x1 cell}
             pcor_all: [24959x24959 double]
    pcor_sampling_num: [24959x24959 int16]
          samples_num: 59097
          RoundNumber: 20000
             SigEdges: [1254559x5 table]
                notes: []



ggm.SigEdges(1:10,:)

ans =

  10x5 table

          GeneA                GeneB            Pcor      SamplingTime       r   
    _________________    _________________    ________    ____________    _______

    'ENSG00000211752'    'ENSG00000211801'     0.11589        124         0.89774
    'ENSG00000211752'    'ENSG00000211799'    0.084472        121         0.87783
    'ENSG00000211752'    'ENSG00000211800'    0.062617        132         0.87172
    'ENSG00000211752'    'ENSG00000211792'     0.12971        145         0.88666
    'ENSG00000211752'    'ENSG00000211794'     0.05501        138         0.88525
    'ENSG00000211752'    'ENSG00000211788'    0.048441        124         0.89315
    'ENSG00000211752'    'ENSG00000211710'     0.09642        132         0.88509
    'ENSG00000211752'    'ENSG00000211734'     0.11832        130         0.89419
    'ENSG00000211752'    'ENSG00000211803'    0.034128        119         0.84301
    'ENSG00000211752'    'ENSG00000211804'    0.022549        128         0.81303
```
<B>ggm.SigEdges</B> contains the gene pairs with |pcors| >= 0.02, listing their gene names, partial correlation coefficients (pcor), the number of iteration they are sampled (SamplingTime), and their Pearson correlation coefficients (r).
```matlab
% save all gene pairs with pcor >= a cutoff value to a file for gene co-expression construction
i = ggm.SigEdges{:,4} >= 0.045
writetable(ggm.SigEdges(i,1:3),'human.ggm.network.txt','Delimiter','tab','WriteVariableNames',FALSE)
```
The network can then be clustered via network clustering algorithm to obtain gene co-expression module and used for down-stream analysis.

## References

Lachmann A, Torre D, Keenan AB, Jagodnik KM, Lee HJ, Wang L, Silverstein MC, and Ma'ayan A. 2018. Massive mining of publicly available RNA-seq data from human and mouse. *Nature Communications* 9:1366.

Ma S, Gong Q, and Bohnert HJ. 2007. An Arabidopsis gene network based on the graphical Gaussian model. *Genome Research* 17:1614-1625.

Wang Y, Zhang Y, Yu N, Li B, Gong J, Mei Y, Bao J, Ma S. 2022. Decoding transcriptional regulation via a human gene expression predictor. *submitted* 

