### Load iRefIndex annotations

``` r
library('iRefR')
```

    ## Loading required package: igraph
    ## Loading required package: graph
    ## 
    ## Attaching package: 'graph'
    ## 
    ## The following objects are masked from 'package:igraph':
    ## 
    ##     degree, edges
    ## 
    ## Loading required package: RBGL
    ## 
    ## Attaching package: 'RBGL'
    ## 
    ## The following object is masked from 'package:igraph':
    ## 
    ##     transitivity

``` r
library('dplyr')
```

    ## 
    ## Attaching package: 'dplyr'
    ## 
    ## The following object is masked from 'package:graph':
    ## 
    ##     union
    ## 
    ## The following object is masked from 'package:stats':
    ## 
    ##     filter
    ## 
    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library('readr')

# Settings
target_species = 'Mus musculus'

if (target_species == 'Homo sapiens') {
    taxon_id    = 9606
} else if (target_species == 'Mus musculus') {
    taxon_id    = 10090
}

# Organism-specific info
orgdb_name = sub(' ', '.', target_species)
taxon_field = sprintf('taxid:%s(%s)', taxon_id, target_species)

# Load human interaction data
input_dir = file.path(Sys.getenv('REF'), 'irefindex', '14.0')

# 2015/05/26 There appears to be a bug which resulted in the compressed
# .txt.zip having a different name from the .txt itself (0407 vs. 0704)
compressed_file = file.path(input_dir, sprintf('%s.mitab.07042015.txt.zip',
                                                taxon_id))
network_filename = sprintf('%s.mitab.04072015.txt', taxon_id)

# Load iRefIndex
iref = unique(read_delim(
    unz(description=compressed_file, filename=network_filename), 
    delim='\t', quote=''
))
```

    ## 
    |===========                                                                    |  13%   10 MB
    |=================                                                              |  21%   15 MB
    |=======================                                                        |  28%   21 MB
    |============================                                                   |  36%   26 MB
    |==================================                                             |  43%   31 MB
    |========================================                                       |  50%   37 MB
    |==============================================                                 |  57%   42 MB
    |===================================================                            |  64%   47 MB
    |=========================================================                      |  71%   52 MB
    |==============================================================                 |  78%   57 MB
    |====================================================================           |  85%   62 MB
    |==========================================================================     |  92%   67 MB
    |===============================================================================|  99%   72 MB
    |================================================================================| 100%   73 MB

``` r
# fix first colname
colnames(iref) = c('uidA', colnames(iref)[-c(1)])

# remove cross-species ixns
iref = iref %>% filter(taxa == taxon_field & taxb == taxon_field)
```

### Convert protein ids to entrez gene ids

**NOTE**

Occasionally icrogids will map to more than one gene identifier. For example, icrogid 4397381 maps to Entrez gene ids "100271849" and "4207", which correspond to "MEF2B myocyte enhancer factor 2B" and "EF2BNB-MEF2B MEF2BNB-MEF2B readthrough".

To keep things simple for now, whenever there is more than one possible mapped id, we will keep the first id matched.

``` r
# load mapping
#mapping_file = file.path(input_dir, sprintf('%s.RData', 
#                                            '9606.mitab.04072015.conversion_table'))

## If mapping file already exists, load it it in
#if (file.exists(mapping_file)) {
#    load(mapping_file)
#    id_mapping = extended_table
#    rm(extended_table)
#} else {
    # Otherwise, generate a mapping
id_mapping = create_id_conversion_table(
    iref, input_dir, sprintf('%s.mitab.04072015.conversion_table', taxon_id))
```

    ## Generating ID Conversion table...
    ## 25% completed...
    ## 50% completed...
    ## 75% completed...
    ## Conversion table has been copied to:
    ## /cbcb/lab/nelsayed/ref_data/irefindex/14.0/10090.mitab.04072015.conversion_table.txt

``` r
#}
```

``` r
# Source protein
gene_ids_a = c()

for (x in iref$icrogida) {
    gene_ids = convert_protein_ID('icrogid', x, 'entrezgene/locuslink', id_mapping)
    gene_ids = gene_ids[1]
    gene_ids_a = append(gene_ids_a, gene_ids)
}

# Target protein
gene_ids_b = c()

for (x in iref$icrogidb) {
    gene_ids = convert_protein_ID('icrogid', x, 'entrezgene/locuslink', id_mapping)
    gene_ids = gene_ids[1]

    gene_ids_b = append(gene_ids_b, gene_ids)
}

# Drop edges for which the source and target proteins could not both be mapped
# to Entrez genes
mask = !(is.na(gene_ids_a) | is.na(gene_ids_b))

gene_ids_a = gene_ids_a[mask]
gene_ids_b = gene_ids_b[mask]
```

### Convert to ENSEMBL IDs

``` r
# Load OrganismDb
library(orgdb_name, character.only=TRUE)
```

    ## Loading required package: AnnotationDbi
    ## Loading required package: stats4
    ## Loading required package: BiocGenerics
    ## Loading required package: parallel
    ## 
    ## Attaching package: 'BiocGenerics'
    ## 
    ## The following objects are masked from 'package:parallel':
    ## 
    ##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    ##     clusterExport, clusterMap, parApply, parCapply, parLapply,
    ##     parLapplyLB, parRapply, parSapply, parSapplyLB
    ## 
    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     combine, intersect, setdiff, union
    ## 
    ## The following object is masked from 'package:stats':
    ## 
    ##     xtabs
    ## 
    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, as.vector, cbind,
    ##     colnames, do.call, duplicated, eval, evalq, Filter, Find, get,
    ##     intersect, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rep.int, rownames, sapply, setdiff, sort,
    ##     table, tapply, union, unique, unlist, unsplit
    ## 
    ## Loading required package: Biobase
    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.
    ## 
    ## Loading required package: GenomeInfoDb
    ## Loading required package: S4Vectors
    ## 
    ## Attaching package: 'S4Vectors'
    ## 
    ## The following object is masked from 'package:dplyr':
    ## 
    ##     rename
    ## 
    ## The following object is masked from 'package:igraph':
    ## 
    ##     compare
    ## 
    ## Loading required package: IRanges
    ## 
    ## Attaching package: 'IRanges'
    ## 
    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     collapse, desc, slice
    ## 
    ## The following object is masked from 'package:igraph':
    ## 
    ##     simplify
    ## 
    ## 
    ## Attaching package: 'AnnotationDbi'
    ## 
    ## The following object is masked from 'package:GenomeInfoDb':
    ## 
    ##     species
    ## 
    ## The following object is masked from 'package:dplyr':
    ## 
    ##     select
    ## 
    ## Loading required package: OrganismDbi
    ## Loading required package: GenomicFeatures
    ## Loading required package: GenomicRanges
    ## Loading required package: GO.db
    ## Loading required package: DBI
    ## 
    ## Loading required package: org.Mm.eg.db
    ## 
    ## Loading required package: TxDb.Mmusculus.UCSC.mm10.knownGene

``` r
orgdb = get(orgdb_name)

# Create a mapping for all of the observed ENTREZ ids
entrez_to_ensembl_mapping = AnnotationDbi::select(
    orgdb, keys=gene_ids_a, keytype='ENTREZID', columns=c('ENSEMBL', 'GENENAME')
)
```

    ## Warning in .generateExtraRows(tab, keys, jointype): 'select' and duplicate query keys resulted in 1:many mapping
    ## between keys and return rows

    ## Warning in .generateExtraRows(tab, keys, jointype): 'select' and duplicate query keys resulted in 1:many mapping
    ## between keys and return rows

``` r
# Convert source and target gene ids
ensembl_ids_a = entrez_to_ensembl_mapping$ENSEMBL[
    match(gene_ids_a, entrez_to_ensembl_mapping$ENTREZID)
]
ensembl_ids_b = entrez_to_ensembl_mapping$ENSEMBL[
    match(gene_ids_b, entrez_to_ensembl_mapping$ENTREZID)
]
```

### Save resulting network

``` r
# create output dataframe
result = data.frame(a=ensembl_ids_a, b=ensembl_ids_b)

# remove any rows containing one or more NAs
result = result[complete.cases(result),]

# remove duplicate entries
result = result[!duplicated(result),]

# sort result
result = result[with(result, order(a,b)),]

# write output
outfile = sub('.mitab', '', sub('txt.zip', 'ensembl.csv', compressed_file))
write.csv(result, outfile, quote=FALSE, row.names=FALSE)

# gzip output
system(sprintf('(gzip -f %s) &', outfile))
```

System information
------------------

``` r
toLatex(sessionInfo())
```

\begin{itemize}\raggedright
  \item R version 3.1.3 (2015-03-09), \verb|x86_64-unknown-linux-gnu|
  \item Locale: \verb|LC_CTYPE=en_US.UTF-8|, \verb|LC_NUMERIC=C|, \verb|LC_TIME=en_US.UTF-8|, \verb|LC_COLLATE=en_US.UTF-8|, \verb|LC_MONETARY=en_US.UTF-8|, \verb|LC_MESSAGES=en_US.UTF-8|, \verb|LC_PAPER=en_US.UTF-8|, \verb|LC_NAME=C|, \verb|LC_ADDRESS=C|, \verb|LC_TELEPHONE=C|, \verb|LC_MEASUREMENT=en_US.UTF-8|, \verb|LC_IDENTIFICATION=C|
  \item Base packages: base, datasets, graphics, grDevices,
    methods, parallel, stats, stats4, utils
  \item Other packages: AnnotationDbi~1.28.2, Biobase~2.26.0,
    BiocGenerics~0.12.1, colorout~1.0-3, DBI~0.3.1, dplyr~0.4.1,
    GenomeInfoDb~1.2.5, GenomicFeatures~1.18.7,
    GenomicRanges~1.18.4, GO.db~3.0.0, graph~1.44.1, igraph~0.7.1,
    IRanges~2.0.1, iRefR~1.13, knitrBootstrap~1.0.0,
    Mus.musculus~1.1.2, OrganismDbi~1.8.1, org.Mm.eg.db~3.0.0,
    RBGL~1.42.0, readr~0.1.1, rmarkdown~0.6.1, RSQLite~1.0.0,
    S4Vectors~0.4.0, setwidth~1.0-3,
    TxDb.Mmusculus.UCSC.mm10.knownGene~3.0.0, vimcom~1.2-5
  \item Loaded via a namespace (and not attached): assertthat~0.1,
    base64enc~0.1-2, BatchJobs~1.6, BBmisc~1.9,
    BiocParallel~1.0.3, biomaRt~2.22.0, Biostrings~2.34.1,
    bitops~1.0-6, brew~1.0-6, checkmate~1.5.3, codetools~0.2-11,
    digest~0.6.8, evaluate~0.7, fail~1.2, foreach~1.4.2,
    formatR~1.2, GenomicAlignments~1.2.2, htmltools~0.2.6,
    iterators~1.0.7, knitr~1.10.5, lazyeval~0.1.10, magrittr~1.5,
    markdown~0.7.7, Rcpp~0.11.6, RCurl~1.95-4.6, Rsamtools~1.18.3,
    rtracklayer~1.26.3, sendmailR~1.2-1, stringi~0.4-1,
    stringr~1.0.0, tools~3.1.3, XML~3.98-1.2, XVector~0.6.0,
    yaml~2.1.13, zlibbioc~1.12.0
\end{itemize}
