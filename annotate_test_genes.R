library(ioutil)
library(argparser)

##### parse arguments ######
args <- arg_parser("program");
args <- add_argument(args, "-genes", help="target genes file", default="/work-zfs/abattle4/ashis/progdata/covid19_eqtl/covid19_eqtl_test_genes.txt")
args <- add_argument(args, "-annot", help="gene annotation file (txt)", default="/work-zfs/abattle4/lab_data/annotation/gencode.v26/gencode.v26.annotation.gene.txt")
args <- add_argument(args, "-map", help="gene mappability file", default="/work-zfs/abattle4/lab_data/annotation/mappability_hg38_gencode26/hg38_gene_mappability.txt")
args <- add_argument(args, "-th", help="minimum threshold for gene mappability", default=0.8)
args <- add_argument(args, "-oann", help="output file with gene annotations", default="results/genes_annotated.txt")
args <- add_argument(args, "-oensembl", help="output file with filtered ensembl gene ids", default="results/genes_ensemblid.txt")

argv = parse_args(args)

genes_fn = argv$genes
gene_annot_fn = argv$annot
gene_mappability_fn = argv$map
min_gene_mappability = argv$th
annotated_genes_fn = argv$oann
ensembl_fn = argv$oensembl

# gene_annot_fn = "/work-zfs/abattle4/lab_data/annotation/gencode.v26/gencode.v26.annotation.gene.txt"
# genes_fn = "/work-zfs/abattle4/ashis/progdata/covid19_eqtl/covid19_eqtl_test_genes.txt"
# annotated_genes_fn = "/work-zfs/abattle4/ashis/progdata/covid19_eqtl/covid19_eqtl_test_genes_annotated.txt"
# ensembl_fn = "/work-zfs/abattle4/ashis/progdata/covid19_eqtl/covid19_eqtl_test_genes_ensembl_id.txt"
# gene_mappability_fn = "/work-zfs/abattle4/lab_data/annotation/mappability_hg38_gencode26/hg38_gene_mappability.txt"
# min_gene_mappability = 0.8

gene_annot_df = read_df(gene_annot_fn, header = T, row.names = F)
genes_df = read_df(genes_fn, header = F, row.names = F)
colnames(genes_df) = "gene_name"
gene_mappability_df = read_df(gene_mappability_fn, header = F, row.names = F)
colnames(gene_mappability_df) = c('gene_id', 'gene_mappability')

annotated_genes_df = merge(genes_df, gene_annot_df)
annotated_genes_df = merge(annotated_genes_df, gene_mappability_df, by = 'gene_id', all.x = T, all.y = F)
annotated_genes_df = annotated_genes_df[annotated_genes_df$gene_mappability >= min_gene_mappability, , drop = F]

write_df(annotated_genes_df, file = annotated_genes_fn, row.names = F, col.names = T)
write_df(annotated_genes_df$gene_id, file = ensembl_fn, row.names = F, col.names = F)
