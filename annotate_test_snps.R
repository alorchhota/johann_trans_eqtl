library(ioutil)
library(genomicsutil)
library(argparser)

##### parse arguments ######
args <- arg_parser("program");
args <- add_argument(args, "-snps", help="target snps file", default="data/lbm_snps_hg38.txt")
args <- add_argument(args, "-annot", help="snp annotation file (txt)", default="data/variant_annot_MAF_05_not_in_repeat.txt.gz")
args <- add_argument(args, "-th", help="minimum threshold for snp mappability", default=1)
args <- add_argument(args, "-oann", help="output file with snp annotations", default="results/lbm_snps_annotated.txt")
args <- add_argument(args, "-osnp", help="output file with filtered snp ids", default="results/lbm_snps_ids.txt")

argv = parse_args(args)

snps_fn = argv$snps
snp_annot_fn = argv$annot
min_snp_mappability = argv$th
annotated_snps_fn = argv$oann
filtered_snps_fn = argv$osnp

### read target snps
snps_df = read_df(snps_fn, header = T, row.names = F)
snps_df$chr = extend_chr(snps_df$chr)
colnames(snps_df) = gsub(pattern = "snp", replacement = 'rsid', x = colnames(snps_df))

### read snp annotations
snp_annot_df = read_df(snp_annot_fn, sep = '\t', header = T, row.names = F)

### get target snp annotations
target_snps_annot = merge(snps_df, snp_annot_df, by = c('chr', 'pos'))
filetered_snps = target_snps_annot$snp[target_snps_annot$mappability >= min_snp_mappability]

### save
write_df(target_snps_annot, file = annotated_snps_fn, row.names = F, col.names = T)
write_df(filetered_snps, file = filtered_snps_fn, row.names = F, col.names = F)
