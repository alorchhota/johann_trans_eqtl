##############################################################################
### compute genome-wide trans-eqtl (inter-chromosomal) stats of given genes
### filtered for cross-mapping 
##############################################################################
library(ioutil)
library(miscutil)
library(MatrixEQTL)
library(argparser)
library(stringr)

##### parse arguments ######
args <- arg_parser("program");
args <- add_argument(args, "-expr", help="expression file", default="/work-zfs/abattle4/lab_data/GTEx_v8_trans_eqtl_data_processed_by_brian/processed_expression/Small_Intestine_Terminal_Ileum.v8.normalized_expression.bed")
args <- add_argument(args, "-cov", help="covariate file", default="/work-zfs/abattle4/lab_data/GTEx_v8_trans_eqtl_data_processed_by_brian/covariates/Small_Intestine_Terminal_Ileum.v8.covariates.txt")
args <- add_argument(args, "-geno_pfx", help="genotype file name prefix with path (before chromosome in the name)", default="/work-zfs/abattle4/lab_data/GTEx_v8_trans_eqtl_data_processed_by_brian/processed_genotypes/GTEx_Analysis_2017-06-05_v9_WholeGenomeSeq_838Indiv_")
args <- add_argument(args, "-geno_chr", help="comma separated chromosomes", default="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX")
args <- add_argument(args, "-geno_sfx", help="genotype file name suffix (after chromosome in the name)", default="_dosage_MAF_05_not_in_repeat.RData")
args <- add_argument(args, "-geno_na", help="character for NA in the genotype file", default="-")
args <- add_argument(args, "-genes", help="target genes file", default="data/lbm_genes_ensemblid.txt")
args <- add_argument(args, "-snps", help="target snps file", default="data/lbm_snps_hg38_snpid.txt")
args <- add_argument(args, "-annot", help="gene annotation file (txt)", default="/work-zfs/abattle4/lab_data/annotation/gencode.v26/gencode.v26.annotation.gene.txt")
args <- add_argument(args, "-crossmap", help="crossmapping file", default='/work-zfs/abattle4/lab_data/annotation/mappability_hg38_gencode26/hg38_cross_mappability_strength.txt')
args <- add_argument(args, "-d", help="distance threshold for cross-mappability", default=1e6)
args <- add_argument(args, "-p", help="pvalue threshold to save test statistic", default=1e-5)
args <- add_argument(args, "-o", help="annotated eqtl file", default="results/target_eqtls")

argv = parse_args(args)
expr_fn = argv$expr
cov_fn = argv$cov
geno_pfx = argv$geno_pfx
geno_chromosomes_input = argv$geno_chr
geno_sfx = argv$geno_sfx
na_genotype_value = argv$geno_na
target_genes_fn = argv$genes
target_snps_fn = argv$snps
gene_annot_fn = argv$annot
crossmap_fn = argv$crossmap
snp_gene_dist_for_crossmap = argv$d
p_threshold = argv$p
out_pfx = argv$o

### process settings and annotation
geno_chromosomes = parse_delimitted_string(geno_chromosomes_input, delim = ',', rm.empty = T)

stopifnot(file.exists(target_genes_fn))
target_genes_df = read_df(target_genes_fn, header = F, row.names = F)
target_genes = unique(target_genes_df[,1])

target_snps_fn = str_trim(target_snps_fn)
target_snps = NULL
if(nchar(target_snps_fn) > 0){
  stopifnot(file.exists(target_snps_fn))
  target_snps_df = read_df(target_snps_fn, header = F, row.names = F)
  target_snps = unique(target_snps_df[,1])
}

gencode_df = read_df(gene_annot_fn, header = T, row.names = F)
target_genes_chromosomes = as.character(sapply(target_genes, function(g) gencode_df[gencode_df$gene_id ==g ,'chr'] ))
target_genes_per_chr = tapply(X = target_genes, INDEX = target_genes_chromosomes, FUN = c)

### read expression
expr_df = read_df(expr_fn, header = T, row.names = F)
rownames(expr_df) = expr_df[,'gene_id']
tested_target_genes = intersect(target_genes, rownames(expr_df))
expr_df = expr_df[tested_target_genes ,5:ncol(expr_df), drop = F]
stopifnot(nrow(expr_df) >=1 )

### read covariate
cov_df = read_df(cov_fn, header = T, row.names = T)
stopifnot(length(intersect(colnames(expr_df), colnames(cov_df))) == ncol(expr_df))
cov_df = cov_df[,colnames(expr_df), drop = F]
n_uniq = sapply(rownames(cov_df), function(idx)length(unique(as.numeric(cov_df[idx,]))))
if(any(n_uniq<=1)){
  const_covariates = names(n_uniq[n_uniq<=1])
  cov_df = cov_df[-which(rownames(cov_df) %in% const_covariates), , drop=F]
}

### function to read genotype file
read_genotype_file <- function(geno_fn){
  stopifnot(endsWith(geno_fn, ".RData"))
  load(geno_fn)
  colnames(genotype_mat_not_in_repeat) = gsub('\\.', '-', colnames(genotype_mat_not_in_repeat))
  genotype_mat_not_in_repeat[genotype_mat_not_in_repeat==na_genotype_value] = NA
  for(cn in colnames(genotype_mat_not_in_repeat))
    genotype_mat_not_in_repeat[,cn] = as.numeric(genotype_mat_not_in_repeat[,cn])
  return(genotype_mat_not_in_repeat)
}


############ filter eqtls for cross-mappability ########################
verbose_print("processing data to filter cross-mappable eqtls ...")
crossmap_df = read_df(crossmap_fn, sep = '\t', header = F, row.names = F)

### crossmap to or from target genes
filter_crossmap_by_target_genes <- function(crossmap_df, incl.genes){
  # crossmap_df : cross-mappability data frame where first two columns are genes
  stopifnot(class(crossmap_df)=='data.frame' && ncol(crossmap_df)>=2)
  stopifnot(class(crossmap_df[,1])=='character' && class(crossmap_df[,2])=='character')
  stopifnot(class(incl.genes) == 'character')
  
  # gene_fac = as.factor(crossmap_df[,1])
  # gene_levels = levels(gene_fac)
  # gene_levels_included = gene_levels %in% incl.genes
  # gene1_included = gene_levels_included[as.integer(gene_fac)]
  
  gene_fac = as.factor(crossmap_df[,2])
  gene_levels = levels(gene_fac)
  gene_levels_included = gene_levels %in% incl.genes
  gene2_included = gene_levels_included[as.integer(gene_fac)]
  
  #return(crossmap_df[gene1_included | gene2_included, , drop=F])
  return(crossmap_df[gene2_included, , drop=F])
}
crossmap_df = filter_crossmap_by_target_genes(crossmap_df, incl.genes = tested_target_genes)


### function: fast acess genes location (crh, tss)
gene_loc_env = new.env(hash = T)
tss_values =  as.integer(apply(gencode_df, MARGIN = 1, FUN = function(row){
  ifelse(row['strand']=="+", row['start_pos'], row['end_pos'])
}))
gencode_df$tss = tss_values
tmp <- mapply(function(g, chr, tss){
  gene_loc_env[[g]] <<- list(chr=chr, tss=tss); 
  return()
}, gencode_df$gene_id, gencode_df$chr, gencode_df$tss)
get_gene_location <- function(g){
  return(list(chr=gene_loc_env[[g]]$chr, tss=gene_loc_env[[g]]$tss))
}

### function: given chr, get annotation of all genes in the given chr
chromosomes = unique(gencode_df$chr)
gene_anot_by_chr = lapply(chromosomes, function(chr) {gencode_df[gencode_df$chr==chr,,drop=F]})
names(gene_anot_by_chr) = chromosomes
get_gene_annot_in_chr <- function(chr){
  return(gene_anot_by_chr[[chr]])
}

### function: given a position, get nearby genes
### note: the genes do not have to be in the list of test genes
get_nearby_genes <- function(chr, pos, d=1e6){
  chr_annot_df = get_gene_annot_in_chr(chr)
  left_pos = pos - d
  right_pos = pos + d
  near_annot_df = chr_annot_df[chr_annot_df$tss >= left_pos & chr_annot_df$tss <=right_pos, ]
  return(near_annot_df$gene_id)
}

### for each gene, get cross-mappable genes from the given gene in different chr
### note: the cross-mappable genes must be in the list of test genes
crossmap_per_gene = tapply(crossmap_df[,2], crossmap_df[,1], c)
trans_crossmap_per_gene_env = new.env(hash = T)
tmp = lapply(names(crossmap_per_gene), function(g){
  # cg = intersect(crossmap_per_gene[[g]], tested_genes) # cross-map gene must be tested
  cg = crossmap_per_gene[[g]]
  g_chr = get_gene_location(g)$chr
  cg_chr = sapply(cg, function(x) get_gene_location(x)$chr)
  trans_cg = cg[cg_chr != g_chr]
  trans_crossmap_per_gene_env[[g]] <<- trans_cg
  return()
})

get_trans_crossmap_genes <- function(g){
  return(trans_crossmap_per_gene_env[[g]])
}

### function: given a snp location, 
### find which cross-mappable trans genes were tests with the snp
get_crossmappable_trans_genes_per_snp <- function(chr, pos, d=1e6){
  genes = get_nearby_genes(chr, pos, d)  
  cross_genes_list <- lapply(genes, get_trans_crossmap_genes) # TODO: p1
  cross_genes = unique(unlist(cross_genes_list))
  return(cross_genes)
}

### function: given a snp id, find it's location
get_snp_loc_from_id <- function(snpid){
  parts = strsplit(snpid, split = '_')[[1]]
  chr = ifelse(startsWith(parts[1], prefix = 'chr'), parts[1], paste0('chr', parts[1]))
  pos = as.integer(parts[2])
  return(list(chr = chr, pos = pos))
}

#######################################
### trans-association per chromosomes
trans_assoc_per_chr = lapply(geno_chromosomes, function(chr){
  verbose_print(sprintf("calling eqtls for chromosome %s ...", chr))
  chr_expr_df = expr_df
  chr_trans_genes = rownames(chr_expr_df)
  if(chr %in% names(target_genes_per_chr)){
    chr_trans_genes = setdiff(rownames(chr_expr_df), target_genes_per_chr[[chr]])
    chr_expr_df = chr_expr_df[chr_trans_genes, , drop = F]
    if(length(chr_trans_genes) <= 0 )
      return(NULL)
  }
  
  ### read genotypes
  geno_fn = sprintf("%s%s%s", geno_pfx, chr, geno_sfx)
  geno_df = read_genotype_file(geno_fn)
  
  if(length(target_snps) > 0){
    chr_snps = intersect(rownames(geno_df), target_snps)
    geno_df = geno_df[chr_snps, , drop = F]
  }
  
  ### prepare data for matrix-eqtl
  common_samples = Reduce(intersect, list(colnames(chr_expr_df), colnames(geno_df), colnames(cov_df)))
  meqtl_expr_slice = SlicedData$new(as.matrix(chr_expr_df[,common_samples,drop=F]))
  meqtl_snp_slice = SlicedData$new(as.matrix(geno_df[,common_samples,drop=F]))
  meqtl_cov_slice = SlicedData$new(as.matrix(cov_df[,common_samples,drop=F]))
  
  ### run matrix-eQTL
  me = Matrix_eQTL_engine(snps = meqtl_snp_slice,
                          gene = meqtl_expr_slice,
                          cvrt = meqtl_cov_slice,
                          output_file_name = NULL,
                          pvOutputThreshold = p_threshold,
                          useModel = modelLINEAR, 
                          verbose = FALSE,
                          pvalue.hist = FALSE,
                          min.pv.by.genesnp = FALSE,
                          noFDRsaveMemory = FALSE)
  
  me$all$eqtls$snps = as.character(me$all$eqtls$snps)  # convert factor to chracter
  me$all$eqtls$gene = as.character(me$all$eqtls$gene)  # convert factor to chracter
  
  ### compute cross-mappability from snp in current chr to target genes
  chr_snps = rownames(geno_df)
  chr_snps_parts = strsplit(chr_snps, split = "_")
  chr_snps_positions = as.integer(sapply(chr_snps_parts, function(x) x[2]))
  chr_crossmaps = lapply(seq_along(chr_snps), function(si){
    cg = get_crossmappable_trans_genes_per_snp(chr = chr, pos = chr_snps_positions[si], d = snp_gene_dist_for_crossmap)
    cg = intersect(cg, chr_trans_genes)  # we are interested only in target genes
    if(length(cg) > 0){
      cm = data.frame(snps = chr_snps[si], gene = cg, crossmap = T, stringsAsFactors = F)
      return(cm)
    } else {
      return(NULL)
    }
  })
  chr_crossmaps_df = do.call(rbind, chr_crossmaps)
  if(is.null(chr_crossmaps_df)){
    chr_crossmaps_df = data.frame(snps = 'dummy', gene = 'dummy', crossmap = NA, stringsAsFactors = F)
    chr_crossmaps_df = chr_crossmaps_df[-1,,drop=F]
  }
  
  ### compute total number of test and total number of cross-mappable tests
  total_snp_gene_pairs = nrow(geno_df) * nrow(chr_expr_df)
  total_crossmappable_snp_gene_pairs = nrow(chr_crossmaps_df)
  total_not_crossmappable_snp_gene_pairs = total_snp_gene_pairs - total_crossmappable_snp_gene_pairs
  
  ### annotate cross-mappability of tests with p<1e-5
  chr_crossmaps_df = chr_crossmaps_df[chr_crossmaps_df$snps %in% unique(me$all$eqtls$snps), , drop = F]
  eqtl_stat_df = merge(me$all$eqtls, chr_crossmaps_df, by = c('snps', 'gene'), all.x = T, all.y = F)

  ### clean memory  
  rm(me, chr_expr_df, geno_df, 
     meqtl_expr_slice, meqtl_snp_slice, meqtl_cov_slice, 
     chr_snps_parts, chr_snps_positions, chr_crossmaps, chr_crossmaps_df)
  gc(reset = T)
  
  return(list(stats = eqtl_stat_df, 
              snps = chr_snps,
              genes = chr_trans_genes,
              total_snps = length(chr_snps),
              total_genes = length(chr_trans_genes),
              total_snp_gene_pairs = total_snp_gene_pairs, 
              total_crossmappable_snp_gene_pairs = total_crossmappable_snp_gene_pairs,
              total_not_crossmappable_snp_gene_pairs = total_not_crossmappable_snp_gene_pairs))
})
names(trans_assoc_per_chr) = geno_chromosomes

### save trans-association stats
trans_res_fn = sprintf("%s_trans_assoc_p_%s.rds", out_pfx, gsub(as.character(p_threshold), pattern = "-", replacement = "_"))
saveRDS(trans_assoc_per_chr, file = trans_res_fn)


### merge all tests
verbose_print("merging all tests ...")
all_eqtls_df = NULL
for(tapc in trans_assoc_per_chr){
  all_eqtls_df = rbind(all_eqtls_df, tapc$stats)
}

### filter cross-mappabile eqtls and compute FDR
not_cross_mappable_eqtls_df = all_eqtls_df[is.na(all_eqtls_df$crossmap) | all_eqtls_df$crossmap == F, ,drop=F]
total_not_crossmappable_snp_gene_pairs = sum(sapply(trans_assoc_per_chr, function(x) x$total_not_crossmappable_snp_gene_pairs))
not_cross_mappable_eqtls_df$FDR = p.adjust(not_cross_mappable_eqtls_df$pvalue, method = 'BH', n = total_not_crossmappable_snp_gene_pairs)
not_cross_mappable_eqtls_df = not_cross_mappable_eqtls_df[order(not_cross_mappable_eqtls_df$pvalue),,drop=F]

### compute gene-level FDR
min_pvalue_per_gene = tapply(not_cross_mappable_eqtls_df$pvalue, not_cross_mappable_eqtls_df$gene, FUN = min)
adjusted_min_pvalue_per_gene = min_pvalue_per_gene * 1e6
adjusted_min_pvalue_per_gene[adjusted_min_pvalue_per_gene >= 1] = 1
genes_tested = Reduce(union, lapply(trans_assoc_per_chr, function(x) x$genes))
n_genes_tested = length(genes_tested)
gene_level_FDR = p.adjust(adjusted_min_pvalue_per_gene, method = "BH", n = n_genes_tested)
not_cross_mappable_eqtls_df$egene_min_pvalue = min_pvalue_per_gene[not_cross_mappable_eqtls_df$gene]
not_cross_mappable_eqtls_df$egene_FDR = gene_level_FDR[not_cross_mappable_eqtls_df$gene]

### save crossmap-filtered eqtls
eqtl_fn = sprintf("%s_crossmap_filtered_trans_eqtls_p_%s.txt.gz", out_pfx, gsub(as.character(p_threshold), pattern = "-", replacement = "_"))
write_df(not_cross_mappable_eqtls_df, file = gzfile(eqtl_fn), row.names = F, col.names = T)
