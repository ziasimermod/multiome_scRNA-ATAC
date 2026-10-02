#!/usr/bin/env Rscript
# Run after Step 06: Rscript analysis/07d_whole_wnn_cluster_rna_atac_baseline.R
# Main teaching pipeline baseline: whole approved WNN clusters, separately in each
# cohort. One pooled CON and one pooled HFD library per cohort => descriptive
# effect sizes and cell/feature counts only, with no mouse-level diet P value.
options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({
 library(SeuratObject); library(Signac); library(Matrix)
 library(dplyr); library(tidyr); library(tibble); library(readr); library(ggplot2)
})
root <- normalizePath('.')
stopifnot(file.exists(file.path(root,'multiome_scRNA-ATAC.Rproj')))
Sys.setenv(MULTIOME_DATASET_VERSION='v2_resequenced',MULTIOME_COHORT_ID='wt')
source(file.path(root,'config','project_config.R'))
helper <- file.path(root, "analysis", "R",
 'wt_crypt_followup_helpers.R')
source(helper)
base <- file.path(RESULTS_ROOT,DATASET_VERSION,'independent')
out <- file.path(base,'whole_wnn_cluster_baseline_2026-10-01')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
MIN_N <- 30L
# Every cohort is analyzed within its own approved Step 06 WNN clusters.
cohorts <- c('wt','ppar','il17')
all_counts <- list(); all_summary <- list(); all_top <- list(); manifest <- list()
for (cohort in cohorts) {
 message('Whole WNN baseline: ',cohort)
 objpath <- file.path(base,cohort,'06_annotation','multiome_annotated.rds')
 decisionpath <- file.path(root,'config','datasets',DATASET_VERSION,'cohorts',cohort,'cell_state_annotations.csv')
 need(file.exists(objpath) && file.exists(decisionpath),paste('Missing Step 06 inputs for',cohort))
 dest <- file.path(out,cohort); dir.create(dest,recursive=TRUE,showWarnings=FALSE)
 decision <- read_csv_checked(decisionpath)
 require_columns(decision,c('cluster','cell_state','approved','annotation_confidence'),'annotation decision')
 need(!anyDuplicated(decision$cluster) && all(decision$approved),paste('Unapproved or duplicated annotation in',cohort))
 obj <- readRDS(objpath)
 need(all(c('RNA','ATAC') %in% names(obj)),paste('Missing RNA/ATAC in',cohort))
 meta <- obj[[]]; meta$cell_barcode <- rownames(meta)
 require_columns(meta,c('wnn_cluster','cell_state','diet','cell_barcode'),'annotated object')
 meta$wnn_cluster <- as.character(meta$wnn_cluster)
 meta$diet <- as.character(meta$diet)
 need(!anyDuplicated(meta$cell_barcode) &&
   setequal(unique(meta$diet),c('CON','HFD')) &&
   setequal(unique(meta$wnn_cluster),as.character(decision$cluster)),
   paste('WNN clusters or diet assignment changed in',cohort))
 key <- match(meta$wnn_cluster,as.character(decision$cluster))
 need(!anyNA(key) && identical(as.character(meta$cell_state),as.character(decision$cell_state[key])),
   paste('Step 06 cell-state assignment differs from approved decisions in',cohort))
 units <- meta |>
   dplyr::count(wnn_cluster,cell_state,diet,name='n') |>
   tidyr::complete(tidyr::nesting(wnn_cluster,cell_state),diet=c('CON','HFD'),fill=list(n=0)) |>
   tidyr::pivot_wider(names_from=diet,values_from=n,names_prefix='n_') |>
   dplyr::mutate(cohort=cohort,eligible=n_CON>=MIN_N & n_HFD>=MIN_N,
     annotation_confidence=decision$annotation_confidence[match(wnn_cluster,as.character(decision$cluster))],
     .before=1)
 need(sum(units$n_CON+units$n_HFD)==nrow(meta),paste('Missing nuclei in',cohort))
 write_table(units,dest,'01_cluster_diet_counts_and_eligibility.csv')
 all_counts[[cohort]] <- units
 manifest[[cohort]] <- tibble(cohort=cohort,source=normalizePath(objpath),
   bytes=file.info(objpath)$size,modified=as.character(file.info(objpath)$mtime),
   n_nuclei=nrow(meta),n_clusters=nrow(units),n_eligible=sum(units$eligible))
 for(modality in c('RNA','ATAC')) {
   assay <- obj[[modality]]
   need('counts' %in% SeuratObject::Layers(assay),paste('Missing joined counts:',cohort,modality))
   counts <- align_columns(SeuratObject::LayerData(assay,layer='counts'),meta$cell_barcode,
     paste(cohort,modality))$matrix
   need(!anyDuplicated(rownames(counts)),paste('Duplicated features:',cohort,modality))
   depth <- Matrix::colSums(counts)
   need(all(is.finite(depth)),paste('Invalid depth:',cohort,modality))
   dir.create(file.path(dest,modality),showWarnings=FALSE)
   stats <- list(); top <- list()
   for(j in seq_len(nrow(units))) {
     u <- units[j,]
     if (!u$eligible) next
     ix <- which(meta$wnn_cluster==u$wnn_cluster)
     message('  ',modality,' WNN ',u$wnn_cluster,': ',u$n_CON,' CON / ',u$n_HFD,' HFD')
     eff <- pooled_effect(counts,meta,ix,modality,total_depth=depth) |>
       dplyr::mutate(cohort=cohort,wnn_cluster=u$wnn_cluster,cell_state=u$cell_state,.before=1)
     # One full, independently indexed table per whole WNN cluster. No
     # subclusters, score gates, or Gourab marker filters enter this baseline.
     write_table(eff,file.path(dest,modality),paste0('WNN_',u$wnn_cluster,'_pooled_effects.csv.gz'))
     stats[[length(stats)+1L]] <- tibble(cohort=cohort,modality=modality,
       wnn_cluster=u$wnn_cluster,cell_state=u$cell_state,n_CON=u$n_CON,n_HFD=u$n_HFD,
       total_features=nrow(eff),evaluable_features=sum(eff$evaluable),
       clear_open_or_up=sum(eff$clear_effect & eff$log2FC_HFD_vs_CON>0),
       clear_close_or_down=sum(eff$clear_effect & eff$log2FC_HFD_vs_CON<0),
       median_abs_log2FC_evaluable=if(any(eff$evaluable)) median(abs(eff$log2FC_HFD_vs_CON[eff$evaluable])) else NA_real_)
     top[[length(top)+1L]] <- eff |>
       dplyr::filter(evaluable,clear_effect) |>
       dplyr::group_by(direction) |>
       dplyr::slice_max(abs(log2FC_HFD_vs_CON),n=100,with_ties=FALSE) |>
       dplyr::ungroup()
     rm(eff); invisible(gc())
   }
   write_table(dplyr::bind_rows(stats),dest,paste0('02_',modality,'_whole_cluster_effect_summary.csv'))
   write_table(dplyr::bind_rows(top),dest,paste0('03_',modality,'_top_descriptive_effects.csv.gz'))
   all_summary[[paste(cohort,modality,sep='_')]] <- dplyr::bind_rows(stats)
   rm(counts,assay,depth,stats,top); invisible(gc())
 }
 rm(obj,meta,units,decision); invisible(gc())
}
write_table(dplyr::bind_rows(manifest),out,'00_input_manifest.csv')
write_table(dplyr::bind_rows(all_counts),out,'01_all_cohort_whole_cluster_counts.csv')
summary <- dplyr::bind_rows(all_summary)
write_table(summary,out,'02_all_cohort_whole_cluster_effect_summary.csv')
# Counts and feature-effect summaries are cohort-specific. Cluster IDs do not
# encode cross-genotype homology and filtered ATAC peak sets need exact interval
# overlap or fragment recount before any interval-level genotype comparison.
p <- summary |>
 tidyr::pivot_longer(c(clear_open_or_up,clear_close_or_down),names_to='direction',values_to='n') |>
 ggplot(aes(x=reorder(paste0('W',wnn_cluster,' ',cell_state),n),y=n,fill=direction)) +
 geom_col() + coord_flip() + facet_grid(modality~cohort,scales='free_y',space='free_y') +
 scale_fill_manual(values=c(clear_open_or_up='#E69F00',clear_close_or_down='#5A3E9A'),
  labels=c('HFD higher','HFD lower')) +
 labs(x='Whole WNN cluster and approved cell state',y='Descriptive clear effects',
  title='Whole annotated WNN cluster baseline: pooled RNA and ATAC effects',
  subtitle='Each cluster compared within cohort; eligible groups have >=30 nuclei per diet',
  caption='Thresholded effect counts depend on coverage, detection and cluster size. No mouse-level replicated diet test.') +
 theme_bw(base_size=10) + theme(legend.title=element_blank(),plot.caption=element_text(hjust=0))
ggplot2::ggsave(file.path(out,'04_whole_cluster_effect_counts.pdf'),p,width=17,height=13,
 device=grDevices::pdf,useDingbats=FALSE)
writeLines(c('Whole annotated WNN cluster RNA/ATAC baseline complete.',
 'One pooled library per diet; no biological-replicate P values.',
 'Whole-cluster IDs are cohort-specific; do not directly join genotype peaks without exact coordinate parity.'),
 file.path(out,'07D_COMPLETE.txt'))
message('Finished: ',out)
