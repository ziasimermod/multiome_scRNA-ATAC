#!/usr/bin/env Rscript
# Run: Rscript /scratch/dsaiz/multiome_scRNA-ATAC/analysis/08m_make_isc_meeting_figures.R --no-coverage
# Run without --no-coverage to add the three-cohort coverage panels.
# Frozen Step 08H/I/L populations and intervals. One pooled library per diet:
# no animal-level P values or biological confidence intervals are available.
options(stringsAsFactors = FALSE)
RUN_COVERAGE <- !('--no-coverage' %in% commandArgs(trailingOnly=TRUE))
suppressPackageStartupMessages({
  library(SeuratObject); library(Signac); library(ggplot2)
  library(dplyr); library(tidyr); library(readr); library(tibble)
  library(GenomicRanges); library(IRanges); library(patchwork)
})
REPO_ROOT <- "/scratch/dsaiz/multiome_scRNA-ATAC"
stopifnot(dir.exists(REPO_ROOT),
  file.exists(file.path(REPO_ROOT, 'multiome_scRNA-ATAC.Rproj')))
root <- normalizePath(REPO_ROOT, mustWork=TRUE)
Sys.setenv(MULTIOME_DATASET_VERSION='v2_resequenced', MULTIOME_COHORT_ID='wt')
source(file.path(root, 'config/project_config.R'))
wt <- OUTPUT_DIR
base <- file.path(RESULTS_ROOT, DATASET_VERSION, 'independent')
h <- file.path(wt, '08h_wt_ISC_region_freeze_sex_audit')
i <- file.path(wt, '08i_frozen_WT_peaks_genotype_context')
l <- file.path(wt, '08l_wt_ISC_open_close_motifs')
g <- file.path(wt, '08g_C0_C2_and_C0W0_ISC')
out <- file.path(wt, '08m_ISC_meeting_figures')
dir.create(out, recursive=TRUE, showWarnings=FALSE)
plots <- file.path(out, 'plots'); dir.create(plots, showWarnings=FALSE)
read <- function(path) { if (!file.exists(path)) stop('Missing: ', path); readr::read_csv(path, show_col_types=FALSE) }
save <- function(p, name, w=12, h=7) ggplot2::ggsave(file.path(plots,name), p, width=w, height=h, device=grDevices::pdf, useDingbats=FALSE)
require_cols <- function(x, cols, label) { z <- setdiff(cols,names(x)); if(length(z)) stop(label, ' missing: ',paste(z,collapse=', ')) }
cohorts <- c('wt','ppar','il17')
parents <- list(wt=c('0','8'),ppar=c('3','6'),il17=c('9','10'))
core <- c('Lgr5','Ascl2','Smoc2','Axin2','Lrig1','Rnf43')
cutoff <- read(file.path(g,'02_C0W0_gate_size_and_score_sensitivity.csv')) |>
  dplyr::filter(gate=='top15') |>
  dplyr::pull(minimum_primary_score)
stopifnot(length(cutoff)==1L,is.finite(cutoff))
regions <- read(file.path(l,'01_all_622_open_close_genotype_and_support.csv.gz'))
require_cols(regions,c('feature','chrom','start_1based','end_1based','WT_direction','C0_supported','count_supported_40','anchor_log2FC_HFD_vs_CON'), '08L regions')
stopifnot(nrow(regions)==622L,!anyDuplicated(regions$feature),sum(regions$C0_supported)==37L)

# Recorded 08H/08I memberships are authoritative: never re-rank by diet.
gate_files <- c(wt=file.path(h,'01_nucleus_level_sex_proxy_and_gate_audit.csv.gz'),
  ppar=file.path(i,'02_ppar_crypt_nuclei_scores.csv.gz'),
  il17=file.path(i,'02_il17_crypt_nuclei_scores.csv.gz'))
audit <- list(); selected <- list()
for (cohort in cohorts) {
  message('UMAP and gate: ',cohort)
  parent_ids <- parents[[cohort]]
  obj <- readRDS(file.path(base,cohort,'06_annotation','multiome_annotated.rds'))
  meta <- obj[[]]; meta$cell_barcode <- rownames(meta)
  require_cols(meta,c('cell_barcode','diet','wnn_cluster','cell_state'),'Step 06 metadata')
  stopifnot('wnn.umap' %in% names(obj),!anyDuplicated(meta$cell_barcode),setequal(unique(as.character(meta$diet)),c('CON','HFD')))
  emb <- SeuratObject::Embeddings(obj[['wnn.umap']])
  stopifnot(all(meta$cell_barcode %in% rownames(emb)))
  score <- read(gate_files[[cohort]])
  require_cols(score,c('cell_barcode','diet'),'gate table')
  stopifnot(!anyDuplicated(score$cell_barcode),all(score$cell_barcode %in% meta$cell_barcode))
  if(cohort=='wt') {
    require_cols(score,c('parent_wnn_cluster','C0W0_top15','ISC_core6_UCell','focused_cluster'),'WT gate')
    parent <- as.character(score$parent_wnn_cluster)
    is_sel <- as.logical(score$C0W0_top15)
    score$parent <- parent
  } else {
    require_cols(score,c('wnn_cluster','top15','ISC_core6_UCell'),'genotype gate')
    score$parent <- as.character(score$wnn_cluster)
    is_sel <- as.logical(score$top15) & score$parent %in% parent_ids
  }
  stopifnot(!anyNA(is_sel),identical(as.character(score$diet),as.character(meta$diet[match(score$cell_barcode,meta$cell_barcode)])),
    identical(score$parent, as.character(meta$wnn_cluster[match(score$cell_barcode,meta$cell_barcode)])))
  score$selected <- is_sel
  display_max <- max(.01, as.numeric(stats::quantile(score$ISC_core6_UCell[score$parent %in% parent_ids],.99,na.rm=TRUE)))
  score$cohort <- cohort
  selected[[cohort]] <- score |>
    dplyr::filter(selected) |>
    dplyr::select(cohort,cell_barcode,diet,parent,ISC_core6_UCell)
  audit[[cohort]] <- score |>
    dplyr::filter(.data$parent %in% .env$parent_ids) |>
    dplyr::group_by(parent,diet) |>
    dplyr::summarise(n_parent=dplyr::n(),n_selected=sum(selected),
      selected_median_score=if(any(selected)) median(ISC_core6_UCell[selected]) else NA_real_,.groups='drop') |>
    dplyr::mutate(cohort=.env$cohort,.before=1)
  xy <- tibble(cell_barcode=meta$cell_barcode,diet=as.character(meta$diet),
    cluster=as.character(meta$wnn_cluster),state=as.character(meta$cell_state),
    x=emb[meta$cell_barcode,1],y=emb[meta$cell_barcode,2]) |>
    dplyr::left_join(score |> dplyr::select(cell_barcode,parent,ISC_core6_UCell,selected),
      by='cell_barcode',relationship='one-to-one') |>
    dplyr::mutate(parent_cell=!is.na(parent) & .data$parent %in% .env$parent_ids,
      selected=replace_na(selected,FALSE),score_plot=if_else(parent_cell,ISC_core6_UCell,NA_real_))
  counts <- audit[[cohort]] |>
    dplyr::group_by(diet) |>
    dplyr::summarise(n=sum(n_selected),.groups='drop')
  label <- paste(paste(counts$diet,counts$n,sep=' n='),collapse='  |  ')
  rule <- if(cohort=='wt') 'WT: W0+C0 focused crypt gate; pooled top 15% across C0/W0 as recorded in 08H' else
    'Top 15% within each approved parent (both diets ranked together; score > 0)'
  note <- paste('Six-gene UCell (maxRank=400):',paste(core,collapse=', '),'\n',rule,'\nSelected:',label,
    '; WT reference minimum score =',sprintf('%.4f',cutoff),
    '; plasma upper limit = cohort parent 99th percentile',sprintf('%.3f',display_max))
  p <- ggplot(xy,aes(x,y)) +
    geom_point(color='#D9DAE0',size=.12,alpha=.35) +
    geom_point(data=xy |> dplyr::filter(parent_cell),aes(color=score_plot),size=.38,alpha=.8) +
    geom_point(data=xy |> dplyr::filter(selected),shape=21,stroke=.16,color='#16151C',aes(fill=score_plot),size=.85) +
    scale_color_viridis_c(option='plasma',limits=c(0,display_max),oob=scales::squish,na.value='grey80',name='UCell') +
    scale_fill_viridis_c(option='plasma',limits=c(0,display_max),oob=scales::squish,guide='none') +
    labs(title=paste(toupper(cohort),'| selected ISC-enriched nuclei'),subtitle='Pale: all other WNN nuclei; plasma: candidate parent nuclei; outlined: selected',
      caption=note,x='WNN UMAP 1',y='WNN UMAP 2') + theme_classic(base_size=12) +
    theme(plot.caption=element_text(size=9,hjust=0)) + coord_equal()
  save(p,paste0('01_',cohort,'_gate_WNN.pdf'))
  p_split <- p + facet_wrap(~diet,nrow=1) + labs(title=paste(toupper(cohort),'| selected nuclei by diet'))
  save(p_split,paste0('02_',cohort,'_gate_WNN_by_diet.pdf'),w=14,h=6.5)
  rm(obj,meta,emb,score,xy); invisible(gc())
}
readr::write_csv(dplyr::bind_rows(audit),file.path(out,'01_parent_gate_counts_by_cohort_and_diet.csv'))
readr::write_csv(dplyr::bind_rows(selected),file.path(out,'01_selected_nuclei_with_frozen_scores.csv.gz'))

# WT selection is an effect-size threshold in one pooled comparison, not a DAR test.
summary <- regions |>
 dplyr::count(WT_direction,name='n') |>
 dplyr::mutate(percent=100*n/sum(n),definition='WT C0/W0 top15; evaluable; abs pooled log2FC >=0.5')
readr::write_csv(summary,file.path(out,'02_WT_622_direction_summary.csv'))
chr <- regions |>
 dplyr::count(chrom,WT_direction,name='n') |>
 tidyr::complete(chrom,WT_direction,fill=list(n=0)) |>
 dplyr::mutate(chrom=factor(chrom,levels=paste0('chr',c(1:19,'X','Y'))))
readr::write_csv(chr,file.path(out,'02_WT_622_chromosome_counts.csv'))
universe <- read(file.path(l,'02_evaluable_WT_universe_sequence_covariates.csv.gz'))
require_cols(universe,c('feature','chrom'),'08L evaluable universe')
stopifnot(nrow(universe)==2028L,all(regions$feature %in% universe$feature))
chr_denom <- universe |> dplyr::count(chrom,name='n_evaluable') |>
  dplyr::left_join(regions |> dplyr::count(chrom,name='n_candidate'),by='chrom',relationship='one-to-one') |>
  dplyr::mutate(n_candidate=tidyr::replace_na(n_candidate,0L),
    fraction=n_candidate/n_evaluable,chrom=factor(chrom,levels=paste0('chr',c(1:19,'X','Y'))))
readr::write_csv(chr_denom,file.path(out,'02_WT_chromosome_evaluable_denominators.csv'))
p <- ggplot(chr,aes(chrom,n,fill=WT_direction)) + geom_col(width=.8) +
 scale_fill_manual(values=c(HFD_open='#E69F00',HFD_closed='#5A3E9A'),name='WT HFD effect') +
 scale_y_continuous(expand=expansion(mult=c(0,.08))) +
 labs(title='622 WT ISC-enriched candidate regions by chromosome',
 subtitle=paste0(sum(regions$WT_direction=='HFD_open'),' HFD-opening; ',sum(regions$WT_direction=='HFD_closed'),' HFD-closing'),
 caption='Counts reflect screened WT consensus peaks; chromosome size and evaluable peak counts differ. No biological replicate P values.',
 x=NULL,y='Candidate regions') + theme_classic(base_size=12) +
 theme(axis.text.x=element_text(angle=50,hjust=1),plot.caption=element_text(hjust=0))
save(p,'03_WT_622_regions_by_chromosome.pdf',w=13,h=7)
p_rate <- ggplot(chr_denom,aes(chrom,fraction)) + geom_col(fill='#7B42A8') +
 geom_text(aes(label=paste0(n_candidate,'/',n_evaluable)),vjust=-.25,size=2.7) +
 scale_y_continuous(labels=scales::percent_format(accuracy=1),expand=expansion(mult=c(0,.13))) +
 labs(title='Fraction of evaluable WT peaks meeting the candidate rule by chromosome',
 subtitle='Numerator = frozen candidate regions; denominator = all evaluable WT C0/W0 top-15 peaks',
 caption='Descriptive screening fraction, not a chromosome enrichment test or a DAR P value.',
 x=NULL,y='Candidates / evaluable peaks') + theme_classic(base_size=11) +
 theme(axis.text.x=element_text(angle=50,hjust=1),plot.caption=element_text(hjust=0))
save(p_rate,'03b_WT_candidate_fraction_by_chromosome.pdf',w=13,h=7)
# The genotype bars recount these *same* 622 coordinates; they are not
# genotype-wide DAR counts or a test of interaction.
cohort_622 <- dplyr::bind_rows(lapply(cohorts,function(cohort) {
  if(cohort=='wt') return(tibble(cohort=cohort,gate='WT C0/W0 top15',
    n_frozen=622L,n_evaluable=622L,n_HFD_open=sum(regions$WT_direction=='HFD_open'),
    n_HFD_closed=sum(regions$WT_direction=='HFD_closed')))
  dplyr::bind_rows(lapply(c('top15','all'),function(gate) {
    key <- paste0(cohort,'_annotated_union_',gate,'_')
    evaluable <- regions[[paste0(key,'evaluable')]]
    clear <- regions[[paste0(key,'clear_effect')]]
    effect <- regions[[paste0(key,'log2FC_HFD_vs_CON')]]
    stopifnot(!anyNA(evaluable),!anyNA(clear))
    tibble(cohort=cohort,gate=if(gate=='top15') 'Within-parent top15' else 'Whole approved parents',
      n_frozen=622L,n_evaluable=sum(evaluable),
      n_HFD_open=sum(evaluable & clear & effect>0),
      n_HFD_closed=sum(evaluable & clear & effect<0))
  }))
}))
cohort_622$panel <- factor(paste(toupper(cohort_622$cohort),cohort_622$gate,'\n',
  cohort_622$n_evaluable,'/622 evaluable'),levels=paste(toupper(cohort_622$cohort),
    cohort_622$gate,'\n',cohort_622$n_evaluable,'/622 evaluable'))
readr::write_csv(cohort_622,file.path(out,'03_exact_WT_622_interval_cohort_descriptive_counts.csv'))
p_cross <- cohort_622 |>
  tidyr::pivot_longer(c(n_HFD_open,n_HFD_closed),names_to='direction',values_to='n') |>
  ggplot(aes(panel,n,fill=direction)) + geom_col() +
  geom_text(aes(label=ifelse(n>0,n,'')),position=position_stack(vjust=.5),color='white',size=3.3) +
  scale_fill_manual(values=c(n_HFD_open='#E69F00',n_HFD_closed='#5A3E9A'),
    labels=c('HFD opening','HFD closing')) +
  labs(title='Pooled effects at the same 622 WT candidate intervals',
   subtitle='Top-15% gates shown beside full knockout parent unions; labels include evaluable denominators',
   caption='WT top-15 selected these loci. Genotype parent identities differ; thresholded counts change with gate size and coverage. No biological replicate P value or interaction test.',
   x=NULL,y='Exact WT intervals with clear pooled effect') +
  theme_classic(base_size=11) + theme(legend.title=element_blank(),plot.caption=element_text(hjust=0),
    axis.text.x=element_text(size=9))
save(p_cross,'03c_same_622_intervals_by_cohort_and_gate.pdf',w=13,h=7)

# Predeclare examples from the reviewed higher-count C0-supported set.
examples <- c('chr9-65121705-65122834','chr4-40142794-40143817','chr15-79166451-79167322')
stopifnot(all(examples %in% regions$feature),all(regions$C0_supported[match(examples,regions$feature)]),
  all(regions$count_supported_40[match(examples,regions$feature)]))
fx <- dplyr::bind_rows(lapply(examples,function(f) {
 r <- regions[match(f,regions$feature),]
 tibble(feature=f,cohort=c('wt','ppar','il17'),
   contrast=c('C0_all','annotated_union_all','annotated_union_all'),
   effect=c(r$C0_all_log2FC_HFD_vs_CON,r$ppar_annotated_union_all_log2FC_HFD_vs_CON,
     r$il17_annotated_union_all_log2FC_HFD_vs_CON),
   evaluable=c(TRUE,r$ppar_annotated_union_all_evaluable,r$il17_annotated_union_all_evaluable),
   CON_counts=c(NA_real_,r$ppar_annotated_union_all_CON_counts,r$il17_annotated_union_all_CON_counts),
   HFD_counts=c(NA_real_,r$ppar_annotated_union_all_HFD_counts,r$il17_annotated_union_all_HFD_counts))
}))
readr::write_csv(fx,file.path(out,'03_prespecified_locus_pooled_effects.csv'))
p <- ggplot(fx |> dplyr::filter(evaluable),aes(cohort,effect,fill=cohort)) +
 geom_hline(yintercept=0,color='grey60') + geom_col(width=.63) +
 geom_text(aes(label=sprintf('%+.2f',effect)),vjust=-.4,size=3.2) +
 facet_wrap(~feature,nrow=1) +
 scale_fill_manual(values=c(wt='#E69F00',ppar='#5A3E9A',il17='#397B89'),guide='none') +
 labs(title='Predeclared WT opening examples: pooled ATAC effects at exact loci',
 subtitle='WT focused C0 versus each knockout whole approved crypt-parent union',
 caption='Log2FC compares one pooled HFD library with one pooled CON library per genotype; no mouse-level inferential P value.',
 x=NULL,y='Pooled ATAC log2FC, HFD / CON') + theme_classic(base_size=11)
save(p,'04_three_loci_pooled_effects.pdf',w=13,h=6)

# Coverage normalization is by group cell count and sequencing depth in Signac.
# Per-locus common numeric ymax is derived from preliminary *normalized* tracks.
# No q95 per-panel clipping; all cohort/diet tracks share the same final y axis.
status <- list()
for (f in if (RUN_COVERAGE) examples else character()) {
 rr <- regions[match(f,regions$feature),]
 gr <- GenomicRanges::GRanges(rr$chrom,IRanges::IRanges(rr$start_1based,rr$end_1based))
 span <- GenomicRanges::GRanges(rr$chrom,IRanges::IRanges(max(1L,rr$start_1based-4000L),rr$end_1based+4000L))
 prelim <- list(); cells_n <- list()
 for (cohort in cohorts) {
   message('Coverage: ',f,' ',cohort)
   obj <- readRDS(file.path(base,cohort,'06_annotation','multiome_annotated.rds'))
   cells <- selected[[cohort]]$cell_barcode
   if(cohort=='wt') {
     wt_cells <- read(gate_files[['wt']])
     cells <- wt_cells$cell_barcode[as.character(wt_cells$focused_cluster)=='0']
   } else {
     scores <- read(gate_files[[cohort]])
     cells <- scores$cell_barcode[as.character(scores$wnn_cluster) %in% parents[[cohort]]]
   }
   m <- obj[[]]; stopifnot(all(cells %in% rownames(m)))
   nn <- table(factor(as.character(m[cells,'diet']),levels=c('CON','HFD')))
   cells_n[[cohort]] <- nn
   if(length(Signac::Fragments(obj[['ATAC']]))==0L || any(nn<30L)) {
     status[[length(status)+1L]] <- tibble(feature=f,cohort=cohort,status='skipped_no_fragments_or_small_group',detail='')
     rm(obj); next
   }
   tryCatch({
     set.seed(1000L+match(f,examples)*10L+match(cohort,cohorts))
     p0 <- Signac::CoveragePlot(obj,region=span,assay='ATAC',cells=cells,group.by='diet',
       annotation=FALSE,peaks=FALSE,links=FALSE,tile=FALSE,region.highlight=gr,
       window=100L,max.downsample=3000L,downsample.rate=1,ymax=NULL)
     # Signac's first component is the normalized coverage ggplot.
     track <- p0[[1]]
     y <- if(inherits(track,'ggplot')) track$data$score else NULL
     if(is.null(y) || !length(y)) {
       built <- ggplot2::ggplot_build(track)$data
       y <- unlist(lapply(built,function(layer) if('y' %in% names(layer)) layer$y else numeric()),use.names=FALSE)
     }
     y <- y[is.finite(y)]
     if(!length(y)) stop('Cannot read normalized coverage track for common ymax')
     prelim[[cohort]] <- list(cells=cells,max=max(y),n=nn)
   },error=function(e) {
     status[[length(status)+1L]] <<- tibble(feature=f,cohort=cohort,status='preflight_failed',detail=conditionMessage(e))
   })
   rm(obj); invisible(gc())
 }
 if(length(prelim)!=length(cohorts)) {
   message('Skipping cross-cohort coverage for ',f,': incomplete preliminary tracks')
   rm(prelim); invisible(gc()); next
 }
 shared_ymax <- max(1e-3,max(vapply(prelim,`[[`,numeric(1),'max'))*1.04)
 pages <- list()
 for(cohort in cohorts) {
   rec <- prelim[[cohort]]
   tryCatch({
     obj <- readRDS(file.path(base,cohort,'06_annotation','multiome_annotated.rds'))
     set.seed(1000L+match(f,examples)*10L+match(cohort,cohorts))
     p <- Signac::CoveragePlot(obj,region=span,assay='ATAC',cells=rec$cells,
       group.by='diet',annotation=FALSE,peaks=FALSE,links=FALSE,tile=FALSE,
       region.highlight=gr,window=100L,max.downsample=3000L,downsample.rate=1,
       ymax=shared_ymax) +
       patchwork::plot_annotation(title=paste(toupper(cohort), 'CON',rec$n[['CON']],'/ HFD',rec$n[['HFD']]),
          subtitle=paste('Fixed normalized y max',signif(shared_ymax,3)))
     pages[[cohort]] <- p
     status[[length(status)+1L]] <- tibble(feature=f,cohort=cohort,status='saved',detail='')
   },error=function(e) {
     status[[length(status)+1L]] <<- tibble(feature=f,cohort=cohort,status='plot_failed',detail=conditionMessage(e))
   })
   if(exists('obj',inherits=FALSE)) rm(obj)
   invisible(gc())
 }
 if(length(pages)==length(cohorts)) {
   assembled <- patchwork::wrap_plots(pages,ncol=1) +
     patchwork::plot_annotation(title=paste('Same GRCm39 locus:',f),
      subtitle='WT C0; PPAR and IL17 full approved ISC parent unions. Highlight = frozen WT interval.',
      caption='Normalized Tn5 insertion frequency per group (cells × mean depth); shared y limit for this locus. One pooled library per diet.')
   save(assembled,paste0('05_',f,'_three_cohort_common_scale_coverage.pdf'),w=14,h=14)
 }
 rm(prelim,pages); invisible(gc())
}
coverage_status <- if(length(status)) dplyr::bind_rows(status) else
  tibble(feature=character(),cohort=character(),status=character(),detail=character())
readr::write_csv(coverage_status,file.path(out,'05_coverage_status.csv'))
writeLines(c('Meeting figure script v2 completed.',
  '622 WT candidates are descriptive pooled effects, not statistically significant DARs.',
  if (RUN_COVERAGE) 'Coverage status must be checked before presenting any locus.' else 'Coverage disabled by --no-coverage.',
  paste('Figures:',plots)),file.path(out,'08M_COMPLETE.txt'))
message('Finished: ',out)
