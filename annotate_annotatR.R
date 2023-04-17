#!/usr/bin/env Rs#cript

# Loading libraries
library(annotatr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(dplyr)
library(GenomicRanges)
library(ggplot2)

# Data already loaded

# Making into GRange object
t_grange=makeGRangesFromDataFrame(t.cpg, 
                                     keep.extra.columns = TRUE,
                                     seqnames.field = "chr",
                                     start.field = "start",
                                     end.field = "end",
                                     strand.field = "strand")

n_grange=makeGRangesFromDataFrame(n.cpg, 
                                     keep.extra.columns = TRUE,
                                     seqnames.field = "chr",
                                     start.field = "start",
                                     end.field = "end",
                                     strand.field = "strand")

#Creating annotation vector
annotate_1=c('hg38_cpgs', 'hg38_genes_introns', 'hg38_genes_cds', 'hg38_genes_exons', 'hg38_genes_promoters', 'hg38_basicgenes', 'hg38_genes_intergenic', 'hg38_cpg_islands', 'hg38_cpg_shelves','hg38_cpg_shores', 'hg38_enhancers_fantom')

# build_annotation() creates a GRange object of all the annotations combined
annotate_gr=build_annotations(genome = 'hg38', annotations = annotate_1)

#Intersecting the regions we read with annotations
annotated_t.cpg=annotate_regions(regions = t_grange, 
                                    annotations = annotate_gr, 
                                    ignore.strand = TRUE,
                                    quiet = FALSE)

annotated_n.cpg=annotate_regions(regions = n_grange,
                                    annotations = annotate_gr,
                                    ignore.strand = TRUE,
                                    quiet = FALSE)

# Converting the output into data frame
df_annotated_t.cpg=data.frame(annotated_t.cpg)
df_annotated_n.cpg=data.frame(annotated_n.cpg)

# Removing NA's from the data frame
df_annotated_t.cpg_2=na.omit(df_annotated_t.cpg)
df_annotated_n.cpg_2=na.omit(df_annotated_n.cpg)


# Summarising Over Annotations
anno_t_summary=summarize_annotations(annotated_regions = annotated_t.cpg,
                                        quiet = TRUE)
anno_n_summary=summarize_annotations(annotated_regions = annotated_n.cpg,
                                        quiet = TRUE)

# Plotting annotation data
#annots_order=c('hg38_cpgs', 'hg38_genes_introns', 'hg38_genes_cds', 'hg38_genes_exons', 'hg38_genes_promoters', 'hg38_basicgenes', 'hg38_genes_intergenic', 'hg38_cpg_islands', 'hg38_cpg_inter', 'hg38_enhancers_fantom')

#t.cpg
#plot_1=plot_annotation(annotated_regions = annotated_t.cpg, annotation_order = annots_order, plot_title = "Annotation plot of t.cpg", x_label="Gene Annotations", y_label="Count")

#n.cpg
#plot_2=plot_annotation(annotated_regions = annotated_n.cpg, annotation_order = annots_order, plot_title = "Annotation plot of n.cpg", x_label="Gene Annotations", y_label="Count")

#Plotting data with ggplot
# t.cpg annotation plot
#plot_1=ggplot(df_annotated_t.cpg_2, aes(x=annot.type))+geom_histogram(stat="count", fill="red")+labs(x="Annotation Type", y="Count", title="Annotation plot of t.cpg")+scale_x_discrete(guide=guide_axis(angle=45))+ylim(0e+00, 4e+07)+theme(axis.text.x=element_text(size=16), axis.text.y=element_text(size=16))
#ggsave(filename="annotation_plot_t.cpg_modified.pdf",width=10, height=8, units="in", dpi=300)

#n.cpg annotation plot
#plot_7=ggplot(df_annotated_n.cpg_2, aes(x=annot.type))+geom_histogram(stat="count", fill="blue")+labs(x="Annotation Type", y="Count", title="Annotation plot of n.cpg")+scale_x_discrete(guide=guide_axis(angle=45))+ylim(0e+00, 4e+07)+theme(axis.text.x=element_text(size=16), axis.text.y=element_text(size=16))
#ggsave(filename="annotation_plot_n.cpg_modified.pdf", width=10, height=8, units="in", dpi=300)

