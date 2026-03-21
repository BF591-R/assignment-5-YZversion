library('tidyverse')
library('SummarizedExperiment')
library('DESeq2')
library('biomaRt')
library('testthat')
library('fgsea')

#' Function to generate a SummarizedExperiment object with counts and coldata
#' to use in DESeq2
#'
#' @param csv_path (str): path to the file verse_counts.tsv
#' @param metafile (str): path to the metadata sample_metadata.csv
#' @param selected_times (list): list of sample timepoints to use
#' 
#'   
#' @return SummarizedExperiment object with subsetted counts matrix
#'   and sample data. Ensure that the timepoints column used as input 
#'   to the model design has 'vP0' set as the reference factor level. Your 
#'   colData dataframe should have columns named samplename and timepoint.
#' @export
#'
#' @examples se <- make_se('verse_counts.tsv', 'sample_metadata.csv', c('vP0', 'vAd'))
make_se <- function(counts_csv, metafile_csv, selected_times) {
    counts_df <- readr::read_tsv(counts_csv, show_col_types = FALSE)
    metadata <- readr::read_csv(metafile_csv, show_col_types = FALSE)

    sample_df <- metadata %>%
        dplyr::filter(timepoint %in% selected_times) %>%
        dplyr::select(samplename, timepoint) %>%
        dplyr::arrange(factor(timepoint, levels = selected_times), samplename)

    if (nrow(sample_df) == 0) {
        stop("No samples found for the requested timepoints.")
    }

    counts_subset <- counts_df %>%
        dplyr::select(gene, dplyr::all_of(sample_df$samplename))

    counts_mat <- counts_subset %>%
        tibble::column_to_rownames("gene") %>%
        as.matrix()

    col_data <- as.data.frame(sample_df)
    rownames(col_data) <- col_data$samplename
    col_data$timepoint <- factor(col_data$timepoint, levels = selected_times)
    if ("vP0" %in% levels(col_data$timepoint)) {
        col_data$timepoint <- stats::relevel(col_data$timepoint, ref = "vP0")
    }

    SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = counts_mat),
        colData = col_data
    )
}

#' Function that runs DESeq2 and returns a named list containing the DESeq2
#' results as a dataframe and the dds object returned by DESeq2
#'
#' @param se (obj): SummarizedExperiment object containing counts matrix and
#' coldata
#' @param design: the design formula to be used in DESeq2
#'
#' @return list with DESeqDataSet object after running DESeq2 and results from
#'   DESeq2 as a dataframe
#' @export
#'
#' @examples results <- return_deseq_res(se, ~ timepoint)
return_deseq_res <- function(se, design) {
    dds <- DESeq2::DESeqDataSet(se, design = design)
    dds <- DESeq2::DESeq(dds)
    res_df <- DESeq2::results(dds) %>% as.data.frame()

    list(
        results = res_df,
        dds = dds
    )
}

#' Function that takes the DESeq2 results dataframe, converts it to a tibble and
#' adds a column to denote plotting status in volcano plot. Column should denote
#' whether gene is either 1. Significant at padj < .10 and has a positive log
#' fold change, 2. Significant at padj < .10 and has a negative log fold change,
#' 3. Not significant at padj < .10. Have the values for these labels be UP,
#' DOWN, NS, respectively. The column should be named `volc_plot_status`. Ensure
#' that the column name for your rownames is called "genes". 
#'
#' @param deseq2_res (df): results from DESeq2 
#' @param padj_threshold (float): threshold for considering significance (padj)
#'
#' @return Tibble with all columns from DESeq2 results and one additional column
#'   labeling genes by significant and up-regulated, significant and
#'   downregulated, and not significant at padj < .10.
#'   
#' @export
#'
#' @examples labeled_results <- label_res(res, .10)
label_res <- function(deseq2_res, padj_threshold) {
    res_df <- as.data.frame(deseq2_res)
    if (!"genes" %in% names(res_df)) {
        res_df <- tibble::rownames_to_column(res_df, var = "genes")
    }

    res_tbl <- tibble::as_tibble(res_df)

    res_tbl %>%
        dplyr::mutate(
            volc_plot_status = dplyr::case_when(
                !is.na(padj) & padj < padj_threshold &
                    !is.na(log2FoldChange) & log2FoldChange > 0 ~ "UP",
                !is.na(padj) & padj < padj_threshold &
                    !is.na(log2FoldChange) & log2FoldChange < 0 ~ "DOWN",
                TRUE ~ "NS"
            )
        )
}

#' Function to plot the unadjusted p-values as a histogram
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one additional
#' column denoting status in volcano plot
#'
#' @return ggplot: a histogram of the raw p-values from the DESeq2 results
#' @export
#'
#' @examples pval_plot <- plot_pvals(labeled_results)
plot_pvals <- function(labeled_results) {
  ggplot2::ggplot(
    dplyr::filter(labeled_results, !is.na(pvalue)),
    ggplot2::aes(x = pvalue)
  ) +
    ggplot2::geom_histogram(binwidth = 0.025) +
    ggplot2::labs(
      x = "Raw p-value",
      y = "Count"
    )
}

#' Function to plot the log2foldchange from DESeq2 results in a histogram
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one additional
#' column denoting status in volcano plot
#' @param padj_threshold (float): threshold for considering significance (padj)
#'
#' @return ggplot: a histogram of log2FC values from genes significant at padj 
#' threshold of 0.1
#' @export
#'
#' @examples log2fc_plot <- plot_log2fc(labeled_results, .10)
plot_log2fc <- function(labeled_results, padj_threshold) {
  sig <- labeled_results %>%
    dplyr::filter(!is.na(padj), padj < padj_threshold, !is.na(log2FoldChange))
  
  ggplot2::ggplot(sig, ggplot2::aes(x = log2FoldChange)) +
    ggplot2::geom_histogram(binwidth = 0.25) +
    ggplot2::labs(
      x = "log2 fold change",
      y = "Count"
    )
}

#' Function to make scatter plot of normalized counts for top ten genes ranked
#' by ascending padj
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#' @param dds_obj (obj): The object returned by running DESeq (dds) containing
#' the updated DESeqDataSet object with test results
#' @param num_genes (int): Number of genes to plot
#'
#' @return ggplot: a scatter plot with the normalized counts for each sample for
#' each of the top ten genes ranked by ascending padj
#' @export
#'
#' @examples norm_counts_plot <- scatter_norm_counts(labeled_results, dds, 10)
scatter_norm_counts <- function(labeled_results, dds_obj, num_genes) {
  res_tbl <- labeled_results
  if (!"genes" %in% names(res_tbl)) {
    res_tbl <- tibble::rownames_to_column(as.data.frame(res_tbl), var = "genes")
  } else {
    res_tbl <- tibble::as_tibble(res_tbl)
  }
  
  top_genes <- res_tbl %>%
    dplyr::filter(!is.na(padj)) %>%
    dplyr::arrange(padj, dplyr::desc(abs(log2FoldChange))) %>%
    dplyr::slice_head(n = num_genes) %>%
    dplyr::pull(genes)
  
  norm_counts <- DESeq2::counts(dds_obj, normalized = TRUE)
  available_genes <- top_genes[top_genes %in% rownames(norm_counts)]
  
  if (length(available_genes) == 0) {
    stop("None of the requested genes are present in the normalized count matrix.")
  }
  
  counts_long <- norm_counts[available_genes, , drop = FALSE] %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "genes") %>%
    tidyr::pivot_longer(
      cols = -genes,
      names_to = "samplename",
      values_to = "normalized_count"
    )
  
  col_data <- as.data.frame(SummarizedExperiment::colData(dds_obj))
  
  plot_data <- counts_long %>%
    dplyr::left_join(col_data, by = "samplename") %>%
    dplyr::filter(!is.na(normalized_count), normalized_count > 0) %>%
    dplyr::mutate(
      genes = factor(genes, levels = available_genes)
    )
  
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = timepoint, y = normalized_count, color = timepoint)
  ) +
    ggplot2::geom_jitter(width = 0.15, height = 0, na.rm = TRUE) +
    ggplot2::facet_wrap(~genes, scales = "free_y") +
    ggplot2::scale_y_continuous(trans = "log10") +
    ggplot2::labs(
      x = "Timepoint",
      y = "Normalized counts"
    )
}
#' Function to generate volcano plot from DESeq2 results
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#'
#' @return ggplot: a scatterplot (volcano plot) that displays log2foldchange vs
#'   -log10(padj) and labeled by status
#' @export
#'
#' @examples volcano_plot <- plot_volcano(labeled_results)
#' 
plot_volcano <- function(labeled_results) {
  plot_data <- labeled_results %>%
    dplyr::filter(!is.na(log2FoldChange), !is.na(padj), padj > 0) %>%
    dplyr::mutate(
      neg_log10_padj = -log10(padj)
    )
  
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = log2FoldChange,
      y = neg_log10_padj,
      color = volc_plot_status
    )
  ) +
    ggplot2::geom_point(na.rm = TRUE) +
    ggplot2::labs(
      x = "log2 fold change",
      y = "-log10(padj)"
    )
}

#' Function to generate a named vector ranked by log2FC descending
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#' @param id2gene_path (str): Path to the file containing the mapping of
#' ensembl IDs to MGI symbols
#'
#' @return Named vector with gene symbols as names, and log2FoldChange as values
#' ranked in descending order
#' @export
#'
#' @examples rnk_list <- make_ranked_log2fc(labeled_results, 'data/id2gene.txt')

make_ranked_log2fc <- function(labeled_results, id2gene_path) {
    res_tbl <- labeled_results
    if (!"genes" %in% names(res_tbl)) {
        res_tbl <- tibble::rownames_to_column(as.data.frame(res_tbl), var = "genes")
    } else {
        res_tbl <- tibble::as_tibble(res_tbl)
    }

    id_map <- readr::read_tsv(
        id2gene_path,
        col_names = c("gene_id", "symbol"),
        show_col_types = FALSE
    )

    ranked_tbl <- res_tbl %>%
        dplyr::filter(!is.na(log2FoldChange)) %>%
        dplyr::left_join(id_map, by = c("genes" = "gene_id")) %>%
        dplyr::filter(!is.na(symbol)) %>%
        dplyr::group_by(symbol) %>%
        dplyr::slice_max(order_by = abs(log2FoldChange), n = 1, with_ties = FALSE) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(dplyr::desc(log2FoldChange))

    stats_vec <- ranked_tbl$log2FoldChange
    stats::setNames(stats_vec, ranked_tbl$symbol)
}

#' Function to run fgsea with arguments for min and max gene set size
#'
#' @param gmt_file_path (str): Path to the gene sets of interest in GMT format
#' @param rnk_list (named vector): Named vector generated previously with gene 
#' symbols and log2Fold Change values in descending order
#' @param min_size (int): Minimum number of genes in gene sets to be allowed
#' @param max_size (int): Maximum number of genes in gene sets to be allowed
#'
#' @return Tibble of results from running fgsea
#' @export
#'
#' @examples fgsea_results <- run_fgsea('data/m2.cp.v2023.1.Mm.symbols.gmt', rnk_list, 15, 500)
run_fgsea <- function(gmt_file_path, rnk_list, min_size, max_size) {
  if (length(rnk_list) == 0) {
    stop("Ranked list is empty; cannot run fgsea.")
  }
  
  gene_sets <- fgsea::gmtPathways(gmt_file_path)
  ordered_stats <- sort(rnk_list, decreasing = TRUE)
  
  fgsea::fgsea(
    pathways = gene_sets,
    stats = ordered_stats,
    minSize = min_size,
    maxSize = max_size,
    eps = 0,
    BPPARAM = BiocParallel::SerialParam()
  ) %>%
    data.table::as.data.table() %>%
    dplyr::arrange(dplyr::desc(NES)) %>%
    tibble::as_tibble()
}

#' Function to plot top ten positive NES and top ten negative NES pathways
#' in a barchart
#'
#' @param fgsea_results (tibble): the fgsea results in tibble format returned by
#'   the previous function
#' @param num_paths (int): the number of pathways for each direction (top or
#'   down) to include in the plot. Set this at 10.
#'
#' @return ggplot with a barchart showing the top twenty pathways ranked by positive
#' and negative NES
#' @export
#'
#' @examples fgsea_plot <- top_pathways(fgsea_results, 10)
top_pathways <- function(fgsea_results, num_paths){
    if (!"pathway" %in% names(fgsea_results) || !"NES" %in% names(fgsea_results)) {
        stop("fgsea_results must contain 'pathway' and 'NES' columns.")
    }

    pos <- fgsea_results %>%
        dplyr::filter(NES > 0) %>%
        dplyr::arrange(dplyr::desc(NES)) %>%
        dplyr::slice_head(n = num_paths)

    neg <- fgsea_results %>%
        dplyr::filter(NES < 0) %>%
        dplyr::arrange(NES) %>%
        dplyr::slice_head(n = num_paths)

    combined <- dplyr::bind_rows(pos, neg) %>%
        dplyr::mutate(
            direction = dplyr::if_else(NES >= 0, "Positive NES", "Negative NES"),
            pathway = forcats::fct_reorder(pathway, NES)
        )

    ggplot2::ggplot(combined, ggplot2::aes(x = pathway, y = NES, fill = direction)) +
        ggplot2::geom_col() +
        ggplot2::coord_flip() +
        ggplot2::scale_fill_manual(
            values = c("Positive NES" = "#009E73", "Negative NES" = "#CC79A7")
        ) +
        ggplot2::labs(
            x = "Pathway",
            y = "Normalized Enrichment Score",
            fill = "Direction",
            title = "Top pathways by NES"
        ) +
        ggplot2::theme_bw()
}

