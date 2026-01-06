#' Create Isoform Similarity Matrix from GTF
#'
#' Parses a GTF file and creates a similarity matrix for isoforms of a specified
#' gene based on shared exon structure. This matrix can be used with
#' \code{\link{multivariate_graph_reg}} for graph-regularized prediction.
#'
#' @param gtf_path character, path to GTF file (can be gzipped)
#' @param gene character, gene name (HGNC symbol) or Ensembl gene ID
#' @param method character, similarity method:
#'   \itemize{
#'     \item "jaccard": Jaccard index based on exon overlap (default)
#'     \item "overlap_coef": Overlap coefficient = overlap / min(length_i, length_j)
#'     \item "binary": 1 if any overlap, 0 otherwise
#'   }
#' @param transcript_ids character vector, optional subset of transcript IDs to include.
#'   If NULL, uses all transcripts for the gene.
#' @param min_transcripts int, minimum number of transcripts required (default 2)
#' @param verbose logical, print progress messages
#'
#' @return A list containing:
#'   \itemize{
#'     \item similarity_matrix: symmetric matrix of pairwise isoform similarities
#'     \item transcript_ids: character vector of transcript IDs (matrix row/col names)
#'     \item gene_id: matched gene identifier
#'     \item n_transcripts: number of transcripts
#'     \item n_exons: named vector of exon counts per transcript
#'   }
#'
#' @details
#' The function computes pairwise similarity between isoforms based on their
#' exon structure. Isoforms that share more exonic sequence are considered
#' more similar, reflecting the biological intuition that they likely share
#' more cis-regulatory effects.
#'
#' The Jaccard similarity is computed as:
#' \deqn{J(i,j) = \frac{|overlap(exons_i, exons_j)|}{|union(exons_i, exons_j)|}}
#'
#' where overlap and union are computed in terms of genomic base pairs.
#'
#' @examples
#' \dontrun{
#' # Create similarity matrix for BRCA1 isoforms
#' sim <- similarity_from_gtf("gencode.v45.annotation.gtf.gz", "BRCA1")
#'
#' # Use with graph-regularized prediction
#' model <- multivariate_graph_reg(X, Y, similarity_matrix = sim$similarity_matrix)
#' }
#'
#' @importFrom rtracklayer import
#' @importFrom GenomicRanges findOverlaps pintersect reduce
#' @importFrom S4Vectors queryHits subjectHits mcols
#' @export
similarity_from_gtf <- function(gtf_path,
                                 gene,
                                 method = c("jaccard", "overlap_coef", "binary"),
                                 transcript_ids = NULL,
                                 min_transcripts = 2,
                                 verbose = TRUE) {

  method <- match.arg(method)

  # Check for required package

  if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    stop("Package 'rtracklayer' is required. Install with: BiocManager::install('rtracklayer')")
  }

  # Import GTF
  if (verbose) message("Reading GTF file...")
  gtf <- rtracklayer::import(gtf_path)

  # Find the gene
  gene_match <- .find_gene_in_gtf(gtf, gene)

  if (is.null(gene_match)) {
    stop("Gene '", gene, "' not found in GTF. ",
         "Try using gene_name (e.g., 'BRCA1') or gene_id (e.g., 'ENSG00000012048').")
  }

  gene_gtf <- gene_match$matches
  matched_id <- gene_match$matched_id

  if (verbose) message("Found gene: ", matched_id, " (matched by ", gene_match$matched_by, ")")


  # Extract exons
  exons <- gene_gtf[gene_gtf$type == "exon"]

  if (length(exons) == 0) {
    stop("No exons found for gene ", matched_id)
  }

  # Get transcript IDs
  all_transcripts <- unique(exons$transcript_id)
  all_transcripts <- all_transcripts[!is.na(all_transcripts)]

  # Subset if requested
  if (!is.null(transcript_ids)) {
    # Try exact match first
    matched <- transcript_ids[transcript_ids %in% all_transcripts]

    # Try without version numbers if no exact matches
    if (length(matched) == 0) {
      tx_no_version <- sub("\\..*$", "", all_transcripts)
      req_no_version <- sub("\\..*$", "", transcript_ids)
      match_idx <- which(tx_no_version %in% req_no_version)
      matched <- all_transcripts[match_idx]
    }

    if (length(matched) == 0) {
      stop("None of the requested transcript_ids found in GTF for this gene")
    }

    all_transcripts <- matched
    if (verbose) message("Using ", length(all_transcripts), " of ", length(transcript_ids), " requested transcripts")
  }

  if (length(all_transcripts) < min_transcripts) {
    stop("Gene ", matched_id, " has only ", length(all_transcripts),
         " transcripts (minimum ", min_transcripts, " required)")
  }

  if (verbose) message("Computing similarity for ", length(all_transcripts), " transcripts...")

  # Compute similarity matrix
  n <- length(all_transcripts)
  sim_matrix <- matrix(0, nrow = n, ncol = n,
                       dimnames = list(all_transcripts, all_transcripts))

  # Count exons per transcript
  n_exons <- sapply(all_transcripts, function(tx) {
    sum(exons$transcript_id == tx, na.rm = TRUE)
  })

  for (i in 1:n) {
    sim_matrix[i, i] <- 1  # Self-similarity

    if (i < n) {
      for (j in (i + 1):n) {
        exons_i <- exons[exons$transcript_id == all_transcripts[i]]
        exons_j <- exons[exons$transcript_id == all_transcripts[j]]

        sim <- .compute_exon_similarity(exons_i, exons_j, method)

        sim_matrix[i, j] <- sim
        sim_matrix[j, i] <- sim
      }
    }
  }

  if (verbose) {
    mean_sim <- mean(sim_matrix[upper.tri(sim_matrix)])
    message("Mean pairwise similarity: ", round(mean_sim, 3))
  }

  return(list(
    similarity_matrix = sim_matrix,
    transcript_ids = all_transcripts,
    gene_id = matched_id,
    n_transcripts = n,
    n_exons = n_exons,
    method = method
  ))
}


#' Find gene in GTF with flexible matching
#' @keywords internal
.find_gene_in_gtf <- function(gtf, target_gene) {
  mcols_gtf <- S4Vectors::mcols(gtf)

  # Try gene_name first
  if ("gene_name" %in% colnames(mcols_gtf)) {
    matches <- gtf[mcols_gtf$gene_name == target_gene & !is.na(mcols_gtf$gene_name)]
    if (length(matches) > 0) {
      return(list(matches = matches, matched_by = "gene_name", matched_id = target_gene))
    }
  }

  # Try gene_id
  if ("gene_id" %in% colnames(mcols_gtf)) {
    matches <- gtf[mcols_gtf$gene_id == target_gene & !is.na(mcols_gtf$gene_id)]
    if (length(matches) > 0) {
      return(list(matches = matches, matched_by = "gene_id", matched_id = target_gene))
    }
  }

  # Try gene_id without version
  if ("gene_id" %in% colnames(mcols_gtf)) {
    gene_ids_no_version <- sub("\\..*$", "", mcols_gtf$gene_id)
    target_no_version <- sub("\\..*$", "", target_gene)

    matches <- gtf[gene_ids_no_version == target_no_version & !is.na(mcols_gtf$gene_id)]
    if (length(matches) > 0) {
      matched_full_id <- unique(mcols_gtf$gene_id[gene_ids_no_version == target_no_version])[1]
      return(list(matches = matches, matched_by = "gene_id_no_version", matched_id = matched_full_id))
    }
  }

  return(NULL)
}


#' Compute similarity between two sets of exons
#' @keywords internal
.compute_exon_similarity <- function(exons_i, exons_j, method) {
  overlaps <- GenomicRanges::findOverlaps(exons_i, exons_j)

  if (length(overlaps) == 0) {
    return(0)
  }

  # Compute overlap width
  overlap_ranges <- GenomicRanges::pintersect(
    exons_i[S4Vectors::queryHits(overlaps)],
    exons_j[S4Vectors::subjectHits(overlaps)]
  )
  overlap_width <- sum(BiocGenerics::width(overlap_ranges))

  if (method == "binary") {
    return(1)
  }

  # Compute total widths
  width_i <- sum(BiocGenerics::width(exons_i))
  width_j <- sum(BiocGenerics::width(exons_j))

  if (method == "overlap_coef") {
    # Overlap coefficient: overlap / min(width_i, width_j)
    return(overlap_width / min(width_i, width_j))
  }

  # Jaccard: overlap / union
  union_ranges <- GenomicRanges::reduce(c(exons_i, exons_j))
  union_width <- sum(BiocGenerics::width(union_ranges))

  return(overlap_width / union_width)
}


#' Create similarity matrix for a subset of transcripts
#'
#' Given a pre-computed similarity result and a subset of transcript IDs,
#' returns the corresponding submatrix. Useful when Y matrix has fewer
#' transcripts than the full gene.
#'
#' @param sim_result result from \code{\link{similarity_from_gtf}}
#' @param transcript_ids character vector of transcript IDs to subset
#' @param match_without_version logical, try matching without version numbers
#'
#' @return subsetted similarity matrix
#'
#' @export
subset_similarity <- function(sim_result, transcript_ids, match_without_version = TRUE) {
  available <- rownames(sim_result$similarity_matrix)

  # Try exact match
  matched <- transcript_ids[transcript_ids %in% available]

  # Try without version if needed
  if (length(matched) < length(transcript_ids) && match_without_version) {
    unmatched <- transcript_ids[!transcript_ids %in% available]
    available_no_ver <- sub("\\..*$", "", available)
    unmatched_no_ver <- sub("\\..*$", "", unmatched)

    for (i in seq_along(unmatched)) {
      idx <- which(available_no_ver == unmatched_no_ver[i])
      if (length(idx) > 0) {
        matched <- c(matched, available[idx[1]])
      }
    }
  }

  if (length(matched) == 0) {
    stop("No matching transcripts found in similarity matrix")
  }

  if (length(matched) < length(transcript_ids)) {
    warning(length(transcript_ids) - length(matched), " transcript(s) not found in similarity matrix")
  }

  return(sim_result$similarity_matrix[matched, matched, drop = FALSE])
}
