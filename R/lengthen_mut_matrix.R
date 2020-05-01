#' Lengthen mutation matrix
#' 
#' A mutation_matrix calculated on a GRangesList or GR object modified by 'split_muts_region()',
#' will contain a column per combination of sample and genomic region. In essence different regions
#' are treated as different samples. This function will transform the matrix, so that these regions
#' are instead treated as different mutation types. For example, instead of 'C[C>T]G', you might have
#' the feature 'C[C>T]G Promoter'. The number of rows in the matrix will thus be
#' multiplied by the number of regions.
#' After using 'split_muts_region()', use 'mut_matrix()' to get a mut_matrix that can be used
#' for this function.
#' The result can be plotted with plot_long_profile, but could also be used for NMF, refitting ect.
#' @param mut_matrix Mutation matrix
#' @return mut_matrix
#' 
#' @export
#' @importFrom magrittr %>% 
#'
#' @seealso \code{\link{lengthen_mut_matrix_single_sample}}
#' @family genomic_regions
#' @examples
#' 
#' ## See the 'split_muts_region()' and 'mut_matrix()' examples for how we obtained the
#' ## mutation matrix information:
#' mut_mat_split_region <- readRDS(system.file("states/mut_mat_data.rds",
#'                                 package="MutationalPatterns"))
#'                                 
#' long_mut_mat = lengthen_mut_matrix(mut_mat_split_region)
#' 
#' 
#' ## This also works on indels:
#' ## See the 'split_muts_region()' and 'count_indels_context()' examples for how we
#' ## obtained the indel counts:
#' indel_counts_split <- readRDS(system.file("states/blood_indels_counts_split_region.rds",
#'                                 package="MutationalPatterns"))
#' 
#' ## Transform the indel counts into a matrix.
#' indel_m_split = as.matrix(dplyr::select(indel_counts_split, -muttype, -muttype_sub))
#' rownames(indel_m_split) = stringr::str_c(indel_counts_split$muttype, 
#'                                         indel_counts_split$muttype_sub, sep = "_")
#' 
#' ## Lengthen the matrix
#' lengthen_mut_matrix(indel_m_split)
#' 
lengthen_mut_matrix = function(mut_matrix){
    col_names = colnames(mut_matrix)
    max_dots_in_name = col_names %>% 
        stringr::str_count("\\.") %>% 
        max()
    if (max_dots_in_name > 1){
        stop("The column names of the mutation matrix contain too many dots. 
             There should only be a dot in between the sample name and the type")
    }
    sample_names = stringr::str_remove(col_names, "\\..*")
    samples = unique(sample_names)
    mut_mat_l = purrr::map(samples, function(sample) mut_matrix[,sample == sample_names, drop = F])
    mut_matrix = purrr::map(mut_mat_l, lengthen_mut_matrix_single_sample) %>% 
        do.call(cbind, .)
    return(mut_matrix)
}

#' Lengthen mutation matrix for a single sample
#' 
#' This function is called by 'lengthen_mut_matrix()' to lengthen a single sample.
#' @param mut_matrix Mutation matrix
#' 
#' @importFrom magrittr %>% 
#' @seealso \code{\link{lengthen_mut_matrix}}
#' @family genomic_regions
#' @return mut_matrix single sample
lengthen_mut_matrix_single_sample = function(mut_matrix){
    col_names = colnames(mut_matrix)
    sample = col_names %>%
        stringr::str_remove("\\..*") %>% 
        unique()
    
    #Determine the new rownames of the mut_mat
    types = stringr::str_remove(col_names, ".*\\.")
    types_rownames = rep(types, each = nrow(mut_matrix))
    new_rownames = stringr::str_c(rownames(mut_matrix), "_", types_rownames)
    
    #Change shape of mutation matrix to a single column
    dim(mut_matrix) = c(length(new_rownames), 1)
    rownames(mut_matrix) = new_rownames
    colnames(mut_matrix) = sample
    return(mut_matrix)
}