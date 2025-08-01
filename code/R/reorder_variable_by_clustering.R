reorder_variable_by_clustering <- function(data, variable_to_order, response_variable, count_variable, value_fill = NULL, return_dendro = FALSE){
  
  if(is.null(value_fill)){
    wide <- data[c(variable_to_order, response_variable, count_variable)] %>%
      pivot_wider(names_from = one_of(response_variable), values_from = one_of(count_variable))
  }else{
    wide <- data[c(variable_to_order, response_variable, count_variable)] %>%
      pivot_wider(names_from = one_of(response_variable), values_from = one_of(count_variable), values_fill = value_fill)
  }
  
  dendrogram <- wide %>%
    dist() %>%
    hclust() %>%
    as.dendrogram() 
  order_indices <- dendrogram %>%
    order.dendrogram()
  
  data[variable_to_order] <- factor(data[[variable_to_order]], levels = wide[[variable_to_order]][order_indices])
  
  if(return_dendro){
    return(list(reord_df = data,
                dendrogram = dendrogram))
  }else{
    return(data)
  }
}

# Consider if you're then going to group by another group (nesting, essentially)
