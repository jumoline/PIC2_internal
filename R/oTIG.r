#' Function to infer relations between terms based on shared members
#'
#' \code{oTIG} is supposed to infer relations between terms based on shared members. It returns an object of class "igraph".
#'
#' @param data a tibble with two columns 'name' and 'members' (each is a vector containing members separated by ', '). Alternatively a data frame (matrix) with rows for members (row names) and columns for terms (column names)
#' @param method how to infer relations between two terms. It can be 'jaccard' (minimum spanning tree minimising the jaccard distance), 'empirical' (the number of the shared members is no less than 'empirical.cutoff' of two-end terms), and 'hybrid' (the union of results from 'jaccard' and 'empirical')
#' @param empirical.cutoff the empirical cutoff used to declare the relation between two terms (the number of the shared members is no less than 'empirical.cutoff' of two-end terms). By default, it is set to 0.5. Only works when the method is 'empirical' or 'hybrid' above
#' @return an object of class "igraph" with node attributes ('name', 'members', 'n_members', 'xcoord' and 'ycoord') and edge attributes ('shared', 'n_shared' and 'weight')
#' @note none
#' @export
#' @seealso \code{\link{oTIG}}
#' @include oTIG.r
#' @examples
#' set.seed(825)
#' mat <- ifelse(matrix(rnorm(1000),200,5)>0, 1, 0)
#' rownames(mat) <- paste0('R', seq(nrow(mat)))
#' colnames(mat) <- paste0('C', seq(ncol(mat)))
#' 
#' \dontrun{
#' tig <- oTIG(mat, method='hybrid')
#' gp <- oGGnetwork(tig, node.label='name', node.label.size=2.5, node.label.color='darkblue', node.label.force=1, node.xcoord='xcoord', node.ycoord='ycoord', colormap="cyan4-cyan4", edge.size='weight', node.shape=18, node.size='n_members', node.size.title="Num of \ngenes", node.size.range=c(2,5), edge.color='cyan4', edge.color.alpha=0.3, edge.curve=0,edge.arrow.gap=0.01)
#' }

oTIG <- function(data, method=c("hybrid","jaccard","empirical"), empirical.cutoff=0.5)
{
    
    method <- match.arg(method)
    
	name <- members <- value <- flag <- from_members <- to_members <- n_from <- n_to <- n_shared <- shared <- from <- to <- jaccard <- small <- large <- empirical <- NULL
    
    df <- NULL
	if(is(data,'data.frame') | is(data,'matrix')){
		if(all(c('name','members') %in% colnames(data))){
			df <- data %>% dplyr::select(name, members)
		}else{
			data %>% tibble::as_tibble(rownames='members') %>% tidyr::pivot_longer(cols=-members,names_to='name', values_to='value') %>% dplyr::filter(value==1) %>% dplyr::select(-value) %>% dplyr::group_by(name) %>% dplyr::summarise(members=stringr::str_c(members,collapse=', ')) %>% dplyr::ungroup() -> df
		}
	}
	
	##########################
	##########################
	if(nrow(df)==0){
		return(NULL)
	}
	##########################
	##########################
    
    ## df_nodes
	df %>% dplyr::mutate(n_members=purrr::map_int(members,~stringr::str_split(.x,', ',simplify=TRUE) %>% length())) -> df_nodes
	
	## big_projected
	df_nodes %>% tidyr::separate_rows(members, sep=", ") %>% dplyr::mutate(flag=1) %>% dplyr::select(members,name,flag) %>% tidyr::pivot_wider(names_from=members, values_from=flag, values_fill=list(flag=0)) %>% tibble::column_to_rownames("name") %>% oBicreate() %>% oBiproject -> big_projected
	
	# df: from  to    n_from  n_to n_shared shared                            exp  dist
	df_e_big <- oIG2TB(big_projected,'edges')
	df_e_big %>% dplyr::inner_join(df_nodes %>% dplyr::transmute(from=name,from_members=members), by="from") %>% dplyr::inner_join(df_nodes %>% transmute(to=name,to_members=members), by="to") %>% dplyr::mutate(shared=purrr::map2_chr(from_members,to_members, function(x,y){
		intersect(stringr::str_split(x,', ',simplify=TRUE), stringr::str_split(y,', ',simplify=TRUE)) %>% stringr::str_c(collapse=', ')
	})) %>% dplyr::mutate(n_from=purrr::map_int(from_members,~stringr::str_split(.x,', ',simplify=TRUE) %>% length())) %>%
	dplyr::mutate(n_to=purrr::map_int(to_members,~str_split(.x,', ',simplify=TRUE) %>% length())) %>% 
	dplyr::mutate(n_shared=purrr::map_int(shared,~str_split(.x,', ',simplify=TRUE) %>% length())) -> df
	
	## jaccard distance and empirical
	df %>% select(from,to,n_from,n_to,n_shared,shared) %>% mutate(empirical=empirical.cutoff*min(n_from,n_to)) %>% mutate(jaccard=1-n_shared/(n_from+n_to-n_shared)) -> df
	
	## jaccard: df_e_mst
	### ig_tmp (updated from big_projected)
	ig_tmp <- oTB2IG(edges=df %>% select(from,to,jaccard,n_shared,shared), nodes=df_nodes)
	### mst minimising jaccard distance
	mst <- igraph::minimum.spanning.tree(ig_tmp, weights=E(ig_tmp)$jaccard)
	### important: always from < to
	df_e_jaccard <- oIG2TB(mst,'edges') %>% dplyr::mutate(small=ifelse(from<=to,from,to), large=ifelse(from>to,from,to)) %>% dplyr::transmute(from=small, to=large, shared=shared, n_shared=n_shared) %>% dplyr::arrange(from)

	## empirical: df_e_empirical
	df_e_empirical <- df %>% dplyr::filter(n_shared>=empirical) %>% dplyr::select(from, to, shared, n_shared)
	
	## tig: a network between terms
	if(method=='hybrid'){
		df_edges <- rbind(df_e_jaccard, df_e_empirical) %>% dplyr::distinct()
		tig <- oTB2IG(edges=df_edges, nodes=df_nodes)
	}else if(method=='jaccard'){
		df_edges <- df_e_jaccard
		tig <- oTB2IG(edges=df_edges, nodes=df_nodes)
	}else if(method=='empirical'){
		df_edges <- df_e_empirical
		tig <- oTB2IG(edges=df_edges, nodes=df_nodes)
	}
	
	## append the edge attribute 'weight' for 'n_shared' rescaled into [1,2]
	w <- E(tig)$n_shared
	E(tig)$weight <- 1 + 2*(w - min(w))/(max(w) - min(w))
	tig %>% oLayout("graphlayouts.layout_with_stress") -> tig
	
	#gp <- oGGnetwork(tig, node.label='name', node.label.size=2.5, node.label.color='darkblue', node.label.force=1, node.xcoord='xcoord', node.ycoord='ycoord', colormap="cyan4-cyan4", edge.size='weight', node.shape=18, node.size='n_members', node.size.title="Num of \ngenes", node.size.range=c(2,5), edge.color='cyan4', edge.color.alpha=0.3, edge.curve=0,edge.arrow.gap=0.01)
	
	tig
}

