#' Function to generate clades from DAG
#'
#' \code{oDAGclade} is supposed to generate clades from a direct acyclic graph (DAG; an ontology), that is, subgraphs induced by each of input nodes.
#'
#' @param ig an object of class "igraph" to represent DAG
#' @param nodes_query nodes in query
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @return 
#' NULL or an object of class "igraph" or a list of "igraph" objects
#' @note none
#' @export
#' @seealso \code{\link{oDAGclade}}
#' @include oDAGclade.r
#' @examples
#' \dontrun{
#' ###########################################################
#' ig <- oRDS('ig.GOMF', placeholder=placeholder)
#' nodes_query <- ig %>% oIG2TB('nodes') %>% filter(distance==1) %>% pull(name)
#' clade <- oDAGclade(ig, nodes_query)
#' }

oDAGclade <- function(ig, nodes_query=NULL, verbose=TRUE)
{
    
    if(!is(ig,"igraph")){
        warnings("The function must apply to the 'igraph' object.\n")
        return(NULL)
    }
    
    if(is(nodes_query,"igraph.vs")){
        nodes_query <- nodes_query$name
    }
    if(is.null(nodes_query)){
    	# if nodes_query is not provided, all tips are used
    	message("Nodes in query are not provided.\n")
    	return(ig)
    }
    
    ## check nodes in query
    ind <- match(nodes_query, V(ig)$name)
    nodes_found <- nodes_query[!is.na(ind)]
    if(length(nodes_found)==0){
    	warnings("Nodes in query cannot be found in the input graph.\n")
        return(NULL)
    }else{
        nodes_query <- V(ig)[nodes_found]$name
    }
	
	node <- NULL
	
	tibble::tibble(node=nodes_query) %>% mutate(ig=map(node, function(x){
		neighs.out <- igraph::neighborhood(ig, order=vcount(ig), nodes=x, mode="out")
		nodeInduced <- neighs.out %>% unlist %>% unique()
		subg <- igraph::induced.subgraph(ig, vids=nodeInduced)
		
		# updata distance
		if(!is.null(V(subg)$distance)){
			V(subg)$distance <- igraph::distances(subg, v=x, to=V(subg), mode="out") %>% as.numeric()
		}
		# updata namespace
		if(!is.null(V(subg)$namespace)){
			V(subg)$namespace <- stringr::str_replace_all(x,"_"," ") %>% stringr::str_to_title()
		}
		subg
	})) -> df
	
    if(verbose){
    	message(sprintf("%d clades created\n", nrow(df)), appendLF=TRUE)
    }
	
	df %>% tibble::deframe() -> clade
	
	if(length(clade)==1){
		clade <- clade[[1]]
	}
    
    return(clade)
}
