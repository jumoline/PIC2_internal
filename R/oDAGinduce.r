#' Function to generate a DAG subgraph induced by input nodes
#'
#' \code{oDAGinduce} is supposed to produce a subgraph induced by input nodes, given a direct acyclic graph (DAG; an ontology). The input is a graph of "igraph", nodes, and the mode defining the paths to the root of DAG. The induced subgraph contains nodes and their ancestors along with the defined paths to the root of DAG.
#'
#' @param ig an object of class "igraph" to represent DAG
#' @param nodes_query nodes in query
#' @param path.mode the mode of paths induced by vertices/nodes with input annotation data. It can be "all_paths" for all possible paths to the root, "shortest_paths" for only one path to the root (for each node in query), "all_shortest_paths" for all shortest paths to the root (i.e. for each node, find all shortest paths with the equal lengths), "shortest_paths_epath" (indeed only edges in the shortest path are kept)
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @return 
#' NULL or an induced subgraph, an object of class "igraph".
#' @note For the mode "shortest_paths", the induced subgraph is the most concise, and thus informative for visualisation when there are many nodes in query, while the mode "all_paths" results in the complete subgraph.
#' @export
#' @seealso \code{\link{oDAGinduce}}
#' @include oDAGinduce.r
#' @examples
#' \dontrun{
#' ###########################################################
#' ig <- oRDS('ig.GOMF', placeholder=placeholder)
#' ig1 <- oDAGinduce(ig, path.mode="all_paths")
#' ig2 <- oDAGinduce(ig, path.mode="shortest_paths")
#' ig3 <- oDAGinduce(ig, path.mode="all_shortest_paths")
#' # very slow below
#' ig4 <- oDAGinduce(ig, path.mode="shortest_paths_epath")
#' }

oDAGinduce <- function(ig, nodes_query=NULL, path.mode=c("all_paths","shortest_paths","all_shortest_paths","shortest_paths_epath"), verbose=TRUE)
{
    
    startT <- Sys.time()
    if(verbose){
    	message(sprintf("Starting ... (at %s)\n", as.character(startT)), appendLF=TRUE)
    }
	
####################################################################################
	
    path.mode <- match.arg(path.mode)
    
    if(!is(ig,"igraph")){
        warnings("The function must apply to the 'igraph' object.\n")
        return(NULL)
    }
    
    if(is(nodes_query,"igraph.vs")){
        nodes_query <- nodes_query$name
    }
    if(is.null(nodes_query)){
    	# if nodes_query is not provided, all tips are used
    	message("Nodes in query are not provided; all tips are used instead.\n")
    	if(path.mode!='all_paths'){
			tips <- which(igraph::degree(ig,mode='out')==0)
			nodes_query <- V(ig)[tips]$name
		}else{
			# itself for 'all_paths'
			return(ig)
		}
    }
    
    ## check nodes in query
    ind <- match(nodes_query, V(ig)$name)
    nodes_found <- nodes_query[!is.na(ind)]
    if(length(nodes_found)==0){
    	warnings("Nodes in query cannot be found in the input graph.\n")
        return(NULL)
    }else{
        nodes_query <- V(ig)[nodes_found]
    }

    ## DAG being induced from nodes in query
    if(path.mode=="all_paths"){
        ## find all ancestors for any node
        neighs.in <- igraph::neighborhood(ig, order=igraph::vcount(ig), nodes=nodes_query$name, mode="in")
        nodeInduced <- neighs.in %>% unlist %>% unique()
        subg <- igraph::induced.subgraph(ig, vids=nodeInduced)
        
    }else if(path.mode=="all_shortest_paths"){
    	# calculates _all_ shortest paths between pairs of vertices
    	# creates a subgraph of a graph, containing specified vertices and all the edges among them
        root <- which(igraph::degree(ig,mode='in')==0)
        res <- igraph::all_shortest_paths(ig, from=root, to=nodes_query, mode='out')
        nodeInduced <- res$res %>% unlist %>% unique()
        subg <- igraph::induced.subgraph(ig, vids=nodeInduced)
        
    }else if(path.mode=="shortest_paths"){
    	# calculates a single shortest path
    	# creates a subgraph of a graph, containing specified vertices and all the edges among them
		root <- which(igraph::degree(ig,mode='in')==0)
		res <- igraph::shortest_paths(ig, from=root, to=nodes_query, mode='out', output="vpath")
		nodeInduced <- res$vpath %>% unlist %>% unique()
		subg <- igraph::induced.subgraph(ig, vids=nodeInduced)
    }else if(path.mode=="shortest_paths_epath"){
    	# calculates a single shortest path
    	# keeps edges in shorted pathes
		root <- which(igraph::degree(ig,mode='in')==0)
		res <- igraph::shortest_paths(ig, from=root, to=nodes_query, mode='out', output="epath")
		subg <- igraph::subgraph.edges(ig, eids=res$epath %>% unlist())
    }

####################################################################################
    endT <- Sys.time()
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    
    if(verbose){
    	message(sprintf("\nEnded (at %s)", as.character(endT)), appendLF=TRUE)
    	message(sprintf("Runtime in total: %d secs\n", runTime), appendLF=TRUE)
    }
    
    return(subg)
}
