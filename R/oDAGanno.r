#' Function to generate a DAG subgraph induced by the input annotation data
#'
#' \code{oDAGanno} is supposed to produce a subgraph induced by the input annotation data, given a direct acyclic graph (DAG; an ontology). The input is a graph of "igraph", a list of the vertices containing annotation data, and the mode defining the paths to the root of DAG. The induced subgraph contains vertices (with annotation data) and their ancestors along with the defined paths to the root of DAG. The annotations at these vertices (including their ancestors) can also be updated according to the true-path rule: those annotated to a term should also be annotated by its all ancestor terms.
#'
#' @param ig an object of class "igraph" to represent DAG
#' @param anno the vertices/nodes for which annotation data are provided. It can be a sparse Matrix of class "dgCMatrix" (with variants/genes as rows and terms as columns), or a list of nodes/terms each containing annotation data, or an object of class 'SET'
#' @param path.mode the mode of paths induced by vertices/nodes with input annotation data. It can be "all_paths" for all possible paths to the root, "shortest_paths" for only one path to the root (for each node in query), "all_shortest_paths" for all shortest paths to the root (i.e. for each node, find all shortest paths with the equal lengths), "shortest_paths_epath" (indeed only edges in the shortest path are kept)
#' @param true.path.rule logical to indicate whether the true-path rule should be applied to propagate annotations. By default, it sets to true
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @return 
#' NULL or an induced subgraph, an object of class "igraph". In addition to the original attributes to nodes and edges, the return subgraph is also appended by two node attributes: 1) "anno" containing a list of variants/genes either as original annotations (and inherited annotations; 2) "IC" standing for information content defined as negative 10-based log-transformed frequency of variants/genes annotated to that term.
#' @note For the mode "shortest_paths", the induced subgraph is the most concise, and thus informative for visualisation when there are many nodes in query, while the mode "all_paths" results in the complete subgraph.
#' @export
#' @seealso \code{\link{oDAGinduce}}
#' @include oDAGanno.r
#' @examples
#' \dontrun{
#' ###########################################################
#' g <- oRDS('ig.GOMF', placeholder=placeholder)
#' anno <- oRDS('org.Hs.egGOMF', placeholder=placeholder)
#' dag <- oDAGanno(g, anno)
#' }

oDAGanno <- function(ig, anno, path.mode=c("all_paths","shortest_paths","all_shortest_paths","shortest_paths_epath"), true.path.rule=TRUE, verbose=TRUE)
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
    
    if(is(anno,"SET")){
    	name <- member <- NULL
        originAnnos <- anno$info %>% dplyr::select(name,member) %>% tibble::deframe()
    }else if(is(anno,"list")){
        originAnnos <- anno
    }else if(is(anno,"dgCMatrix")){
		D <- anno
		originAnnos <- sapply(1:ncol(D), function(j){
			names(which(D[,j]!=0))
		})
		names(originAnnos) <- colnames(anno)
    }else{
    	warnings("The input annotation must be either 'SET' or 'list' or 'dgCMatrix' object.\n")
    	return(NULL)
    }
    
    ## check nodes in annotation
    if (is.list(originAnnos)){
        originNodes <- names(originAnnos)
        
        ind <- match(originNodes, V(ig)$name)
        nodes_mapped <- originNodes[!is.na(ind)]
        if(length(nodes_mapped)==0){
            warnings("The input annotation data do not contain terms matched to the nodes/terms in the input graph.\n")
            return(NULL)
        }
    }
    
	if(verbose){
		message(sprintf("Generate a DAG subgraph via `%s` (%s) ...", path.mode, as.character(Sys.time())), appendLF=TRUE)
	}
    
    ## generate a subgraph of a direct acyclic graph (DAG) induced by terms from input annotations
    dag <- oDAGinduce(ig, originNodes, path.mode=path.mode)
    allNodes <- V(dag)$name
    
	## create a new (empty) hash environment
	## node2domain.HoH: 1st key (node/term), 2nd key (domain), value (origin/inherit)
	node2domain.HoH <- new.env(hash=TRUE, parent=emptyenv())
	
	## assigin original annotations to "node2domain.HoH"
	lapply(allNodes, function(node){
		e <- new.env(hash=TRUE, parent=emptyenv())
		if(node %in% originNodes){
			sapply(originAnnos[[node]], function(x){
				assign(as.character(x), "o", envir=e)
			})
		}
		assign(node, e, envir=node2domain.HoH)
	})

    ## whether true-path rule will be applied
    if(true.path.rule){
    
		if(verbose){
			message(sprintf("Apply true-path rule (%s) ...", as.character(Sys.time())), appendLF=TRUE)
		}
    
		## get the levels list
		level2node <- oDAGlevel(dag, level.mode="longest_path", return.mode="level2node")
		## build a hash environment from the named list "level2node"
		## level2node.Hash: key (level), value (a list of nodes/terms)
		level2node.Hash <- list2env(level2node)
		nLevels <- length(level2node)
		for(i in nLevels:1) {
			currNodes <- get(as.character(i), envir=level2node.Hash, mode='character')
			
			#########################################
			currNodes <- currNodes[!is.na(currNodes)]
			#########################################
			
			## get the incoming neighbors (excluding self) that are reachable (i.e. nodes from i-1 level)
			adjNodesList <- lapply(currNodes, function(node){
				neighs.in <- igraph::neighborhood(dag, order=1, nodes=node, mode="in")
				setdiff(V(dag)[unlist(neighs.in)]$name, node)
			})
			names(adjNodesList) <- currNodes

			## push the domains from level i to level i - 1
			lapply(currNodes, function(node){
				## get the domains from this node
				domainsID <- ls(get(node, envir=node2domain.HoH, mode='environment'))

				## assigin inherit annotations to "node2domain.HoH"
				lapply(adjNodesList[[node]], function(adjNode){
					adjEnv <- get(adjNode, envir=node2domain.HoH, mode='environment')
					sapply(domainsID, function(domainID){
						assign(domainID, "i", envir=adjEnv)
					})
				})
			})
		
			if(verbose){
				message(sprintf("\tAt level %d, there are %d nodes, and %d incoming neighbors.", i, length(currNodes), length(unique(unlist(adjNodesList)))), appendLF=TRUE)
			}
		
		}
	}
	
	node2domains <- as.list(node2domain.HoH)[allNodes]
	domain_annotation <- lapply(node2domains, function(node){
		vec <- unlist(as.list(node))
		res <- names(vec)
		names(res) <- vec
		sort(res)
	})
	
    ## append 'anno' attributes to the graph
    V(dag)$anno <- domain_annotation

    ## append 'IC' attributes to the graph
    counts <- sapply(V(dag)$anno, length)
    IC <- -1*log10(counts/max(counts))
    ### force those 'Inf' to be NA
    if(1){
    	IC[is.infinite(IC)] <- NA
    }
    V(dag)$IC <- IC
   ####################################################################################
    endT <- Sys.time()
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    
    if(verbose){
    	message(sprintf("\nEnded (at %s)", as.character(endT)), appendLF=TRUE)
    	message(sprintf("Runtime in total: %d secs\n", runTime), appendLF=TRUE)
    }
    
    return(dag)
}
