#' Function to create a bipartitle graph
#'
#' \code{oBicreate} is supposed to create a bipartitle graph.
#'
#' @param data a data frame/matrix to create a bipartitle graph
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @return It returns an igraph object.
#' @note none
#' @export
#' @seealso \code{\link{oBicreate}}
#' @include oBicreate.r
#' @examples
#' data <- cbind(c(1,1,0,1), c(0,1,2,0))
#' ig <- oBicreate(data)

oBicreate <- function(data, verbose=TRUE)
{
    
    if(is(data,"matrix") | is(data,"data.frame") | is(data,"dgCMatrix")){
    	if(is(data,"dgCMatrix")){
    		data <- as.matrix(data)
    	}

        if(nrow(data) < ncol(data)){
        	data <- t(data)
        }
        # all zeros replaced wtih NA (no edge)
        data[data==0] <- NA
        
        if(is.null(colnames(data))){
        	colnames(data) <- paste0('C',1:ncol(data))
        }
        if(is.null(rownames(data))){
        	rownames(data) <- paste0('R',1:nrow(data))
        }
    }else{
    	return(NULL)
    }
    
    # df_xnodes
    df_xnodes <- data.frame(name=colnames(data), type='xnode', stringsAsFactors=F)
    # df_ynodes
    df_ynodes <- data.frame(name=rownames(data), type='ynode', stringsAsFactors=F)
	## df_nodes
	df_nodes <- rbind(df_xnodes, df_ynodes)

	## df_relations
	df_relations <- oSM2DF(data, verbose=verbose)
	colnames(df_relations) <- c("from","to","weight")
	
	## create a bipartitle graph
	ig <- igraph::graph.data.frame(d=df_relations, directed=T, vertices=df_nodes)
	
	vec_type <- sort(table(V(ig)$type))
	if(verbose){
		message(sprintf("The igraph object has %d nodes (%d xnode '%s' type and %d ynode '%s' type) and %d edges (%s) ...", vcount(ig), vec_type[1],names(vec_type)[1], vec_type[2],names(vec_type)[2], ecount(ig), as.character(Sys.time())), appendLF=T)
	}

    invisible(ig)
}


