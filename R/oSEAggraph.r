#' Function to visualise enrichment results using a ggraph-like lauout
#'
#' \code{oSEAggraph} is supposed to visualise enrichment results using a ggraph-like lauout.
#'
#' @param obj an object of class "eset" or "ls_eset". Alterntively, it can be a data frame having all these columns ('name','adjp','or','zscore'; 'group' optionally)
#' @param ig an object of class "igraph" with node attribute 'name'. Note: the node labels would be the node attribute 'name' unless the node attribute 'label' is explicitely provided. If provided, only those terms within it will be visualised. By default, it is NULL meaning no surch restriction
#' @param fixed logical to indicate whether all terms in ig will be visualised. By default, it is TURE; otherwise only overlapped terms from eset will be visualised
#' @param node.color which statistics will be used for node coloring. It can be "or" for the odds ratio, "adjp" for adjusted p value (FDR) and "zscore" for enrichment z-score
#' @param colormap short name for the colormap. It can be one of "jet" (jet colormap), "bwr" (blue-white-red colormap), "gbr" (green-black-red colormap), "wyr" (white-yellow-red colormap), "br" (black-red colormap), "yr" (yellow-red colormap), "wb" (white-black colormap), "rainbow" (rainbow colormap, that is, red-yellow-green-cyan-blue-magenta), and "ggplot2" (emulating ggplot2 default color palette). Alternatively, any hyphen-separated HTML color names, e.g. "lightyellow-orange" (by default), "blue-black-yellow", "royalblue-white-sandybrown", "darkgreen-white-darkviolet". A list of standard color names can be found in \url{https://html-color-codes.info/color-names/index.html}
#' @param zlim the minimum and maximum values for which colors should be plotted
#' @param node.size which statistics will be used for node size. It can be "or" for the odds ratio, "adjp" for adjusted p value (FDR) and "zscore" for enrichment z-score
#' @param slim the minimum and maximum values for which sizes should be plotted
#' @param node.size.range the range of actual node size
#' @param node.label.size the text size of the node labelings. By default, it is 2. If 0, all labellings will be disabled
#' @param leave the logic specifying whether or not only leaves (nodes/labellings) shown. This can be disenabled if the layout does not support tips
#' @param ncolumns an integer specifying the number of columns for facet_wrap. By defaul, it is NULL (decided on according to the number of groups that will be visualised)
#' @param ... additional graphic parameters used in oGGraph
#' @return
#' a ggplot2 object appended with 'ig', 'data' which should contain columns 'x','y', 'name' (the same as V(ig)$name), 'label' (if not given in ig, a 'name' varient), 'data_enrichment' (enrichment results), and 'gp_template' with labelling if multiple groups (together with no labelling for the colored plots).
#' @note none
#' @export
#' @seealso \code{\link{oGGraph}}
#' @include oSEAggraph.r
#' @examples
#' \dontrun{
#' # circular visualisation of enriched AA pathways
#' gp <- oSEAggraph(eset, ig)
#' 
#' ###############################
#' # advanced use: multiple groups
#' 
#' # circular visualisation of enriched AA pathways
#' gp <- oSEAggraph(ls_eset, ig)
#' gp
#' gp$gp_template
#' }

oSEAggraph <- function(obj, ig=NULL, fixed=FALSE, node.color=c("zscore","adjp","or"), colormap="grey90-orange-darkred", zlim=NULL, node.size=c("adjp","zscore","or"), slim=NULL, node.size.range=c(0.5,4), node.label.size=2, leave=TRUE, ncolumns=NULL, ...)
{
	## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    node.color <- match.arg(node.color)
    node.size <- match.arg(node.size)
    
    if(is.null(obj)){
        warnings("There is no enrichment in the 'obj' object.\n")
        return(NULL)
    }
    
    if(is(obj,'eSET')){
    	df_enrichment_group <- obj$info
		df_enrichment_group$group <- 'group'
	}else if(is(obj,'eSAD')){	
		obj %>% oSEAextract() -> df_enrichment_group
	}else{
		df_enrichment_group <- obj
	}
	
	if(is(df_enrichment_group,'data.frame')){
		
		if(all(c('name','adjp','or','zscore') %in% colnames(df_enrichment_group))){
			if(!('group' %in% colnames(df_enrichment_group))){
				df_enrichment_group %>% dplyr::mutate(group='group') -> df_enrichment_group
			}
		}
	}
	
	if(is(df_enrichment_group$group,'factor')){
		if(length(unique(df_enrichment_group$group)) != length(levels(df_enrichment_group$group))){
			df_enrichment_group$group <- factor(df_enrichment_group$group, levels=sort(unique(df_enrichment_group$group)))
		}
	}

	##########################################################
	# restrict those nodes provided in 'ig'
	if(!is(ig,"igraph")){
		if(is(obj,'eSET') & is(obj$ig,'igraph')){
			ig <- obj$ig
		}else{
			return(NULL)
		}
	}
	
	if(!fixed){
		ind <- match(V(ig)$name, df_enrichment_group$name)
		nodes_query <- V(ig)$name[!is.na(ind)]
		if(is(suppressWarnings(try(ig <- oDAGinduce(ig, nodes_query, path.mode="all_paths"), TRUE)),"try-error")){
			ig <- NULL
		}
	}
	##########################################################

	gp <- NULL
	
	if(is(ig,"igraph")){
		
		#########################
		## replace those infinite
		df_enrichment_group$or[is.infinite(df_enrichment_group$or)] <- max(df_enrichment_group$or[!is.infinite(df_enrichment_group$or)])
		#########################
		
		# convert 'ig' into 'ls_igg' by group
		ls_df <- split(x=df_enrichment_group[,c("name","zscore","adjp","or")], f=df_enrichment_group$group)
		ls_igg <- lapply(ls_df, function(df_enrichment){
			igg <- ig
			V(igg)$zscore <- 0
			ind <- match(V(igg)$name, df_enrichment$name)
			V(igg)$zscore[!is.na(ind)] <- df_enrichment$zscore[ind[!is.na(ind)]]

			V(igg)$adjp <- 0
			ind <- match(V(igg)$name, df_enrichment$name)
			V(igg)$adjp[!is.na(ind)] <- -1*log10(df_enrichment$adjp)[ind[!is.na(ind)]]
	
			V(igg)$or <- 0
			ind <- match(V(igg)$name, df_enrichment$name)
			V(igg)$or[!is.na(ind)] <- log2(df_enrichment$or)[ind[!is.na(ind)]]
			igg
		})
		
		if(node.color=="or"){
			node.color.title <- expression(log[2](OR))
			if(is.null(zlim)){
				zlim <- c(0, ceiling(max(log2(df_enrichment_group$or))))
			}
		}else if(node.color=="adjp"){
			node.color.title <- expression(-log[10](adjP))
			if(is.null(zlim)){
				zlim <- c(0, ceiling(max(-1*log10(df_enrichment_group$adjp))))
			}
		}else if(node.color=="zscore"){
			node.color.title <- "Z-score"
			if(is.null(zlim)){
				zlim <- c(0, ceiling(max(df_enrichment_group$zscore)))
			}
		}
		
		if(node.size=="or"){
			node.size.title <- expression(log[2](OR))
			if(is.null(slim)){
				slim <- c(0, ceiling(max(log2(df_enrichment_group$or))))
			}
		}else if(node.size=="adjp"){
			node.size.title <- expression(-log[10](adjP))
			if(is.null(slim)){
				slim <- c(0, ceiling(max(-1*log10(df_enrichment_group$adjp))))
			}
		}else if(node.size=="zscore"){
			node.size.title <- "Z-score"
			if(is.null(slim)){
				slim <- c(0, ceiling(max(df_enrichment_group$zscore)))
			}
		}
		
		x <- y <- leaf <- x_group <- y_group <- group <- group_id <- color <- size <- NULL
		
		if(length(ls_igg)>1){
			gp_template <- oGGraph(ig, leave=leave, node.label.size=node.label.size, ...)
			gp <- oGGraph(ig, leave=leave, node.label.size=0, node.size=0, ...)
		}else{
			gp <- oGGraph(ig, leave=leave, node.label.size=node.label.size, node.size=0, ...)
		}
		
		## node/leave labels
		if(leave & !is.null(gp$data$leaf)){
			df_data <- subset(gp$data, leaf==T)
		}else{
			df_data <- gp$data
		}
		
		ls_df <- lapply(1:length(ls_igg), function(i){
			g <- ls_igg[[i]]

			group <- names(ls_igg)[i]	

			ind <- match(V(g)$name, df_data$name)
			if(node.color=="or"){
				color <- V(g)$or[!is.na(ind)]
			}else if(node.color=="adjp"){
				color <- V(g)$adjp[!is.na(ind)]
			}else if(node.color=="zscore"){
				color <- V(g)$zscore[!is.na(ind)]
			}
			color[color<=zlim[1]] <- zlim[1]
			color[color>=zlim[2]] <- zlim[2]
		
			if(node.size=="or"){
				size <- V(g)$or[!is.na(ind)]
			}else if(node.size=="adjp"){
				size <- V(g)$adjp[!is.na(ind)]
			}else if(node.size=="zscore"){
				size <- V(g)$zscore[!is.na(ind)]
			}
			size[size<=slim[1]] <- slim[1]
			size[size>=slim[2]] <- slim[2]
		
			df <- data.frame(df_data[ind[!is.na(ind)],], color=color, size=size, group=group, group_id=i, stringsAsFactors=F)
		})
		df <- do.call(rbind, ls_df)
		
		gp <- gp + geom_point(data=df, aes(x=x*1, y=y*1, color=color, size=size))
		
		gp <- gp + scale_size_continuous(limits=slim, range=node.size.range, guide=guide_legend(node.size.title,title.position="top",ncol=1)) + scale_colour_gradientn(colors=oColormap(colormap)(64), limits=zlim, guide=guide_colorbar(title=node.color.title,title.position="top",barwidth=0.4))
		
		## facet by 'group' artificially added 
		if(length(ls_igg)>1){
			if(is.null(ncolumns)){
				ncolumns <- ceiling(sqrt(length(ls_igg)))
			}
			group <- NULL
			gp <- gp + facet_wrap(~group, ncol=ncolumns)
			gp <- gp + theme(strip.background=element_rect(fill="transparent",color="transparent"), strip.text=element_text(size=12,face="bold"), strip.placement="inside", panel.spacing=unit(0,"lines"))
			
			gp$gp_template <- gp_template
		}
		
		gp$data_enrichment <- df
	}
	
    return(gp)
}
