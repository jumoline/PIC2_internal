#' Function to visualise GSEA results using dot plot
#'
#' \code{oGSEAdotplot} is supposed to visualise GSEA results using dot plot. It returns an object of class "ggplot" or a list of "ggplot" objects.
#'
#' @param eGSEA an object of class "eGSEA"
#' @param top the number of the top enrichments to be visualised. Alternatively, the gene set names can be queried
#' @param colormap short name for the colormap. It can be one of "jet" (jet colormap), "bwr" (blue-white-red colormap), "gbr" (green-black-red colormap), "wyr" (white-yellow-red colormap), "br" (black-red colormap), "yr" (yellow-red colormap), "wb" (white-black colormap), and "rainbow" (rainbow colormap, that is, red-yellow-green-cyan-blue-magenta). Alternatively, any hyphen-separated HTML color names, e.g. "blue-black-yellow", "royalblue-white-sandybrown", "darkgreen-white-darkviolet". A list of standard color names can be found in \url{https://html-color-codes.info/color-names/index.html}
#' @param zlim the minimum and maximum z values for which colors should be plotted
#' @param ncolors the number of colors specified over the colormap
#' @param xlab the label for x-axis. If NULL, it is 'Target ranks'
#' @param title the title. If NULL, it is term name followed by the number of its annotations
#' @param subtitle the subtitle. It can be used to show 'leading' info, 'enrichment' info or 'both'
#' @param clab the label for colorbar. By default, it is '5-star ratings'
#' @param x.scale how to transform the x scale. It can be "normal" for no transformation, "sqrt" for square root transformation, and "log" for log-based transformation
#' @param peak logical to indicate whether the peak location is shown
#' @param peak.color the peak color
#' @param leading logical to indicate whether the leading targets are texted. Alterntively, leading can be numeric to restict the top targets displayed
#' @param leading.size the size of leading targets' texts. It only works when the parameter 'leading' is enabled
#' @param leading.color the label color of leading targets' texts
#' @param leading.alpha the 0-1 value specifying transparency of leading targets' texts
#' @param leading.padding the padding around the leading targets' texts
#' @param leading.arrow the arrow pointing to the leading targets
#' @param leading.force the repelling force between leading targets' texts
#' @param leading.query which genes in query will be labelled. By default, it sets to NULL meaning all genes will be displayed. If labels in query can not be found, then all will be displayed
#' @param leading.query.only logical to indicate whether only genes in query will be displayed. By default, it sets to FALSE. It only works when labels in query are enabled/found
#' @param leading.edge.only logical to indicate whether only the leading edge will be shown. By default, it sets to FALSE
#' @param leading.label.direction the leading label direction. It can be "none", "left" (aligned to the left-most edge)
#' @param compact logical to indicate whether the compact/void theme is used. If TRUE, axes and legend info will be hidden
#' @param font.family the font family for texts
#' @param signature logical to indicate whether the signature is assigned to the plot caption. By default, it sets TRUE showing which function is used to draw this graph
#' @param ... additional paramters associated with ggrepel::geom_text_repel. If queried, it has high priority (eg, color='darkred',size=2,alpha=0.6,fontface='bold')
#' @return an object of class "ggplot" or a list of "ggplot" objects.
#' @note none
#' @export
#' @seealso \code{\link{oGSEAdotplot}}
#' @include oGSEAdotplot.r
#' @examples
#' \dontrun{
#' gp <- oGSEAdotplot(eGSEA, top=1)
#' #gp <- oGSEAdotplot(eGSEA, top=1, peak=FALSE, compact=TRUE, signature=FALSE)
#' gp
#' 
#' ls_gp <- oGSEAdotplot(eGSEA, top=1:4, signature=FALSE)
#' library(gridExtra)
#' grid.arrange(grobs=ls_gp, ncol=2)
#' }

oGSEAdotplot <- function(eGSEA, top=1, colormap="spectral", zlim=NULL, ncolors=64, xlab=NULL, title=c('setID','none'), subtitle=c('leading','enrichment','both','none'), clab='Priority\nrating', x.scale=c("normal","sqrt","log"), peak=TRUE, peak.color='black', leading=FALSE, leading.size=2.5, leading.color='steelblue4', leading.alpha=1, leading.padding=0.2, leading.arrow=0.01, leading.force=0.01, leading.query=NULL, leading.query.only=FALSE, leading.edge.only=FALSE, leading.label.direction=c("none","left"), compact=FALSE, font.family="sans", ...)
{3
	
	x.scale <- match.arg(x.scale)
	title <- match.arg(title)
	subtitle <- match.arg(subtitle)
	leading.label.direction <- match.arg(leading.label.direction)
	
    if(!is(eGSEA,"eGSEA")){
    	stop("The function must apply to a 'eGSEA' object.\n")
    }
    
    df_summary <- eGSEA$df_summary %>% as.data.frame()
    nSet <- nrow(df_summary)
    
    ## determine which gene set
    if(is(top,"integer") | is(top,"numeric")){
     	top <- as.integer(top)
    	ind <- which((top <= nSet) & (top >= 1))
        if(length(ind)>0){
        	which.terms <- top[ind]
        }else{
        	which.terms <- NULL
        }
        
    }else{
        ind <- which(df_summary$setID %in% top)
        if(length(ind)>0){
        	which.terms <- ind
        }else{
        	which.terms <- NULL
        }
        
    }
    
    if(is.null(which.terms)){
    	return(NULL)
    }
    
    Hits <- Rank <- RES <- Score <- x <- GeneID <- NULL
    ls_gp <- lapply(which.terms, function(which.term){

		df_full <- eGSEA$full[[which.term]]
		df_leading <- subset(df_full, Hits==3)
		
		nLead <- df_summary[which.term, "nLead"]
		nAnno <- df_summary[which.term, "nAnno"]
		nes <- df_summary[which.term, "nes"]
		pvalue <- df_summary[which.term, "pvalue"]
		adjp <- df_summary[which.term, "adjp"]
		
		leading_info <- paste("Peak (rank=", df_leading$Rank, ")",
						 "\nLeading genes (n=", nLead, ")",
						 "\nSignificance (NES=", nes,
						 ", P=", pvalue,
						 ", FDR=", adjp,")",
						 sep="",collapse="")
		
		###############
		if(leading.edge.only){
			df_full <- subset(df_full,Rank<=df_leading$Rank)
		}
		###############
				
		bp <- ggplot(df_full, aes(x=Rank, y=RES, colour=Score))
		bp <- bp + geom_point(size=0.5)
		bp <- bp + geom_hline(yintercept=0, color="grey")
		bp <- bp + geom_segment(data=subset(df_full,Hits>=1), aes(xend=Rank, yend=0), size=0.4)
		if(is.null(zlim)){
			zlim[1] <- min(df_full$Score)
			zlim[2] <- max(df_full$Score)
		}
		bp <- bp + scale_colour_gradientn(colors=oColormap(colormap)(ncolors), limits=zlim, guide=guide_colorbar(title=clab,title.position="top",barwidth=0.5,nbin=5,draw.ulim=FALSE,draw.llim=FALSE))
		
		if(leading | leading>0){
			if(leading>1){
				df_genes <- subset(df_full,Hits>=1 & Rank<=leading)
			}else{
				df_genes <- subset(df_full,Hits>=1)
			}
			
			vec <- eGSEA$leading[[which.term]]
			if(0){
				# why this?
				ind <- match(vec, df_genes$Rank)
				df_genes <- df_genes[ind,]
				df_genes$GeneID <- names(vec)
			}else{
				ind <- match(df_genes$GeneID, names(vec))
				df_genes <- df_genes[!-is.na(ind),]
			}
			
			df_genes_query <- NULL
			df_genes_noquery <- NULL
			if(!is.null(leading.query)){
				ind <- match(df_genes$GeneID, leading.query)
				if(sum(!is.na(ind))!=0){
					df_genes_query <- df_genes[!is.na(ind),]
					df_genes_noquery <- df_genes[is.na(ind),]
				}
			}
			
			if(is.null(df_genes_query)){
				if(leading.label.direction=='none'){
					bp <- bp + ggrepel::geom_text_repel(data=df_genes, aes(x=Rank,y=RES,label=GeneID), lineheight=0.8, size=leading.size, color=leading.color, alpha=leading.alpha, fontface='italic', box.padding=unit(0.5,"lines"), point.padding=unit(leading.padding,"lines"), segment.color='grey80', segment.alpha=0.5, segment.size=0.2, arrow=arrow(length=unit(leading.arrow,'npc')), force=leading.force, max.overlaps=Inf, ...)
				}else{
					bp <- bp + ggrepel::geom_text_repel(data=df_genes, aes(x=Rank,y=RES,label=GeneID), lineheight=0.8, size=leading.size, color=leading.color, alpha=leading.alpha, fontface='italic', box.padding=unit(0.5,"lines"), point.padding=unit(leading.padding,"lines"), segment.color='grey80', segment.alpha=0.5, segment.size=0.2, arrow=arrow(length=unit(leading.arrow,'npc')), max.overlaps=Inf, direction="y", hjust=0, nudge_x=-1*df_genes$Rank-0.05*df_leading$Rank, ...)
				}
				
			}else{
				if(leading.query.only){
					if(leading.label.direction=='none'){
						bp <- bp + ggrepel::geom_text_repel(data=df_genes_query, aes(x=Rank,y=RES,label=GeneID), lineheight=0.8, size=leading.size, color=leading.color, alpha=leading.alpha, fontface='italic', box.padding=unit(0.5,"lines"), point.padding=unit(leading.padding,"lines"), segment.color='grey80', segment.alpha=0.5, segment.size=0.2, arrow=arrow(length=unit(leading.arrow,'npc')), force=leading.force, max.overlaps=Inf, ...)
					}else{
						bp <- bp + ggrepel::geom_text_repel(data=df_genes_query, aes(x=Rank,y=RES,label=GeneID), lineheight=0.8, size=leading.size, color=leading.color, alpha=leading.alpha, fontface='italic', box.padding=unit(0.5,"lines"), point.padding=unit(leading.padding,"lines"), segment.color='grey80', segment.alpha=0.5, segment.size=0.2, arrow=arrow(length=unit(leading.arrow,'npc')), max.overlaps=Inf, direction="y", hjust=0, nudge_x=-1*df_genes_query$Rank-0.05*df_leading$Rank, ...)
					}
					
				}else{
					bp <- bp + ggrepel::geom_text_repel(data=df_genes_noquery, aes(x=Rank,y=RES,label=GeneID), lineheight=0.8, size=leading.size, color=leading.color, alpha=leading.alpha, fontface='italic', box.padding=unit(0.5,"lines"), point.padding=unit(leading.padding,"lines"), segment.color='grey50', segment.alpha=0.5, segment.size=0.2, arrow=arrow(length=unit(leading.arrow,'npc')), force=leading.force, max.overlaps=Inf)
					bp <- bp + ggrepel::geom_text_repel(data=df_genes_query, aes(x=Rank,y=RES,label=GeneID), lineheight=0.8, box.padding=unit(0.5,"lines"), point.padding=unit(leading.padding,"lines"), segment.color='grey50', segment.alpha=0.5, segment.size=0.2, arrow=arrow(length=unit(leading.arrow,'npc')), force=leading.force, max.overlaps=Inf, ...)
					
				}
			}

		}else{
			df_genes <- subset(df_full,Hits>=1)
		}
		
		if(peak){
			bp <- bp + geom_point(data=df_leading, aes(x=Rank, y=RES), colour=peak.color, alpha=1, size=1) + geom_segment(data=df_leading, aes(xend=Rank,yend=0), linewidth=0.5, colour=peak.color, linetype="solid") 
			#bp <- bp + ggrepel::geom_text_repel(data=df_leading, aes(x=Rank,y=RES,label=leading_info), size=2, color='blue', alpha=0.8, fontface='bold.italic')
		}
		
		bp <- bp  + theme_bw() + theme(legend.position="right", legend.title=element_text(color="black",face="bold",size=9), axis.title.y=element_text(color="black"), axis.title.x=element_text(color="black"), panel.border=element_blank(), panel.grid.minor=element_blank(), panel.grid.major=element_blank())
		if(is.null(xlab)){
			if(leading.edge.only){
				if(x.scale=='log'){
					xlab <- "Target ranks (log-scaled) at the leading prioritisation"
				}else{
					xlab <- "Target ranks at the leading prioritisation"
				}
				
			}else{
				xlab <- paste0("Target ranks (from 1 to ", nrow(df_full), ")")	
			}
			
		}
		bp <- bp + xlab(xlab) + ylab("Running enrichment score")

		## title
		if(title=='none'){
			title <- NA
		}else if(title=='setID'){
			title <- paste0(df_summary[which.term,"setID"], " (n=", nAnno, ")")
		}
		
		if(is.null(title)){
			title <- paste0(df_summary[which.term,"setID"], " (n=", nAnno, ")")
		}
		if(subtitle=='both'){
			subtitle <- paste("Peak (rank=", df_leading$Rank, "), ",
							 "Leading targets (n=", nLead, " out of ", nAnno,")\n",
							 "Enrichment (NES=", nes,
							 ", P-value=", pvalue,
							 ", FDR=", adjp,")",
							 sep="",collapse="")
		}else if(subtitle=='leading'){
			subtitle <- paste("Peak (rank=", df_leading$Rank, "), ",
							 "Leading targets (n=", nLead, " out of ", nAnno,")",
							 sep="",collapse="")
		}else if(subtitle=='enrichment'){
			subtitle <- paste("Enrichment (NES=", nes,
							 ", P-value=", pvalue,
							 ", FDR=", adjp,")",
							 sep="",collapse="")
		}else{
			subtitle <- ''
		}
		if(subtitle!=''){
			bp <- bp + labs(title=title, subtitle=subtitle) + theme(plot.title=element_text(hjust=0.5,size=10), plot.subtitle=element_text(hjust=0.5,size=8))
		}else{
			if(!is.na(title)){
				bp <- bp + labs(title=title) + theme(plot.title=element_text(hjust=0.5,size=12))
			}
		}
	
		## x scale
		if(x.scale=="sqrt"){
			x <- NULL
			bp <- bp + scale_x_continuous(trans=scales::sqrt_trans(), breaks=scales::trans_breaks("log10", function(x) 10^x, n=4))
		}else if(x.scale=="log"){
			x <- .x <- NULL
			#bp <- bp + scale_x_continuous(trans=scales::log_trans(), breaks=scales::trans_breaks("log10", function(x) 10^x, n=4)) + annotation_logticks(sides='b')
			bp <- bp + scale_x_log10(breaks=scales::trans_breaks("log10", function(x) 10^x, n=4), labels = scales::trans_format("log10", scales::math_format(10^.x))) + annotation_logticks(sides='b')
			
		}

		## change font family to 'Arial'
		bp <- bp + theme(text=element_text(family=font.family))

		## put arrows on x- and y-axis
		#gp <- bp + theme(axis.line.x=element_line(arrow=arrow(angle=30,length=unit(0.25,"cm"), type="open")), axis.line.y=element_line(arrow=arrow(angle=30,length=unit(0.25,"cm"), type="open")))
		gp <- bp + theme(axis.line.x=element_line(), axis.line.y=element_line())
		
		# whether is compact
		if(compact){
			gp <- gp + theme_void() + theme(legend.position="none")
			if(!is.na(title)){
				gp <- gp + labs(title=title) + theme(plot.title=element_text(hjust=0.5,size=8),plot.margin=unit(rep(0,4),rep("lines",4)))
			}
			
		}
		
		gp$leading <- df_genes
		
		invisible(gp)
    })
    names(ls_gp) <- which.terms
    
    if(length(ls_gp)==1){
    	invisible(ls_gp[[1]])
    }else{
    	invisible(ls_gp)
    }
}
