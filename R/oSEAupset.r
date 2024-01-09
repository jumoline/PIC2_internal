#' Function to visualise enrichment results using a upset plot
#'
#' \code{oSEAupset} is supposed to visualise enrichment results using a upset plot. The top is the kite plot, and visualised below is the combination matrix for overlapped genes. It returns an object of class "ggplot".
#'
#' @param obj an object of class "eSET". Alterntively, it can be a tibble having all these columns (named as 'name','adjp','nO','overlap')
#' @param top the number of the top terms (sorted according to adjp). For the "eSET" object, if it is 'auto' (for eSET), only the significant terms (see below adjp.cutoff) will be displayed
#' @param adjp.cutoff FDR cutoff used to declare the significant terms. By default, it is set to 0.05. Only works when top is 'auto' above
#' @param color the point color
#' @param shape the point shape
#' @param size.range the point size range
#' @param slim the minimum and maximum point size values, defaulting to the range of the nO
#' @param size.title the size legend title
#' @param wrap.width a positive integer specifying wrap width of name
#' @param label.height.unit the unit specifying the per-row height of the combination matrix. By default, it is NULL. If specified (such as 8), it will be used to decide 'combmatrix.label.height' and 'combmatrix.label.text'
#' @param member.levels how to define the levels of members in the combination matri. It can be 'num' (the number of non-zeros for a member) or 'customised' (see 'levels.customised' below) or 'auto' (alphabeta order; by default)
#' @param levels.customised the customised levels of members in the combination matrix. Only works when member.levels is 'customised' above
#' @param sortBy how to sort bars/combinations. It can be 'pvalue' or 'name'
#' @return an object of class "ggplot"
#' @note none
#' @export
#' @seealso \code{\link{oSEAupset}}
#' @include oSEAupset.r
#' @examples
#' \dontrun{
#' obj %>% oSEAupset()
#' obj %>% oSEAextract() %>% filter(onto=='PSG') %>% oSEAupset()
#' }

oSEAupset <- function(obj, top=10, adjp.cutoff=0.05, color="steelblue4", shape=18, size.range=c(2,5), slim=NULL, size.title="Overlap", wrap.width=NULL, label.height.unit=5.5, member.levels=c("auto","num","customised"), levels.customised=NULL, sortBy=c("pvalue","name"))
{
    
    if(is.null(obj)){
        warnings("There is no enrichment in the 'obj' object.\n")
        return(NULL)
    }
    
    member.levels <- match.arg(member.levels)
    sortBy <- match.arg(sortBy)
    
	name <- nO <- adjp <- overlap <- n <- NULL
    name_overlap <- NULL
    
    adjp.cutoff <- as.numeric(adjp.cutoff)
    if(adjp.cutoff==1){
    	adjp.cutoff <- adjp.cutoff + 0.01
    }
    
    df <- NULL
    if(is(obj,'eSET')){
    	df <- obj %>% oSEAextract()
	}else{
		df <- obj
	}
	
	if(is(df,'data.frame')){
		
		if(all(c('name','adjp','overlap') %in% colnames(df))){
		
			if(!('nO' %in% colnames(df))){
				df %>% dplyr::mutate(nO=1) -> df
			}
			
			df <- df %>% mutate(name_overlap=str_c(name,',',overlap))
		}
	}
	
	if(is(top,'numeric')){
		top <- as.integer(top)
		df %>% dplyr::top_n(top,-adjp) -> df
	}
	
	##########################
	## important: the column 'overlap' is not necessarily unique
	
	if(sortBy=='pvalue'){
		pvalue <- name <- NULL
		df %>% dplyr::filter(adjp<adjp.cutoff) %>% dplyr::arrange(pvalue,name) %>% distinct(name_overlap,.keep_all=TRUE) %>% dplyr::mutate(name_overlap=forcats::fct_inorder(name_overlap)) -> df
	}else if(sortBy=='name'){
		name_overlap <- name <- NULL
		df %>% dplyr::filter(adjp<adjp.cutoff) %>% dplyr::arrange(name) %>% dplyr::distinct(name_overlap,.keep_all=TRUE) %>%  dplyr::mutate(name_overlap=forcats::fct_inorder(name_overlap)) -> df
	}
	#df %>% dplyr::filter(adjp<adjp.cutoff) %>% dplyr::arrange(adjp,name) -> df
	
	##########################
	##########################
	if(nrow(df)==0){
		return(NULL)
	}
	##########################
	##########################

	## text wrap
	if(!is.null(wrap.width)){
		width <- as.integer(wrap.width)
		res_list <- lapply(df$name, function(x){
			x <- gsub('_', ' ', x)
			y <- strwrap(x, width=width)
			if(length(y)>1){
				paste0(y[1], '...')
			}else{
				y
			}
		})
		df$name <- unlist(res_list)
	}
	
	size <- NULL
	
	df$size <- df$nO
	if(is.null(slim)){
		tmp <- df$size
		slim <- c(0, ceiling(max(tmp[!is.infinite(tmp)])))
	}
	df$size[df$size<=slim[1]] <- slim[1]
	df$size[df$size>=slim[2]] <- slim[2]
	
	###########################################
	
	df %>% ggplot(aes(x=name_overlap,y=-log10(adjp))) + geom_col(fill="steelblue", color='transparent', width=0.1, alpha=0.8) + geom_point(aes(size=nO), shape=shape, color=color, fill='white', alpha=0.8) + theme_classic() + theme(axis.title.x=element_blank()) + ggrepel::geom_text_repel(aes(label=name), color='black', size=2, segment.color='grey80', segment.alpha=0.5, max.overlaps=Inf) -> gp
	
	########################
	levels <- NULL
	if(member.levels=='num'){
		#df %>% tidyr::separate_rows(overlap, sep=', |,') %>% dplyr::count(overlap) %>% dplyr::arrange(-n,overlap) %>% dplyr::pull(overlap) -> levels
		
		tmp_overlap <- NULL
		df %>% dplyr::mutate(tmp_overlap=as.character(overlap)) %>% tidyr::separate_rows(tmp_overlap, sep=', |,') %>% dplyr::count(tmp_overlap) %>% dplyr::arrange(-n,tmp_overlap) %>% dplyr::pull(tmp_overlap) -> levels
		
	}else if(member.levels=='auto'){
		#df %>% tidyr::separate_rows(overlap, sep=', |,') %>% dplyr::count(overlap) %>% dplyr::arrange(overlap) %>% dplyr::pull(overlap) -> levels
	
		tmp_overlap <- NULL
		df %>% dplyr::mutate(tmp_overlap=as.character(overlap)) %>% tidyr::separate_rows(tmp_overlap, sep=', |,') %>% dplyr::count(tmp_overlap) %>% dplyr::arrange(tmp_overlap) %>% dplyr::pull(tmp_overlap) -> levels
		
	}else if(member.levels=='customised'){
		levels <- levels.customised
	}
	########################
	gp + ggupset::axis_combmatrix(sep=", |,", levels=levels) + ggupset::theme_combmatrix(combmatrix.panel.point.color.fill="steelblue", combmatrix.panel.point.color.empty="grey90", combmatrix.panel.point.size=1.5, combmatrix.panel.line.size=0.1, combmatrix.label.height=unit(label.height.unit*(length(levels)+1),"pt"), combmatrix.label.text=element_text(size=label.height.unit)) -> gp
	#gp + facet_grid(.~namespace,scales="free_x",space="free_x")
	gp <- gp + theme(axis.title.x=element_blank()) + ylab(expression(-log[10]("FDR")))
	
	## size
	gp <- gp + scale_size_continuous(limits=slim, range=size.range, guide=guide_legend(size.title,title.position="top",keywidth=0.4,ncol=1,order=1))
	if(length(size.range)==1 | df %>% dplyr::count(nO) %>% nrow() == 1){
		gp <- gp + guides(size=FALSE)
	}
	
	gp
}

