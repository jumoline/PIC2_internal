#' Function to prioritise pathways based on GSEA analysis of prioritised genes
#'
#' \code{oPierGSEA} is supposed to prioritise pathways given prioritised genes and the ontology in query. It is done via gene set enrichment analysis (GSEA). It returns an object of class "eGSEA". 
#
#' @param pNode an object of class "pNode" (or "sTarget" or "dTarget"). Alternatively, it can be a data frame with two columns ('priority' and 'rank')
#' @param priority.top the number of the top targets used for GSEA. By default, it is NULL meaning all targets are used
#' @param customised.genesets a list each containing gene symbols. By default, it is NULL. If the list provided, it will overtake the previous parameter "ontology"
#' @param size.range the minimum and maximum size of members of each term in consideration. By default, it sets to a minimum of 10 but no more than 500
#' @param type It can be 'simple' or 'multilevel'
#' @param weight an integer specifying score weight. It can be "0" for unweighted (an equivalent to Kolmogorov-Smirnov, only considering the rank), "1" for weighted by input gene score (by default), and "2" for over-weighted, and so on
#' @param nperm the number of random permutations. For each permutation, gene-score associations will be permutated so that permutation of gene-term associations is realised
#' @param p.adjust.method the method used to adjust p-values. It can be one of "BH", "BY", "bonferroni", "holm", "hochberg" and "hommel". The first two methods "BH" (widely used) and "BY" control the false discovery rate (FDR: the expected proportion of false discoveries amongst the rejected hypotheses); the last four methods "bonferroni", "holm", "hochberg" and "hommel" are designed to give strong control of the family-wise error rate (FWER). Notes: FDR is a less stringent condition than FWER
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true
#' @param silent logical to indicate whether the messages will be silent completely. By default, it sets to false. If true, verbose will be forced to be false
#' @return
#' an object of class "eGSEA", a list with following components:
#' \itemize{
#'  \item{\code{df_summary}: a data frame of nTerm x 12 containing gene set enrichment analysis result, where nTerm is the number of terms/genesets, and the 9 columns are "setID" (i.e. "Term ID"), "nAnno" (i.e. number in members annotated by a term), "nLead" (i.e. number in members as leading genes), "peak" (i.e. the rank at peak), "total" (i.e. the total number of genes analysed), "es" (i.e. enrichment score), "nes" (i.e. normalised enrichment score; enrichment score but after being normalised by gene set size), "pvalue" (i.e. nominal p value), "adjp" (i.e. adjusted p value; p value but after being adjusted for multiple comparisons), "frac" (nLead/nAnno), "member" (members at the leading edge), and "member_rank" (members and their ranks)}
#'  \item{\code{leading}: a list of gene sets, each storing leading gene info (i.e. the named vector with names for gene symbols and elements for priority rank). Always, gene sets are identified by "setID"}
#'  \item{\code{full}: a list of gene sets, each storing full info on gene set enrichment analysis result (i.e. a data frame of nGene x 6, where nGene is the number of genes, and the 6 columns are "GeneID", "Rank" for priority rank, "Score" for priority score, "RES" for running enrichment score,  "Hits" for gene set hits info with 1 for gene hit, 2 for leading gene hit, 3 for the point defining leading genes, 0 for no hit). Always, gene sets are identified by "setID"}
#'  \item{\code{cross}: a matrix of nTerm X nTerm, with an on-diagnal cell for the leading genes observed in an individaul term, and off-diagnal cell for the overlapped leading genes shared between two terms}
#' }
#' @note none
#' @export
#' @seealso \code{\link{oPierGSEA}}
#' @include oPierGSEA.r
#' @examples
#' \dontrun{
#' # a) provide the seed nodes/genes with the weight info
#' ## get genes within 500kb away from AS GWAS lead SNPs
#' seeds.genes <- ImmunoBase$AS$genes_variants
#' ## seeds weighted according to distance away from lead SNPs
#' data <- 1- seeds.genes/500000
#'
#' # b) perform priority analysis
#' pNode <- oPierGenes(data=data, network="PCommonsDN_medium",restart=0.7)
#' 
#' # c) do pathway-level priority using GSEA
#' eGSEA <- oPierGSEA(pNode=pNode, ontology="DGIdb", nperm=2000, placeholder=placeholder)
#' bp <- oGSEAbarplot(eGSEA, top_num="auto", displayBy="nes")
#' gp <- oGSEAdotplot(eGSEA, top=1)
#' }

oPierGSEA <- function(pNode, priority.top=NULL, customised.genesets, size.range=c(10,500), type=c("simple","multilevel"), weight=1, nperm=NULL, p.adjust.method=c("BH","BY","bonferroni","holm","hochberg","hommel"), verbose=TRUE, silent=FALSE)
{
    startT <- Sys.time()
    if(!silent){
    	message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=TRUE)
    	message("", appendLF=TRUE)
    }else{
    	verbose <- FALSE
    }
    ####################################################################################
    
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    type <- match.arg(type)
    p.adjust.method <- match.arg(p.adjust.method)

    weight <- as.integer(weight)
    
    if(is(pNode,"pNode")){
    	if(is(pNode$priority,"tbl")){
    		name1 <- NULL
    		df_priority <- pNode$priority %>% dplyr::mutate(name1=name) %>% tibble::column_to_rownames('name1')
    	}else{
    		df_priority <- pNode$priority
    	}
        df_priority <- df_priority[, c("name","seed","weight","priority","rank")]      
        
    }else if(is(pNode,"sTarget") | is(pNode,"dTarget")){
    	if(is(pNode$priority,"tbl")){
    		name1 <- NULL
    		df_priority <- pNode$priority %>% dplyr::mutate(name1=name) %>% tibble::column_to_rownames('name1')
    	}else{
    		df_priority <- pNode$priority
    	}
    	df_priority <- df_priority[, c("name","rank","rating")]
    	df_priority$priority <- df_priority$rating

	}else if(is(pNode,"data.frame")){
    	df_priority <- pNode[,c(1:2)] %>% as.data.frame()
    	colnames(df_priority) <- c("priority","rank")
    }else{
    	stop("The function must apply to a 'pNode' or 'sTarget' or 'dTarget' object.\n")
    }
    
    ###############
	## priority top
	if(!is.null(priority.top)){
		priority.top <- as.integer(priority.top)
		if(priority.top > nrow(df_priority)){
			priority.top <- nrow(df_priority)
		}else if(priority.top <= 1){
			priority.top <- nrow(df_priority)
		}
	}else{
		priority.top <- nrow(df_priority)
	}
    df_priority <- df_priority[1:priority.top,]
    ###############
    
    data <- df_priority
    
  #############################################################################################

    if(is(customised.genesets,"SET")){
    	name <- member <- NULL
        anno <- customised.genesets$info %>% dplyr::mutate(member=purrr::map(member,~sort(.x))) %>% dplyr::select(name,member) %>% tibble::deframe()
    }else{
        if(is.list(customised.genesets)){
            if(is.null(names(customised.genesets))){
                names(customised.genesets) <- paste("C", 1:length(customised.genesets), sep="")
            }
            #######
            customised.genesets <- customised.genesets[lapply(customised.genesets,length)>0]
            if(length(customised.genesets)==0){
            	return(NULL)
            }	
            anno <- customised.genesets	
        }else{
			stop("There is no input for the ontology.\n")
		}
	}
	
	#####
	ind <- which(sapply(anno,length) >= size.range[1] & sapply(anno,length) <= size.range[2])
	if(length(ind)>0){
		anno <- anno[ind]
	}else{
		return(NULL)
	}
	#####
	
    #############################################################################################
    
    flag_fgsea <- FALSE
    pkgs <- c("fgsea")
    if(all(pkgs %in% rownames(utils::installed.packages()))){
        tmp <- sapply(pkgs, function(pkg) {
            requireNamespace(pkg, quietly=TRUE)
        })
        if(all(tmp)){
        	flag_fgsea <- TRUE
        }
    }
    
    if(flag_fgsea){
    
		if(verbose){
			now <- Sys.time()
			message(sprintf("\n#######################################################"))
			message(sprintf("'fgsea' from the fgsea package is being called (%s):", as.character(now)), appendLF=TRUE)
			message(sprintf("#######################################################"))
		}
		
		stats <- data$priority
		names(stats) <- rownames(data)
		if(type=='simple'){
			## fgseaSimple
			suppressWarnings(fgseaRes <- fgsea::fgseaSimple(pathway=anno, stats=stats, minSize=size.range[1], maxSize=size.range[2], gseaParam=weight, scoreType="pos", nperm=nperm))
		}else if(type=='multilevel'){
			## fgseaMultilevel: based on the adaptive multilevel splitting Monte Carlo approach. This allows us to exceed the results of simple sampling and calculate arbitrarily small P-values
			suppressWarnings(fgseaRes <- fgsea::fgseaMultilevel(pathway=anno, stats=stats, minSize=size.range[1], maxSize=size.range[2], gseaParam=weight, scoreType="pos", nPermSimple=nperm))
		}
		
		tab <- tibble::tibble(setID         = fgseaRes$pathway,
						   ES           = fgseaRes$ES,
						   nES          = fgseaRes$NES,
						   pvalue       = fgseaRes$pval,
						   adjp         = fgseaRes$padj,
						   setSize      = fgseaRes$size
						  )
		pvalue <- nES <- ES <- NULL
		res <- tab %>% dplyr::arrange(pvalue,-nES,-ES)
	
		if(verbose){
			message(sprintf("#######################################################"))
			message(sprintf("'fgsea' has been finished (%s)!", as.character(Sys.time())), appendLF=TRUE)
			message(sprintf("#######################################################\n"))
		}
	
	}else{
	
		return(NULL)
	}

    #############################################################################################
	### get Leading genes
    if(TRUE){

		########################
		geneid <- rownames(data)
		nGene <- nrow(data)
		## score rank
		rank.score <- data$priority
		ind <- order(rank.score, decreasing=TRUE)
		rank.score.sorted <- rank.score[ind]
		geneid.sorted <- geneid[ind]
    
		ls_df_leading <- lapply(res$setID, function(x){
			## initialisation
			nHit <- length(anno[[x]])
			nMiss <- nGene - nHit
			## observed
			observed.point <- rep(-1/nMiss, nGene)
			flag <- match(anno[[x]], geneid.sorted)
			###### remove NA
			flag <- flag[!is.na(flag)]
			######
			if(weight==0) {
				observed.point[flag] <- 1/nHit
			}else if(weight==1){
				hit_tmp <- abs(rank.score.sorted[flag])
				observed.point[flag] <- hit_tmp/sum(hit_tmp)
			}else{
				hit_tmp <- abs(rank.score.sorted[flag] ** weight)
				observed.point[flag] <- hit_tmp/sum(hit_tmp)      
			}
			RES <- cumsum(observed.point)
			max.RES <- max(RES)
			min.RES <- min(RES)
			es.observed <- signif(ifelse(max.RES>abs(min.RES), max.RES, min.RES), digits=5)
			###########################################################
			if(0){
				## TODO: based on max.RES and min.RES
				es.position <- ifelse(max.RES>abs(min.RES), which.max(RES), which.min(RES))
			}else{
				## TODO: based on sign of ES
				if(res$ES[res$setID==x] > 0){
					es.position <- which.max(RES)
				}else{
					es.position <- which.min(RES)
				}
			}
			###########################################################
			## for leading genes
			if(RES[es.position]<0){
				ind <- which(flag >= es.position)
			}else{
				ind <- which(flag <= es.position)
			}
			hits <- rep(0, length(RES))
			hits[flag] <- 1
			hits[flag[ind]] <- 2
			hits[es.position] <- 3
			
			GeneID <- geneid.sorted
			df_leading <- tibble::tibble(GeneID=GeneID, Rank=1:length(RES), Score=rank.score.sorted, RES=RES, Hits=hits)
			
		})
		names(ls_df_leading) <- res$setID
		
		ls_leadingGenes <- lapply(ls_df_leading, function(x){
			x$GeneID[x$Hits>=2]
		})
		
    }
	
	# leadingGenes
	if(1){
		
		leadingGenes <- lapply(ls_leadingGenes,function(x){
			ind <- match(x, rownames(df_priority))
			rank <- df_priority$rank[ind]
			names(rank) <- x
			return(rank)
		})
		
		## append leading genes
		res$nLead <- sapply(leadingGenes,length)
		## append peak
		res$peak <- sapply(leadingGenes,max)
	}
	
	## append "term_name" and "term_distance"
	summary <- tibble::tibble(setID=res$setID, nAnno=res$setSize, nLead=res$nLead, peak=res$peak, total=rep(nrow(df_priority),length(res$peak)), es=res$ES, nes=res$nES, pvalue=res$pvalue, adjp=res$adjp)
	
	## adjust p-values
	pvalue <- NULL
	summary <- summary %>% dplyr::mutate(adjp=stats::p.adjust(pvalue, method=p.adjust.method))
	
	## scientific notation
	summary$es <- signif(summary$es, digits=3)
	summary$nes <- signif(summary$nes, digits=3)
	summary$pvalue <- signif(summary$pvalue, digits=2)
	summary$pvalue <- ifelse(summary$pvalue<0.01 & summary$pvalue!=0, as.numeric(format(summary$pvalue,scientific=TRUE)), summary$pvalue)
	summary$adjp <- signif(summary$adjp, digits=2)
	summary$adjp <- ifelse(summary$adjp<0.01 & summary$adjp!=0, as.numeric(format(summary$adjp,scientific=TRUE)), summary$adjp)
	
    ################################
    cross <- matrix(0, nrow=length(leadingGenes), ncol=length(leadingGenes))
    if(length(leadingGenes)>=2){
		for(i in seq(1,length(leadingGenes)-1)){
			x1 <- names(leadingGenes[[i]])
			for(j in seq(i+1,length(leadingGenes))){
				x2 <- names(leadingGenes[[j]])
				cross[i,j] <- length(intersect(x1, x2))
				cross[j,i] <- length(intersect(x1, x2))
			}
		}
		colnames(cross) <- rownames(cross) <- names(leadingGenes)
		diag(cross) <- sapply(leadingGenes, length)
    }
    ####################################################################################
	
	# summary_leading
	if(1){
		## df_summary
		nLead <- nAnno <- setID <- name <- value <- lead_member <- NULL
		df_summary <- summary %>% dplyr::mutate(frac=nLead/nAnno)
		
		## df_leading
		df_leading <- leadingGenes[df_summary %>% pull(setID)] %>% tibble::enframe() %>% dplyr::transmute(setID=name, member=value) %>% dplyr::mutate(member_rank=purrr::map_chr(member, function(x) stringr::str_c(names(x),' (',x,')', collapse=', ')))
		
		## df_summary_leading
		df_summary_leading <- df_summary %>% dplyr::inner_join(df_leading, by='setID')
		summary <- df_summary_leading
	}
	
    eGSEA <- list(df_summary = summary,
    			  leading = leadingGenes,
    			  full = ls_df_leading,
    			  cross = cross
                 )
    class(eGSEA) <- "eGSEA"
    
####################################################################################
    endT <- Sys.time()
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    
    if(!silent){
    	message(paste(c("\nEnd at ",as.character(endT)), collapse=""), appendLF=TRUE)
    	message(paste(c("Runtime in total is: ",runTime," secs\n"), collapse=""), appendLF=TRUE)
    }
    
    invisible(eGSEA)
}
