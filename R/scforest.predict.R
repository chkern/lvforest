#' Estimating Valid Latent Variable Scores from Conditionally Causal Models
#'
#' @param trained Trained model created by scforest.train(). 
#' @param data Dataset to be used for latent variable score prediction. Note that the observed responses of your model cannot contain missing data and that it must contain all partitioning variables from training.
#' @param idvar ID-Variable to distinguish individual data points in the dataset.
#' @param conf Character vector of potential confounding variables. Default = 'Partitioning Variables from training'.
#' @param latvar Latent variable of interest to estimate latent variable scores. 
#' @param indmethod Method for testing for significance of dHSIC. Default = gamma.
#' @param exclude_unconf Only use latent variable scores if there is unconfoundedness. Default = TRUE.
#' @param stdscores Latent variable scores as z-scores.
#' @examples 
#' \dontrun{
#' predicted <- scforest.predict(trained = trained,
#' data = simu, idvar = "id", latvar = "LatVar4", exclude_unconf = T, stdscores = F)
#' }

#' @export
scforest.predict <- function(trained,data,idvar,conf=trained$Info$input,latvar=latvars[1],indmethod="gamma",exclude_unconf=TRUE,stdscores=FALSE){
  ### Warnung
  if(nrow(trained$Bestnodes) == 0){stop("scforest error: no suitable nodes produced in training.")}
  
  #### Gather Information
  data = as.data.frame(data)
  model = trained$Info$model
  ordered = trained$Info$ordered
  manifs = trained[[3]][[1]][[2]][[  which(sapply(trained[[3]][[1]][[2]], function(y) class(y)) == "lavaan")[1]  ]]@Model@dimNames[[1]][[1]]
  latvars = trained[[3]][[1]][[2]][[ which(sapply(trained[[3]][[1]][[2]], function(y) class(y)) == "lavaan")[1]  ]]@Model@dimNames[[1]][[2]]
  bestnodes = trained$Bestnodes[order(as.numeric(trained$Bestnodes$tree),as.numeric(trained$Bestnodes$node)),]
  tree = unique(bestnodes$tree)
  
  #### Pick models used for prediction out of trained model
  bstfts = sapply(trained[[3]], function(y) y[[2]]);bestfits = list()
  for(i in 1:length(unique(bestnodes$tree))){  bestfits[[i]] = bstfts[[i]][  which( partykit::nodeids(trained[[3]][[i]][[1]], terminal = TRUE)  %in%  bestnodes[bestnodes$tree==tree[i],2])  ] }
  names(bestfits) = names(trained[[3]])
  
  ### Warnungen
  if(length(bestfits) != length(tree)  ){stop("scforest fatal error: length of best_fits in trained model does not correspond with information in bestnodes.")} #weiter zum nächsten tree wenn liste mit best fits und bestnodes nicht übereinstimmen!
  if( is.null(data) ){stop("scforest error: please define data set for prediction including same manifest variables and partitioning variables as data set for training.")}
  if( is.null(latvar) ){stop("scforest error: please define latent variable of interest for LV-scores.")}
  
  ### Prediction & Independence tests...
  predis <- scpredtree(data,bestnodes,bestfits,tree,idvar,manifs,latvar,conf,indmethod,exclude_unconf,stdscores)
  
  ### Compile results
  scfor <- scpredcompile(predis,data,idvar,latvar,tree)
  return(scfor)
}
################################################################################
### Support functions

`%>%` <- dplyr::`%>%`
`%dopar%` <- foreach::`%dopar%`

scpredconf <- function(i,j,preds,latvar,conf,bestnodes,tree,indmethod){ #irgendwelche Probleme mit dopar... alles exportiert?
  try({
    types <- sapply(conf, function(y) ifelse(class(data[,y]) == "factor","discrete","gaussian")  )
    typenames <- names(types)
    
    if(nrow(preds)>0){
      varlist <- lapply(names(types),function(y){as.double(preds[,y])})
      ind_tests <- mapply(function(x,y) { dHSIC::dhsic.test(Y=as.matrix(preds[,latvar]),X=as.matrix(varlist[[which(names(types)==y)]]),matrix.input=F,method=indmethod,pairwise = F,kernel=c("gaussian",x))$p.value }, types, typenames  )
    } else {ind_tests <- NA}
    unconf_info<- c()
    unconf_info[1] <- tree[j]
    unconf_info[2] <- bestnodes[bestnodes$tree==tree[j],"node"][i]
    if(nrow(preds)>0){unconf_info[3] <- if(any( ind_tests < (0.05/length(ind_tests)) )){"confounded"} else {"not confounded"} } else {unconf_info[3] <- "no predictions"}
    unconf_info[4] <- bestnodes[bestnodes$tree==tree[j],"decision_rule"][i]
  },silent=T )#unabhängigkeit von allen input variablen (multivariate distribution)
  return(list(ind_tests,unconf_info))
}

scpredtree <- function(data,bestnodes,bestfits,tree,idvar,manifs,latvar,conf,indmethod,exclude_unconf,stdscores){
  ncores <- parallel::detectCores()-1
  cl <- parallel::makeCluster(spec=ncores) 
  doParallel::registerDoParallel(cl)
  predis <- foreach::foreach(j=1:length(tree), .packages=c("lavaan","dHSIC","dplyr"), .export=c("scpredconf")) %dopar% {
    ind_tree=list()
    totpreds=data.frame()
    for(i in 1:length(bestnodes[bestnodes$tree==tree[j],2]) ){
      data_node <- data %>% dplyr::filter(eval(parse(text=  bestnodes[bestnodes$tree==tree[j],"decision_rule"][i]  )) ) %>% dplyr::select(all_of(c(idvar,manifs,conf)))
      preds <- tryCatch({ 
        scores <- lavaan::lavPredict(bestfits[[j]][[i]],newdata = data_node)
        if(stdscores){scoressd <- apply(scores, 2, sd);scores <- t(  apply(scores, 1, function(x) x/scoressd)  )} 
        cbind(data_node,scores)   
      }, error=function(cond){return(data.frame())})
      
      ###Independence-test für jeden Node
      ind_tree[[i]] <- scpredconf(i,j,preds,latvar,conf,bestnodes,tree,indmethod)
      
      ###Preds nur verwenden wenn parameter unconfounded
      if(exclude_unconf){
        if( !(ind_tree[[i]][[2]][3] %in%c("no predictions","confounded"))){ totpreds <- rbind(totpreds,preds) }   
      } else { totpreds <- rbind(totpreds,preds) }
    }
    list(totpreds,ind_tree)
  }
  parallel::stopCluster(cl)
  return(predis)
}

scpredcompile <- function(predis,data,idvar,latvar,tree){
  goodpreds <- as.data.frame(data[,idvar]);colnames(goodpreds)=idvar;v=1
  for (i in 1:length(predis) ){ #every tree
    if(nrow(predis[[i]][[1]])>0) {preds <- predis[[i]][[1]][,c(idvar,latvar)]} else {next};v=v+1
    diffid <- setdiff(data[,idvar],preds[,idvar])
    diffid <- as.data.frame(cbind(diffid,rep(NA, length(diffid)) ))
    colnames(diffid) <- colnames(preds)
    preds <- as.data.frame(rbind(preds,diffid))
    if(nrow(preds) == nrow(data)  ){ #error handler
      preds <- preds[order(match(preds[,idvar], data[,idvar] )), ]
      goodpreds <- cbind(goodpreds,preds[,latvar])
      colnames(goodpreds)[v]<-paste0("tree",tree[i])
    }
  }
  goodpreds <- as.data.frame(goodpreds)
  if(ncol(goodpreds)>2){goodpreds$mean <- rowMeans(goodpreds[,2:ncol(goodpreds)], na.rm = TRUE)}
  
  unconf_table <-  data.frame()  
  for(i in 1:length(tree)){
    for(j in 1:length(predis[[i]][[2]])){
      unconf_var <- names(predis[[i]][[2]][[j]][[1]])[which(as.numeric(predis[[i]][[2]][[j]][[1]]) < (0.05/length(predis[[i]][[2]][[j]][[1]]))  )]
      if( !(identical(unconf_var, character(0)) |  is.null(unconf_var))  ) {unconf_row <- c(predis[[i]][[2]][[j]][[2]], unconf_var )} else {unconf_row <- c(predis[[i]][[2]][[j]][[2]], NA )}
      unconf_table <- rbind(unconf_table, unconf_row) 
    }
  }
  colnames(unconf_table) <- c("tree","node","confounded params","decision_rule","confounder")
  
  gathered <- list()
  gathered[[1]] <- goodpreds
  names(gathered)[[1]] <- "Goodpreds"
  gathered[[2]] <- unconf_table
  names(gathered)[[2]] <- "Unconfoundedness_table"
  return(gathered)
}







