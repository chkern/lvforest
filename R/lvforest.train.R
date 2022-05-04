#' Estimating Valid Latent Variable Scores from Conditionally Causal Models
#'
#' @param data Dataset to be analyzed. Note that the observed responses of your model cannot contain missing data. 
#' @param model Model description in lavaan-terminology.
#' @param input Character vector of partitioning variables. Note that every partitioning variable must either be defined as factor (if it should be treated as categorical) or numeric.
#' @param ordered Character vector of ordinal observed variables of the model. Default = NULL.
#' @param ntrees Number of trees to be computed. Default = 100.
#' @param split Number of partitioning variables to be selected for random split selection at every split point within a tree ('mtry' argument in ctree function). Default = 2.
#' @param minsize Minimum size of terminal nodes of trees. Needs to be large enough to estimate the model parameters. Default = 300. 
#' @param cutoff_rmsea Cutoff value for model fit evaluation. Default = .05.
#' @param cutoff_loading Cutoff value for multiplicative parameters in model. Note that models with parameter estimates not significantly different form 0 are excluded per default. Default = 0.
#' @param dbsamp Double sampling. Default = TRUE.
#' @param bagging Numeric vector defining the proportion of original data used for tree growing (first number in vector) and re-fitting terminal nodes (second number in vector). Default = NULL.
#' @param std.lv Standardization of latent variable variances in model. Default = FALSE.
#' @param ctree_control Control settings for ctree function. Default = ctree_control(minbucket = minsize, mtry = split).
#' 
#' @examples 
#' \dontrun{
#' trained <- lvforest.train(
#' input = c("conf","dicho1","cat1","num1","cat2","rand1","rand2","rand3","rand4","rand5"),
#' model <- 'LatVar1 =~ simuvar1 + beta2*simuvar2 + beta3*simuvar3 
#'  LatVar2 =~ simuvar4 + beta5*simuvar5 + beta6*simuvar6
#'  LatVar3 =~ simuvar7 + beta8*simuvar8 + beta9*simuvar9
#'  LatVar4 =~ LatVar1 + lambda2*LatVar2 + lambda3*LatVar3 + delta*simuvar_effect',
#' ordered = c("simuvar1","simuvar2","simuvar3","simuvar4",
#' "simuvar5","simuvar6","simuvar7","simuvar8","simuvar9"),
#' data = simu,
#' cutoff_rmsea = .03,
#' ctree_control = ctree_control(minbucket = 200, mtry = 1, testtype = "Teststatistic") 
#' #Bonferroni correction off 
#' )
#' }

#' @export
lvforest.train <- function(data,model,input,ntrees=100,split=2,minsize=300,cutoff_rmsea=.05,cutoff_loading=.2,dbsamp=T,bagging=NULL,ordered=NULL,std.lv=FALSE,ctree_control=partykit::ctree_control(minbucket=minsize, mtry= split)){
  
  #### Trees berechnen
  trees=lvtrees(data,model,input,ntrees,cutoff_rmsea,cutoff_loading,dbsamp,bagging,ordered,std.lv,ctree_control)
  
  #### Liste benennen
  for(i in 1:ntrees){
    names(trees)[[i]] <- paste0("iteration",i)
    if(length(trees[[i]])==1){names(trees[[i]]) <- paste0("tree",i)}
    if(length(trees[[i]])>1){names(trees[[i]]) <- c(paste0("tree",i),paste0("fit",i),paste0("modelfit",i),paste0("DataGrow",i),paste0("DataRefit",i))}
  }
  
  #### Bestnodes erstellen 
  mofi <- data.frame()
  for(i in 1:length(trees)){
    if ( length(trees[[i]])>1 && !is.character(trees[[i]][[3]]) ) mofi <- rbind(mofi,trees[[i]][[3]])
  }
  mofi = na.exclude(mofi[colSums(!is.na(mofi)) > 0]) 
  mofi_stable = mofi[mofi$RMSEA<cutoff_rmsea & mofi$`parameter stability`=="stable",]
  
  if(!(ncol(mofi_stable)==0)){
    #### trees "bereinigen"
    pred_trees <- trees[ as.numeric( unique(mofi_stable$tree))  ]
    names(pred_trees) <- names( trees[ as.numeric( unique(mofi_stable$tree)) ] )
    mofi_stable = mofi_stable[order(mofi_stable$RMSEA),]
  }
  suc_trees <- trees[sapply(trees, function(y) length(y)==5)]
  names(suc_trees) <- names(trees[sapply(trees, function(y) length(y)==5)])
  
  ##### Information
  if(!dbsamp & is.null(bagging)){direct=TRUE} else {direct=FALSE}
  info=list()
  info[[1]] <- ntrees; names(info)[[1]] <- "ntrees"
  info[[2]] <- length(suc_trees); names(info)[[2]] <- "successful_iterations"
  if(!(ncol(mofi)==0)){info[[3]] <- length(pred_trees)} else {info[[3]] <- 0}; names(info)[[3]] <- "pred_iterations" 
  info[[4]] <- input; names(info)[[4]] <- "input"
  info[[5]] <- direct; names(info)[[5]] <- "direct"
  info[[6]] <- dbsamp; names(info)[[6]] <- "dbsamp"
  if(is.null(bagging)){info[[7]] <- "NULL"} else {info[[7]] <- bagging}; names(info)[[7]] <- "bagging"
  info[[8]] <- cutoff_rmsea; names(info)[[8]] <- "cutoff_rmsea"
  info[[9]] <- cutoff_rmsea; names(info)[[9]] <- "cutoff_loading"
  info[[10]] <- ordered; names(info)[[10]] <- "ordered"
  info[[11]] <- model; names(info)[[11]] <- "model"
  
  #### Results erstellen
  lvfor <- list()
  lvfor[[1]] <- info; names(lvfor)[[1]] <- "Info"
  lvfor[[2]] <- mofi; names(lvfor)[[2]] <- "Nodes"
  if(!(ncol(mofi_stable)==0)){lvfor[[3]] <- pred_trees; names(lvfor)[[3]] <- "Pred_trees"}
  if(!(ncol(mofi_stable)==0)){lvfor[[4]] <- mofi_stable; names(lvfor)[[4]] <- "Bestnodes"}
  lvfor[[5]] <- suc_trees; names(lvfor)[[5]] <- "Suc_trees"
  
  return(lvfor)
}


################################################################################
### Support functions

`%>%` <- dplyr::`%>%`
`%dorng%` <- doRNG::`%dorng%`
`%dopar%` <- foreach::`%dopar%`

lvdata <- function(data,model,input){
  fitsi <- try(lavaan::cfa(model = model, data = data, estimator = "ML",do.fit=FALSE))
  manifs <- fitsi@Model@dimNames[[1]][[1]] #manifest variables
  latvars <- fitsi@Model@dimNames[[1]][[2]] #latent variables
  data = as.data.frame(data)
  data = data[,c(input,manifs)] %>% tidyr::replace_na(list("-99"))
  return(list(manifs,latvars,data))
}

lvform <- function(scores,input){
  formy = "";  for (i in 1:length(scores)){formy=paste(formy, paste(scores[i],"+"))  }
  formy = substr(formy,1,nchar(formy)-2); form = paste(formy,"~")
  for (i in 1:length(input)){ form=paste(form,paste(input[i],"+")) }
  form = substr(form,1,nchar(form)-2)
  return(form)
}

lvsampling <- function(data,dbsamp,bagging){
  if (dbsamp){
    folds <- caret::createFolds(y = rownames(data), k= 2)
    folds1 <- folds$Fold1
    folds2 <- folds$Fold2
  } else if (!is.null(bagging) && bagging < 1) {
    folds1=folds2 <- sample(x=rownames(data),size=bagging*nrow(data))
  } else{folds1=folds2=rownames(data)} 
  return(list(folds1,folds2))
}

lvscores <- function(data,model,std.lv){
  csresult <- tryCatch({
    fit_num = lavaan::cfa(model = model, data = data,  estimator = "ML",std.lv = std.lv)
    fit_num_scores <- lavaan::lavScores(fit_num)
    colnames(fit_num_scores) = stringr::str_replace_all(colnames(fit_num_scores), "[^[:alnum:]]", "")
    scores = colnames(fit_num_scores)
    data = cbind(data,fit_num_scores)
    list(data,scores)},error=function(e){warning("lvforest warning: Your model does not converge with ML estimator. Try 'dbsamp=TRUE'.");stop(e)})
  return(csresult)
}

lvtrees <- function(data,model,input,ntrees,cutoff_rmsea,cutoff_loading,dbsamp,bagging,ordered,std.lv,ctree_control){
  #Preprocess & gather Info
  direct=FALSE
  dt=lvdata(data,model,input);manifs=dt[[1]];latvars=dt[[2]];data=dt[[3]]
  
  #Scores ausrechnen wenn "direct"
  if(dbsamp) bagging=NULL
  if(!dbsamp & is.null(bagging)){direct=T}
  if(direct){sc=lvscores(data,model,std.lv);data=sc[[1]];scores=sc[[2]]}
  
  #All Iterations Multicore
  ncores <- parallel::detectCores()-1
  cl <- parallel::makeCluster(spec=ncores) 
  doParallel::registerDoParallel(cl)
  trees <- foreach::foreach(j=1:ntrees, .packages=c("stringr","lavaan","caret","partykit","strucchange"),.export=c("lvdata","lvform","lvsampling","lvscores")) %dorng% { 
    
    #### Data Partitioning
    sp = lvsampling(data,dbsamp,bagging)
    treedata = data[which(as.numeric(rownames(data)) %in%  sp[[1]]),]
    datafit = data[which(as.numeric(rownames(data)) %in%  sp[[2]]),]
    
    #### MOB mit der ersten Haelfte der Daten
    #### Scores ausrechnen wenn nicht "direct"
    try({
      if(!direct){sc=lvscores(data=treedata,model,std.lv);treedata=sc[[1]];scores=sc[[2]]}
      tree <- partykit::ctree(as.formula(lvform(scores,input)),treedata, control = ctree_control  )
    },silent = T) 
    ####
    
    if(exists("tree")){ #erster error handler
      
      #### Re-Fit & Modelfit table
      fit_ord <- list()
      ni <- partykit::nodeids(tree, terminal = TRUE)
      modelfit <- data.frame()
      rls <- partykit:::.list.rules.party(tree)
      types <- sapply(input, function(y) ifelse(class(data[,y]) == "factor","LMuo","maxLM")  )#welcher parameter-stablilitaetstest
      typenames <- names(types)
      
      if(!(sum(rls=="")>0) & length(ni)!=0){ #zweiter Error hanlder
        
        ###Nodes refitten und Ergebnisse in Tabelle schreiben
        for(i in 1:length(ni)){
          mf<-c()
          data_refit <- subset(datafit,eval(parse(text=rls[i])) ) #richtiger Terminal node
          fit_ord[[i]] <- tryCatch({lavaan::cfa(model = model, data = data_refit, start = start, ordered = ordered, estimator = "WLS", std.lv = std.lv, control=list(iter.max=100))},error=function(e){return("n.c.")},warning=function(e){return("n.c.")})
          names(fit_ord)[[i]] <- paste0("node",ni[i])
          if(is.character(fit_ord[[i]])){next}
          pars <- tryCatch({cbind(fit_ord[[i]]@ParTable$lhs,fit_ord[[i]]@ParTable$op,fit_ord[[i]]@ParTable$rhs,fit_ord[[i]]@ParTable$est,fit_ord[[i]]@ParTable$se)},error=function(e){return("n.c.")})
          if(any(pars!="n.c.")) {
            parvars <- as.numeric(pars[pars[,2]=="~~",4]) 
            parload <- pars[pars[,2]=="=~",c(4,5)]
            zv <- as.numeric(parload[,1])/as.numeric(parload[,2]) #z-values 
            zvp <- sapply(zv[zv!=Inf],function(y) pnorm(q=y, lower.tail=FALSE)) #pvalues
          } else {next} #nur Varianzen & factor loadings
          if(any(parvars<0 | parvars>30) | any(zvp>0.05 & zvp!=1) | any(abs(as.numeric(parload[,1]))<cutoff_loading) ){next} #negative Varianzen oder zu kleine factor loadings --> next iteration
          
          mf[1] <- j; mf[2] <- ni[i]; if(dbsamp){mf[3] <- nrow(subset(treedata,eval(parse(text=rls[i])) ));mf[4] <- fit_ord[[i]]@Data@nobs[[1]]} else { mf[3] <- nrow(data_refit); mf[4] <- NA}
          mf[5] <- lavaan::fitMeasures(fit_ord[[i]],"rmsea");mf[6] <- lavaan::fitMeasures(fit_ord[[i]],"rmsea.ci.lower");mf[7] <- lavaan::fitMeasures(fit_ord[[i]],"rmsea.ci.upper");mf[8] <- lavaan::fitMeasures(fit_ord[[i]],"pvalue"); mf[9] <- rls[which(names(rls)==ni[i])]
          
          ### Parameter stability 
          if(lavaan::fitMeasures(fit_ord[[i]],"rmsea")<cutoff_rmsea){  #stability test nur für models mit gutem fit
            datastable <- subset(data,eval(parse(text=rls[i])) ) #richtiger Terminal node aber mit gesamten Daten!
            fit_num <- tryCatch({lavaan::cfa(model = model, data = datastable, estimator = "ML", std.lv = std.lv, control=list(iter.max=100))},error=function(e){return("n.c.")}) 
            if(!is.character(fit_num)) {
              fluc_tests <- tryCatch({mapply(function(x,y) { strucchange::sctest(fit_num, order.by = datastable[,y],   functional = x)$p.value }, types, typenames  )},error=function(e){return("n.c.")}) #was wenn "solution has not been found"?
              if(is.character(fluc_tests)){mf[10] <- "n.c."} else {mf[10] <- tryCatch({ifelse(any(fluc_tests < (0.05/length(fluc_tests)) & fluc_tests != 0),"unstable","stable")},error=function(e){return("n.c.")})  }
            } else {mf[10] <- "n.c."} #"not converged"
          } else {mf[10] <- "n.e."} #"not execuded"
          
          modelfit <- rbind(modelfit,mf)
        }
        
        if(nrow(modelfit)==0){modelfit<-"none converged"} else {
          modelfit <- as.data.frame(modelfit)
          colnames(modelfit) <- c("tree","node","n_fit","n_refit","RMSEA","RMSEA C.I. lower","RMSEA C.I. upper","p-value_chisq","decision_rule","parameter stability")
          modelfit <- modelfit[order(modelfit[,5],modelfit[,6]),] 
        }
        list(tree,fit_ord,modelfit,sp[[1]],sp[[2]]) 
      } else {tree}
    } else {NA}
  }
  parallel::stopCluster(cl)
  return(trees)
}

