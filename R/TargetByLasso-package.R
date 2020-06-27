#' TargetByLasso
#' @description This package is used for key target prediction of a disease by employing the results of drug test matrix
#' @param single.effect Single drugs with effects that can be binomial or continue values
#' @param combo.effect Combination drugs with effects, only allowing binomial here
#' @param ref.matrix A matrix from STITCH, containing all drug-target information
#' @param single.select.times Repeat times for single feature selection, the default value is 100
#' @param single.nfold Fold number for single selection cross-validation, the default value is 10
#' @param combo.select.times Repeat times for combination features selection, the default value is 10
#' @param combo.nfold Fold number for combination features selection, the default value is 10
#' @param model The function has two models for combo selection 1 = glmnet, 2 = biglasso, the default is glmnet
#' @param n.cores Number of cores of CPU asking for tasks running
#' @importFrom reshape2 acast
#' @importFrom glmnet cv.glmnet
#' @importFrom parallel detectCores
#' @importFrom doParallel registerDoParallel
#' @importFrom doParallel stopImplicitCluster
#' @importFrom biglasso cv.biglasso
#' @keywords systemtic diseases, key target prediction, target relationship prediction
#' @export
#"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
TargetPre <- function(single.effect,             # single drugs with effects
                             combo.effect,              # combination drugs with effects
                             ref.matrix,                # a matrix containing all drug-target information
                             single.select.times = 100, # repeat times for single feature selection
                             single.nfold        = 10,  # fold number for single selection cross-validation
                             combo.select.times  = 10,  # repeat times for combination features selection
                             combo.nfold         = 10,  # fold number for combination features selection
                             model               = 1,   # we have to models for combo selection (1 = glmnet, 2 = biglasso)
                             n.cores = detectCores()    # number of cores of CPU asking for tasks running
){

  # Assigning proper names for the read-in data frame.
  colnames(single.effect) <- c("Drug","effect")
  colnames(combo.effect)   <- c("Drug.1","Drug.2","effect")
  names(ref.matrix)[1]    <- "Drug"

  # Organizing single effects (d) and reference matrix (dt) to generate a matrix for single feature selection (d.dt).
  model.sin           <- merge(single.effect, ref.matrix, by = "Drug")
  rownames(model.sin) <- model.sin$Drug
  model.sin           <- model.sin[,-1]
  model.sin           <- model.sin[,colSums(abs(model.sin)) != 0]
  d.dt                <- cbind(model.sin[1],data.frame(t(unique(t(model.sin[,-1])))))

  # Single feature selection
  df.single.target <- NULL

  registerDoParallel(n.cores)
  for (i in 1:single.select.times) {

    if(length(table(single.effect$effect)) == 2){
      # For the single effects are binomial: effective or not.
      cvfit.single <- cv.glmnet(as.matrix(d.dt[,-1]),
                                d.dt$effect,
                                family = "binomial",
                                type.measure = "class",
                                parallel = T,
                                nfolds = single.nfold)
    } else {
      # For the single effects are continuous variables.
      cvfit.single <- cv.glmnet(as.matrix(d.dt[,-1]),
                                d.dt$effect,
                                parallel = T,
                                nfolds = single.nfold)}

    t = as.matrix(coef(cvfit.single,s="lambda.min"))
    df.single.target = cbind(df.single.target, t)
  }
  stopImplicitCluster()

  # Modifying the outcome, and removing duplicated results.
  df.single.target <- df.single.target[-1,]
  df.single.target <- df.single.target[rowSums(df.single.target) != 0,]
  df.single.target <- df.single.target[,colSums(df.single.target) != 0]
  df.single.target <- t(unique(t(df.single.target)))

  # Aggregating the features sharing the same pattern in the original matrix, and getting all the features' names
  same.single.list = list()
  for (i in 1:nrow(df.single.target)) {
    same.single.list[[i]] <- which(apply(model.sin[,-1], 2, function(x) return(all(x == model.sin[,colnames(model.sin) ==  rownames(df.single.target)[i]])
    )))
  }
  target.all <- names(unlist(same.single.list))

  # Making the combination drug effects into a symmetric matrix, removing unnecessary drugs,
  # and reformatting it into drug-drug effects.
  com.mt      <- acast(combo.effect, Drug.1 ~ Drug.2)
  com.mt      <- com.mt[rownames(com.mt) %in%
                          intersect(intersect(as.character(combo.effect$Drug.1),
                                              as.character(combo.effect$Drug.2)),
                                    as.character(ref.matrix$Drug)),
                        colnames(com.mt) %in%
                          intersect(intersect(as.character(combo.effect$Drug.1),
                                              as.character(combo.effect$Drug.2)),
                                    as.character(ref.matrix$Drug))]
  dd          <- melt(com.mt)
  dd          <- na.omit(dd)
  dd$dd.names <- paste(dd$Var1,dd$Var2, sep = ".")

  # According to the newly generated dd names and selected features, screening the corresponding reference
  # matrix content and getting the Kronecker product
  ref.df           <- ref.matrix[ref.matrix$Drug %in% unique(as.character(dd$Var2,dd$Var1)),]
  ref.df           <- ref.df[,colnames(ref.matrix) %in% c("Drug",target.all)] # target.all gives best accuracy
  rownames(ref.df) <- ref.df$Drug
  ref.df           <- ref.df[-1]
  ref.mt.kronecker <- kronecker(as.matrix(ref.df), as.matrix(ref.df))

  # Naming rows and select proper dd pairs from the screened matrix.
  rown                       <- as.character(rownames(ref.df))
  drug.combn                 <- data.frame(v1 = rep(rown, length(rown)), v2 = rep(rown, each = length(rown)))
  rownames(ref.mt.kronecker) <- paste(drug.combn$v1, drug.combn$v2, sep = ".")
  ref.mt                     <- ref.mt.kronecker[rownames(ref.mt.kronecker) %in% dd$dd.names,]

  # Naming columns and removing duplicated content, like a-b, and b-a.
  coln                        <- as.character(colnames(ref.df))
  target.combn                <- data.frame(v1 = rep(coln, length(coln)),v2 = rep(coln, each = length(coln)))
  target.combn$tt.name        <- paste(target.combn$v1, target.combn$v2, sep = ".")
  colnames(ref.mt)            <- target.combn$tt.name
  target.combn.unique         <- data.frame(t(combn(coln,2)))
  target.combn.unique$tt.name <- paste(target.combn.unique$X1, target.combn.unique$X2, sep = ".")
  ref.mt                      <- ref.mt[,colnames(ref.mt) %in% target.combn.unique$tt.name]

  # Removing columns with contents are all 0
  ref.mt <- ref.mt[,colSums(abs(ref.mt)) != 0]

  # Organizing combo effects (dd) and reference matrix (ddtt) to generate a matrix for single feature selection (dd.ddtt).
  ref.mt               <- data.frame(ref.mt)
  ref.mt$dd.names      <- rownames(ref.mt)
  ref.mt.sub           <- merge(dd[c(3,4)], ref.mt, by = "dd.names")
  rownames(ref.mt.sub) <- ref.mt.sub$dd.names
  ref.mt.sub           <- ref.mt.sub[,-1]
  ref.mt.sub           <- ref.mt.sub[,colSums(abs(ref.mt.sub)) != 0]
  dd.ddtt              <- cbind(ref.mt.sub[1],data.frame(t(unique(t(ref.mt.sub[-1])))))

  df.target.pair = NULL

  if(model == 1){
    registerDoParallel(n.cores)
    for (i in 1:combo.select.times) {
      cvfit.combo <- cv.glmnet(as.matrix(dd.ddtt[-1]),
                               dd.ddtt$value,
                               family = "binomial",
                               type.measure = "class",
                               nfolds = combo.nfold,
                               parallel=TRUE)
      t = as.matrix(coef(cvfit.combo, s = "lambda.min"))
      df.target.pair = cbind(df.target.pair, t)}
    stopImplicitCluster()
  } else if (model == 2) {
    for (i in 1:combo.select.times) {
      cvfit.combo <- cv.biglasso(as.big.matrix(as.matrix(dd.ddtt[-1])),
                                 dd.ddtt$value,
                                 nfolds = combo.nfold,
                                 ncores = n.cores)
      t = as.matrix(coef(cvfit.combo, s = "lambda.min"))
      df.target.pair = cbind(df.target.pair, t)}}


  # Modifying the outcome, and removing duplicated results.
  df.target.pair <- df.target.pair[-1,]
  df.target.pair <- df.target.pair[rowSums(df.target.pair) != 0,]
  df.target.pair <- df.target.pair[,colSums(df.target.pair) != 0]
  df.target.pair <- t(unique(t(df.target.pair)))

  # Aggregating the feature pairs sharing the same pattern in the original matrix,
  # and getting all the feature pairs' names
  same.pair.list = list()
  if (length(df.target.pair) ==0){
    cat("No combo features selected, please try another model. \n")
  } else {
    for (i in 1:nrow(df.target.pair)) {
      same.pair.list[[i]] <- which(apply(ref.mt.sub[,-1], 2, function(x) return(all(x == ref.mt.sub[,colnames(ref.mt.sub) == rownames(df.target.pair)[i]]))))}
  }



  # Gathering all the expecting outcomes as a result.
  result = list(single.selected.targets                 = df.single.target,
                single.targets.with.duplicated.patterns = same.single.list,
                combo.selected.targets                  = df.target.pair,
                combo.targets.with.duplicated.patterns  = same.pair.list)
  return(result)
}
