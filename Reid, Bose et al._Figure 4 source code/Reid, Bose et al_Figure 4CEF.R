#=======================================================================================================================================================#
#                                                            CPI Project Figure 4                                                                       #
#=======================================================================================================================================================#

library(ggplot2)
library(ggthemes)
theme_lasso_mse <- function(base_size=14, base_family="arial") {
  (theme_foundation(base_size=base_size, base_family=base_family)+
     theme_classic() +
     theme(  plot.margin=unit(c(5,5,5,5),"mm"),
             plot.title = element_text(face = "bold", size = rel(2), hjust = 0.5),
             # panel.border = element_rect(colour = NA)
             text = element_text(),
             axis.title = element_text( size = rel(2)),
             axis.title.y = element_text(angle=90,vjust =0.5, margin=margin(0,10,0,0,"mm")),
             axis.title.x = element_text(vjust = -0.2, margin=margin(10,0,0,0,"mm")),
             axis.text.x = element_text(angle = 0, hjust = 0.5, size = rel(2)),
             axis.text.y = element_text(angle = 0, hjust = 0.5, size = rel(2)),
             axis.line = element_line(colour="black", size = rel(2)),
             axis.ticks = element_line(size = rel(2)),
             axis.ticks.length = unit(0.5, "cm"),
             # panel.grid.major = element_line(colour="#f0f0f0"),
             # panel.grid.minor = element_blank(),
             legend.key = element_rect(colour = NA),
             legend.position = "right",
             legend.direction = "vertical",
             legend.text = element_text(size = rel(1.75)),
             legend.key.size= unit(0.2, "cm"),
             legend.margin = margin(0,0,10,0, "mm"),
             legend.title = element_blank(),
             strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
             strip.text = element_text(face="bold")
     ))
}  

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# reorder matrix according the values of one column i, in decreasing order
max2min <- function(X, by_i, names_i, decreasing = T){
  K <- X
  if(decreasing){
    K <- K[order(K[,by_i]),]
  }
  else if(!decreasing){
    K <- K[order(-K[,by_i]),]
  } 
  K[,names_i] <-  factor(K[,names_i], levels=unique(K[,names_i]))
  
  return(K)
}


# 1. Read Data ---------------------------------------------------------------------------------------------------------------------------
Baseline_Marrow_Revised <- read.csv('Revised_Baseline_Marrow.csv')


# 2. MIC Calculation (Fig C) -------------------------------------------------------------------------------------------------------------

# Mutual Information Coefficient (Reshef et al. 2011)
# obtain the MIC between the random variable of metabolite rel. intensity values
# and the response vector (1 = CR)

library(minerva)

transp_Baseline_Marrow_Revised <- Baseline_Marrow_Revised
transp_Baseline_Marrow_Revised$Name <- as.character(transp_Baseline_Marrow_Revised$Name)

which(transp_Baseline_Marrow_Revised$Name == "S-Adenosyl-L-Homocysteine") #SAH is identified twice in data
transp_Baseline_Marrow_Revised$Name[358] <- c('S-Adenosyl-L-homocysteine')

rownames(transp_Baseline_Marrow_Revised) <- transp_Baseline_Marrow_Revised$Name
transp_Baseline_Marrow_Revised <- transp_Baseline_Marrow_Revised[,-(1:2)]
transp_Baseline_Marrow_Revised <- t(transp_Baseline_Marrow_Revised)
transp_Baseline_Marrow_Revised <- cbind(transp_Baseline_Marrow_Revised, as.numeric(grepl("CR", rownames(transp_Baseline_Marrow_Revised))))

# calculate MIC values
MIC_BaselineMarrowRevised <- mine(transp_Baseline_Marrow_Revised)

# dataframe for plotting MIC of top MIC(metabolites, Tx)
MIC_BaselineMarrowRevised_Tx <- data.frame(matrix(nrow = nrow(MIC_BaselineMarrowRevised$MIC)-1, ncol = 2))
colnames(MIC_BaselineMarrowRevised_Tx) <- c("Metabolite", "MIC")
MIC_BaselineMarrowRevised_Tx$Metabolite <- as.character(rownames( MIC_BaselineMarrowRevised$MIC[-nrow(MIC_BaselineMarrowRevised$MIC),]))

#reformat metabolite names
MIC_BaselineMarrowRevised_Tx$Metabolite <- stringr::str_to_title(MIC_BaselineMarrowRevised_Tx$Metabolite)

# copy over MIC values for metabolite with response vector
MIC_BaselineMarrowRevised_Tx$MIC <- MIC_BaselineMarrowRevised$MIC[-nrow(MIC_BaselineMarrowRevised$MIC),ncol(MIC_BaselineMarrowRevised$MIC)]

# return ordered dataframe
MIC_BaselineMarrowRevised_Tx <-  max2min(MIC_BaselineMarrowRevised_Tx, 2, 1, decreasing = F)



# 3. Predicting Response from Metabolites (Fig E,F) --------------------------------------------------------------------------------------

library(glmnet)
library(fscaret)
library(pROC)
library(fbroc)
library(sets)

# 3.1. Format 

lasso_Baseline_MarrowRevised <- Baseline_Marrow_Revised
lasso_Baseline_MarrowRevised$Name <- as.character(lasso_Baseline_MarrowRevised$Name)
which(lasso_Baseline_MarrowRevised$Name == "S-Adenosyl-L-Homocysteine") #SAH is identified twice in data
lasso_Baseline_MarrowRevised$Name[358] <- c('S-Adenosyl-L-homocysteine')
rownames(lasso_Baseline_MarrowRevised) <- lasso_Baseline_MarrowRevised$Name
lasso_Baseline_MarrowRevised <- lasso_Baseline_MarrowRevised[,-(1:2)]
lasso_Baseline_MarrowRevised <- t(lasso_Baseline_MarrowRevised)


# 3.2 Build Predictor Sets 

s1 <- set("Succinyl Carnitine", "aspartate",  "(S)-dihydroorotate")
s2 <- set("L-Phenylalanine", "L-glutamate", "3-Hydroxybutyrylcarnitine","3-hydroxy-isovaleryl carnitine", "Glutaryl Carnitine", "L-Asparagine")

combinations <- set_power(s2)
metabolite_sets <- vector("list", 1+set_cardinality(combinations)) # 2 + |combinations| - 1 empty set

metabolite_sets[[1]] <- s1
metabolite_sets[[2]] <- s2

i <- 3
for(s in combinations){
  if(set_cardinality(s) > 0 ){
    metabolite_sets[[i]] <- set_union(s1, s)
    i <- i + 1
  }  
}


# 3.3 Logistic Regression Model 
# Refer to 'Elements of Statistical Learning 2013 - Section 6.6.1 (p. 266)' for methodology

# Leave-One-Out Predictions, with ridge regression since we don't want to shrink variables 

get_metabolite_set_columns <- function(metabolite_combinations, df = lasso_Baseline_MarrowRevised){
  # returns column indices for metabolites in [metabolite_combinations]
  metabolite_columns <- c()
  for(metabolite in metabolite_combinations){
    metabolite_columns <- c(metabolite_columns, which(colnames(lasso_Baseline_MarrowRevised) == metabolite) )
  }
  return(metabolite_columns)
}

# Binary response vector
response_vector <- ifelse(grepl("CR", rownames(lasso_Baseline_MarrowRevised)), 1, 0)

# set seed for reproducibility
set.seed(110)

roc_curves_by_combination <- data.frame()
for(s in metabolite_sets){
  if(set_cardinality(s) > 1){
    metabolite_cols <- get_metabolite_set_columns(s, df = lasso_Baseline_MarrowRevised)
    
    k_fold_preds <- lapply(1:nrow(lasso_Baseline_MarrowRevised), function(x){
      # remove one sample, then fit log reg model with LOO fold cross validation
      fit <- cv.glmnet(lasso_Baseline_MarrowRevised[-x, metabolite_cols, drop = F] ,
                       response_vector[-x],
                       alpha = 0,                # ridge penalty
                       family="binomial",        # logistic regression model
                       nfolds = 17,              # equivalent to leave-one-out
                       type.measure = "class",   # mis-classification error
                       intercept = F,            # model uses only metabolites as predictors
                       standardize = T)          # standardizes predictors before fitting (but coefficients returned on original scale)
      
      # predict holdout sample with fit, using lambda corresponding to min classification error
      pred <- predict(fit$glmnet.fit, lasso_Baseline_MarrowRevised[x, metabolite_cols, drop = F], s= fit$lambda.min ,type="response")
      
      return(data.frame(pred, true = response_vector[x]))
    })
    
    k_fold_preds <- do.call(rbind, k_fold_preds)
    
    roc_metrics <-  roc(response = k_fold_preds$true, predictor = k_fold_preds$X1) 
    temp_df <- data.frame('TPR' =  roc_metrics$sensitivities, 'FPR'  = roc_metrics$specificities,
                          'AUC' = rep(roc_metrics$auc, times = length(roc_metrics$sensitivities) ),
                          'Combination' = rep(paste(as.character(s), collapse=" & "), times=length(roc_metrics$sensitivities) ) )
    roc_curves_by_combination <- rbind(roc_curves_by_combination, temp_df)

    # IF BOOTSTRAPPING:
      # boostrap_roc <- boot.roc(pred = k_fold_preds$X1, true.class = as.logical(k_fold_preds$true), n.boot = 50, use.cache = T )
      # rows are single bootstrap sample
      # boostrap_roc$boot.tpr 
      # boostrap_roc$boot.fpr
      # boostrap_roc$auc
      # boostrap_roc$roc
      
  }
}
write.csv(roc_curves_by_combination, 'ROC_Curves_all_combos.csv')
# ----------------------------------------------------------------------------------------------------------------------------------------