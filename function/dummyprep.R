
library(fastDummies)

dummyprep <- function(dat, numeric.break = 5, remove_first_dummy=T){
    
    new.dat <- dat
    
    for (i in 1:ncol(dat)){
        
        if (length(unique(unlist(dat[,i])))>12){
            
            bb <-  unique(quantile(unlist(dat[,i]), probs = seq(0,1,1/numeric.break), na.rm = T))
            j <- cut(unlist(dat[,i]), breaks = bb,
                     labels = 1:(length(bb)-1))
            new.dat[,i] <- j
            
        }else{
            
            new.dat[,i] <- factor(unlist(dat[,i]))
        }
        
    }
    
    new.dat <- dummy_columns(new.dat,remove_selected_columns = TRUE, 
                             remove_first_dummy = remove_first_dummy,
                             ignore_na = TRUE)
    for (k in 1:ncol(new.dat)){
        new.dat[which(is.na(new.dat[,k])),k]<- 0
    }
    return(new.dat)
}

dummyprep.factor <- function(dat, numeric.break = 5, remove_first_dummy=T){
     
    factor_var <- NULL
    
    for (i in 1:ncol(dat)){
        
        if (length(unique(unlist(dat[,i])))<=12){
             
            factor_var <- c(factor_var,i)
        }
    }
    
    new.dat <- dat[,factor_var]
    new.dat2 <- dat[,-factor_var]
    for (j in 1:length(factor_var)){
        new.dat[,j] <- factor(new.dat[,j])
    }
    
    new.dat <- dummy_columns(new.dat,remove_selected_columns = TRUE, 
                             remove_first_dummy = remove_first_dummy,
                             ignore_na = TRUE)
    for (k in 1:ncol(new.dat)){
        new.dat[which(is.na(new.dat[,k])),k]<- 0
    }
    new.dat <- cbind(new.dat2,new.dat)
    return(new.dat)
}
