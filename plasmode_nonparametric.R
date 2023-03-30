#### Arguments
# data:         A dataframe containing the real data to be used to generate simulated datasets (data.frame)
# size:         The number of observations to be generated (numeric)
# n:            Sample size of the original data (numeric)
# nsim:         The number of simulated datasets to be generated
# timeobs:      Number of time-points of data (numeric)
# Y:            Name of the outcome variable (character)
# id:           Name of the id variable (character)
# A1, A2, ...:  Names of the binary treatment variables (character)
# L1, L2, ...:  Names of the covariates (vector of character)
# silent:       TRUE/FALSE to indicate whether or not the function should provide update on the number of datasets generated

#### Notes
# Categorical variables should be factors to be treated as such by the function

#### Output
# A list of nsim simulated datasets


############################################Non parametric plasmode###########################################

library(randomForest);
plasmode_longitudinal_non_param<- function(data, size, n, nsim=1, timeobs, Y, id, A1, A2, A3,
                                           L1, L2, L3, silent = FALSE,...){
  
  ### Arguments check
  if(!is.data.frame(data)) stop("data needs to be a data frame")
  if(!is.numeric(n)) stop("n needs to be a numeric variable")
  if(length(n) != 1) stop("n needs to be of length 1")
  if(n < 0) stop("The sample size should be a positive number")
  if((n %% 1) != 0){warning("The sample size was not an integer and has been rounded down"); n = n - n %% 1;}
  if(!is.numeric(timeobs)){stop("timeobs needs to be a numeric variable")}
  if(length(timeobs) != 1){stop("timeobs needs to be of length 1")}
  if(timeobs<=1) stop("timeobs should be greater to 1")
  if((timeobs %% 1) != 0){warning("timeobs was not an integer and has been rounded down"); timeobs = timeobs - timeobs %% 1;}
  if (is.null(id)) stop("An identifier variable 'id' needs to be provided")
  if (is.null(size)) stop("the size of the simulated sample (size) needs to be provided")
  if (nsim<0 || is.null(nsim)) stop("the number of simulated samples needs to be provided and to be a positive number")

  
  ids= data.frame(matrix(nrow = size, ncol = nsim))
  list_data_plasmode <- list() ## ouput of the function
  for(i in 1:nsim){
    idxs <- sample(n, size, replace = TRUE); ### Generate the variable of identifiers of the simulated data
    
    ids[,i] <- data[[id]][idxs]; 
    
    dat.plasmode= data[idxs,]; ### create the data frame of the simulated data
    
    
    ## This is the simulation of A1 (first treatment which is a binomial variable) from the regression A1~L1
    groupvar=c(get(apply(expand.grid("L",1), 1, paste0, collapse=""))) ##this is L1: independent variable used in the model of regression
    
    formula0=as.formula(paste(mget(paste0("A",1)),paste(groupvar, collapse=" + "), sep=" ~ "))  ## the formula of regression A1~L1
    
    
    ## assign the formula to a variable called modA.1
    assign(paste("modA",1,sep="."), randomForest(formula0, data = data, importance = TRUE))
    ##generate the variable A1star by using predict on dat.plasmode
    assign(paste("A",1,".star",sep=""),predict(get(paste("modA",1,sep=".")),newdata=dat.plasmode));
    
    ##include the generated variable in the simulated data dat.plasmode
    dat.plasmode[[eval(parse(text=paste0("A",1)))]] <- get(paste("A",1,".star",sep=""))
    dat.plasmode[[eval(parse(text=paste0("A",1)))]]<-factor(dat.plasmode[[eval(parse(text=paste0("A",1)))]])
    
    
    for (j in 2:timeobs){ ### j is the variable that indicates the time of observation
      
      
      ## getting the index of the elements of our covariates L1, L2, L3....
      for (u in 1:length(get(paste0("L",j)))){
        
        ## creating the dependent variables which are the elements of L1, L2, L3,...
        measurevar2 <- eval(parse(text=paste(paste0("L",j,"[",u,"]"),"",sep="")));
        
        ## creating the independent variables which are the treatments and the covariates from 1 to j-1 
        groupvars2<-c(mget(paste0("A",1:c(j-1))), unlist(mget(paste("L", 1:c(j-1), sep=""))));
        
        ## creating the formula of the model for the covariates( for example: L2~ A1+L1)
        formula2=as.formula(paste(measurevar2, paste(groupvars2, collapse=" + "), sep=" ~ "));
        
        ## assign the formula to a variable called modL.u
        assign(paste("modL",u,j,sep="."), randomForest(formula2, data = data, importance = TRUE))
        ## creating the new variable called Lj.star[u],using predict
        assign(paste0("L",j,".star","[",u,"]"),predict(get(paste("modL",u,j,sep=".")),newdata=dat.plasmode));
        
      }
      ### this is the dependent variable of treatment for each regression for each time
      measurevar <- mget(paste0("A",j)); 
      
      ### creating the independent variables(treatment at t-1 plus the covariates from t=1 to t=j)
      groupvars<-c(mget(paste0("A",1:c(j-1))),unlist(mget(paste("L", 1:j, sep=""))));
      
      ### creating the formula of regression, for example(A2~ A1+L1+L2)
      formula=as.formula(paste(measurevar, paste(groupvars, collapse=" + "), sep=" ~ "));
      
      
      ## assign the formula to a variable called modA.1
      assign(paste("modA",j,sep="."), randomForest(formula, data = data, importance = TRUE))
      ##generate the variable A1star by using predict on validation data
      assign(paste("A",j,".star",sep=""),predict(get(paste("modA",j,sep=".")),newdata=dat.plasmode));
      
      
      ## including the variable A.star in the simulated data
      dat.plasmode[[eval(parse(text=paste0("A",j)))]] <- get(paste("A",j,".star",sep=""))
      dat.plasmode[[eval(parse(text=paste0("A",j)))]] <-factor(dat.plasmode[[eval(parse(text=paste0("A",j)))]])
    }
    ### Simulation of the response variable, this is a function for a continuous outcome
    
    groupvars0<-c(mget(paste0("A",1:timeobs)),unlist(mget(paste("L", 1:timeobs, sep=""))));
    formulaY=as.formula(paste(Y, paste(groupvars0, collapse=" + "), sep=" ~ "));
    # model random Forest
    modelY <- randomForest(formulaY, 
                           data = data)
    
    ## creating the new variable
    Y.star=predict(modelY,newdata=dat.plasmode)
    ## inserting the new variable in the simulated data
    dat.plasmode[[Y]]=Y.star
    
    list_data_plasmode[[i]] <- dat.plasmode
    if(silent == FALSE) cat("plasmode dataset", i, "has been generated.\n"); 
  }
  ## return the list of simulated data
  return(list_data_plasmode)
  
  ## remove all the objects in the function besides the output
  rm(list=ls())
}

