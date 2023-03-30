############################################Non parametric contrefactual plasmode###########################################

#### Arguments
# data:        A dataframe containing the real data to be used to generate simulated datasets (data.frame)
# n:           Sample size of the original data (numeric)
# timeobs:     Number of time-points of data (numeric)
# Y:           Name of the outcome variable (character)
# id:          Name of the id variable (character)
# A1, A2, ...: Names of the binary treatment variables (character)
# L1, L2, ...: Names of the covariates (vector of character)
# A.fixed:     Numeric matrix of dimension n.regimes x timeobs, where each row represents a treatment regime for which
#              counterfactual outcomes should be generated
# silent:      TRUE/FALSE to indicate whether or not the function should provide update on the number of datasets generated

#### Notes
# Categorical variables should be factors in data to be considered as such by the function

#### Output
# A list of n.regimes simulated datasets

       

library(randomForest);
plasmode_longitudinal_non_param1<- function(data, n, timeobs, Y, id, A1, A2, A3, L1, L2, L3,
                                            A.fixed, silent = FALSE, ...){
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


  list_data_plasmode_conterfactual<- list() ## ouput of the function
  ### create the data frame of the simulation
  data1=data
  for (w in 1:nrow(A.fixed)) { ## to select each conterfactual treatment (each row)
    
    for (j in 2:timeobs){ ### j is the variable that indicates the time of observation
      
      ## getting the index of the elements of the covariates L1, L2, L3....
      for (u in 1:length(get(paste0("L",j)))){ 
        
        ## creating the dependent variables which are the elements of L1, L2, L3,...
        measurevar2 <- eval(parse(text=paste(paste0("L",j,"[",u,"]"),"",sep="")));
        
        ## creating the independent variables which are the treatments and the covariates from 1 to j-1 
        groupvars2<-c(mget(paste0("A",1:c(j-1))), unlist(mget(paste("L", 1:c(j-1), sep=""))));
        
        ## creating the formula of the model for the covariates (for example: L2~ A1+L1)
        formula2=as.formula(paste(measurevar2, paste(groupvars2, collapse=" + "), sep=" ~ "));
        
        
        assign(paste("modL",u,j,sep="."), randomForest(formula2, data = data))
        
        # Prediction 
        newdat = data1; ## create the prediction data
        for (l in 1:timeobs) { 
          newdat[[paste0("A",l)]]<-rep(A.fixed[w,l],each=n);
          newdat[[paste0("A",l)]]<-factor(newdat[[eval(parse(text=paste0("A",l)))]], levels=levels(data[,eval(parse(text=paste0("A",l)))]))
        }
        
        assign(paste0("L",j,".star","[",u,"]"),predict(get(paste("modL",u,j,sep=".")),newdata=newdat,type="res"));
        
        newdat=data
        
        ## stock the new variable in the simulated data
        data1[[eval(parse(text=paste(paste0("L",j,"[",u,"]"),"",sep="")))]] <- get(paste0("L",j,".star","[",u,"]"))
      }
    }
    
    
    ### Simulation of the response variable 
    groupvars0<-c(mget(paste0("A",1:timeobs)),unlist(mget(paste("L", 1:timeobs, sep=""))));
    formulaY=as.formula(paste(Y, paste(groupvars0, collapse=" + "), sep=" ~ "));
    
    # model random Forest
    modelY <- randomForest(formulaY, 
                           data = data) 
    
    # 5. Prediction 
    newdat = data1;
    for (l in 1:timeobs) {
      newdat[[paste0("A",l)]]<- rep(A.fixed[w,l],each=n);
      newdat[[paste0("A",l)]]<-factor(newdat[[eval(parse(text=paste0("A",l)))]], levels=levels(data[,eval(parse(text=paste0("A",l)))]))
    }
    ## creating the new variable
    Y.star=predict(modelY,newdata=newdat)
    newdat=data
    
    ## inserting the new variable in the simulated data
    data1[[Y]]=Y.star
    
    ## change the values of the treatment at each time based on A.fixed 
    for (d in 1:timeobs) {
      eval( parse(text=paste("data1","[,",paste0("A",d),"] <- rep(A.fixed[w,d],each=n)",sep="")))
      data1[[paste0("A",d)]]<-factor(data1[[eval(parse(text=paste0("A",d)))]], levels=levels(data[,eval(parse(text=paste0("A",d)))]))
    }
    ## get the list of plasmode simulations
    list_data_plasmode_conterfactual[[w]] <- data1
    if(silent == FALSE) cat("plasmode dataset", w, "has been generated.\n"); 
  }
  ## return the list of simulated data
  return(list_data_plasmode_conterfactual)
}

