#### Arguments
# data:         A dataframe containing the real data to be used to generate simulated datasets (data.frame)
# n:            Sample size of the original data (numeric)
# size:         The number of observations to be generated (numeric)
# nsim:         The number of simulated datasets to be generated
# timeobs:      Number of time-points of data (numeric)
# Y:            Name of the outcome variable (character)
# id:           Name of the id variable (character)
# A1, A2, ...:  Names of the binary treatment variables (character)
# L1, L2, ...:  Names of the covariates (vector of character)
# distribution: Names of the distribution of each variable in data (character vector,
#               either "gaussian", "binomial" or "multinomial")
# silent:       TRUE/FALSE to indicate whether or not the function should provide update on the number of datasets generated

#### Notes
# If distribution is not provided, the function tries to guess it based on whether the variables
# are factors or not and their number of levels

#### Output
# A list of nsim simulated datasets


library(mlogit);
plasmode_longitudinal_param<- function(data, size, n, nsim=1, timeobs, Y, id, A1, A2, A3,
                                       L1, L2, L3, distribution = NULL, silent = FALSE,...){
  expit=plogis;

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

 
  # If the distribution of the variables was not provided, determine it
  if(is.null(distribution)){
    distribution = rep("gaussian", ncol(data));
    for(j in 1:ncol(data)){
      if(is.factor(data[,j]) == TRUE){
        if(nlevels(data[,j]) == 2) distribution[j] = "binomial";
        if(nlevels(data[,j]) > 2) distribution[j] = "multinomial";
      } 
    }
  }

  ids= data.frame(matrix(nrow = size, ncol = nsim))
  list_data_plasmode <- list()
  for(i in 1:nsim){
    idxs <- sample(n, size, replace = TRUE); ### Generate the id variable of the simulated data
    
    ids[,i] <- data[[id]][idxs]; 
    
    dat.plasmode= data[idxs,]; ### create the data frame of the simulation
    
    ## This is the simulation of A1 from the regression A1~L1
    groupvar=c(get(apply(expand.grid("L",1), 1, paste0, collapse=""))) ## this is L1: independent variable used in the regression model
    
    formula0=as.formula(paste(mget(paste0("A",1)),paste(groupvar, collapse=" + "), sep=" ~ "))  ## the regression formula A~L1
    
    ## assign the formula to a variable called modA.1
    assign(paste("modA",1,sep="."), glm(formula0, data = data, family = "binomial")); 
    
    ##here we get the coefficients of our model modA.1
    assign(paste("coefs",1,sep=""), coef(get(paste("modA",1,sep=".")))); 
    
    ## getting the model matrix of the model: XA1 at the first time of observation
    assign(paste("XA",1,sep=""), model.matrix(get(paste("modA",1,sep=".")),data = dat.plasmode));
    
    ##generate the variable A1star by using rbinom with probs=XA1%*%coefs
    assign(paste("A",1,".star",sep=""),rbinom(size, 1, expit(as.matrix(get(paste("XA",1,sep=""))%*%get(paste("coefs",1,sep=""))))));
    
    ##include the generated variable in the simulated data dat.plasmode
    dat.plasmode[[eval(parse(text=paste0("A",1)))]] <- get(paste("A",1,".star",sep=""))
    dat.plasmode[[eval(parse(text=paste0("A",1)))]] <-factor(dat.plasmode[[eval(parse(text=paste0("A",1)))]])
    
    for (j in 2:timeobs){ ### j is the variable that indicates the time of observation
      
      
      ## getting the index of the elements of our covariables L1, L2, L3....
      for (u in 1:length(get(paste0("L",j)))){
        
        ## creating the dependent variables which are the elements of L1, L2, L3,...
        measurevar2 <- eval(parse(text=paste(paste0("L",j,"[",u,"]"),"",sep="")));
        
        ## creating the independent variables which are the treatments and the covariates from 1 to j-1 
        groupvars2<-c(mget(paste0("A",1:c(j-1))), unlist(mget(paste("L", 1:c(j-1), sep=""))));
        
        
        ## indicate the type of the distribution
        if(distribution[which(colnames(dat.plasmode) == measurevar2)] == "gaussian"){
          ## creating the formula of the model for covariates (for example: L2~ A1+L1)
          formula2=as.formula(paste(measurevar2, paste(groupvars2, collapse=" + "), sep=" ~ "));
          
          ##assign the formula to modL.u.j
          assign(paste("modL",u,j,sep="."), lm(formula2, data = data)); 
          
          ##get the coefficients of the model
          assign(paste("coefsL",u,j,sep=""), coef(get(paste("modL",u,j,sep="."))));
          
          ## get the sigma of the model
          assign(paste("sigmaL",u,j,sep=""), sigma(eval(parse(text=paste("modL",u,j,sep=".")))));
          
          ## get the model matrix XLuj
          assign(paste("XL",u,j,sep=""), model.matrix(get(paste("modL",u,j,sep=".")),data = dat.plasmode));
          
          ## creating the new variable called Lj.star[u],using rnorm
          assign(paste0("L",j,".star","[",u,"]"),rnorm(size,
                                                       mean=get(paste("XL",u,j,sep=""))%*%get(paste("coefsL",u,j,sep="")),sd= get(paste("sigmaL",u,j,sep=""))));
          
          ## stocking the variable in the simulated data
          dat.plasmode[[eval(parse(text=paste(paste0("L",j,"[",u,"]"),"",sep="")))]] <- get(paste0("L",j,".star","[",u,"]"))
          
          ## in case of binomial distribution 
        } else if(distribution[which(colnames(dat.plasmode) == measurevar2)] == "binomial"){
          
          ## creating the formula of the model for covariates (for example: L2~ A1+L1)
          formula2=as.formula(paste(measurevar2, paste(groupvars2, collapse=" + "), sep=" ~ "));
          ## assign the formula to modL.u.j
          assign(paste("modL",u,j,sep="."), glm(formula2, data = data, family = "binomial"));
          ## get coefficients
          assign(paste("coefsL",u,j,sep=""), coef(get(paste("modL",u,j,sep="."))));
          
          ## get the model matrix XLuj
          assign(paste("XL",u,j,sep=""), model.matrix(get(paste("modL",u,j,sep=".")),data = dat.plasmode));
          
          ## create new variable Lj.star[u]
          assign(paste0("L",j,".star","[",u,"]"), rbinom(size, 1, expit(as.matrix(get(paste("XL",u,j,sep=""))%*%get(paste("coefsL",u,j,sep=""))))));
          ## stock the new variable in the simulated data
          dat.plasmode[[eval(parse(text=paste(paste0("L",j,"[",u,"]"),"",sep="")))]] <- get(paste0("L",j,".star","[",u,"]"))
       
          dat.plasmode[[eval(parse(text=paste(paste0("L",j,"[",u,"]"),"",sep="")))]] <- factor(dat.plasmode[[eval(parse(text=paste(paste0("L",j,"[",u,"]"),"",sep="")))]])
          
           }
        ## in the case of multinomial distribution
        if(distribution[which(colnames(dat.plasmode) == measurevar2)] == "multinomial"){
          ## create data
          ds = mlogit.data(data = data, shape = "wide", choice = measurevar2, varying = NULL, idvar = id);
          ## formula of mlogit
          formula3=as.formula(paste(measurevar2, paste(groupvars2, collapse=" + "), sep="~ 1 |"))
          ## assign formula to modL.u.j
          assign(paste("modL",u,j,sep="."), mlogit(formula3, data = ds));
          ## get coefficients
          assign(paste("coefsL",u,j,sep=""), matrix(coef(get(paste("modL",u,j,sep="."))),
                                                    nrow = nlevels(dat.plasmode[[measurevar2]]) - 1, byrow = FALSE));
          
          ## model matrix XLuj
          formula3b=as.formula(paste(measurevar2, paste(groupvars2, collapse=" + "), sep="~"))
          assign(paste("XL",u,j,sep=""), model.matrix(formula3b, data = dat.plasmode))
          
          Pa=as.matrix(rep(0,size))
          
          ## to generalize for all possible levels
          for (w in 1:(nlevels(dat.plasmode[[measurevar2]])-1)) {
            Paa<-Pa+exp(get(paste("XL",u,j,sep=""))%*%get(paste("coefsL",u,j,sep=""))[w,])
          }
          
          ## get the first probability 
          P1 = (1 + Paa)^(-1);
          
          ## get the other probabilities 
          for(k in 2:nlevels(dat.plasmode[[measurevar2]])){
            assign(paste0("P", k), exp(get(paste("XL",u,j,sep=""))%*%get(paste("coefsL",u,j,sep=""))[k-1,])*P1)
          } 
          ## create new variable Lj.star[u]
          assign(paste0("L",j,".star","[",u,"]"), t(apply(data.frame(mget(paste("P", 1:nlevels(dat.plasmode[[measurevar2]]), sep=""))), 1, rmultinom, n = 1, size = 1))%*%c(1:nlevels(dat.plasmode[[measurevar2]])));
          
          ## stock the new variable in the simulated data
          dat.plasmode[[measurevar2]] <- get(paste0("L",j,".star","[",u,"]"))
          ### i modified here !
          dat.plasmode[[measurevar2]]<-factor(dat.plasmode[[measurevar2]])
        }}
      
      ### this is the dependent variable of treatment for each regression for each time
      measurevar <- mget(paste0("A",j)); 
      
      ### creating the independent variables(treatment at t-1 plus the covariables from t=1 to t=j)
      groupvars<-c(mget(paste0("A",1:c(j-1))),unlist(mget(paste("L", 1:j, sep=""))));
      
      ### creating the regression formula, for example(A2~ A1+L1+L2)
      formula=as.formula(paste(measurevar, paste(groupvars, collapse=" + "), sep=" ~ "));
      
      ## assign the formula to an output called modA.j depending on the time of observation
      assign(paste("modA",j,sep="."), glm(formula, data = data, family = "binomial"));
      
      ## get the coefficients of the model under the name of coefs,j: j is the time of observation
      assign(paste("coefs",j,sep="") ,coef(get(paste("modA",j,sep="."))));
      
      ## getting the model matrix of the model: XAj at each time of observation
      assign(paste("XA",j,sep=""),model.matrix(get(paste("modA",j,sep=".")),data = dat.plasmode));
      
      ## update the variable of treatment at each t by generating a binomial variable called A.star
      assign(paste("A",j,".star",sep=""),rbinom(size, 1, expit(as.matrix(get(paste("XA",j,sep=""))%*%get(paste("coefs",j,sep=""))))));
      
      ## including the variable A.star in the simulated data
      dat.plasmode[[eval(parse(text=paste0("A",j)))]] <- get(paste("A",j,".star",sep=""))
      dat.plasmode[[eval(parse(text=paste0("A",j)))]]<-factor(dat.plasmode[[eval(parse(text=paste0("A",j)))]])}
    
  
    groupvars0<-c(mget(paste0("A",1:timeobs)),unlist(mget(paste("L", 1:timeobs, sep=""))));
    formulaY=as.formula(paste(Y, paste(groupvars0, collapse=" + "), sep=" ~ "));
    ## indicate the type of the distribution
    if(distribution[which(colnames(dat.plasmode) == Y)] == "gaussian"){
      ## formula of Y in function of the treatments and covariables at all time (the time of observation is indicated by the user in the arguments of the function)
      ## creating the model
      modY=lm(formulaY, data = data)
      ## coefs of the model
      coefsY=coef(modY)
      ## sigma of the model
      sigmaY=sigma(modY)
      ##model matrix
      Xy=model.matrix(modY,data=dat.plasmode)
      ## creating the new variable
      Y.star=rnorm(size,mean=Xy%*%coefsY,sd=sigmaY)
      ## inserting the new variable in the simulated data
      dat.plasmode[[Y]]=Y.star}
    
 
    if(distribution[which(colnames(dat.plasmode) == Y)] == "binomial"){
      ## model 
      modY=glm(formulaY, data = data,family="binomial")
      ## coefs and sigma
      coefsY=coef(modY) 
      sigmaY=sigma(modY)
      ## model matrix
      Xy=model.matrix(modY,data=dat.plasmode)
      ## generate the new variable by rbinom
      Y.star=rbinom(size,1,expit(as.matrix(Xy%*%coefsY)))
      dat.plasmode[[Y]]=Y.star
      dat.plasmode[[Y]]=factor(dat.plasmode[[Y]])}
    
    
    if(distribution[which(colnames(dat.plasmode) == Y)] == "multinomial"){
      ## create data
      ds = mlogit.data(data = data, shape = "wide", choice = Y, varying = NULL, idvar = id);
      ## model
      formula2=as.formula(paste(Y, paste(groupvars0, collapse=" + "), sep="~ 1 |"))
      modY=mlogit(formula2,data=ds)
      ## get coefficients
      coefsY=matrix(coef(modY),nrow = nlevels(dat.plasmode[[Y]]) - 1, byrow = FALSE) 
      ## model matrix XLuj
      formula2b=as.formula(paste(Y, paste(groupvars0, collapse=" + "), sep="~"))
      Xy=model.matrix(formula2b, data = dat.plasmode)
      Pa=as.matrix(rep(0,size))
      ## to generalize for all possible levels
      for (w in 1:(nlevels(dat.plasmode[[Y]])-1)) {
        Pa=Pa+exp(Xy%*%coefsY[w,])
      }
      
      ## get the first probability 
      P1 = (1 + Pa)^(-1);
      
      ## get the other probabilities 
      for(k in 2:nlevels(dat.plasmode[[Y]])){
        assign(paste0("P", k), exp(Xy%*%coefsY[k-1,])*P1)
      } 
      ## generate the new variable by rmultinom
      Y.star=t(apply(data.frame(mget(paste("P", 1:nlevels(dat.plasmode[[Y]]), sep=""))), 1, rmultinom, n = 1, size = 1))%*%c(1:nlevels(dat.plasmode[[Y]]))
      dat.plasmode[[Y]]=Y.star
      dat.plasmode[[Y]]=factor(dat.plasmode[[Y]])}
    
    list_data_plasmode[[i]] <- dat.plasmode
    if(silent == FALSE) cat("plasmode dataset", i, "has been generated.\n"); 
  }
  ## return the list of simulated data
  return(list_data_plasmode)
  
  ## remove all the objects in the function besides the output
  rm(list=ls())
}


