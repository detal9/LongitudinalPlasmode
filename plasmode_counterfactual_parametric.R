#### Arguments
# data:         A dataframe containing the real data to be used to generate simulated datasets (data.frame)
# n:            Sample size of the original data (numeric)
# timeobs:      Number of time-points of data (numeric)
# Y:            Name of the outcome variable (character)
# id:           Name of the id variable (character)
# A1, A2, ...:  Names of the binary treatment variables (character)
# L1, L2, ...:  Names of the covariates (vector of character)
# A.fixed:      Numeric matrix of dimension n.regimes x timeobs, where each row represents a treatment regime for which
#               counterfactual outcomes should be generated
# distribution: Names of the distribution of each variable in data (character vector,
#               either "gaussian", "binomial" or "multinomial")
# silent:       TRUE/FALSE to indicate whether or not the function should provide update on the number of datasets generated

#### Notes
# If distribution is not provided, the function tries to guess it based on whether the variables
# are factors or not and their number of levels

#### Output
# A list of n.regimes simulated datasets


plasmode_longitudinal_param1<- function(data, n, timeobs, Y, id, A1, A2, A3, 
                                        L1, L2, L3, A.fixed, distribution = NULL, silent = FALSE,...){
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
  
  ### If the distribution of the variables was not provided, determine it
  if(is.null(distribution)){
    distribution = rep("gaussian", ncol(data));
    for(j in 1:ncol(data)){
      if(is.factor(data[,j]) == TRUE){
        if(nlevels(data[,j]) == 2) distribution[j] = "binomial";
        if(nlevels(data[,j]) > 2) distribution[j] = "multinomial";
      } 
    }
  }

  list_data_plasmode_conterfactual <- list()
  ### create the data frame of the simulation
  data1=data
  for (w in 1:nrow(A.fixed)) { ## to select each conterfactual treatment (each row)
    for (j in 2:timeobs){ ### j is the variable that indicates the time of observation
      ## getting the index of the elements of our covariates L1, L2, L3....
      for (u in 1:length(get(paste0("L",j)))){
        
        ## creating the dependent variables which are the elements of L1, L2, L3,...
        measurevar2 <- eval(parse(text=paste(paste0("L",j,"[",u,"]"),"",sep="")));
        
        ## creating the independent variables which are the treatments and the covariates from 1 to j-1 
        groupvars2<-c(mget(paste0("A",1:c(j-1))), unlist(mget(paste("L", 1:c(j-1), sep=""))));
        
        ## creating the formula of the model of covariables( for example: L2~ A1+L1)
        formula2=as.formula(paste(measurevar2, paste(groupvars2, collapse=" + "), sep=" ~ "));
        
        ## indicate the type of the distribution
        if(distribution[which(colnames(data) == measurevar2)] == "gaussian"){
          ##assign the formula to modL.u.j
          assign(paste("modL",u,j,sep="."), lm(formula2, data = data)); 
          
          ##get the coefficients of the model
          assign(paste("coefsL",u,j,sep=""), coef(get(paste("modL",u,j,sep="."))));
          
          ## get the sigma of the model
          assign(paste("sigmaL",u,j,sep=""), sigma(eval(parse(text=paste("modL",u,j,sep=".")))));
          
          ## get the model matrix XLuj
          assign(paste("XL",u,j,sep=""), model.matrix(get(paste("modL",u,j,sep="."))));
          
          ## set the treatment to the A.fixed matrix chosen by the user (a loop for each treatment)
          for (l in 1:c(j-1)) {
            eval(parse(text=paste(paste("XL",u,j,sep=""),"[,",grep(paste0("A",l),colnames(get(paste("XL",u,j,sep="")))),"]<- rep(A.fixed[w,l],each=n)",sep="")));
          }
          
          ## creating the new variable called Lj.star[u],using rnorm
          assign(paste0("L",j,".star","[",u,"]"),rnorm(n,
                                                       mean=get(paste("XL",u,j,sep=""))%*%get(paste("coefsL",u,j,sep="")),sd= get(paste("sigmaL",u,j,sep=""))));
          
          ## stocking the variable in the simulated data
          data1[[eval(parse(text=paste(paste0("L",j,"[",u,"]"),"",sep="")))]] <- get(paste0("L",j,".star","[",u,"]"))
          
          
          ## in case of binomial distribution 
        } else if(distribution[which(colnames(data) == measurevar2)] == "binomial"){
          ## assign the formula to modL.u.j
          assign(paste("modL",u,j,sep="."), glm(formula2, data = data, family = "binomial"));
          ## get coefficients
          assign(paste("coefsL",u,j,sep="") ,coef(get(paste("modL",u,j,sep="."))));
          
          ## get the model matrix XLuj
          assign(paste("XL",u,j,sep="") ,model.matrix(get(paste("modL",u,j,sep="."))));
          
          ## set A to A.fixed
          for (l in 1:c(j-1)) {
            eval(parse(text=paste(paste("XL",u,j,sep=""),"[,",grep(paste0("A",l),colnames(get(paste("XL",u,j,sep="")))),"]<- rep(A.fixed[w,l],each=n)",sep="")));
          }
          ## create new variable Lj.star[u]
          assign(paste0("L",j,".star","[",u,"]"),rbinom(n, 1, expit(as.matrix(get(paste("XL",u,j,sep=""))%*%get(paste("coefsL",u,j,sep=""))))));
          ## stock the new variable in the simulated data
          data1[[eval(parse(text=paste(paste0("L",j,"[",u,"]"),"",sep="")))]] <- get(paste0("L",j,".star","[",u,"]"))
        }
        ## in the case of multinomial distribution
        if(distribution[which(colnames(data) == measurevar2)] == "multinomial"){
          ## create data
          ds = mlogit.data(data = data, shape = "wide", choice = measurevar2, varying = NULL, idvar = id);
          ## formula of mlogit
          formula3=as.formula(paste(measurevar2, paste(groupvars2, collapse=" + "), sep="~ 1 |"))
          ## assign formula to modL.u.j
          assign(paste("modL",u,j,sep="."), mlogit(formula3, data = ds));
          ## get coefficients
          assign(paste("coefsL",u,j,sep=""), matrix(coef(get(paste("modL",u,j,sep="."))),
                                                    nrow = nlevels(data[[measurevar2]]) - 1, byrow = FALSE));
          
      
          ## model matrix XLuj
          formula3b=as.formula(paste(measurevar2, paste(groupvars2, collapse=" + "), sep="~"))
          assign(paste("XL",u,j,sep=""), model.matrix(formula3b, data = data))
          
          ## set the treatment to the A.fixed matrix chosen by the user 
          for (l in 1:c(j-1)) {
            eval(parse(text=paste(paste("XL",u,j,sep=""),"[,",grep(paste0("A",l),colnames(get(paste("XL",u,j,sep="")))),"]<- rep(A.fixed[w,l],each=n)",sep="")));
          }
          Pa=as.matrix(rep(0,n))
          
          ## to generalize for all possible levels
          for (z in 1:(nlevels(data[[measurevar2]])-1)) {
            Pa=Pa+exp(get(paste("XL",u,j,sep=""))%*%get(paste("coefsL",u,j,sep=""))[z,])
          }
          
          
          ## get the first probability 
          P1 = (1 + Pa)^(-1);
          
          ## get the other probabilities 
          for(k in 2:nlevels(data[[eval(parse(text=paste(paste0("L",j,"[",u,"]"),"",sep="")))]])){
            assign(paste0("P", k), exp(get(paste("XL",u,j,sep=""))%*%get(paste("coefsL",u,j,sep=""))[k-1,])*P1)
          } 
          ## create new variable Lj.star[u]
          assign(paste0("L",j,".star","[",u,"]"), t(apply(data.frame(mget(paste("P", 1:nlevels(data[[eval(parse(text=paste(paste0("L",j,"[",u,"]"),"",sep="")))]]), sep=""))), 1, rmultinom, n = 1, size = 1))%*%c(1:nlevels(data[[eval(parse(text=paste(paste0("L",j,"[",u,"]"),"",sep="")))]])));
          
          ## stock the new variable in the simulated data
          data1[[eval(parse(text=paste(paste0("L",j,"[",u,"]"),"",sep="")))]] <- get(paste0("L",j,".star","[",u,"]"))}
          data1[[eval(parse(text=paste(paste0("L",j,"[",u,"]"),"",sep="")))]]<-factor(data1[[eval(parse(text=paste(paste0("L",j,"[",u,"]"),"",sep="")))]])
        }}
    
    ### Simulation of the response variable, this is a function for continuous outcome!
    groupvars0<-c(mget(paste0("A",1:timeobs)), unlist(mget(paste("L", 1:timeobs, sep=""))));
    formulaY=as.formula(paste(Y, paste(groupvars0, collapse=" + "), sep=" ~ "));
    ## indicate the type of the distribution
    if(distribution[which(colnames(data) == Y)] == "gaussian"){
      ## formula of Y in function of the treatments and covariables at all time (the time of observation is indicated by the user in the arguments of the function)
      ## creating the model
      modY=lm(formulaY, data = data)
      ## coefs of the model
      coefsY=coef(modY)
      ## sigma of the model
      sigmaY=sigma(modY)
      ##model matrix
      Xy=model.matrix(modY)
      
      ## set the treatment to the A.fixed matrix chosen by the user 
      Xy[,grep("A",colnames(Xy))]<-matrix(rep(A.fixed[w,1:timeobs],each=n),ncol=timeobs)
      ## creating the new variable
      Y.star=rnorm(n,mean=Xy%*%coefsY,sd=sigmaY)
      ## inserting the new variable in the simulated data
      data1[[Y]]=Y.star}
    
    
    if(distribution[which(colnames(data) == Y)] == "binomial"){
      ## model 
      modY=glm(formulaY, data = data,family="binomial")
      ## coefs and sigma
      coefsY=coef(modY) 
      sigmaY=sigma(modY)
      ## model matrix
      Xy=model.matrix(modY)
      ## set the treatment to the A.fixed matrix chosen by the user 
      Xy[,grep("A",colnames(Xy))]<-matrix(rep(A.fixed[w,1:timeobs],each=n), ncol=timeobs)
      ## generate the new variable by rbinom
      Y.star=rbinom(size,1,expit(as.matrix(Xy%*%coefsY)))
      data1[[Y]]=Y.star}
    
    
    if(distribution[which(colnames(data) == Y)] == "multinomial"){
      ## create data
      ds = mlogit.data(data = data, shape = "wide", choice = Y, varying = NULL, idvar = id);
      ## model
      formula2=as.formula(paste(Y, paste(groupvars0, collapse=" + "), sep="~ 1 |"))
      modY=mlogit(formula2,data=ds)
      ## get coefficients
      coefsY=matrix(coef(modY),nrow = nlevels(dat.plasmode[[Y]]) - 1, byrow = FALSE) 
      ## model matrix XLuj
      formula2b=as.formula(paste(Y, paste(groupvars0, collapse=" + "), sep="~"))
      Xy=model.matrix(formula2b, data = data)
      ## set the treatment to the A.fixed matrix choosed by the user 
      Xy[,grep("A",colnames(Xy))]<-matrix(rep(A.fixed[w,1:timeobs],each=n),ncol=timeobs)
      
      Pa=as.matrix(rep(0,n))
      ## to generalize for all possible levels
      for (w in 1:(nlevels(data[[Y]])-1)) {
        Pa=Pa+exp(Xy%*%coefsY[w,])
      }
      
      ## get the first probability 
      P1 = (1 + Pa)^(-1);
      
      ## get the other probabilities
      for(k in 2:nlevels(data[[Y]])){
        assign(paste0("P", k), exp(Xy%*%coefsY[k-1,])*P1)
      } 
      ## generate the new variable by rmultinom
      Y.star=t(apply(data.frame(mget(paste("P", 1:nlevels(data[[Y]]), sep=""))), 1, rmultinom, n = 1, size = 1))%*%c(1:nlevels(data[[Y]]))
      data1[[Y]]=Y.star
      data1[[Y]]=factor(data1[[Y]])}
    
    ## change the values of the treatment at each time based on A.fixed 
    for (d in 1:timeobs) {
      eval( parse(text=paste("data1","[,",paste0("A",d),"] <- rep(A.fixed[w,d],each=n)",sep="")))
    }
    ## get the list of plasmode simulations
    list_data_plasmode_conterfactual[[w]] <- data1
    if(silent == FALSE) cat("plasmode dataset", w, "has been generated.\n"); 
  }
  ## return the list of simulated data
  return(list_data_plasmode_conterfactual)
}
