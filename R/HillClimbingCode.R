#object of class "hillclimbing" contains the slots approach (the ascent method we select - SA, MA or LA), Budget, Number of Inputs, 
#Number of Connections, and the bounds of the different coefficients.
setClass(
  Class="hillclimbing", 
  representation=representation(    
    approach = "character",
    Budget="numeric", 
    Num_Inputs = "numeric", 
    Num_Connections = "numeric",
    Num_Landscapes ="numeric",
    Upperbound_inputvars = "numeric",
    Lowerbound_inputvars = "numeric",
    Upperbound_squares = "numeric",
    Lowerbound_squares = "numeric",
    Upperbound_crossprod = "numeric",
    Lowerbound_crossprod = "numeric"    
  )
)


# validity condition to ensure that the number of connections is less than the number of inputs, and that the approach name is correct 
setValidity("hillclimbing", function(object){
  if(object@Num_Connections>=object@Num_Inputs){
    return("number of connections must be lower than number of inputs")
  }
  else if(object@approach!="steepascent"&&object@approach!="medianascent"&&object@approach!="leastascent"){
    return("incorrect ascent method")
  }
  else{
    return(TRUE)}}
)


#class "climbingTools". This will contain arguments that will be used by 
#functions that are instrumental to the different hill climbing models, 
#such as the function to calculate profit or to find the neighborhoods.

setClass(
  Class="climbingTools", 
  representation=representation(
    Coeff_inputvars = "numeric",
    Coeff_squares = "numeric",
    Coeff_crossprod = "numeric",
    X1="numeric",
    Num_Connections="numeric"
  )
)

#function for standard error
#for methods with a name that doesn't correspond to any existing function
#we need to first define a 'generic' function before defining the method
setGeneric(
  name="standardError",
  def=function(object){standardGeneric("standardError")
  }
)
setMethod(
  f="standardError", 
  signature="climbingTools", 
  def=function(object){
  X1 <- object@X1
  return(sqrt(var(X1,na.rm=TRUE)/length(na.omit(X1))))
  }
)



#profit function
setGeneric(
  name="allocationProfit",
  def=function(object){standardGeneric("allocationProfit")
  }
)
setMethod(
  f="allocationProfit", 
  signature="climbingTools", 
  def=function(object){
    X1 <- object@X1
    Coeff_inputvars <- object@Coeff_inputvars
    Coeff_squares <- object@Coeff_squares
    Coeff_crossprod <- object@Coeff_crossprod
    X3=NULL		
    
    #mat1 is a special matrix which when applied to an allocation, gives the cross product of its elements
    mat1 <- multcomp:::contrMat(1:length(X1), "Tukey")
    
    #X3 contains the cross product
    X3=c(X3, apply(abs(mat1), 1, function(row){
      #product of row of mat1 and row of X1
      prod(X1[row == 1])}))
    
    #profit using the equation given in the paper
    profit <- sum(Coeff_inputvars*X1)+sum(Coeff_squares*((X1)^2))+sum((Coeff_crossprod)*X3)
    return(profit)
  }
)



#calculates the neighborhoods of an allocation by accepting an allocation, X, 
#and the number of connections, c, from the user. This will be used when we evaluate the various ascent methods
setGeneric(
  name="findNeighborhood",
  def=function(object){standardGeneric("findNeighborhood")
  }
)
setMethod(
  f="findNeighborhood", 
  signature="climbingTools", 
  def=function(object){    
    c <- object@Num_Connections
    X <- object@X1
    
    mat1 <- multcomp:::contrMat(1:(c+1), "Tukey")
    mat1 <- rbind(mat1,-1*mat1 )
    
    n.inputs <- length(X)
    n.groups <- n.inputs%/%(c+1)
    g1 <- 1:(c+1)
    groups <- sapply(0:(n.groups-1), function(i){
      g1 + (c+1)*i
    })
    
    index <- 1:ncol(groups)
    nbds <- lapply(1:ncol(groups), function(i){
      g <- groups[,i]
      temp.mat <- apply(mat1, 1, function(m){
        X.temp <- X
        X.temp[g] <- X.temp[g] + m
        return(X.temp)  
      })	
      return(temp.mat)
    })
    
    nbds <- t(do.call(cbind, nbds))
    neighborhoods <- NULL
    
    for(i in 1:nrow(nbds))
    {
      if(min(nbds[i,])>=0){
        neighborhoods <- rbind(neighborhoods, nbds[i,])	
      }		 
    }	
    
    #if there is no neighborhood, we return the argument X
    if(is.null(neighborhoods)){
      return(matrix(X, nrow=1))
    }
    
    else{
      return(neighborhoods)
    }	
  }  
)



#hillclimb() method that will climb the landscape depending on the approach selected by 
setGeneric(
  name="hillclimb",
  def=function(object){standardGeneric("hillclimb")
  }
)
setMethod(
  f="hillclimb", 
  signature="hillclimbing", 
  def=function(object){
    Num_Inputs <- object@Num_Inputs
    Budget <- object@Budget
    Num_Connections <- object@Num_Connections
    Upperbound_inputvars <- object@Upperbound_inputvars
    Lowerbound_inputvars <- object@Lowerbound_inputvars
    Upperbound_squares <- object@Upperbound_squares
    Lowerbound_squares <- object@Lowerbound_squares
    Upperbound_crossprod <- object@Upperbound_crossprod
    Lowerbound_crossprod <- object@Lowerbound_crossprod
    approach <- object@approach
    Num_Landscapes <- object@Num_Landscapes
    
    #landscape_Profit will contain local maxima of each landscape
    landscape_Profit = NULL
    
    #steps will store the number of steps it takes to reach the local maxima on a landscape
    steps <- NULL  
     
    #loop for 100 landscapes
    for(j in 1:Num_Landscapes){
      #coefficient vectors for the input variables, squares and cross product
      c1=c(runif(Num_Inputs, Lowerbound_inputvars, Upperbound_inputvars))
      c2=c(runif(Num_Inputs, Lowerbound_squares, Upperbound_squares))    
      c3=c(runif(choose(Num_Inputs, 2), Lowerbound_crossprod, Upperbound_crossprod))
      
      #X1 is an arbitrary allocation where we start the climb; randomly selected allocation
      X1 <- rmultinom(n=1, size=Budget, prob=rep(1/Num_Inputs, Num_Inputs))[,1]
      
      #neighborhoods of X1
      nbds <- findNeighborhood(new("climbingTools", X1=X1, Num_Connections = Num_Connections))
      
      #current_profit contains the profit from our current allocation
      current_Profit <- allocationProfit(new("climbingTools", X1=X1, Coeff_inputvars=c1, Coeff_squares=c2, Coeff_crossprod=c3))
      
      #Prof_nbds vector contains profit of each neighborhood
      Prof_nbds <- apply(nbds, 1, function(i){allocationProfit(new("climbingTools", X1=i, Coeff_inputvars=c1, Coeff_squares=c2, Coeff_crossprod=c3))})
      
      #subset of neighborhoods higher than the current profit level
      Prof2 <- Prof_nbds[Prof_nbds>current_Profit]
      
      #will store the number of steps for this landscape
      count <- 0
      
      #loop breaks when we there are no neighborhoods higher than our current location
      while(length(Prof2)!=0){ 
        #if statements that moves us to the next allocation depending on the ascent method we choose
        if(approach=="steepascent"){
          #m is the index of the allocation which we must move on to
          m <- which(Prof_nbds==max(Prof2))
        }
        else if(approach=="medianascent"){
          #medianProf2 is the index of the median of Prof2, whether its length is odd or even
          medianProf2 <- sapply(median(Prof2), function(x){which.min(abs(x - Prof2))})
          m <-which(Prof_nbds==Prof2[medianProf2])
        }
        else if(approach=="leastascent"){
          m <- which(Prof_nbds==min(Prof2)) 
        }    
        
        #move us to the allocation which is the lowest of the nbds higher than our current position
        X1 <- nbds[m,]
        
        current_Profit <- allocationProfit(new("climbingTools", X1=X1, Coeff_inputvars=c1, Coeff_squares=c2, Coeff_crossprod=c3)) 
        
        #neighborhoods of the new allocation we climb on to
        nbds <- findNeighborhood(new("climbingTools", X1=X1, Num_Connections = Num_Connections))
        
        #profits of the new neighborhoods we climb on to
        Prof_nbds <- apply(nbds, 1, function(i){allocationProfit(new("climbingTools", X1=i, Coeff_inputvars=c1, Coeff_squares=c2, Coeff_crossprod=c3))}) 	
        
        Prof2 <- Prof_nbds[Prof_nbds>current_Profit]
        
        count <- count+1 										 
      }
      #store the normalized profit by dividing the maxima by the budget^2+budget
      normalizedProfit <- (current_Profit)/(Budget+(Budget^2))*100
      
      steps <- c(steps, count)
      
      #landscape_Profit vector contains the local maxima of each landscape
      landscape_Profit=c(landscape_Profit, normalizedProfit)
    } 
    #returns a list with the standard error and mean of local maxima of 100 landscapes
    return(list(allProfits=landscape_Profit, StandardError = standardError(new("climbingTools", X1=landscape_Profit)), 
    MeanNormalizedProfit = mean(landscape_Profit), MeanSteps = mean(steps), 
    StandardErrorSteps = standardError(new("climbingTools", X1=steps))))
  }
) 


#method to return table with the mean normalized profits for all three methods for all connections from 1 to 5
setGeneric(
  name="hillclimbtable",
  def=function(object){standardGeneric("hillclimbtable")
  }
)
setMethod(
  f="hillclimbtable", 
  signature="hillclimbing", 
  def=function(object){
    #accessor methods to extract relevant data
    Budget <- object@Budget
    Num_Inputs <- object@Num_Inputs
    Upperbound_inputvars <- object@Upperbound_inputvars
    Lowerbound_inputvars <- object@Lowerbound_inputvars
    Upperbound_squares <- object@Upperbound_squares
    Lowerbound_squares <- object@Lowerbound_squares
    Upperbound_crossprod <- object@Upperbound_crossprod
    Lowerbound_crossprod <- object@Lowerbound_crossprod
    Num_Landscapes <- object@Num_Landscapes
    steepascent <- NULL
    medianascent <- NULL
    leastascent <- NULL
    for (i in 1:5){
      object <- new("hillclimbing", approach="steepascent", 
                    Budget=Budget, Num_Inputs=Num_Inputs, Num_Connections=i, Num_Landscapes=Num_Landscapes,
                    Upperbound_inputvars=Upperbound_inputvars, Lowerbound_inputvars=Lowerbound_inputvars, 
                    Upperbound_squares=Upperbound_squares, Lowerbound_squares=Lowerbound_squares, 
                    Upperbound_crossprod=Upperbound_crossprod, Lowerbound_crossprod=Lowerbound_crossprod)
      steepascent <- c(steepascent, hillclimb(object)$MeanNormalizedProfit)
    }
    for (i in 1:5){
      object <- new("hillclimbing", approach="medianascent", 
                    Budget=Budget, Num_Inputs=Num_Inputs, Num_Connections=i, Num_Landscapes=Num_Landscapes,
                    Upperbound_inputvars=Upperbound_inputvars, Lowerbound_inputvars=Lowerbound_inputvars, 
                    Upperbound_squares=Upperbound_squares, Lowerbound_squares=Lowerbound_squares, 
                    Upperbound_crossprod=Upperbound_crossprod, Lowerbound_crossprod=Lowerbound_crossprod)
      medianascent <- c(medianascent, hillclimb(object)$MeanNormalizedProfit)
    }
    for (i in 1:5){
      object <- new("hillclimbing", approach="leastascent", 
                    Budget=Budget, Num_Inputs=Num_Inputs, Num_Connections=i, Num_Landscapes=Num_Landscapes,
                    Upperbound_inputvars=Upperbound_inputvars, Lowerbound_inputvars=Lowerbound_inputvars, 
                    Upperbound_squares=Upperbound_squares, Lowerbound_squares=Lowerbound_squares, 
                    Upperbound_crossprod=Upperbound_crossprod, Lowerbound_crossprod=Lowerbound_crossprod)
      leastascent <- c(leastascent, hillclimb(object)$MeanNormalizedProfit)
    }
    approachName <- c("1", "2", "3", "4", "5")
    data <- rbind(approachName, steepascent, medianascent, leastascent)
    return(data)
  }
)



#method to boxplots with the mean normalized profits for all three methods and for all connections from 1 to 5
setGeneric(
  name="hillclimbboxplot",
  def=function(object){standardGeneric("hillclimbboxplot")
  }
)
setMethod(
  f="hillclimbboxplot", 
  signature="hillclimbing", 
  def=function(object){
    Budget <- object@Budget
    Num_Inputs <- object@Num_Inputs
    Upperbound_inputvars <- object@Upperbound_inputvars
    Lowerbound_inputvars <- object@Lowerbound_inputvars
    Upperbound_squares <- object@Upperbound_squares
    Lowerbound_squares <- object@Lowerbound_squares
    Upperbound_crossprod <- object@Upperbound_crossprod
    Lowerbound_crossprod <- object@Lowerbound_crossprod
    Num_Landscapes <- object@Num_Landscapes
    d1 <- NULL
    d2 <- NULL
    d3 <- NULL
    for (i in 1:5){
      object <- new("hillclimbing", approach="steepascent", 
                    Budget=Budget, Num_Inputs=Num_Inputs, Num_Connections=i, Num_Landscapes=Num_Landscapes,
                    Upperbound_inputvars=Upperbound_inputvars, Lowerbound_inputvars=Lowerbound_inputvars, 
                    Upperbound_squares=Upperbound_squares, Lowerbound_squares=Lowerbound_squares, 
                    Upperbound_crossprod=Upperbound_crossprod, Lowerbound_crossprod=Lowerbound_crossprod)
      d1 <- rbind(d1, data.frame(y=hillclimb(object)$allProfits, connections=factor(i), methodtype="steep"))
    }
    for (i in 1:5){
      object <- new("hillclimbing", approach="medianascent", 
                    Budget=Budget, Num_Inputs=Num_Inputs, Num_Connections=i, Num_Landscapes=Num_Landscapes,
                    Upperbound_inputvars=Upperbound_inputvars, Lowerbound_inputvars=Lowerbound_inputvars, 
                    Upperbound_squares=Upperbound_squares, Lowerbound_squares=Lowerbound_squares, 
                    Upperbound_crossprod=Upperbound_crossprod, Lowerbound_crossprod=Lowerbound_crossprod)
      d2 <- rbind(d2, data.frame(y=hillclimb(object)$allProfits, connections=factor(i), methodtype="median"))
    }
    for (i in 1:5){
      object <- new("hillclimbing", approach="leastascent", 
                    Budget=Budget, Num_Inputs=Num_Inputs, Num_Connections=i, Num_Landscapes=Num_Landscapes,
                    Upperbound_inputvars=Upperbound_inputvars, Lowerbound_inputvars=Lowerbound_inputvars, 
                    Upperbound_squares=Upperbound_squares, Lowerbound_squares=Lowerbound_squares, 
                    Upperbound_crossprod=Upperbound_crossprod, Lowerbound_crossprod=Lowerbound_crossprod)
      d3 <- rbind(d3, data.frame(y=hillclimb(object)$allProfits, connections=factor(i), methodtype="least"))
    }
    data.full <- rbind(d1, d2, d3)
    return(qplot(x=methodtype,y=y, data=data.full,geom="boxplot") + facet_grid(. ~connections, labeller=label_both))
  }
)



#twoD.simplex plots all combinations of three inputs that satisfy the budget constraint on a 2-D grid on a simplex, scaled to 100.
setGeneric(
  name="twoD.simplex",
  def=function(object){
    standardGeneric("twoD.simplex")
  }
)
setMethod(
  f="twoD.simplex", 
  signature="hillclimbing", 
  def=function(object){
    #findVectors() finds every possible allocation that meets the budget constraint for a given number of inputs
    V <- findVectors(3, object@Budget)
    X1=V[,1]
    X2=V[,2]
    X3=V[,3]
    X=data.frame(X1, X2, X3)
    ggtern(data=X,aes(X1, X2, X3), scale=1) + geom_point(fill="red", shape=21,size=2) + theme_bw()
  }
)



#method to produce randomly generated landscape
setGeneric(
  name="landscapePlot",
  def=function(object){standardGeneric("landscapePlot")
  }
)
setMethod(
  f="landscapePlot", 
  signature="hillclimbing", 
  def=function(object){
    Budget <- object@Budget
    Upperbound_inputvars <- object@Upperbound_inputvars
    Lowerbound_inputvars <- object@Lowerbound_inputvars
    Upperbound_squares <- object@Upperbound_squares
    Lowerbound_squares <- object@Lowerbound_squares
    Upperbound_crossprod <- object@Upperbound_crossprod
    Lowerbound_crossprod <- object@Lowerbound_crossprod
    Num_Inputs <- 3
    c1=c(runif(Num_Inputs, Lowerbound_inputvars, Upperbound_inputvars))
    c2=c(runif(Num_Inputs, Lowerbound_squares, Upperbound_squares))    
    c3=c(runif(choose(Num_Inputs, 2), Lowerbound_crossprod, Upperbound_crossprod))
    all_allocations <- findVectors(groups= Num_Inputs,size=Budget)
    profits <- apply(all_allocations, 1, function(i){
      allocationProfit(new("climbingTools", X1=i, 
        Coeff_inputvars=c1, Coeff_squares=c2, Coeff_crossprod=c3))})
    plot3d(x=all_allocations[,1],y=all_allocations[,2],z=profits, 
           col="red", xlab="Input 1", ylab="Input 2")
  }
)

