as.ggsList <- function(jmfit){
  if(inherits(jmfit, "JMbayes")|inherits(jmfit, "mvJMbayes")){
    lapply(names(jmfit$mcmc), function(memberName){
      member = jmfit$mcmc[[memberName]]
      if(!is.null(member)){
        if(memberName=="inv_D" | memberName=="D"){
          nsubmember = dim(member)[3]
          subMembers = vector("list", nsubmember)
          for(i in 1:nsubmember){
            subMembers[[i]] = ggs(as.mcmc(as.matrix(member[,,i])))
          }
          return(subMembers)
        }else if(memberName=="b"){
          nsubmember = dim(member)[1]
          subMembers = vector("list", nsubmember)
          for(i in 1:nsubmember){
            subMembers[[i]] = ggs(as.mcmc(t(member[i,,])))
          }
          return(subMembers)
        }else{
          return(ggs(as.mcmc(as.matrix(member))))
        }
      }else{
        return(NULL)
      }
    })
  }else{
    stop("The object you passed is not a JMbayes/mvJMbayes object")
  }
}