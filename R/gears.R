#' Personal greeting
#'
#' @description Greet a person and appropriately capitalize their name.
#'
#' @param name Your name (character string; e.g. "john doe").
#'
#' @return A character string, capitalized to title case.
#' @export
#'
#' @examples
#' hello("james bond")
gears<-function()
{
  # revised MHStep
  basefn<-"50nodes20and10edges200sampUneqVarUpdate1to10"
  #basefn<-"20nodes20and10edges100sampNopenal1to10"

  starts<-1
  ends<-10

  library(pscl)
  dataAll<-read.table("simuData10nodes20edges200samp.txt",header=F)
  gamAll<-read.table("gams10nodes20edges200samp.txt",header=F)
  betaAll<-read.table("Beta10nodes20edges200samp.txt",header=F)


  nodes<-ncol(dataAll)
  numData<-100
  sampsize<-nrow(dataAll)/numData

  possParents<-seq(0,(nodes-1))

  # sample gamCom, gamX, or gamY (the gam in the function)
  # beta is the coefficients in the regressions
  # sigma2 is the pre-specified large variances in the mixed distribution of beta_ij
  # sig2 is for the variances in the regression residuals.
  # ComXY is the data, either Z or X or Y

  sampGam<-function(gam,beta,sigma2,sig2,n,data)
  {
    #pp is the probability to be included
    pp<-0.1
    sig22<-rep(sig2,nodes)
    for (kk in 2:nodes)
    {
      loc<-sum(possParents[1:(kk-1)])+1
      for (j in 1:(kk-1))
      {
        gamTmp<-gam
        betaTmp<-beta[loc:(loc+possParents[kk]-1)]*gam[loc:(loc+possParents[kk]-1)]
        a<--(beta[(loc+j-1)])^2/(2*sigma2)-n/2*log(sig22[kk])-sum((data[,kk]-as.matrix(data[,1:(kk-1)])%*%betaTmp)^2)/(2*sig22[kk])
        gamTmp[(loc+j-1)]<-1-gam[(loc+j-1)]
        betaTmp<-beta[loc:(loc+possParents[kk]-1)]*gamTmp[loc:(loc+possParents[kk]-1)]
        b<- -n/2*log(sig22[kk])-sum((data[,kk]-as.matrix(data[,1:(kk-1)])%*%betaTmp)^2)/(2*sig22[kk])
        r<-1/(1+exp(b-a)*(1-pp)/pp)*gam[(loc+j-1)]+1/(1+(1-pp)/pp*exp(a-b+(beta[(loc+j-1)])^2/(2*sigma2)))*(1-gam[(loc+j-1)])
        #			r<-1/(1+exp(b-a))*gam[(loc+j-1)]+1/(1+exp(a-b))*(1-gam[(loc+j-1)])

        judge<-runif(1)
        gam[(loc+j-1)]<-sum(r>=judge)
      }
    }
    return(gam)
  }

  sampGamCom<-function(gam,beta,sigma2,sigX2,sigY2,nx,n,data)
  {
    #pp is the probability to be included
    pp<-0.1
    sigX22<-rep(sigX2,nodes)
    sigY22<-rep(sigY2,nodes)
    for (kk in 2:nodes)
    {
      loc<-sum(possParents[1:(kk-1)])+1
      for (j in 1:(kk-1))
      {
        gamTmp<-gam
        betaTmp<-beta[loc:(loc+possParents[kk]-1)]*gam[loc:(loc+possParents[kk]-1)]
        a1<--(beta[(loc+j-1)])^2/(2*sigma2)-nx/2*log(sigX22[kk])-sum((data[1:nx,kk]-as.matrix(data[1:nx,1:(kk-1)])%*%betaTmp)^2)/(2*sigX22[kk])
        a2<--(beta[(loc+j-1)])^2/(2*sigma2)-(n-nx)/2*log(sigY22[kk])-sum((data[(nx+1):n,kk]-as.matrix(data[(nx+1):n,1:(kk-1)])%*%betaTmp)^2)/(2*sigY22[kk])
        a<-a1+a2

        gamTmp[(loc+j-1)]<-1-gam[(loc+j-1)]
        betaTmp<-beta[loc:(loc+possParents[kk]-1)]*gamTmp[loc:(loc+possParents[kk]-1)]
        b1<- -nx/2*log(sigX22[kk])-sum((data[1:nx,kk]-as.matrix(data[1:nx,1:(kk-1)])%*%betaTmp)^2)/(2*sigX22[kk])
        b2<- -(n-nx)/2*log(sigY22[kk])-sum((data[(nx+1):n,kk]-as.matrix(data[(nx+1):n,1:(kk-1)])%*%betaTmp)^2)/(2*sigY22[kk])
        b<-b1+b2

        r<-1/(1+exp(b-a)*(1-pp)/pp)*gam[(loc+j-1)]+1/(1+(1-pp)/pp*exp(a-b+(beta[(loc+j-1)])^2/(2*sigma2)))*(1-gam[(loc+j-1)])
        #			r<-1/(1+exp(b-a))*gam[(loc+j-1)]+1/(1+exp(a-b))*(1-gam[(loc+j-1)])

        judge<-runif(1)
        gam[(loc+j-1)]<-sum(r>=judge)
      }
    }
    return(gam)
  }


  # Calculate probabilities of a given graph. probGraph(gamX,betaX,sigma2,sig2,nx,dataX)

  probGraph<-function(gam,beta,sigma2,sig2,n,data)
  {
    prob<-0
    sig22<-rep(sig2,nodes)
    for (kk in 2:nodes)
    {
      loc<-sum(possParents[1:(kk-1)])+1
      gamTmp<-gam
      betaTmp<-beta[loc:(loc+possParents[kk]-1)]*gam[loc:(loc+possParents[kk]-1)]
      prob<-prob-n/2*log(sig22[kk])-sum((data[,kk]-as.matrix(data[,1:(kk-1)])%*%betaTmp)^2)/(2*sig22[kk])
    }
    return(prob)
  }

  probGraphCom<-function(gam,beta,sigma2,sigX2,sigY2,n,nx,data)
  {
    prob<-0
    sigX22<-rep(sigX2,nodes)
    sigY22<-rep(sigY2,nodes)

    for (kk in 2:nodes)
    {
      loc<-sum(possParents[1:(kk-1)])+1
      gamTmp<-gam
      betaTmp<-beta[loc:(loc+possParents[kk]-1)]*gam[loc:(loc+possParents[kk]-1)]
      part1<- -nx/2*log(sigX22[kk])-sum((data[1:nx,kk]-as.matrix(data[1:nx,1:(kk-1)])%*%betaTmp)^2)/(2*sigX22[kk])
      part2<- -(n-nx)/2*log(sigY22[kk])-sum((data[(nx+1):n,kk]-as.matrix(data[(nx+1):n,1:(kk-1)])%*%betaTmp)^2)/(2*sigY22[kk])
      prob<-prob+part1+part2
    }
    return(prob)
  }

  # sample eta -- the indicator variable for common or differential networks
  sampEta<-function(data,eta,betaComTmp,betaXTmp,betaYTmp,sigX2,sigY2,nx,gamCom,gamX,gamY)
  {
    n<-nrow(data)
    ny<-n-nx
    #	sigC2<-(sigX2+sigY2)/2

    # the common network
    #	sumCom<--n/2*sum(log(sigC2))
    sumCom<-0
    betaCom<-betaComTmp*gamCom
    for (i in 2:nodes)
    {
      loc<-sum(possParents[1:(i-1)])+1
      #		betaCom<-beta[loc:(loc+possParents[i]-1)]*gam[loc:(loc+possParents[i]-1)]
      part1<- -sum((data[1:nx,i]-as.matrix(data[1:nx,1:(i-1)])%*%betaCom[loc:(loc+possParents[i]-1)])^2)/(2*sigX2[i])
      part2<- -sum((data[(nx+1):n,i]-as.matrix(data[(nx+1):n,1:(i-1)])%*%betaCom[loc:(loc+possParents[i]-1)])^2)/(2*sigY2[i])
      sumCom<-sumCom+part1+part2
    }
    sumCom<- sumCom-nx/2*sum(log(sigX2))-(n-nx)/2*sum(log(sigY2))

    # the X network
    dataX<-data[1:nx,]
    sumX<--nx/2*sum(log(sigX2))
    #	sumX<-0
    betaX<-betaXTmp*gamX
    for (i in 2:nodes)
    {
      loc<-sum(possParents[1:(i-1)])+1
      sumX<- sumX-sum((dataX[,i]-as.matrix(dataX[,1:(i-1)])%*%betaX[loc:(loc+possParents[i]-1)])^2)/(2*sigX2[i])
    }

    # the Y network
    dataY<-data[(nx+1):n,]
    sumY<--ny/2*sum(log(sigY2))
    #	sumY<-0
    betaY<-betaYTmp*gamY
    for (i in 2:nodes)
    {
      loc<-sum(possParents[1:(i-1)])+1
      sumY<- sumY-sum((dataY[,i]-as.matrix(dataY[,1:(i-1)])%*%betaY[loc:(loc+possParents[i]-1)])^2)/(2*sigY2[i])
    }
    disagrY<-sum(abs(gamY-gamCom))/length(gamCom)
    disagrX<-sum(abs(gamX-gamCom))/length(gamCom)
    disagrXY<-sum(abs(gamX-gamY))/length(gamCom)
    #	cat("disagreeY ",disagrY," disagreeX ",disagrX," disagreeXY ",disagrXY,"\n")
    #	if ((disagrY<=0.05 && disagrX<=0.05)||disagrXY<=0.05)
    #	{
    #		eta<-1
    #	}
    #	else
    #	{
    #		r<-exp(sumCom-sumX-sumY)
    sumGamCom<-sum(gamCom)
    sumGamX<-sum(gamX)
    sumGamY<-sum(gamY)
    disXY<-length(gamX)-sum((gamX&gamY)|((1-gamX)&(1-gamY)))
    #		cat("disXY ",disXY," sumGamCom ",sumGamCom," sumGamX ",sumGamX," sumGamY ",sumGamY,"\n")
    dis<-0
    #		dis<-log(n)*log(disXY/sumGamCom)
    #		dis<-n*log(disXY/sumGamCom)

    # 	      r<-1+exp(sumX+sumY-sumCom-sumGamX*log(nx)-sumGamY*log(ny)+sumGamCom*log(n))
    #The following is what has been used
    # 	      r<-1+exp(sumX+sumY-sumCom-sumGamX*log(nx)-sumGamY*log(ny)+sumGamCom*log(n)+dis)
    # The following was the one being checked, and now it is the final one. It takes dis=0, and use half of the original (sumGamX*log(nx)-sumGamY*log(ny)+sumGamCom*log(n))
    r<-1+exp(sumX+sumY-sumCom-0.5*sumGamX*log(nx)-0.5*sumGamY*log(ny)+0.5*sumGamCom*log(n)+dis)
    # Taking out the penalty. This is to check the need of penalty
    # 	      r<-1+exp(sumX+sumY-sumCom)

    # 	      r<-1+exp(sumX+sumY-sumCom)
    ## 	      r<-1+exp(sumX+sumY-sumCom-mean(log(sumGamX)*log(nx),log(sumGamY)*log(ny))+log(sumGamCom)*log(n))
    # 	      r<-1+exp(sumX+sumY-sumCom-log(sumGamX)*log(nx)-log(sumGamY)*log(ny)+log(sumGamCom)*log(n))
    # 	      r<-1+exp(sumX+sumY-sumCom-sumGamX/nodes*log(nx)-sumGamY/nodes*log(ny)+sumGamCom/nodes*log(n))
    # 	      r<-1+exp(sumX+sumY-sumCom-sumGamX/sqrt(nodes)*log(nx)-sumGamY/sqrt(nodes)*log(ny)+sumGamCom/sqrt(nodes)*log(n))
    # 	      r<-1+exp(sumX+sumY-sumCom-sumGamX/sqrt(nodes)*log(nx)-sumGamY/sqrt(nodes)*log(ny)+sumGamCom/sqrt(nodes)*log(n)+log(dis))

    r<-r^(-1)
    eta<-sum(round(r,digits=1)>=0.5)
    #		cat("dis ",dis," sumX+sumY ",sumX+sumY," sumCom ",sumCom," sumXY-sumCom ",sumX+sumY-sumCom,"\n")
    #		cat("r ",r," penalty ",0.5*(-sumGamX*log(nx)-sumGamY*log(ny)+sumGamCom*log(n))+dis,"\n")
    #		judge<-runif(1)
    #		eta<-sum(r>judge)
    #      	cat("r ",r," eta ",eta,"\n")
    #	}
    return(eta)
  }


  # sample beta -- the coefficients in common or differential networks. i is for node, k is for the coefficient to be updated

  sampBeta<-function(beta,data, i,k,gam,sig2)
  {
    sig22<-rep(sig2,nodes)
    loc<-sum(possParents[1:(i-1)])+1
    gamTmp<-gam
    gamTmp[loc+k-1]<-0
    betaTmp<-beta*gamTmp
    betaTmp<-betaTmp[loc:(loc+possParents[i]-1)]
    betaTmp<-betaTmp[-k]
    dataTmp<-data[,1:(i-1)]
    if (i>2)
    {
      dataTmp<-dataTmp[,-k]
      CZ<-t(data[,i]-as.matrix(dataTmp)%*%betaTmp)%*%as.matrix(data[,k])
    }
    else
    {
      CZ<-t(data[,i])%*%as.matrix(data[,k])
    }
    denom<-sum(data[,k]^2)*sigma2+sig22[i]
    meanBeta<-sigma2*CZ/denom
    stdBeta<-sqrt(sig22[i]*sigma2/denom)
    temp<-rnorm(1,meanBeta,stdBeta)
    beta[loc+k-1]<-temp
    return(beta)
  }

  sampBetaCom<-function(beta,data,nx,i,k,gam,sigX2,sigY2)
  {
    sigX22<-rep(sigX2,nodes)
    sigY22<-rep(sigY2,nodes)

    loc<-sum(possParents[1:(i-1)])+1
    gamTmp<-gam
    gamTmp[loc+k-1]<-0
    betaTmp<-beta*gamTmp
    betaTmp<-betaTmp[loc:(loc+possParents[i]-1)]
    betaTmp<-betaTmp[-k]
    dataTmp<-data[,1:(i-1)]
    dataTmpX<-data[1:nx,1:(i-1)]
    dataTmpY<-data[(nx+1):n,1:(i-1)]


    if (i>2)
    {
      dataTmpX<-dataTmpX[,-k]
      dataTmpY<-dataTmpY[,-k]

      CZ1<-sigY2*t(data[1:nx,i]-as.matrix(dataTmpX)%*%betaTmp)%*%as.matrix(data[1:nx,k])
      CZ2<-sigX2*t(data[(nx+1):n,i]-as.matrix(dataTmpY)%*%betaTmp)%*%as.matrix(data[(nx+1):n,k])
      CZ=CZ1+CZ2
    }
    else
    {
      CZ<-sigY2*t(data[1:nx,i])%*%as.matrix(data[1:nx,k])+sigX2*t(data[(nx+1):n,i])%*%as.matrix(data[(nx+1):n,k])
    }
    denom<-sum(data[1:nx,k]^2)*sigma2*sigY2+sum(data[(nx+1):n,k]^2)*sigma2*sigX2+sigX2*sigY2
    meanBeta<-sigma2*CZ/denom
    stdBeta<-sqrt(sigX2*sigY2*sigma2/denom)
    temp<-rnorm(1,meanBeta,stdBeta)
    beta[loc+k-1]<-temp
    return(beta)
  }


  # sample sigC2, sigX2, and sigY2. Note at this stage these sigma-squares are assumed common for
  # all the nodes. This may have to be revised to assume node-specific variances.
  # i is for node

  #sampSig<-function(beta,gam, data,i)
  #{
  #	loc<-sum(possParents[1:(i-1)])+1
  #	shape<-nrow(data)/2+a0
  #	betaTmp<-beta*gam
  #	scale<-sum((data[,i]-as.matrix(data[,1:(i-1)])%*%betaTmp[loc:(loc+possParents[i]-1)])^2)/2+b0
  #	return(rigamma(1,shape, scale))
  #}

  MHStep<-function(order,k,iii,nodes,etaSampled,betaSampledX,betaSampledUniqueX,gamSampledX,gamSampledUniqueX,
                   betaSampledY,betaSampledUniqueY,gamSampledY,gamSampledUniqueY,
                   betaSampledCom,betaSampledUniqueCom,gamSampledCom,gamSampledUniqueCom,
                   sigX2SampledUnique,sigY2SampledUnique,nx,data,probTmp)
  {
    n<-nrow(data)
    #	cat("iii ",iii,"\n")
    #	cat("order in ",order[k,(iii-1),],"\n")
    if (iii==((K-k)*(B+N)+1))
    {
      orderTmp<-order[k,1,]
      gamX<-gamSampledX[k,1,2:ncol(gamSampledX[k,,])]
      betaX<-betaSampledX[k,1,2:ncol(betaSampledX[k,,])]

      gamY<-gamSampledY[k,1,2:ncol(gamSampledY[k,,])]
      betaY<-betaSampledY[k,1,2:ncol(betaSampledY[k,,])]

      betaCom<-betaSampledCom[k,1,2:ncol(betaSampledCom[k,,])]

      sigX2<-sigX2SampledUnique[1]
      sigY2<-sigY2SampledUnique[1]
    }
    else
    {
      orderTmp<-order[k,(iii-1),]
      gamX<-gamSampledX[k,(iii-1),2:ncol(gamSampledX[k,,])]
      betaX<-betaSampledX[k,(iii-1),2:ncol(betaSampledX[k,,])]

      gamY<-gamSampledY[k,(iii-1),2:ncol(gamSampledY[k,,])]
      betaY<-betaSampledY[k,(iii-1),2:ncol(betaSampledY[k,,])]

      betaCom<-betaSampledCom[k,(iii-1),2:ncol(betaSampledCom[k,,])]

      sigX2<-sigX2SampledUnique[order[k,(iii-1),1]]
      sigY2<-sigY2SampledUnique[order[k,(iii-1),1]]
    }
    #	cat("length(orderTmp) ",length(orderTmp)," order ",orderTmp[2:length(orderTmp)],"\n")
    dataTmp<-data[,orderTmp[2:length(orderTmp)]]

    probO<-OrderProb(gamX,betaX,gamY,betaY,betaCom,sigma2,sigX2,sigY2,n,dataTmp,nx)

    probOld<-probO$prob
    #	cat("probOld ",probOld,"\n")
    gamXOld<-probO$avgGamX
    betaXOld<-probO$avgBetaX

    gamYOld<-probO$avgGamY
    betaYOld<-probO$avgBetaY

    gamComOld<-probO$avgGamCom
    betaComOld<-probO$avgBetaCom

    sigX2Old<-probO$sigX2
    sigY2Old<-probO$sigY2
    etaOld<-probO$eta

    #	loc<-sort(sample(seq(1,nodes), floor(nodes/2), replace = FALSE, prob = NULL))
    loc<-sort(sample(seq(1,nodes), 2, replace = FALSE, prob = NULL))

    #	cat("loc ",loc,"\n")
    orderTmp1<-orderTmp
    tmp<-orderTmp[(loc[1]+1)]
    for (mm in 2:length(loc))
    {
      orderTmp1[(loc[mm-1]+1)]<-orderTmp[(loc[mm]+1)]
    }
    orderTmp1[(loc[mm]+1)]<-tmp
    if (iii==((K-k)*(B+N)+20))
    {
      orderTmp1[2:length(orderTmp1)]<-seq(1,nodes)
    }
    #	cat("length(orderTmp1) ",length(orderTmp1)," orderTmp1 ",orderTmp1,"\n")
    dataTmp<-data[,orderTmp1[2:length(orderTmp1)]]

    for (kk in 1:K)
    {
      if (kk==1)
      {
        uniqueOrder<-unique(order[kk,1:(iii-1),1])
      }
      else
      {
        uniqueOrder<-unique(c(uniqueOrder,unique(order[kk,1:(iii-1),1])))
      }
    }
    NAloc<-which(is.na(uniqueOrder)==1)
    if (length(NAloc!=0))	uniqueOrder<-uniqueOrder[-NAloc]

    for (kk in 1:length(uniqueOrder))
    {
      check<-0
      for (t in 1:K)
      {
        loc<-max(which(order[t,,1]==uniqueOrder[kk]))
        if (loc>0)
        {
          diff<-orderTmp1[2:length(orderTmp1)]-order[t,loc,2:(nodes+1)]
          if ((sum(abs(diff)))==0)
          {
            check<-1
            locPick<-loc
            kChose<-t
            break
          }
        }
        if (check==1) break
      }
      if (check==1) break
    }# end for (kk in 1:length(uniqueOrder))
    if (check==1)
    {
      betaX<-betaSampledUniqueX[order[kChose,loc,1],]
      gamX<-gamSampledUniqueX[order[kChose,loc,1],]

      betaY<-betaSampledUniqueY[order[kChose,loc,1],]
      gamY<-gamSampledUniqueY[order[kChose,loc,1],]

      betaCom<-betaSampledUniqueCom[order[kChose,loc,1],]
      gamCom<-gamSampledUniqueCom[order[kChose,loc,1],]

      sigX2<-sigX2SampledUnique[order[kChose,loc,1]]
      sigY2<-sigY2SampledUnique[order[kChose,loc,1]]
    }
    else
    {
      # X network
      gamX<-NULL
      for (i in 2:nodes)
      {
        if (i==2)
        {
          gamX<-rbinom((i-1),1,prob=0.5)
        }
        else
        {
          gamX<-c(gamX,rbinom((i-1),1,prob=0.5))
        }
      } # end for (i in 2:nodes)
      betaX<-rep(1,length(gamX))

      # Y network
      gamY<-NULL
      for (i in 2:nodes)
      {
        if (i==2)
        {
          gamY<-rbinom((i-1),1,prob=0.5)
        }
        else
        {
          gamY<-c(gamY,rbinom((i-1),1,prob=0.5))
        }
      } # end for (i in 2:nodes)
      betaY<-rep(1,length(gamY))

      # identical network
      gamCom<-as.numeric(gamX|gamY)
      betaCom<-rep(1,length(gamCom))


      sigX2<-sigY2<-1
    } # end else
    probO<-OrderProb(gamX,betaX,gamY,betaY,betaCom,sigma2,sigX2,sigY2,n,dataTmp,nx)

    #	cat("probO$prob ",probO$prob," probOld ",probOld,"\n")

    lpprobk<- -max(-probO$prob,H[k])/T[k]
    lpprobOldk<- -max(-probOld,H[k])/T[k]
    ratio<-min(1,exp(lpprobk-lpprobOldk))
    #	cat("ratio ",ratio,"\n")
    if (ratio>runif(1))
    {
      # assign the new order to the iteration iii.
      order[k,iii,]<-orderTmp1
      if (check==1)
      {
        order[k,iii,1]<-order[kChose,locPick,1]
        betaSampledUniqueX[order[kChose,locPick,1],]<-probO$avgBetaX
        gamSampledUniqueX[order[kChose,locPick,1],]<-probO$avgGamX

        betaSampledUniqueY[order[kChose,locPick,1],]<-probO$avgBetaY
        gamSampledUniqueY[order[kChose,locPick,1],]<-probO$avgGamY

        betaSampledUniqueCom[order[kChose,locPick,1],]<-probO$avgBetaCom
        gamSampledUniqueCom[order[kChose,locPick,1],]<-probO$avgGamCom

        sigX2SampledUnique[order[kChose,locPick,1]]<-probO$sigX2
        sigY2SampledUnique[order[kChose,locPick,1]]<-probO$sigY2


        betaSampledX[k,iii,1]<-order[k,iii,1]
        betaSampledX[k,iii,2:ncol(betaSampledX[k,,])]<-probO$avgBetaX
        gamSampledX[k,iii,1]<-order[k,iii,1]
        gamSampledX[k,iii,2:ncol(gamSampledX[k,,])]<-probO$avgGamX

        betaSampledY[k,iii,1]<-order[k,iii,1]
        betaSampledY[k,iii,2:ncol(betaSampledY[k,,])]<-probO$avgBetaY
        gamSampledY[k,iii,1]<-order[k,iii,1]
        gamSampledY[k,iii,2:ncol(gamSampledY[k,,])]<-probO$avgGamY

        betaSampledCom[k,iii,1]<-order[k,iii,1]
        betaSampledCom[k,iii,2:ncol(betaSampledCom[k,,])]<-probO$avgBetaCom
        gamSampledCom[k,iii,1]<-order[k,iii,1]
        gamSampledCom[k,iii,2:ncol(gamSampledCom[k,,])]<-probO$avgGamCom

        etaSampled[k,iii,1]<-order[k,iii,1]
        etaSampled[k,iii,2]<-probO$eta

        probTmp<-probO$prob
      }
      else
      {
        order[k,iii,1]<-max(uniqueOrder)+1
        betaSampledUniqueX<-rbind(betaSampledUniqueX,probO$avgBetaX)
        gamSampledUniqueX<-rbind(gamSampledUniqueX,probO$avgGamX)

        betaSampledUniqueY<-rbind(betaSampledUniqueY,probO$avgBetaY)
        gamSampledUniqueY<-rbind(gamSampledUniqueY,probO$avgGamY)

        betaSampledUniqueCom<-rbind(betaSampledUniqueCom,probO$avgBetaCom)
        gamSampledUniqueCom<-rbind(gamSampledUniqueCom,probO$avgGamCom)

        sigX2SampledUnique<-c(sigX2SampledUnique,probO$sigX2)
        sigY2SampledUnique<-c(sigY2SampledUnique,probO$sigY2)

        betaSampledX[k,iii,1]<-order[k,iii,1]
        betaSampledX[k,iii,2:ncol(betaSampledX[k,,])]<-probO$avgBetaX
        gamSampledX[k,iii,1]<-order[k,iii,1]
        gamSampledX[k,iii,2:ncol(gamSampledX[k,,])]<-probO$avgGamX

        betaSampledY[k,iii,1]<-order[k,iii,1]
        betaSampledY[k,iii,2:ncol(betaSampledY[k,,])]<-probO$avgBetaY
        gamSampledY[k,iii,1]<-order[k,iii,1]
        gamSampledY[k,iii,2:ncol(gamSampledY[k,,])]<-probO$avgGamY

        betaSampledCom[k,iii,1]<-order[k,iii,1]
        betaSampledCom[k,iii,2:ncol(betaSampledCom[k,,])]<-probO$avgBetaCom
        gamSampledCom[k,iii,1]<-order[k,iii,1]
        gamSampledCom[k,iii,2:ncol(gamSampledCom[k,,])]<-probO$avgGamCom

        etaSampled[k,iii,1]<-order[k,iii,1]
        etaSampled[k,iii,2]<-probO$eta

        probTmp<-probO$prob
      }
    } # end if (ratio>runif(1))
    else
    {
      if (iii==((K-k)*(B+N)+1))
      {
        order[k,iii,]<-order[k,1,]
      }
      else
      {
        order[k,iii,]<-order[k,(iii-1),]
      }
      betaSampledX[k,iii,1]<-order[k,iii,1]
      #		cat("length(betaXOld) ",length(betaXOld)," length(2:ncol(betaSampledX[k,,])) ",
      #				length(2:ncol(betaSampledX[k,,])),"\n")
      betaSampledX[k,iii,2:ncol(betaSampledX[k,,])]<-betaXOld
      gamSampledX[k,iii,1]<-order[k,iii,1]
      gamSampledX[k,iii,2:ncol(gamSampledX[k,,])]<-gamXOld


      betaSampledY[k,iii,1]<-order[k,iii,1]
      #		cat("length(betaYOld) ",length(betaYOld)," length(2:ncol(betaSampledY[k,,])) ",
      #				length(2:ncol(betaSampledY[k,,])),"\n")
      betaSampledY[k,iii,2:ncol(betaSampledY[k,,])]<-betaYOld
      gamSampledY[k,iii,1]<-order[k,iii,1]
      gamSampledY[k,iii,2:ncol(gamSampledY[k,,])]<-gamYOld

      betaSampledCom[k,iii,1]<-order[k,iii,1]
      betaSampledCom[k,iii,2:ncol(betaSampledCom[k,,])]<-betaComOld
      gamSampledCom[k,iii,1]<-order[k,iii,1]
      gamSampledCom[k,iii,2:ncol(gamSampledCom[k,,])]<-gamComOld

      etaSampled[k,iii,1]<-order[k,iii,1]
      etaSampled[k,iii,2]<-etaOld

      probTmp<-probOld
    }
    #	cat("iii ",iii," k ",k," eta ",etaSampled[k,iii,2]," order out ",order[k,iii,],"\n")
    return(list(order=order,betaSampledX=betaSampledX,betaSampledUniqueX=betaSampledUniqueX,gamSampledX=gamSampledX,
                gamSampledUniqueX=gamSampledUniqueX,
                betaSampledY=betaSampledY,betaSampledUniqueY=betaSampledUniqueY,gamSampledY=gamSampledY,
                gamSampledUniqueY=gamSampledUniqueY,
                betaSampledCom=betaSampledCom,betaSampledUniqueCom=betaSampledUniqueCom,gamSampledCom=gamSampledCom,
                gamSampledUniqueCom=gamSampledUniqueCom,
                etaSampled=etaSampled,
                sigX2SampledUnique=sigX2SampledUnique,sigY2SampledUnique=sigY2SampledUnique,probTmp=probTmp))
  }



  OrderProb<-function(gamX,betaX,gamY,betaY,betaCom,sigma2,sigX2,sigY2,n,data,nx)
  {
    gamXSum<-0
    sumBetaX<-0

    gamYSum<-0
    sumBetaY<-0

    gamComSum<-0
    sumBetaCom<-0

    eta<-0

    dataX<-data[1:nx,]
    dataY<-data[(nx+1):n,]
    ny<-n-nx

    for (ii in 1:loops)
    {
      gamX<-sampGam(gamX,betaX,sigma2,sigX2,nx,dataX)
      for (i in 2:nodes)
      {
        loc<-sum(possParents[1:(i-1)])+1
        nonZero<-which(gamX[loc:(loc+possParents[i]-1)]!=0)
        for (k1 in nonZero)
        {
          betaX<-sampBeta(betaX,dataX, i,k1,gamX,sigX2)
        }
      }

      gamY<-sampGam(gamY,betaY,sigma2,sigY2,nx,dataY)
      for (i in 2:nodes)
      {
        loc<-sum(possParents[1:(i-1)])+1
        nonZero<-which(gamY[loc:(loc+possParents[i]-1)]!=0)
        for (k1 in nonZero)
        {
          betaY<-sampBeta(betaY,dataY, i,k1,gamY,sigY2)
        }
      }

      gamCom<-as.numeric(gamX|gamY)

      sigC2<-(sigY2+sigX2)/2

      gamCom<-sampGamCom(gamCom,betaCom,sigma2,sigX2,sigY2,nx,n,data)

      for (i in 2:nodes)
      {
        loc<-sum(possParents[1:(i-1)])+1
        nonZero<-which(gamCom[loc:(loc+possParents[i]-1)]!=0)
        for (k1 in nonZero)
        {
          betaCom<-sampBetaCom(betaCom,data,nx,i,k1,gamCom,sigX2,sigY2)
        }
      }

      shape1<-nrow(dataX)/2+a0
      scale1<-sum((dataX[,1])^2)/2+b0
      sigX2<-rigamma(1,shape1,scale1)

      shape1<-nrow(dataY)/2+a0
      scale1<-sum((dataY[,1])^2)/2+b0
      sigY2<-rigamma(1,shape1,scale1)

      #		sigC2<-((n-nx-1)*sigY2+(nx-1)*sigX2)/(n-2)

      if (ii>(loops/2))
      {
        gamXSum<-gamXSum+gamX
        sumBetaX<-sumBetaX+betaX

        gamYSum<-gamYSum+gamY
        sumBetaY<-sumBetaY+betaY

        gamComSum<-gamComSum+gamCom
        sumBetaCom<-sumBetaCom+betaCom
      }

    } # end number of loops
    avgGamX<-round(gamXSum/(loops/2))
    avgBetaX<-sumBetaX/(loops/2)*avgGamX

    avgGamY<-round(gamYSum/(loops/2))
    avgBetaY<-sumBetaY/(loops/2)*avgGamY

    avgGamCom<-round(gamComSum/(loops/2))
    avgBetaCom<-sumBetaCom/(loops/2)*avgGamCom

    #	eta<-1

    sigX2Tmp<-rep(sigX2,nodes)
    sigY2Tmp<-rep(sigY2,nodes)
    eta<-sampEta(data,eta,avgBetaCom,avgBetaX,avgBetaY,sigX2Tmp,sigY2Tmp,nx,avgGamCom,avgGamX,avgGamY)


    #	eta<-0
    if (eta<=0.5)
    {
      probX<-probGraph(avgGamX,avgBetaX,sigma2,sigX2,nx,dataX)
      probY<-probGraph(avgGamY,avgBetaY,sigma2,sigY2,ny,dataY)
      prob<-probX+probY
    }
    else
    {
      prob<-probGraphCom(avgGamCom,avgBetaCom,sigma2,sigX2,sigY2,n,nx,data)
    }
    #	cat("avgGamX ",avgGamX,"\n")
    #	cat("avgBetaX ",avgBetaX,"\n")
    return(list(eta=eta,avgGamX=avgGamX,avgBetaX=avgBetaX,avgGamY=avgGamY,avgBetaY=avgBetaY,
                avgGamCom=avgGamCom,avgBetaCom=avgBetaCom,sigX2=sigX2,sigY2=sigY2,sigC2=sigC2,prob=prob))
  }






  # number of MCMC iterations
  loops<-50
  K<-5
  # N is the number of iterations to be used to collect results
  B<-10
  N<-10
  orderloops<-(B+N)*K+B+N



  agreement<-rep(0,numData)
  truePos<-rep(0,numData)
  truePosBeta<-rep(0,numData)
  trueNeg<-rep(0,numData)
  falsePos<-rep(0,numData)

  agreementX<-rep(0,numData)
  truePosX<-rep(0,numData)
  truePosXBeta<-rep(0,numData)
  trueNegX<-rep(0,numData)
  falsePosX<-rep(0,numData)

  agreementY<-rep(0,numData)
  truePosY<-rep(0,numData)
  truePosYBeta<-rep(0,numData)
  trueNegY<-rep(0,numData)
  falsePosY<-rep(0,numData)

  choice<-rep(0,numData)
  jj1<-1
  jj2<-1

  jj<-1
  pee<-0.1

  accurGX<-accurGY<-accurGCom<-rep(0,numData)
  fpGX<-fpGY<-fpGCom<-rep(0,numData)
  tpGX<-tpGY<-tpGCom<-rep(0,numData)
  tpBX<-tpBY<-tpBCom<-rep(0,numData)

  avgEtaEst<-rep(0,numData)

  #for (jj in 1:numData)
  for (jj in starts:ends)
  {
    cat("data ",jj,"\n")
    #   set.seed(2000*jj)
    # use a different seed to check influence from sampling error
    set.seed(12345*jj)

    data<-dataAll[(1+(jj-1)*sampsize):(jj*sampsize),]
    n<-nrow(data)
    nx<-nrow(data)/2

    #  this is for separated networks
    gamTrueX<-gamAll[(2*(jj-1)+1),]
    gamTrueY<-gamAll[(2*(jj-1)+2),]

    trueBetaX<-betaAll[(2*(jj-1)+1),]
    trueBetaY<-betaAll[(2*(jj-1)+2),]

    #  this will be used for combined networks
    #   gamTrueX<-gamAll[jj,]
    #   gamTrueY<-gamAll[jj,]
    gamTrueCom<-gamAll[jj,]

    #   trueBetaX<-betaAll[jj,]
    #   trueBetaY<-betaAll[jj,]
    trueBetaCom<-betaAll[jj,]

    #  Initial values

    # pre-specified large variances in the mixed distribution of beta_ij
    sigma2<-100

    #initial values
    #   eta<-1
    a0<-0.001
    b0<-0.001

    gamCom<-NULL
    for (i in 2:nodes)
    {
      if (i==2)
      {
        gamCom<-rbinom((i-1),1,prob=0.5)
      }
      else
      {
        gamCom<-c(gamCom,rbinom((i-1),1,prob=0.5))
      }
    }

    # X network
    gamX<-NULL
    for (i in 2:nodes)
    {
      if (i==2)
      {
        gamX<-rbinom((i-1),1,prob=0.5)
      }
      else
      {
        gamX<-c(gamX,rbinom((i-1),1,prob=0.5))
      }
    }

    # Y network
    gamY<-NULL
    for (i in 2:nodes)
    {
      if (i==2)
      {
        gamY<-rbinom((i-1),1,prob=0.5)
      }
      else
      {
        gamY<-c(gamY,rbinom((i-1),1,prob=0.5))
      }
    }

    # initinal values for beta
    betaCom<-rep(1,length(gamCom))
    betaX<-rep(1,length(gamX))
    betaY<-rep(1,length(gamCom))
    betaJoint<-pmax(betaCom,betaX,betaY)

    # initial values for the variances in the regression residuals
    sigC2<-sigX2<-sigY2<-1

    gamComSum<-rep(0,length(gamCom))
    gamXSum<-rep(0,length(gamX))
    gamYSum<-rep(0,length(gamY))
    sumBetaCom<-sumBetaX<-sumBetaY<-0

    betaXTruePos<-betaYTruePos<-betaComTruePos<-0

    #  set H0, H1, ..., HK

    # D is the rings containing samples of orders. In total, K rings.
    # the first element is for sample index, the second is for order of nodes (nodes), the order index (1), and the chain index (1),
    # and the third is for ring index

    D<-array(NA,dim=c((K*orderloops),(nodes+2),K))
    Dindex<-rep(0,K)

    # order is for the ordering collected in each chain from each iteration. The first element is the chain index, the second is for
    # the iteration index, and the third is the ordering.
    order<-array(NA,dim=c(K, K*orderloops,(nodes+1)))

    # a beta matrix corresponding to each order with the first column an index for order
    # this index should corresponding to the order sampled.
    OrderIndex<-1
    betaSampledX<-array(NA,dim=c(K, K*orderloops,(nodes*(nodes-1)/2+1)))
    betaSampledUniqueX<-matrix(betaX,nrow=1)
    gamSampledX<-array(NA,dim=c(K, K*orderloops,(nodes*(nodes-1)/2+1)))
    gamSampledUniqueX<-matrix(gamX,nrow=1)

    betaSampledY<-array(NA,dim=c(K, K*orderloops,(nodes*(nodes-1)/2+1)))
    betaSampledUniqueY<-matrix(betaY,nrow=1)
    gamSampledY<-array(NA,dim=c(K, K*orderloops,(nodes*(nodes-1)/2+1)))
    gamSampledUniqueY<-matrix(gamY,nrow=1)

    betaSampledCom<-array(NA,dim=c(K, K*orderloops,(nodes*(nodes-1)/2+1)))
    betaSampledUniqueCom<-matrix(betaCom,nrow=1)
    gamSampledCom<-array(NA,dim=c(K, K*orderloops,(nodes*(nodes-1)/2+1)))
    gamSampledUniqueCom<-matrix(gamCom,nrow=1)

    etaSampled<-array(NA,dim=c(K, K*orderloops,2))

    #   dataX<-data[1:nx,]
    # probably only use the following to get an idea of the choice of H0
    #   TestprobOrder<-probGraph(gamX,betaX,sigma2,sigX2,nx,dataX)
    # calculated -probOrder divided by 8000
    H0<-10000/8000
    #   H0<-10000/200

    #   H5<--TestprobOrder
    geo<-2
    H1<-H0*geo
    H2<-H1*geo
    H3<-H2*geo
    H4<-H3*geo
    H<-c(H0,H1,H2,H3,H4)

    #  set T0 to TK
    T0<-1
    T1<-T0*geo
    T2<-T1*geo
    T3<-T2*geo
    T4<-T3*geo
    T<-c(T0,T1,T2,T3,T4)


    for (iii in 1:K)
    {
      order[iii,1,]<-c(OrderIndex,seq(1,nodes))
      #	order[iii,1,]<-c(OrderIndex,sample(seq(1,nodes),nodes))

      betaSampledX[iii,1,]<-c(OrderIndex,betaX)
      gamSampledX[iii,1,]<-c(OrderIndex,gamX)

      betaSampledY[iii,1,]<-c(OrderIndex,betaY)
      gamSampledY[iii,1,]<-c(OrderIndex,gamY)

      betaSampledCom[iii,1,]<-c(OrderIndex,betaCom)
      gamSampledCom[iii,1,]<-c(OrderIndex,gamCom)

      etaSampled[iii,1,]<-c(OrderIndex,1)
    }
    #   orderTmp<-order[1,1,]
    #   cat("start order ",orderTmp,"\n")
    #   dataXtmp<-dataX[,orderTmp[2:length(orderTmp)]]

    #   probO<-OrderProb(gamX,betaX,sigma2,sigX2,nx,dataXtmp)
    #   probTmp<-probO$prob
    probTmp<-0

    sigX2SampledUnique<-sigX2
    sigY2SampledUnique<-sigY2

    for (iii in 2:orderloops)
    {
      for (k in K:1)
      {
        #	cat("iii ",iii," k ",k,"\n")

        if (iii>((K-k)*(B+N)))
        {
          if (k==K)
          {
            out<-MHStep(order,k,iii,nodes,etaSampled,betaSampledX,betaSampledUniqueX,gamSampledX,gamSampledUniqueX,
                        betaSampledY,betaSampledUniqueY,gamSampledY,gamSampledUniqueY,
                        betaSampledCom,betaSampledUniqueCom,gamSampledCom,gamSampledUniqueCom,
                        sigX2SampledUnique,sigY2SampledUnique,nx,data,probTmp)
            #			cat("k=K ",k,"\n")
            order<-out$order
            betaSampledX<-out$betaSampledX
            betaSampledUniqueX<-out$betaSampledUniqueX
            gamSampledX<-out$gamSampledX
            gamSampledUniqueX<-out$gamSampledUniqueX

            betaSampledY<-out$betaSampledY
            betaSampledUniqueY<-out$betaSampledUniqueY
            gamSampledY<-out$gamSampledY
            gamSampledUniqueY<-out$gamSampledUniqueY

            betaSampledCom<-out$betaSampledCom
            betaSampledUniqueCom<-out$betaSampledUniqueCom
            gamSampledCom<-out$gamSampledCom
            gamSampledUniqueCom<-out$gamSampledUniqueCom

            sigX2SampledUnique<-out$sigX2SampledUnique
            sigY2SampledUnique<-out$sigY2SampledUnique
            probTmp<-out$probTmp

            etaSampled<-out$etaSampled
          } # end if (k==K)
          else
          {
            judge<-runif(1)
            if (pee<judge)
            {
              #				cat("pee ",pee, " runif(1) ",judge,"\n")
              out<-MHStep(order,k,iii,nodes,etaSampled,betaSampledX,betaSampledUniqueX,gamSampledX,gamSampledUniqueX,
                          betaSampledY,betaSampledUniqueY,gamSampledY,gamSampledUniqueY,
                          betaSampledCom,betaSampledUniqueCom,gamSampledCom,gamSampledUniqueCom,
                          sigX2SampledUnique,sigY2SampledUnique,nx,data,probTmp)
              order<-out$order
              betaSampledX<-out$betaSampledX
              betaSampledUniqueX<-out$betaSampledUniqueX
              gamSampledX<-out$gamSampledX
              gamSampledUniqueX<-out$gamSampledUniqueX

              betaSampledY<-out$betaSampledY
              betaSampledUniqueY<-out$betaSampledUniqueY
              gamSampledY<-out$gamSampledY
              gamSampledUniqueY<-out$gamSampledUniqueY

              betaSampledCom<-out$betaSampledCom
              betaSampledUniqueCom<-out$betaSampledUniqueCom
              gamSampledCom<-out$gamSampledCom
              gamSampledUniqueCom<-out$gamSampledUniqueCom

              etaSampled<-out$etaSampled

              sigX2SampledUnique<-out$sigX2SampledUnique
              sigY2SampledUnique<-out$sigY2SampledUnique
              probTmp<-out$probTmp
            } # end if (pee<runif(1))
            else
            {
              #				cat("Not pee<runif(1)","\n")

              if (iii==((K-k)*(B+N)+1))
              {
                orderTmp<-order[k,1,]
                gamX<-gamSampledUniqueX[1,]
                betaX<-betaSampledUniqueX[1,]

                gamY<-gamSampledUniqueY[1,]
                betaY<-betaSampledUniqueY[1,]

                gamCom<-gamSampledUniqueCom[1,]
                betaCom<-betaSampledUniqueCom[1,]

                sigX2<-sigX2SampledUnique[1]
                sigY2<-sigY2SampledUnique[1]
              }
              else
              {
                orderTmp<-order[k,(iii-1),]
                gamX<-gamSampledX[k,(iii-1),2:ncol(gamSampledX[k,,])]
                betaX<-betaSampledX[k,(iii-1),2:ncol(betaSampledX[k,,])]

                gamY<-gamSampledY[k,(iii-1),2:ncol(gamSampledY[k,,])]
                betaY<-betaSampledY[k,(iii-1),2:ncol(betaSampledY[k,,])]

                gamCom<-gamSampledCom[k,(iii-1),2:ncol(gamSampledCom[k,,])]
                betaCom<-betaSampledCom[k,(iii-1),2:ncol(betaSampledCom[k,,])]

                sigX2<-sigX2SampledUnique[order[k,(iii-1),1]]
                sigY2<-sigY2SampledUnique[order[k,(iii-1),1]]
              }
              dataTmp<-data[,orderTmp[2:length(orderTmp)]]
              probO<-OrderProb(gamX,betaX,gamY,betaY,betaCom,sigma2,sigX2,sigY2,n,dataTmp,nx)


              for (j in 1:(K-1))
              {
                if (H[j]<(-probO$prob) && (-probO$prob)<H[j+1])
                {
                  ring<-j
                  break
                }
              }
              if ((-probO$prob)>H[K])
              {
                ring<-K
              }

              # sample from rings
              selectRing<-sample(seq(ring,K),1)
              len<-length(D[,2,selectRing])-sum(is.na(D[,2,selectRing]))
              if (len==0)
              {
                out<-MHStep(order,k,iii,nodes,etaSampled,betaSampledX,betaSampledUniqueX,gamSampledX,gamSampledUniqueX,
                            betaSampledY,betaSampledUniqueY,gamSampledY,gamSampledUniqueY,
                            betaSampledCom,betaSampledUniqueCom,gamSampledCom,gamSampledUniqueCom,
                            sigX2SampledUnique,sigY2SampledUnique,nx,data,probTmp)

                order<-out$order
                betaSampledX<-out$betaSampledX
                betaSampledUniqueX<-out$betaSampledUniqueX
                gamSampledX<-out$gamSampledX
                gamSampledUniqueX<-out$gamSampledUniqueX

                betaSampledY<-out$betaSampledY
                betaSampledUniqueY<-out$betaSampledUniqueY
                gamSampledY<-out$gamSampledY
                gamSampledUniqueY<-out$gamSampledUniqueY

                betaSampledCom<-out$betaSampledCom
                betaSampledUniqueCom<-out$betaSampledUniqueCom
                gamSampledCom<-out$gamSampledCom
                gamSampledUniqueCom<-out$gamSampledUniqueCom

                etaSampled<-out$etaSampled

                sigX2SampledUnique<-out$sigX2SampledUnique
                sigY2SampledUnique<-out$sigY2SampledUnique
              }
              else
              {

                loc1<-which(is.na(D[,2,selectRing]))
                #				loc<-sample((D[-loc1,2,selectRing]),1)

                cand<-D[-loc1,2,selectRing]
                if (length(cand)==1)
                {
                  loc<-which((D[,2,selectRing])==cand[1])
                }
                else
                {
                  sampleIndex<-sample(x=cand,size=1)
                  loc<-which((D[,2,selectRing])==sampleIndex)
                }


                orderIndex<-D[loc[1],2,selectRing]
                sigX2Tmp<-sigX2SampledUnique[orderIndex]
                sigY2Tmp<-sigY2SampledUnique[orderIndex]

                gamXTmp<-gamSampledUniqueX[orderIndex,]
                betaXTmp<-betaSampledUniqueX[orderIndex,]

                gamYTmp<-gamSampledUniqueY[orderIndex,]
                betaYTmp<-betaSampledUniqueY[orderIndex,]

                gamComTmp<-gamSampledUniqueCom[orderIndex,]
                betaComTmp<-betaSampledUniqueCom[orderIndex,]

                orderTmp<-D[loc[1],2:(nodes+2),selectRing]
                #				cat("Here length(orderTmp) ",length(orderTmp)," ",orderTmp,"\n")
                dataTmp<-data[,orderTmp[2:length(orderTmp)]]

                probTmp<-OrderProb(gamXTmp,betaXTmp,gamYTmp,betaYTmp,betaComTmp,sigma2,sigX2Tmp,sigY2Tmp,n,dataTmp,nx)

                lpyCurrent<--max(-probTmp$prob,H[(k)])/T[(k)]
                lpyHigh<--max(-probTmp$prob,H[(k+1)])/T[(k+1)]
                lpxCurrent<--max(-probO$prob,H[(k)])/T[(k)]
                lpxHigh<--max(-probO$prob,H[(k+1)])/T[(k+1)]
                lratio<-lpyCurrent+lpxHigh-lpxCurrent-lpyHigh
                if (exp(lratio)>runif(1))
                {
                  order[k,iii,]<-D[loc[1],2:(nodes+2),selectRing]
                  betaSampledUniqueX[order[k,iii,1],]<-probTmp$avgBetaX
                  gamSampledUniqueX[order[k,iii,1],]<-probTmp$avgGamX

                  betaSampledUniqueY[order[k,iii,1],]<-probTmp$avgBetaY
                  gamSampledUniqueY[order[k,iii,1],]<-probTmp$avgGamY

                  betaSampledUniqueCom[order[k,iii,1],]<-probTmp$avgBetaCom
                  gamSampledUniqueCom[order[k,iii,1],]<-probTmp$avgGamCom

                  sigX2SampledUnique[order[k,iii,1]]<-probTmp$sigX2
                  sigY2SampledUnique[order[k,iii,1]]<-probTmp$sigY2

                  betaSampledX[k,iii,2:ncol(betaSampledX[k,,])]<-probTmp$avgBetaX
                  betaSampledX[k,iii,1]<-order[k,iii,1]
                  gamSampledX[k,iii,2:ncol(gamSampledX[k,,])]<-probTmp$avgGamX
                  gamSampledX[k,iii,1]<-order[k,iii,1]

                  betaSampledY[k,iii,2:ncol(betaSampledY[k,,])]<-probTmp$avgBetaY
                  betaSampledY[k,iii,1]<-order[k,iii,1]
                  gamSampledY[k,iii,2:ncol(gamSampledY[k,,])]<-probTmp$avgGamY
                  gamSampledY[k,iii,1]<-order[k,iii,1]

                  betaSampledCom[k,iii,2:ncol(betaSampledCom[k,,])]<-probTmp$avgBetaCom
                  betaSampledCom[k,iii,1]<-order[k,iii,1]
                  gamSampledCom[k,iii,2:ncol(gamSampledCom[k,,])]<-probTmp$avgGamCom
                  gamSampledCom[k,iii,1]<-order[k,iii,1]

                  etaSampled[k,iii,2]<-probTmp$eta
                  etaSampled[k,iii,1]<-order[k,iii,1]
                }
                else
                {
                  if (iii==((K-k)*(B+N)+1))
                  {
                    order[k,iii,]<-order[k,1,]
                  }
                  else
                  {
                    order[k,iii,]<-order[k,iii-1,]
                  }
                  betaSampledUniqueX[order[k,iii,1],]<-probO$avgBetaX
                  gamSampledUniqueX[order[k,iii,1],]<-probO$avgGamX

                  betaSampledUniqueY[order[k,iii,1],]<-probO$avgBetaY
                  gamSampledUniqueY[order[k,iii,1],]<-probO$avgGamY

                  betaSampledUniqueCom[order[k,iii,1],]<-probO$avgBetaCom
                  gamSampledUniqueCom[order[k,iii,1],]<-probO$avgGamCom

                  sigX2SampledUnique[order[k,iii,1]]<-probO$sigX2
                  sigY2SampledUnique[order[k,iii,1]]<-probO$sigY2

                  betaSampledX[k,iii,2:ncol(betaSampledX[k,,])]<-probO$avgBetaX
                  betaSampledX[k,iii,1]<-order[k,iii,1]
                  gamSampledX[k,iii,2:ncol(gamSampledX[k,,])]<-probO$avgGamX
                  gamSampledX[k,iii,1]<-order[k,iii,1]

                  betaSampledY[k,iii,2:ncol(betaSampledY[k,,])]<-probO$avgBetaY
                  betaSampledY[k,iii,1]<-order[k,iii,1]
                  gamSampledY[k,iii,2:ncol(gamSampledY[k,,])]<-probO$avgGamY
                  gamSampledY[k,iii,1]<-order[k,iii,1]

                  betaSampledCom[k,iii,2:ncol(betaSampledCom[k,,])]<-probO$avgBetaCom
                  betaSampledCom[k,iii,1]<-order[k,iii,1]
                  gamSampledCom[k,iii,2:ncol(gamSampledCom[k,,])]<-probO$avgGamCom
                  gamSampledCom[k,iii,1]<-order[k,iii,1]

                  etaSampled[k,iii,2]<-probO$eta
                  etaSampled[k,iii,1]<-order[k,iii,1]
                }
              } # end else
            } # end else
            #		cat("k ",k," iii ",iii," etaSampled ",etaSampled[k,iii,2],"\n")
          } # end else corresponding to if (k==K)
        } # end if (iii>((K-k)*(B+N))
        if (iii>((K-k)*(B+N)+B))
        {
          sigX2<-sigX2SampledUnique[order[k,iii,1]]
          sigY2<-sigY2SampledUnique[order[k,iii,1]]

          gamX<-gamSampledX[k,iii,2:ncol(gamSampledX[k,,])]
          betaX<-betaSampledX[k,iii,2:ncol(betaSampledX[k,,])]

          gamY<-gamSampledY[k,iii,2:ncol(gamSampledY[k,,])]
          betaY<-betaSampledY[k,iii,2:ncol(betaSampledY[k,,])]

          gamCom<-gamSampledCom[k,iii,2:ncol(gamSampledCom[k,,])]
          betaCom<-betaSampledCom[k,iii,2:ncol(betaSampledCom[k,,])]

          orderTmp<-order[k,iii,]
          dataTmp<-data[,orderTmp[2:length(orderTmp)]]

          probO<-OrderProb(gamX,betaX,gamY,betaY,betaCom,sigma2,sigX2,sigY2,n,dataTmp,nx)

          for (j in 1:(K-1))
          {
            if (H[j]<(-probO$prob)&& (-probO$prob)<H[(j+1)])
            {
              D[Dindex[j],1,j]<-k
              D[Dindex[j],2,j]<-order[k,iii,1]
              D[Dindex[j],3:(nodes+2),j]<-order[k,iii,2:(nodes+1)]
              Dindex[j]<-Dindex[j]+1
              break
            }
          }
          if ((-probO$prob)>H[K])
          {
            #			cat("Dindex[K] ",Dindex[K],"\n")
            D[Dindex[K],1,K]<-k
            D[Dindex[K],2,K]<-order[k,iii,1]
            #			cat("here ",D[Dindex[K],1,K]," K ",K,"\n")
            D[Dindex[K],3:(nodes+2),K]<-order[k,iii,2:(nodes+1)]
            Dindex[K]<-Dindex[K]+1
          }
        } # end if (iii>((K-k)*(B+N)+B)
      } # end for (k in K:0)


    } # end for (iii in 2:orderloops)
    #}


    parentsTrueGX<-parentsTrueGY<-parentsTrueGCom<-matrix(rep(0,nodes*nodes),nrow=nodes,ncol=nodes)

    parentsEstGX<-parentsEstGY<-parentsEstGCom<-matrix(rep(0,nodes*nodes),nrow=nodes,ncol=nodes)

    parentsTrueBX<-parentsTrueBY<-parentsTrueBCom<-matrix(rep(0,nodes*nodes),nrow=nodes,ncol=nodes)
    parentsEstBX<-parentsEstBY<-parentsEstBCom<-matrix(rep(0,nodes*nodes),nrow=nodes,ncol=nodes)

    #avgEtaEst<-1

    orderEst<-rep(0,nodes)




    for (kk in 0:(nodes-2))
    {
      start<-sum(0:kk)+1
      end<-start+kk
      t<-1
      for (kkk in start:end)
      {
        parentsTrueGX[(kk+2),t]<-gamTrueX[1,kkk]
        parentsTrueGY[(kk+2),t]<-gamTrueY[1,kkk]
        parentsTrueGCom[(kk+2),t]<-gamTrueCom[1,kkk]

        parentsTrueBX[(kk+2),t]<-trueBetaX[1,kkk]
        parentsTrueBY[(kk+2),t]<-trueBetaY[1,kkk]
        parentsTrueBCom[(kk+2),t]<-trueBetaCom[1,kkk]

        t<-t+1
      }
    }

    Sumeta<-sumEtaEst<-0
    counts<-0



    for (iii in ((K-1)*(B+N)+B):((K-1)*(B+N)+B+B+2*N))
    {
      for (t in 1:K)
      {
        if (counts==0)
        {
          sumEtaEst<-etaSampled[t,iii,2]
          counts<-counts+1
        }
        else
        {
          sumEtaEst<-etaSampled[t,iii,2]+sumEtaEst
          counts<-counts+1
        }
      }
    }
    avgEtaEst[jj]<-sumEtaEst/counts
    cat("avgEtaEst jj ",jj," ",avgEtaEst[jj],"\n")



    counts<-0

    for (iii in ((K-1)*(B+N)+B):((K-1)*(B+N)+B+B+2*N))
    {

      for (t in 1:K)
      {
        for (kk in 0:(nodes-2))
        {
          orderEst[1:(kk+2)]<-order[t,iii, 2:(kk+3)]
          start<-sum(0:kk)+1+1
          end<-start+kk
          tt<-1
          for (kkk in start:end)
          {
            #				cat("orderEst[kk+2] ",orderEst[kk+2]," orderEst[tt] ",orderEst[tt],"\n")
            #				cat("gamSampledX[t,iii,kkk] ",gamSampledX[t,iii,kkk],"\n")

            if (avgEtaEst[jj]<=0.5)
            {
              if (etaSampled[t,iii,2]==0)
              {
                parentsEstGX[orderEst[kk+2],orderEst[tt]]<-gamSampledX[t,iii,kkk]
                parentsEstGY[orderEst[kk+2],orderEst[tt]]<-gamSampledY[t,iii,kkk]
                parentsEstBX[orderEst[kk+2],orderEst[tt]]<-betaSampledX[t,iii,kkk]
                parentsEstBY[orderEst[kk+2],orderEst[tt]]<-betaSampledY[t,iii,kkk]
                tt<-tt+1
              }
            }
            else
            {
              if (etaSampled[t,iii,2]==1)
              {
                parentsEstGCom[orderEst[kk+2],orderEst[tt]]<-gamSampledCom[t,iii,kkk]
                parentsEstBCom[orderEst[kk+2],orderEst[tt]]<-betaSampledCom[t,iii,kkk]
                tt<-tt+1
              }
            }


          }
          #			print(parentsEstGX)
        }
        if (counts==0)
        {
          if (avgEtaEst[jj]<=0.5 && etaSampled[t,iii,2]==0)
          {
            sumParentsEstGX<-parentsEstGX
            sumParentsEstGY<-parentsEstGY
            sumParentsEstBX<-parentsEstBX
            sumParentsEstBY<-parentsEstBY
            counts<-counts+1
          }
          if (avgEtaEst[jj]>0.5 && etaSampled[t,iii,2]==1)
          {
            sumParentsEstGCom<-parentsEstGCom
            sumParentsEstBCom<-parentsEstBCom
            counts<-counts+1
          }
        }
        else
        {
          if (avgEtaEst[jj]<=0.5 && etaSampled[t,iii,2]==0)
          {
            sumParentsEstGX<-parentsEstGX+sumParentsEstGX
            sumParentsEstGY<-parentsEstGY+sumParentsEstGY
            sumParentsEstBX<-parentsEstGX+sumParentsEstBX
            sumParentsEstBY<-parentsEstGY+sumParentsEstBY
            counts<-counts+1
          }
          if (avgEtaEst[jj]>0.5 && etaSampled[t,iii,2]==1)
          {
            sumParentsEstGCom<-parentsEstGCom+sumParentsEstGCom
            sumParentsEstBCom<-parentsEstGCom+sumParentsEstBCom
            counts<-counts+1
          }
        }
      }
    }

    if (avgEtaEst[jj]<=0.5)
    {
      avgParentsEstGX<-round(sumParentsEstGX/counts,0)
      diffGX<-avgParentsEstGX-parentsTrueGX
      accurGX[jj]<-1-sum(abs(diffGX))/(nodes*(nodes-1)/2)
      fpGX[jj]<-sum(diffGX>0)/(nodes*(nodes-1)/2-sum(parentsTrueGX))
      tpGX[jj]<-1-abs(sum(diffGX<0))/sum(parentsTrueGX)
      avgParentsEstBX<-round(sumParentsEstBX/counts,0)
      diffBX<-avgParentsEstBX+(parentsTrueGX-1)
      tpBX[jj]<-length(which(diffBX==1))/sum(parentsTrueGX)

      avgParentsEstGY<-round(sumParentsEstGY/counts,0)
      diffGY<-avgParentsEstGY-parentsTrueGY
      accurGY[jj]<-1-sum(abs(diffGY))/(nodes*(nodes-1)/2)
      fpGY[jj]<-sum(diffGY>0)/(nodes*(nodes-1)/2-sum(parentsTrueGY))
      tpGY[jj]<-1-abs(sum(diffGY<0))/sum(parentsTrueGY)
      avgParentsEstBY<-round(sumParentsEstBY/counts,0)
      diffBY<-avgParentsEstBY+(parentsTrueGY-1)
      tpBY[jj]<-length(which(diffBY==1))/sum(parentsTrueGY)
    }
    else
    {
      avgParentsEstGCom<-round(sumParentsEstGCom/counts,0)
      diffGCom<-avgParentsEstGCom-parentsTrueGCom
      accurGCom[jj]<-1-sum(abs(diffGCom))/(nodes*(nodes-1)/2)
      fpGCom[jj]<-sum(diffGCom>0)/(nodes*(nodes-1)/2-sum(parentsTrueGCom))
      tpGCom[jj]<-1-abs(sum(diffGCom<0))/sum(parentsTrueGCom)
      avgParentsEstBCom<-round(sumParentsEstBCom/counts,0)
      diffBCom<-avgParentsEstBCom+(parentsTrueGCom-1)
      tpBCom[jj]<-length(which(diffBCom==1))/sum(parentsTrueGCom)
    }


  } # end data=jj


  save(avgEtaEst,accurGCom,fpGCom,tpGCom,tpBCom,accurGX,fpGX,tpGX,tpBX,accurGY,fpGY,tpGY,tpBY,file=paste(basefn,".Rdata",sep=""))

  # summary for beta
  cat("summary for eta \n")
  quantile(avgEtaEst,prob=c(0,0.025,0.05,0.5,0.95,0.975,1))
  mean(avgEtaEst)
  sqrt(var(avgEtaEst))

  #cat("identical networks \n")
  # false positives edges
  cat("summary for fpGCom \n")
  quantile(fpGCom[which(avgEtaEst>0.5)],prob=c(0,0.025,0.05,0.5,0.95,0.975,1))
  mean(fpGCom[which(avgEtaEst>0.5)])
  sqrt(var(fpGCom[which(avgEtaEst>0.5)]))

  # true positives edges
  cat("summary for tpGCom \n")
  quantile(tpGCom[which(avgEtaEst>0.5)],prob=c(0,0.025,0.05,0.5,0.95,0.975,1))
  mean(tpGCom[which(avgEtaEst>0.5)])
  sqrt(var(tpGCom[which(avgEtaEst>0.5)]))

  # true positives edges and directions
  cat("summary for tpBCom \n")
  quantile(tpBCom[which(avgEtaEst>0.5)],prob=c(0,0.025,0.05,0.5,0.95,0.975,1))
  mean(tpBCom[which(avgEtaEst>0.5)])
  sqrt(var(tpBCom[which(avgEtaEst>0.5)]))

  cat("differential networks \n")
  # differential networks
  # false positives edges X
  cat("summary for fpGX \n")
  quantile(fpGX[which(avgEtaEst<=0.5)],prob=c(0,0.025,0.05,0.5,0.95,0.975,1))
  mean(fpGX[which(avgEtaEst<=0.5)])
  sqrt(var(fpGX[which(avgEtaEst<=0.5)]))

  # true positives edges X
  cat("summary for tpGX \n")
  quantile(tpGX[which(avgEtaEst<=0.5)],prob=c(0,0.025,0.05,0.5,0.95,0.975,1))
  mean(tpGX[which(avgEtaEst<=0.5)])
  sqrt(var(tpGX[which(avgEtaEst<=0.5)]))

  # true positives edges and directions X
  cat("summary for tpBX \n")
  quantile(tpBX[which(avgEtaEst<=0.5)],prob=c(0,0.025,0.05,0.5,0.95,0.975,1))
  mean(tpBX[which(avgEtaEst<=0.5)])
  sqrt(var(tpBX[which(avgEtaEst<=0.5)]))

  # false positives edges Y
  cat("summary for fpGY \n")
  quantile(fpGY[which(avgEtaEst<=0.5)],prob=c(0,0.025,0.05,0.5,0.95,0.975,1))
  mean(fpGY[which(avgEtaEst<=0.5)])
  sqrt(var(fpGY[which(avgEtaEst<=0.5)]))

  # true positives edges Y
  cat("summary for tpGY \n")
  quantile(tpGY[which(avgEtaEst<=0.5)],prob=c(0,0.025,0.05,0.5,0.95,0.975,1))
  mean(tpGY[which(avgEtaEst<=0.5)])
  sqrt(var(tpGY[which(avgEtaEst<=0.5)]))

  # true positives edges and directions Y
  cat("summary for tpBY \n")
  quantile(tpBY[which(avgEtaEst<=0.5)],prob=c(0,0.025,0.05,0.5,0.95,0.975,1))
  mean(tpBY[which(avgEtaEst<=0.5)])
  sqrt(var(tpBY[which(avgEtaEst<=0.5)]))


  #sink()



}
