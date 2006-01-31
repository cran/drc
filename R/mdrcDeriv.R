"mdrcDeriv" <- function()
{
    ## Handling derivatives with and without Box-Cox transformation ... experimental at the moment
    if (derFlag)
    {
        if (!boxcox)
        {
            ## Handling the first derivative
            opfctDer<-function(parm)
            {
                #parmVal <- parm2mat(parm)
                resVec<-resp-fctEval
                fctDer<-drcDer(dose,parmVal)

                retList<-list()
                for (i in 1:numNames)
                {
                    retList[[i]]<-tapply(resVec*fctDer[,i],collapse[,i],sum, na.rm=T)  # na.rm=T is used to remove zero term appearing as NaN in expLog
                }
                -2*c(unlist(retList))
            }

            ## Handling the second derivative
            opfctDer2<-function(parmVal)
            {
                #parmVal <- parm2mat(parm)
                resVec<-resp-fctEval
                fctDer<-drcDer(dose,parmVal)
                fctDer2<-drcDer2(dose,parmVal)

                retMat<-matrix(0,csVec[numNames+1],csVec[numNames+1])
                k<-1
                for (i in 1:numNames)
                {
                    for (j in i:numNames)
                    {
                        tempVec<-tapply(-fctDer[,i]*fctDer[,j]+resVec*fctDer2[,k],factor(paste(collapse[,i],"x",collapse[,j])),sum,na.rm=T)

                        if (length(c((csVec[i]+1):csVec[i+1]))==1)
                        {
                            retMat[c((csVec[i]+1):csVec[i+1]),c((csVec[j]+1):(csVec[j+1]))]<-tempVec
                        } else {
                            diag(retMat[c((csVec[i]+1):csVec[i+1]),c((csVec[j]+1):(csVec[j+1]))])<-tempVec}
                            # special behaviour of 'diag' for vectors of length one
                        k<-k+1
                    }
                }
                retMat[lower.tri(retMat)]<-t(retMat)[lower.tri(retMat)]

                -2*retMat
            }
        } else { # in case Box-Cox transformation is applied

            if (abs(lambda)>bcTol) {fctEval2<-(fctEval^lambda-1)/lambda} else {fctEval2<-log(fctEval)}

            ## Handling the first derivative
            opfctDer<-function(parm)
            {
                #parmVal <- parm2mat(parm)
                resVec<-resp-fctEval2
                fctDer<-drcDer(dose,parmVal)

                retList<-list()
                for (i in 1:numNames)
                {
                    retList[[i]]<-tapply(resVec*fctDer[,i]*(fctEval^(lambda-1)),collapse[,i],sum, na.rm=T)
                }
                -2*c(unlist(retList))
            }

            ## Handling the second derivative
            opfctDer2<-function(parmVal)
            {
                #parmVal <- parm2mat(parm)
                resVec<-resp-fctEval2
                fctDer<-drcDer(dose,parmVal)
                fctDer2<-drcDer2(dose,parmVal)

                retMat<-matrix(0,csVec[numNames+1],csVec[numNames+1])
                k<-1
                for (i in 1:numNames)
                {
                    for (j in i:numNames)
                    {
                        tempVec<-tapply(-fctDer[,i]*fctDer[,j]*(fctEval^(lambda-2))+resVec*((lambda-1)*(fctEval^(lambda-2))*fctDer[,i]*fctDer[,j]
                                        +fctDer2[,k]*(fctEval^(lambda-1))),factor(paste(collapse[,i],"x",collapse[,j])),sum,na.rm=T)

                        if (length(c((csVec[i]+1):csVec[i+1]))==1)
                        {
                            retMat[c((csVec[i]+1):csVec[i+1]),c((csVec[j]+1):(csVec[j+1]))]<-tempVec
                        } else {
                            diag(retMat[c((csVec[i]+1):csVec[i+1]),c((csVec[j]+1):(csVec[j+1]))])<-tempVec}
                        k<-k+1
                    }
                }
                retMat[lower.tri(retMat)]<-t(retMat)[lower.tri(retMat)]
                -2*retMat
            }
        }
    }  # end of handling derivatives
}
