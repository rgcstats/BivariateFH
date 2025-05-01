
library(rpart)
library(rpart.plot)


rho.e.values=c(0,0.1,0.3,0.5,0.7,0.9) 
phi.11.values<-1
phi.22.values<-c(1,1.5,3,5,10)

rho.v.values=c(0,0.1,0.3,0.5,0.7,0.9)
r.values=c(1/10,1/5,1/3,2/3,1,1.5,3,5,10)     # r is the ration of two variances i.e r=sigma22/sigma11
sigma.11.values<-c(1/10,1/5,1/3,2/3,1,1.5,3,5,10)


MSE.BLUP.UFH<-NULL 
MSE.BLUP.BFH<-NULL 
MSE.BLUP.UFH.determinant <- MSE.BLUP.BFH.determinant <- NULL

for(i6 in c(1:length(sigma.11.values))){  
  sigma.11.values.1<-sigma.11.values[i6]       
  
  for(i5 in c(1:length(r.values))){
    r.values.1<- r.values[i5]
    
    for(i4 in c(1:length(rho.v.values))){
      rho.v.values.1<-rho.v.values[i4]
      
      for(i3 in c(1:length(phi.22.values))){
        phi.22.values.1<-  phi.22.values[i3]
        
        for(i2 in c(1:length(phi.11.values))){
          phi.11.values.1<-  phi.11.values[i2]
          
          for(i1 in c(1:length(rho.e.values))){
            
            phi.12.values<- rho.e.values[i1]*sqrt(phi.11.values.1*phi.22.values.1)
            phi.21.values<-phi.12.values
            phi.mat<-matrix(c(phi.11.values.1,phi.12.values,phi.21.values, phi.22.values.1),nrow=2,ncol=2)
            sigma.22.values<-r.values.1*sigma.11.values.1
            sigma.12.values<-rho.v.values.1*sqrt(sigma.11.values.1*sigma.22.values)
            sigma.21.values<-sigma.12.values
            sigma.mat<-matrix(data=c(sigma.11.values.1,sigma.12.values,sigma.21.values, sigma.22.values),nrow=2,ncol=2) 
            
            MSE.BLUP.UFH.1<-(diag(sigma.mat)*diag(phi.mat))/(diag(sigma.mat)+diag(phi.mat))      # formula of UFH
            MSE.BLUP.BFH.1<-diag(sigma.mat-t(sigma.mat)%*%solve(phi.mat+sigma.mat)%*%sigma.mat)  # formula of BFH
            
            MSE.BLUP.UFH.crossprod12 <- (phi.mat[1,1]*phi.mat[2,2]*sigma.mat[1,2]+sigma.mat[1,1]*sigma.mat[2,2]*phi.mat[1,2])/
              (sigma.mat[1,1]+phi.mat[1,1])/(sigma.mat[2,2]+phi.mat[2,2])
            MSE.BLUP.UFH.determinant.1 <- MSE.BLUP.UFH.1[1]*MSE.BLUP.UFH.1[2]-MSE.BLUP.UFH.crossprod12^2
            
            MSE.BLUP.BFH.matrix <- sigma.mat-t(sigma.mat)%*%solve(phi.mat+sigma.mat)%*%sigma.mat
            MSE.BLUP.BFH.determinant.1 <- det(MSE.BLUP.BFH.matrix)
           
            # MSE without determinant
            
            MSE.BLUP.UFH<-rbind(MSE.BLUP.UFH,MSE.BLUP.UFH.1)   
            MSE.BLUP.BFH<-rbind(MSE.BLUP.BFH,MSE.BLUP.BFH.1)
            
            # MSE with determinant
            
            MSE.BLUP.UFH.determinant <- c(MSE.BLUP.UFH.determinant,MSE.BLUP.UFH.determinant.1)
            MSE.BLUP.BFH.determinant <- c(MSE.BLUP.BFH.determinant,MSE.BLUP.BFH.determinant.1)
           
          }
        }
      }
    }
  }
}

results<-data.frame(rho.e=rep(rho.e.values,length(phi.11.values)*length(phi.22.values)*length(rho.v.values)*length(r.values)*length(sigma.11.values)),
                      phi11=rep(rep(phi.11.values,each=length(rho.e.values)),length(phi.22.values)*length(rho.v.values)*length(r.values)*length(sigma.11.values)),
                      phi22=rep(rep(phi.22.values,each=length(rho.e.values)*length(phi.11.values)),length(rho.v.values)*length(r.values)*length(sigma.11.values)),
                      rho.v=rep(rep(rho.v.values,each=length(rho.e.values)*length(phi.11.values)*length(phi.22.values)),length(r.values)*length(sigma.11.values)),  
                      r=rep(rep(r.values,each=length(rho.e.values)*length(phi.11.values)*length(phi.22.values)*length(rho.v.values)),length(sigma.11.values)),
                      sigma11=rep(sigma.11.values,each=length(rho.e.values)*length(phi.11.values)*length(phi.22.values)*length(rho.v.values)*length(r.values)),
                      UFH.1=MSE.BLUP.UFH[,1],UFH.2=MSE.BLUP.UFH[,2],BFH.1=MSE.BLUP.BFH[,1],BFH.2=MSE.BLUP.BFH[,2],MSE.BLUP.UFH.determinant=MSE.BLUP.UFH.determinant,MSE.BLUP.BFH.determinant=MSE.BLUP.BFH.determinant)

# Graph from Numerical study
results$gain<-results$MSE.BLUP.BFH.determinant/results$MSE.BLUP.UFH.determinant
results$r.s<-results$phi22/results$phi11
results$ratio<-round(results$r.s/results$r,1)
results$sigma22<-results$r*results$sigma11
results$gain.1<-results$BFH.1/results$UFH.1
results$gain.2<-results$BFH.2/results$UFH.2

pdf("REN01.pdf", width = 16, height = 8) 
par(mfrow=c(1,1),mar=c(5,7.2,4.1,0.5),mgp=c(4, 1, 0))

boxplot(results[,"MSE.BLUP.BFH.determinant"]/results[,"MSE.BLUP.UFH.determinant"]~results[,"rho.v"],ylab="RE",ylim=c(0,3),xlab=expression(rho[v]),outline=FALSE,
        main="First Variable",cex.axis=1.5, cex.lab = 1.5, cex.main=1.5)
dev.off()

pdf("REN02.pdf", width = 16, height = 8)
par(mfrow=c(1,1),mar=c(5,7.2,4.1,0.5),mgp=c(4, 1, 0))

boxplot(results[,"MSE.BLUP.BFH.determinant"]/results[,"MSE.BLUP.UFH.determinant"]~results[,"rho.e"],ylab="RE",xlab= expression(rho[e]),ylim=c(0.2,1),outline=FALSE,
        main="First Variable",cex.axis=1.5, cex.lab = 1.5, cex.main=1.5)
dev.off()


pdf("GainScatterplot.pdf", width = 12, height = 10) 
par(mfrow=c(1,1), mar=c(5,7.2,4.1,0.5),mgp=c(4, 1, 0)) 
plot(results$gain.1,results$gain.2, ylab="RE of Second Variable", xlab=" RE of First Variable",
     cex.axis=1.7, cex.lab = 1.7)
abline(0,1)
dev.off()


pdf("RENRatioR.pdf", width = 20, height = 10)
par(mfrow=c(1,1), mar=c(5,7.2,4.1,0.5),mgp=c(4, 1, 0),las=2) 

boxplot(results[,"MSE.BLUP.BFH.determinant"]/results[,"MSE.BLUP.UFH.determinant"]~results[,"ratio"], ylab="RE",xlab="R",ylim=c(0.2,1),outline=FALSE,
        main="First Variable",cex.axis=1.5, cex.lab = 1.5, cex.main=1.7)

dev.off()

pdf("REN07.pdf", width = 20, height = 10) 
par(mfrow=c(1,1))

sub.results<-(results$rho.e==0.7) & (results$rho.v==0.7)  # subsection of results.2
plot(results$ratio[sub.results],results$gain[sub.results],log="x",xlab="R",ylab="RE",main="First Variable",cex.axis=1.5, cex.lab = 1.6, cex.main=1.7) 

dev.off()


pdf("REN08.pdf", width = 20, height = 10) 
par(mfrow=c(1,1))

sub.results<-(results$rho.e==0.1) & (results$rho.v==0.7)  
plot(results$ratio[sub.results],results$gain[sub.results],log="x",xlab="R",ylim=c(0.2,1),ylab="RE",main="First Variable",cex.axis=1.5, cex.lab = 1.5, cex.main=1.7) 

dev.off()

pdf("REN09.pdf", width = 20, height = 10) 
par(mfrow=c(1,1))

sub.results<-(results$rho.e==0.7) & (results$rho.v==0.1)  
plot(results$ratio[sub.results],results$gain.1[sub.results],log="x",xlab="R",ylim=c(0.2,1),ylab="RE",main="First Variable",cex.axis=1.5, cex.lab = 1.5, cex.main=1.7) 

dev.off()

# Regression Tree Analysis

fit.final.1 <- rpart(gain.1~ratio+rho.e+rho.v+sigma11+sigma22+r+r.s+phi11+phi22,data=results,method="anova",control=list(maxdepth=4)) 

pdf("REGTreeFirstVar.pdf", width=14,height=12) 
par(mfrow=c(1,1))
prp(fit.final.1, main="Regression Tree of relative rfficiency (RE) for First Variable",
    extra=100, # display number and percent of obs
    compress=TRUE,
    split.round=1,
    type = 1, # Label all nodes, not just leaves
    nn=TRUE, # display the node numbers
    nn.font=4,  #Font for the node numbers. Default 3, italic
    branch=1,# change angle of branch lines 
    Margin=0, #Extra white space around the tree, as a fraction of the graph width. Default 0,meaning no extra space. To add say 10% space around the tree use Margin=0.1.
    cex=1.5, #display size of the split
    cex.main=1.5,
    #trace=1, # print the automatically calculated cex
    fallen.leaves=TRUE,#Default FALSE. If TRUE, display the leaves at the bottom of the graph.
    yesno=TRUE, #Default TRUE, meaning write yes and no on the appropriate sides of the top split.
    xsep="/",#separates the factor levels in split labels
    split.cex=0.9, # make the split text larger than the node text
    uniform=T) #If TRUE (the default), the vertical spacing of the nodes is uniform
dev.off()


fit.final.2 <- rpart(gain.2~ratio+rho.e+rho.v+sigma11+sigma22+r+r.s+phi11+phi22,data=results,method="anova",control=list(maxdepth=4)) 
pdf("REGTreeSecondVar.pdf", width=12,height=10) 
par(mfrow=c(1,1))
prp(fit.final.2, main="Regression Tree of relative rfficiency (RE) for Second Variable",
    extra=100, # display number and percent of obs
    compress=TRUE,
    split.round=1,
    type = 1, # Label all nodes, not just leaves
    nn=TRUE, # display the node numbers
    nn.font=4,  #Font for the node numbers. Default 3, italic
    branch=1,# change angle of branch lines 
    Margin=0, #Extra white space around the tree, as a fraction of the graph width. Default 0,meaning no extra space. To add say 10% space around the tree use Margin=0.1.
    cex=1.5, #display size of the split
    cex.main=1.5,
    #trace=1, # print the automatically calculated cex
    fallen.leaves=TRUE,#Default FALSE. If TRUE, display the leaves at the bottom of the graph.
    yesno=TRUE, #Default TRUE, meaning write yes and no on the appropriate sides of the top split.
    xsep="/",#separates the factor levels in split labels
    split.cex=0.9, # make the split text larger than the node text
    uniform=T) #If TRUE (the default), the vertical spacing of the nodes is uniform
dev.off()

