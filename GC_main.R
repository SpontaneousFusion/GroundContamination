# if(dev.list())
# cat(1)
# dev.off()
# install.packages("kriging")
rm(list=ls())
gc()
graphics.off()
### INITIALIZATION LOADING LIBRARIES ###
# install.packages("numDeriv")
# install.packages("matlib")
library(mvtnorm)
library(RColorBrewer)



# Total fluence vs frid distance
# 
# 
# 

fCalcLg=function(approx.det.h=5,middle=F,points.to.calculate=500,x.sc=1,y.sc=1)
{
  
  lg.sizes=c(1,1.5,2,2.5,3,3.5,4,4.5,seq(5,100))
  lg.sizes
  tot.fluence.vector=c()
  ui=1
  for(ui in 1:length(lg.sizes))
  {
    
    extent.dst=points.to.calculate*lg.sizes[ui]
    dist.points=seq(-extent.dst,extent.dst,lg.sizes[ui])
    
    single.square=(dist.points[2]-dist.points[1])^2
    single.square
    
    if(middle==F)
    {
      mat.x.app=matrix(data=dist.points,
                       nrow=length(dist.points),ncol=length(dist.points))
      
      mat.y.app=matrix(data=dist.points,
                       nrow=length(dist.points),ncol=length(dist.points),byrow=T)
      
    }else{
      mat.x.app=matrix(data=dist.points,
                       nrow=length(dist.points),ncol=length(dist.points))+(lg.sizes[ui]/2)*x.sc
      
      mat.y.app=matrix(data=dist.points,
                       nrow=length(dist.points),ncol=length(dist.points),byrow=T)+(lg.sizes[ui]/2)*y.sc
      
    }
    
    approx.det.dist=sqrt(mat.x.app^2+mat.y.app^2+approx.det.h^2)
    
    approx.total.fluence= sum(((single.square)*   exp(-0.0103*approx.det.dist))  /(4*pi*approx.det.dist^2) )
    
    tot.fluence.vector[ui]=approx.total.fluence
    
    
  }    
  
  
  return(data.frame(lg.sizes,tot.fluence.vector))
  
}




############################ global parameters #####################
############################ global parameters #####################


# sim=100000
sim=30000



sigma_MH = 0.1
sigma.reduction.factor=0.8
# sigma_MH = 1

infinite.extent=10


mu.soil=0.00813*1520

mu.air=0.00949484


deposition.depth=0.05



ngamma=1

scale.factor=10


detector.height=15



    
    
  


    ### INITIALIZATION OF DATA FROM WILLIAM ###
    
   
      col.pal=brewer.pal(9,"YlGnBu")
      col.pal=rev(colorRampPalette(col.pal)(1000))
    

      

    
    
    
    ########## position matrix generation ################
    
area.matrix.bounds=c(-100,100)
      measurement.matrix.bounds=c(-100,100)
      
      
      area.matrix.size=c(scale.factor,scale.factor)

      measurement.matrix.size=c(scale.factor,scale.factor)
      scale.factor
      # measurement.matrix.bounds=c(-100,100)
      
      
      
      ### Surface matrix
      area.matrix.x=matrix(data=seq(area.matrix.bounds[1],area.matrix.bounds[2],length.out=area.matrix.size[1]),nrow=area.matrix.size[1],ncol=area.matrix.size[2],byrow = T)
      area.matrix.y=matrix(data=seq(area.matrix.bounds[1],area.matrix.bounds[2],length.out=area.matrix.size[2]),nrow=area.matrix.size[2],ncol=area.matrix.size[1],byrow = F)
      area.matrix.z=matrix(data=0,nrow=area.matrix.size[2],ncol=area.matrix.size[1],byrow = F)
      area.matrix=list(area.matrix.x,area.matrix.y,area.matrix.z)
      area.matrix
      
      
      
      measurement.matrix.x=matrix(data=seq(measurement.matrix.bounds[1],measurement.matrix.bounds[2],length.out=measurement.matrix.size[1]),nrow=measurement.matrix.size[2],ncol=measurement.matrix.size[1],byrow = T)
      measurement.matrix.y=matrix(data=seq(measurement.matrix.bounds[1],measurement.matrix.bounds[2],length.out=measurement.matrix.size[2]),nrow=measurement.matrix.size[2],ncol=measurement.matrix.size[1],byrow = F)
      measurement.matrix.z=matrix(data=detector.height,nrow=measurement.matrix.size[2],ncol=measurement.matrix.size[1],byrow = F)
      measurement.matrix=list(measurement.matrix.x,measurement.matrix.y,measurement.matrix.z)
      measurement.matrix
      

    
   
    
    ### Defining defaults in case of no specification ###
    area.list=area.matrix
    measurement.list=measurement.matrix
    
    
    # measurement.matrix.size[2]
    list.size=area.matrix.size[1]*area.matrix.size[2]
    
    area.size=area.matrix.size[1]*area.matrix.size[2]
    # area.size
    measurement.size=measurement.matrix.size[1]*measurement.matrix.size[2]
    # measurement.size
    
    
    
    
    geom.eff=matrix(nrow=measurement.size,ncol=measurement.size)
    # geom.eff
    
    ##################### Establishing activity conversion ratio ########################
    
    # total.side.length=abs(area.matrix.bounds[1]-area.matrix.bounds[2])
    # total.area=total.side.length^2
    

    activity.normalization=((area.matrix.bounds[2]-area.matrix.bounds[1])/scale.factor)^2
    source.area=((area.matrix.bounds[2]-area.matrix.bounds[1])/scale.factor)^2
    
    
    



    
    ####################### Fluence correction #########################
    detector.height
    fluence.grid=fCalcLg(detector.height)
    fluence.grid[,2]=fluence.grid[,2]/fluence.grid[1,2]
    fluence.grid
    lg.id=abs(round(measurement.matrix.y[1]-measurement.matrix.y[2],0))
    
    fluence.id=which(fluence.grid[,1]==lg.id)
    
    fluence.correction=fluence.grid[fluence.id,2]
    fluence.correction
    
    ######################### Calculating GEOM EFF ################################
    cat("\nEvaluating detector response function...")
    distance=list()
    measurement.list
    # deposition.depth=0
    i=1
    for(i in 1:measurement.size)
    {

      dst=sqrt(
        (area.list[[1]]-measurement.list[[1]][i])^2+
          
          (area.list[[2]]-measurement.list[[2]][i])^2+
          
          (area.list[[3]]-measurement.list[[3]][i])^2)
      
      
      dst.2d=sqrt(
        (area.list[[1]]-measurement.list[[1]][i])^2+
          
          (area.list[[2]]-measurement.list[[2]][i])^2)
      

      ######## Deposition Depth #####################

      l=sqrt( deposition.depth^2+  (  ( deposition.depth*as.vector(dst.2d)  )/(detector.height+deposition.depth)       )^2            )
      
      
      
      
      geom.eff[i,]=fluence.correction*1e3*(   (exp(-mu.air*(as.vector(dst)-as.vector(l) )     )*exp(-mu.soil*as.vector(l)))/
                                                (4*pi*as.vector(dst^2)))
      
      
    }
    
    

    

    
    
    ########################## FUNCTIONS ######################
    

    fCalcInfiniteResponse=function(extent,area.matrix.size,area.matrix.bounds,height,activity=1,ngamma=1,efficiency=1)
    {
      if(extent<1)
      {return(-Inf)}
      
      
      matrix.step.size=abs((area.matrix.bounds[1]-area.matrix.bounds[2]))/(area.matrix.size[1]-1)
      cat("\nNumber of points for infinite plane extent: ", extent,"")
      cat("\nInfinite plane extent to be calculated: ", round(extent*matrix.step.size)," m")
      
      coordinate.vector=seq(area.matrix.bounds[1]-extent*matrix.step.size,
                            area.matrix.bounds[2]+extent*matrix.step.size,
                            matrix.step.size)
      
      meas.coordinate.vector=seq(area.matrix.bounds[1],
                                 area.matrix.bounds[2],
                                 matrix.step.size)


      # activity.z=matrix(height,nrow=2*(limits+exclusion)+1,ncol=2*(limits+exclusion)+1)
      activity.x=matrix(data=coordinate.vector,nrow=length(coordinate.vector),ncol=length(coordinate.vector),byrow=T)

      
      activity.y=matrix(data=coordinate.vector,nrow=length(coordinate.vector),ncol=length(coordinate.vector),byrow=F)

      
      activity.z=matrix(data=height,nrow=length(coordinate.vector),ncol=length(coordinate.vector),byrow=F)

      
      activity.matrix=matrix(data=1,nrow=length(coordinate.vector),ncol=length(coordinate.vector),byrow=F)

      
      activity.matrix[which(activity.x>=area.matrix.bounds[1]&activity.x<=area.matrix.bounds[2]&activity.y>=area.matrix.bounds[1]&activity.y<=area.matrix.bounds[2])]=0

      
      to.return=matrix(nrow=area.matrix.size[1],ncol=area.matrix.size[2])


      for(i in 1:area.matrix.size[1])
      {

        for(j in 1:area.matrix.size[2])
        {

          
          distances=sqrt(
            (activity.x-meas.coordinate.vector[i])^2+
              (activity.y-meas.coordinate.vector[j])^2+
              (activity.z)^2
          )
          
          distances.2d=sqrt(
            (activity.x-meas.coordinate.vector[i])^2+
              (activity.y-meas.coordinate.vector[j])^2
          )

          l=sqrt( deposition.depth^2+  (  ( deposition.depth*distances.2d  )/(detector.height+deposition.depth)       )^2            )
          

          
          to.return[j,i]=fluence.correction*1e3*sum ( activity.matrix* as.vector   (exp(-mu.air*(distances-l))  *exp(-mu.soil*l))          /(4*pi*distances^2)  ) 
          
          
        }
        
        
      }

      return(to.return)
    }
    

    
    infinite.response=fCalcInfiniteResponse(100,area.matrix.size,area.matrix.bounds,height=detector.height)
    infinite.response
    

    fLikelihood=function(detector.data,geom.eff,lambda)
    {
      
      if(min(lambda)<0)
        return(-Inf)
      
      
      return.likelihood=sum(dpois(detector.data, (geom.eff%*% lambda)+mean(lambda)*as.vector(infinite.response),log=T))
      
      
      return(return.likelihood)
    }
    
    

    
    
    
    ################## OTHER CODE ##############################
    
    
    
  
    ################## MCMC PART ##############################
    ### Starting parameters ##########
    
    

      
      
      accept <- 0
      count <- 0
      i=1
      count=0
      


      

      beta_sigma = 10^-4
      sim
      
      number_of_parameters=  area.matrix.size[1]*area.matrix.size[2]
      Sigma = diag(  number_of_parameters)*10^2
      Sigma2 = diag(  number_of_parameters)*10^2*0.5
      mean_est   = matrix(0, nrow=1, ncol = number_of_parameters  )
      Sigma_est  = matrix(0, nrow=number_of_parameters, ncol=   number_of_parameters)
      lambda_vec = matrix(0, nrow= sim, ncol =   number_of_parameters)
      # lambda_star_vec=matrix(0, nrow= sim, ncol =   number_of_parameters)
      sigma_mh_vec=matrix(0,nrow=sim,ncol=  number_of_parameters)
      
      
      # surface.activity.matrix%*%geom.eff.matrix
      # image(surface.activity.matrix,col=col.pal)
      # surf.act=as.vector(activity.matrix)
      # surf.actfCal
      
      # fCalcResponseMatrixOptim(1)
      
      ################## Generate data ###################
    
      activity.matrix=matrix(100,nrow=area.matrix.size[1],ncol=area.matrix.size[2])
      activity.matrix[35]=350
      activity.matrix=activity.matrix*activity.normalization
   scale.factor
      # geom.eff%*%as.vector(activity.matrix)
      mcmc.data=round( matrix(      geom.eff%*%as.vector(activity.matrix),nrow=scale.factor,ncol=scale.factor)+mean(activity.matrix)*infinite.response ,0)


      
      

  
  
      lik.vals=c()
      # stop()
      it.accept=0
      
      cat("\nStarting MCMC...")

      time.segment.start=Sys.time()
      
      
      i=1    
      k=1
      
      lambda=rep(100*activity.normalization,scale.factor^2)
      
      count=0
      
      for(i in 1:sim)
      {
        
        count  = count + 1
        
        lambda_star = lambda
        
        lambda_star[1:number_of_parameters] = rmvnorm(1,mean  = lambda[1:number_of_parameters], sigma = sigma_MH^2 * Sigma, method="chol")
        lambda_star
        U = runif(1)
        acc = min(1,  exp( fLikelihood(mcmc.data,geom.eff,lambda=lambda_star) - fLikelihood(mcmc.data,geom.eff,lambda=lambda)))
        acc
        
        if(U < acc )
        {
          lambda = lambda_star
          accept = accept + 1 #counting number of acceptance
          it.accept=it.accept+1
        }
        lambda_vec[i,] = lambda
        
        
        ################# ADAPTIVE MCMC ####################
        
        #adaptive algoritm
        w = 1/i
        mean_est = (1-w) *  + w*lambda[1:number_of_parameters,drop=F]
        
        Sigma_est = (1-w) * Sigma_est + w *(lambda[1:number_of_parameters,drop=F] - mean_est)%*% t(lambda[1:number_of_parameters,drop=F] - mean_est)
        if(i > 100)
        {
          Sigma = diag(diag(Sigma_est)) + diag(length(mean_est)) * beta_sigma
        }
        
        
        
        if(i%%100 == 0)
        {
          if(accept/count < 0.24 )
            sigma_MH = max(0.95,(1-5/sqrt(i))) * sigma_MH
          
          else
            sigma_MH = min(1.05, 1+5/sqrt(i)) * sigma_MH
          
          accept = 0
          count  = 0
          
          
          
          
        }
        
        
        if(i%%1000 == 0)
        {
          
          # cat("loooooooooooooooooooooooooooooooooooooooooooooooool")
          
          cat("\n",i/sim*100," % done...")
          cat("\nAcceptance: ",it.accept/1000)
          

          
          
          
          
          
          
          
          
          par(mfrow=c(2,1))
          plot(lambda_vec[seq(0,sim,100),1],xlab="MCMC iteration",ylim=c(min(lambda_vec),max(lambda_vec)),type="l")
          for(i in 2:(length(lambda_vec[1,])))
          {
            points(lambda_vec[seq(0,sim,100),i],col=i,type="l")
          }
          
          plot(lik.vals[max(1,(length(lik.vals)-50)):(length(lik.vals))],type="l",ylab="",main="Likelihood of latest 500 iterations")
          
        }
        
        
        if(i%%100 == 0)
        {

          lik.vals[k]=fLikelihood(mcmc.data,geom.eff,lambda=lambda)

          k=k+1

        }
        
        
        
        
    
    
        
  }
      
      
      time.segment.end=Sys.time()

      time.difference=as.numeric(difftime(time.segment.end,time.segment.start,units="secs"))

      

      cat("\n\nMCMC time: ",time.difference," seconds, acceptance: ",accept/sim)
      

      map.matrix=matrix(nrow=scale.factor,ncol=scale.factor)
      for(i in 1:length(lambda_vec[1,]))
      {

        
        ev.density=density(lambda_vec[,i])
        MAP.value=round( ev.density$x[  which(ev.density$y==max(ev.density$y))       ] )

        map.matrix[i]=MAP.value
        
        
      }
      

      
    
    ########################### Fluence ########################
    {
      # plotname=paste0("~/Dropbox/Bayesiansk statistikprojekt/R scripts/Project_folder/JW/Article_4/FIGS/FLU_h",detector.height,"_c",cleanup.scenario,".png")
      # png(filename = plotname,width = 1000,height=1000)
      par(mfrow=c(1,1))
      # par(mfrow=c(1,2))
      par(mar=c(2.3,2.3,4.5,1),mgp=c(1.3,0.45,0),cex=1,cex.axis=1,cex.lab=1,cex.main=1.2,cex.sub=1,family="serif",lwd=2)
      
      
      
      
      par(mfrow=c(1,1),mar=c(3,3,2,1))
      layout(matrix(c(1,2),nrow=1),widths=c(0.8,0.2))
      par(mar=c(2.3,4.3,3.5,1),mgp=c(1.3,0.45,0),cex=1,cex.axis=1,cex.lab=1,cex.main=1.2,cex.sub=1,family="serif",lwd=2)
      
      image(activity.matrix,col=col.pal,axes=F,xlab="x position (m)",ylab="y position (m)")
      # title(main=paste0("Detector height ",detector.height,", cleanup scenario ",cleanup.scenario))
      box()
      axis(1,at=seq(0,1,0.25),labels=seq(0,1,0.25)*140-70)
      axis(2,at=seq(0,1,0.25),labels=seq(0,1,0.25)*140-70)
      
      par(mar=c(2.3,0,3.5,4))
      plot.new()
      plot.window(xlim=c(-0.1,0.1),ylim=c(0,1),xaxs="i",yaxs="i")
      # box()
      n_of_p=1000
      for(i in 1:100)
      {
        points(x=rep(-0.5+(i/100),n_of_p),y=seq(0,1,length.out=n_of_p),col=col.pal,pch=15,cex=0.5)   
      }
      
      # points(x=rep(0,n_of_p),y=seq(0,1,length.out=n_of_p),col=col.pal,pch=15,cex=0.2) 
      axis(4,at=seq(0,1,length.out=10),las=1,labels=round(seq(min(activity.matrix),max(activity.matrix),length.out=10)/activity.normalization),cex.axis=0.7)
      # mtext("Fluence (1/m^2)",side=4,line=2.2)
      mtext(expression("Activity (kBq/m"^2*")"),side=4,line=2.2,cex=1)
      
      mtext(paste0("Simulated activity"), side = 3, line = -2, outer = TRUE,cex=1,font=2)
      
      
      
      # dev.off()
    }
    
    # col.pal=hcl.colors(1000)
    ########################### Reconstruction ########################
    {
      
      # zlims=c(0,500*1)
      
      # plotname=paste0("~/Dropbox/Bayesiansk statistikprojekt/R scripts/Project_folder/JW/Article_4/FIGS/REC_h",detector.height,"_c",cleanup.scenario,".png")
      # png(filename = plotname,width = 1000,height=1000)
      par(mfrow=c(1,1))
      # par(mfrow=c(1,2))
      par(mar=c(2.3,2.3,4.5,1),mgp=c(1.3,0.45,0),cex=1,cex.axis=1,cex.lab=1,cex.main=1.2,cex.sub=1,family="serif",lwd=2)
      # col.pal=RColorBrewer::brewer.pal(9,"Blues")
      # col.pal = rev(colorRampPalette(col.pal)(1000))
      # col.pal=hcl.colors(1000,"YlGnBu")
      
      
      # par(mfrow=c(1,1),mar=c(3,3,2,1))
      layout(matrix(c(1,2),nrow=1),widths=c(0.8,0.2))
      par(mar=c(2.3,2.3,3.5,1),mgp=c(1.3,0.45,0),cex=1,cex.axis=1,cex.lab=1,cex.main=1.2,cex.sub=1,family="serif",lwd=2)
      
      image(map.matrix,col=col.pal,axes=F,xlab="x position (m)",ylab="y position (m)")
      # title(main=paste0("Reconstructed activity distribution, h=",detector.height,", cl=",cleanup.scenario))
      box()
      axis(1,at=seq(0,1,0.25),labels=seq(0,1,0.25)*140-70)
      axis(2,at=seq(0,1,0.25),labels=seq(0,1,0.25)*140-70)
      
      par(mar=c(2.3,0,3.5,4))
      plot.new()
      plot.window(xlim=c(-0.1,0.1),ylim=c(0,1),xaxs="i",yaxs="i")
      box()
      n_of_p=1000
      for(i in 1:100)
      {
        points(x=rep(-0.5+(i/100),n_of_p),y=seq(0,1,length.out=n_of_p),col=col.pal,pch=15,cex=0.5)   
      }
      
      # points(x=rep(0,n_of_p),y=seq(0,1,length.out=n_of_p),col=col.pal,pch=15,cex=0.2) 
      axis(4,at=seq(0,1,length.out=10),las=1,labels=round(seq(min(map.matrix)[1],max(map.matrix)[1],length.out=10)/activity.normalization)  )
      # axis(4,at=seq(0,1,length.out=10),las=1,labels=round(seq(min(zlims),max(zlims),length.out=10)))
      
      
      mtext(expression("Activity (kBq/m"^2*")"),side=4,line=2.2,cex=1)
      
      mtext(paste0("Reconstructed activity distribution, h=",detector.height), side = 3, line = -2, outer = TRUE,cex=1,font=2)
      
      # dev.off()
    }
    

