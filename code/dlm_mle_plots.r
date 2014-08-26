# Set data and graphics path
dpath = "../data/"
gpath = "../graphs/"

dlm.mle.plots <- function(nt, nsims, mod, x, beta0, beta1, phi, rho, sigma2s, sigma2b, sigma2m)
{
  # only first element of rho, sigma2b, and sigma2m used
  rho = rho[1]; sigma2b = sigma2b[1]; sigma2m = sigma2m[1]
  
  # Length of time, x, and mod
  Nt = length(nt)
  Nx = length(x)
  Nm = length(mod)
  
  # Plot time/SNR/phi-rho ratio versus parameter for every mod and x
  for(i in 1:Nx)
  {
    for(k in 1:Nm)
    { 
      # Length of SNR and phi-rho ratio
      if(mod[k] == "M111")
      {
        snr = sigma2s / sigma2b
        prr = phi / rho
      } else {
        snr = sigma2s / sigma2m
        prr = phi
      }
      Ns = length(snr)
      Nr = length(prr)
      
      # Time by parameter plots
      for(l in 1:Ns)
      {
        for(m in 1:Nr)
        {
          pdf(file=paste(gpath,paste(nsims,"TR",mod[k],x[i],beta0,beta1,phi[m],rho,sigma2s[l],sigma2b,sigma2m,sep="-"),".pdf",sep=""),width=ifelse(mod[k] == "M111",20,15),height=5*Nt)
          par(mfrow=c(Nt,3+(mod[k]=="M111")),mar=c(7,7,6,2)+0.1,mgp=c(4,1,0))
          for(j in 1:Nt)
          {
            # Load MLEs and find error/convergence rate
            load(paste(dpath,paste(nt[j],nsims,mod[k],x[i],beta0,beta1,phi[m],rho,sigma2s[l],sigma2b,sigma2m,sep="-"),".rdata",sep=""))
            error.rate = mean(out.mle$error)
            conv.rate = mean(out.mle$converge[which(out.mle$error == 0)] == 0)
            ind = which(out.mle$error == 0 & out.mle$converge == 0)
            
            # beta
            require(KernSmooth)
            b = 1.06*c(sd(out.mle$theta.mle[ind,1]),sd(out.mle$theta.mle[ind,2]))*nsims^(-.2)
            tmp = bkde2D(cbind(out.mle$theta.mle[ind,1],out.mle$theta.mle[ind,2]), b)
            if(j == 1) xlab=expression(beta[0]) else xlab=""
            image(tmp$x1,tmp$x2,-tmp$fhat,xlab=xlab,ylab=paste(nt[j],"TRs"),cex.lab=2.5,cex.axis=1.25)
            contour(tmp$x1,tmp$x2,tmp$fhat,add=TRUE,drawlabels=F)
            points(beta0,beta1,pch='x',col=4,lwd=3,cex=2.5)
            abline(h = 0, col='gray',lwd=3)
            if(j == 1)
            {
              mtext(expression(beta[1]),side=2,cex=1.7,line=2)
              title(eval(bquote(expression(M[.(paste(strsplit(mod[k],"")[[1]][2:4],sep="",collapse=""))]))),cex.main=2,line=3)
              mtext(eval(bquote(expression(paste(beta,"= (",.(beta0),",",.(beta1),")")))),side=3,cex=1.2)
              mtext(paste("Error rate: ",round(error.rate,3),", Conv. rate: ",round(conv.rate,3),sep=""),side=3,line=-1)
            }
            
            # phi
            if(mod[k] == "M111")
            {
              b = 1.06*c(sd(out.mle$theta.mle[ind,3]),sd(out.mle$theta.mle[ind,4]))*nsims^(-.2)
              tmp = bkde2D(cbind(out.mle$theta.mle[ind,3],out.mle$theta.mle[ind,4]), b)
              if(j == 1)
              {
                xlab=expression(phi)
                ylab=expression(rho)
              } else {
                xlab=""
                ylab=""
              }
              image(tmp$x1,tmp$x2,-tmp$fhat,xlab=xlab,ylab=ylab,cex.lab=2.5,cex.axis=1.25)
              contour(tmp$x1,tmp$x2,tmp$fhat,add=TRUE,drawlabels=F)
              points(phi[m],rho,pch='x',col=4,lwd=3,cex=2.5)
              if(j == 1) mtext(eval(bquote(expression(paste("(",phi,", ",rho,")"," = (",.(phi[m]),", ",.(rho),")",sep="")))),side=3,cex=1.2)        
            } else {
              if(j == 1) xlab=expression(phi) else xlab=""
              hist(out.mle$theta.mle[ind,3],xlab=xlab,ylab="",main="",cex.lab=2.5,cex.axis=1.25)
              abline(v=phi[m],col=4,lwd=3)
              if(j == 1) mtext(eval(bquote(expression(paste(phi," = ",.(phi[m]),sep="")))),side=3,cex=1.2)
            }
            
            # sigma2s
            if(mod[k] == "M111")
            {
              b = 1.06*c(sd(out.mle$theta.mle[ind,5]),sd(out.mle$theta.mle[ind,6]))*nsims^(-.2)
              tmp = bkde2D(cbind(out.mle$theta.mle[ind,5],out.mle$theta.mle[ind,6]), b)
              if(j == 1)
              {
                xlab=expression(sigma[s]^2)
                ylab=expression(sigma[b]^2)
              } else {
                xlab=""
                ylab=""
              }
              image(tmp$x1,tmp$x2,-tmp$fhat,xlab=xlab,ylab=ylab,cex.lab=2.5,cex.axis=1.25)
              contour(tmp$x1,tmp$x2,tmp$fhat,add=TRUE,drawlabels=F)
              points(sigma2s[l],sigma2b,pch='x',col=4,lwd=3,cex=3)
              if(j == 1) mtext(eval(bquote(expression(paste("(",sigma[s]^2,", ",sigma[b]^2,")"," = (",.(sigma2s[l]),", ",.(sigma2b),")",sep="")))),side=3,cex=1.2)        
            } else {
              b = 1.06*c(sd(out.mle$theta.mle[ind,5]),sd(out.mle$theta.mle[ind,7]))*nsims^(-.2)
              tmp = bkde2D(cbind(out.mle$theta.mle[ind,5],out.mle$theta.mle[ind,7]), b)
              if(j == 1)
              {
                xlab=expression(sigma[s]^2)
                ylab=expression(sigma[m]^2)
              } else {
                xlab=""
                ylab=""
              }
              image(tmp$x1,tmp$x2,-tmp$fhat,xlab=xlab,ylab=ylab,cex.lab=2.5,cex.axis=1.25)
              contour(tmp$x1,tmp$x2,tmp$fhat,add=TRUE,drawlabels=F)
              points(sigma2s[l],sigma2m,pch='x',col=4,lwd=3,cex=3)
              if(j == 1) mtext(eval(bquote(expression(paste("(",sigma[s]^2,", ",sigma[m]^2,")"," = (",.(sigma2s[l]),", ",.(sigma2m),")",sep="")))),side=3,cex=1.2)
            }
            
            # sigma2m
            if(mod[k] == "M111")
            {
              if(j == 1) xlab=expression(sigma[m]^2) else xlab=""
              hist(out.mle$theta.mle[ind,7],xlab=xlab,ylab="",main="",cex.lab=2.5,cex.axis=1.25)
              abline(v=sigma2m,col=4,lwd=3)
              if(j == 1) mtext(eval(bquote(expression(paste(sigma[m]^2," = ",.(sigma2m),sep="")))),side=3,cex=1.2)
            }
          }
          dev.off()
        }
      }

      # SNR by parameter plots
      for(m in 1:Nr)
      {
        for(j in 1:Nt)
        {
          # Load MLEs and find maximum and minimums
          ymax.beta0 = ymax.beta1 = ymax.phi = ymax.rho = ymax.sigma2m.x = ymax.sigma2m.y = -Inf
          ymin.beta0 = ymin.beta1 = ymin.phi = ymin.rho = ymin.sigma2m.x = ymin.sigma2m.y = Inf
          for(l in 1:Ns)
          {
            load(paste(dpath,paste(nt[j],nsims,mod[k],x[i],beta0,beta1,phi[m],rho,sigma2s[l],sigma2b,sigma2m,sep="-"),".rdata",sep=""))          
            ind = which(out.mle$error == 0 & out.mle$converge == 0)
            require(KernSmooth)
            b = 1.06*c(sd(out.mle$theta.mle[ind,1]),sd(out.mle$theta.mle[ind,2]))*nsims^(-.2)
            tmp = bkde2D(cbind(out.mle$theta.mle[ind,1],out.mle$theta.mle[ind,2]), b)
            ymax.beta0 = max(ymax.beta0,tmp$x1); ymax.beta1 = max(ymax.beta1,tmp$x2)
            ymin.beta0 = min(ymin.beta0,tmp$x1); ymin.beta1 = min(ymin.beta1,tmp$x2)
            if(mod[k] == "M111")
            {
              b = 1.06*c(sd(out.mle$theta.mle[ind,3]),sd(out.mle$theta.mle[ind,4]))*nsims^(-.2)
              tmp = bkde2D(cbind(out.mle$theta.mle[ind,3],out.mle$theta.mle[ind,4]), b)
              ymax.phi = max(ymax.phi,tmp$x1); ymax.rho = max(ymax.rho,tmp$x2)
              ymin.phi = min(ymin.phi,tmp$x1); ymin.rho = min(ymin.rho,tmp$x2)
              tmp = hist(out.mle$theta.mle[ind,7],plot=FALSE)
              ymax.sigma2m.x = max(ymax.sigma2m.x,tmp$breaks); ymin.sigma2m.x = min(ymin.sigma2m.x,tmp$breaks)
              ymax.sigma2m.y = max(ymax.sigma2m.y,tmp$counts); ymin.sigma2m.y = min(ymin.sigma2m.y,tmp$counts)
            } else {
              tmp = hist(out.mle$theta.mle[ind,3],plot=FALSE)
              ymax.phi = max(ymax.phi,tmp$breaks); ymin.phi = min(ymin.phi,tmp$breaks)
              ymax.rho = max(ymax.rho,tmp$counts); ymin.rho = min(ymin.rho,tmp$counts)
            }
          }

          pdf(file=paste(gpath,paste(nsims,nt[j],mod[k],x[i],beta0,beta1,phi[m],rho,"SNR",sigma2b,sigma2m,sep="-"),".pdf",sep=""),width=ifelse(mod[k] == "M111",20,15),height=5*Ns)
          par(mfrow=c(Ns,3+(mod[k]=="M111")),mar=c(7,12,9,2)+0.1,mgp=c(6.5,1,0))
          for(l in 1:Ns)
          {
            # Load MLEs and find error/convergence rate
            load(paste(dpath,paste(nt[j],nsims,mod[k],x[i],beta0,beta1,phi[m],rho,sigma2s[l],sigma2b,sigma2m,sep="-"),".rdata",sep=""))
            error.rate = mean(out.mle$error)
            conv.rate = mean(out.mle$converge[which(out.mle$error == 0)] == 0)
            ind = which(out.mle$error == 0 & out.mle$converge == 0)
            
            # SNR label
            if(mod[k] == "M111") snr.lab = bquote(expression(paste(sigma[s]^2/sigma[b]^2,"=",.(snr[l])))) else snr.lab = bquote(expression(paste(sigma[s]^2/sigma[m]^2,"=",.(snr[l])))) 
            
            # beta
            require(KernSmooth)
            b = 1.06*c(sd(out.mle$theta.mle[ind,1]),sd(out.mle$theta.mle[ind,2]))*nsims^(-.2)
            tmp = bkde2D(cbind(out.mle$theta.mle[ind,1],out.mle$theta.mle[ind,2]), b)
            if(l == 1) xlab=expression(beta[0]) else xlab=""
            if(l == 1) axes = TRUE else axes = FALSE
            plot(tmp$x1,tmp$x2,axes=axes,xlim=c(ymin.beta0,ymax.beta0),ylim=c(ymin.beta1,ymax.beta1),xlab=xlab,ylab=eval(snr.lab),cex.lab=4,cex.axis=1.9)
            u = par("usr")
            rect(u[1], u[3], u[2], u[4], col = heat.colors(12)[12])
            image(tmp$x1,tmp$x2,-tmp$fhat,add=TRUE)
            contour(tmp$x1,tmp$x2,tmp$fhat,add=TRUE,drawlabels=F)
            points(beta0,beta1,pch='x',col=4,lwd=3,cex=3.5)
            abline(h = 0, col='gray',lwd=3)
            if(l == 1)
            {
              mtext(expression(beta[1]),side=2,cex=3,line=2.5)
              title(eval(bquote(expression(M[.(paste(strsplit(mod[k],"")[[1]][2:4],sep="",collapse=""))]))),cex.main=4,line=5)
              mtext(eval(bquote(expression(paste(beta,"= (",.(beta0),",",.(beta1),")")))),side=3,cex=2)
            }
            mtext(paste("Conv. rate: ",round(conv.rate,3),sep=""),side=3,line=-1.5,cex=1.4)
            box()
            #if(error.rate != 0 | conv.rate != 1) print(c(round(error.rate,3),round(conv.rate,3),paste(nsims,nt[j],mod[k],x[i],beta0,beta1,phi[m],rho,"SNR",sigma2b,sigma2m,sep="-")))
            
            
            # phi
            if(mod[k] == "M111")
            {
              b = 1.06*c(sd(out.mle$theta.mle[ind,3]),sd(out.mle$theta.mle[ind,4]))*nsims^(-.2)
              tmp = bkde2D(cbind(out.mle$theta.mle[ind,3],out.mle$theta.mle[ind,4]), b)
              if(l == 1)
              {
                xlab=expression(phi)
                ylab=expression(rho)
                axes=TRUE
              } else {
                xlab=""
                ylab=""
                axes=FALSE
              }
              plot(tmp$x1,tmp$x2,axes=axes,xlim=c(ymin.phi,ymax.phi),ylim=c(ymin.rho,ymax.rho),xlab=xlab,ylab=ylab,cex.lab=4,cex.axis=1.9)
              u = par("usr")
              rect(u[1], u[3], u[2], u[4], col = heat.colors(12)[12])              
              image(tmp$x1,tmp$x2,-tmp$fhat,add=TRUE)
              contour(tmp$x1,tmp$x2,tmp$fhat,add=TRUE,drawlabels=F)
              points(phi[m],rho,pch='x',col=4,lwd=3,cex=3.5)
              if(l == 1) mtext(eval(bquote(expression(paste("(",phi,", ",rho,")"," = (",.(phi[m]),", ",.(rho),")",sep="")))),side=3,cex=2)        
            } else {
              if(l == 1) xlab=expression(phi) else xlab=""
              hist(out.mle$theta.mle[ind,3],axes=(l==1),xlim=c(ymin.phi,ymax.phi),ylim=c(ymin.rho,ymax.rho),xlab=xlab,ylab="",main="",cex.lab=4,cex.axis=1.9)
              abline(v=phi[m],col=4,lwd=4)
              if(l == 1) mtext(eval(bquote(expression(paste(phi," = ",.(phi[m]),sep="")))),side=3,cex=2)
            }
            box()
            
            # sigma2s
            if(mod[k] == "M111")
            {
              b = 1.06*c(sd(out.mle$theta.mle[ind,5]),sd(out.mle$theta.mle[ind,6]))*nsims^(-.2)
              tmp = bkde2D(cbind(out.mle$theta.mle[ind,5],out.mle$theta.mle[ind,6]), b)
              if(l == 1)
              {
                xlab=expression(sigma[s]^2)
                ylab=expression(sigma[b]^2)
              } else {
                xlab=""
                ylab=""
              }
              image(tmp$x1,tmp$x2,-tmp$fhat,xlab=xlab,ylab=ylab,cex.lab=4,cex.axis=1.9)
              contour(tmp$x1,tmp$x2,tmp$fhat,add=TRUE,drawlabels=F)
              points(sigma2s[l],sigma2b,pch='x',col=4,lwd=3,cex=3.5)
              mtext(eval(bquote(expression(paste("(",sigma[s]^2,", ",sigma[b]^2,")"," = (",.(sigma2s[l]),", ",.(sigma2b),")",sep="")))),side=3,cex=2)        
            } else {
              b = 1.06*c(sd(out.mle$theta.mle[ind,5]),sd(out.mle$theta.mle[ind,7]))*nsims^(-.2)
              tmp = bkde2D(cbind(out.mle$theta.mle[ind,5],out.mle$theta.mle[ind,7]), b)
              if(l == 1)
              {
                xlab=expression(sigma[s]^2)
                ylab=expression(sigma[m]^2)
              } else {
                xlab=""
                ylab=""
              }
              image(tmp$x1,tmp$x2,-tmp$fhat,xlab=xlab,ylab=ylab,cex.lab=4,cex.axis=1.9)
              contour(tmp$x1,tmp$x2,tmp$fhat,add=TRUE,drawlabels=F)
              points(sigma2s[l],sigma2m,pch='x',col=4,lwd=3,cex=3.5)
              mtext(eval(bquote(expression(paste("(",sigma[s]^2,", ",sigma[m]^2,")"," = (",.(sigma2s[l]),", ",.(sigma2m),")",sep="")))),side=3,cex=2)
            }
            box()
            
            # sigma2m
            if(mod[k] == "M111")
            {
              if(l == 1) xlab=expression(sigma[m]^2) else xlab=""
              hist(out.mle$theta.mle[ind,7],axes=(l==1),xlim=c(ymin.sigma2m.x,ymax.sigma2m.x),ylim=c(ymin.sigma2m.y,ymax.sigma2m.y),xlab=xlab,ylab="",main="",cex.lab=4,cex.axis=1.9)
              abline(v=sigma2m,col=4,lwd=4)
              if(l == 1) mtext(eval(bquote(expression(paste(sigma[m]^2," = ",.(sigma2m),sep="")))),side=3,cex=2)
              box()
            }
          }
          dev.off()
        }
      }
      
      # phi-rho ratio by parameter plots
      for(l in 1:Ns)
      {
        for(j in 1:Nt)
        {
          # Load MLEs and find maximum and minimums
          ymax.beta0 = ymax.beta1 = ymax.s2s = ymax.s2b = ymax.sigma2m.x = ymax.sigma2m.y = -Inf
          ymin.beta0 = ymin.beta1 = ymin.s2s = ymin.s2b = ymin.sigma2m.x = ymin.sigma2m.y = Inf
          for(m in 1:Nr)
          {
            load(paste(dpath,paste(nt[j],nsims,mod[k],x[i],beta0,beta1,phi[m],rho,sigma2s[l],sigma2b,sigma2m,sep="-"),".rdata",sep=""))
            ind = which(out.mle$error == 0 & out.mle$converge == 0)
            require(KernSmooth)
            b = 1.06*c(sd(out.mle$theta.mle[ind,1]),sd(out.mle$theta.mle[ind,2]))*nsims^(-.2)
            tmp = bkde2D(cbind(out.mle$theta.mle[ind,1],out.mle$theta.mle[ind,2]), b)
            ymax.beta0 = max(ymax.beta0,tmp$x1); ymax.beta1 = max(ymax.beta1,tmp$x2)
            ymin.beta0 = min(ymin.beta0,tmp$x1); ymin.beta1 = min(ymin.beta1,tmp$x2)
            if(mod[k] == "M111")
            {
              b = 1.06*c(sd(out.mle$theta.mle[ind,5]),sd(out.mle$theta.mle[ind,6]))*nsims^(-.2)
              tmp = bkde2D(cbind(out.mle$theta.mle[ind,5],out.mle$theta.mle[ind,6]), b)
              ymax.s2s = max(ymax.s2s,tmp$x1); ymax.s2b = max(ymax.s2b,tmp$x2)
              ymin.s2s = min(ymin.s2s,tmp$x1); ymin.s2b = min(ymin.s2b,tmp$x2)
              tmp = hist(out.mle$theta.mle[ind,7],plot=FALSE)
              ymax.sigma2m.x = max(ymax.sigma2m.x,tmp$breaks); ymin.sigma2m.x = min(ymin.sigma2m.x,tmp$breaks)
              ymax.sigma2m.y = max(ymax.sigma2m.y,tmp$counts); ymin.sigma2m.y = min(ymin.sigma2m.y,tmp$counts)
            } else {
              b = 1.06*c(sd(out.mle$theta.mle[ind,5]),sd(out.mle$theta.mle[ind,7]))*nsims^(-.2)
              tmp = bkde2D(cbind(out.mle$theta.mle[ind,5],out.mle$theta.mle[ind,7]), b)
              ymax.s2s = max(ymax.s2s,tmp$x1); ymax.s2b = max(ymax.s2b,tmp$x2)
              ymin.s2s = min(ymin.s2s,tmp$x1); ymin.s2b = min(ymin.s2b,tmp$x2)
            }
          }
          
          pdf(file=paste(gpath,paste(nsims,nt[j],mod[k],x[i],beta0,beta1,"PRR",rho,sigma2s[l],sigma2b,sigma2m,sep="-"),".pdf",sep=""),width=ifelse(mod[k] == "M111",20,15),height=5*Nr)
          par(mfrow=c(Nr,3+(mod[k]=="M111")),mar=c(7,12,9,2)+0.1,mgp=c(6.5,1,0))
          for(m in 1:Nr)
          {
            # Load MLEs and find error/convergence rate
            load(paste(dpath,paste(nt[j],nsims,mod[k],x[i],beta0,beta1,phi[m],rho,sigma2s[l],sigma2b,sigma2m,sep="-"),".rdata",sep=""))
            error.rate = mean(out.mle$error)
            conv.rate = mean(out.mle$converge[which(out.mle$error == 0)] == 0)
            ind = which(out.mle$error == 0 & out.mle$converge == 0)
            
            # SNR label
            if(mod[k] == "M111") prr.lab = bquote(expression(paste(phi/rho,"=",.(prr[m])))) else prr.lab = bquote(expression(paste(phi,"=",.(phi[m]))))
            
            # beta
            require(KernSmooth)
            b = 1.06*c(sd(out.mle$theta.mle[ind,1]),sd(out.mle$theta.mle[ind,2]))*nsims^(-.2)
            tmp = bkde2D(cbind(out.mle$theta.mle[ind,1],out.mle$theta.mle[ind,2]), b)
            if(m == 1) xlab=expression(beta[0]) else xlab=""            
            if(m == 1) axes = TRUE else axes = FALSE
            plot(tmp$x1,tmp$x2,axes=axes,xlim=c(ymin.beta0,ymax.beta0),ylim=c(ymin.beta1,ymax.beta1),xlab=xlab,ylab=eval(prr.lab),cex.lab=4,cex.axis=1.9)
            u = par("usr")
            rect(u[1], u[3], u[2], u[4], col = heat.colors(12)[12])            
            image(tmp$x1,tmp$x2,-tmp$fhat,add=TRUE)
            contour(tmp$x1,tmp$x2,tmp$fhat,add=TRUE,drawlabels=F)
            points(beta0,beta1,pch='x',col=4,lwd=3,cex=3.5)
            abline(h = 0, col='gray',lwd=3)
            if(m == 1)
            {
              mtext(expression(beta[1]),side=2,cex=3,line=2.5)
              title(eval(bquote(expression(M[.(paste(strsplit(mod[k],"")[[1]][2:4],sep="",collapse=""))]))),cex.main=4,line=5)
              mtext(eval(bquote(expression(paste(beta,"= (",.(beta0),",",.(beta1),")")))),side=3,cex=2)
            }
            mtext(paste("Conv. rate: ",round(conv.rate,3),sep=""),side=3,line=-1.5,cex=1.4)
            box()
            #print(c(round(error.rate,3),round(conv.rate,3),paste(nsims,nt[j],mod[k],x[i],beta0,beta1,"PRR",rho,sigma2s[l],sigma2b,sigma2m,sep="-")))
            
            # phi
            if(mod[k] == "M111")
            {
              b = 1.06*c(sd(out.mle$theta.mle[ind,3]),sd(out.mle$theta.mle[ind,4]))*nsims^(-.2)
              tmp = bkde2D(cbind(out.mle$theta.mle[ind,3],out.mle$theta.mle[ind,4]), b)
              if(m == 1)
              {
                xlab=expression(phi)
                ylab=expression(rho)
              } else {
                xlab=""
                ylab=""
              }
              image(tmp$x1,tmp$x2,-tmp$fhat,xlab=xlab,ylab=ylab,cex.lab=4,cex.axis=1.9)
              contour(tmp$x1,tmp$x2,tmp$fhat,add=TRUE,drawlabels=F)
              points(phi[m],rho,pch='x',col=4,lwd=3,cex=3.5)
              mtext(eval(bquote(expression(paste("(",phi,", ",rho,")"," = (",.(phi[m]),", ",.(rho),")",sep="")))),side=3,cex=2)        
            } else {
              if(m == 1) xlab=expression(phi) else xlab=""
              hist(out.mle$theta.mle[ind,3],xlab=xlab,ylab="",main="",cex.lab=4,cex.axis=1.9)
              abline(v=phi[m],col=4,lwd=4)
              mtext(eval(bquote(expression(paste(phi," = ",.(phi[m]),sep="")))),side=3,cex=2)
            }
            
            # sigma2s
            if(mod[k] == "M111")
            {
              b = 1.06*c(sd(out.mle$theta.mle[ind,5]),sd(out.mle$theta.mle[ind,6]))*nsims^(-.2)
              tmp = bkde2D(cbind(out.mle$theta.mle[ind,5],out.mle$theta.mle[ind,6]), b)
              if(m == 1)
              {
                xlab=expression(sigma[s]^2)
                ylab=expression(sigma[b]^2)
                axes=TRUE
              } else {
                xlab=""
                ylab=""
                axes=FALSE
              }
              plot(tmp$x1,tmp$x2,axes=axes,xlim=c(ymin.s2s,ymax.s2s),ylim=c(ymin.s2b,ymax.s2b),xlab=xlab,ylab=ylab,cex.lab=4,cex.axis=1.9)
              u = par("usr")
              rect(u[1], u[3], u[2], u[4], col = heat.colors(12)[12])              
              image(tmp$x1,tmp$x2,-tmp$fhat,add=TRUE)
              contour(tmp$x1,tmp$x2,tmp$fhat,add=TRUE,drawlabels=F)
              points(sigma2s[l],sigma2b,pch='x',col=4,lwd=3,cex=3.5)
              if(m == 1) mtext(eval(bquote(expression(paste("(",sigma[s]^2,", ",sigma[b]^2,")"," = (",.(sigma2s[l]),", ",.(sigma2b),")",sep="")))),side=3,cex=2)        
            } else {
              b = 1.06*c(sd(out.mle$theta.mle[ind,5]),sd(out.mle$theta.mle[ind,7]))*nsims^(-.2)
              tmp = bkde2D(cbind(out.mle$theta.mle[ind,5],out.mle$theta.mle[ind,7]), b)
              if(m == 1)
              {
                xlab=expression(sigma[s]^2)
                ylab=expression(sigma[m]^2)
                axes=TRUE
              } else {
                xlab=""
                ylab=""
                axes=FALSE
              }              
              plot(tmp$x1,tmp$x2,axes=axes,xlim=c(ymin.s2s,ymax.s2s),ylim=c(ymin.s2b,ymax.s2b),xlab=xlab,ylab=ylab,cex.lab=4,cex.axis=1.9)
              u = par("usr")
              rect(u[1], u[3], u[2], u[4], col = heat.colors(12)[12])               
              image(tmp$x1,tmp$x2,-tmp$fhat,add=TRUE)
              contour(tmp$x1,tmp$x2,tmp$fhat,add=TRUE,drawlabels=F)
              points(sigma2s[l],sigma2m,pch='x',col=4,lwd=3,cex=3.5)
              if(m == 1) mtext(eval(bquote(expression(paste("(",sigma[s]^2,", ",sigma[m]^2,")"," = (",.(sigma2s[l]),", ",.(sigma2m),")",sep="")))),side=3,cex=2)
            }
            box()
            
            # sigma2m
            if(mod[k] == "M111")
            {
              if(m == 1) xlab=expression(sigma[m]^2) else xlab=""
              hist(out.mle$theta.mle[ind,7],axes=(m==1),xlim=c(ymin.sigma2m.x,ymax.sigma2m.x),ylim=c(ymin.sigma2m.y,ymax.sigma2m.y),xlab=xlab,ylab="",main="",cex.lab=4,cex.axis=1.9)
              abline(v=sigma2m,col=4,lwd=4)
              if(m == 1) mtext(eval(bquote(expression(paste(sigma[m]^2," = ",.(sigma2m),sep="")))),side=3,cex=2)
            }
          }
          dev.off()
        }
      }      
    }
  }
}

dlm.mle.plots(c(250,500,1000),1000,c("M101","M011"),"conv",750,15,c(.1,.5,.9),.5,c(1,5,10,15,20),10,10)
dlm.mle.plots(c(250,500,1000),1000,c("M101","M011"),"norm",750,15,.5,.5,5,10,10)
dlm.mle.plots(c(250,500,1000),1000,"M101","none",750,15,.5,.5,5,10,10)

dlm.mle.plots(c(250,500,1000),1000,"M111","norm",750,15,.6,.6,10,10,10)

dlm.mle.plots(c(250,500,1000),1000,"M111","conv",750,15,c(.3,.6,.9),.3,c(1,5,10,15,20),1,10)
dlm.mle.plots(c(250,500,1000),1000,"M111","conv",750,15,c(.3,.6,.9),.6,c(1,5,10,15,20),1,10)
dlm.mle.plots(c(250,500,1000),1000,"M111","conv",750,15,c(.3,.6,.9),.9,c(1,5,10,15,20),1,10)
dlm.mle.plots(c(250,500,1000),1000,"M111","conv",750,15,c(.3,.6,.9),.3,c(1,5,10,15,20),5,10)
dlm.mle.plots(c(250,500,1000),1000,"M111","conv",750,15,c(.3,.6,.9),.6,c(1,5,10,15,20),5,10)
dlm.mle.plots(c(250,500,1000),1000,"M111","conv",750,15,c(.3,.6,.9),.9,c(1,5,10,15,20),5,10)
dlm.mle.plots(c(250,500,1000),1000,"M111","conv",750,15,c(.3,.6,.9),.3,c(1,5,10,15,20),10,10)
dlm.mle.plots(c(250,500,1000),1000,"M111","conv",750,15,c(.3,.6,.9),.6,c(1,5,10,15,20),10,10)
dlm.mle.plots(c(250,500,1000),1000,"M111","conv",750,15,c(.3,.6,.9),.9,c(1,5,10,15,20),10,10)
dlm.mle.plots(c(250,500,1000),1000,"M111","conv",750,15,c(.3,.6,.9),.3,c(1,5,10,15,20),15,10)
dlm.mle.plots(c(250,500,1000),1000,"M111","conv",750,15,c(.3,.6,.9),.6,c(1,5,10,15,20),15,10)
dlm.mle.plots(c(250,500,1000),1000,"M111","conv",750,15,c(.3,.6,.9),.9,c(1,5,10,15,20),15,10)
dlm.mle.plots(c(250,500,1000),1000,"M111","conv",750,15,c(.3,.6,.9),.3,c(1,5,10,15,20),20,10)
dlm.mle.plots(c(250,500,1000),1000,"M111","conv",750,15,c(.3,.6,.9),.6,c(1,5,10,15,20),20,10)
dlm.mle.plots(c(250,500,1000),1000,"M111","conv",750,15,c(.3,.6,.9),.9,c(1,5,10,15,20),20,10)


