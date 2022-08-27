require(caper)
require(geiger)

diagnostics.plot<-function(mod.res){
  old.par = par(no.readonly = TRUE)
  par(mfrow=c(2, 2))
  par(mar=c(3, 3, 1, 0.5))
  hist(mod.res$phyres, probability=T, xlab="", ylab="", main="")
  mtext(text="histogram of residuals", side=3, line=0)
  x=seq(min(mod.res$phyres), max(mod.res$phyres), length.out=100)
  lines(x, dnorm(x, mean=0, sd=sd(mod.res$phyres)))
  qqnorm(mod.res$phyres, main="", pch=19)
  qqline(mod.res$phyres)
  mtext(text="qq-plot of residuals", side=3, line=0)
  plot(fitted(mod.res), mod.res$phyres, pch=19, col=grey(level=0.25, alpha=0.75))
  abline(h=0, lty=2)
  mtext(text="residuals against fitted values", side=3, line=0)
  par(old.par)
}

influence.pgls<-function(pgls.res){
	xterms=as.character(pgls.res$call)[-1]
  xnames=names(pgls.res$call)[-1]
  orig.data=pgls.res$data$data
  orig.phy=pgls.res$data$phy
  orig.data=data.frame(orig.data, taxon=orig.phy$tip.label)
  orig.all.data=pgls.res$data
  xterms[xnames=="data"]="sel.test.data"
  xterms[grepl(x=xterms, pattern="ML")]=paste("\"", xterms[grepl(x=xterms, pattern="ML")], "\"", sep="")
  xcall=paste(paste(xnames, xterms, sep="="), collapse=", ")
  xcall=paste(c("pgls(", xcall, ")"), collapse="")
  taxon.names=orig.phy$tip.label
  estimates=matrix(NA, nrow=length(taxon.names), ncol=length(coefficients(pgls.res)))
  SEs=matrix(NA, nrow=length(taxon.names), ncol=length(coefficients(pgls.res)))
  params=matrix(NA, nrow=length(taxon.names), ncol=3)
  old.sem=options()$show.error.messages
	options(show.error.messages=F)
  for(i in 1:length(taxon.names)){
    sel.orig.data=orig.data[-i, ]
    sel.orig.phy=drop.tip(phy=orig.phy, tip=taxon.names[i], rooted=T)
    sel.test.data=comparative.data(phy = sel.orig.phy, data = sel.orig.data, names.col = taxon, vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE)
    sel.pgls=try(eval(parse(text=xcall)))
    if(class(sel.pgls)!="try-error"){
      params[i,]=sel.pgls$param
      sel.pgls=summary(sel.pgls)$coefficients
      estimates[i,]=sel.pgls[,1]
      SEs[i,]=sel.pgls[,2]
    }
  }
	options(show.error.messages=old.sem)
  rownames(params)=taxon.names
  colnames(params)=names(pgls.res$param)
  rownames(estimates)=taxon.names
  colnames(estimates)=names(coefficients(pgls.res))
  rownames(SEs)=taxon.names
  colnames(SEs)=names(coefficients(pgls.res))
	xx.sum=cbind(summary(pgls.res)$coefficients[, "Estimate"], t(apply(estimates, 2, range, na.rm=T)))
	colnames(xx.sum)=c("orig", "min", "max")
  return(list(summary=xx.sum, taxa.failed=sum(is.na(estimates[,1])), estimates=estimates, SEs=SEs, params=params))
}

