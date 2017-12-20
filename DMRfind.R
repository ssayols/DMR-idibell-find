##################################
##
## Metode de cerca de DMR basat en la correlacio de bases properes. Es basa en que els canvis significatius
## que trobem en una CpG els trobem de manera similar a CpG del voltant
##
##    1-calcular la correlacio amb les bases veines (+/-Wkb) de totes les CpG i fer una curva
##      smooth Rt=f(x)=smooth(x=dist,y=corr) on Rt seria la correlacio teorica a una distancia 'x'
##    2-buscar DM CpG amb un t-test
##    3-mentre R(w1) dins els CI de Rt(w1) i mateixa direccio de la diferencia:
##      3.1-buscar la seguent CpG a la dreta i calcular la correlacio R(w1)=corr(x=M(0),y=M(w1))
##    4-mentre R(w2) dins els CI de Rt(w2) i mateixa direccio de la diferencia:
##      4.1-buscar la seguent CpG a la esquerra i calcular la correlacio R(w2)=corr(x=M(0),y=M(w2))
##    5-reportar la dmr de [-v,+w]
##
##################################

##
## getRt: model teoric de correlacio. Agafa totes les CpG i mostres disponibles a la taula betas, en calcula
##	la correlacio de la mitja de metilacio de totes les mostres en una CpG (coord x) respecte les seves veines
##  a distancia 'abs(d)' (coord y).
##
##	IN:
##		-betas: taula de betas de n files (CpG) i m columnes (mostres)
##		-fData450k: posicions de les CpG c("TargetID","CHR","MAPINFO")
##		-W: distancia maxima entre CpG
##		-cores: nombre de CPUs a utilitzar per paralelitzar
##
##	OUT: objecte de la classe S3
##		-df=dades calculades: distancia, correlacio, diferencia respecte fit
##		-fit=model polinomial loess
##		-se.up, se.down=error estandard del model polinomial
##		-disp.up,disp.down=dispersio calculada com un model polinomial de les distancies dels punts respecte
##			el model polinomial loess
##
getRt <- function(betas,fData450k,W=1000,cores=4) {

	# calcular les distancies entre les N CpG
	require(GenomicRanges)
	cat("Get distance between CpG sites... ")
	rd1 <- with(fData450k[fData450k$TargetID %in% rownames(betas),c("TargetID","CHR","MAPINFO")],
				GRanges(CHR,IRanges(MAPINFO,MAPINFO),TargetID=TargetID))
	mm <- findOverlaps(rd1,rd1,maxgap=W)
	sel <- data.frame(as.data.frame(rd1[as.matrix(mm)[,1],]),as.data.frame(rd1[as.matrix(mm)[,2],]))
	sel$dist <- abs(sel$start - sel$start.1) # no tenim en compte l'strand, ni si la CpG esta a prom-body-inter
	sel <- sel[sel$dist > 0,]

	# calcular les mitges de metilacio i per utilitzar-les com a coordenades x,y
	require(parallel)
	cat("done!\nCorrelate distance vs. meth average... ")
	m <- unlist(mclapply(1:nrow(betas),function(i) mean(betas[i,],na.rm=T),mc.cores=cores))
	names(m) <- rownames(betas)
	x <- by(sel,sel$dist,function(sel) {
				cor(x=m[sel$TargetID],y=m[sel$TargetID.1]) })
	cat("done!\n")

	# fer el la curva smooth + CI al voltant (CI95 == 2*SE)
	df <- data.frame(d=as.numeric(names(as.list(x))),c=unlist(as.list(x)))
	l <- loess(c ~ d,df)	# model loess corr~dist
	p <- predict(l,df$d,se=TRUE)	# calcular els SE que ens serviran per calcular el CI95
	
	# calcular dispersio dels punts respecte el model l
	df$e <- df$c - p$fit
	d1 <- loess(e ~ d,df[df$e > 0,])
	p1 <- predict(d1,df$d)
	d2 <- loess(e ~ d,df[!(df$e > 0),])
	p2 <- predict(d2,df$d)
	
	# retorn
	return(structure(list(df=df,
						  fit=l,
						  predict=p$fit,
						  se.up=p$se.fit,
						  se.down=p$se.fit,
						  disp.up=ifelse(is.na(p1),0,p1),
						  disp.down=ifelse(is.na(p2),0,p2)),
					 class="Rt"))
}

##
## funcio generica plot per la classe Rt
##
plot.Rt <- function(x,n=1,...) {

	# curva smooth + CI al voltant
	plot(x$df$d,x$df$c,cex=.2,xlab="dist",ylab="corr")			# dibuixar els punts
	lines(x$df$d,x$predict,col="red",lwd=1.5)		# dibuixar la curva smooth
	lines(x$df$d,x$predict + n * x$se.up  ,col="blue",lty=2,lwd=1.5)	# dibuixar el CI95 de dalt
	lines(x$df$d,x$predict - n * x$se.down,col="blue",lty=2,lwd=1.5)	# dibuixar el CI95 de baix

	# dibuixar dispersio dels punts respecte el model teoric
	lines(x$df$d,x$predict + n * x$disp.up  ,col="green",lty=2,lwd=1.5)
	lines(x$df$d,x$predict + n * x$disp.down,col="green",lty=2,lwd=1.5)
}

##
## findDMR: funcio que busca una regio de metilacio diferencial a partir d'una CpG diferencialment metilada
##	expandint a CpG a dreta i esquerra mentre la correlacio entre la CpG original i les del voltant sigui
##	semblant a l'esperada (Rt entre disp.up i disp.down):
##
##    2-buscar DM CpG amb un t-test (pval corregit per FDR)
##    3-mentre R(w1) dins els CI de Rt(w1):
##      3.1-buscar la seguent CpG a la dreta i calcular la correlacio R(w1)=corr(x=M(0),y=M(w1))
##    3-mentre R(w2) dins els CI de Rt(w2):
##      4.1-buscar la seguent CpG a la esquerra i calcular la correlacio R(w2)=corr(x=M(0),y=M(w2))
##    5-reportar la dmr de [-v,+w]
##
##	IN:
##		-betas: taula de betas de n files (CpG) i m columnes (mostres)
##		-A: samples' classes, as in t.test(Y ~ A)
##		-fData450k: posicions de les CpG c("TargetID","CHR","MAPINFO")
##		-Rt: correlacio teorica
##		-sig: pvalor de tall del test
##		-W: distancia maxima entre CpG
##		-cores: nombre de CPUs a utilitzar per paralelitzar
##
##	OUT:
##		un objecte de la classe dmr amb coordenades chr-ini-fin
##
findDMR <- function(betas,A,fData450k,Rt,sig=.01,W=1000,cores=4) {

	fData450k <- fData450k[fData450k$TargetID %in% rownames(betas),c("TargetID","CHR","MAPINFO")]
	fData450k <- fData450k[order(fData450k$CHR,fData450k$MAPINFO),]

	## Buscar DM CpG amb un t-test
	cat("Get differentially methylated sites (t-test)... ")
	require(parallel)
	p <- unlist(mclapply(1:nrow(betas),function(i) t.test(betas[i,] ~ A,na.rm=T)$p.value,mc.cores=cores))
	p <- p.adjust(p,method="fdr")
	names(p) <- rownames(betas)

	## Entrenar HMM amb les correlacions de les cpg diferencialment metilades
	## Suposem que una CpG no DMP correlaciona com Rt fins a +/-W, mentre que una DMP no
#	require(HMM)
#	hmm <- estimarHMM(betas,A,fData450k,p,Rt,sig,W)

	## Expansio per la dreta
	cat("done!\nExtending rightwards case 1... ")
	w1 <- mclapply(names(p)[p < sig],expandir,betas[,A == levels(A)[1]],fData450k,dir=+1,mc.cores=cores)
	cat("done!\nExtending rightwards case 2... ")
	w2 <- mclapply(names(p)[p < sig],expandir,betas[,A == levels(A)[2]],fData450k,dir=+1,mc.cores=cores)

	## Expansio per l'esquerra
	cat("done!\nExtending leftwards case 1... ")
	v1 <- mclapply(names(p)[p < sig],expandir,betas[,A == levels(A)[1]],fData450k,dir=-1,mc.cores=cores)
	cat("done!\nExtending leftwards case 2... ")
	v2 <- mclapply(names(p)[p < sig],expandir,betas[,A == levels(A)[2]],fData450k,dir=-1,mc.cores=cores)
	cat("done!\n")

	## reportar la dmr de [-v,+w]
	dmr <- sapply(1:sum(p < sig),function(i) rbind(w1[[i]],w2[[i]],v1[[i]],v2[[i]]))
	dmr.ini.fin <- lapply(dmr,function(x) {
		if(nrow(x) > 0) 
			data.frame(chr=unique(x$chr),ini=min(c(x$p1,x$p2)),fin=max(c(x$p1,x$p2))) 
		else
			data.frame()
	})
	return(structure(list(dmr=Reduce(rbind,dmr.ini.fin)),class="dmr"))
}

##
## Expansio: funcio interna utilitzada per findDMR
##
expandir <- function(TargetID,betas,fData450k,dir=1,n=1,max.depth=10,i=0,W=1000) {

	chr <- fData450k$CHR[fData450k$TargetID == TargetID]
	pos <- fData450k$MAPINFO[fData450k$TargetID == TargetID]
	if(dir > 0)
		sel <- fData450k[fData450k$CHR == chr & fData450k$MAPINFO > pos & fData450k$MAPINFO < pos + W,]
	else
		sel <- fData450k[fData450k$CHR == chr & fData450k$MAPINFO > pos - W & fData450k$MAPINFO < pos,]

	# reportar la correlacio amb les seves veines i expandir-les si correlacionen
	r <- data.frame()
	if(nrow(sel) > 0) {	# si te veins
		sel$dist <- abs(sel$MAPINFO - pos)
		sel <- sel[order(sel$dist),]

		for(i in 1:nrow(sel)) {
			corr <- abs(cor(x=betas[TargetID,],y=betas[sel$TargetID[i],]))
			down <- Rt$predict[sel$dist[i]] + n * Rt$disp.down[sel$dist[i]]
#			up   <- Rt$predict[sel$dist[i]] + n * Rt$disp.up  [sel$dist[i]]
			if(corr > down) {	# expandir la DMR si la correlacio es mes alta del que esperem
				r <- rbind(r,data.frame(t1=TargetID,t2=sel$TargetID[i],chr=chr,p1=pos,p2=sel$MAPINFO[i],c=corr))
				if(i < max.depth) {	# continuar l'expansio si no hem passat el nombre maxim de nivells
					r <- rbind(r,expandir(sel$TargetID[i],betas,fData450k,dir,n,max.depth,i=i+1,W))
				}
			}
			else {	# no correlaciona, s'ha trencat la DMR i ja no tractem d'expandir mes per aquest costat
				break
			}
		}
	}

	return(r)
}

##
## plot.drm: dibuixar la dmr utilitzant boxplots a les posicions on hi ha les CpG
##
plot.dmr <- function(x,A,W=1000,f="plot.dmr.pdf",cores=4,...) {

	require(ggplot2)
	require(reshape2)
	require(parallel)

	names(A) <- colnames(betas)

	# recorrer les DMR i dibuixar totes les CpG com a boxplots a +/-W parells de base
	p <- mclapply(1:nrow(x$dmr),function(i) {

		# dades de metilacio de les CpG dins la finestra
		cg <- fData450k[fData450k$CHR == x$dmr$chr[i] & 
						fData450k$MAPINFO >= x$dmr$ini[i] - W & 
						fData450k$MAPINFO <= x$dmr$fin[i] + W,c("TargetID","MAPINFO")]
		df <- melt(betas[cg$TargetID,])
		df <- merge(df,as.data.frame(A),by.x="Var2",by.y=0)	# afegir el grup
		df <- merge(df,cg,by.x="Var1",by.y="TargetID")

		# dmr
		df2 <- data.frame(xini=x$dmr$ini[i],xfin=x$dmr$fin[i],yini=-0.01,yfin=-0.01)

		p <- ggplot(df,aes(x=MAPINFO,y=value,colour=A)) + 
				geom_boxplot(aes(group=paste(MAPINFO,A)),position="identity") +
				geom_smooth() +
				geom_segment(data=df2,aes(x=xini,y=yini,xend=xfin,yend=yfin),colour="red",size=2) +
				geom_point(shape=19) +
				ggtitle(paste0(x$dmr$chr[i],":",x$dmr$ini[i] - W,"-",x$dmr$fin[i] + W))

		return(p)
	},mc.cores=cores)

	# fer un PDF
	pdf(f)
	lapply(p,print)
	dev.off()
}
