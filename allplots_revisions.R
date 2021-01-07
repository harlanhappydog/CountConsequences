
alpha_sig_level <- 0.05
lvls <- list()

exp(c(0.5, 1, 1.5, 2, 2.5))

lvls[[1]] <- c(50, 100, 250, 500, 
                  1000, 2000, 5000, 10000)		    # n
lvls[[2]] <- c(0.5, 1, 1.5, 2, 2.5)					# beta_0
lvls[[3]] <- c(Inf, 4/2, 1, 1/2, 1/3)	 			# phi
lvls[[4]] <- c(0, 0.05, 0.1, 0.2, 0.5)				# omega

dsgn <- as.matrix(expand.grid(lvls[[1]], lvls[[2]], lvls[[3]], lvls[[4]]))

colnames(dsgn)<-c("n","beta_0","phi","omega")

dim(dsgn)


resultsmat <- data.frame(matrix(0,dim(dsgn)[1], 41))
colnames(resultsmat)<- c("iii",
"n", "beta_0", "phi", "omega", 
"LRT_pval", 
"vuong_P_zip_pval", 
"vuong_NB_zinb_pval", 
"Poisson_pval",
"nb_pval",
"zip_pval",
"zinb_pval",
"pval1",
"pvalAIC",
"pvalAICc",
"pvalBIC",
"correct_pval",
"correct_selection",
"correct_selectionAIC",
"correct_selectionAICc",
"correct_selectionBIC",
"selectionLRT",
"selectionvuongP",
"selectionvuongNB",
"propPoisson",
"propNB",
"propZIP",
"propZINB",
"propPoissonAIC",
"propPoissonAICc",
"propNBAIC",
"propNBAICc",
"propZIPAIC",
"propZIPAICc",
"propZINBAIC",
"propZINBAICc",
"propPoissonBIC",
"propNBBIC",
"propZIPBIC",
"propZINBBIC",
"nSim")



#for(iii in 1:10){
for(iii in 1:dim(dsgn)[1]){
	
	tryCatch({

results_iii <- readRDS(paste(paste("~/cluster_results/results_Nov2020_15000",iii, sep="_"),".rds",sep=""))	


LRT_ <- mean(results_iii[,"LRT_pval"]<0.05)
A_ <- mean(results_iii[results_iii[,"LRT_pval"]>=0.05, "vuong_P_zip_pval"]<0.05)
B_ <- mean(results_iii[results_iii[,"LRT_pval"]<0.05, "vuong_NB_zinb_pval"]<0.05)
C_ <- mean(results_iii[results_iii[,"LRT_pval"]>=0.05 & results_iii[,"vuong_P_zip_pval"]>=0.05, "poisson_pval"]<0.05)
D_ <- mean(results_iii[results_iii[,"LRT_pval"]>=0.05 & results_iii[,"vuong_P_zip_pval"]<0.05, "zip_pval"]<0.05)
E_ <- mean(results_iii[results_iii[,"LRT_pval"]<0.05 & results_iii[,"vuong_NB_zinb_pval"]>=0.05, "nb_pval"]<0.05)
F_ <- mean(results_iii[results_iii[,"LRT_pval"]<0.05 & results_iii[,"vuong_NB_zinb_pval"]<0.05, "zinb_pval"]<0.05)


C_all<-mean(results_iii[, "poisson_pval"]<0.05)
D_all<-mean(results_iii[, "zip_pval"]<0.05)
E_all<-mean(results_iii[, "nb_pval"]<0.05)
F_all<-mean(results_iii[, "zinb_pval"]<0.05)


rejections <- apply(results_iii[,6:16],2, function(x){mean(na.omit(x)<= alpha_sig_level)})


if(results_iii[1,"p_ZI"]==0 & results_iii[1,"phi"]==Inf){
	rejections["correct_pval"] <- rejections["poisson_pval"]
	rejections["correct_selection"] <- mean(results_iii[,"choice"]==1)
	rejections["correct_selectionAIC"] <- mean(results_iii[,"choiceAIC"]==1)
	rejections["correct_selectionAICc"] <- mean(results_iii[,"choiceAICc"]==1)	
	rejections["correct_selectionBIC"] <- mean(results_iii[,"choiceBIC"]==1)
	}

if(results_iii[1,"p_ZI"]==0 & results_iii[1,"phi"]<Inf){
	rejections["correct_pval"]<-rejections["nb_pval"]
	rejections["correct_selection"] <- mean(results_iii[,"choice"]==2)
	rejections["correct_selectionAIC"] <- mean(results_iii[,"choiceAIC"]==2)
	rejections["correct_selectionAICc"] <- mean(results_iii[,"choiceAICc"]==2)	
	rejections["correct_selectionBIC"] <- mean(results_iii[,"choiceBIC"]==2)
	}

if(results_iii[1,"p_ZI"]>0 & results_iii[1,"phi"]==Inf){
	rejections["correct_pval"] <- rejections["zip_pval"]
	rejections["correct_selection"] <- mean(results_iii[,"choice"]==3)
	rejections["correct_selectionAIC"] <- mean(results_iii[,"choiceAIC"]==3)
	rejections["correct_selectionAICc"] <- mean(results_iii[,"choiceAICc"]==3)	
	rejections["correct_selectionBIC"] <- mean(results_iii[,"choiceBIC"]==3)
	}

 if(results_iii[1,"p_ZI"]>0 & results_iii[1,"phi"]<Inf){
 	rejections["correct_pval"]<-rejections["zinb_pval"]
	rejections["correct_selection"] <- mean(results_iii[,"choice"]==4)
	rejections["correct_selectionAIC"] <- mean(results_iii[,"choiceAIC"]==4)
	rejections["correct_selectionAICc"] <- mean(results_iii[,"choiceAICc"]==4)	
	rejections["correct_selectionBIC"] <- mean(results_iii[,"choiceBIC"]==4)
 	}


rejections["selectionLRT"]<-mean((results_iii[,"LRT_pval"]<=alpha_sig_level), na.rm=TRUE)		
rejections["selectionvuongP"]<-mean((results_iii[,"vuong_P_zip_pval"]<=alpha_sig_level), na.rm=TRUE)
rejections["selectionvuongNB"]<-mean((results_iii[,"vuong_NB_zinb_pval"]<=alpha_sig_level), na.rm=TRUE)
	
rejections["propPoisson"] <- mean(results_iii[,"choice"]==1)
rejections["propNB"] <- mean(results_iii[,"choice"]==2)
rejections["propZIP"] <- mean(results_iii[,"choice"]==3)
rejections["propZINB"] <- mean(results_iii[,"choice"]==4)

rejections["propPoissonAIC"] <- mean(results_iii[,"choiceAIC"]==1)
rejections["propNBAIC"] <- mean(results_iii[,"choiceAIC"]==2)
rejections["propZIPAIC"] <- mean(results_iii[,"choiceAIC"]==3)
rejections["propZINBAIC"] <- mean(results_iii[,"choiceAIC"]==4)

rejections["propPoissonAICc"] <- mean(results_iii[,"choiceAICc"]==1)
rejections["propNBAICc"] <- mean(results_iii[,"choiceAICc"]==2)
rejections["propZIPAICc"] <- mean(results_iii[,"choiceAICc"]==3)
rejections["propZINBAICc"] <- mean(results_iii[,"choiceAICc"]==4)

rejections["propPoissonBIC"] <- mean(results_iii[,"choiceBIC"]==1)
rejections["propNBBIC"] <- mean(results_iii[,"choiceBIC"]==2)
rejections["propZIPBIC"] <- mean(results_iii[,"choiceBIC"]==3)
rejections["propZINBBIC"] <- mean(results_iii[,"choiceBIC"]==4)

rejections
allresults <- c(rejections)


	resultsmat
	simdat <-as.data.frame(c(iii=iii,unlist(c(dsgn[iii,], allresults, results_iii[1,"nSim"]))))

colnames(resultsmat) <- rownames(as.data.frame(simdat))
	resultsmat[iii,] <- unlist(as.data.frame(simdat))
	print(resultsmat[iii,])
	allresults<-NULL
	results_iii<-NULL
	
}

, 
		 error=function(e){})
}


# should be integer(0)
c(1:dim(resultsmat)[1])[resultsmat$n==0]
# should be 1
mean(rowSums((resultsmat[,c("propPoisson", "propZIP", "propNB", "propZINB")])))


head(resultsmat);tail(resultsmat)
resultsmat <- resultsmat[resultsmat$n>0,]
dim(resultsmat)
library(ggplot2)


dim(resultsmat)
resultsmat[,"prop_inflate"]<-	0

for(jjj in 1:dim(resultsmat)[1]){
	if(resultsmat[jjj,"Poisson_pval"] > resultsmat[jjj,"correct_pval"]){
	resultsmat[jjj,"prop_inflate"]<-	resultsmat[jjj,"prop_inflate"]+resultsmat[jjj,"propPoisson"]
	}

	if(resultsmat[jjj,"nb_pval"] > resultsmat[jjj,"correct_pval"]){
	resultsmat[jjj,"prop_inflate"]<-	resultsmat[jjj,"prop_inflate"]+resultsmat[jjj,"propNB"]
	}
	
	if(resultsmat[jjj,"zip_pval"] > resultsmat[jjj,"correct_pval"]){
	resultsmat[jjj,"prop_inflate"]<-	resultsmat[jjj,"prop_inflate"]+resultsmat[jjj,"propZIP"]
	}
	
	if(resultsmat[jjj,"zinb_pval"] > resultsmat[jjj,"correct_pval"]){
	resultsmat[jjj,"prop_inflate"]<-	resultsmat[jjj,"prop_inflate"]+resultsmat[jjj,"propZINB"]
	}
}
tail(resultsmat)


dim(resultsmat)
resultsmat[,"prop_inflate2"] <-	0

for(jjj in 1:dim(resultsmat)[1]){
	if((resultsmat[jjj,"Poisson_pval"] > resultsmat[jjj,"correct_pval"]) & (resultsmat[jjj,"Poisson_pval"] > 0.05)){
	resultsmat[jjj,"prop_inflate2"]<-	resultsmat[jjj,"prop_inflate2"]+resultsmat[jjj,"propPoisson"]
	}

	if((resultsmat[jjj,"nb_pval"] > resultsmat[jjj,"correct_pval"]) & (resultsmat[jjj,"nb_pval"]  > 0.05)){
	resultsmat[jjj,"prop_inflate2"]<-	resultsmat[jjj,"prop_inflate2"] +resultsmat[jjj,"propNB"]
	}
	
	if((resultsmat[jjj,"zip_pval"] > resultsmat[jjj,"correct_pval"]) & (resultsmat[jjj,"zip_pval"] > 0.05) ){
	resultsmat[jjj,"prop_inflate2"]<-	resultsmat[jjj,"prop_inflate2"]+resultsmat[jjj,"propZIP"]
	}
	
	if((resultsmat[jjj,"zinb_pval"] > resultsmat[jjj,"correct_pval"]) & (resultsmat[jjj,"zinb_pval"] >0.05)){
	resultsmat[jjj,"prop_inflate2"]<-	resultsmat[jjj,"prop_inflate2"]+resultsmat[jjj,"propZINB"]
	}
}
head(resultsmat)
resultsmatsmall <- resultsmat


head(resultsmatsmall)
  resultsmatsmall[resultsmatsmall[,"phi"]==Inf & resultsmatsmall[,"omega"]==0,]  
  
  
resultsmatsmall<-  resultsmatsmall[resultsmatsmall$n>0 ,]
resultsmatsmall$type1 <- resultsmatsmall$pval1
resultsmatsmall$typeAIC <- resultsmatsmall$pvalAIC
resultsmatsmall$typeAICc <- resultsmatsmall$pvalAICc
resultsmatsmall$typeBIC <- resultsmatsmall$pvalBIC
resultsmatsmall$correct <- resultsmatsmall$correct_pval 


resultsmatsmall$g <- as.factor(paste(as.character(resultsmatsmall$phi),as.character(resultsmatsmall$omega), as.character(resultsmatsmall$beta_0), sep="_"))


resultsmatsmall$diff <- resultsmatsmall$type1-resultsmatsmall$correct


max(resultsmatsmall[resultsmatsmall[,"phi"]<Inf & resultsmatsmall[,"omega"]==0,][,"zip_pval"] )
max(resultsmatsmall[resultsmatsmall[,"phi"]<Inf & resultsmatsmall[,"omega"]==0,][,"Poisson_pval"] )

unique(resultsmatsmall[,"omega"])
unique(resultsmatsmall[,"phi"])

resultsmatsmall2 <-transform(resultsmatsmall, phi =factor(phi, levels=sort(unique(resultsmatsmall$phi)), c(
 expression(paste(phi," = 1/3")), 
 expression(paste(phi," = 1/2")), 
 expression(paste(phi," = ",1)), 
 expression(paste(phi," = ",2)),  
 expression(paste(phi," = ",Inf)) ) ))

resultsmatsmall<-resultsmatsmall2

unique(resultsmatsmall$phi)

unique(resultsmatsmall$omega)

resultsmatsmall <-transform(resultsmatsmall, omega =factor(omega, levels=sort(unique(resultsmatsmall$omega)), c(
 expression(paste(omega," = ",0.0)),
 expression(paste(omega," = ",0.05)), 
 expression(paste(omega," = ",0.10)), 
 expression(paste(omega," = ",0.20)), 
 expression(paste(omega," = ",0.50)) ) ))


resultsmatsmall$lll <- (as.numeric(as.factor(paste(as.character(resultsmatsmall$omega),as.character(resultsmatsmall$phi), sep="_"))))

qqq<-1
resultsmatsmall$lab<-NULL

for(i in 1:(length(resultsmatsmall$lll))){
resultsmatsmall$lab[i] <- qqq
if(!resultsmatsmall$lll[i]==resultsmatsmall$lll[i+1]){qqq<-qqq+1}
}
resultsmatsmall$lab
resultsmatsmall$lab[resultsmatsmall$lab[-length(resultsmatsmall$lab)]==resultsmatsmall$lab[-1]]<-""
resultsmatsmall$lab[length(resultsmatsmall$lab)]<-"25"


#######  correct
p <-ggplot(resultsmatsmall, aes(x = n, y = correct, group = g, colour=as.factor(beta_0))) + geom_point() + geom_line()  + 

facet_grid(phi ~omega, labeller=label_parsed) + 
labs(x = paste("sample size (n)"), y = expression("probability of p" < alpha)) + 
coord_cartesian(xlim = c(40, 15000), ylim = c(0, 0.1)) + 
guides(color=guide_legend(title = expression(paste("intercept parameter (", beta[0], ") : "))))  +  
theme(axis.title = element_text(size = 15), axis.text.x = element_text(vjust = 3), legend.text=element_text(size=15), axis.ticks = element_blank(), legend.title=element_text(size=16), axis.text = element_text(size=10.5),   panel.grid.minor.x = element_blank(), plot.title=element_text(size=17), strip.text.x = element_text(size = 13, colour = "black", angle = 0), strip.text.y = element_text(size = 12, colour = "black", angle = 270))+ 

geom_hline(yintercept = alpha_sig_level) + labs(title="Type 1 error obtained with correct model")

pdf("~/plots/correct.pdf", width=9.5, height=7)


p +  geom_text(
    size    = 4,
    data    = resultsmatsmall,
    mapping = aes(x = Inf, y = Inf, label = as.numeric(lab)),
    hjust   = 1.05,
    vjust   = 1.5, col="black"
  ) + 
  scale_y_continuous(breaks = c(0, 0.025, 0.05, 0.075), labels=c("0", "0.025", "0.050", "0.075"), sec.axis = sec_axis(trans ~ ., labels = NULL,  name = expression(" " %<-% "direction of increasing overdispersion" %<-% "\n"), t))  + theme(legend.position="bottom") +
  scale_x_continuous(trans='log10', breaks = c(50, 250, 1000, 10000), sec.axis = dup_axis(trans ~ ., labels = NULL, name = expression(" " %->% "direction of increasing zero-inflation" %->% "\n"))) 

dev.off()
#######  correct_selection
p <-ggplot(resultsmatsmall, aes(x = n, y = correct_selection, group = g, colour=as.factor(beta_0))) + geom_point() + geom_line()  + 

facet_grid(phi ~omega, labeller=label_parsed) + 
labs(x = paste("sample size (n)"), y = expression("probability of selecting correct model")) + 
coord_cartesian(xlim = c(40, 15000), ylim = c(0, 1)) + 
guides(color=guide_legend(title = expression(paste("intercept parameter (", beta[0], ") : "))))  +  
theme(axis.title = element_text(size = 15), axis.text.x = element_text(vjust = 3), legend.text=element_text(size=15), axis.ticks = element_blank(), legend.title=element_text(size=16), axis.text = element_text(size=10.5),   panel.grid.minor.x = element_blank(), plot.title=element_text(size=17), strip.text.x = element_text(size = 13, colour = "black", angle = 0), strip.text.y = element_text(size = 12, colour = "black", angle = 270))+ 

 labs(title="Probability of selecting the correct model following the 7-step testing scheme")

pdf("~/plots/correct_selection.pdf", width=9.5, height=7)


p +  geom_text(
    size    = 4,
    data    = resultsmatsmall,
    mapping = aes(x = Inf, y = Inf, label = as.numeric(lab)),
    hjust   = 1.05,
    vjust   = 1.5, col="black"
  ) + 
  scale_y_continuous(sec.axis = sec_axis(trans ~ ., labels = NULL,  name = expression(" " %<-% "direction of increasing overdispersion" %<-% "\n"), t))  + theme(legend.position="bottom") +
  scale_x_continuous(trans='log10', breaks = c(50, 250, 1000, 10000), sec.axis = dup_axis(trans ~ ., labels = NULL, name = expression(" " %->% "direction of increasing zero-inflation" %->% "\n"))) 

dev.off()


#######  type1
p <-ggplot(resultsmatsmall, aes(x = n, y = type1, group = g, colour=as.factor(beta_0))) + geom_point() + geom_line()  + 

facet_grid(phi ~omega, labeller=label_parsed) + 
labs(x = paste("sample size (n)"), y = expression("probability of p" < alpha)) + 
coord_cartesian(xlim = c(40, 15000), ylim = c(0, 0.1)) + 
guides(color=guide_legend(title = expression(paste("intercept parameter (", beta[0], ") : "))))  +  
theme(axis.title = element_text(size = 15), axis.text.x = element_text(vjust = 3), legend.text=element_text(size=15), axis.ticks = element_blank(), legend.title=element_text(size=16), axis.text = element_text(size=10.5),   panel.grid.minor.x = element_blank(), plot.title=element_text(size=17), strip.text.x = element_text(size = 13, colour = "black", angle = 0), strip.text.y = element_text(size = 12, colour = "black", angle = 270))+ 

geom_hline(yintercept = alpha_sig_level) + labs(title="Type 1 error obtained following the 7-step testing scheme")

pdf("~/plots/type1.pdf", width=9.5, height=7)


p +  geom_text(
    size    = 4,
    data    = resultsmatsmall,
    mapping = aes(x = Inf, y = Inf, label = as.numeric(lab)),
    hjust   = 1.05,
    vjust   = 1.5, col="black"
  ) + 
  scale_y_continuous(breaks = c(0, 0.025, 0.05, 0.075), labels=c("0", "0.025", "0.050", "0.075"), sec.axis = sec_axis(trans ~ ., labels = NULL,  name = expression(" " %<-% "direction of increasing overdispersion" %<-% "\n"), t))  + theme(legend.position="bottom") +
  scale_x_continuous(trans='log10', breaks = c(50, 250, 1000, 10000), sec.axis = dup_axis(trans ~ ., labels = NULL, name = expression(" " %->% "direction of increasing zero-inflation" %->% "\n"))) 

dev.off()


#######  correct_selectionAIC
p <-ggplot(resultsmatsmall, aes(x = n, y = correct_selectionAIC, group = g, colour=as.factor(beta_0))) + geom_point() + geom_line()  + 

facet_grid(phi ~omega, labeller=label_parsed) + 
labs(x = paste("sample size (n)"), y = expression("probability of selecting correct model")) + 
coord_cartesian(xlim = c(40, 15000), ylim = c(0, 1)) + 
guides(color=guide_legend(title = expression(paste("intercept parameter (", beta[0], ") : "))))  +  
theme(axis.title = element_text(size = 15), axis.text.x = element_text(vjust = 3), legend.text=element_text(size=15), axis.ticks = element_blank(), legend.title=element_text(size=16), axis.text = element_text(size=10.5),   panel.grid.minor.x = element_blank(), plot.title=element_text(size=17), strip.text.x = element_text(size = 13, colour = "black", angle = 0), strip.text.y = element_text(size = 12, colour = "black", angle = 270))+ 

labs(title="Probability that the correct model is the one with the lowest AIC")

pdf("~/plots/correct_selectionAIC.pdf", width=9.5, height=7)

p +  geom_text(
    size    = 4,
    data    = resultsmatsmall,
    mapping = aes(x = Inf, y = Inf, label = as.numeric(lab)),
    hjust   = 1.05,
    vjust   = 1.5, col="black"
  ) + 
  scale_y_continuous(sec.axis = sec_axis(trans ~ ., labels = NULL,  name = expression(" " %<-% "direction of increasing overdispersion" %<-% "\n"), t))  + theme(legend.position="bottom") +
  scale_x_continuous(trans='log10', breaks = c(50, 250, 1000, 10000), sec.axis = dup_axis(trans ~ ., labels = NULL, name = expression(" " %->% "direction of increasing zero-inflation" %->% "\n"))) 
dev.off()


#######  correct_selectionBIC
p <-ggplot(resultsmatsmall, aes(x = n, y = correct_selectionBIC, group = g, colour=as.factor(beta_0))) + geom_point() + geom_line()  + 

facet_grid(phi ~omega, labeller=label_parsed) + 
labs(x = paste("sample size (n)"), y = expression("probability of selecting correct model")) + 
coord_cartesian(xlim = c(40, 15000), ylim = c(0, 1)) + 
guides(color=guide_legend(title = expression(paste("intercept parameter (", beta[0], ") : "))))  +  
theme(axis.title = element_text(size = 15), axis.text.x = element_text(vjust = 3), legend.text=element_text(size=15), axis.ticks = element_blank(), legend.title=element_text(size=16), axis.text = element_text(size=10.5),   panel.grid.minor.x = element_blank(), plot.title=element_text(size=17), strip.text.x = element_text(size = 13, colour = "black", angle = 0), strip.text.y = element_text(size = 12, colour = "black", angle = 270))+ 

labs(title="Probability that the correct model is the one with the lowest BIC")

pdf("~/plots/correct_selectionBIC.pdf", width=9.5, height=7)

p +  geom_text(
    size    = 4,
    data    = resultsmatsmall,
    mapping = aes(x = Inf, y = Inf, label = as.numeric(lab)),
    hjust   = 1.05,
    vjust   = 1.5, col="black"
  ) + 
  scale_y_continuous(sec.axis = sec_axis(trans ~ ., labels = NULL,  name = expression(" " %<-% "direction of increasing overdispersion" %<-% "\n"), t))  + theme(legend.position="bottom") +
  scale_x_continuous(trans='log10', breaks = c(50, 250, 1000, 10000), sec.axis = dup_axis(trans ~ ., labels = NULL, name = expression(" " %->% "direction of increasing zero-inflation" %->% "\n"))) 

dev.off()


#######  correct_selectionAICc
p <-ggplot(resultsmatsmall, aes(x = n, y = correct_selectionAICc, group = g, colour=as.factor(beta_0))) + geom_point() + geom_line()  + 

facet_grid(phi ~omega, labeller=label_parsed) + 
labs(x = paste("sample size (n)"), y = expression("probability of selecting correct model")) + 
coord_cartesian(xlim = c(40, 15000), ylim = c(0, 1)) + 
guides(color=guide_legend(title = expression(paste("intercept parameter (", beta[0], ") : "))))  +  
theme(axis.title = element_text(size = 15), axis.text.x = element_text(vjust = 3), legend.text=element_text(size=15), axis.ticks = element_blank(), legend.title=element_text(size=16), axis.text = element_text(size=10.5),   panel.grid.minor.x = element_blank(), plot.title=element_text(size=17), strip.text.x = element_text(size = 13, colour = "black", angle = 0), strip.text.y = element_text(size = 12, colour = "black", angle = 270))+ 

 labs(title="Probability that the correct model is the one with the lowest AICc")

pdf("~/plots/correct_selectionAICc.pdf", width=9.5, height=7)


p +  geom_text(
    size    = 4,
    data    = resultsmatsmall,
    mapping = aes(x = Inf, y = Inf, label = as.numeric(lab)),
    hjust   = 1.05,
    vjust   = 1.5, col="black"
  ) + 
  scale_y_continuous(sec.axis = sec_axis(trans ~ ., labels = NULL,  name = expression(" " %<-% "direction of increasing overdispersion" %<-% "\n"), t))  + theme(legend.position="bottom") +
  scale_x_continuous(trans='log10', breaks = c(50, 250, 1000, 10000), sec.axis = dup_axis(trans ~ ., labels = NULL, name = expression(" " %->% "direction of increasing zero-inflation" %->% "\n"))) 

dev.off()


#######  typeAIC
p <-ggplot(resultsmatsmall, aes(x = n, y = typeAIC, group = g, colour=as.factor(beta_0))) + geom_point() + geom_line()  + 

facet_grid(phi ~omega, labeller=label_parsed) + 
labs(x = paste("sample size (n)"), y = expression("probability of p" < alpha)) + 
coord_cartesian(xlim = c(40, 15000), ylim = c(0, 0.1)) + 
guides(color=guide_legend(title = expression(paste("intercept parameter (", beta[0], ") : "))))  +  
theme(axis.title = element_text(size = 15), axis.text.x = element_text(vjust = 3), legend.text=element_text(size=15), axis.ticks = element_blank(), legend.title=element_text(size=16), axis.text = element_text(size=10.5),   panel.grid.minor.x = element_blank(), plot.title=element_text(size=17), strip.text.x = element_text(size = 13, colour = "black", angle = 0), strip.text.y = element_text(size = 12, colour = "black", angle = 270))+ 

geom_hline(yintercept = alpha_sig_level) + labs(title="Type  1  error  obtained from model with the lowest AIC")


pdf("~/plots/typeAIC.pdf", width=9.5, height=7)

p +  geom_text(
    size    = 4,
    data    = resultsmatsmall,
    mapping = aes(x = Inf, y = Inf, label = as.numeric(lab)),
    hjust   = 1.05,
    vjust   = 1.5, col="black"
  ) + 
  scale_y_continuous(breaks = c(0, 0.025, 0.05, 0.075), labels=c("0", "0.025", "0.050", "0.075"), sec.axis = sec_axis(trans ~ ., labels = NULL,  name = expression(" " %<-% "direction of increasing overdispersion" %<-% "\n"), t))  + theme(legend.position="bottom") +
  scale_x_continuous(trans='log10', breaks = c(50, 250, 1000, 10000), sec.axis = dup_axis(trans ~ ., labels = NULL, name = expression(" " %->% "direction of increasing zero-inflation" %->% "\n"))) 
dev.off()



#######  typeBIC
p <-ggplot(resultsmatsmall, aes(x = n, y = typeBIC, group = g, colour=as.factor(beta_0))) + geom_point() + geom_line()  + 

facet_grid(phi ~omega, labeller=label_parsed) + 
labs(x = paste("sample size (n)"), y = expression("probability of p" < alpha)) + 
coord_cartesian(xlim = c(40, 15000), ylim = c(0, 0.1)) + 
guides(color=guide_legend(title = expression(paste("intercept parameter (", beta[0], ") : "))))  +  
theme(axis.title = element_text(size = 15), axis.text.x = element_text(vjust = 3), legend.text=element_text(size=15), axis.ticks = element_blank(), legend.title=element_text(size=16), axis.text = element_text(size=10.5),   panel.grid.minor.x = element_blank(), plot.title=element_text(size=17), strip.text.x = element_text(size = 13, colour = "black", angle = 0), strip.text.y = element_text(size = 12, colour = "black", angle = 270))+ 

geom_hline(yintercept = alpha_sig_level) + labs(title="Type  1  error  obtained from model with the lowest BIC")
pdf("~/plots/typeBIC.pdf", width=9.5, height=7)

p +  geom_text(
    size    = 4,
    data    = resultsmatsmall,
    mapping = aes(x = Inf, y = Inf, label = as.numeric(lab)),
    hjust   = 1.05,
    vjust   = 1.5, col="black"
  ) + 
  scale_y_continuous(breaks = c(0, 0.025, 0.05, 0.075), labels=c("0", "0.025", "0.050", "0.075"), sec.axis = sec_axis(trans ~ ., labels = NULL,  name = expression(" " %<-% "direction of increasing overdispersion" %<-% "\n"), t))  + theme(legend.position="bottom") +
  scale_x_continuous(trans='log10', breaks = c(50, 250, 1000, 10000), sec.axis = dup_axis(trans ~ ., labels = NULL, name = expression(" " %->% "direction of increasing zero-inflation" %->% "\n"))) 

dev.off() 



#######  typeAICc
p <-ggplot(resultsmatsmall, aes(x = n, y = typeAICc, group = g, colour=as.factor(beta_0))) + geom_point() + geom_line()  + 

facet_grid(phi ~omega, labeller=label_parsed) + 
labs(x = paste("sample size (n)"), y = expression("probability of p" < alpha)) + 
coord_cartesian(xlim = c(40, 15000), ylim = c(0, 0.1)) + 
guides(color=guide_legend(title = expression(paste("intercept parameter (", beta[0], ") : "))))  +  
theme(axis.title = element_text(size = 15), axis.text.x = element_text(vjust = 3), legend.text=element_text(size=15), axis.ticks = element_blank(), legend.title=element_text(size=16), axis.text = element_text(size=10.5),   panel.grid.minor.x = element_blank(), plot.title=element_text(size=17), strip.text.x = element_text(size = 13, colour = "black", angle = 0), strip.text.y = element_text(size = 12, colour = "black", angle = 270))+ 

geom_hline(yintercept = alpha_sig_level) + labs(title="Type  1  error  obtained from model with the lowest AICc")

pdf("~/plots/typeAICc.pdf", width=9.5, height=7)

p +  geom_text(
    size    = 4,
    data    = resultsmatsmall,
    mapping = aes(x = Inf, y = Inf, label = as.numeric(lab)),
    hjust   = 1.05,
    vjust   = 1.5, col="black"
  ) + 
  scale_y_continuous(breaks = c(0, 0.025, 0.05, 0.075), labels=c("0", "0.025", "0.050", "0.075"), sec.axis = sec_axis(trans ~ ., labels = NULL,  name = expression(" " %<-% "direction of increasing overdispersion" %<-% "\n"), t))  + theme(legend.position="bottom") +
  scale_x_continuous(trans='log10', breaks = c(50, 250, 1000, 10000), sec.axis = dup_axis(trans ~ ., labels = NULL, name = expression(" " %->% "direction of increasing zero-inflation" %->% "\n"))) 

dev.off()


#######  selectionLRT
p <-ggplot(resultsmatsmall, aes(x = n, y = selectionLRT, group = g, colour=as.factor(beta_0))) + geom_point() + geom_line()  + 

facet_grid(phi ~omega, labeller=label_parsed) + 
labs(x = paste("sample size (n)"), y = expression("probability of rejecting the null")) + 
coord_cartesian(xlim = c(40, 15000), ylim = c(0, 1)) + 
guides(color=guide_legend(title = expression(paste("intercept parameter (", beta[0], ") : "))))  +  
theme(axis.title = element_text(size = 15), axis.text.x = element_text(vjust = 3), legend.text=element_text(size=15), axis.ticks = element_blank(), legend.title=element_text(size=16), axis.text = element_text(size=10.5),   panel.grid.minor.x = element_blank(), plot.title=element_text(size=17), strip.text.x = element_text(size = 13, colour = "black", angle = 0), strip.text.y = element_text(size = 12, colour = "black", angle = 270))+ 

labs(title="Probability  that  the  D&L  test  rejects  the  null  of no overdispersion")

pdf("~/plots/selectionLRT.pdf", width=9.5, height=7)

p +  geom_text(
    size    = 4,
    data    = resultsmatsmall,
    mapping = aes(x = Inf, y = Inf, label = as.numeric(lab)),
    hjust   = 1.05,
    vjust   = 1.5, col="black"
  ) + 
  scale_y_continuous(sec.axis = sec_axis(trans ~ ., labels = NULL,  name = expression(" " %<-% "direction of increasing overdispersion" %<-% "\n"), t))  + theme(legend.position="bottom") +
  scale_x_continuous(trans='log10', breaks = c(50, 250, 1000, 10000), sec.axis = dup_axis(trans ~ ., labels = NULL, name = expression(" " %->% "direction of increasing zero-inflation" %->% "\n"))) 

dev.off()


#######  selectionvuongP
p <-ggplot(resultsmatsmall, aes(x = n, y = selectionvuongP, group = g, colour=as.factor(beta_0))) + geom_point() + geom_line()  + 

facet_grid(phi ~omega, labeller=label_parsed) + 
labs(x = paste("sample size (n)"), y = expression("probability of rejecting the null")) + 
coord_cartesian(xlim = c(40, 15000), ylim = c(0, 1)) + 
guides(color=guide_legend(title = expression(paste("intercept parameter (", beta[0], ") : "))))  +  
theme(axis.title = element_text(size = 15), axis.text.x = element_text(vjust = 3), legend.text=element_text(size=15), axis.ticks = element_blank(), legend.title=element_text(size=16), axis.text = element_text(size=10.5),   panel.grid.minor.x = element_blank(), plot.title=element_text(size=17), strip.text.x = element_text(size = 13, colour = "black", angle = 0), strip.text.y = element_text(size = 12, colour = "black", angle = 270))+ 

labs(title="Prob. that the Vuong test (Poisson vs. ZIP) rejects the null of no zero-inflation")

pdf("~/plots/selectionvuongP.pdf", width=9.5, height=7)

p +  geom_text(
    size    = 4,
    data    = resultsmatsmall,
    mapping = aes(x = Inf, y = Inf, label = as.numeric(lab)),
    hjust   = 1.05,
    vjust   = 1.5, col="black"
  ) + 
  scale_y_continuous(sec.axis = sec_axis(trans ~ ., labels = NULL,  name = expression(" " %<-% "direction of increasing overdispersion" %<-% "\n"), t))  + theme(legend.position="bottom") +
  scale_x_continuous(trans='log10', breaks = c(50, 250, 1000, 10000), sec.axis = dup_axis(trans ~ ., labels = NULL, name = expression(" " %->% "direction of increasing zero-inflation" %->% "\n"))) 
dev.off()


#######  selectionvuongNB
p <-ggplot(resultsmatsmall, aes(x = n, y = selectionvuongNB, group = g, colour=as.factor(beta_0))) + geom_point() + geom_line()  + 

facet_grid(phi ~omega, labeller=label_parsed) + 
labs(x = paste("sample size (n)"), y = expression("probability of rejecting the null")) + 
coord_cartesian(xlim = c(40, 15000), ylim = c(0, 1)) + 
guides(color=guide_legend(title = expression(paste("intercept parameter (", beta[0], ") : "))))  +  
theme(axis.title = element_text(size = 15), axis.text.x = element_text(vjust = 3), legend.text=element_text(size=15), axis.ticks = element_blank(), legend.title=element_text(size=16), axis.text = element_text(size=10.5),   panel.grid.minor.x = element_blank(), plot.title=element_text(size=17), strip.text.x = element_text(size = 13, colour = "black", angle = 0), strip.text.y = element_text(size = 12, colour = "black", angle = 270))+ 

labs(title="Prob. that the Vuong test (NB vs. ZINB) rejects the null of no zero-inflation")

pdf("~/plots/selectionvuongP.pdf", width=9.5, height=7)

p +  geom_text(
    size    = 4,
    data    = resultsmatsmall,
    mapping = aes(x = Inf, y = Inf, label = as.numeric(lab)),
    hjust   = 1.05,
    vjust   = 1.5, col="black"
  ) + 
  scale_y_continuous(sec.axis = sec_axis(trans ~ ., labels = NULL,  name = expression(" " %<-% "direction of increasing overdispersion" %<-% "\n"), t))  + theme(legend.position="bottom") +
  scale_x_continuous(trans='log10', breaks = c(50, 250, 1000, 10000), sec.axis = dup_axis(trans ~ ., labels = NULL, name = expression(" " %->% "direction of increasing zero-inflation" %->% "\n"))) 

dev.off()


#######  Poisson_pval
p <-ggplot(resultsmatsmall, aes(x = n, y = Poisson_pval, group = g, colour=as.factor(beta_0))) + geom_point() + geom_line()  + 

facet_grid(phi ~omega, labeller=label_parsed) + 
labs(x = paste("sample size (n)"), y = expression("probability of p" < alpha)) + 
coord_cartesian(xlim = c(40, 15000), ylim = c(0, 0.2)) + 
guides(color=guide_legend(title = expression(paste("intercept parameter (", beta[0], ") : "))))  +  
theme(axis.title = element_text(size = 15), axis.text.x = element_text(vjust = 3), legend.text=element_text(size=15), axis.ticks = element_blank(), legend.title=element_text(size=16), axis.text = element_text(size=10.5),   panel.grid.minor.x = element_blank(), plot.title=element_text(size=17), strip.text.x = element_text(size = 13, colour = "black", angle = 0), strip.text.y = element_text(size = 12, colour = "black", angle = 270))+ 

geom_hline(yintercept = alpha_sig_level) + labs(title=expression(paste("Probability that the Poisson model rejects the null of ", beta[X], " = 0")))


pdf("~/plots/Poisson_pval.pdf", width=9.5, height=7)

p +  geom_text(
    size    = 4,
    data    = resultsmatsmall,
    mapping = aes(x = Inf, y = Inf, label = as.numeric(lab)),
    hjust   = 1.05,
    vjust   = 1.5, col="black"
  ) + 
  scale_y_continuous(breaks = c(0, 0.05, 0.10, 0.15), labels=c("0", "0.05", "0.10", "0.15"), sec.axis = sec_axis(trans ~ ., labels = NULL,  name = expression(" " %<-% "direction of increasing overdispersion" %<-% "\n"), t))  + theme(legend.position="bottom") +
  scale_x_continuous(trans='log10', breaks = c(50, 250, 1000, 10000), sec.axis = dup_axis(trans ~ ., labels = NULL, name = expression(" " %->% "direction of increasing zero-inflation" %->% "\n"))) 

dev.off()



#######  zip_pval
p <-ggplot(resultsmatsmall, aes(x = n, y = zip_pval, group = g, colour=as.factor(beta_0))) + geom_point() + geom_line()  + 

facet_grid(phi ~omega, labeller=label_parsed) + 
labs(x = paste("sample size (n)"), y = expression("probability of p" < alpha)) + 
coord_cartesian(xlim = c(40, 15000), ylim = c(0, 0.2)) + 
guides(color=guide_legend(title = expression(paste("intercept parameter (", beta[0], ") : "))))  +  
theme(axis.title = element_text(size = 15), axis.text.x = element_text(vjust = 3), legend.text=element_text(size=15), axis.ticks = element_blank(), legend.title=element_text(size=16), axis.text = element_text(size=10.5),   panel.grid.minor.x = element_blank(), plot.title=element_text(size=17), strip.text.x = element_text(size = 13, colour = "black", angle = 0), strip.text.y = element_text(size = 12, colour = "black", angle = 270))+ 

geom_hline(yintercept = alpha_sig_level) + labs(title=expression(paste("Probability that the ZIP model rejects the null of ", beta[X], " = ", gamma[X], " = 0")))


pdf("~/plots/zip_pval.pdf", width=9.5, height=7)

p +  geom_text(
    size    = 4,
    data    = resultsmatsmall,
    mapping = aes(x = Inf, y = Inf, label = as.numeric(lab)),
    hjust   = 1.05,
    vjust   = 1.5, col="black"
  ) + 
  scale_y_continuous(breaks = c(0, 0.05, 0.10, 0.15), labels=c("0", "0.05", "0.10", "0.15"), sec.axis = sec_axis(trans ~ ., labels = NULL,  name = expression(" " %<-% "direction of increasing overdispersion" %<-% "\n"), t))  + theme(legend.position="bottom") +
  scale_x_continuous(trans='log10', breaks = c(50, 250, 1000, 10000), sec.axis = dup_axis(trans ~ ., labels = NULL, name = expression(" " %->% "direction of increasing zero-inflation" %->% "\n"))) 

dev.off()




#######  nb_pval
p <-ggplot(resultsmatsmall, aes(x = n, y = nb_pval, group = g, colour=as.factor(beta_0))) + geom_point() + geom_line()  + 

facet_grid(phi ~omega, labeller=label_parsed) + 
labs(x = paste("sample size (n)"), y = expression("probability of p" < alpha)) + 
coord_cartesian(xlim = c(40, 15000), ylim = c(0, 0.2)) + 
guides(color=guide_legend(title = expression(paste("intercept parameter (", beta[0], ") : "))))  +  
theme(axis.title = element_text(size = 15), axis.text.x = element_text(vjust = 3), legend.text=element_text(size=15), axis.ticks = element_blank(), legend.title=element_text(size=16), axis.text = element_text(size=10.5),   panel.grid.minor.x = element_blank(), plot.title=element_text(size=17), strip.text.x = element_text(size = 13, colour = "black", angle = 0), strip.text.y = element_text(size = 12, colour = "black", angle = 270))+ 

geom_hline(yintercept = alpha_sig_level) + labs(title=expression(paste("Probability that the NB model rejects the null of ", beta[X], " = 0")))

pdf("~/plots/nb_pval.pdf", width=9.5, height=7)

p +  geom_text(
    size    = 4,
    data    = resultsmatsmall,
    mapping = aes(x = Inf, y = Inf, label = as.numeric(lab)),
    hjust   = 1.05,
    vjust   = 1.5, col="black"
  ) + 
  scale_y_continuous(breaks = c(0, 0.05, 0.10, 0.15), labels=c("0", "0.05", "0.10", "0.15"), sec.axis = sec_axis(trans ~ ., labels = NULL,  name = expression(" " %<-% "direction of increasing overdispersion" %<-% "\n"), t))  + theme(legend.position="bottom") +
  scale_x_continuous(trans='log10', breaks = c(50, 250, 1000, 10000), sec.axis = dup_axis(trans ~ ., labels = NULL, name = expression(" " %->% "direction of increasing zero-inflation" %->% "\n"))) 

dev.off()  


#######  zinb_pval
p <-ggplot(resultsmatsmall, aes(x = n, y = zinb_pval, group = g, colour=as.factor(beta_0))) + geom_point() + geom_line()  + 

facet_grid(phi ~omega, labeller=label_parsed) + 
labs(x = paste("sample size (n)"), y = expression("probability of p" < alpha)) + 
coord_cartesian(xlim = c(40, 15000), ylim = c(0, 0.2)) + 
guides(color=guide_legend(title = expression(paste("intercept parameter (", beta[0], ") : "))))  +  
theme(axis.title = element_text(size = 15), axis.text.x = element_text(vjust = 3), legend.text=element_text(size=15), axis.ticks = element_blank(), legend.title=element_text(size=16), axis.text = element_text(size=10.5),   panel.grid.minor.x = element_blank(), plot.title=element_text(size=17), strip.text.x = element_text(size = 13, colour = "black", angle = 0), strip.text.y = element_text(size = 12, colour = "black", angle = 270))+ 

geom_hline(yintercept = alpha_sig_level) + labs(title=expression(paste("Probability that the ZINB model rejects the null of ", beta[X], " = ", gamma[X], " = 0")))

pdf("~/plots/zinb_pval.pdf", width=9.5, height=7)

p +  geom_text(
    size    = 4,
    data    = resultsmatsmall,
    mapping = aes(x = Inf, y = Inf, label = as.numeric(lab)),
    hjust   = 1.05,
    vjust   = 1.5, col="black"
  ) + 
  scale_y_continuous(breaks = c(0, 0.05, 0.10, 0.15), labels=c("0", "0.05", "0.10", "0.15"), sec.axis = sec_axis(trans ~ ., labels = NULL,  name = expression(" " %<-% "direction of increasing overdispersion" %<-% "\n"), t))  + theme(legend.position="bottom") +
  scale_x_continuous(trans='log10', breaks = c(50, 250, 1000, 10000), sec.axis = dup_axis(trans ~ ., labels = NULL, name = expression(" " %->% "direction of increasing zero-inflation" %->% "\n"))) 

dev.off() #######  correct
p <-ggplot(resultsmatsmall, aes(x = n, y = correct, group = g, colour=as.factor(beta_0))) + geom_point() + geom_line()  + 

facet_grid(phi ~omega, labeller=label_parsed) + 
labs(x = paste("sample size (n)"), y = expression("probability of p" < alpha)) + 
coord_cartesian(xlim = c(40, 15000), ylim = c(0, 0.1)) + 
guides(color=guide_legend(title = expression(paste("intercept parameter (", beta[0], ") : "))))  +  
theme(axis.title = element_text(size = 15), axis.text.x = element_text(vjust = 3), legend.text=element_text(size=15), axis.ticks = element_blank(), legend.title=element_text(size=16), axis.text = element_text(size=10.5),   panel.grid.minor.x = element_blank(), plot.title=element_text(size=17), strip.text.x = element_text(size = 13, colour = "black", angle = 0), strip.text.y = element_text(size = 12, colour = "black", angle = 270))+ 

geom_hline(yintercept = alpha_sig_level) + labs(title="Type 1 error obtained with correct model")

pdf("~/plots/correct.pdf", width=9.5, height=7)


p +  geom_text(
    size    = 4,
    data    = resultsmatsmall,
    mapping = aes(x = Inf, y = Inf, label = as.numeric(lab)),
    hjust   = 1.05,
    vjust   = 1.5, col="black"
  ) + 
  scale_y_continuous(breaks = c(0, 0.025, 0.05, 0.075), labels=c("0", "0.025", "0.050", "0.075"), sec.axis = sec_axis(trans ~ ., labels = NULL,  name = expression(" " %<-% "direction of increasing overdispersion" %<-% "\n"), t))  + theme(legend.position="bottom") +
  scale_x_continuous(trans='log10', breaks = c(50, 250, 1000, 10000), sec.axis = dup_axis(trans ~ ., labels = NULL, name = expression(" " %->% "direction of increasing zero-inflation" %->% "\n"))) 

dev.off()
#######  correct_selection
p <-ggplot(resultsmatsmall, aes(x = n, y = correct_selection, group = g, colour=as.factor(beta_0))) + geom_point() + geom_line()  + 

facet_grid(phi ~omega, labeller=label_parsed) + 
labs(x = paste("sample size (n)"), y = expression("probability of selecting correct model")) + 
coord_cartesian(xlim = c(40, 15000), ylim = c(0, 1)) + 
guides(color=guide_legend(title = expression(paste("intercept parameter (", beta[0], ") : "))))  +  
theme(axis.title = element_text(size = 15), axis.text.x = element_text(vjust = 3), legend.text=element_text(size=15), axis.ticks = element_blank(), legend.title=element_text(size=16), axis.text = element_text(size=10.5),   panel.grid.minor.x = element_blank(), plot.title=element_text(size=17), strip.text.x = element_text(size = 13, colour = "black", angle = 0), strip.text.y = element_text(size = 12, colour = "black", angle = 270))+ 

 labs(title="Probability of selecting the correct model following the 7-step testing scheme")

pdf("~/plots/correct_selection.pdf", width=9.5, height=7)


p +  geom_text(
    size    = 4,
    data    = resultsmatsmall,
    mapping = aes(x = Inf, y = Inf, label = as.numeric(lab)),
    hjust   = 1.05,
    vjust   = 1.5, col="black"
  ) + 
  scale_y_continuous(sec.axis = sec_axis(trans ~ ., labels = NULL,  name = expression(" " %<-% "direction of increasing overdispersion" %<-% "\n"), t))  + theme(legend.position="bottom") +
  scale_x_continuous(trans='log10', breaks = c(50, 250, 1000, 10000), sec.axis = dup_axis(trans ~ ., labels = NULL, name = expression(" " %->% "direction of increasing zero-inflation" %->% "\n"))) 

dev.off()


#######  type1
p <-ggplot(resultsmatsmall, aes(x = n, y = type1, group = g, colour=as.factor(beta_0))) + geom_point() + geom_line()  + 

facet_grid(phi ~omega, labeller=label_parsed) + 
labs(x = paste("sample size (n)"), y = expression("probability of p" < alpha)) + 
coord_cartesian(xlim = c(40, 15000), ylim = c(0, 0.1)) + 
guides(color=guide_legend(title = expression(paste("intercept parameter (", beta[0], ") : "))))  +  
theme(axis.title = element_text(size = 15), axis.text.x = element_text(vjust = 3), legend.text=element_text(size=15), axis.ticks = element_blank(), legend.title=element_text(size=16), axis.text = element_text(size=10.5),   panel.grid.minor.x = element_blank(), plot.title=element_text(size=17), strip.text.x = element_text(size = 13, colour = "black", angle = 0), strip.text.y = element_text(size = 12, colour = "black", angle = 270))+ 

geom_hline(yintercept = alpha_sig_level) + labs(title="Type 1 error obtained following the 7-step testing scheme")

pdf("~/plots/type1.pdf", width=9.5, height=7)


p +  geom_text(
    size    = 4,
    data    = resultsmatsmall,
    mapping = aes(x = Inf, y = Inf, label = as.numeric(lab)),
    hjust   = 1.05,
    vjust   = 1.5, col="black"
  ) + 
  scale_y_continuous(breaks = c(0, 0.025, 0.05, 0.075), labels=c("0", "0.025", "0.050", "0.075"), sec.axis = sec_axis(trans ~ ., labels = NULL,  name = expression(" " %<-% "direction of increasing overdispersion" %<-% "\n"), t))  + theme(legend.position="bottom") +
  scale_x_continuous(trans='log10', breaks = c(50, 250, 1000, 10000), sec.axis = dup_axis(trans ~ ., labels = NULL, name = expression(" " %->% "direction of increasing zero-inflation" %->% "\n"))) 

dev.off()


#######  correct_selectionAIC
p <-ggplot(resultsmatsmall, aes(x = n, y = correct_selectionAIC, group = g, colour=as.factor(beta_0))) + geom_point() + geom_line()  + 

facet_grid(phi ~omega, labeller=label_parsed) + 
labs(x = paste("sample size (n)"), y = expression("probability of selecting correct model")) + 
coord_cartesian(xlim = c(40, 15000), ylim = c(0, 1)) + 
guides(color=guide_legend(title = expression(paste("intercept parameter (", beta[0], ") : "))))  +  
theme(axis.title = element_text(size = 15), axis.text.x = element_text(vjust = 3), legend.text=element_text(size=15), axis.ticks = element_blank(), legend.title=element_text(size=16), axis.text = element_text(size=10.5),   panel.grid.minor.x = element_blank(), plot.title=element_text(size=17), strip.text.x = element_text(size = 13, colour = "black", angle = 0), strip.text.y = element_text(size = 12, colour = "black", angle = 270))+ 

labs(title="Probability that the correct model is the one with the lowest AIC")

pdf("~/plots/correct_selectionAIC.pdf", width=9.5, height=7)

p +  geom_text(
    size    = 4,
    data    = resultsmatsmall,
    mapping = aes(x = Inf, y = Inf, label = as.numeric(lab)),
    hjust   = 1.05,
    vjust   = 1.5, col="black"
  ) + 
  scale_y_continuous(sec.axis = sec_axis(trans ~ ., labels = NULL,  name = expression(" " %<-% "direction of increasing overdispersion" %<-% "\n"), t))  + theme(legend.position="bottom") +
  scale_x_continuous(trans='log10', breaks = c(50, 250, 1000, 10000), sec.axis = dup_axis(trans ~ ., labels = NULL, name = expression(" " %->% "direction of increasing zero-inflation" %->% "\n"))) 
dev.off()


#######  correct_selectionBIC
p <-ggplot(resultsmatsmall, aes(x = n, y = correct_selectionBIC, group = g, colour=as.factor(beta_0))) + geom_point() + geom_line()  + 

facet_grid(phi ~omega, labeller=label_parsed) + 
labs(x = paste("sample size (n)"), y = expression("probability of selecting correct model")) + 
coord_cartesian(xlim = c(40, 15000), ylim = c(0, 1)) + 
guides(color=guide_legend(title = expression(paste("intercept parameter (", beta[0], ") : "))))  +  
theme(axis.title = element_text(size = 15), axis.text.x = element_text(vjust = 3), legend.text=element_text(size=15), axis.ticks = element_blank(), legend.title=element_text(size=16), axis.text = element_text(size=10.5),   panel.grid.minor.x = element_blank(), plot.title=element_text(size=17), strip.text.x = element_text(size = 13, colour = "black", angle = 0), strip.text.y = element_text(size = 12, colour = "black", angle = 270))+ 

labs(title="Probability that the correct model is the one with the lowest BIC")

pdf("~/plots/correct_selectionBIC.pdf", width=9.5, height=7)

p +  geom_text(
    size    = 4,
    data    = resultsmatsmall,
    mapping = aes(x = Inf, y = Inf, label = as.numeric(lab)),
    hjust   = 1.05,
    vjust   = 1.5, col="black"
  ) + 
  scale_y_continuous(sec.axis = sec_axis(trans ~ ., labels = NULL,  name = expression(" " %<-% "direction of increasing overdispersion" %<-% "\n"), t))  + theme(legend.position="bottom") +
  scale_x_continuous(trans='log10', breaks = c(50, 250, 1000, 10000), sec.axis = dup_axis(trans ~ ., labels = NULL, name = expression(" " %->% "direction of increasing zero-inflation" %->% "\n"))) 

dev.off()


#######  correct_selectionAICc
p <-ggplot(resultsmatsmall, aes(x = n, y = correct_selectionAICc, group = g, colour=as.factor(beta_0))) + geom_point() + geom_line()  + 

facet_grid(phi ~omega, labeller=label_parsed) + 
labs(x = paste("sample size (n)"), y = expression("probability of selecting correct model")) + 
coord_cartesian(xlim = c(40, 15000), ylim = c(0, 1)) + 
guides(color=guide_legend(title = expression(paste("intercept parameter (", beta[0], ") : "))))  +  
theme(axis.title = element_text(size = 15), axis.text.x = element_text(vjust = 3), legend.text=element_text(size=15), axis.ticks = element_blank(), legend.title=element_text(size=16), axis.text = element_text(size=10.5),   panel.grid.minor.x = element_blank(), plot.title=element_text(size=17), strip.text.x = element_text(size = 13, colour = "black", angle = 0), strip.text.y = element_text(size = 12, colour = "black", angle = 270))+ 

 labs(title="Probability that the correct model is the one with the lowest AICc")

pdf("~/plots/correct_selectionAICc.pdf", width=9.5, height=7)


p +  geom_text(
    size    = 4,
    data    = resultsmatsmall,
    mapping = aes(x = Inf, y = Inf, label = as.numeric(lab)),
    hjust   = 1.05,
    vjust   = 1.5, col="black"
  ) + 
  scale_y_continuous(sec.axis = sec_axis(trans ~ ., labels = NULL,  name = expression(" " %<-% "direction of increasing overdispersion" %<-% "\n"), t))  + theme(legend.position="bottom") +
  scale_x_continuous(trans='log10', breaks = c(50, 250, 1000, 10000), sec.axis = dup_axis(trans ~ ., labels = NULL, name = expression(" " %->% "direction of increasing zero-inflation" %->% "\n"))) 

dev.off()


#######  typeAIC
p <-ggplot(resultsmatsmall, aes(x = n, y = typeAIC, group = g, colour=as.factor(beta_0))) + geom_point() + geom_line()  + 

facet_grid(phi ~omega, labeller=label_parsed) + 
labs(x = paste("sample size (n)"), y = expression("probability of p" < alpha)) + 
coord_cartesian(xlim = c(40, 15000), ylim = c(0, 0.1)) + 
guides(color=guide_legend(title = expression(paste("intercept parameter (", beta[0], ") : "))))  +  
theme(axis.title = element_text(size = 15), axis.text.x = element_text(vjust = 3), legend.text=element_text(size=15), axis.ticks = element_blank(), legend.title=element_text(size=16), axis.text = element_text(size=10.5),   panel.grid.minor.x = element_blank(), plot.title=element_text(size=17), strip.text.x = element_text(size = 13, colour = "black", angle = 0), strip.text.y = element_text(size = 12, colour = "black", angle = 270))+ 

geom_hline(yintercept = alpha_sig_level) + labs(title="Type  1  error  obtained from model with the lowest AIC")


pdf("~/plots/typeAIC.pdf", width=9.5, height=7)

p +  geom_text(
    size    = 4,
    data    = resultsmatsmall,
    mapping = aes(x = Inf, y = Inf, label = as.numeric(lab)),
    hjust   = 1.05,
    vjust   = 1.5, col="black"
  ) + 
  scale_y_continuous(breaks = c(0, 0.025, 0.05, 0.075), labels=c("0", "0.025", "0.050", "0.075"), sec.axis = sec_axis(trans ~ ., labels = NULL,  name = expression(" " %<-% "direction of increasing overdispersion" %<-% "\n"), t))  + theme(legend.position="bottom") +
  scale_x_continuous(trans='log10', breaks = c(50, 250, 1000, 10000), sec.axis = dup_axis(trans ~ ., labels = NULL, name = expression(" " %->% "direction of increasing zero-inflation" %->% "\n"))) 
dev.off()



#######  typeBIC
p <-ggplot(resultsmatsmall, aes(x = n, y = typeBIC, group = g, colour=as.factor(beta_0))) + geom_point() + geom_line()  + 

facet_grid(phi ~omega, labeller=label_parsed) + 
labs(x = paste("sample size (n)"), y = expression("probability of p" < alpha)) + 
coord_cartesian(xlim = c(40, 15000), ylim = c(0, 0.1)) + 
guides(color=guide_legend(title = expression(paste("intercept parameter (", beta[0], ") : "))))  +  
theme(axis.title = element_text(size = 15), axis.text.x = element_text(vjust = 3), legend.text=element_text(size=15), axis.ticks = element_blank(), legend.title=element_text(size=16), axis.text = element_text(size=10.5),   panel.grid.minor.x = element_blank(), plot.title=element_text(size=17), strip.text.x = element_text(size = 13, colour = "black", angle = 0), strip.text.y = element_text(size = 12, colour = "black", angle = 270))+ 

geom_hline(yintercept = alpha_sig_level) + labs(title="Type  1  error  obtained from model with the lowest BIC")
pdf("~/plots/typeBIC.pdf", width=9.5, height=7)

p +  geom_text(
    size    = 4,
    data    = resultsmatsmall,
    mapping = aes(x = Inf, y = Inf, label = as.numeric(lab)),
    hjust   = 1.05,
    vjust   = 1.5, col="black"
  ) + 
  scale_y_continuous(breaks = c(0, 0.025, 0.05, 0.075), labels=c("0", "0.025", "0.050", "0.075"), sec.axis = sec_axis(trans ~ ., labels = NULL,  name = expression(" " %<-% "direction of increasing overdispersion" %<-% "\n"), t))  + theme(legend.position="bottom") +
  scale_x_continuous(trans='log10', breaks = c(50, 250, 1000, 10000), sec.axis = dup_axis(trans ~ ., labels = NULL, name = expression(" " %->% "direction of increasing zero-inflation" %->% "\n"))) 

dev.off() 



#######  typeAICc
p <-ggplot(resultsmatsmall, aes(x = n, y = typeAICc, group = g, colour=as.factor(beta_0))) + geom_point() + geom_line()  + 

facet_grid(phi ~omega, labeller=label_parsed) + 
labs(x = paste("sample size (n)"), y = expression("probability of p" < alpha)) + 
coord_cartesian(xlim = c(40, 15000), ylim = c(0, 0.1)) + 
guides(color=guide_legend(title = expression(paste("intercept parameter (", beta[0], ") : "))))  +  
theme(axis.title = element_text(size = 15), axis.text.x = element_text(vjust = 3), legend.text=element_text(size=15), axis.ticks = element_blank(), legend.title=element_text(size=16), axis.text = element_text(size=10.5),   panel.grid.minor.x = element_blank(), plot.title=element_text(size=17), strip.text.x = element_text(size = 13, colour = "black", angle = 0), strip.text.y = element_text(size = 12, colour = "black", angle = 270))+ 

geom_hline(yintercept = alpha_sig_level) + labs(title="Type  1  error  obtained from model with the lowest AICc")

pdf("~/plots/typeAICc.pdf", width=9.5, height=7)

p +  geom_text(
    size    = 4,
    data    = resultsmatsmall,
    mapping = aes(x = Inf, y = Inf, label = as.numeric(lab)),
    hjust   = 1.05,
    vjust   = 1.5, col="black"
  ) + 
  scale_y_continuous(breaks = c(0, 0.025, 0.05, 0.075), labels=c("0", "0.025", "0.050", "0.075"), sec.axis = sec_axis(trans ~ ., labels = NULL,  name = expression(" " %<-% "direction of increasing overdispersion" %<-% "\n"), t))  + theme(legend.position="bottom") +
  scale_x_continuous(trans='log10', breaks = c(50, 250, 1000, 10000), sec.axis = dup_axis(trans ~ ., labels = NULL, name = expression(" " %->% "direction of increasing zero-inflation" %->% "\n"))) 

dev.off()


#######  selectionLRT
p <-ggplot(resultsmatsmall, aes(x = n, y = selectionLRT, group = g, colour=as.factor(beta_0))) + geom_point() + geom_line()  + 

facet_grid(phi ~omega, labeller=label_parsed) + 
labs(x = paste("sample size (n)"), y = expression("probability of rejecting the null")) + 
coord_cartesian(xlim = c(40, 15000), ylim = c(0, 1)) + 
guides(color=guide_legend(title = expression(paste("intercept parameter (", beta[0], ") : "))))  +  
theme(axis.title = element_text(size = 15), axis.text.x = element_text(vjust = 3), legend.text=element_text(size=15), axis.ticks = element_blank(), legend.title=element_text(size=16), axis.text = element_text(size=10.5),   panel.grid.minor.x = element_blank(), plot.title=element_text(size=17), strip.text.x = element_text(size = 13, colour = "black", angle = 0), strip.text.y = element_text(size = 12, colour = "black", angle = 270))+ 

labs(title="Probability  that  the  D&L  test  rejects  the  null  of no overdispersion")

pdf("~/plots/selectionLRT.pdf", width=9.5, height=7)

p +  geom_text(
    size    = 4,
    data    = resultsmatsmall,
    mapping = aes(x = Inf, y = Inf, label = as.numeric(lab)),
    hjust   = 1.05,
    vjust   = 1.5, col="black"
  ) + 
  scale_y_continuous(sec.axis = sec_axis(trans ~ ., labels = NULL,  name = expression(" " %<-% "direction of increasing overdispersion" %<-% "\n"), t))  + theme(legend.position="bottom") +
  scale_x_continuous(trans='log10', breaks = c(50, 250, 1000, 10000), sec.axis = dup_axis(trans ~ ., labels = NULL, name = expression(" " %->% "direction of increasing zero-inflation" %->% "\n"))) 

dev.off()


#######  selectionvuongP
p <-ggplot(resultsmatsmall, aes(x = n, y = selectionvuongP, group = g, colour=as.factor(beta_0))) + geom_point() + geom_line()  + 

facet_grid(phi ~omega, labeller=label_parsed) + 
labs(x = paste("sample size (n)"), y = expression("probability of rejecting the null")) + 
coord_cartesian(xlim = c(40, 15000), ylim = c(0, 1)) + 
guides(color=guide_legend(title = expression(paste("intercept parameter (", beta[0], ") : "))))  +  
theme(axis.title = element_text(size = 15), axis.text.x = element_text(vjust = 3), legend.text=element_text(size=15), axis.ticks = element_blank(), legend.title=element_text(size=16), axis.text = element_text(size=10.5),   panel.grid.minor.x = element_blank(), plot.title=element_text(size=17), strip.text.x = element_text(size = 13, colour = "black", angle = 0), strip.text.y = element_text(size = 12, colour = "black", angle = 270))+ 

labs(title="Prob. that the Vuong test (Poisson vs. ZIP) rejects the null of no zero-inflation")

pdf("~/plots/selectionvuongP.pdf", width=9.5, height=7)

p +  geom_text(
    size    = 4,
    data    = resultsmatsmall,
    mapping = aes(x = Inf, y = Inf, label = as.numeric(lab)),
    hjust   = 1.05,
    vjust   = 1.5, col="black"
  ) + 
  scale_y_continuous(sec.axis = sec_axis(trans ~ ., labels = NULL,  name = expression(" " %<-% "direction of increasing overdispersion" %<-% "\n"), t))  + theme(legend.position="bottom") +
  scale_x_continuous(trans='log10', breaks = c(50, 250, 1000, 10000), sec.axis = dup_axis(trans ~ ., labels = NULL, name = expression(" " %->% "direction of increasing zero-inflation" %->% "\n"))) 
dev.off()


#######  selectionvuongNB
p <-ggplot(resultsmatsmall, aes(x = n, y = selectionvuongNB, group = g, colour=as.factor(beta_0))) + geom_point() + geom_line()  + 

facet_grid(phi ~omega, labeller=label_parsed) + 
labs(x = paste("sample size (n)"), y = expression("probability of rejecting the null")) + 
coord_cartesian(xlim = c(40, 15000), ylim = c(0, 1)) + 
guides(color=guide_legend(title = expression(paste("intercept parameter (", beta[0], ") : "))))  +  
theme(axis.title = element_text(size = 15), axis.text.x = element_text(vjust = 3), legend.text=element_text(size=15), axis.ticks = element_blank(), legend.title=element_text(size=16), axis.text = element_text(size=10.5),   panel.grid.minor.x = element_blank(), plot.title=element_text(size=17), strip.text.x = element_text(size = 13, colour = "black", angle = 0), strip.text.y = element_text(size = 12, colour = "black", angle = 270))+ 

labs(title="Prob. that the Vuong test (NB vs. ZINB) rejects the null of no zero-inflation")

pdf("~/plots/selectionvuongNB.pdf", width=9.5, height=7)

p +  geom_text(
    size    = 4,
    data    = resultsmatsmall,
    mapping = aes(x = Inf, y = Inf, label = as.numeric(lab)),
    hjust   = 1.05,
    vjust   = 1.5, col="black"
  ) + 
  scale_y_continuous(sec.axis = sec_axis(trans ~ ., labels = NULL,  name = expression(" " %<-% "direction of increasing overdispersion" %<-% "\n"), t))  + theme(legend.position="bottom") +
  scale_x_continuous(trans='log10', breaks = c(50, 250, 1000, 10000), sec.axis = dup_axis(trans ~ ., labels = NULL, name = expression(" " %->% "direction of increasing zero-inflation" %->% "\n"))) 

dev.off()


#######  Poisson_pval
p <-ggplot(resultsmatsmall, aes(x = n, y = poisson_pval, group = g, colour=as.factor(beta_0))) + geom_point() + geom_line()  + 

facet_grid(phi ~omega, labeller=label_parsed) + 
labs(x = paste("sample size (n)"), y = expression("probability of p" < alpha)) + 
coord_cartesian(xlim = c(40, 15000), ylim = c(0, 0.2)) + 
guides(color=guide_legend(title = expression(paste("intercept parameter (", beta[0], ") : "))))  +  
theme(axis.title = element_text(size = 15), axis.text.x = element_text(vjust = 3), legend.text=element_text(size=15), axis.ticks = element_blank(), legend.title=element_text(size=16), axis.text = element_text(size=10.5),   panel.grid.minor.x = element_blank(), plot.title=element_text(size=17), strip.text.x = element_text(size = 13, colour = "black", angle = 0), strip.text.y = element_text(size = 12, colour = "black", angle = 270))+ 

geom_hline(yintercept = alpha_sig_level) + labs(title=expression(paste("Probability that the Poisson model rejects the null of ", beta[X], " = 0")))


pdf("~/plots/Poisson_pval.pdf", width=9.5, height=7)

p +  geom_text(
    size    = 4,
    data    = resultsmatsmall,
    mapping = aes(x = Inf, y = Inf, label = as.numeric(lab)),
    hjust   = 1.05,
    vjust   = 1.5, col="black"
  ) + 
  scale_y_continuous(breaks = c(0, 0.05, 0.10, 0.15), labels=c("0", "0.05", "0.10", "0.15"), sec.axis = sec_axis(trans ~ ., labels = NULL,  name = expression(" " %<-% "direction of increasing overdispersion" %<-% "\n"), t))  + theme(legend.position="bottom") +
  scale_x_continuous(trans='log10', breaks = c(50, 250, 1000, 10000), sec.axis = dup_axis(trans ~ ., labels = NULL, name = expression(" " %->% "direction of increasing zero-inflation" %->% "\n"))) 

dev.off()



#######  zip_pval
p <-ggplot(resultsmatsmall, aes(x = n, y = zip_pval, group = g, colour=as.factor(beta_0))) + geom_point() + geom_line()  + 

facet_grid(phi ~omega, labeller=label_parsed) + 
labs(x = paste("sample size (n)"), y = expression("probability of p" < alpha)) + 
coord_cartesian(xlim = c(40, 15000), ylim = c(0, 0.2)) + 
guides(color=guide_legend(title = expression(paste("intercept parameter (", beta[0], ") : "))))  +  
theme(axis.title = element_text(size = 15), axis.text.x = element_text(vjust = 3), legend.text=element_text(size=15), axis.ticks = element_blank(), legend.title=element_text(size=16), axis.text = element_text(size=10.5),   panel.grid.minor.x = element_blank(), plot.title=element_text(size=17), strip.text.x = element_text(size = 13, colour = "black", angle = 0), strip.text.y = element_text(size = 12, colour = "black", angle = 270))+ 

geom_hline(yintercept = alpha_sig_level) + labs(title=expression(paste("Probability that the ZIP model rejects the null of ", beta[X], " = ", gamma[X], " = 0")))


pdf("~/plots/zip_pval.pdf", width=9.5, height=7)

p +  geom_text(
    size    = 4,
    data    = resultsmatsmall,
    mapping = aes(x = Inf, y = Inf, label = as.numeric(lab)),
    hjust   = 1.05,
    vjust   = 1.5, col="black"
  ) + 
  scale_y_continuous(breaks = c(0, 0.05, 0.10, 0.15), labels=c("0", "0.05", "0.10", "0.15"), sec.axis = sec_axis(trans ~ ., labels = NULL,  name = expression(" " %<-% "direction of increasing overdispersion" %<-% "\n"), t))  + theme(legend.position="bottom") +
  scale_x_continuous(trans='log10', breaks = c(50, 250, 1000, 10000), sec.axis = dup_axis(trans ~ ., labels = NULL, name = expression(" " %->% "direction of increasing zero-inflation" %->% "\n"))) 

dev.off()




#######  nb_pval
p <-ggplot(resultsmatsmall, aes(x = n, y = nb_pval, group = g, colour=as.factor(beta_0))) + geom_point() + geom_line()  + 

facet_grid(phi ~omega, labeller=label_parsed) + 
labs(x = paste("sample size (n)"), y = expression("probability of p" < alpha)) + 
coord_cartesian(xlim = c(40, 15000), ylim = c(0, 0.2)) + 
guides(color=guide_legend(title = expression(paste("intercept parameter (", beta[0], ") : "))))  +  
theme(axis.title = element_text(size = 15), axis.text.x = element_text(vjust = 3), legend.text=element_text(size=15), axis.ticks = element_blank(), legend.title=element_text(size=16), axis.text = element_text(size=10.5),   panel.grid.minor.x = element_blank(), plot.title=element_text(size=17), strip.text.x = element_text(size = 13, colour = "black", angle = 0), strip.text.y = element_text(size = 12, colour = "black", angle = 270))+ 

geom_hline(yintercept = alpha_sig_level) + labs(title=expression(paste("Probability that the NB model rejects the null of ", beta[X], " = 0")))

pdf("~/plots/nb_pval.pdf", width=9.5, height=7)

p +  geom_text(
    size    = 4,
    data    = resultsmatsmall,
    mapping = aes(x = Inf, y = Inf, label = as.numeric(lab)),
    hjust   = 1.05,
    vjust   = 1.5, col="black"
  ) + 
  scale_y_continuous(breaks = c(0, 0.05, 0.10, 0.15), labels=c("0", "0.05", "0.10", "0.15"), sec.axis = sec_axis(trans ~ ., labels = NULL,  name = expression(" " %<-% "direction of increasing overdispersion" %<-% "\n"), t))  + theme(legend.position="bottom") +
  scale_x_continuous(trans='log10', breaks = c(50, 250, 1000, 10000), sec.axis = dup_axis(trans ~ ., labels = NULL, name = expression(" " %->% "direction of increasing zero-inflation" %->% "\n"))) 

dev.off()  


#######  zinb_pval
p <-ggplot(resultsmatsmall, aes(x = n, y = zinb_pval, group = g, colour=as.factor(beta_0))) + geom_point() + geom_line()  + 

facet_grid(phi ~omega, labeller=label_parsed) + 
labs(x = paste("sample size (n)"), y = expression("probability of p" < alpha)) + 
coord_cartesian(xlim = c(40, 15000), ylim = c(0, 0.2)) + 
guides(color=guide_legend(title = expression(paste("intercept parameter (", beta[0], ") : "))))  +  
theme(axis.title = element_text(size = 15), axis.text.x = element_text(vjust = 3), legend.text=element_text(size=15), axis.ticks = element_blank(), legend.title=element_text(size=16), axis.text = element_text(size=10.5),   panel.grid.minor.x = element_blank(), plot.title=element_text(size=17), strip.text.x = element_text(size = 13, colour = "black", angle = 0), strip.text.y = element_text(size = 12, colour = "black", angle = 270))+ 

geom_hline(yintercept = alpha_sig_level) + labs(title=expression(paste("Probability that the ZINB model rejects the null of ", beta[X], " = ", gamma[X], " = 0")))

pdf("~/plots/zinb_pval.pdf", width=9.5, height=7)

p +  geom_text(
    size    = 4,
    data    = resultsmatsmall,
    mapping = aes(x = Inf, y = Inf, label = as.numeric(lab)),
    hjust   = 1.05,
    vjust   = 1.5, col="black"
  ) + 
  scale_y_continuous(breaks = c(0, 0.05, 0.10, 0.15), labels=c("0", "0.05", "0.10", "0.15"), sec.axis = sec_axis(trans ~ ., labels = NULL,  name = expression(" " %<-% "direction of increasing overdispersion" %<-% "\n"), t))  + theme(legend.position="bottom") +
  scale_x_continuous(trans='log10', breaks = c(50, 250, 1000, 10000), sec.axis = dup_axis(trans ~ ., labels = NULL, name = expression(" " %->% "direction of increasing zero-inflation" %->% "\n"))) 

dev.off() 
  
  
########  

resultsmatsmall[resultsmatsmall[,"n"]>250 & resultsmatsmall[,"phi"]%in%unique(resultsmatsmall[,"phi"])[2] & resultsmatsmall[,"omega"]%in%unique(resultsmatsmall[,"omega"])[1], "nb_pval"]

resultsmatsmall[resultsmatsmall[,"n"]==100 
& resultsmatsmall[,"beta_0"] == 0.5
& resultsmatsmall[,"phi"]%in%unique(resultsmatsmall[,"phi"])[5] 
& resultsmatsmall[,"omega"]%in%unique(resultsmatsmall[,"omega"])[5], "zinb_pval"]

resultsmatsmall[resultsmatsmall[,"n"]==100 
& resultsmatsmall[,"beta_0"] == 2.5
& resultsmatsmall[,"phi"]%in%unique(resultsmatsmall[,"phi"])[2] 
& resultsmatsmall[,"omega"]%in%unique(resultsmatsmall[,"omega"])[2], "zinb_pval"]

resultsmatsmall[
 resultsmatsmall[,"phi"]%in%unique(resultsmatsmall[,"phi"])[1] 
& resultsmatsmall[,"omega"]%in%unique(resultsmatsmall[,"omega"])[1], "selectionLRT"]


range(resultsmatsmall[resultsmatsmall[,"phi"]%in%unique(resultsmatsmall[,"phi"])[1] & resultsmatsmall[,"omega"]%in%unique(resultsmatsmall[,"omega"])[1], "LRT_pval"])

(resultsmatsmall[resultsmatsmall[,c("n")]%in%100 & resultsmatsmall[,c("beta_0")]%in%0.5 ,c("beta_0","n","omega", "phi","correct")])


(resultsmatsmall[resultsmatsmall[,c("n")]%in%100 & resultsmatsmall[,c("beta_0")]%in%2.5 ,c("beta_0","n","omega", "phi","correct")])

print(resultsmatsmall[resultsmatsmall[,"phi"]%in%resultsmatsmall[,"phi"][1] & resultsmatsmall[,"omega"]%in%resultsmatsmall[,"omega"][1], c("n", "beta_0", "phi", "omega", "typeAIC")][, "typeAIC"],digits=3)

range(round(resultsmatsmall[resultsmatsmall[,"phi"]%in%resultsmatsmall[,"phi"][1] & resultsmatsmall[,"omega"]%in%resultsmatsmall[,"omega"][1], c("n", "beta_0", "phi", "omega", "selectionLRT")][,"selectionLRT"],3))


range(round(resultsmatsmall[resultsmatsmall[,"phi"]%in%resultsmatsmall[,"phi"][1] & resultsmatsmall[,"omega"]%in%resultsmatsmall[,"omega"][1], c("n", "beta_0", "phi", "omega", "correct_selection")][,"correct_selection"],3))


resultsmatsmall[resultsmatsmall[,"beta_0"]%in%0.5 & resultsmatsmall[,"phi"]%in%resultsmatsmall[,"phi"][1] & resultsmatsmall[,"omega"]%in%resultsmatsmall[,"omega"][1], c("n", "beta_0", "phi", "omega", "typeAIC")][,]

#"scenario 3"
print(resultsmatsmall[resultsmatsmall[,"n"]%in%250  & resultsmatsmall[,"beta_0"]%in%0.5  & resultsmatsmall[,"phi"]%in%unique(resultsmatsmall[,"phi"])[1] & resultsmatsmall[,"omega"]%in%resultsmatsmall[,"omega"][1],c("poisson_pval",  "zip_pval", "nb_pval", "zinb_pval")], digits=2)

print(resultsmatsmall[resultsmatsmall[,"n"]%in%250  & resultsmatsmall[,"beta_0"]%in%0.5  & resultsmatsmall[,"phi"]%in%unique(resultsmatsmall[,"phi"])[1] & resultsmatsmall[,"omega"]%in%resultsmatsmall[,"omega"][1],], digits=2)


#"scenario 6"
print(resultsmatsmall[resultsmatsmall[,"n"]%in%2000  & resultsmatsmall[,"beta_0"]%in%0.5  & resultsmatsmall[,"phi"]%in%unique(resultsmatsmall[,"phi"])[1] & resultsmatsmall[,"omega"]%in%unique(resultsmatsmall[,"omega"])[1],c("type1", "typeAIC", "typeBIC", "typeAICc", "poisson_pval",  "zip_pval", "nb_pval", "zinb_pval")], digits=2)


#"scenario 46"
print(resultsmatsmall[resultsmatsmall[,"n"]%in%2000  & resultsmatsmall[,"beta_0"]%in%0.5  & resultsmatsmall[,"phi"]%in%unique(resultsmatsmall[,"phi"])[2] & resultsmatsmall[,"omega"]%in%unique(resultsmatsmall[,"omega"])[1],c("type1", "typeAIC", "typeBIC", "typeAICc", "poisson_pval",  "zip_pval", "nb_pval", "zinb_pval")], digits=2)


#"scenario 57"
print(resultsmatsmall[resultsmatsmall[,"n"]%in%50  & resultsmatsmall[,"beta_0"]%in%1.5  & resultsmatsmall[,"phi"]%in%unique(resultsmatsmall[,"phi"])[2] & resultsmatsmall[,"omega"]%in%resultsmatsmall[,"omega"][1],c("Poisson_pval",  "zip_pval", "nb_pval", "zinb_pval")], digits=2)

print(resultsmatsmall[resultsmatsmall[,"n"]%in%50  & resultsmatsmall[,"beta_0"]%in%1.5  & resultsmatsmall[,"phi"]%in%unique(resultsmatsmall[,"phi"])[2] & resultsmatsmall[,"omega"]%in%resultsmatsmall[,"omega"][1],], digits=2)



#"scenario 242"
print(resultsmatsmall[resultsmatsmall[,"n"]%in%100  & resultsmatsmall[,"beta_0"]%in%0.5  & resultsmatsmall[,"phi"]%in%unique(resultsmatsmall[,"phi"])[2] & resultsmatsmall[,"omega"]%in%unique(resultsmatsmall[,"omega"])[2],c("Poisson_pval",  "zip_pval", "nb_pval", "zinb_pval")], digits=2)


#"scenario 251"
print(resultsmatsmall[resultsmatsmall[,"n"]%in%250  & resultsmatsmall[,"beta_0"]%in%1  & resultsmatsmall[,"phi"]%in%unique(resultsmatsmall[,"phi"])[2] & resultsmatsmall[,"omega"]%in%unique(resultsmatsmall[,"omega"])[2],c("propNBBIC", "propZINBBIC", "type1", "typeAIC","typeBIC", "poisson_pval",  "zip_pval", "nb_pval", "zinb_pval")], digits=2)

print(resultsmatsmall[resultsmatsmall[,"n"]%in%250  & resultsmatsmall[,"beta_0"]%in%1  & resultsmatsmall[,"phi"]%in%unique(resultsmatsmall[,"phi"])[2] & resultsmatsmall[,"omega"]%in%unique(resultsmatsmall[,"omega"])[2],], digits=2)


print(resultsmatsmall[resultsmatsmall[,"n"]%in%250  & resultsmatsmall[,"beta_0"]%in%1  & resultsmatsmall[,"phi"]%in%unique(resultsmatsmall[,"phi"])[2] & resultsmatsmall[,"omega"]%in%unique(resultsmatsmall[,"omega"])[2],c("propZINBAICc")], digits=2)


#"scenario 402"
print(resultsmatsmall[resultsmatsmall[,"n"]%in%100  & resultsmatsmall[,"beta_0"]%in%0.5  & resultsmatsmall[,"phi"]%in%unique(resultsmatsmall[,"phi"])[1] & resultsmatsmall[,"omega"]%in%unique(resultsmatsmall[,"omega"])[3],c("Poisson_pval",  "zip_pval", "nb_pval", "zinb_pval")], digits=2)


print(resultsmatsmall[resultsmatsmall[,"n"]%in%100  & resultsmatsmall[,"beta_0"]%in%0.5  & resultsmatsmall[,"phi"]%in%unique(resultsmatsmall[,"phi"])[1] & resultsmatsmall[,"omega"]%in%unique(resultsmatsmall[,"omega"])[3],], digits=2)





maketable1<-function(iii){
	
results_iii <- readRDS(paste(paste("~/cluster_results/results_Nov2020_15000",iii, sep="_"),".rds",sep=""))	


LRT_ <- mean(na.rm = FALSE, results_iii[,"LRT_pval"]<0.05)
A_<-mean(na.rm=TRUE, results_iii[results_iii[,"LRT_pval"]>=0.05, "vuong_P_zip_pval"]<0.05)
B_<-mean(na.rm=TRUE, results_iii[results_iii[,"LRT_pval"]<0.05, "vuong_NB_zinb_pval"]<0.05)
C_<-mean(na.rm=TRUE, results_iii[results_iii[,"LRT_pval"]>=0.05 & results_iii[,"vuong_P_zip_pval"]>=0.05, "poisson_pval"]<0.05)
D_<-mean(na.rm=TRUE, results_iii[results_iii[,"LRT_pval"]>=0.05 & results_iii[,"vuong_P_zip_pval"]<0.05, "zip_pval"]<0.05)
E_<-mean(na.rm=TRUE, results_iii[results_iii[,"LRT_pval"]<0.05 & results_iii[,"vuong_NB_zinb_pval"]>=0.05, "nb_pval"]<0.05)
F_<-mean(na.rm=TRUE, results_iii[results_iii[,"LRT_pval"]<0.05 & results_iii[,"vuong_NB_zinb_pval"]<0.05, "zinb_pval"]<0.05)




A_n<-length(na.omit(results_iii[results_iii[,"LRT_pval"]>=0.05, "vuong_P_zip_pval"]))
B_n<-length(na.omit(results_iii[results_iii[,"LRT_pval"]<0.05, "vuong_NB_zinb_pval"]))
C_n<-length(na.omit(results_iii[results_iii[,"LRT_pval"]>=0.05 & results_iii[,"vuong_P_zip_pval"]>=0.05, "poisson_pval"]))
D_n<-length(na.omit(results_iii[results_iii[,"LRT_pval"]>=0.05 & results_iii[,"vuong_P_zip_pval"]<0.05, "zip_pval"]))
E_n<-length(na.omit(results_iii[results_iii[,"LRT_pval"]<0.05 & results_iii[,"vuong_NB_zinb_pval"]>=0.05, "nb_pval"]))
F_n<-length(na.omit(results_iii[results_iii[,"LRT_pval"]<0.05 & results_iii[,"vuong_NB_zinb_pval"]<0.05, "zinb_pval"]))
c(A_n, B_n, C_n, D_n, E_n, F_n)



C_AIC<-mean(na.rm=TRUE, results_iii[results_iii[,"choiceAIC"]==1, "poisson_pval"]<0.05)
D_AIC<-mean(na.rm=TRUE, results_iii[results_iii[,"choiceAIC"]==3, "zip_pval"]<0.05)
E_AIC<-mean(na.rm=TRUE, results_iii[results_iii[,"choiceAIC"]==2, "nb_pval"]<0.05)
F_AIC<-mean(na.rm=TRUE, results_iii[results_iii[,"choiceAIC"]==4, "zinb_pval"]<0.05)

C_BIC<-mean(na.rm=TRUE, results_iii[results_iii[,"choiceBIC"]==1, "poisson_pval"]<0.05)
D_BIC<-mean(na.rm=TRUE, results_iii[results_iii[,"choiceBIC"]==3, "zip_pval"]<0.05)
E_BIC<-mean(na.rm=TRUE, results_iii[results_iii[,"choiceBIC"]==2, "nb_pval"]<0.05)
F_BIC<-mean(na.rm=TRUE, results_iii[results_iii[,"choiceBIC"]==4, "zinb_pval"]<0.05)

# Figure 1 numbers:
print(cbind(c(
100*round(mean(na.rm=TRUE, results_iii[,"LRT_pval"]>=0.05),4),
100*round(mean(na.rm=TRUE, results_iii[,"LRT_pval"]<0.05),4),

100*round(mean(na.rm=TRUE, results_iii[,"LRT_pval"]>=0.05 & results_iii[,"vuong_P_zip_pval"]>=0.05),4),
100*round(mean(na.rm=TRUE, results_iii[,"LRT_pval"]>=0.05 & results_iii[,"vuong_P_zip_pval"]<0.05),4),
100*round(mean(na.rm=TRUE, results_iii[,"LRT_pval"]<0.05 & results_iii[,"vuong_NB_zinb_pval"]>=0.05),4),
100*round(mean(na.rm=TRUE, results_iii[,"LRT_pval"]<0.05 & results_iii[,"vuong_NB_zinb_pval"]<0.05),4),



round(sum(na.rm=TRUE, results_iii[results_iii[,"choice"]==1, "poisson_pval"]>=0.05)/150,2),
round(sum(na.rm=TRUE, results_iii[results_iii[,"choice"]==1, "poisson_pval"]<0.05)/150,2),
round(sum(na.rm=TRUE, results_iii[results_iii[,"choice"]==3, "zip_pval"]>=0.05)/150,2),
round(sum(na.rm=TRUE, results_iii[results_iii[,"choice"]==3, "zip_pval"]<0.05)/150,2),
round(sum(na.rm=TRUE, results_iii[results_iii[,"choice"]==2, "nb_pval"]>=0.05)/150,2),
round(sum(na.rm=TRUE, results_iii[results_iii[,"choice"]==2, "nb_pval"]<0.05)/150,2),
round(sum(na.rm=TRUE, results_iii[results_iii[,"choice"]==4, "zinb_pval"]>=0.05)/150,2),
round(sum(na.rm=TRUE, results_iii[results_iii[,"choice"]==4, "zinb_pval"]<0.05)/150,2),


round(sum(na.rm=TRUE, results_iii[results_iii[,"choiceAIC"]==1, "poisson_pval"]>=0.05)/150,2),
round(sum(na.rm=TRUE, results_iii[results_iii[,"choiceAIC"]==1, "poisson_pval"]<0.05)/150,2),
round(sum(na.rm=TRUE, results_iii[results_iii[,"choiceAIC"]==3, "zip_pval"]>=0.05)/150,2),
round(sum(na.rm=TRUE, results_iii[results_iii[,"choiceAIC"]==3, "zip_pval"]<0.05)/150,2),
round(sum(na.rm=TRUE, results_iii[results_iii[,"choiceAIC"]==2, "nb_pval"]>=0.05)/150,2),
round(sum(na.rm=TRUE, results_iii[results_iii[,"choiceAIC"]==2, "nb_pval"]<0.05)/150,2),
round(sum(na.rm=TRUE, results_iii[results_iii[,"choiceAIC"]==4, "zinb_pval"]>=0.05)/150,2),
round(sum(na.rm=TRUE, results_iii[results_iii[,"choiceAIC"]==4, "zinb_pval"]<0.05)/150,2),


round(sum(na.rm=TRUE, results_iii[results_iii[,"choiceBIC"]==1, "poisson_pval"]>=0.05)/150,2),
round(sum(na.rm=TRUE, results_iii[results_iii[,"choiceBIC"]==1, "poisson_pval"]<0.05)/150,2),
round(sum(na.rm=TRUE, results_iii[results_iii[,"choiceBIC"]==3, "zip_pval"]>=0.05)/150,2),
round(sum(na.rm=TRUE, results_iii[results_iii[,"choiceBIC"]==3, "zip_pval"]<0.05)/150,2),
round(sum(na.rm=TRUE, results_iii[results_iii[,"choiceBIC"]==2, "nb_pval"]>=0.05)/150,2),
round(sum(na.rm=TRUE, results_iii[results_iii[,"choiceBIC"]==2, "nb_pval"]<0.05)/150,2),
round(sum(na.rm=TRUE, results_iii[results_iii[,"choiceBIC"]==4, "zinb_pval"]>=0.05)/150,2),
round(sum(na.rm=TRUE, results_iii[results_iii[,"choiceBIC"]==4, "zinb_pval"]<0.05)/150,2))))




C_all<-mean(na.rm=TRUE, results_iii[, "poisson_pval"]<0.05)
D_all<-mean(na.rm=TRUE, results_iii[, "zip_pval"]<0.05)
E_all<-mean(na.rm=TRUE, results_iii[, "nb_pval"]<0.05)
F_all<-mean(na.rm=TRUE, results_iii[, "zinb_pval"]<0.05)




cbind(c(C_BIC, D_BIC, E_BIC, F_BIC), c(C_AIC, D_AIC, E_AIC, F_AIC),c(C_, D_, E_, F_),c(C_all, D_all, E_all, F_all))

sevenstepselect<-as.numeric(resultsmatsmall[iii,c("propPoisson",      "propZIP",    "propNB" ,"propZINB")])

AICselect<-as.numeric(c(resultsmatsmall[iii,c("propPoissonAIC",      "propZIPAIC",    "propNBAIC" ,"propZINBAIC")]))

BICselect<-as.numeric(c(resultsmatsmall[iii,c("propPoissonBIC",      "propZIPBIC",    "propNBBIC" ,"propZINBBIC")]))


print(rbind(
c(C_all, D_all, E_all, F_all),
c(C_, D_, E_, F_),
sevenstepselect,
c(C_AIC, D_AIC, E_AIC, F_AIC),
AICselect,
c(C_BIC, D_BIC, E_BIC, F_BIC),
BICselect))

table1<-print(xtable(rbind(
c(C_all, D_all, E_all, F_all),
c(C_, D_, E_, F_),
sevenstepselect,
c(C_AIC, D_AIC, E_AIC, F_AIC),
AICselect,
c(C_BIC, D_BIC, E_BIC, F_BIC),
BICselect)), include.rownames=FALSE)

return(table1)

}
library(xtable)
maketable1(3)
maketable1(6)
maketable1(46)
maketable1(57)
maketable1(242)
maketable1(402)
maketable1(243)
maketable1(251)


mean(resultsmatsmall[, "correct_selection"])
mean(resultsmatsmall[, "correct_selectionAIC"])
mean(resultsmatsmall[, "correct_selectionAICc"])
mean(resultsmatsmall[, "correct_selectionBIC"])

mean(resultsmatsmall[, "type1"])
mean(resultsmatsmall[, "typeAIC"])
mean(resultsmatsmall[, "typeAICc"])
mean(resultsmatsmall[, "typeBIC"])

