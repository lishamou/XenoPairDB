args <- commandArgs(trailingOnly=TRUE)
B <- if(length(args)) normalizePath(args[1]) else normalizePath('.')
O <- file.path(B,'figures');dir.create(O,showWarnings=FALSE)
dir.create(file.path(B,'logs'),showWarnings=FALSE)
blue <- '#41789C'; orange <- '#C87942'; gray <- '#B9C2CB'
read <- function(n) read.delim(file.path(O,paste0('figure',n,'_source.tsv')),check.names=FALSE,stringsAsFactors=FALSE)
basepar <- function() par(family='sans',cex=.8,las=1,bty='l',mgp=c(2.1,.5,0),tcl=-.25)
plotfig <- function(n) {
 d <- read(n)
 if(n==1){
  par(mfrow=c(1,2),mar=c(4,8,3,1));basepar()
  evidence<-factor(d$evidence,levels=c('direct','disputed','control'));counts<-table(evidence,d$system)
  barplot(counts,horiz=TRUE,col=c(blue,orange,gray),border=NA,xlab='Molecular axes',main='a  Catalogue composition',cex.names=.8)
  legend('bottomright',legend=c('Direct','Disputed','Control'),fill=c(blue,orange,gray),bty='n',cex=.7)
  par(mar=c(4,4,3,1));z<-table(factor(d$sequence_available,levels=c('True','False')),factor(ifelse(d$evidence=='control','Control','Non-control'),levels=c('Control','Non-control')))
  xx<-barplot(z,col=c(blue,gray),border=NA,ylim=c(0,72),ylab='Molecular axes',main='b  Sequence availability');text(rep(xx,each=2),c(z[1,1]/2,z[1,1]+z[2,1]/2,z[1,2]/2,z[1,2]+z[2,2]/2),labels=c(z))
  legend('topleft',legend=c('Available','Unavailable'),fill=c(blue,gray),bty='n',cex=.7)
 }
 if(n==2){
  par(mfrow=c(1,1),mar=c(4,9,3,1));basepar()
  z<-rbind(d$A,d$R);colnames(z)<-d$donor_short
  barplot(z,horiz=TRUE,col=c(blue,gray),border=NA,xlab='Mixed-divergence scenario score',cex.names=.8,ylim=c(0,9),xlim=c(0,1200),xaxt='n',main='Coverage within the sequence-evaluable subset')
  axis(1,at=seq(0,1200,200))
  legend('top',legend=c('Genetically addressed','Residual'),fill=c(blue,gray),bty='n',horiz=TRUE,cex=.75)
 }
 if(n==3){
  par(mfrow=c(1,2),mar=c(4,4,3,1));basepar();a<-d[d$panel=='a',];ref<-as.character(a$historical_positive)=='True'
  plot(ecdf(a$combined[!ref]),col=gray,lwd=2,verticals=TRUE,do.points=FALSE,xlim=c(0,1),xlab='Recalculated reference-set score',ylab='Cumulative fraction',main='a  Reference-set distributions')
  plot(ecdf(a$combined[ref]),add=TRUE,col=blue,lwd=2,verticals=TRUE,do.points=FALSE)
  legend(x=.02,y=.93,legend=c('Untested background (n=4,114)','Literature reference (n=8)'),col=c(gray,blue),lwd=2,bty='n',cex=.62)
  c<-d[d$panel=='b',];p<-as.character(c$primary)=='True';plot(c$antigenicity_v2,c$log2fc,type='n',xlab='Catalogue antigen-priority score',ylab='Median peptide log2 fold change',main='b  External peptide comparison',ylim=range(c$log2fc)+c(-.3,.6))
  points(c$antigenicity_v2[p],c$log2fc[p],pch=16,col=blue);points(c$antigenicity_v2[!p],c$log2fc[!p],pch=1,col=orange,cex=1.2);text(c$antigenicity_v2[!p],c$log2fc[!p],labels='RPN1',pos=3,cex=.75)
  legend('bottomleft',legend=c('Primary cohort (n=6)','Post hoc RPN1 addition'),col=c(blue,orange),pch=c(16,1),bty='n',cex=.65)
 }
 if(n==4){
  par(mfrow=c(1,2),mar=c(5,4,3,1));basepar();genes<-c('ICAM1','ULBP1','CD86','CD274','THBD');raw<-d[d$mode=='raw_all_CA',];raw<-raw[match(genes,raw$gene),];filtered<-d[d$mode=='filtered_CA',];filtered<-filtered[match(genes,filtered$gene),]
  x<-rbind(raw$rmsd_fit,filtered$rmsd_fit);colnames(x)<-genes;barplot(x,beside=TRUE,col=c(blue,orange),border=NA,ylim=c(0,max(x)*1.15),ylab='C-alpha RMSD (Å)',main='a  Alignment sensitivity',las=2);legend('topright',legend=c('Raw fit','Filtered fit'),fill=c(blue,orange),bty='n',cex=.7)
  x<-rbind(raw$residue_pairs,filtered$residue_pairs);colnames(x)<-genes;barplot(x,beside=TRUE,col=c(blue,orange),border=NA,ylim=c(0,max(x)*1.2),ylab='Residue-pair contacts at 5 Å',main='b  Contact summaries',las=2);legend('topright',legend=c('Raw fit','Filtered fit'),fill=c(blue,orange),bty='n',cex=.7)
 }
}
for(n in 1:4){
 png(file.path(O,paste0('figure',n,'.png')),width=7.1,height=3.3,units='in',res=300,type='quartz');plotfig(n);dev.off()
 pdf(file.path(O,paste0('figure',n,'.pdf')),width=7.1,height=3.3,useDingbats=FALSE);plotfig(n);dev.off()
 tiff(file.path(O,paste0('figure',n,'.tiff')),width=7.1,height=3.3,units='in',res=600,type='quartz');plotfig(n);dev.off()
}
writeLines(capture.output(sessionInfo()),file.path(B,'logs','R_SESSION_INFO.txt'))
cat('Four R figures exported from frozen source tables\n')
