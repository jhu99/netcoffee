# R code: format function changing various dataformat to standard format
format <- function(filelist)#filename1 inputfile filename2 outputfile
{
  for(filename in filelist)
  {
    originData <- as.matrix(read.csv(filename, sep = "\t", quote="", header=FALSE, comment="!"));
    simpleData <- originData;
    nrow <- dim(simpleData)[1];
    for(i in 1:nrow)
    {
      protein1 <- simpleData[i,1];
      protein2 <- simpleData[i,2];
      part1 <- strsplit(protein1,split="|",fixed=TRUE)[[1]][2];
      part2 <- strsplit(protein2,split="|",fixed=TRUE)[[1]][2];
      simpleData[i,1] <- part1;
      simpleData[i,2] <- part2;
    }
    filename2 <- paste(filename,'data',sep=".");
    write.table(simpleData,file=filename2,sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE);
  }
}
#format test cf file
format_testfile <- function(filelist)#filename1 inputfile filename2 outputfile
{
  for(filename in filelist)
  {
    originData <- as.matrix(read.csv(filename, sep = " ", quote="", header=TRUE, row.names=NULL,comment="!"));
    simpleData <- originData;
    nrow <- dim(simpleData)[1];
    s <- rep(10^-50,nrow);
    simpleData <- cbind(simpleData,s);
    filename2 <- paste(filename,'data',sep=".");
    write.table(simpleData,file=filename,sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE);
  }
}

format_dip <- function(filelist)#filename1 inputfile filename2 outputfile
{
  for(filename in filelist)
  {
    originData <- as.matrix(read.csv(filename, sep = "\t", quote="", header=TRUE, row.names=NULL,comment="!"));
    simpleData <- originData;
    nrow <- dim(simpleData)[1];
    for(i in 1:nrow)
    {
      protein1 <- simpleData[i,1];
      protein2 <- simpleData[i,2];
      part1 <- strsplit(protein1,split="|",fixed=TRUE)[[1]][1];
      part2 <- strsplit(protein2,split="|",fixed=TRUE)[[1]][1];
      simpleData[i,1] <- part1;
      simpleData[i,2] <- part2;
    }
    filename2 <- paste(filename,'data',sep=".");
    write.table(simpleData[,c(1,2,10,11)],file=filename2,sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE);
  }
}
format_dip_blast <- function(filelist)#filename1 inputfile filename2 outputfile
{
  for(filename in filelist)
  {
    originData <- as.matrix(read.csv(filename, sep = "\t", quote="", header=TRUE, row.names=NULL,comment="!"));
    simpleData <- originData;
    nrow <- dim(simpleData)[1];
    for(i in 1:nrow)
    {
      protein1 <- simpleData[i,1];
      protein2 <- simpleData[i,2];
      part1 <- strsplit(protein1,split="|",fixed=TRUE)[[1]][3];
      part2 <- strsplit(protein2,split="|",fixed=TRUE)[[1]][3];
      simpleData[i,1] <- part1;
      simpleData[i,2] <- part2;
    }
    filename2 <- paste(filename,'data',sep=".");
    write.table(simpleData,file=filename2,sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE);
  }
}
format_dip_int2 <- function(filelist)
{
  for(filename in filelist)
  {
    originData <- as.matrix(read.csv(filename, sep = "\t", quote="", header=TRUE, row.names=NULL,comment="!"));
    simpleData <- originData;
    nrow <- dim(simpleData)[1];
    ones <- rep(0.9,nrow);
    simpleData <- cbind(simpleData,ones);
    filename2 <- paste(filename,'data',sep=".");
    write.table(simpleData,file=filename2,sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE);
  }
}
formatSelection <- function(filename,filename2,cols)#filename1 inputfile filename2 outputfile
{
  originData <- as.matrix(read.csv(filename, sep = "\t", quote="", header=FALSE, comment="!"));
  simpleData <- originData[,cols];
  write.table(simpleData,file=filename2,sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE);
}

formatNet2Tab <- function(filelist)
{
	for (filename in filelist)
	{
		filename2=paste(filename,".tab",sep="");
		originData <- read.csv(filename, sep = "\t", quote="", header=FALSE, comment="!", skip=2);
		simpleData <- originData[,c(1,2)];
    simpleData=rbind(c("INTERACTOR_A","INTERACTOR_B"),simpleData);
		write.table(simpleData,file=filename2,sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE);
	}
}
formatModel <- function(fileslist)
{
  browser();
  for (filename in fileslist)
  {
    originData <- as.matrix(read.csv(filename, sep = "\t", quote="", header=FALSE));
    rowNum <- dim(originData)[1];
    nonreduncancy <- hash();
    requestData <- c();
    for ( row in 1:rowNum)
    {
      protein1 <- originData[row,1];
      protein2 <- originData[row,2];
      if (protein1 > protein2)
      {
        temp = protein1;
        protein1=protein2;
        protein2=temp;
      }
      key <- paste(protein1,protein2,sep="");
      if(protein1 == protein2)
      {
        ;
      }
      else
      {
        if (is.null(nonreduncancy[[key]]))
        {
          .set(nonreduncancy,key,1);
          entries <- originData[row,];
          requestData <- rbind(requestData,entries);
        }else
        {
          nonreduncancy[[key]]=nonreduncancy[[key]]+1;
        }
      }
    }
    fileparts <- strsplit(filename,split="NA",fixed=TRUE);
    outfilename <- paste(fileparts[[1]][1],fileparts[[1]][2],fileparts[[1]][3], sep="");
    write.table(requestData,file=outfilename,sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE);
  }
}
formatModelUp <- function(fileslist)
{
  for (filename in fileslist)
  {
    originData <- as.matrix(read.csv(filename, sep = "\t", quote="", header=FALSE));
    rowNum <- dim(originData)[1];
    requestData <- c();
    for ( row in 1:rowNum)
    {
      protein1 <- originData[row,1];
      protein2 <- originData[row,2];
      if(row ==1 )
      {
        entries <- originData[row,];
        idlabel1 <- entries[1];
        idlabel2 <- entries[2];
        idac1 <- strsplit(idlabel1,split="|",fixed=TRUE)[[1]][2];
        idac2 <- strsplit(idlabel2,split="|",fixed=TRUE)[[1]][2];
        entries[1] <- idac1; entries[2] <- idac2;
        requestData <- rbind(requestData,entries);
      }
      else if (protein1 != originData[row-1,1] ||
        protein2 != originData[row-1,2])
      {
        entries <- originData[row,];
        idlabel1 <- entries[1];
        idlabel2 <- entries[2];
        idac1 <- strsplit(idlabel1,split="|",fixed=TRUE)[[1]][2];
        idac2 <- strsplit(idlabel2,split="|",fixed=TRUE)[[1]][2];
        entries[1] <- idac1; entries[2] <- idac2;
        requestData <- rbind(requestData,entries);
      }else
      {
      }
    }
    fileparts <- strsplit(filename,split=".homology",fixed=TRUE);
    outfilename <- paste(fileparts[[1]][1],fileparts[[1]][2],fileparts[[1]][3], sep="");
    write.table(requestData,file=outfilename,sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE);
  }
}
addEntries <- function(object,value,keys)
{
  if(keys[1]=="")
  {
  }
  else{
    keylist <- strsplit(keys,"?;? ");#split string to vector
    for (i in keylist)
    {
      .set(object,keys=i,values=value);
    }
  }
}
empiricalDistribution <- function(filelist)
{
  i = 1;
  density <- c();
  browser();
  for (filename in filelist)
  {
    originData <- as.matrix(read.csv(filename, sep = "\t", quote="", header=FALSE));
    x <- as.numeric(originData[,11]);   #e-value 
    d <- density(x, bw= "nrd0", kernel= "epanechnikov");
    density <- c(density,d);
    outfile <- paste(filename,"distribution", "data", sep=".");
    #save(d,file=outfile);
    #write.table(d,file=outfile,sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE);
    hist(x);
     pdf(file=paste(filename,"histgram", "pdf",sep="."))
     hist(x,prob=TRUE, breaks=400,col = "grey", plot=TRUE,xlab="e-value",main="");
     lines(d,col = "red");
     dev.off();
     pdf(paste(filename,"distribution","pdf",sep="."))
     plot(d,col ="red",xlab="e-value", main="");
     dev.off();
  } 
  save(density,file="distribution.Rdata");
}
loadIdMaps <- function(idfiles)
{
  library("hash");
  ref2Uni <- hash();#refSeq2UniprotKB-AC
  entrez2Uni <- hash();#entrez2UniprotKB-AC
  gi2Uni <- hash();#gi2UniprotKB-AC
  #worm2Uni <- hash();#worm2UniprotKB-AC
  uni2Goterms <- c();#UniprotKB-AC2goterms
  browser();
  for (filename in idfiles)
  {
    idData <- as.matrix(read.csv(filename, sep = "\t", quote="", header=FALSE));
    rowNum <- dim(idData)[1];
    for (row in 1:rowNum)
    {
#       entrez <- idData[row,3];
#       addEntries(entrez2Uni,idData[row,1],entrez);
#       refSeq <- idData[row,4];
#       addEntries(ref2Uni,idData[row,1],refSeq);
#       gi <- idData[row,5];
#       addEntries(gi2Uni,idData[row,1],gi);
      goterm <- idData[row,7];
      x<-strsplit(goterm,"?;? ");
      uni2Goterms <- rbind(uni2Goterms,c(idData[row,1],x[[1]]));
#      uni2Goterms[idData[row,1]] <- strsplit(goterm,"?;? ");
    }
  }
  browser();
  save(goterm,gi2Uni,ref2Uni,entrez2Uni,file="./UniprotKBIdMap.RData")
}
uniprotId <- function(names)
{
  for (name in names)
  {
    nameData <- strsplit(name,"[:]");
    type <- strsplit(nameData[[1]][1],"[ /]");
    id <- nameData[[1]][2];
    UniprotKBId = switch(type[[1]][1],
                         genbank = gi2Uni[[id]],
                         refseq = ref2Uni[[id]],
                         entrezgene = entrez2Uni[[id]],
                         uniprotkb = id,
                         NULL
    )
    if ( is.null(UniprotKBId) == FALSE ){
      return(UniprotKBId);
    }
  }
  return(NULL);
}
#Format raw interaction to a simple framework.
formatInteraction <- function(files)
{
  for (filename in files)
  {
    originData <- as.matrix(read.csv(filename, sep = "\t", quote="", header=FALSE, comment="#"));
    simpleData <- originData[,c(1,2,15)];
    nrow <- dim(simpleData)[1];
    for(i in 1:nrow)
    {
      protein1 <- simpleData[i,1];
      protein2 <- simpleData[i,2];
      score <- simpleData[i,3];
      part1 <- strsplit(protein1,split=":",fixed=TRUE)[[1]][2];
      part2 <- strsplit(protein2,split=":",fixed=TRUE)[[1]][2];
      part3 <- strsplit(score,split="miscore:",fixed=TRUE)[[1]][2];
      simpleData[i,1] <- part1;
      simpleData[i,2] <- part2;
      simpleData[i,3] <- part3;
    }
    filename2 <- paste(filename,'si',sep=".");#simple frame work
    write.table(simpleData,file=filename2,sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE);
  }
}
networkUni <- function(species)
{
  library("hash");
  networks <- c();
  browser();
  for (sp in species)
  {
    files <- list.files(path="./",pattern=sp);
    network <- hash();
    outputfile <- paste(sp,".int", sep ="");
    for (filename in files)
    {
      simpleData <- as.matrix(read.csv(filename, sep = "\t", quote="", header=TRUE));
      rowNum <- dim(simpleData)[1];
      columnNum <- dim(simpleData)[2];
      for (row in 1:rowNum)
      {
        proteinA <- simpleData[row,1];
        namesA <- strsplit(proteinA,"|", fixed = TRUE);
        nameA <- uniprotId(namesA);
        proteinB <- simpleData[row,2];
        namesB <- strsplit(proteinB,"|", fixed=TRUE);
        nameB <- uniprotId(namesB);
        if( !is.null(nameA) && !is.null(nameB))
        {
          if(nameA > nameB)
          {
            temp = nameA;
            nameA = nameB;
            nameB =temp;
          }
          key <- paste(nameA,nameB,sep=";");
          if (is.null(network[[key]]))
          {
            .set(network,key,1);
          }else
          {
            values(network, keys=key) <-  values(network, keys=key)+1;
          }
        }
      }
    }
    #networks <- c(networks,network);
  }
#  save(networks,file="./networks.RData");
}

connectMySQL <- function()
{
  library("RMySQL");
  con <- dbConnect(MySQL(),user="root",password="dlh4721rt",dbname="homology_list");
  return(con);
}

calculateScore <- function(con)
{
  score <- dbGetQuery(con,"select bitscore1, bitscore2, bitscore3 from table_123");
  row <- dim(score)[1];
  column <- dim(score)[2];
  nodeScores <- c();
  for (i in 1:row)
  {
    if ( i%%10000 ==0)
    {
      print(i);
    }
    nodeScore = 0;
    for (j in 1:column)
    {
      x <- score[i,j];
      currScore = dn(x,j);
      nodeScore = nodeScore + currScore;
    }
    nodeScores <- rbind(nodeScores,nodeScore);
  }
  nodeScores
  #dbWriteTable(con, "table_test", value=nodeScores, overwrite=FALSE, append = TRUE);
}

dn <- function(x,j)
{
  index1 <- ceiling((x-density[j*14-13][[1]][1])/(density[j*14-13][[1]][2]-density[j*14-13][[1]][1]));
  prob1 <- density[j*14-12][[1]][index1];
  if( is.na(prob1) )
  {
    prob1 = 1.0e-100;
  }
  if ( prob1 < 1.0e-100 )
  {
    prob1 = 1.0e-100;
  }
  index2 <- ceiling((x-density[j*14-6][[1]][1])/(density[j*14-6][[1]][2]-density[j*14-6][[1]][1]));
  prob2 <- density[j*14-5][[1]][index2];
  if( is.na(prob2) )
  {
    prob2 = 1.0e-100;
  }
  if ( prob2 < 1.0e-100 )
  {
    prob2 = 1.0e-100;
  }
  score = log10(prob1/prob2);
  score
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%extract uniprotId from sequence id in "ur" format files%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
getHomologylist <- function(filelist)
{
  for(filename in filelist)
  {
    originData <- as.matrix(read.csv(filename, sep = "\t", quote="", header=FALSE));
    rowNum <- dim(originData)[1];
    #browser();
    outfilename <- paste(filename,"data",sep=".");
    for(i in 1:rowNum)
    {
      protein1 <- originData[i,1];
      protein2 <- originData[i,2];
      originData[i,1] <- strsplit(protein1,"|", fixed=TRUE)[[1]][2];
      originData[i,2] <- strsplit(protein2,"|", fixed=TRUE)[[1]][2];
    }
    write.table(originData,file=outfilename,sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE);
  }
}

drawdistribution<-function()
{
  xoriginData1 <- as.matrix(read.csv("./u_opmodel/ce_dm.ur", sep = "\t", quote="", header=TRUE, comment="!"));
  xoriginData2 <- as.matrix(read.csv("./u_opmodel/ce_hs.ur", sep = "\t", quote="", header=TRUE, comment="!"));
  xoriginData3 <- as.matrix(read.csv("./u_opmodel/dm_hs.ur", sep = "\t", quote="", header=TRUE, comment="!"));
  yoriginData1 <- as.matrix(read.csv("./u_nullmodel/data1/ce_dm_null.ur", sep = "\t", quote="", header=TRUE, comment="!"));
  yoriginData2 <- as.matrix(read.csv("./u_nullmodel/data1/ce_hs_null.ur", sep = "\t", quote="", header=TRUE, comment="!"));
  yoriginData3 <- as.matrix(read.csv("./u_nullmodel/data1/dm_hs_null.ur", sep = "\t", quote="", header=TRUE, comment="!"));
  x <- xoriginData1[,3];
  x <- c(x,xoriginData2[,3]);
  x <- c(x,xoriginData3[,3]);
  x <- as.numeric(x);
  y <- yoriginData1;
  y <- c(y,yoriginData2);
  y <- c(y,yoriginData3);
  load("./model_distri.RData");
  browser();
  plot(ne,np,type="o",col="red",log="xy",xlab="E-Value",ylab="Probability",main="Distribution of H&N model");
  lines(he,hp,type="o",col="blue");
  legend("bottomleft",legend=c("N Model","H Model"),lwd=2.5,text.col=c("red","blue"),col=c("red","blue"),pch=c("o","o"));
}

drawblastdistribitscore <- function(filelist)
{
  allscore <- c();
  i=0;
  x <- seq(from=0,by=10,length=95);
  
  browser();
  for (filename in filelist)
  {
    i=i+1;
    score <- as.matrix(read.csv(filename, sep = " ", quote="", header=FALSE, comment="#"));
    logratio <- log10(score[,2]/score[,1]);
    filename2 = paste(filename,"pdf", sep="-bit.");
    filename3 = paste(filename,"logratio.pdf", sep="-bit-");
    pdf(file=filename2)
    par(mar=c(5,5,4,2)+0.3)
    plot(x,score[,1],type="o",,pch=1,col="red",log="y",xlab="Bitscore",ylab="Probability",cex.lab=1.5,cex.axis=1.2,cex.main=2.0);
    lines(x,score[,2],type="o",col="blue",pch=4);
    legend("topright",legend=c("N Model","H Model"),lwd=2.5,text.col=c("red","blue"),col=c("red","blue"),pch=c(1,4));
    dev.off();
    pdf(file=filename3)
    par(mar=c(5,5,4,2)+0.3)
    plot(x[1:88],logratio[1:88],type="o",pch=1,col="red",log="x",xlab="E-Value",ylab="Score",cex.lab=1.5,cex.axis=1.2,cex.main=2.0);
    dev.off();
    if(i==1)
    {
      allscore = score;
    }
    else 
    {allscore <- score + allscore;}
    print(i);
  }
  browser();
  
  allscore = allscore/i;
  logratio <- log10(allscore[,2]/allscore[,1]);
  write.table(logratio,file="./dataset/bldata/blastdatae-7/score_composit.bmodel",sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE);
  
  pdf(file="./dataset/bldata/blastdatae-7/composite-bit.pdf")
  par(mar=c(5,5,4,2)+0.3)
  plot(x,allscore[,1],type="o",pch=1,log="y",xlab="bitscore",ylab="Probability",cex.lab=2.0,cex.axis=1.5);
  lines(x,allscore[,2],type="o",pch=4);
  legend("topright",legend=c("N Model","H Model"),lwd=2.5,pch=c(1,4));
  dev.off();
  
  pdf(file="./dataset/bldata/blastdatae-7/logratio-bit.pdf")
  par(mar=c(5,5,4,2)+0.3)
  plot(x[1:86],logratio[1:86],type="o",pch=1,col="red",log="x",xlab="bitscore",ylab="Score",cex.lab=2.0,cex.axis=1.5,cex.main=2.0);
  dev.off();
  
#   pdf(file="logratio_regression.pdf")
#   par(mar=c(5,5,4,2)+0.1)
#   beta1 =  -1.176e-07;
#   beta2 = -5.576e-05;
#   beta3 = -8.437e-03;
#   beta4 = 2.378e+00 
#   ylogratio <- beta1*(log10(x))^3+beta2*(log10(x))^2+beta3*log10(x)+beta4;
#   plot(x,ylogratio,type="o",pch=1,col="red",log="x",xlab="e-Value",ylab="Score",cex.lab=2.0,cex.axis=1.5,cex.main=2.0);
#   dev.off();
#   write.table(ylogratio,file="./regression_score_composit_4s.model",sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE);
}

drawblastdistri <- function(filelist)
{
  allscore <- c();
  i=0;
  x <- 10^c(seq(from=-180,by=2,length=95));
  browser();
  for (filename in filelist)
  {
    i=i+1;
    score <- as.matrix(read.csv(filename, sep = " ", quote="", header=FALSE, comment="#"));
    logratio <- log10(score[,2]/score[,1]);
    filename2 = paste(filename,"pdf", sep=".");
    filename3 = paste(filename,"logratio.pdf", sep="-");
    pdf(file=filename2)
    par(mar=c(5,5,4,2)+0.3)
    plot(x,score[,1],type="o",,pch=1,col="red",log="xy",xlab="E-Value",ylab="Probability",cex.lab=1.5,cex.axis=1.2,cex.main=2.0);
    lines(x,score[,2],type="o",col="blue",pch=4);
    legend("topleft",legend=c("N Model","H Model"),lwd=2.5,text.col=c("red","blue"),col=c("red","blue"),pch=c(1,4));
    dev.off();
    pdf(file=filename3)
    par(mar=c(5,5,4,2)+0.3)
    plot(x[1:88],logratio[1:88],type="o",pch=1,col="red",log="x",xlab="E-Value",ylab="Score",cex.lab=1.5,cex.axis=1.2,cex.main=2.0);
    dev.off();
    if(i==1)
    {
      allscore = score;
    }
    else 
    {allscore <- score + allscore;}
    print(i);
  }
  browser();
  
  allscore = allscore/i;
  logratio <- log10(allscore[,2]/allscore[,1]);
  write.table(logratio,file="./dataset/bldata/blastdatae-7/score_composit.model",sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE);
  
  pdf(file="./dataset/bldata/blastdatae-7/composite.pdf")
  par(mar=c(5,5,4,2)+0.3)
  plot(x,allscore[,1],type="o",pch=1,log="xy",xlab="E-Value",ylab="Probability",cex.lab=1.5,cex.axis=1.2);
  lines(x,allscore[,2],type="o",pch=4);
  legend("topleft",legend=c("N Model","H Model"),lwd=2.5,pch=c(1,4));
  dev.off();
  
  pdf(file="./dataset/bldata/blastdatae-7/logratio.pdf")
  par(mar=c(5,5,4,2)+0.3)
  plot(x[1:88],logratio[1:88],type="o",pch=1,col="red",log="x",xlab="E-Value",ylab="Score",cex.lab=1.5,cex.axis=1.2,cex.main=2.0);
  dev.off();
  
#   pdf(file="logratio_regression.pdf")
#   par(mar=c(5,5,4,2)+0.1)
#   beta1 =  -1.176e-07;
#   beta2 = -5.576e-05;
#   beta3 = -8.437e-03;
#   beta4 = 2.378e+00 
#   ylogratio <- beta1*(log10(x))^3+beta2*(log10(x))^2+beta3*log10(x)+beta4;
#   plot(x,ylogratio,type="o",pch=1,col="red",log="x",xlab="e-Value",ylab="Score",cex.lab=2.0,cex.axis=1.5,cex.main=2.0);
#   dev.off();
#   write.table(ylogratio,file="./regression_score_composit_4s.model",sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE);
}

convertGOA2fsst <- function(filelist)
{
  for(filename in filelist)
  {
    outfilename=paste(filename,"fsst",sep=".");
    #outfilename=paste(filename,"evals",sep=".");
    originData <- as.matrix(read.csv(filename, sep = "\t", quote="", header=FALSE, comment="!"));
    DataView <- originData[,c(2,5,7,9)];
    #DataView <- originData[,c(1,2,12)]
#    numrow <- dim(DataView)[1];
#    newData <- DataView;
#    j=0;
#    #browser();
#    for(i in 1:numrow)
#    {
#       if(DataView[i,3]=="IEA" || DataView[i,3]=="ISS")
#		next;
#	   j=j+1;
#       newData[j,]=DataView[i,];
#    }
#    write.table(newData[1:j,],file=outfilename,sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE);
    write.table(DataView,file=outfilename,sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE);
  }
}
clearzero <-function(mat)
{
  len <- dim(mat,1);
  newmat <- c();
  for(i in 1:len)
  {
    x <- mat[i,];
  }
}
extractMulFunsim <- function(files)
{
  mycol <- c("grey","gold");
  boxcol <- c("blue","red");
  originData1 <- as.matrix(read.csv(files[1], sep = "\t", quote="", header=TRUE, comment="#"));
  #originData2 <- as.matrix(read.csv(files[2], sep = "\t", quote="", header=TRUE, comment="#"));
  originData3 <- as.matrix(read.csv(files[3], sep = "\t", quote="", header=TRUE, comment="#"));
  rfunsim1 <- as.numeric(originData1[,1]);
  funsim1 <- as.numeric(originData1[,2]);
  rfunsimall1 <- as.numeric(originData1[,3]);
  funsimall1 <- as.numeric(originData1[,4]);
  mf1 <- as.numeric(originData1[,5]);
  bp1 <- as.numeric(originData1[,6]);
  cc1 <- as.numeric(originData1[,7]);
  #rfunsim2 <- as.numeric(originData2[,1]);
  #funsim2 <- as.numeric(originData2[,2]);
  #rfunsimall2 <- as.numeric(originData2[,3]);
  #funsimall2 <- as.numeric(originData2[,4]);
  #mf2 <- as.numeric(originData2[,5]);
  #bp2 <- as.numeric(originData2[,6]);
  #cc2 <- as.numeric(originData2[,7]);
  rfunsim3 <- as.numeric(originData3[,1]);
  funsim3 <- as.numeric(originData3[,2]);
  rfunsimall3 <- as.numeric(originData3[,3]);
  funsimall3 <- as.numeric(originData3[,4]);
  mf3 <- as.numeric(originData3[,5]);
  bp3 <- as.numeric(originData3[,6]);
  cc3 <- as.numeric(originData3[,7]);
  mybreak=c(0.0,0.6,1.0);
  #par(oma=c(5,5,5,5))
  #par(mar=c(6,7,5,3)+0.1)
  
  #pdf(file="rfunsim_record.pdf");
  jpeg(filename = "rfunsim_record.jpeg")
  par(mar=c(5,5,4,2)+0.1)
  boxplot(rfunsim1,rfunsim3,names=c("IsoRank-N","NetCoffee"),main="(c)",cex.lab=2.0,cex.axis=1.5,cex.main=2.0,ylab=expression(bar("rfunSim")))
  dev.off();
  #pdf(file="mf_record.pdf");
  jpeg(filename = "mf_record.jpeg")
  par(mar=c(5,5,4,2)+0.1)
  boxplot(mf1,mf3,names=c("IsoRank-N","NetCoffee"),main="(b)",cex.lab=2.0,cex.axis=1.5,cex.main=2.0,ylab=expression(bar("MFscore")))
  dev.off();
  #pdf(file="bp_record.pdf");
  jpeg(filename = "bp_record.jpeg")
  par(mar=c(5,5,4,2)+0.1)
  boxplot(bp1,bp3,names=c("IsoRank-N","NetCoffee"),main="(a)",cex.lab=2.0,cex.axis=1.5,cex.main=2.0,ylab=expression(bar("BPscore")))
  dev.off();
  #pdf(file="comparison_record.pdf");
  jpeg(filename = "comparison_record.jpeg")
  f1 <- hist(rfunsim1,mybreak,plot=FALSE);
  #f2 <- hist(rfunsim2,mybreak,plot=FALSE);
  f3 <- hist(rfunsim3,mybreak,plot=FALSE);
  d1 <- f1$density;
  #d2 <- f2$density;
  d3 <- f3$density;
  x1 <- c(d1[1]*60,d1[2]*40);
  #x2 <- c(d2[1]*60,d2[2]*40);
  x3 <- c(d3[1]*60,d3[2]*40);
  x <- cbind(x1,x3);
  par(mar=c(5,5,4,2)+0.1);
  barplot(x,beside=TRUE,names.arg=c("IsoRank-N","NetCoffee"),legend.text=c("FP","TP"), ylab="Percentage(%)",cex.names=1.5,cex.lab=2.0,,cex.main=2.0,cex.axis=1.5,main="(d)");
  dev.off();
  
#   pdf(file="isorankn_record.pdf");
#   browser();
#   f <- hist(rfunsim1,mybreak, freq=FALSE, plot=FALSE, col=mycol, xlab="average of rfunSim for each alignment record",main="");
#   barplot(f$density*20, col=mycol, names.arg=c("[0.0,0.6)","[0.6,1.0)"),ylim=c(0,50),ylab="Percetage(%)",xlab="rfunSim");
#   dev.off();
#   pdf(file="mnetaligner_record.pdf");
#   f <- hist(rfunsim2,mybreak, freq=FALSE, plot=FALSE, col=mycol, xlab="average of rfunSim for each alignment record",main="");
#   barplot(f$density*20, col=mycol, names.arg=c("[0.0,0.2)","[0.2,0.4)","[0.4,0.6)","[0.6,0.8)","[0.8,1.0)"),ylim=c(0,50),ylab="Percetage(%)",xlab="rfunSim");
#   dev.off();
#  pdf(file="tcoffee_record.pdf");
#  hist(rfunsim3,mybreak, freq=FALSE, xlab="average of rfunSim for each alignment record",main="");
#  dev.off();
}

extractrFunsim <- function(files)
{
  originData1 <- as.matrix(read.csv(files[1], sep = "\t", quote="", header=TRUE, comment="#"));
  originData2 <- as.matrix(read.csv(files[2], sep = "\t", quote="", header=TRUE, comment="#"));
  originData3 <- as.matrix(read.csv(files[3], sep = "\t", quote="", header=TRUE, comment="#"));
  rfunsim1 <- as.numeric(originData1[,2]);
  funsim1 <- as.numeric(originData1[,3]);
  rfunsimall1 <- as.numeric(originData1[,4]);
  funsimall1 <- as.numeric(originData1[,5]);
  mf1 <- as.numeric(originData1[,6]);
  bp1 <- as.numeric(originData1[,7]);
  cc1 <- as.numeric(originData1[,8]);
  rfunsim2 <- as.numeric(originData2[,2]);
  funsim2 <- as.numeric(originData2[,3]);
  rfunsimall2 <- as.numeric(originData2[,4]);
  funsimall2 <- as.numeric(originData2[,5]);
  mf2 <- as.numeric(originData2[,6]);
  bp2 <- as.numeric(originData2[,7]);
  cc2 <- as.numeric(originData2[,8]);
  rfunsim3 <- as.numeric(originData3[,2]);
  funsim3 <- as.numeric(originData3[,3]);
  rfunsimall3 <- as.numeric(originData3[,4]);
  funsimall3 <- as.numeric(originData3[,5]);
  mf3 <- as.numeric(originData3[,6]);
  bp3 <- as.numeric(originData3[,7]);
  cc3 <- as.numeric(originData3[,8]);
  mybreak=seq(from=0.0,to=1.0,by=0.2);  
  
  pdf(file="rfunsim.pdf");
  boxplot(rfunsim1,rfunsim2,rfunsim3,names=c("IsoRankN","M-NetAligner","tCoffee"),main="rfunsim")
  dev.off();
  pdf(file="funsim.pdf");
  boxplot(funsim1,funsim2,funsim3,names=c("IsoRankN","M-NetAligner","tCoffee"),main="funsim")
  dev.off();
  pdf(file="mf.pdf");
  boxplot(mf1,mf2,mf3,names=c("IsoRankN","M-NetAligner","tCoffee"),main="MF")
  dev.off();
  pdf(file="bp.pdf");
  boxplot(bp1,bp2,bp3,names=c("IsoRankN","M-NetAligner","tCoffee"),main="BP")
  dev.off();
  pdf(file="isorankn.pdf");
  hist(rfunsim1,mybreak, freq=FALSE,main="isorankn")
  dev.off();
  pdf(file="mnetaligner.pdf");
  hist(rfunsim2,mybreak, freq=FALSE,main="mnetaligner")
  dev.off();
  pdf(file="tcoffee.pdf");
  hist(rfunsim3,mybreak, freq=FALSE,main="tcoffee")
  dev.off();
}

plotConvergency <- function(filename, para)
{
  originData <- as.matrix(read.csv(filename, sep = "\t", quote="", header=FALSE, comment="!"));
  y <- originData[,1];
  size <- nrow(originData);
  x <- c(1:size);
  boxtext=paste(expression(alpha),"=",para);
  jpeg(file="./result/images/convergence.jpeg",quality=100,width=1200,height=1200,pointsize=34);
  plot(x,y,type="o",pch=1,col="red",cex.lab=3.0,cex.axis=2.5,cex.main=2.0);
  text(80000,500,boxtext,cex=2);
  dev.off();
}

scorefunction <- function(filename)
{
  browser();
  library("nls2");
  originData <- as.matrix(read.csv(filename, sep = "\t", quote="", header=FALSE, comment="!"));
  y <- originData[1:86,];
  x <- seq(from=-180,to=-10,by=2);
  fmode <- nls(y ~ beta1*x^3+beta2*x^2+beta3*x+beta4,start=list(beta1=0.0,beta2=0.0,beta3=0.0,beta4=3.0));
  fmode
}

format2table <- function(files)#filename1 inputfile filename2 outputfile
{
  j=1;
  precision = mat.or.vec(1,3)
  for(filename in files)
  {
    originData <- as.matrix(read.csv(filename, sep = "\t", quote="", header=FALSE, comment="#"));
    simpleData <- originData[,c(1,5,6)];
    nrow <- dim(simpleData)[1];
    tp  <- 0;
    metaData <- mat.or.vec(10,10);
    for(i in 1:nrow)
    {
      m = floor(simpleData[i,2]*10)+1;#mf
      n = floor(simpleData[i,3]*10)+1;#bp
      if(m>8 || n>6) tp= tp+1;
      #metaData[m,n]=metaData[m,n]+1;
    }
    precision[j] = (tp*100)/nrow;
    j=j+1;
  }
  precision = precision[c(1,3)]
  jpeg(filename="precision.jpeg");
  par(mar=c(5,5,4,2)+0.1);
  barplot(precision,beside=TRUE,names.arg=c("IsoRank-N","NetCoffee"), ylab="Percentage(%)",cex.names=1.5,cex.lab=2.0,,cex.main=2.0,cex.axis=1.5,main="(d)");
  dev.off();
}

eraseGoterm <- function(filelist)#erase IEA and ISS
{
  for(filename in filelist)
  {
    outfilename=paste(filename,"reduced",sep=".");
    #outfilename=paste(filename,"evals",sep=".");
    originData <- as.matrix(read.csv(filename, sep = "\t", quote="", header=FALSE, comment="!"));
    DataView <- originData[,c(2,5,7,9)];
    nrow <- dim(originData)[1];
    ncol <- dim(originData)[2];
    browser()
    newData <- originData;
    j=1;
    for(i in 1:nrow)
    {
      if(originData[i,7]=="IEA" || originData[i,7]=="ISS")
      {
      }else
      {
        newData[j,] <- originData[i,];
        j=j+1;
      }
    } 
    browser()
    #DataView <- originData[,c(1,2,12)]
    write.table(newData[1:j,],file=outfilename,sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE);
  }
}

drawPValue <- function()
{
  y1 <- as.matrix(read.csv("./isorankn_pvalue.txt", sep = "\t", quote="", header=FALSE, comment="!"));
  y2 <- as.matrix(read.csv("./tcoffee_pvalue.txt", sep = "\t", quote="", header=FALSE, comment="!"));
  x <- 1:2000;
  boxcol=c("red","blue")
  jpeg("pvalue.jpeg");
  #par(mar=c(5,5,4,2)+0.1);
  Y1 <- sort(y1,decreasing=FALSE);
  Y2 <- sort(y2,decreasing=FALSE);
  plot(x,Y2[1:2000],log="y",type="o",pch=1,col="blue")
  lines(x,Y1[1:2000],type="o",pch=4,col="red")
  legend("topleft",legend=c("isorankn","T-Coffee"),lwd=2.5,text.col=c("red","blue"),col=c("red","blue"),pch=c(4,1))
  dev.off();
  jpeg("pvalue_box.jpeg");
  boxplot(Y1,Y2,col=boxcol,names=c("M-NetAligner","T-Coffee"),cex.lab=2.0,cex.axis=1.5,cex.main=2.0,log="y")
  dev.off();
}

calculateMatchSet_i <- function(filename)
{
    browser();
	originData <- as.matrix(read.csv(filename, sep = "\t", quote="", header=FALSE, comment="#"));
	data1_qualified <- originData[1:6,];
	data1_all <- originData[7:12,];
	data2_qualified <- originData[13:20,];
	data2_all <- originData[21:28,];
	ratio1 = data1_qualified/data1_all;
	ratio2 = data2_qualified/data2_all;
	print(ratio1)
	print(ratio2)
}

drawMatchSet_i <- function(filename,numspecies)
# This function doesn't work now. Because the match-set-i.data was changed.
# The percentage = # of qualified match-sets/ all match-sets including match-sets containing uncharacterized proteins
{
  originData <- as.matrix(read.csv(filename,sep="\t",quote="",header=FALSE,comment="#"));
  browser();
  result=originData[,1]/(originData[,2]);
  net=result[seq(1,14,by=2)];
  iso=result[seq(2,14,by=2)];
  x <- net;
  x <- rbind(x,iso);
  y <- x*100;
  browser();
  pdf(file="./result/images/precision.pdf");
  par(mar=c(5,5,4,2)+0.1)
  barplot(y,beside=TRUE,names.arg=c(0.0,0.3,0.4,0.5,0.6,0.7,1.0), ylab="Percentage(%)",col=c("gray","white"),cex.names=1.5,cex.lab=2.0,,cex.main=2.0,cex.axis=1.5,ylim=c(0,60));
  #legend("topleft",legend=c("NetCoffee","IsoRank-N"),lwd=2.5,col=c("gray"));
  dev.off()
}

drawBoxplotMulFunsim <- function()
{
	ALPHA <- c("0.0",0.3,0.4,0.5,0.6,0.7,"1.0");
	filelist <- c();
	for (i in ALPHA)
	{ 
	    filename1 =paste("./result/","five_species",sep="");
		filename =paste(filename1,"/netcoffee/alpha_",sep="");
		filename2 =paste(filename1,"/isorankn/alpha_",sep="");
		filename1 =paste(filename,i,sep="");
		filename1 =paste(filename1,"/aveFunSim.result",sep="");
		filename2 =paste(filename2,i,sep="");
		filename2 =paste(filename2,"/aveFunSim.result",sep="");
		filelist <- rbind(filelist,filename1);
		filelist <- rbind(filelist,filename2);
	}
	browser();
	originData1 <- as.matrix(read.csv(filelist[1],sep="\t",quote="",header=FALSE,comment="#"));
	originData2 <- as.matrix(read.csv(filelist[2],sep="\t",quote="",header=FALSE,comment="#"));
	originData3 <- as.matrix(read.csv(filelist[3],sep="\t",quote="",header=FALSE,comment="#"));
	originData4 <- as.matrix(read.csv(filelist[4],sep="\t",quote="",header=FALSE,comment="#"));
	originData5 <- as.matrix(read.csv(filelist[5],sep="\t",quote="",header=FALSE,comment="#"));
	originData6 <- as.matrix(read.csv(filelist[6],sep="\t",quote="",header=FALSE,comment="#"));
	originData7 <- as.matrix(read.csv(filelist[7],sep="\t",quote="",header=FALSE,comment="#"));
	originData8 <- as.matrix(read.csv(filelist[8],sep="\t",quote="",header=FALSE,comment="#"));
	originData9 <- as.matrix(read.csv(filelist[9],sep="\t",quote="",header=FALSE,comment="#"));
	originData10 <- as.matrix(read.csv(filelist[10],sep="\t",quote="",header=FALSE,comment="#"));
	originData11 <- as.matrix(read.csv(filelist[11],sep="\t",quote="",header=FALSE,comment="#"));
	originData12 <- as.matrix(read.csv(filelist[12],sep="\t",quote="",header=FALSE,comment="#"));
	originData13 <- as.matrix(read.csv(filelist[13],sep="\t",quote="",header=FALSE,comment="#"));
	originData14 <- as.matrix(read.csv(filelist[14],sep="\t",quote="",header=FALSE,comment="#"));
	
	data1 <- originData1[,c(1,5,6)];
	data2 <- originData2[,c(1,5,6)];
	data3 <- originData3[,c(1,5,6)];
	data4 <- originData4[,c(1,5,6)];
	data5 <- originData5[,c(1,5,6)];
	data6 <- originData6[,c(1,5,6)];
	data7 <- originData7[,c(1,5,6)];
	data8 <- originData8[,c(1,5,6)];
	data9 <- originData9[,c(1,5,6)];
	data10 <- originData10[,c(1,5,6)];
	data11 <- originData11[,c(1,5,6)];
	data12 <- originData12[,c(1,5,6)];
	data13 <- originData13[,c(1,5,6)];
	data14 <- originData14[,c(1,5,6)];
  
   browser();

	pdf(filename="./result/images/rfunsim_record.pdf");
	par(mar=c(5,5,4,2)+0.1)
	boxplot(data1[,1],data2[,1],
			data3[,1],data4[,1],
			data5[,1],data6[,1],
			data7[,1],data8[,1],
			data9[,1],data10[,1],
			ylim=c(0,1.2),
			names=c(0.3,0.3,0.4,0.4,0.5,0.5,0.6,0.6,0.7,0.7),
			border = c("black","darkgrey","black","darkgrey","black","darkgrey","black","darkgrey","black","darkgrey"),
			cex.lab=2.0,cex.axis=1.5,cex.main=2.0,ylab=expression(bar("rfunSim")));
			legend("topleft",legend=c("NetCoffee","IsoRank-N"),lwd=2.5,col=c("red","blue"),pch=c(22,22));
	dev.off()
  
  pdf(filename="./result/images/mf_record.pdf");
  par(mar=c(5,5,4,2)+0.1)
	boxplot(data1[,2],data2[,2],
			data3[,2],data4[,2],
			data5[,2],data6[,2],
			data7[,2],data8[,2],
			data9[,2],data10[,2],
			ylim=c(0,1.2),
			names=c(0.3,0.3,0.4,0.4,0.5,0.5,0.6,0.6,0.7,0.7),
			border = c("red","blue","red","blue","red","blue","red","blue","red","blue"),
			cex.lab=2.0,cex.axis=1.5,cex.main=2.0,ylab=expression(bar("MFscore")));
			legend("topleft",legend=c("NetCoffee","IsoRank-N"),lwd=2.5,col=c("red","blue"),pch=c(22,22));
	dev.off()
  
  pdf(filename="./result/images/bp_record.pdf");
  par(mar=c(5,5,4,2)+0.1)
	boxplot(data1[,3],data2[,3],
			data3[,3],data4[,3],
			data5[,3],data6[,3],
			data7[,3],data8[,3],
			data9[,3],data10[,3],
			ylim=c(0,1.2),
			names=c(0.3,0.3,0.4,0.4,0.5,0.5,0.6,0.6,0.7,0.7),
			border = c("red","blue","red","blue","red","blue","red","blue","red","blue"),
			cex.lab=2.0,cex.axis=1.5,cex.main=2.0,ylab=expression(bar("BPscore")));
			legend("topleft",legend=c("NetCoffee","IsoRank-N"),lwd=2.5,col=c("red","blue"),pch=c(22,22));
	dev.off()
	print(filelist);
}

drawEdgeScoreDistri <- function()
{
	originData <- as.matrix(read.csv("./result/five_species/netcoffee/alpha_0.0/scoreRecords.txt",sep="\t",quote="",header=FALSE,comment="#"));
  browser();
  myseq <- as.numeric(originData[,2]);
  top <- as.numeric(originData[,5]);
  amyseq <- ( myseq -min(myseq))/(max(myseq)-min(myseq));
  atop <- (top-min(top))/(max(top)-min(top));
  pdf("./result/images/edgescore_hist1.pdf")
  h1<-hist(amyseq,breaks=10,plot=FALSE);
  count1 <-h1$counts;
  count1 <- 100*count1/sum(count1);
  barplot(count1,names.arg=c(seq(from=0.1,to=1.0,by=0.1)),xlab="Sequence-based score",ylab="Probability(%)",main="(a)")
  dev.off();
  pdf("./result/images/edgescore_hist2.pdf")
  h2<-hist(atop,breaks=10,plot=FALSE);
  count2 <-h2$counts;
  count2 <- 100*count2/sum(count2);
  barplot(count2,names.arg=c(seq(from=0.1,to=1.0,by=0.1)),xlab="Topology-based score",ylab="Probability(%)",main="(b)")  
  dev.off();
  pdf("./result/images/edgescore_box1.pdf");
  boxplot(amyseq,atop,names=c("seqence-based score","topology-based score"));
  dev.off();
  pdf("./result/images/edgescore_hist3.pdf");
  atop <- atop^0.10;
  h3 <- hist(atop,breaks=10,plot=FALSE);
  count3 <- h3$counts;
  count3 <- 100*count3/sum(count3);
   barplot(count3,names.arg=c(seq(from=0.1,to=1.0,by=0.1)),xlab="Rescaled topology-based score",ylab="Probability(%)",main="(c)")
  dev.off();
   pdf("./result/images/edgescore_box2.pdf");
  boxplot(amyseq,atop,names=c("seqence-based score","topology-based score"));
  dev.off();
}

drawBoxplotZeroMatchSet <- function()
{
  matchset <- c("matchsets-three.txt", "matchsets-four.txt","matchsets-five.txt");
  figlabel <- c("(a) i=","(b) i=","(c) i=","(d) i=","(e) i=","(f) i=","(g) i=","(h) i=","(i) i=");
  #for graemlin data
  species=3;
  #for graemlin data
  fignum=1;
  for (mat in matchset )
  {
    ALPHA <- c(0.0,0.3,0.4,0.5,0.6,0.7);
    filelist <- c();
	  for (i in ALPHA)
	  { 
	      filename1 =paste("./result/","five_species",sep="");
		  filename =paste(filename1,"/netcoffee/alpha_",sep="");
		  filename1 =paste(filename,i,sep="");
		  filename1 =paste(filename1,mat,sep="/");
		  filelist <- rbind(filelist,filename1);
	  }
    print(filelist);
    originData1 <- as.matrix(read.csv(filelist[1],sep="\t",quote="",header=FALSE,comment="#"));
  	originData2 <- as.matrix(read.csv(filelist[2],sep="\t",quote="",header=FALSE,comment="#"));
  	originData3 <- as.matrix(read.csv(filelist[3],sep="\t",quote="",header=FALSE,comment="#"));
  	originData4 <- as.matrix(read.csv(filelist[4],sep="\t",quote="",header=FALSE,comment="#"));
  	originData5 <- as.matrix(read.csv(filelist[5],sep="\t",quote="",header=FALSE,comment="#"));
  	originData6 <- as.matrix(read.csv(filelist[6],sep="\t",quote="",header=FALSE,comment="#"));
    
    data1 <- originData1[,c(1,5,6)];
    data2 <- originData2[,c(1,5,6)];
  	data3 <- originData3[,c(1,5,6)];
  	data4 <- originData4[,c(1,5,6)];
  	data5 <- originData5[,c(1,5,6)];
  	data6 <- originData6[,c(1,5,6)];

  	xcolnames = c(0.0,0.3,0.4,0.5,0.6,0.7);

    boxtext=paste(figlabel[fignum], species, sep="");
  	fignum=fignum+1;
  	output3=paste("./result/images/bp_net_",mat,sep="");
    output3=paste(output3,"jpeg",sep=".")
    jpeg(filename=output3,quality=100,width=1200,height=1200,pointsize=34);
    par(mar=c(5,5,4,2)+0.1);
	boxplot(data1[,3],data2[,3],
			data3[,3],data4[,3],
			data5[,3],data6[,3],
			ylim=c(0,1.2),
			#xaxt="n",
			names=c(0.0, 0.3,0.4,0.5,0.6,0.7),
			border=c("black","antiquewhite4","black","antiquewhite4","black","antiquewhite4"),
			col= c("gray","gray","gray","gray","gray","gray"),
			cex.lab=2.0,cex.axis=1.5,cex.main=2.0,ylab=expression(bar("BPscore")));
			#legend("topleft",legend=c("NetCoffee","IsoRank-N"),lwd=2.5,col=c("gray","black"),pch=c(15,22));
			text(5,1.1,boxtext,cex=2);
			#axis(1, at = c(1.5,3.5,5.5,7.5,9.5), labels = xcolnames, cex.axis=1.5);
	dev.off()
	
	boxtext=paste(figlabel[fignum], species, sep="");
  	fignum=fignum+1;
    
    browser();
  
    output2=paste("./result/images/mf_net_",mat,sep="");
    output2=paste(output2,"jpeg",sep=".")
  jpeg(filename=output2,quality=100,width=1200,height=1200,pointsize=34);
  par(mar=c(5,5,4,2)+0.1)
	boxplot(data1[,2],data2[,2],
			data3[,2],data4[,2],
			data5[,2],data6[,2],
			ylim=c(0,1.2),
			#xaxt="n",
			names=c(0.0,0.3,0.4,0.5,0.6,0.7),
			border = c("black","antiquewhite4","black","antiquewhite4","black","antiquewhite4"),
			col= c("gray","gray","gray","gray","gray","gray"),
			cex.lab=2.0,cex.axis=1.5,cex.main=2.0,ylab=expression(bar("MFscore")));
			#legend("topleft",legend=c("NetCoffee","IsoRank-N"),lwd=2.5,col=c("black","antiquewhite4"),pch=c(15,22));
			text(5,1.1,boxtext,cex=2);
			#axis(1, at = c(1.5,3.5,5.5,7.5,9.5), labels = xcolnames, cex.axis=1.5);
	dev.off()
	
	boxtext=paste(figlabel[fignum], species, sep="");
  	fignum=fignum+1;    
    output1=paste("./result/images/rfunsim_net_",mat,sep="");
    output1=paste(output1,"jpeg",sep=".")
    jpeg(filename=output1,quality=100,width=1200,height=1200,pointsize=34);
  par(mar=c(5,5,4,2)+0.1)
	boxplot(data1[,1],data2[,1],
			data3[,1],data4[,1],
			data5[,1],data6[,1],
			ylim=c(0,1.2),
			#xaxt="n",
			names=c(0.0,0.3,0.4,0.5,0.6,0.7),
			border = c("black","antiquewhite4","black","antiquewhite4","black","antiquewhite4"),
			col= c("gray","gray","gray","gray","gray","gray"),
			cex.lab=2.0,cex.axis=1.5,cex.main=2.0,ylab=expression(bar("rfunSim")));
			#legend("topleft",legend=c("NetCoffee","IsoRank-N"),lwd=2.5,col=c("black","antiquewhite4"),pch=c(15,22));
			text(5,1.1,boxtext,cex=2);
	#axis(1, at = c(1.5,3.5,5.5,7.5,9.5), labels = xcolnames, cex.axis=1.5);
	dev.off()

	species=species+1;
    
  }
}

drawBoxplotTestLogModel <- function()
{
  filename1="./result/investigation/toplog/aveFunSim.result";
  filename2="./result/investigation/toplogratio/aveFunSim.result";
  originData1 <- as.matrix(read.csv(filename1,sep="\t",quote="",header=FALSE,comment="#"));
  originData2 <- as.matrix(read.csv(filename2,sep="\t",quote="",header=FALSE,comment="#"));
  browser();
  nrow1=dim(originData1)[1];
  nrow2=dim(originData2)[1];
  data1 <- as.numeric(originData1[,c(1,5,6)]);
  data2 <- as.numeric(originData2[,c(1,5,6)]);  
  dim(data1) <- c(nrow1,3);
  dim(data2) <- c(nrow2,3);
  xcolnames = c("-Log","LogRatio");
  figlabel <- c("(a)","(b)","(c)");
  fignum=1;
  browser();
  
  boxtext=figlabel[fignum];
    fignum=fignum+1;
  	output3=paste("./result/investigation/bp","pdf",sep=".");
  pdf(file=output3);
  par(mar=c(5,5,4,2)+0.1)
	boxplot(data1[,3],data2[,3],
			ylim=c(0,1.2),
			xaxt="n",
			#names=c(0.3,0.4,0.5,0.6,0.7),
			border=c("black","antiquewhite4"),
			col= c("gray","white"),
			cex.lab=2.0,cex.axis=1.5,cex.main=2.0,ylab=expression(bar("BPscore")));
			#legend("topleft",legend=c("NetCoffee","IsoRank-N"),lwd=2.5,col=c("gray","black"),pch=c(15,22));
			text(5,1.1,boxtext,cex=2);
			axis(1, at = c(1.0,2.0), labels = xcolnames, cex.axis=1.5);
	dev.off()
browser();
	
	boxtext=figlabel[fignum];
  	fignum=fignum+1;
  
    output2=paste("./result/investigation/mf","pdf",sep=".");
  pdf(file=output2);
  par(mar=c(5,5,4,2)+0.1)
	boxplot(data1[,2],data2[,2],
			ylim=c(0,1.2),
			xaxt="n",
			#names=c(0.3,0.4,0.5,0.6,0.7),
			border = c("black","antiquewhite4"),
			col= c("gray","white"),
			cex.lab=2.0,cex.axis=1.5,cex.main=2.0,ylab=expression(bar("MFscore")));
			#legend("topleft",legend=c("NetCoffee","IsoRank-N"),lwd=2.5,col=c("black","antiquewhite4"),pch=c(15,22));
			text(5,1.1,boxtext,cex=2);
			axis(1, at = c(1.0,2.0), labels = xcolnames, cex.axis=1.5);
	dev.off()
browser();
	
	boxtext=figlabel[fignum];
  	fignum=fignum+1;    
    output1=paste("./result/investigation/rfunsim","pdf",sep=".");
    pdf(file=output1);
  par(mar=c(5,5,4,2)+0.1)
	boxplot(data1[,1],data2[,1],
			ylim=c(0,1.2),
			xaxt="n",
			#names=c(0.3,0.4,0.5,0.6,0.7),
			border = c("black","antiquewhite4"),
			col= c("gray","white"),
			cex.lab=2.0,cex.axis=1.5,cex.main=2.0,ylab=expression(bar("rfunSim")));
			#legend("topleft",legend=c("NetCoffee","IsoRank-N"),lwd=2.5,col=c("black","antiquewhite4"),pch=c(15,22));
			text(5,1.1,boxtext,cex=2);
	axis(1, at = c(1.0,2.0), labels = xcolnames, cex.axis=1.5);
	dev.off();
}

drawBoxplotMatchSet_i <- function()
{
  matchset <- c("matchsets-three.txt", "matchsets-four.txt","matchsets-five.txt");
  #
  figlabel <- c("(a) i=","(b) i=","(c) i=","(d) i=","(e) i=","(f) i=","(g) i=","(h) i=","(i) i=");
  #for graemlin data
  species=3;
  #for graemlin data
  fignum=1;
  for (mat in matchset )
  {
    ALPHA <- c("0.0",0.3,0.4,0.5,0.6,0.7,"1.0");
    filelist <- c();
	  for (i in ALPHA)
	  { 
	    filename1 =paste("./result/","five_species",sep="");
		  filename =paste(filename1,"/netcoffee/alpha_",sep="");
		  filename2 =paste(filename1,"/isorankn/alpha_",sep="");
		  filename1 =paste(filename,i,sep="");
		  filename1 =paste(filename1,mat,sep="/");
		  filename2 =paste(filename2,i,sep="");
		  #for graemlin data
		  #filename2="./result/six_species/graemlin";
		  #for graemlin data
		  filename2 =paste(filename2,mat,sep="/");
		  filelist <- rbind(filelist,filename1);
		  filelist <- rbind(filelist,filename2);
	  }
    print(filelist);
    browser();
    originData1 <- as.matrix(read.csv(filelist[1],sep="\t",quote="",header=FALSE,comment="#"));
  	originData2 <- as.matrix(read.csv(filelist[2],sep="\t",quote="",header=FALSE,comment="#"));
  	originData3 <- as.matrix(read.csv(filelist[3],sep="\t",quote="",header=FALSE,comment="#"));
  	originData4 <- as.matrix(read.csv(filelist[4],sep="\t",quote="",header=FALSE,comment="#"));
  	originData5 <- as.matrix(read.csv(filelist[5],sep="\t",quote="",header=FALSE,comment="#"));
  	originData6 <- as.matrix(read.csv(filelist[6],sep="\t",quote="",header=FALSE,comment="#"));
  	originData7 <- as.matrix(read.csv(filelist[7],sep="\t",quote="",header=FALSE,comment="#"));
  	originData8 <- as.matrix(read.csv(filelist[8],sep="\t",quote="",header=FALSE,comment="#"));
  	originData9 <- as.matrix(read.csv(filelist[9],sep="\t",quote="",header=FALSE,comment="#"));
  	originData10 <- as.matrix(read.csv(filelist[10],sep="\t",quote="",header=FALSE,comment="#"));
    originData11 <- as.matrix(read.csv(filelist[11],sep="\t",quote="",header=FALSE,comment="#"));
    originData12 <- as.matrix(read.csv(filelist[12],sep="\t",quote="",header=FALSE,comment="#"));
    originData13 <- as.matrix(read.csv(filelist[13],sep="\t",quote="",header=FALSE,comment="#"));
    originData14 <- as.matrix(read.csv(filelist[14],sep="\t",quote="",header=FALSE,comment="#"));
    
    data1 <- originData1[,c(1,5,6)];
    data2 <- originData2[,c(1,5,6)];
  	data3 <- originData3[,c(1,5,6)];
  	data4 <- originData4[,c(1,5,6)];
  	data5 <- originData5[,c(1,5,6)];
  	data6 <- originData6[,c(1,5,6)];
  	data7 <- originData7[,c(1,5,6)];
  	data8 <- originData8[,c(1,5,6)];
  	data9 <- originData9[,c(1,5,6)];
  	data10 <- originData10[,c(1,5,6)];
    data11 <- originData11[,c(1,5,6)];
    data12 <- originData12[,c(1,5,6)];
    data13 <- originData13[,c(1,5,6)];
    data14 <- originData14[,c(1,5,6)];

  	#xcolnames = seq(from=0.3,to=0.7,by=0.1);
    xcolnames = c(0.0, 0.3, 0.4, 0.5, 0.6, 0.7, 1.0);

    boxtext=paste(figlabel[fignum], species, sep="");
  	fignum=fignum+1;
  	output3=paste("./result/images/bp_",mat,sep="");
    output3=paste(output3,"pdf",sep=".")
  pdf(file=output3);
  par(mar=c(5,5,4,2)+0.1)
	boxplot(data1[,3],data2[,3],
			data3[,3],data4[,3],
			data5[,3],data6[,3],
			data7[,3],data8[,3],
			data9[,3],data10[,3],
      data11[,3],data12[,3],
      data13[,3],data14[,3],
			ylim=c(0,1.2),
			xaxt="n",
			#names=c(0.0,0.3,0.4,0.5,0.6,0.7,1.0),
			border=c("black","antiquewhite4","black","antiquewhite4","black","antiquewhite4","black","antiquewhite4","black","antiquewhite4","black","antiquewhite4","black","antiquewhite4"),
			col= c("gray","white","gray","white","gray","white","gray","white","gray","white","gray","white","gray","white"),
			cex.lab=2.0,cex.axis=1.5,cex.main=2.0,ylab=expression(bar("BPscore")));
			#legend("topleft",legend=c("NetCoffee","IsoRank-N"),lwd=2.5,col=c("gray","black"),pch=c(15,22));
			text(5,1.1,boxtext,cex=2);
			axis(1, at = c(1.5,3.5,5.5,7.5,9.5,11.5,13.5), labels = xcolnames, cex.axis=1.5);
	dev.off()
    browser();
	
	boxtext=paste(figlabel[fignum], species, sep="");
  	fignum=fignum+1;
  
    output2=paste("./result/images/mf_",mat,sep="");
    output2=paste(output2,"pdf",sep=".")
  pdf(file=output2);
  par(mar=c(5,5,4,2)+0.1)
	boxplot(data1[,2],data2[,2],
			data3[,2],data4[,2],
			data5[,2],data6[,2],
			data7[,2],data8[,2],
			data9[,2],data10[,2],
      data11[,2],data12[,2],
        data13[,2],data14[,2],
			ylim=c(0,1.2),
			xaxt="n",
			#names=c(0.3,0.4,0.5,0.6,0.7),
			border = c("black","antiquewhite4","black","antiquewhite4","black","antiquewhite4","black","antiquewhite4","black","antiquewhite4","black","antiquewhite4","black","antiquewhite4","black","antiquewhite4","black","antiquewhite4"),
			col= c("gray","white","gray","white","gray","white","gray","white","gray","white","gray","white","gray","white"),
			cex.lab=2.0,cex.axis=1.5,cex.main=2.0,ylab=expression(bar("MFscore")));
			#legend("topleft",legend=c("NetCoffee","IsoRank-N"),lwd=2.5,col=c("black","antiquewhite4"),pch=c(15,22));
			text(5,1.1,boxtext,cex=2);
			axis(1, at = c(1.5,3.5,5.5,7.5,9.5,11.5,13.5), labels = xcolnames, cex.axis=1.5);
	dev.off()
	
	boxtext=paste(figlabel[fignum], species, sep="");
  	fignum=fignum+1;    
    output1=paste("./result/images/rfunsim_",mat,sep="");
    output1=paste(output1,"pdf",sep=".")
    pdf(file=output1);
  par(mar=c(5,5,4,2)+0.1)
	boxplot(data1[,1],data2[,1],
			data3[,1],data4[,1],
			data5[,1],data6[,1],
			data7[,1],data8[,1],
			data9[,1],data10[,1],
      data11[,1],data12[,1],
      data13[,1],data14[,1],
			ylim=c(0,1.2),
			xaxt="n",
			#names=c(0.3,0.4,0.5,0.6,0.7),
			border = c("black","antiquewhite4","black","antiquewhite4","black","antiquewhite4","black","antiquewhite4","black","antiquewhite4","black","antiquewhite4","black","antiquewhite4"),
			col= c("gray","white","gray","white","gray","white","gray","white","gray","white","gray","white","gray","white"),
			cex.lab=2.0,cex.axis=1.5,cex.main=2.0,ylab=expression(bar("rfunSim")));
			#legend("topleft",legend=c("NetCoffee","IsoRank-N"),lwd=2.5,col=c("black","antiquewhite4"),pch=c(15,22));
			text(5,1.1,boxtext,cex=2);
	axis(1, at = c(1.5,3.5,5.5,7.5,9.5,11.5,13.5), labels = xcolnames, cex.axis=1.5);
	dev.off()

	species=species+1;
    
  }
}

