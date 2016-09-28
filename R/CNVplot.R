CNVplot<- function (df,Start,End,copyNumber,genome,title,yLabel) {
  opar<-par()
  if(missing(copyNumber)){ copyNumber<-FALSE}
  if(missing(Start)){ Start<-FALSE}
  if(missing(End)){End<-FALSE}
  if(missing(yLabel)){yLabel=""}
  if(missing(genome)){ genome="hg19"}
  if (genome != 19 & genome != 38 & genome != "hg19" & genome != "hg38"){ stop("FATAL ERROR: genome must be either hg19 or hg38")}
  columns=ncol(df)
  trio=FALSE
  single=FALSE
  if(columns==3){ single=TRUE}
  else if(columns==5) { trio=TRUE}
  else { stop("FATAL ERROR: input data frame must have 3 or 5 columns")}
  df<-df[order(df[,2],decreasing = F),]
  if(nrow(df)==0){ stop("FATAL ERROR: input data frame has no rows") }
  offset<-NULL
  offsetTitle<-NULL
  if(df[1,2] > 1000000) {
    offset<-1000000
    offsetTitle<-"(Mb)"
  }
  if(df[1,2] < 1000000 & df[1,2] > 1000) {
    offset<-1000
    offsetTitle<-"(kb)"
  }
  if(df[1,2] <= 1000) {
    offset<-1
    offsetTitle<-"(bp)"
  }
  chr <- df[1,1]
  start<-min(df[,2],na.rm = T) # no offset
  end<-max(df[,2],na.rm = T)   # no offset
  df[,2]<-df[,2]/offset # offset the df
  #	Parameters for plotting
  plot.start<-min(df[,2],na.rm = T) #offset start
  plot.end<-max(df[,2],na.rm = T) #offset end
  xpos<-seq(plot.start,plot.end,length.out = 4)
  y.plot<-NULL
  maxim<-max(df[,3],na.rm=T)
  minin<-min(df[,3],na.rm=T)
  if(trio==TRUE){
    maxim<-max(df[,c(3,4,5)],na.rm = T)
    minin<-min(df[,c(3,4,5)],na.rm = T)
  }
  minin<-abs(minin)
  if(minin >= maxim){ y.plot<- minin}
  if(maxim > minin){ y.plot <- maxim}
  y.plot = ceiling(y.plot)+1
  mirror.y<- -1*y.plot
  # colors for copy number
  if(copyNumber != FALSE){
    if(copyNumber == "del" | copyNumber == "DEL" | (copyNumber < 2 & is.numeric(copyNumber))){cols<-"#ff6666"}
    if(copyNumber == "dup" | copyNumber == "DUP" | (copyNumber >= 2 & is.numeric(copyNumber))){cols<-"#99ccff"}
  }
  else { cols<-"#000000"}
  if (Start != FALSE & End != FALSE & copyNumber !=FALSE) {
    Start <- Start/offset
    End <- End/offset
  }
  else { cnv.points <- FALSE }
  rug <- df[, 2]
  temp.rug <- seq(1:length(rug))
  rug <- as.data.frame(cbind(rug, temp.rug))
  rug$temp.rug <- 2
  if (nrow(df) == 1) {
    plot.start <- plot.start - (1/offset)
    plot.end <- plot.end + (1/offset)
    xpos <- plot.start + (1/offset)
  }
  genomic.space <- NULL
  x <- seq(plot.start, plot.end)
  y <- seq(1:length(x))
  genomic.space <- as.data.frame(cbind(x, y))
  userGenes=NULL
  if(genome == 19 | genome == "hg19"){
    data(hg19_exons)
    chrTemp <- chr
    chrTemp <- gsub("chr", "", chrTemp)
    chrQuery <- paste("chr", chrTemp, sep = "")
    userGenes <- subset(hg19_exons, chr == chrQuery & ((txStart > start & txEnd < end) | (txStart < start & txEnd < end & txEnd > start) | (txStart > start & txEnd > end & txStart < end) | (txStart < start & txStart < end & txEnd > start & txEnd > end)))
    hg19_exons<-NULL
    rm(hg19_exons)
  }
  else if (genome == 38 | genome == "hg38"){
    data(hg38_exons)
    chrTemp <- chr
    chrTemp <- gsub("chr", "", chrTemp)
    chrQuery <- paste("chr", chrTemp, sep = "")
    userGenes <- subset(hg38_exons, chr == chrQuery & ((txStart > start & txEnd < end) | (txStart < start & txEnd < end & txEnd > start) | (txStart > start & txEnd > end & txStart < end) | (txStart < start & txStart < end & txEnd > start & txEnd > end)))
    hg38_exons<-NULL
    rm(hg38_exons)
  }
  if (nrow(userGenes) == 0) {
    if (nrow(df) == 1) {
      plot.start <- plot.start - (1000/offset)
      plot.end <- plot.end + (1000/offset)
    }
    par(mai = c(0.1, 0.7, 0.1, 0.1), omi = c(1, 0.5, 1, 1))
    layout(matrix(c(1, 2, 3), 3, 1, byrow = T), heights = c(0.5,8.5, 1.5))
    if(trio==TRUE){ trioPlot(df, rug, plot.start, plot.end, mirror.y,y.plot, title,yLabel) }
    else if(single==TRUE){ singlePlot(df,rug,plot.start,plot.end,mirror.y,y.plot,title,yLabel)}
    xtitle <- paste("genome position", chrQuery, offsetTitle,sep = " ")
    plot(genomic.space, type = "n", xlab = "", ylim = c(-35,35), ylab = "", yaxt = "n", xaxt = "n", xlim = c(plot.start,plot.end))
    segments(x0 = Start, y0 = 0, x1 = End, y1 = 0,col = cols, lwd = 10, lend = 2)
    axis(1, at = xpos, labels = formatC(xpos, digits = 3,format = "f"))
    mtext(text = xtitle, SOUTH <- 1, line = 3)
    layout(matrix(c(1), 1, 1, byrow = T))
    suppressWarnings(par(opar))
  }
  else {
    userGenes <- as.data.frame(userGenes)
    userGenes <- userGenes[order(userGenes$txStart), ]
    userGenes$offset <- NULL
    counter <- 0
    for (i in 1:nrow(userGenes)) {
      if (counter == 0) { userGenes$offset[i] <- -5 }
      if (counter == 1) { userGenes$offset[i] <- 5 }
      if (counter == 2) { userGenes$offset[i] <- -10.67}
      if (counter == 3) { userGenes$offset[i] <- 10.67 }
      if (counter == 4) { userGenes$offset[i] <- -17.77 }
      if (counter == 5) { userGenes$offset[i] <- 17.77 }
      if (counter == 6) { userGenes$offset[i] <- -24.89 }
      if (counter == 7) { userGenes$offset[i] <- 24.89 }
      if (counter == 8) { userGenes$offset[i] <- -32 }
      if (counter == 9) { userGenes$offset[i] <- 32 }
      counter <- counter + 1
      if (counter == 9) { counter = 0 }
    }
    ySeg <- 0
    par(mai = c(0.1, 0.7, 0.1, 0.1), omi = c(1, 0.5, 1, 1))
    if (nrow(userGenes) > 3) {
      layout(matrix(c(1, 2, 3), 3, 1, byrow = T), heights = c(0.5,5, 4.5))
      ySeg <- 27
    }
    else {
      layout(matrix(c(1, 2, 3), 3, 1, byrow = T), heights = c(0.5,7, 3))
      ySeg <- 20
    }
    if (nrow(userGenes) > 6) {
      layout(matrix(c(1, 2, 3), 3, 1, byrow = T), heights = c(0.5,4, 5.5))
      ySeg <- 33.5
    }
    if(trio==TRUE){ trioPlot(df, rug, plot.start, plot.end, mirror.y,y.plot, title,yLabel) }
    else if(single==TRUE){ singlePlot(df,rug,plot.start,plot.end,mirror.y,y.plot,title,yLabel)}
    xtitle <- paste("genome position", chrQuery, offsetTitle,sep = " ")
    plot(genomic.space, type = "n", xlab = "", ylim = c(-35,35), ylab = "", yaxt = "n", xaxt = "n", xlim = c(plot.start,plot.end))
    genePlotter(userGenes, start, end, xpos, offset)
    segments(x0 = Start, y0 = ySeg, x1 = End, y1 = ySeg,col = cols, lwd = 10, lend = 2)
    mtext(text = xtitle, SOUTH <- 1, line = 3)
    axis(1, at = xpos, labels = formatC(xpos, digits = 3,format = "f"))
    layout(matrix(c(1), 1, 1, byrow = T))
    suppressWarnings(par(opar))
  }
}
singlePlot <- function (df, rug, plot.start, plot.end, mirror.y, y.plot,title,yLabel) {
  plot(rug, type = "h", xlab = "", ylab = "", xaxt = "n", yaxt = "n", cex = 2, xlim = c(plot.start, plot.end), ylim = c(0,1.6), frame.plot = F)
  if (missing(title)) { title <- "" }
  else {
    mtext(text = title, NORTH <- 3, line = 2)
  }
  plot(x = df[, 2], y = df[, 3],  xlab = "", ylab = yLabel,pch = 16, xaxt = "n", ylim = c(mirror.y, y.plot), xlim = c(plot.start,plot.end), cex = 1.2, col = "#ff7400")
}
trioPlot <- function (df, rug, plot.start, plot.end, mirror.y, y.plot,title,yLabel) {
  plot(rug, type = "h", xlab = "", ylab = "", xaxt = "n", yaxt = "n", cex = 2, xlim = c(plot.start, plot.end), ylim = c(0,1.6), frame.plot = F)
  if (missing(title)) { title <- "" }
  else {
    mtext(text = title, NORTH <- 3, line = 2)
  }
  plot(x = df[, 2], y = df[, 5], xlab = "", ylab = yLabel,pch = 16, xaxt = "n", ylim = c(mirror.y, y.plot), xlim = c(plot.start,plot.end), cex = 1.2, col = "#009999")
  points(x = df[, 2], y = df[, 4], pch = 16, cex = 1.2,col = "#00cc00")
  points(x = df[, 2], y = df[, 3], pch = 16, cex = 1.2,col = "#ff7400")
  legend("topleft", legend = c("Proband", "Mother", "Father"),col = c("#ff7400", "#00cc00", "#009999"), pch = 16, cex = 1,bty = "n", pt.cex = 1.2)
}
genePlotter <- function (genes,start,end,xpos,offset) {

  for(i in 1:nrow(genes)){

    row<- genes[i,]
    plotSet<-row$offset
    tx<-seq(row$txStart,row$txEnd)
    temp<-seq(row$txStart,row$txEnd)
    tx<-as.data.frame(cbind(tx,temp))
    tx$temp<-plotSet
    tx$tx<-tx$tx/offset
    ########

    ########
    start.exons<-strsplit(as.character(row$exonStarts),",")
    end.exons<-strsplit(as.character(row$exonEnds),",")

    start.exons<-as.numeric(unlist(start.exons))
    end.exons<-as.numeric(unlist(end.exons))

    start.exons<-start.exons[order(start.exons)]
    end.exons<-end.exons[order(end.exons)]

    exons<-NULL

    for (i in 1:length(end.exons)) {
      temp<-c(seq(start.exons[i],end.exons[i]))
      exons<-c(exons,temp)
    }

    temp<-seq(1,length(exons))
    exons<-as.data.frame(cbind(exons,temp))
    exons$temp<-plotSet

    #######

    cds.exons<-subset(exons, exons >= row$cdsStart & exons <= row$cdsEnd)
    utr.exons<-subset(exons, exons <= row$cdsStart | exons >= row$cdsEnd)
    cds.exons$exons<-cds.exons$exons/offset
    utr.exons$exons<-utr.exons$exons/offset

    #######

    gene.label<-NULL
    geneTitlePos<-NULL
    genomic.space<-seq(start,end)
    genomic.space<-genomic.space/offset

    if(row$txStart/offset < median(genomic.space)) {
      geneTitlePos<-row$txEnd/offset
    } else {
      geneTitlePos<-row$txStart/offset
    }


    if(row$txStart<start & row$txStart<end & row$txEnd>start & row$txEnd>end){
      geneTitlePos<-median(xpos)
    }

    points(tx,type="l",lwd=1.33)
    points(utr.exons,pch="|",cex=1)
    points(cds.exons,pch="|",col="#000072",cex=1.2)
    if (row$strand == "-") { gene.label <- paste("<- ",row$gene)}
    if (row$strand == "+") { gene.label <- paste(row$gene," ->")}

    textSize<-0.87
    geneOff<-6
    if(nrow(genes)>=6){
      textSize<-0.66
      geneOff<-3
    }
    text(x=geneTitlePos,y=(plotSet-geneOff),labels=gene.label,font=4,cex=textSize)

  }
}
