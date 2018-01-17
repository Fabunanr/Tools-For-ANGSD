library(shiny)
library(genomeIntervals)
library(lattice)
library(Hmisc)
library(ape)
library(data.table)
options(shiny.maxRequestSize = -1)
thetas.headers <- c("(indexStart,indexStop)(firstPos_withData,lastPos_withData)(WinStart,WinStop)","Chr","WinCenter","tW","tP","tF","tH","tL","Tajima","fuf","fud","fayh","zeng","nSites")

not.loaded <- TRUE

# Define server logic required to draw a histogram
shinyServer(


  function(input, output) {

    dataInputThetas = reactive({
      data <- input$userThetas
      path <- as.character(data$datapath)
      thetas <- fread(input=path, sep="\t")
      setnames(thetas,thetas.headers)
      return(thetas)
    })

    dataInputThetasH = reactive({
      data <- input$userThetasH
      path <- as.character(data$datapath)
      thetasH <- fread(input=path, sep="\t")
      setnames(thetasH,thetas.headers)
      return(thetasH)
    })
    
    dataInputThetasL = reactive({
      data <- input$userThetasL
      path <- as.character(data$datapath)
      thetasL <- fread(input=path, sep ="\t")
      setnames(thetasL,thetas.headers)
      return(thetasL)
    })

    gffInput = reactive({
      data <- input$userAnnotations
      path <- as.character(data$datapath)
      gff <- readGff3(path)

      return(gff)

    })

    output$thetaChroms = renderUI({
      if(is.null(input$userThetas)){
        choices <- 10
      }
      else{
      thetas <- dataInputThetas()
      choices <- unique(thetas$Chr)
      }
      selectInput('thetaChrom', 'Chromosome to plot', choices)
    })

    output$thetaChromsH = renderUI({
      if(is.null(input$userThetasH)){
        choices <- 10
      }
      else{
      thetasH <- dataInputThetasH()
      choices <- unique(thetasH$Chr)
      }
      selectInput('thetaChromH', 'Chromosome to plot', choices)
    })
    
    output$thetaChromsL = renderUI({
      if(is.null(input$userThetasL)){
        choices <- 10
      }
      else{
      thetasL <- dataInputThetasL()
      choices <- unique(thetasL$Chr)
      }
      selectInput('thetaChromL', 'Chromosome to plot', choices)
    })
	
    output$thetaPlot <- renderPlot({
#       error handling code to provide a default dataset to graph
      thetas <- tryCatch({
        dataInputThetas()
      }, error = function(err) {
        thetas <- read.table(file=data$datapath,
                             sep="\t",
                             col.names=thetas.headers)
      }
      )

      thetas <- subset(thetas,Chr==input$thetaChrom)
      if(input$annotations){
        validate(need(input$userAnnotations, 'Need GFF file before clicking checkbox!'))
        gff <- gffInput()
        gff.gene <- subset(gff, type="gene")
        gff.df <- data.frame(gff.gene,annotation(gff))
        gff.df.chr <- subset(gff.df, seq_name==thetas$Chr[1])
        if(length(gff.df.chr$seq_name)==0){
          stop("Annotation does not match graphed region. Please make sure the first column of your GFF file matches the Chr column of the .pestPG file.")
        }
        gff.df.gene <- subset(gff.df.chr, type=="gene")
      }

      if(input$subset) {
        thetas.plot <- subset(thetas, WinCenter > input$WinCenterLow & WinCenter < input$WinCenterHigh)
      }
      else {
        thetas.plot <- thetas
      }

      # remove nsites=0
      thetas.plot <- subset(thetas.plot, nSites != 0)
      # remove data points with less than 50 sites. Calculate minimum from data?
      if(input$rm.nsites) {
        thetas.plot <- subset(thetas.plot, nSites > input$nsites)
      }
      #divide thetas by the number of sites in each window
      thetas.plot$tW <- thetas.plot$tW/thetas.plot$nSites
      thetas.plot$tP <- thetas.plot$tP/thetas.plot$nSites
      thetas.plot$tF <- thetas.plot$tF/thetas.plot$nSites
      thetas.plot$tH <- thetas.plot$tH/thetas.plot$nSites
      thetas.plot$tL <- thetas.plot$tW/thetas.plot$nSites

      data <- switch(input$thetaChoice,
                     "Watterson's Theta" = thetas.plot$tW,
                     "Pairwise Theta" = thetas.plot$tP,
                     "Fu and Li's Theta" = thetas.plot$tF,
                     "Fay's Theta" = thetas.plot$tH,
                     "Maximum likelihood (L) Theta" = thetas.plot$tL
      )
      if(input$annotations) {
        plot(thetas.plot$WinCenter,
             data, t="p", pch=19,col=rgb(0,0,0,0.5),
             xlab="Position (bp)",
             ylab=paste(input$thetaChoice,"Estimator Value"),
             main=paste("HIGH Pop Estimators of theta along chromosome", thetas$Chr[1])
        )

        rug(rect(gff.df.gene$X1, -1e2, gff.df.gene$X2, 0, col=rgb(0.18,0.55,0.8,0.75), border=NA))
        if(input$thetaLowess){lines(lowess(thetas.plot$WinCenter,data, f=0.1), col="red")}
      }
      else {
        plot(thetas.plot$WinCenter,
             data, t="p", pch=19,col=rgb(0,0,0,0.5),
             xlab="Position (bp)",
             ylab=paste(input$thetaChoice,"Estimator Value"),
             main=paste("HIGH Pop Estimators of theta along chromosome", thetas$Chr[1])
        )
        if(input$thetaLowess){lines(lowess(thetas.plot$WinCenter,data, f=0.1), col="red")}
      }
    })

    output$thetaPlotH <- renderPlot({
#       error handling code to provide a default dataset to graph
      thetasH <- tryCatch({
        dataInputThetasH()
      }, error = function(err) {
        thetasH <- read.table(file=data$datapath,
                             sep="\t",
                             col.names=thetas.headers)
      }
      )

      thetasH <- subset(thetasH,Chr==input$thetaChromH)
      if(input$annotations){
        validate(need(input$userAnnotations, 'Need GFF file before clicking checkbox!'))
        gff <- gffInput()
        gff.gene <- subset(gff, type="gene")
        gff.df <- data.frame(gff.gene,annotation(gff))
        gff.df.chr <- subset(gff.df, seq_name==thetasH$Chr[1])
        if(length(gff.df.chr$seq_name)==0){
          stop("Annotation does not match graphed region. Please make sure the first column of your GFF file matches the Chr column of the .pestPG file.")
        }
        gff.df.gene <- subset(gff.df.chr, type=="gene")
      }

      if(input$subset) {
        thetasH.plot <- subset(thetasH, WinCenter > input$WinCenterLow & WinCenter < input$WinCenterHigh)
      }
      else {
        thetasH.plot <- thetasH
      }

      # remove nsites=0
      thetasH.plot <- subset(thetasH.plot, nSites != 0)
      # remove data points with less than 50 sites. Calculate minimum from data?
      if(input$rm.nsites) {
        thetasH.plot <- subset(thetasH.plot, nSites > input$nsites)
      }
      #divide thetas by the number of sites in each window
      thetasH.plot$tW <- thetasH.plot$tW/thetasH.plot$nSites
      thetasH.plot$tP <- thetasH.plot$tP/thetasH.plot$nSites
      thetasH.plot$tF <- thetasH.plot$tF/thetasH.plot$nSites
      thetasH.plot$tH <- thetasH.plot$tH/thetasH.plot$nSites
      thetasH.plot$tL <- thetasH.plot$tW/thetasH.plot$nSites

      dataH <- switch(input$thetaChoice,
                     "Watterson's Theta" = thetasH.plot$tW,
                     "Pairwise Theta" = thetasH.plot$tP,
                     "Fu and Li's Theta" = thetasH.plot$tF,
                     "Fay's Theta" = thetasH.plot$tH,
                     "Maximum likelihood (L) Theta" = thetasH.plot$tL
      )
      if(input$annotations) {
        plot(thetasH.plot$WinCenter,
             dataH, t="p", pch=19,col=rgb(0,0,0,0.5),
             xlab="Position (bp)",
             ylab=paste(input$thetaChoice,"Estimator Value"),
             main=paste("MID Pop Estimators of theta along chromosome", thetasH$Chr[1])
        )

        rug(rect(gff.df.gene$X1, -1e2, gff.df.gene$X2, 0, col=rgb(0.18,0.55,0.8,0.75), border=NA))
        if(input$thetaLowess){lines(lowess(thetasH.plot$WinCenter,dataH, f=0.1), col="red")}
      }
      else {
        plot(thetasH.plot$WinCenter,
             dataH, t="p", pch=19,col=rgb(0,0,0,0.5),
             xlab="Position (bp)",
             ylab=paste(input$thetaChoice,"Estimator Value"),
             main=paste("MID Pop Estimators of theta along chromosome", thetasH$Chr[1])
        )
        if(input$thetaLowess){lines(lowess(thetasH.plot$WinCenter,dataH, f=0.1), col="red")}
      }
    })

    output$thetaPlotL <- renderPlot({
#       error handling code to provide a default dataset to graph
      thetasL <- tryCatch({
        dataInputThetasL()
      }, error = function(err) {
        thetasL <- read.table(file=input$thetasL,
                             sep="\t",
                             col.names=thetas.headers)
      }
      )

      thetasL <- subset(thetasL,Chr==input$thetaChromL)
      if(input$annotations){
        validate(need(input$userAnnotations, 'Need GFF file before clicking checkbox!'))
        gff <- gffInput()
        gff.gene <- subset(gff, type="gene")
        gff.df <- data.frame(gff.gene,annotation(gff))
        gff.df.chr <- subset(gff.df, seq_name==thetasL$Chr[1])
        if(length(gff.df.chr$seq_name)==0){
          stop("Annotation does not match graphed region. Please make sure the first column of your GFF file matches the Chr column of the .pestPG file.")
        }
        gff.df.gene <- subset(gff.df.chr, type=="gene")
      }

      if(input$subset) {
        thetasL.plot <- subset(thetasL, WinCenter > input$WinCenterLow & WinCenter < input$WinCenterHigh)
      }
      else {
        thetasL.plot <- thetasL
      }

      # remove nsites=0
      thetasL.plot <- subset(thetasL.plot, nSites != 0)
      # remove data points with less than 50 sites. Calculate minimum from data?
      if(input$rm.nsites) {
        thetasL.plot <- subset(thetasL.plot, nSites > input$nsites)
      }
      #divide thetas by the number of sites in each window
      thetasL.plot$tW <- thetasL.plot$tW/thetasL.plot$nSites
      thetasL.plot$tP <- thetasL.plot$tP/thetasL.plot$nSites
      thetasL.plot$tF <- thetasL.plot$tF/thetasL.plot$nSites
      thetasL.plot$tH <- thetasL.plot$tH/thetasL.plot$nSites
      thetasL.plot$tL <- thetasL.plot$tW/thetasL.plot$nSites

      dataL <- switch(input$thetaChoice,
                     "Watterson's Theta" = thetasL.plot$tW,
                     "Pairwise Theta" = thetasL.plot$tP,
                     "Fu and Li's Theta" = thetasL.plot$tF,
                     "Fay's Theta" = thetasL.plot$tH,
                     "Maximum likelihood (L) Theta" = thetasL.plot$tL
      )
      if(input$annotations) {
        plot(thetasL.plot$WinCenter,
             dataL, t="p", pch=19,col=rgb(0,0,0,0.5),
             xlab="Position (bp)",
             ylab=paste(input$thetaChoice,"Estimator Value"),
             main=paste("LOW Pop Estimators of theta along chromosome", thetasL$Chr[1])
        )

        rug(rect(gff.df.gene$X1, -1e2, gff.df.gene$X2, 0, col=rgb(0.18,0.55,0.8,0.75), border=NA))
        if(input$thetaLowess){lines(lowess(thetasL.plot$WinCenter,dataL, f=0.1), col="red")}
      }
      else {
        plot(thetasL.plot$WinCenter,
             dataL, t="p", pch=19,col=rgb(0,0,0,0.5),
             xlab="Position (bp)",
             ylab=paste(input$thetaChoice,"Estimator Value"),
             main=paste("LOW Pop Estimators of theta along chromosome", thetasL$Chr[1])
        )
        if(input$thetaLowess){lines(lowess(thetasL.plot$WinCenter,dataL, f=0.1), col="red")}
      }
    })
    
    output$selectionPlot <- renderPlot({
      # error handling code to provide a default dataset to graph
      thetas <- tryCatch({
        dataInputThetas()
      }, error = function(err) {
        thetas <- read.table(file=data$datapath,
                             sep="\t",
                             col.names=thetas.headers
        )
      }
      )


      thetas <- subset(thetas,Chr==input$thetaChrom)
      if(input$annotations){
        validate(need(input$userAnnotations, 'Need GFF file before clicking checkbox!'))
        gff <- gffInput()
        gff.gene <- subset(gff, type="gene")
        gff.df <- data.frame(gff.gene,annotation(gff))
        gff.df.chr <- subset(gff.df, seq_name==thetas$Chr[1])
        if(length(gff.df.chr$seq_name)==0){
          stop("Annotation does not match graphed region. Please make sure the first column of your GFF file matches the Chr column of the .pestPG file.")
        }
        gff.df.gene <- subset(gff.df.chr, type=="gene")
      }
      if(input$subset) {
        thetas.plot <- subset(thetas, WinCenter > input$WinCenterLow & WinCenter < input$WinCenterHigh)
      }
      else {
        thetas.plot <- thetas
      }

      # remove nsites=0
      thetas.plot <- subset(thetas.plot, nSites != 0)

      data <- switch(input$selectionChoice,
                     "Tajima's D" = thetas.plot$Tajima,
                     "Fu and Li's F" = thetas.plot$fuf,
                     "Fu and Li's D" = thetas.plot$fud,
                     "Fay and Wu's H" = thetas.plot$fayh,
                     "Zeng's E" = thetas.plot$zeng
      )
      if(input$annotations) {
        plot(thetas.plot$WinCenter,
             data, t="p", pch=19,col=rgb(0,0,0,0.5),
             xlab="Position (bp)",
             ylab=paste(input$selectionChoice,"Estimator Value"),
             main=paste("HIGH Pop Neutrality test statistics along chromosome", thetas$Chr[1])
        )

        rug(rect(gff.df.gene$X1, -1e2, gff.df.gene$X2, 0, col=rgb(0.18,0.55,0.8,0.75), border=NA))
        if(input$thetaLowess){lines(lowess(thetas.plot$WinCenter,data, f=0.1), col="red")}
      }
      else{      
      plot(thetas.plot$WinCenter,
           data, t="p", pch=19,col=rgb(0,0,0,0.5),
           xlab="Position (bp)",
           ylab=input$selectionChoice,
           main=paste("HIGH Pop Neutrality test statistics along chromosome", thetas$Chr[1])
      )
      if(input$selectionLowess){lines(lowess(thetas.plot$WinCenter,data, f=0.1), col="red")}
  }
      })

    output$selectionPlotH <- renderPlot({
      # error handling code to provide a default dataset to graph
      thetasH <- tryCatch({
        dataInputThetasH()
      }, error = function(err) {
        thetasH <- read.table(file=data$datapath,
                             sep="\t",
                             col.names=thetas.headers
        )
      }
      )


      thetasH <- subset(thetasH,Chr==input$thetaChromH)

      if(input$subset) {
        thetasH.plot <- subset(thetasH, WinCenter > input$WinCenterLow & WinCenter < input$WinCenterHigh)
      }
      else {
        thetasH.plot <- thetasH
      }

      # remove nsites=0
      thetasH.plot <- subset(thetasH.plot, nSites != 0)
      thetasH <- subset(thetasH,Chr==input$thetaChromH)
      if(input$annotations){
        validate(need(input$userAnnotations, 'Need GFF file before clicking checkbox!'))
        gff <- gffInput()
        gff.gene <- subset(gff, type="gene")
        gff.df <- data.frame(gff.gene,annotation(gff))
        gff.df.chr <- subset(gff.df, seq_name==thetasH$Chr[1])
        if(length(gff.df.chr$seq_name)==0){
          stop("Annotation does not match graphed region. Please make sure the first column of your GFF file matches the Chr column of the .pestPG file.")
        }
        gff.df.gene <- subset(gff.df.chr, type=="gene")
      }
      dataH <- switch(input$selectionChoice,
                     "Tajima's D" = thetasH.plot$Tajima,
                     "Fu and Li's F" = thetasH.plot$fuf,
                     "Fu and Li's D" = thetasH.plot$fud,
                     "Fay and Wu's H" = thetasH.plot$fayh,
                     "Zeng's E" = thetasH.plot$zeng
      )
      if(input$annotations) {
        plot(thetasH.plot$WinCenter,
             dataH, t="p", pch=19,col=rgb(0,0,0,0.5),
             xlab="Position (bp)",
             ylab=paste(input$selectionChoice,"Estimator Value"),
             main=paste("MID Pop Neutrality test statistics along chromosome", thetasH$Chr[1])
        )

        rug(rect(gff.df.gene$X1, -1e2, gff.df.gene$X2, 0, col=rgb(0.18,0.55,0.8,0.75), border=NA))
        if(input$thetaLowess){lines(lowess(thetasH.plot$WinCenter,dataH, f=0.1), col="red")}
      }
      else{      
      plot(thetasH.plot$WinCenter,
           dataH, t="p", pch=19,col=rgb(0,0,0,0.5),
           xlab="Position (bp)",
           ylab=input$selectionChoice,
           main=paste("MID Pop Neutrality test statistics along chromosome", thetasH$Chr[1])
      )
      if(input$selectionLowess){lines(lowess(thetasH.plot$WinCenter,dataH, f=0.1), col="red")}
  }
    })
    
    output$selectionPlotL <- renderPlot({
      # error handling code to provide a default dataset to graph
      thetasL <- tryCatch({
        dataInputThetasL()
      }, error = function(err) {
        thetasH <- read.table(file=data$datapath,
                             sep="\t",
                             col.names=thetas.headers
        )
      }
      )


      thetasL <- subset(thetasL,Chr==input$thetaChromL)

      if(input$subset) {
        thetasL.plot <- subset(thetasL, WinCenter > input$WinCenterLow & WinCenter < input$WinCenterHigh)
      }
      else {
        thetasL.plot <- thetasL
      }

      # remove nsites=0
      thetasL.plot <- subset(thetasL.plot, nSites != 0)
      if(input$annotations){
        validate(need(input$userAnnotations, 'Need GFF file before clicking checkbox!'))
        gff <- gffInput()
        gff.gene <- subset(gff, type="gene")
        gff.df <- data.frame(gff.gene,annotation(gff))
        gff.df.chr <- subset(gff.df, seq_name==thetasL$Chr[1])
        if(length(gff.df.chr$seq_name)==0){
          stop("Annotation does not match graphed region. Please make sure the first column of your GFF file matches the Chr column of the .pestPG file.")
        }
        gff.df.gene <- subset(gff.df.chr, type=="gene")
      }
      dataL <- switch(input$selectionChoice,
                     "Tajima's D" = thetasL.plot$Tajima,
                     "Fu and Li's F" = thetasL.plot$fuf,
                     "Fu and Li's D" = thetasL.plot$fud,
                     "Fay and Wu's H" = thetasL.plot$fayh,
                     "Zeng's E" = thetasL.plot$zeng
      )
      if(input$annotations) {
        plot(thetasL.plot$WinCenter,
             dataL, t="p", pch=19,col=rgb(0,0,0,0.5),
             xlab="Position (bp)",
             ylab=paste(input$selectionChoice,"Estimator Value"),
             main=paste("LOW Pop Neutrality test statistics along chromosome", thetasL$Chr[1])
        )

        rug(rect(gff.df.gene$X1, -1e2, gff.df.gene$X2, 0, col=rgb(0.18,0.55,0.8,0.75), border=NA))
        if(input$thetaLowess){lines(lowess(thetasL.plot$WinCenter,dataL, f=0.1), col="red")}
      }
      else{      
      plot(thetasL.plot$WinCenter,
           dataL, t="p", pch=19,col=rgb(0,0,0,0.5),
           xlab="Position (bp)",
           ylab=input$selectionChoice,
           main=paste("LOW Pop Neutrality test statistics along chromosome", thetasL$Chr[1])
      )
      if(input$selectionLowess){lines(lowess(thetasL.plot$WinCenter,dataL, f=0.1), col="red")}
  }
    })
})
