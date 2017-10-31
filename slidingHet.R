##This script was written in R version 3.2.3


he <- function(data) {
  nLoci <- (length(names(data))-1)/2
  preSum <- 0
  missingLoci <- 0
  for (i in 1:nLoci){
    alleles <- c(data[,i*2],data[,i*2+1])
    noMissing <- alleles[alleles!=-9]
    uniqueAll <- unique(alleles)
    if (length(noMissing)==0){
      missingLoci <- missingLoci +1
    }
    else {
      heLocal <- 0
      for (j in uniqueAll) {
        
        heLocal <- heLocal + (sum(noMissing==j)/(length(noMissing)))^2
      }
    }
    preSum <- preSum + (heLocal)
  }
  He <- 1-(1/(nLoci-missingLoci))*preSum
  return(He)
}  

nDiv <- function(dataSeq) {
  Seq <- dataSeq
  names(Seq) <- c("Sample","Sequence")
  ##Add a exception to the windows that will have the sequenceless sample
  if (length(Seq[Seq$Sequence=="",2]) == 0){
    
    Seq <- Seq[Seq$Sequence!="",]
    
  }
  
  ##Calculates each haplotype frequency
  Seq["Freq"]<-as.numeric(0)
  #to calculate the haplotype frequency, for each sequence we match it against all sequences (including itself) in order to calculate the number of times that haplotype occurs in one group
  for (m in 1:length(Seq$Sequence)){
    ntimes = 0
    
    for (n in 1:length(Seq$Sequence)){
      if (Seq$Sequence[m] == Seq$Sequence[n]) {      
        
        ntimes = ntimes + 1      
        
      }    
      
    }
    #For each sequence the haplotype frequency is then added (it will include duplicates in this phase)
    Seq$Freq[m] <- ntimes/(length(Seq$Sequence))  
  }
  ##Creates a new dataframe with unique haplotypes
  Seqs = data.frame()
  #Unique haplotypes are extracted from the total list by creating a new empty dataset which is only fed with a new sequence if this sequence is different from previously added ones
  for (s in 1:length(Seq$Sequence)) {
    if (Seq$Sequence[s] %in% Seqs$Sequence) {
      
      next
      
    }
    else  {
      
      rbind(Seqs,data.frame(Sample = (Seq$Sample[s]), Sequence = (Seq$Sequence[s]),Freq= (Seq$Freq[s]))) -> Seqs
      
    }	
    
  }    
  
  Seqs$Freq = as.numeric(Seqs$Freq)
  Seqs$Sequence = as.character(Seqs$Sequence)
  Seqs$Sample = as.character(Seqs$Sample)
  ################################################
  ## If less then two unique haplotypes were sampled in a group, its nucleotide diversity was 0
  if (length (Seqs$Sequence)< 2){
    return (0)
  }
  ## If more than two unique haplotypes were found we calculated the nucleotide diversity by following Nei & Kumar (2000) formula
  
  else if (length(Seqs$Sequence) >= 2) {
    
    
    ##Calculates the ndiversity
    
    n = as.numeric(length(Seq$Sequence))
    
    nDiversity = 0
    
    
    for (i in 1:(length(Seqs$Sequence)-1)){
      
      for (j in i+1:(length(Seqs$Sequence)-i)){
        
        s1=strsplit(Seqs$Sequence[i], "")[[1]]
        s2=strsplit(Seqs$Sequence[j], "")[[1]]
        
        diff = 0
        total = 0    
        
        for (k in 1:length(s1)) {
          
          if (s1[k] == "-" | s2[k]== "-" | s1[k] == "N" | s2[k] == "N") {
            
            next
          } 
          else if(s1[k]!=s2[k]){
            
            diff = diff + 1
            total = total +1
            
          }
          else if(s1[k]==s2[k]){
            
            total = total +1
            
          }
          
        }
        
        nDiversity = nDiversity + ((2*(n/(n-1)))*((Seqs$Freq[i])*(Seqs$Freq[j])*(diff/total)))
        
      }
    }  
    
    
    #Data is then added to a data.frame and saved in the output directory
    return (nDiversity)
  }
  
}

linearCoordinates <<- function(data, transformAxis = "Longitude") {
  #########################
  ####  Linear model  #####
  #########################
  ##  Define the linear model that best represents the data to have the windwow going through it
  tmpData <- data
  
  lModel <- lm(tmpData$LATITUDE~tmpData$LONGITUDE)
  wVal <- coef(lModel)
  ##  Get the line module on the form of slope intercept (y=mx+b)
  intercept <- wVal[[1]]
  slope <- wVal[[2]]
  ##  Transform the latitude into the new linear model predications
  if (tolower(transformAxis) == "longitude") {
      tmpData$transformed<-(tmpData$LATITUDE-intercept)/slope
      return (list(tmpData,wVal))
  } 
  if (tolower(transformAxis) == "latitude") {
    tmpData$transformed <- slope*tmpData$LONGITUDE+intercept
    return (list(tmpData,wVal))
  }
  else{
    return ("Error invalid axis was selected, choices are Latitude or Longitude")
  }
  #######################
  
}

slidingWindow <- function(lineardata,window_width= 0.5, slide_by = 1, seq_included = TRUE) {
  ########################################
  ########## Sliding window ##############
  ########################################
  data <- lineardata[[1]]
  
  slope = lineardata[[2]][[2]]
  intercept = lineardata[[2]][[1]]
  ##  Test variables for easy initializing
  mi <- min(data$transformed)
  ma <- max(data$transformed)
  totalDistnace <- ma - mi
  
  #######################################
  ########## Sliding window #############
  ##########    variables   #############
  #######################################
  
  
  ##  Initializing variables
  startCenter <- mi               ## Where to start the window 
  wWidth <- window_width          ## Width of the window divided by 2 (width from the center to the edge)
  sDistance <- slide_by           ## How much the window should move per iteration
  endCenter <- ma                 ## Last center
  nIts <- (endCenter-startCenter)/sDistance
  
  ## Results storage
  
  fHe <- list ()
  fNd <- list ()
  fCd <- list ()
  fN  <- list ()
  ## For ArcGis Dataframe
  fData <-data[,c("CODE","LATITUDE","LONGITUDE","transformed")]
  
  
  for (i in 1:(nIts+1)) {
    tmpData <- data[data$transformed>=(startCenter-wWidth) & data$transformed <= (startCenter+wWidth),]
    if (dim(tmpData)[1] == 0) {
      fHe <- append(fHe, NA)
      fNd <- append(fNd, NA)
      fCd <- append(fCd, startCenter)
      fN <-  append(fN, 0)
      startCenter <- startCenter + sDistance
      fData[fData[,1] %in% tmpData[,1],i+4] <- i
      next()
    }
    if (seq_included){
      tmpMicros <- tmpData[,c(1,5:dim(tmpData)[2]-2)]
      tmpSeqs <- tmpData[tmpData[,dim(tmpData)[2]-1]!="",c(1,dim(tmpData)[2]-1)]
      fHe <- append(fHe, he(tmpMicros))
      fNd <- append(fNd, nDiv(tmpSeqs))
      fCd <- append(fCd, startCenter)
      fN  <- append(fN, dim(tmpData)[1])
      startCenter <- startCenter + sDistance
    
      ## Build the final dataset
      fData[fData[,1] %in% tmpData[,1],i+4] <- i
    ##
    } else {
      tmpMicros <- tmpData[,c(1,5:dim(tmpData)[2]-1)]
      fHe <- append(fHe, he(tmpMicros))
      fCd <- append(fCd, startCenter)
      fN  <- append(fN, dim(tmpData)[1])
      startCenter <- startCenter + sDistance
      
      ## Build the final dataset
      fData[fData[,1] %in% tmpData[,1],i+4] <- i
      ##  
      
      
      
      
      
    }
  }
  
  ##  Build the dataset per window 
  Windows <- data.frame(Window=c(1:length(fHe)))
  Windows$tLongitude <- unlist(fCd)
  ##  Recover center longitude
  fLg <- list()
  for (i in Windows[,2]){
    fLg <- append(fLg,(slope*i)+intercept)
    
  }
  ##
  Windows$tLatitude <- unlist(fLg)
  Windows$N <- unlist(fN)
  Windows$He <- unlist(fHe)
  if (seq_included) {
    Windows$nDiv <- unlist(fNd)
  }
  
  lB <- NULL
  uB <- NULL
  
  for (i in 1:length(fCd)) {
    lB<-c(lB,(((slope*(fCd[[i]]-wWidth))+intercept)))
    uB<-c(uB,(((slope*(fCd[[i]]+wWidth))+intercept)))
  }
  
  Windows$lowerBound <- lB
  Windows$upperBound <- uB
  
  write.csv(Windows,"perWindowResults.csv",row.names = F)
  write.csv(fData,"sampleGroups.csv",row.names= F)
  return(list(fData,Windows))
}

dataIntegrity  <- function (data,seq_included=TRUE) {
  
  header <- names(data)
  necessaryNames <- c("CODE","SPECIES","LATITUDE","LONGITUDE")
  
  if (all(is.element(necessaryNames,header))) {
    if (seq_included){
      if ((length(header)-5)%%2!=0){
      
        print ("After accounting for CODE,SPECIES,LATITUDE,LONGITUDE and SEQUENCE collumns, number of collumns is not pair check microsats table")
      }else{
      print ("Header seems to be correct")
      print ("Make sure the last collumn of your dataset contains the sequence data")
      }
    }else {
      
      if ((length(header)-4)%%2!=0){
        
        print ("After accounting for CODE,SPECIES,LATITUDE,LONGITUDE and SEQUENCE collumns, number of collumns is not pair check microsats table")
      }else {
      print ("Header seems to be correct")
      }
    }
    
  }else {
    
    print ("ERROR: Data header should contain:")
    print ("a CODE collumn (sample ids)")
    print ("a SPECIES collumn (species name)")
    print ("a LATITUDE collumn (sample Latitude)")
    print ("and LONGITUDE (sample Longitude)")
    print ("REMEMBER: All collumn names should be in all caps")
    
  }
  
  
}
  
dataVisualizer <- function(lineardata, slidingWindowList)  {
  require(rworldmap)
  windows <- slidingWindowList[[2]]
  #####  Map Maker########
  map <- getMap(resolution = "low")
  fData <- slidingWindowList[[1]]
  slope = lineardata[[2]][[2]]
  intercept = lineardata[[2]][[1]]
  ##All Data
  plot (map, xlim= c(min(fData$LONGITUDE)-1,max(fData$LONGITUDE)+1), ylim = c(min(fData$LATITUDE)-1,max(fData$LATITUDE)+1),asp=1)
  points(fData$LONGITUDE,fData$LATITUDE,col = "red", cex = 0.8, pch= 16)
  abline(a= intercept, b = slope)
  
  # 1DWindow Attribution
  pointsMap <- fData[,c(1,2,3,4)]
  for (i in 1:dim(fData)[1]) {
    number <- fData[i,c(5:dim(fData)[2])]
    
    pointsMap$Value[i] <- sum(number, na.rm = T)#/length(number)-sum(is.na(number))
    
  }
  
  colPoints <- rainbow(max(pointsMap$Value))
  plot (map, xlim= c(min(LONGITUDE)-0.5,max(LONGITUDE)+0.5), ylim = c(min(LATITUDE)-1,max(LATITUDE)+1),asp=1)
  points(pointsMap$LONGITUDE, pointsMap$LATITUDE, col=colPoints[pointsMap$Value], pch=16, cex=1)
  
  colWindow <- rainbow(dim(windows)[1])
  
  for (i in 1:dim(windows)[1]) {
    abline(h = windows[i,"lowerBound"],col=colWindow[i])
    abline(h = windows[i,"upperBound"],col=colWindow[i])
    
    
    
    
    
  }
}
