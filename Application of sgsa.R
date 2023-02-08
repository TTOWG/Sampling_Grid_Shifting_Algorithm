## To the Only Wise God - TTOWG

library(gstat)
library(sp)
library(dbscan)
library(tidyverse)

#Non-ergodic consideration

source("function_sgsa.R")

sampledata = read.csv(file.choose(),header=T)
zonalanisomodel = vgm(psill = 0.0005, "Sph", range = 1000000000, anis = c(0,90,0,3.5e-6,1.5e-6))
Integrated3Dmodel = vgm(nugget = 0.0024, psill = 0.0016, range = 3500, model = "Sph", anis=c(90,0,0,0.428571428,4.28571428e-3), add.to = zonalanisomodel)
#variogram_model = vgm(nugget = 0.001, psill = 0.0034, range = 3000, model = "Sph", anis=c(90,0,0,0.5,0.1))
set.seed(123)
ApplicationShifts = sgsa(sampledata =  sampledata, x_origin =  700000, y_origin = 732500, z_origin = 0, nx = 40, ny = 13, nz = 100, deltaX = 400, deltaY =  400, deltaZ = 1, vargmodel =Integrated3Dmodel, nmin =  20, nmax = 50)

TotalPossibleshifts = 78
AllShiftingsmeanporoHolder= data.frame(matrix(0,TotalPossibleshifts,2))
names(AllShiftingsmeanporoHolder) = c("Shifting", "MeanPoro")
AllShiftingsVariogramHolder = data.frame(matrix(0,6,3*TotalPossibleshifts))

# Accessing the sample files (.csv) and computing the statistics of intrest
    # Mean
for (mo in 1:TotalPossibleshifts){
  ShiftSamplePointsandData = read.csv(paste0("Sample_",mo,".csv"),header=T)
  # Mean computations
  Shift_meanporo = mean(ShiftSamplePointsandData$attribute)
  
  #Send the output to the overall holder
  AllShiftingsmeanporoHolder[mo,1]=mo
  AllShiftingsmeanporoHolder[mo,2]=Shift_meanporo
  
  # empirical variogram computation
  coordinates(ShiftSamplePointsandData) = ~x_coord+y_coord+z_coord
  ShiftVariogramCloud = variogram(attribute~1,ShiftSamplePointsandData, alpha = 90, beta = 0, tol.hor = 22.5, cutoff = 8000, cloud = TRUE)
  Shiftdbscan = dbscan(ShiftVariogramCloud[,c(2,3)], eps = 120, minPts = 600)
  ShiftdbscanVariogramCloud = data.frame(ShiftVariogramCloud$dist, ShiftVariogramCloud$gamma, Shiftdbscan$cluster, Direction = 90)    
  ShiftdbscanVariogramCloud$Shiftdbscan.cluster = as.factor(ShiftdbscanVariogramCloud$Shiftdbscan.cluster) #Ensuring R sees the cluster numbers as classification ID and not as numerals
  
  NumberofClusters = length(levels(ShiftdbscanVariogramCloud$Shiftdbscan.cluster))
  ShiftdbscanEmpVariog = data.frame(matrix(0, NumberofClusters-1, 3, dimnames = list(c("1","2","3","4","5","6"),c("np", "dist", "gamma"))))
  for(mm in 1:NumberofClusters-1){
    ShiftdbscanEmpVariog[mm,1]=length(ShiftdbscanVariogramCloud$ShiftVariogramCloud.dist[ShiftdbscanVariogramCloud$Shiftdbscan.cluster==mm]); 
    ShiftdbscanEmpVariog[mm,2]=mean(ShiftdbscanVariogramCloud$ShiftVariogramCloud.dist[ShiftdbscanVariogramCloud$Shiftdbscan.cluster==mm]);
    ShiftdbscanEmpVariog[mm,3]=mean(ShiftdbscanVariogramCloud$ShiftVariogramCloud.gamma[ShiftdbscanVariogramCloud$Shiftdbscan.cluster==mm])
  }
  
  #Send the output to the overall holder
  
  AllShiftingsVariogramHolder[((3*mo)-2):(3*mo)] = data.frame(ShiftdbscanEmpVariog$np,ShiftdbscanEmpVariog$dist,ShiftdbscanEmpVariog$gamma)
  names(AllShiftingsVariogramHolder)[((3*mo)-2):(3*mo)] = c(paste("np",mo,sep = ""),paste("dist",mo,sep = ""),paste("gamma",mo,sep = ""))
  
  
}

    # uncertainty (Histogram) of mean
ggplot(AllShiftingsmeanporoHolder, aes(x= MeanPoro ))+
  geom_histogram(bins = 15,color = "blue", fill = "lightgreen")+
  labs(x = "Mean Porosity", y = "Frequency")
    

    #Variogram Uncertainty (mean and Variance): computations and plots
MeanandVarianceofClusterVariogram = data.frame(matrix(0,nrow(AllShiftingsVariogramHolder),3))
names(MeanandVarianceofClusterVariogram) = c("dist", "meanofgamma", "varianceofgamma")
for(mu in 1:nrow(AllShiftingsVariogramHolder)){
  MeanandVarianceofClusterVariogram[mu,1] = AllShiftingsVariogramHolder[mu,2] 
  MeanandVarianceofClusterVariogram[mu,2] = mean(t(AllShiftingsVariogramHolder[mu,seq(from=3, to = 3*TotalPossibleshifts, by = 3)])) 
  MeanandVarianceofClusterVariogram[mu,3] = var(t(AllShiftingsVariogramHolder[mu,seq(from=3, to = 3*TotalPossibleshifts, by = 3)]))
}    

VarModelValues = variogramLine(Integrated3Dmodel, 8000, 80, dir=c(1,0,0), dist_vector = seq(1,8000,1))
wantedplot = ggplot(data = AllShiftingsVariogramHolder)+
  geom_point(data = MeanandVarianceofClusterVariogram, aes(x = dist, y = meanofgamma, shape = "dc", size = "bi"),color = "blue")+
  geom_line(data = VarModelValues,aes(x = dist, y = gamma, color = "cu"), size = 0.75)+
  geom_point(aes(x = dist1, y = gamma1, shape = "pe", size = "sm"),color = "red")+
  labs(x = "Lag Distance, m", y = "Porosity Variogram")+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 0.006))+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 0.006, by = 0.001))+
  scale_colour_manual(name = '', values =c("cu" = "black"), labels = c("cu" = "Model"))+
  scale_shape_manual(name = '', values = c('dc' = 8,'pe' = 16), labels = c('dc' = 'Mean','pe' = 'Estimates' ))+
  scale_size_manual(name = '', values = c('bi' = 4,'sm' = 1), labels = c('bi' = 'Mean','sm' = 'Estimates' ))+
  theme(legend.position=c(0.6, 0.3))+
  theme(legend.title = element_blank())+
  theme(legend.background = element_rect(fill="lightgrey", size=0.5, linetype="solid"))
for(ss in 1:TotalPossibleshifts){
  wantedplot = wantedplot+geom_point(aes_string(x = AllShiftingsVariogramHolder[,((3*ss)-1)], y = AllShiftingsVariogramHolder[,3*ss]), color = "red")
}
wantedplot
ggsave("VargFamilyPlot.jpg", dpi = 96)

  

#Ergodic considerations
  #Simulations
simpar = gstat(formula = attribute~1, locations = ~x_coord+y_coord+z_coord, data = sampledata, model = Integrated3Dmodel, nmin = 20, nmax = 50)
SamplingPointsCoord = read.csv("Sample Block Centers.csv",header=T)
set.seed(123)
ErgodicRealizations = predict(simpar, newdata = SamplingPointsCoord, nsim = TotalPossibleshifts, debug=-1)

AllRealizationsmeanporoHolder= data.frame(matrix(0,TotalPossibleshifts,2))
names(AllRealizationsmeanporoHolder) = c("Realization", "MeanPoro")
AllRealizationsVariogramHolder = data.frame(matrix(0,6,3*TotalPossibleshifts))

# Accessing the columns of simulation output files and computing the statistics of intrest
for (yo in 1:TotalPossibleshifts){
  RealizationSamplePointsandData = ErgodicRealizations[,c(1:3,yo+3)]
  names(RealizationSamplePointsandData)[4] = c("attribute")
  # Mean computations
  Realization_meanporo = mean(RealizationSamplePointsandData$attribute)
  
  #Send the output to the overall holder
  AllRealizationsmeanporoHolder[yo,1]=yo
  AllRealizationsmeanporoHolder[yo,2]=Realization_meanporo
  
  # empirical variogram computation
  coordinates(RealizationSamplePointsandData) = ~x_coord+y_coord+z_coord
  RealizationVariogramCloud = variogram(attribute~1,RealizationSamplePointsandData, alpha = 90, beta = 0, tol.hor = 22.5, cutoff = 8000, cloud = TRUE)
  Realizationdbscan = dbscan(RealizationVariogramCloud[,c(2,3)], eps = 120, minPts = 600)
  RealizationdbscanVariogramCloud = data.frame(RealizationVariogramCloud$dist, RealizationVariogramCloud$gamma, Realizationdbscan$cluster, Direction = 90)    
  RealizationdbscanVariogramCloud$Realizationdbscan.cluster = as.factor(RealizationdbscanVariogramCloud$Realizationdbscan.cluster) #Ensuring R sees the cluster numbers as classification ID and not as numerals
  
  NumberofClusters = length(levels(RealizationdbscanVariogramCloud$Realizationdbscan.cluster))
  RealizationdbscanEmpVariog = data.frame(matrix(0, NumberofClusters-1, 3, dimnames = list(c("1","2","3","4","5","6"),c("np", "dist", "gamma"))))
  for(rr in 1:NumberofClusters-1){
    RealizationdbscanEmpVariog[rr,1]=length(RealizationdbscanVariogramCloud$RealizationVariogramCloud.dist[RealizationdbscanVariogramCloud$Realizationdbscan.cluster==rr]); 
    RealizationdbscanEmpVariog[rr,2]=mean(RealizationdbscanVariogramCloud$RealizationVariogramCloud.dist[RealizationdbscanVariogramCloud$Realizationdbscan.cluster==rr]);
    RealizationdbscanEmpVariog[rr,3]=mean(RealizationdbscanVariogramCloud$RealizationVariogramCloud.gamma[RealizationdbscanVariogramCloud$Realizationdbscan.cluster==rr])
  }
  
  #Send the output to the overall holder
  
  AllRealizationsVariogramHolder[((3*yo)-2):(3*yo)] = data.frame(RealizationdbscanEmpVariog$np,RealizationdbscanEmpVariog$dist,RealizationdbscanEmpVariog$gamma)
  names(AllRealizationsVariogramHolder)[((3*yo)-2):(3*yo)] = c(paste("np",yo,sep = ""),paste("dist",yo,sep = ""),paste("gamma",yo,sep = ""))  
}

# uncertainty (Histogram) of mean
ggplot(AllRealizationsmeanporoHolder, aes(x= MeanPoro ))+
  geom_histogram(bins = 15,color = "blue", fill = "lightgreen")+
  labs(x = "Mean Porosity", y = "Frequency")

#Variogram Uncertainty (mean and Variance): computations and plots
MeanandVarianceofClusterVariogram_erg = data.frame(matrix(0,nrow(AllRealizationsVariogramHolder),3))
names(MeanandVarianceofClusterVariogram_erg) = c("dist", "meanofgamma", "varianceofgamma")
for(ru in 1:nrow(AllRealizationsVariogramHolder)){
  MeanandVarianceofClusterVariogram_erg[ru,1] = AllRealizationsVariogramHolder[ru,2] 
  MeanandVarianceofClusterVariogram_erg[ru,2] = mean(t(AllRealizationsVariogramHolder[ru,seq(from=3, to = 3*TotalPossibleshifts, by = 3)])) 
  MeanandVarianceofClusterVariogram_erg[ru,3] = var(t(AllRealizationsVariogramHolder[ru,seq(from=3, to = 3*TotalPossibleshifts, by = 3)]))
}    
####################################################################
VarModelValues = variogramLine(Integrated3Dmodel, 8000, 80, dir=c(1,0,0), dist_vector = seq(1,8000,1))
wantedplot = ggplot(data = AllShiftingsVariogramHolder)+
  geom_point(data = MeanandVarianceofClusterVariogram, aes(x = dist, y = meanofgamma, shape = "dc", size = "bi"),color = "blue")+
  geom_line(data = VarModelValues,aes(x = dist, y = gamma, color = "cu"), size = 0.75)+
  geom_point(aes(x = dist1, y = gamma1, shape = "pe", size = "sm"),color = "red")+
  labs(x = "Lag Distance, m", y = "Porosity Variogram")+
  coord_cartesian(xlim = c(0, 8000), ylim = c(0, 0.006))+
  scale_x_continuous(breaks = seq(0, 8000, by = 1000))+
  scale_y_continuous(breaks = seq(0, 0.006, by = 0.001))+
  scale_colour_manual(name = '', values =c("cu" = "black"), labels = c("cu" = "Model"))+
  scale_shape_manual(name = '', values = c('dc' = 8,'pe' = 16), labels = c('dc' = 'Mean','pe' = 'Estimates' ))+
  scale_size_manual(name = '', values = c('bi' = 4,'sm' = 1), labels = c('bi' = 'Mean','sm' = 'Estimates' ))+
  theme(legend.position=c(0.6, 0.3))+
  theme(legend.title = element_blank())+
  theme(legend.background = element_rect(fill="lightgrey", size=0.5, linetype="solid"))
for(ss in 1:TotalPossibleshifts){
  wantedplot = wantedplot+geom_point(aes_string(x = AllShiftingsVariogramHolder[,((3*ss)-1)], y = AllShiftingsVariogramHolder[,3*ss]), color = "red")
}
wantedplot
ggsave("VargFamilyPlot.jpg", dpi = 96)