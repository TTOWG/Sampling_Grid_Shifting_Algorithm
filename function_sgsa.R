## To the Only Wise God - TTOWG

sgsa = function(sampledata, x_origin = 0, y_origin = 0, z_origin = 0, nx = 1, ny = 1, nz = 1, deltaX = 0, deltaY = 0, deltaZ = 0, vargmodel, beta = NULL, nmin, nmax){
  
library(gstat)
library(sp)

samplesize = nrow(sampledata)
dimensionality = ncol(sampledata)-1

basicgrid=expand.grid(1:nx,1:ny,1:nz)
names(basicgrid) = c("x","y","z")

if(identical(c(nx,deltaX),c(1,0))){
  x_coord = c()                                 # deactivates unwanted dimension
  } else {
  x_coord = (deltaX*(basicgrid$x-0.5))+x_origin # x-coordinate of the block center
}
if(identical(c(ny,deltaY),c(1,0))){
  y_coord = c()                                 # deactivates unwanted dimension
  } else {
  y_coord = (deltaY*(basicgrid$y-0.5))+y_origin # y-coordinate of the block center
}
if(identical(c(nz,deltaZ),c(1,0))){
  z_coord = c()                                  # deactivates unwanted dimension
  } else {
  z_coord = (deltaZ*(basicgrid$z-0.5))+z_origin # z-coordinate of the block center
}

scaledgrid = data.frame(cbind(x_coord, y_coord, z_coord))

if (identical(c(deltaY,deltaZ),c(0,0))){
  simpar = gstat(formula = attribute~1, locations = ~x_coord, data = sampledata, model = vargmodel, beta = beta, nmin = nmin, nmax = nmax)
  } else {
  if (identical(c(deltaX,deltaZ),c(0,0))){
    simpar = gstat(formula = attribute~1, locations = ~y_coord, data = sampledata, model = vargmodel, beta = beta, nmin = nmin, nmax = nmax)
  } else {
    if (identical(c(deltaX,deltaY),c(0,0))){
      simpar = gstat(formula = attribute~1, locations = ~z_coord, data = sampledata, model = vargmodel, beta = beta, nmin = nmin, nmax = nmax)
    } else {
    if (identical(deltaZ,0)) {
      simpar = gstat(formula = attribute~1, locations = ~x_coord+y_coord, data = sampledata, model = vargmodel, beta = beta, nmin = nmin, nmax = nmax)
    } else {
      if (identical(deltaY,0)) {
        simpar = gstat(formula = attribute~1, locations = ~x_coord+z_coord, data = sampledata, model = vargmodel, beta = beta, nmin = nmin, nmax = nmax)
      } else {
        if (identical(deltaX,0)) {
          simpar = gstat(formula = attribute~1, locations = ~y_coord+z_coord, data = sampledata, model = vargmodel, beta = beta, nmin = nmin, nmax = nmax)
        } else {
          simpar = gstat(formula = attribute~1, locations = ~x_coord+y_coord+z_coord, data = sampledata, model = vargmodel, beta = beta, nmin = nmin, nmax = nmax)
        } 
        }
      }  
    }
  }
}

simulatedfullgrid = predict(simpar, newdata = scaledgrid, nsim = 1,debug=-1)
write.csv(simulatedfullgrid,"simulated_fullgrid.csv", row.names = F)

fullsampledata = data.frame(matrix(0,samplesize,8))
names(fullsampledata) = c("x_index","y_index","z_index","block_natord","x_coord","y_coord", "z_coord", "attribute")
for(i in 1:samplesize){
  fullsampledata$x_index[i] = if (deltaX>0) ceiling((sampledata$x_coord[i] - x_origin)/deltaX) else 1
  fullsampledata$y_index[i] = if (deltaY>0) ceiling((sampledata$y_coord[i] - y_origin)/deltaY) else 1
  fullsampledata$z_index[i] = if (deltaZ>0) ceiling((sampledata$z_coord[i] - z_origin)/deltaZ) else 1
  fullsampledata$block_natord[i] = ((fullsampledata$z_index[i]-1)*nx*ny)+((fullsampledata$y_index[i]-1)*nx)+fullsampledata$x_index[i]
  fullsampledata$x_coord[i] = if (deltaX>0) sampledata$x_coord[i] else 0
  fullsampledata$y_coord[i] = if (deltaY>0) sampledata$y_coord[i] else 0
  fullsampledata$z_coord[i] = if (deltaZ>0) sampledata$z_coord[i] else 0
  fullsampledata$attribute[i] = sampledata$attribute[i]
}


ES = nx-max(fullsampledata$x_index)
NS = ny-max(fullsampledata$y_index)
DS = nz-max(fullsampledata$z_index)
WS = min(fullsampledata$x_index)
SS = min(fullsampledata$y_index)
US = min(fullsampledata$z_index)
S_lat = WS+ES
S_longit = SS+NS
S_planar = S_longit*S_lat
S_total =S_planar*(US+DS)


#Placeholders
Planarshifts = data.frame(matrix(0,samplesize,(7*S_planar)))
LatitudinalShifts = data.frame(matrix(0,samplesize,(7*S_lat)))

# Original sample Location and westward shifts re-samples

# Determine the range of positions for westward shiftings
West_ends = ((US-1)*S_planar)+((SS-1)*S_lat)+WS 
West_begins = West_ends-WS+1
for(w in 1:WS){
  # Determine position for wth shifting
  Position_west = West_ends-(w-1)
  # Generate the wth coordinate suite
  W_th_CoordSuite = data.frame(fullsampledata$x_index-(w-1),fullsampledata$y_index,fullsampledata$z_index,fullsampledata$block_natord-(w-1),fullsampledata$x_coord-((w-1)*deltaX),fullsampledata$y_coord,fullsampledata$z_coord)
  names(W_th_CoordSuite) = c("x_index","y_index","z_index","block_natord","x_coord","y_coord", "z_coord")
  # Assign generated coordinates to placeholders
  Planarshifts[,((7*(SS-1)*S_lat)+(7*(WS-w))+1):((7*(SS-1)*S_lat)+(7*(WS-(w-1))))] = W_th_CoordSuite
  names(Planarshifts)[((7*(SS-1)*S_lat)+(7*(WS-w))+1):((7*(SS-1)*S_lat)+(7*(WS-(w-1))))] = c(paste("x_index",Position_west, sep = ""),paste("y_index",Position_west, sep = ""),paste("z_index",Position_west, sep = ""),paste("block_natord",Position_west, sep = ""),paste("x_coord",Position_west, sep = ""),paste("y_coord",Position_west, sep = ""),paste("z_coord",Position_west, sep = ""))
  LatitudinalShifts[,(7*(WS-w)+1):(7*(WS-(w-1)))] = W_th_CoordSuite 
  # Sample the full-grid realization at generated coordinates
  W_th_data = data.frame(matrix(0,samplesize,1)) 
  for (s in 1:samplesize){
    W_th_data[s,] = simulatedfullgrid[as.numeric(rownames(simulatedfullgrid))==W_th_CoordSuite$block_natord[s],dimensionality+1]
  }
  names(W_th_data) = c("attribute")
  W_th_Sample_v1 = if (deltaX>0) data.frame(W_th_CoordSuite,W_th_data) else data.frame(W_th_CoordSuite[,-c(which(names(W_th_CoordSuite)=="x_index"),which(names(W_th_CoordSuite)=="x_coord"))],W_th_data)
  W_th_Sample_v2 = if (deltaY>0) W_th_Sample_v1 else W_th_Sample_v1[,-c(which(names(W_th_Sample_v1)=="y_index"),which(names(W_th_Sample_v1)=="y_coord"))]
  W_th_Sample_final = if (deltaZ>0) W_th_Sample_v2 else W_th_Sample_v2[,-c(which(names(W_th_Sample_v2)=="z_index"),which(names(W_th_Sample_v2)=="z_coord"))]
  
  # Export sample data
  write.csv(W_th_Sample_final,paste0("Sample_",Position_west,".csv"), row.names = F)
}

# Eastward shifts re-samples

if(ES>0){
  # Determine the range of positions for eastward shiftings
  East_begins = ((US-1)*S_planar)+((SS-1)*S_lat)+WS+1 
  East_ends = East_begins+ES-1
  for(e in 1:ES){
    # Determine position for eth shifting
    Position_east = East_begins+(e-1)
    # Generate the eth coordinate suite
    e_th_CoordSuite = data.frame(fullsampledata$x_index+e,fullsampledata$y_index,fullsampledata$z_index,fullsampledata$block_natord+e,fullsampledata$x_coord+(e*deltaX),fullsampledata$y_coord,fullsampledata$z_coord)
    names(e_th_CoordSuite) = c("x_index","y_index","z_index","block_natord","x_coord","y_coord", "z_coord")
    # Assign generated coordinates to placeholders
    Planarshifts[,((7*(SS-1)*S_lat)+(7*(WS+(e-1)))+1):((7*(SS-1)*S_lat)+(7*(WS+e)))] = e_th_CoordSuite
    names(Planarshifts)[((7*(SS-1)*S_lat)+(7*(WS+(e-1)))+1):((7*(SS-1)*S_lat)+(7*(WS+e)))] = c(paste("x_index",Position_east, sep = ""),paste("y_index",Position_east, sep = ""),paste("z_index",Position_east, sep = ""),paste("block_natord",Position_east, sep = ""),paste("x_coord",Position_east, sep = ""),paste("y_coord",Position_east, sep = ""),paste("z_coord",Position_east, sep = ""))
    LatitudinalShifts[,(7*(WS+e-1)+1):(7*(WS+e))] = e_th_CoordSuite
    # Sample the full-grid realization at generated coordinates
    e_th_data = data.frame(matrix(0,samplesize,1)) 
    for (t in 1:samplesize){
      e_th_data[t,] = simulatedfullgrid[as.numeric(rownames(simulatedfullgrid))==e_th_CoordSuite$block_natord[t],dimensionality+1]
    }
    names(e_th_data) = c("attribute")
    e_th_Sample_v1 = if (deltaX>0) data.frame(e_th_CoordSuite,e_th_data) else data.frame(e_th_CoordSuite[,-c(which(names(e_th_CoordSuite)=="x_index"),which(names(e_th_CoordSuite)=="x_coord"))],e_th_data)
    e_th_Sample_v2 = if (deltaY>0) e_th_Sample_v1 else e_th_Sample_v1[,-c(which(names(e_th_Sample_v1)=="y_index"),which(names(e_th_Sample_v1)=="y_coord"))]
    e_th_Sample_final = if (deltaZ>0) e_th_Sample_v2 else e_th_Sample_v2[,-c(which(names(e_th_Sample_v2)=="z_index"),which(names(e_th_Sample_v2)=="z_coord"))]
   
    
    # Export sample data
    write.csv(e_th_Sample_final,paste0("Sample_",Position_east,".csv"), row.names = F)
    }
}

# Southward shifts of all latitudinal shifts
if(SS>1){
  # Determine the range of positions for southward shiftings
  South_ends = ((US-1)*S_planar)+((SS-1)*S_lat) 
  South_begins = South_ends-((SS-1)*S_lat)+1
  for(s in 1:(SS-1)){
    for(l in 1:S_lat){
      # Determine position for wth shifting
      Position_south = South_ends-((s-1)*S_lat)-(l-1)
      # Extract the lth latitudinal shifting to be shifted southward 
      l_th_latitudinalshift = LatitudinalShifts[,(7*(S_lat-l)+1):(7*(S_lat-(l-1)))] 
      names(l_th_latitudinalshift) = c("x_index","y_index","z_index","block_natord","x_coord","y_coord", "z_coord")
      # Generate the sth coordinate suite
      s_th_CoordSuite = data.frame(l_th_latitudinalshift$x_index,l_th_latitudinalshift$y_index-s,l_th_latitudinalshift$z_index,l_th_latitudinalshift$block_natord-(s*nx),l_th_latitudinalshift$x_coord,l_th_latitudinalshift$y_coord-(s*deltaY),l_th_latitudinalshift$z_coord)
      names(s_th_CoordSuite) = c("x_index","y_index","z_index","block_natord","x_coord","y_coord", "z_coord")
      # Assign generated coordinates to placeholder
      Planarshifts[,((7*(((SS-1-(s-1))*S_lat)-l)+1)):(7*(((SS-1-(s-1))*S_lat)-(l-1)))] = s_th_CoordSuite
      names(Planarshifts)[((7*(((SS-1-(s-1))*S_lat)-l)+1)):(7*(((SS-1-(s-1))*S_lat)-(l-1)))] = c(paste("x_index",Position_south, sep = ""),paste("y_index",Position_south, sep = ""),paste("z_index",Position_south, sep = ""),paste("block_natord",Position_south, sep = ""),paste("x_coord",Position_south, sep = ""),paste("y_coord",Position_south, sep = ""),paste("z_coord",Position_south, sep = ""))
      # Sample the full-grid realization at generated coordinates
      s_th_data = data.frame(matrix(0,samplesize,1)) 
      for (v in 1:samplesize){
        s_th_data[v,] = simulatedfullgrid[as.numeric(rownames(simulatedfullgrid))==s_th_CoordSuite$block_natord[v],dimensionality+1]
      }
      names(s_th_data) = c("attribute")
      s_th_Sample_v1 = if (deltaX>0) data.frame(s_th_CoordSuite,s_th_data) else data.frame(s_th_CoordSuite[,-c(which(names(s_th_CoordSuite)=="x_index"),which(names(s_th_CoordSuite)=="x_coord"))],s_th_data)
      s_th_Sample_v2 = if (deltaY>0) s_th_Sample_v1 else s_th_Sample_v1[,-c(which(names(s_th_Sample_v1)=="y_index"),which(names(s_th_Sample_v1)=="y_coord"))]
      s_th_Sample_final = if (deltaZ>0) s_th_Sample_v2 else s_th_Sample_v2[,-c(which(names(s_th_Sample_v2)=="z_index"),which(names(s_th_Sample_v2)=="z_coord"))]
     
      
      # Export sample data
      write.csv(s_th_Sample_final,paste0("Sample_",Position_south,".csv"), row.names = F)
      }
  }
}


# Northward shifts of all latitudinal shifts
if(NS>0){
  # Determine the range of positions for northward shiftings
  North_begins = ((US-1)*S_planar)+(SS*S_lat)+1 
  North_ends = North_begins+(NS*S_lat)-1
  for(n in 1:NS){
    for(l in 1:S_lat){
      # Determine position for nth shifting
      Position_north = North_begins+((n-1)*S_lat)+(l-1)
      # Extract the lth latitudinal shifting to be shifted northward
      l_th_latitudinalshift = LatitudinalShifts[,(7*(l-1)+1):(7*l)] 
      names(l_th_latitudinalshift) = c("x_index","y_index","z_index","block_natord","x_coord","y_coord", "z_coord")
      # Generate the nth coordinate suite
      n_th_CoordSuite = data.frame(l_th_latitudinalshift$x_index,l_th_latitudinalshift$y_index+n,l_th_latitudinalshift$z_index,l_th_latitudinalshift$block_natord+(n*nx),l_th_latitudinalshift$x_coord,l_th_latitudinalshift$y_coord+(n*deltaY),l_th_latitudinalshift$z_coord)
      names(n_th_CoordSuite) = c("x_index","y_index","z_index","block_natord","x_coord","y_coord", "z_coord")
      # Assign generated coordinates to placeholder
      Planarshifts[,((7*((SS-1)   +1   +(n-1))*S_lat)+    (7*(l-1))   +1):((7*((SS-1)   +1   +(n-1))*S_lat)+    (7*l))] = n_th_CoordSuite
      names(Planarshifts)[((7*((SS-1)   +1   +(n-1))*S_lat)+    (7*(l-1))   +1):((7*((SS-1)   +1   +(n-1))*S_lat)+    (7*l))] = c(paste("x_index",Position_north, sep = ""),paste("y_index",Position_north, sep = ""),paste("z_index",Position_north, sep = ""),paste("block_natord",Position_north, sep = ""),paste("x_coord",Position_north, sep = ""),paste("y_coord",Position_north, sep = ""),paste("z_coord",Position_north, sep = ""))
      # Sample the full-grid realization at generated coordinates
      n_th_data = data.frame(matrix(0,samplesize,1)) 
      for (a in 1:samplesize){
        n_th_data[a,] = simulatedfullgrid[as.numeric(rownames(simulatedfullgrid))==n_th_CoordSuite$block_natord[a],dimensionality+1]
      }
      names(n_th_data) = c("attribute")
      n_th_Sample_v1 = if (deltaX>0) data.frame(n_th_CoordSuite,n_th_data) else data.frame(n_th_CoordSuite[,-c(which(names(n_th_CoordSuite)=="x_index"),which(names(n_th_CoordSuite)=="x_coord"))],n_th_data)
      n_th_Sample_v2 = if (deltaY>0) n_th_Sample_v1 else n_th_Sample_v1[,-c(which(names(n_th_Sample_v1)=="y_index"),which(names(n_th_Sample_v1)=="y_coord"))]
      n_th_Sample_final = if (deltaZ>0) n_th_Sample_v2 else n_th_Sample_v2[,-c(which(names(n_th_Sample_v2)=="z_index"),which(names(n_th_Sample_v2)=="z_coord"))]
      
      
      # Export sample data
      write.csv(n_th_Sample_final,paste0("Sample_",Position_north,".csv"), row.names = F)
      }
  }
}

# Upward shifts of all Planar shifts
if(US>1){
  # Placeholder
  UpShifts = data.frame(matrix(0,samplesize,(7*(US-1)*S_planar)))
  # Determine the range of positions for upward shiftings
  Up_ends = (US-1)*S_planar 
  Up_begins = 1
  for(u in 1:(US-1)){
    for(p in 1:S_longit){
      for(l in 1:S_lat){
        # Determine position for uth shifting
        Position_up = Up_ends-((u-1)*S_planar)-((p-1)*S_lat)-(l-1)
        # Extract the pth planar shifting to be shifted upward
        p_th_planarshift = Planarshifts[,(7*((S_planar)-((p-1)*S_lat)-l)+1):(7*((S_planar)-((p-1)*S_lat)-(l-1)))] 
        names(p_th_planarshift) = c("x_index","y_index","z_index","block_natord","x_coord","y_coord", "z_coord")
        # Generate the uth coordinate suite
        u_th_CoordSuite = data.frame(p_th_planarshift$x_index,p_th_planarshift$y_index,p_th_planarshift$z_index-u,p_th_planarshift$block_natord-(u*nx*ny),p_th_planarshift$x_coord,p_th_planarshift$y_coord,p_th_planarshift$z_coord-(u*deltaZ))
        names(u_th_CoordSuite) = c("x_index","y_index","z_index","block_natord","x_coord","y_coord", "z_coord")
        # Assign generated coordinates to placeholder
        UpShifts[,(7*(((US-1)*S_planar)-((u-1)*S_planar)-((p-1)*S_lat)-l)+1):(7*(((US-1)*S_planar)-((u-1)*S_planar)-((p-1)*S_lat)-(l-1)))] = u_th_CoordSuite
        names(UpShifts)[(7*(((US-1)*S_planar)-((u-1)*S_planar)-((p-1)*S_lat)-l)+1):(7*(((US-1)*S_planar)-((u-1)*S_planar)-((p-1)*S_lat)-(l-1)))] = c(paste("x_index",Position_up, sep = ""),paste("y_index",Position_up, sep = ""),paste("z_index",Position_up, sep = ""),paste("block_natord",Position_up, sep = ""),paste("x_coord",Position_up, sep = ""),paste("y_coord",Position_up, sep = ""),paste("z_coord",Position_up, sep = ""))
        # Sample the full-grid realization at generated coordinates
        u_th_data = data.frame(matrix(0,samplesize,1)) 
        for (b in 1:samplesize){
          u_th_data[b,] = simulatedfullgrid[as.numeric(rownames(simulatedfullgrid))==u_th_CoordSuite$block_natord[b],dimensionality+1]
        }
        names(u_th_data) = c("attribute")
        u_th_Sample_v1 = if (deltaX>0) data.frame(u_th_CoordSuite,u_th_data) else data.frame(u_th_CoordSuite[,-c(which(names(u_th_CoordSuite)=="x_index"),which(names(u_th_CoordSuite)=="x_coord"))],u_th_data)
        u_th_Sample_v2 = if (deltaY>0) u_th_Sample_v1 else u_th_Sample_v1[,-c(which(names(u_th_Sample_v1)=="y_index"),which(names(u_th_Sample_v1)=="y_coord"))]
        u_th_Sample_final = if (deltaZ>0) u_th_Sample_v2 else u_th_Sample_v2[,-c(which(names(u_th_Sample_v2)=="z_index"),which(names(u_th_Sample_v2)=="z_coord"))]
        
        
        # Export sample data
        write.csv(u_th_Sample_final,paste0("Sample_",Position_up,".csv"), row.names = F)
        }
    }
  }
} else {
  UpShifts = c()
} 


# Downward shifts of all Planar shifts
if(DS>0){
  # Placeholder
  DownShifts = data.frame(matrix(0,samplesize,(7*DS*S_planar)))
  # Determine the range of positions for downward shiftings
  Down_begins = (US*S_planar)+1 
  Down_ends = S_total
  for(d in 1:DS){
    for(p in 1:S_longit){
      for(l in 1:S_lat){
        # Determine position for dth shifting
        Position_down = Down_begins+((d-1)*S_planar)+((p-1)*S_lat)  +(l-1)
        # Extract the pth planar shifting to be shifted downward
        p_th_planarshift = Planarshifts[,(7*(((p-1)*S_lat)+(l-1))+1):(7*(((p-1)*S_lat)+l))] 
        names(p_th_planarshift) = c("x_index","y_index","z_index","block_natord","x_coord","y_coord", "z_coord")
        # Generate the dth coordinate suite
        d_th_CoordSuite = data.frame(p_th_planarshift$x_index,p_th_planarshift$y_index,p_th_planarshift$z_index+d,p_th_planarshift$block_natord+(d*nx*ny),p_th_planarshift$x_coord,p_th_planarshift$y_coord,p_th_planarshift$z_coord+(d*deltaZ))
        names(d_th_CoordSuite) = c("x_index","y_index","z_index","block_natord","x_coord","y_coord", "z_coord")
        # Assign generated coordinates to placeholder
        DownShifts[,(7*(((d-1)*S_planar)+((p-1)*S_lat)+(l-1))+1):(7*(((d-1)*S_planar)+((p-1)*S_lat)+l))] = d_th_CoordSuite
        names(DownShifts)[(7*(((d-1)*S_planar)+((p-1)*S_lat)+(l-1))+1):(7*(((d-1)*S_planar)+((p-1)*S_lat)+l))] = c(paste("x_index",Position_down, sep = ""),paste("y_index",Position_down, sep = ""),paste("z_index",Position_down, sep = ""),paste("block_natord",Position_down, sep = ""),paste("x_coord",Position_down, sep = ""),paste("y_coord",Position_down, sep = ""),paste("z_coord",Position_down, sep = ""))
        # Sample the full-grid realization at generated coordinates
        d_th_data = data.frame(matrix(0,samplesize,1)) 
        for (c in 1:samplesize){
          d_th_data[c,] = simulatedfullgrid[as.numeric(rownames(simulatedfullgrid))==d_th_CoordSuite$block_natord[c],dimensionality+1]
        }
        names(d_th_data) = c("attribute")
        d_th_Sample_v1 = if (deltaX>0) data.frame(d_th_CoordSuite,d_th_data) else data.frame(d_th_CoordSuite[,-c(which(names(d_th_CoordSuite)=="x_index"),which(names(d_th_CoordSuite)=="x_coord"))],d_th_data)
        d_th_Sample_v2 = if (deltaY>0) d_th_Sample_v1 else d_th_Sample_v1[,-c(which(names(d_th_Sample_v1)=="y_index"),which(names(d_th_Sample_v1)=="y_coord"))]
        d_th_Sample_final = if (deltaZ>0) d_th_Sample_v2 else d_th_Sample_v2[,-c(which(names(d_th_Sample_v2)=="z_index"),which(names(d_th_Sample_v2)=="z_coord"))]
        
        
        # Export sample data
        write.csv(d_th_Sample_final,paste0("Sample_",Position_down,".csv"), row.names = F)
        }
    }
  }
} else {
  DownShifts = c()
}

#Aggregating all coordinates suites

GrandShifts = if (identical(nz,1)) data.frame(cbind(UpShifts,data.matrix(Planarshifts),DownShifts)) else data.frame(UpShifts,Planarshifts,DownShifts)
write.csv(GrandShifts, "Grand_Shifts.csv", row.names = F)

return(GrandShifts)

}
