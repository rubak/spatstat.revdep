#' Generates points of natural origin.
#'
#' This function generates points of natural origin around existent points by sampling their distance from the original propagule from a user-defined dispersal kernel (as set in function SDS(), the output of which is to be entered in this function as argument 'Delta'.
#'
#' @param PresLoc data frame with coordinates of sighting locations. Column 2 must be 'y' (latitude on a projected coordinate system with meters as distance unit of measure) and column 3 must be 'x' (longitude on a projected coordinate system with meters as distance unit of measure).
#' @param N total number of points to generate.
#' @param Deltas y and x shifts as generated from function SDS(). This function samples distance and angle of new points form existing ones from this argument.
#' @param year the year lable to attach to the generated points
#'
#' @return data frame of combined source ('PresLoc' argument) and new locations.
#'
#' @author Luca Butikofer
#'
#' @export


NPG<-function(PresLoc,N,Deltas,year){


  if(N>nrow(PresLoc)){
    NEWLOC<-list(NA)
    for(i in 1:nrow(PresLoc)){
      new<-sample(rownames(Deltas),N)
      new<-Deltas[new,]
      new2<- data.frame(matrix(nrow=N, ncol=2))
      for(j in 1:N){
        new2[j,]<-PresLoc[i,2:3]+new[j,1:2]
      }
      NEWLOC[[i]]<-new2
    }
    newloc<-do.call("rbind",NEWLOC)
    newloc<-newloc[sample(row.names(newloc),N),]
    newloc$year<-rep(year,nrow(newloc))
    newloc$Pnat<-rep(1,nrow(newloc))
    newloc$species<-rep(PresLoc$species[1],nrow(newloc))
    newloc$Dist<-rep("NPG",nrow(newloc))
    colnames(newloc)<-c("y","x","year","Pnat","species","Dist")
    newloc<-newloc[,c(3,1,2,5,4,6)]
    return(newloc)

  }else{

    PresLoc2<-PresLoc[sample(row.names(PresLoc),N),]
    new<-Deltas[sample(row.names(Deltas),N),]
    PresLoc2$x<-PresLoc2$x+new[,1]
    PresLoc2$y<-PresLoc2$y+new[,2]
    PresLoc2$year<-rep(year,nrow(PresLoc2))
    PresLoc2$Pnat<-rep(1,nrow(PresLoc2))
    PresLoc2$Dist<-rep('NPG',nrow(PresLoc2))
    return(PresLoc2)

  }
}
