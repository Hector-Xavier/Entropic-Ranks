library("RankProd")
library("entropy")
library("factoextra")

entropic_analysis <- function(ordered_vector,step_up=1,window_size,bins,verbose=FALSE)
{
  differences <- ordered_vector[seq(2,length(ordered_vector))]-ordered_vector[seq(1,length(ordered_vector)-1)]
  
  entropy_plotter <- vector(length=floor((length(differences)-window_size)/step_up))
  for (i in 0:(length(entropy_plotter)-1))
  {
    entropy_plotter[i+1] <- entropy(discretize(differences[(i*step_up+1):(i*step_up+window_size)],numBins=bins),method="Laplace")
  }
  
  temp <- eclust(entropy_plotter, "kmeans", k = 2, nstart = 200, graph = FALSE)
  
  if (verbose)
  {
    message("Calculating entropies of ",length(entropy_plotter)," ovelapping windows on a list of length: ",length(ordered_vector))
    message("Suggested cutoff at element no ",seq(length(temp$cluster))[temp$cluster!=temp$cluster[1]][1]-1," with a mean entropy of: ",mean(entropy_plotter[1:seq(length(temp$cluster))[temp$cluster!=temp$cluster[1]][1]]))
    message("Last 1/3 infimum: ",min(entropy_plotter[floor(length(entropy_plotter)*2/3):length(entropy_plotter)]))
    barplot(entropy_plotter,xlab=c("Observation step = ",step_up),ylab=c("Window size = ",window_size),main=c(bins," bins"),names.arg=seq(length(entropy_plotter))*step_up)
    barplot(temp$silinfo$widths$sil_width,ylim=c(-0.3,1))
  }
  
  return((seq(length(temp$cluster))[temp$cluster!=temp$cluster[1]][1]-1)*step_up)
}

isolate_significant_elements <- function(list_under_analysis,resolution=1,process_log=FALSE)
{
  suggested_surfaces <- c()
  
  if (process_log)
    par(mfrow=c(2,1))
  
  for (i in 1:5)
  {
    for (j in 5:20)
    {
      if (process_log)
      {
        message("Testing for ",i*5," bins with window size equal to ",j*10)
      }
      suggested_surfaces <- c(suggested_surfaces,entropic_analysis(ordered_vector=list_under_analysis,step_up=resolution,bins=i*5,window_size=j*10,verbose=process_log))
    }
  }
  
  if (process_log)
  {
    message("Suggested cutoff surfaces:")
    print(table(suggested_surfaces))
    par(mfrow=c(1,1))
  }
  return(as.integer(rownames(table(suggested_surfaces))[table(suggested_surfaces) == max(table(suggested_surfaces))])[1])
}

#Code written by Hector-Xavier de Lastic
#Development & testing by Hector-Xavier de Lastic & Irene Liampa
#Contact:
#hector.xavier.de.lastic@gmail.com
#irini.liampa@gmail.com
