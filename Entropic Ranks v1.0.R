require("RankProd")
require("entropy")
require("factoextra")

entropic_analysis <- function(ordered_vector,step_up=1,window_size,bins,verbose=FALSE)
{
  differences <- ordered_vector[seq(2,length(ordered_vector))]-ordered_vector[seq(1,length(ordered_vector)-1)]
  
  entropy_plotter <- vector(length=floor((length(differences)-window_size)/step_up))
  for (i in 0:(length(entropy_plotter)-1))
  {
    entropy_plotter[i+1] <- entropy(discretize(differences[(i*step_up+1):(i*step_up+window_size)],numBins=bins),method="Laplace")
  }
  
  temp <- eclust(entropy_plotter, "kmeans", k=2, nstart=200, graph=FALSE)
  
  if (verbose)
  {
    message("Calculating entropies of ",length(entropy_plotter)," ovelapping windows")
    message("Suggested cutoff at feature no ",seq(length(temp$cluster))[temp$cluster!=temp$cluster[1]][1]-1,", at a mean entropy of ",mean(entropy_plotter[1:seq(length(temp$cluster))[temp$cluster!=temp$cluster[1]][1]]))
    message("Last 1/3 entropy infimum: ",min(entropy_plotter[floor(length(entropy_plotter)*2/3):length(entropy_plotter)]))
    barplot(entropy_plotter,xlab=c("Granularity = ",step_up),ylab=c("Window size = ",window_size),main=c(bins," bins"),names.arg=seq(length(entropy_plotter))*step_up)
    barplot(temp$silinfo$widths$sil_width,ylim=c(-0.3,1),ylab="Silhouette width",xlab="Clustered elements",main="K-means clustering quality")
  }
  
  return((seq(length(temp$cluster))[temp$cluster!=temp$cluster[1]][1]-1)*step_up)
}

isolate_significant_elements <- function(data_under_analysis,granularity=1,process_log=FALSE)
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
        message("Reiterating for ",i*5," bins and a sliding window of ",j*10," features")
      }
      suggested_surfaces <- c(suggested_surfaces,entropic_analysis(ordered_vector=data_under_analysis,step_up=granularity,bins=i*5,window_size=j*10,verbose=process_log))
    }
  }
  
  if (process_log)
  {
    print(table(suggested_surfaces))
    message("Most consistent cutoff point: feature no ",as.integer(rownames(table(suggested_surfaces))[table(suggested_surfaces) == max(table(suggested_surfaces))])[1],".")
    par(mfrow=c(1,1))
  }
  return(as.integer(rownames(table(suggested_surfaces))[table(suggested_surfaces) == max(table(suggested_surfaces))])[1])
}

entropic_ranks <- function(data_under_analysis,population_vector,data_origin=NULL,granularity=1,process_log=FALSE,create_output_files=FALSE,is_logged=TRUE,logbase=2,huge_feature_list=FALSE)
{
  if (is.null(data_origin))
  {
    data_origin <- rep(1,length(population_vector))
  }
  
  comparison <- RPadvance(data_under_analysis,cl=population_vector,origin=data_origin,logged=is_logged,na.rm=FALSE,gene.names=rownames(data_under_analysis),plot=process_log,huge=TRUE)

  if (huge_feature_list)
  {
    rank_product_lists <- topGene(comparison,method="pfp",num.gene=5000,logged=is_logged,logbase=logbase,gene.names=rownames(data_under_analysis))
  }else{
    rank_product_lists <- topGene(comparison,cutoff=0.99,method="pfp",num.gene=NULL,logged=is_logged,logbase=logbase,gene.names=rownames(data_under_analysis))
  }
  
  rank_product_lists$Table1 <- rank_product_lists$Table1[1:isolate_significant_elements(rank_product_lists$Table1[,2],granularity,process_log),]
  rank_product_lists$Table2 <- rank_product_lists$Table2[1:isolate_significant_elements(rank_product_lists$Table2[,2],granularity,process_log),]
  
  if (create_output_files)
  {
    if (dim(rank_product_lists$Table1)[1]>0)
    {
      write.table(file="Downregulated list.txt",rank_product_lists$Table1[,c(3:5)],sep="\t",quote=FALSE,row.names=TRUE,col.names=TRUE)
    }
    
    if (dim(rank_product_lists$Table2)[1]>0)
    {
      write.table(file="Upregulated list.txt",rank_product_lists$Table2[,3:5],sep="\t",quote=FALSE,row.names=TRUE,col.names=TRUE)
    }
  }
  return(rank_product_lists)
}
