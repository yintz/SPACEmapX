
ShowTwigSpotsList <- function(...){
  # Capture the arguments as a list
  args <- list(...)
  
  # Check if all arguments are numeric
  if(!all(sapply(args, is.numeric))) {
    stop("All arguments must be numeric.")
  }
  
  if(length(args) == 0) {
    stop("At least one twig number must be provided.")
  }
  
  # Convert the list of arguments to a data frame
  # Convert the list of arguments to a data frame
  df <- data.frame(Number = unlist(args))
  #df <- data.frame(Number = unlist(args))
  return(df)
  i<-1
  while i<-nrow(df) {
  Node<- SelectingSubTreeData(my.subtrees, df[i,1])
  Node$Barcode<-gsub("\\.","-",Node$Barcode)
  Node2<-as.data.frame(str_split_fixed(Node$Barcode, "_", 2))
  Node3<-cbind(Node2,Node$Node)
  Node4<-Node3
  length(unique(Node4$V1))
  
  for (i in c(1:length(unique(Node4$V1)))){
    print(i)
    print(unique(Node4$V1)[i])
    Node5<-subset(Node4,Node4$V1==unique(Node4$V1)[i])
    colnames(Node5)<-c("no","Barcode","Twig")
    Node6<-Node5[,c(2,3)]
    write.csv(Node6, paste("FFPE_V2_Paitent_",unique(Node4$V1)[i],".csv", sep = ""),row.names = FALSE)
    i=i+1
  }
  }
}
