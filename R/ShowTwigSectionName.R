ShowTwigSectionName <- function(twigID){
  
  Node<- SelectingSubTreeData(my.subtrees, twigID)
  #Node<- SelectingSubTreeData(my.subtrees, 2404)
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
    Node5[,3]<-unique(Node4$V1)[i]
    Node6<-Node5[,c(2,3)]
    #write.csv(Node5[,c(2,3)], paste("FFPE_V2_Paitent10_",unique(Node4$V1)[i],".csv", sep = ""))
    write.csv(Node6, paste("FFPE_V2_Paitent10_",unique(Node4$V1)[i],".csv", sep = ""),row.names = FALSE)
  }
}
