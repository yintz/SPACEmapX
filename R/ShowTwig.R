# before you run, you must install the SptialinferCNV along with all the relavent package
# and fully understand what you are doing in order to use this package.
# this function is designed to select the FFPE samples based on the dendrogram tree's trig.
# Once you select the the trig, type the trig ID in the function, it will generate the CSV 
# file with Trig ID for you to load back to the Loupe Brower to view and check.
# the other function called ShowTrigSectionName stores the section name into CSV list.
# it is used to locate the section, which more used for further SptialInferCNV to refine.
# this is the common function to use.


ShowTwig <- function(twigID){
  
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
    Node6<-Node5[,c(2,3)]
    #write.csv(Node5[,c(2,3)], paste("FFPE_V2_Paitent10_",unique(Node4$V1)[i],".csv", sep = ""))
    write.csv(Node6, paste("FFPE_V2_Paitent10_",unique(Node4$V1)[i],".csv", sep = ""),row.names = FALSE)
  }
}
