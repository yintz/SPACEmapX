LoadSTinfo <- function(csv_file_path) {
# check if the file exits
  if (!file.exists(file_path)) {
    stop("Error: the file doesnt exist, please check the path.")
  }
  
  # read csv file and first row as column names
  data <- read.csv(file_path, header = TRUE, stringsAsFactors = FALSE)
  
  # return
  return(data)
}

# df <- LoadSTinfo("your_file.csv")

# check first few lines.
# head(df)



  
#   allSTtable<-read.csv(csvfile)
#   allSTtable
#   n=1
#   while(n<dim(allSTable)[1,])
#   {
#     ENSBMLID_Counts[n]<-paste(allSTtable$SampleID[n],'ENSBMLID_Counts',sep = "_")
#     ENSBMLID_Counts[n]<-ImportCountData("allSTtable$SampleID[n]", allSTtable$SampleGeneExpressionFile[n]")
# Histology[n]<-paste(allSTtable$SampleID[n],'Histology',sep = "_")
# Histology[n]<- ImportHistologicalAnnotations(allSTtable$SampleID[n],  allSTtable$annotation[n])
# n=n+1
# }
# 
# A1_Joined_Counts <- MergingCountAndAnnotationData("A1",A1_Histology, A1_ENSBMLID_Counts)
# A2_Joined_Counts <- MergingCountAndAnnotationData("A2",A2_Histology, A2_ENSBMLID_Counts)
# A3_Joined_Counts <- MergingCountAndAnnotationData("A3",A3_Histology, A3_ENSBMLID_Counts)
# A4_Joined_Counts <- MergingCountAndAnnotationData("A4",A4_Histology, A4_ENSBMLID_Counts)
# A5_Joined_Counts <- MergingCountAndAnnotationData("A5",A5_Histology, A5_ENSBMLID_Counts)
# P1_Joined_Counts <- MergingCountAndAnnotationData("P1",P1_Histology, P1_ENSBMLID_Counts)
# P2_Joined_Counts <- MergingCountAndAnnotationData("P2",P2_Histology, P2_ENSBMLID_Counts)
# P3_Joined_Counts <- MergingCountAndAnnotationData("P3",P3_Histology, P3_ENSBMLID_Counts)
# P4_Joined_Counts <- MergingCountAndAnnotationData("P4",P4_Histology, P4_ENSBMLID_Counts)
# LCI1_Joined_Counts <- MergingCountAndAnnotationData("LCI1",LCI1_Histology, LCI1_ENSBMLID_Counts)
# LII1_Joined_Counts <- MergingCountAndAnnotationData("LII1",LII1_Histology, LII1_ENSBMLID_Counts)
# 
# 
# 
# Counts_joined <- A1_Joined_Counts %>% replace(., is.na(.), 0)
# Counts_joined <- Counts_joined %>% column_to_rownames(., var = "Genes")
# write.table(Counts_joined, "countMatrix.tsv", sep = "\t")
# MergeAllAnnotation<-A1_Joined_Counts
# MergedAll_Final <- FinalAnnotations(MergeAllAnnotation, Counts_joined)
# write.table(MergedAll_Final, "Annotations.tsv", 
#             sep = "\t",
#             quote = FALSE, 
#             col.names = FALSE, 
#             row.names = FALSE)
# 
# 
# Ref_infCNV <- infercnv::CreateInfercnvObject(raw_counts_matrix="countMatrix.tsv", 
#                                              gene_order_file="./siCNV_GeneOrderFile.tsv",
#                                              annotations_file="./Annotations.tsv",
#                                              delim="\t",
#                                              ref_group_names= NULL,
#                                              chr_exclude = c("ChrM"))
# 
# Ref_infCNVNoGroup = infercnv::run(Ref_infCNVNoGroup,
#                                   cutoff=0.1,
#                                   out_dir="./Clone", 
#                                   num_threads = 10,
#                                   denoise=T,
#                                   output_format = "png",
#                                   hclust_method='ward.D2',
#                                   cluster_by_groups=F,
#                                   analysis_mode = "samples",
#                                   tumor_subcluster_partition_method = "qnorm",
#                                   HMM=F)
# 
# 
# 
# 
# 
# 
#   
# A1_ENSBMLID_Counts<-ImportCountData("A1", "./A1/filtered_feature_bc_matrix.h5")
# A2_ENSBMLID_Counts<-ImportCountData("A2", "./A2/filtered_feature_bc_matrix.h5")
# A3_ENSBMLID_Counts<-ImportCountData("A3", "./A3/filtered_feature_bc_matrix.h5")
# A4_ENSBMLID_Counts<-ImportCountData("A4", "./A4/filtered_feature_bc_matrix.h5")
# A5_ENSBMLID_Counts<-ImportCountData("A5", "./A5/filtered_feature_bc_matrix.h5")
# P1_ENSBMLID_Counts<-ImportCountData("P1", "./P1/filtered_feature_bc_matrix.h5")
# P2_ENSBMLID_Counts<-ImportCountData("P2", "./P2/filtered_feature_bc_matrix.h5")
# P3_ENSBMLID_Counts<-ImportCountData("P3", "./P3/filtered_feature_bc_matrix.h5")
# P4_ENSBMLID_Counts<-ImportCountData("P4", "./P4/filtered_feature_bc_matrix.h5")
# LCI1_ENSBMLID_Counts<-ImportCountData("LCI1", "./LCI1/filtered_feature_bc_matrix.h5")
# LII1_ENSBMLID_Counts<-ImportCountData("LII1", "./LII1/filtered_feature_bc_matrix.h5")
# 
# 
# A1_Histology <- ImportHistologicalAnnotations("A1", "./A1/A1.csv")
# A2_Histology <- ImportHistologicalAnnotations("A2", "./A2/A2.csv")
# A3_Histology <- ImportHistologicalAnnotations("A3", "./A3/A3.csv")
# A4_Histology <- ImportHistologicalAnnotations("A4", "./A4/A4.csv")
# A5_Histology <- ImportHistologicalAnnotations("A5", "./A5/A5.csv")
# P1_Histology <- ImportHistologicalAnnotations("P1", "./P1/P1.csv")
# P2_Histology <- ImportHistologicalAnnotations("P2", "./P2/P2.csv")
# P3_Histology <- ImportHistologicalAnnotations("P3", "./P3/P3.csv")
# P4_Histology <- ImportHistologicalAnnotations("P4", "./P4/P4.csv")
# LCI1_Histology <- ImportHistologicalAnnotations("LCI1", "./LCI1/LCI1.csv")
# LII1_Histology <- ImportHistologicalAnnotations("LII1", "./LII1/LII1.csv")
# 
# 
# 
# 
#   
#   
#     # 函数体
#     # 执行的代码块
#     return(output)
# }
