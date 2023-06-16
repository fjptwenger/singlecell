library(Seurat)
library(ggplot2)
library(scales)
library(patchwork)
library(grid)
library(dplyr)
library(SingleR)
library(celldex)
library(scRNAseq)
library(scater)
library(clustree)

path<-'scRNA-seq/'
setwd(path)

#The data were downloaded from the webpage:https://ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE200997
#dinfo<-read.csv('GSE200997_GEO_processed_CRC_10X_cell_annotation.csv')
#dat<-read.csv('GSE200997_GEO_processed_CRC_10X_raw_UMI_count_matrix.csv')

load('data/GSE200997.Rdata')
dat<-GSE200997_expr
dinfo<-GSE200997_saminfo

#Statistical information
patients<-as.character(regmatches(colnames(GSE200997_expr),gregexpr('^[BT]_[a-zA-Z1-9]+',colnames(GSE200997_expr))))
table(substr(unique(patients),1,1)) # 7 Normal 15 Tumor 
table(dinfo$Location) # 28,499 LCC cells, 21,360 RCC cells
table(dinfo$Condition)# 18,273 cells from Normal samples, 31,586 cells from Tumorl samples
table(dinfo[dinfo$Condition=='Tumor',]$Location) #16448 Left, 15138 Right


#seurat object
tumor_id<-dinfo[dinfo$Condition=='Tumor',"X"]
left<-dinfo[dinfo$Condition=='Tumor' & dinfo$Location=='Left' ,"X"]
right<-dinfo[dinfo$Condition=='Tumor' & dinfo$Location=='Right' ,"X"]
adat<-dat[,tumor_id]
dim(adat) #23828 31586
tumor<-CreateSeuratObject(counts = adat, min.cells = 3, project = "GSE200997")

tumor[["percent.mt"]] <- PercentageFeatureSet(tumor, pattern = "^MT-")
#add location
rownames(dinfo)<-dinfo$X
dinfo_L<-dinfo[rownames(tumor@meta.data),]
tumor@meta.data$group<-ifelse(dinfo_L$Location=='Right','RCC','LCC')
dim(tumor@meta.data)#31586     5

#save(tumor,file='tumor.Rdata')
#load('tumor.Rdata')
dat <- subset(tumor, subset=nFeature_RNA > 50 & nFeature_RNA < 6000 &  percent.mt < 20)
dim(dat@meta.data) #30411, 5

sample<-as.character(regmatches(rownames(dat@meta.data),gregexpr('^[BT]_[a-zA-Z1-9]+',rownames(dat@meta.data))))
dat <- AddMetaData(object = dat, metadata = sample, col.name = 'Sample')

dir.create('Figures/')
tiff(file="Figures/Figure 1A.group_feature.tiff",width=15,height=10,units="cm",compression="lzw",bg="white",res=600)
p=VlnPlot(object = dat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),group.by = "group",pt.size=0)
print(p)
dev.off()

dat <- NormalizeData(dat, normalization.method = "LogNormalize", scale.factor = 10000)

dat <- FindVariableFeatures(dat, selection.method = "vst", nfeatures = 3000)
#scale for PCA
dat <- ScaleData(dat,vars.to.regress = c('nFeature_RNA', 'percent.mt'))

dat <- RunPCA(dat, npcs=50)

dat<- FindNeighbors(dat, dims = 1:20)
#dat <- FindClusters(dat, resolution = 0.5) #resolution == (0.4-1.2)
dat <- FindClusters(dat, resolution = c(seq(0,1.6,0.2)))  

p<-clustree(dat@meta.data,prefix="RNA_snn_res.") 
tiff(file='Figures/Figure1B.clustree.tiff',width=20,height=20,units="cm",compression="lzw",bg="white",res=600)
print(p)
dev.off()

Idents(object = dat)<-"RNA_snn_res.0.6"
dat<- RunUMAP(dat, dims = 1:20)

tiff(file='Figure/Figure1C.umap_cluster.tiff',width=25,height=15,units="cm",compression="lzw",bg="white",res=600)
p <- DimPlot(dat, reduction = "umap", label = TRUE, repel = TRUE) 
print(p)
dev.off()
save(dat,file = 'data/dat.Rdata')

numd<-data.frame(dat@meta.data$group,dat@meta.data$RNA_snn_res.0.6)
colnames(numd)<-c('group','cluster')

# cell_annotation: SingleR 
#hpca.se <- celldex::HumanPrimaryCellAtlasData() #download online
#save(hpca.se,file='data/singleR_anno.Rdata')
load('data/singleR_anno.Rdata')
anno_dat<-SingleR(test = dat@assays$RNA@data, ref = hpca.se, labels = hpca.se$label.main)
#save(anno_dat, file = "data/anno_dat.RData")

#Statistical annotation results
#load('data/anno_dat.RData')
table(anno_dat$pruned.labels)
#Astrocyte               B_cell           BM & Prog.         Chondrocytes                  CMP 
#1                 5088                   11                   34                   91 
#DC Embryonic_stem_cells    Endothelial_cells     Epithelial_cells         Erythroblast 
#113                   38                  377                 5955                    6 
#Fibroblasts          Gametocytes                  GMP          Hepatocytes            HSC_CD34+ 
#  194                    3                   24                   59                   12 
#iPS_cells        Keratinocytes           Macrophage                  MEP             Monocyte 
#33                    4                  205                    3                  489 
#MSC            Myelocyte              Neurons          Neutrophils              NK_cell 
#4                    3                   19                   37                 1075 
#Pre-B_cell_CD34-     Pro-B_cell_CD34+        Pro-Myelocyte  Smooth_muscle_cells              T_cells 
#32                   37                    8                   30                15988 
#Tissue_stem_cells 
#364 

#the distribution of clusters
table(dat@meta.data$RNA_snn_res.0.6) #optimal RNA_snn_res.=0.6
#0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20 
#3929 3882 2931 2864 2297 2268 2241 2147 2077  841  764  751  738  502  436  390  374  352  262  183  182 

#cell and cluster
tmp1<-data.frame(rownames(dat@meta.data),dat@meta.data$RNA_snn_res.0.6)

colnames(tmp1)<-c('cell_ID','clusters')
#cell and annotation
tmp2<-data.frame(anno_dat@rownames, anno_dat$pruned.labels)

colnames(tmp2)<-c('cell_ID','cell_name')

cell_cluster_name<-merge(tmp1,tmp2,by='cell_ID')

write.table(table(cell_cluster_name[,c('cell_name','clusters')]),file='cell_cluster_sta.xls',quote = F,sep='\t')

#according to the table(cell_cluster_sta.xls), the  clusters(0-20) could be annotated as follow:
anno<-c('T cells','T cells','B cells','T cells','Epithelial cells',
	        'B cells','T cells','T cells','Epithelial cells','Macrophage',
		        'Epithelial cells','T cells','Epithelial cells','NK cells','Fibroblasts',
		        'Endothelial cells','T cells','Epithelial cells','T cells','Tissue stem cells','B cells')
names(anno) <- levels(dat)

cluster_anno<-data.frame(0:20,anno)
colnames(cluster_anno)<-c('cluster','cell_types')
group_celltypes<-merge(numd,cluster_anno,by='cluster')


dat2<-RenameIdents(dat, anno)
tiff(file='Figures/Figure1D.anno_cluster_umap.tiff',width=25,height=15,units="cm",compression="lzw",bg="white",res=600)
p<-DimPlot(dat2, reduction = "umap", label = TRUE, pt.size = 0.5)  #不删除rmcellid绘制
print(p)
dev.off()


tiff(file='Figures/Figure1E.celltypes_location_ratio.tiff',width=10,height=15,units="cm",compression="lzw",bg="white",res=600)
p<-ggplot(data = group_celltypes, aes(x = cell_types, fill = factor(group))) + geom_bar( position = "fill")+
	  scale_fill_manual(values = c('#fcaf17','#009ad6'))+
	    geom_hline(aes(yintercept = 0.45), colour = "red", size=0.5,linetype="dashed") +
	      geom_hline(aes(yintercept = 0.55), colour = "red", size=0.5,linetype="dashed") +
	        theme(panel.background = element_blank(),panel.grid.major= element_line(color = "white"),panel.grid.minor =element_line(color= "white"),legend.title = element_blank(),legend.position='top')+
		  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
	  print(p)
	  dev.off()

	  #obtained the DEGs of each annotated cell type by comparing the LCC and RCC cells
	  cluster0.markers <- FindMarkers(dat2, ident.1="LCC",ident.2="RCC", subset.ident = c("T cells"),
					                                        logfc.threshold = 0.585, # FC = 0.0625
										                                      group.by = 'group',
										                                      test.use = "roc")
	  cluster1.markers <- FindMarkers(dat2, ident.1="LCC",ident.2="RCC", subset.ident =c('B cells'),
					                                        logfc.threshold = 0.585, # FC = 0.0625
										                                      group.by = 'group',
										                                      test.use = "roc")
	  cluster2.markers <- FindMarkers(dat2, ident.1="LCC", ident.2="RCC",subset.ident =c('Epithelial cells'),
					                                        logfc.threshold = 0.585, # FC = 0.0625
										                                      group.by = 'group',
										                                      test.use = "roc")
	  cluster3.markers <- FindMarkers(dat2, ident.1="L", subset.ident =c('NK cells'),
					                                        logfc.threshold = 0.585, # FC = 0.0625
										                                      group.by = 'group',
										                                      test.use = "roc")
	  cluster4.markers <- FindMarkers(dat2, ident.1="LCC", ident.2="RCC",subset.ident =c('Macrophage'),
					                                        logfc.threshold = 0.585, # FC = 0.0625
										                                      group.by = 'group',
										                                      test.use = "roc")
	  cluster5.markers <- FindMarkers(dat2, ident.1="LCC",ident.2="RCC", subset.ident =c('Fibroblasts'),
					                                  logfc.threshold = 0.585, # FC = 0.0625
									                                  group.by = 'group',
									                                  test.use = "roc")
	  cluster6.markers <- FindMarkers(dat2, ident.1="LCC", ident.2="RCC",subset.ident =c('Endothelial cells'),
					                                  logfc.threshold = 0.585, # FC = 0.0625
									                                  group.by = 'group',
									                                  test.use = "roc")
	  cluster7.markers <- FindMarkers(dat2, ident.1="L", subset.ident =c('Tissue stem cells'),
					                                  logfc.threshold = 0.585, # FC = 0.0625
									                                  group.by = 'group',
									                                  test.use = "roc")
	  dir.create('Tables/')
	  write.table(cluster0.markers,file='Tables/Table_S2.tcell_DEG.xls',quote = F,sep = '\t')
	  write.table(cluster1.markers,file='Tables/Table_S2.bcell_DEG.xls',quote = F,sep = '\t')
	  write.table(cluster2.markers,file='Tables/Table_S2.Epithelial_cells_DEG.xls',quote = F,sep = '\t')
	  write.table(cluster3.markers,file='Tables/Table_S2.NK_cell_DEG.xls',quote = F,sep = '\t')
	  write.table(cluster4.markers,file='Tables/Table_S2.Macrophage_DEG.xls',quote = F,sep = '\t')
	  write.table(cluster5.markers,file='Tables/Table_S2.Fibroblasts_DEG.xls',quote = F,sep = '\t')
	  write.table(cluster6.markers,file='Tables/Table_S2.Endothelial_cells.xls',quote = F,sep = '\t')
	  write.table(cluster7.markers,file='Tables/Table_S2.Tissue_stem_cells_DEG.xls',quote = F,sep = '\t')

	  #Merge to obtain the DEGs
	  tcell_diffgene<-rownames(cluster0.markers)
	  bcell_diffgene<-rownames(cluster1.markers)
	  epithelial_diffgene<-rownames(cluster2.markers)
	  nkcell_diffgene<-rownames(cluster3.markers)
	  macrophage_diffgene<-rownames(cluster4.markers)
	  fibroblasts_diffgene<-rownames(cluster5.markers)
	  Endothelial_cells_diffgene<-rownames(cluster6.markers)
	  stemcells_diffgene<-rownames(cluster7.markers)

	  all_diffgene<-c(tcell_diffgene,bcell_diffgene,epithelial_diffgene,
			                 nkcell_diffgene,macrophage_diffgene,fibroblasts_diffgene,Endothelial_cells_diffgene,stemcells_diffgene
					               
					              )

	  all_diffgene_u<-unique(all_diffgene) #the DEGs list

	  #GO&KEGG analysis
	  library(clusterProfiler)
	  library(org.Hs.eg.db)
	  targetgenelist<-all_diffgene_u
	  addID<-bitr(targetgenelist, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
	  targetgenelist=addID[addID$SYMBOL %in% targetgenelist,"ENTREZID"]
	  go <- enrichGO(gene=targetgenelist, OrgDb = org.Hs.eg.db, ont='ALL',pvalueCutoff = 0.05, readable = T,qvalueCutoff = 0.05,keyType = 'ENTREZID')	
	  tiff(file='Figures/Figure1G.go.tiff',width=30,height=20,units="cm",compression="lzw",bg="white",res=600)
	  p=cnetplot(go, foldChange=targetgenelist,colorEdge=TRUE)
	  print(p)
	  dev.off()

	  kegg <- enrichKEGG(targetgenelist,organism = 'hsa',keyType = 'kegg',
			                        pvalueCutoff = 0.05,use_internal_data = T)

	  tiff(file='Figures/Figure1F.kegg.tiff',width=20,height=20,units="cm",compression="lzw",bg="white",res=600)
	  p=barplot(kegg,showCategory=20,drop=T)
	  print(p)
	  dev.off()

	  #Analyses of the prognostic genes for LCC and RCC patients based on the above DEGs list and TCGA-COAD cohort
	  #survival information of TCGA-COAD
	  #download website: https://xenabrowser.net/datapages/?cohort=GDC%20TCGA%20Colon%20Cancer%20(COAD)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
	  TCGA_survival<-read.table('data/TCGA-COAD.survival.tsv',sep='\t',header = T)
	  #phenotype information
	  TCGA_pheno<-read.table('data/TCGA-COAD.GDC_phenotype.tsv',sep='\t',header = T,quote = '"')
	  colnames(TCGA_pheno)[1]<-'sample'

	  TCGA_pheno_survival<-merge(TCGA_survival,TCGA_pheno,by='sample')

	  #According to the published paper[PMID:33980174]:
	  #the cecum, ascending, and hepatic belong to the right side, 
	  #and the splenic, descending, sigmoid, rectum, and junction belong to the left side.
	  table(TCGA_pheno_survival$site_of_resection_or_biopsy.diagnoses)
	  TCGA_pheno_survival$site<-
		    ifelse(TCGA_pheno_survival$site_of_resection_or_biopsy.diagnoses=='Ascending colon' | 
			              TCGA_pheno_survival$site_of_resection_or_biopsy.diagnoses=='Hepatic flexure of colon'|
				                 TCGA_pheno_survival$site_of_resection_or_biopsy.diagnoses=='Cecum','Right',
					          ifelse(
							            TCGA_pheno_survival$site_of_resection_or_biopsy.diagnoses=='Descending colon' |
									                 TCGA_pheno_survival$site_of_resection_or_biopsy.diagnoses=='Sigmoid colon','Left',
										            ifelse(
												                TCGA_pheno_survival$site_of_resection_or_biopsy.diagnoses=='Transverse colon','Transverse','unknown'
														             
														           )
								             )
				        )                               
	  TCGA_pheno_survival<-TCGA_pheno_survival[TCGA_pheno_survival$site!='unknown',]   
	  TCGA_pheno_survival<-TCGA_pheno_survival[TCGA_pheno_survival$OS.time >0,]

	  #Table 1
	  clinical_data<-read.table('data/COAD_clinicalMatrix.xls',header = T,sep='\t')
	  clinical_data$MSI<-ifelse(clinical_data$MSI_updated_Oct6201=="",clinical_data$CDE_ID_3226963,clinical_data$MSI_updated_Oct62011)
	  clinical_data$sample<-clinical_data$sampleID
	  TCGA_pheno_survival$sample<-substr(TCGA_pheno_survival$sample,1,15)
	  TCGA_pheno_survival_table<-merge(TCGA_pheno_survival,clinical_data[,c('sample','MSI')],by='sample')

	  TCGA_pheno_survival_table$Gender<-TCGA_pheno_survival_table$gender.demographic
	  TCGA_pheno_survival_table$Age<-ifelse(TCGA_pheno_survival_table$age_at_initial_pathologic_diagnosis >60,'>60','<=60')
	  library(stringr)
	  TCGA_pheno_survival_table$T<-str_extract(TCGA_pheno_survival_table$pathologic_T, "T[0-9]+")
	  TCGA_pheno_survival_table$M <- str_extract(TCGA_pheno_survival_table$pathologic_M, "M[0-9]*X*")
	  TCGA_pheno_survival_table$N<-str_extract(TCGA_pheno_survival_table$pathologic_N, "N[0-9]")
	  TCGA_pheno_survival_table$stage<-str_extract(TCGA_pheno_survival_table$tumor_stage.diagnoses, "i+v*")
	  TCGA_pheno_survival_table$microsatellite_status<-TCGA_pheno_survival_table$MSI
	  TCGA_pheno_survival_table$OS<-factor(TCGA_pheno_survival_table$OS,
					                                            levels=c(0,1),
										                                         labels=c("Alive", # 第一个作为参考组
																                                                "Death"))

	  #Table1
	  library(table1)
	  table1(~ factor(Gender) + factor(Age)+factor(microsatellite_status) +factor(stage) +
		          factor(T) +  factor(N) +  factor(M) | site, data=TCGA_pheno_survival_table)

	  #comparing the differece: left VS right vs Transverse
	  library(cowplot)
	  library(tidyverse)
	  library(ggplot2)
	  library(ggsci)
	  library(ggpubr)
	  TCGA_survival_compare<-data.frame(TCGA_pheno_survival$sample,TCGA_pheno_survival$OS.time,TCGA_pheno_survival$site)
	  colnames(TCGA_survival_compare)<-c('sample','OS.time','Group')

	  library(tidyr)
	  rownames(TCGA_survival_compare)<-TCGA_survival_compare$sample
	  TCGA_survival_compare$sample<-NULL
	  TCGA_survival_compare$Group<-ifelse(TCGA_survival_compare$Group=='Left','LCC',ifelse(TCGA_survival_compare$Group=='Right','RCC','Transverse'))
	    
	  TCGA_survival_compare$Group<-as.factor(TCGA_survival_compare$Group)

	  mydata<-TCGA_survival_compare %>% gather(key = "Group",value = "OS.time") %>% na.omit() 
	  my_comparisons <- list( c("LCC", "RCC"),c("RCC", "Transverse"),c("LCC", "Transverse"))

	  tiff(file='Figures/Figure S1.LOcation_OS.tiff',width=18,height=15,units="cm",compression="lzw",bg="white",res=600)
	  p<-ggboxplot(mydata,x='Group',y='OS.time',color='Group',palette = 'jama',add = 'jitter')
	  p+stat_compare_means(comparisons = my_comparisons)
	  dev.off()

	  #for LCC
	  TCGA_pheno_survival$sample<-substr(TCGA_pheno_survival$sample,1,15)
	  TCGA_pheno_survival<-dplyr::distinct(TCGA_pheno_survival,sample,.keep_all=TRUE)
	  rownames(TCGA_pheno_survival)<-TCGA_pheno_survival$sample

	  TCGA_pheno_survival_left<-TCGA_pheno_survival[TCGA_pheno_survival$site=='Left',]
	  group_list<-ifelse(substr(TCGA_pheno_survival_left$sample,14,15)<10,'Tumor','Normal')

	  #load('data/TCGA_colon_expr_filt.Rdata')
	  commone_sample_left<-intersect(TCGA_pheno_survival_left$sample,colnames(TCGA_colon_expr_filt))
	  tumor_sample_left<-commone_sample_left[as.numeric(substr(commone_sample_left,14,15))<10]

	  common_gene_left<-intersect(all_diffgene_u,rownames(TCGA_colon_expr_filt))
	  gene_survival_data_left<-cbind(TCGA_pheno_survival[tumor_sample_left,c('sample','OS','OS.time')],t(TCGA_colon_expr_filt[common_gene_left,tumor_sample_left]))

	  #save(gene_survival_data_left,file = 'gene_survival_data_left.Rdata')

	  load('data/gene_survival_data_left.Rdata')

	  dir.create('Figures/KM_left')
	  library(survminer)
	  library(survival)
	  surv_gene_left<-c()
	  for (i in 4:(ncol(gene_survival_data_left)-1)){
		    gene_survival_data_left$group<-ifelse(gene_survival_data_left[,i]>median(gene_survival_data_left[,i]),'high','low')
	    fit <- survfit(Surv(as.numeric(OS.time), as.numeric(OS)) ~gene_survival_data_left$group, data = gene_survival_data_left)
	      #print(c(colnames(gene_survival_data_left)[i],as.character(surv_pvalue(fit)[4])))
	      a<-as.character(surv_pvalue(fit)[4])
	      if(as.numeric(gsub('p = ','',a))<0.05){
		          surv_gene_left<-append(surv_gene_left,colnames(gene_survival_data_left)[i])
	          #print(paste(colnames(gene_survival_data_left)[i],as.character(surv_pvalue(fit)[4]),sep = '\t'))
	          #write.table(paste(colnames(gene_survival_data_left)[i],as.character(surv_pvalue(fit)[4]),sep = '\t'),file = survfile,append = T,quote = F,sep = '\t',row.names = F,col.names = F)
	          p=ggsurvplot(fit, 
			                        #pval = TRUE, 
			                        pval = toupper(surv_pvalue(fit)[4]),#p改为大写
						                 conf.int = TRUE,
						                 risk.table = TRUE, # Add risk table
								                  risk.table.col = "strata", # Change risk table color by groups
								                  linetype = "strata", # Change line type by groups
										                   surv.median.line = "hv", # Specify median survival
										                   ggtheme = theme_bw(), # Change ggplot2 theme
												                    legend.labs = c("High", "Low"), 
												                    legend.title = paste0("Expression of ",colnames(gene_survival_data_left)[i]), 
														                     palette = c( "#f03b20","#2c7fb8"))
		      
		      filename = paste0('Figures/KM_left/',colnames(gene_survival_data_left)[i],'.tiff')
		      tiff(file=filename,width=15,height=15,units="cm",compression="lzw",bg="white",res=600)
		          print(p)
		          dev.off()
			    }
	  }

	  #define the function for universal  Cox regression analysis
	  processUniCOX <- function(genelist, survival_gene_data) {
		    library(survival)
	    library(survminer)
	      genelist<-gsub(genelist,pattern = '-',replacement = '_')
	      uni_cox<-function(single_gene){
		          
		          formula<-as.formula(paste0('Surv(OS.time, OS) ~ ', single_gene))
	          surv_uni_cox<-summary(coxph(formula,data = survival_gene_data))
		      
		      ph_hypothesis_p<-cox.zph(coxph(formula,data = survival_gene_data))$table[1,3]
		      if(surv_uni_cox$coefficients[,5]<0.05 & ph_hypothesis_p>0.05){
			            single_cox_report<-data.frame(
								          'uni_cox_sig_genes'=single_gene,
									          'beta'=surv_uni_cox$coefficients[,1],
									          'Hazard_Ratio'=exp(surv_uni_cox$coefficients[,1]),
										          'z_pvalue'=surv_uni_cox$coefficients[,5],
										          'wald_pvalue'=as.numeric(surv_uni_cox$waldtest[3]),
											          'Likelihood_pvalue'=as.numeric(surv_uni_cox$logtest[3])
											        )
		            single_cox_report
			        }
		          
		        }
	        uni_cox_list<-lapply(genelist, uni_cox)
	        do.call(rbind,uni_cox_list)
	  }

	  uni_cox_df_left<-processUniCOX(surv_gene_left,gene_survival_data_left)

	  save(uni_cox_df_left,file = 'data/uni_cox_df_left.Rdata')

	  write.table(uni_cox_df_left,file = 'Tables/Table3.LCC.xls',quote = F,sep='\t')

	  #LASSO analysis: define the function
	  dolasso <- function(uni_cox_df,survival_gene_data) {
		    library(glmnet)
	    survival_gene_data<-survival_gene_data[survival_gene_data$OS.time>0,]
	      x<-as.matrix(survival_gene_data[,uni_cox_df$uni_cox_sig_genes])
	      y<-survival_gene_data[,c('OS.time','OS')]
	        names(y)<-c('time','status')
	        y$time<-as.double(y$time)
		  y$status<-as.double(y$status)
		  y<-as.matrix(survival::Surv(y$time,y$status))
		    lasso_fit<-cv.glmnet(x,y,family='cox',type.measure = 'deviance')
		    coefficient<-coef(lasso_fit,s=lasso_fit$lambda.min)	
		      return(coefficient)
	  }


	  coefficient_left<-dolasso(uni_cox_df_left,gene_survival_data_left)
	  active.Index<-which(as.numeric(coefficient_left)!=0)
	  active.coefficients<-as.numeric(coefficient_left)[active.Index]
	  sig_gene_multi_cox_left<-rownames(coefficient_left)[active.Index]

	  #multivariate Cox regression analysis: define the function
	  doPH_hypothesis<-function(sig_gene_multi_cox, survival_gene_data){
		    formula_for_multivariate<-as.formula(paste0('Surv(OS.time, OS) ~ ', paste(sig_gene_multi_cox,sep = '',collapse = '+')))
	    multi_variate_cox<-coxph(formula_for_multivariate,data=survival_gene_data)
	      #check if variances are supported by PH hypothesis
	      ph_hypothesis_multi<-cox.zph(multi_variate_cox)
	      #the last row of the table records the test results on the Global model, delete it.
	      ph_hypo_table<-ph_hypothesis_multi$table[-nrow(ph_hypothesis_multi$table),]	
	        return(ph_hypo_table)
	        
	  }

	  processMultiCOX <- function(sig_gene_multi_cox, survival_gene_data) {
		    ph_hypo_table=doPH_hypothesis(sig_gene_multi_cox, survival_gene_data)	
	    #remove variances not supported by ph hypothesis and perform the 2nd regression
	    formula_for_multivariate_2<-as.formula(paste0('Surv(OS.time, OS) ~ ', paste(rownames(ph_hypo_table)[ph_hypo_table[,3]>0.05],sep = '',collapse = '+')))
	      multi_variate_cox_2<-coxph(formula_for_multivariate_2,data = survival_gene_data)
	      return(multi_variate_cox_2)
	  }

	  multi_variate_cox_left<-processMultiCOX(sig_gene_multi_cox_left, gene_survival_data_left)

	  C_index_left<-multi_variate_cox_left$concordance['concordance']

	  save(sig_gene_multi_cox_left,multi_variate_cox_left,file='data/model_left_TCGA.Rdata')

	  #forest plot:Figure 2A
	  ggforest(model = multi_variate_cox_left,data = gene_survival_data_left,main = 'Hazard ratios of candidate gene',fontsize = 1)
	  ggsave(filename = 'Figures/Figure 2A.multi_cox_forestplot_left.png',width = 7,height = 4)
	  dev.off()

	  #Figure 2B
	  library(ggrisk)
	  tiff(file='Figures/Figure2B.heatmap_risk.tiff',width=15,height=15,units="cm",compression="lzw",bg="white",res=600)
	  p=ggrisk(multi_variate_cox_left,cutoff.value = 'roc',
		     color.A=c(low='blue',high='red'),
		       color.B=c(code.0='#238443',code.1='#D7301F'), 
		       color.C=c(low='green',median='white',high='red'), 
		       )
	  print(p)
	  dev.off()

	  #calculate the risk score of each patient:define the function
	  riskscore<-function(survival_gene_data,candidate_genes_for_cox,cox_report){
		    suppressMessages(library(dplyr))
	    samplenames<-survival_gene_data$sample %>% na.omit()
	      
	      survival_gene_data<-survival_gene_data[survival_gene_data$sample %in% samplenames,]
	      rownames(survival_gene_data)<-survival_gene_data$sample
	        survival_gene_data<-survival_gene_data[samplenames,]
	        
	        risk_score_table<-survival_gene_data[samplenames,candidate_genes_for_cox]
		  
		  rownames(risk_score_table)<-survival_gene_data$sample 
		  
		  all_riskscore=list()
		    for(each_sample in rownames(risk_score_table)){
			        each_total_score=0
		      for(each_sig_gene in colnames(risk_score_table)){
			            
			            each_total_score=each_total_score+risk_score_table[each_sample,each_sig_gene]*(summary(cox_report)$coefficients[each_sig_gene,1])
		            
		          }
		          all_riskscore=append(all_riskscore,each_total_score)
		          
		        }
		    risk_score_table$total_risk_score<-as.numeric(all_riskscore)
		      risk_score_table<-cbind(survival_gene_data[,c('sample','OS.time','OS')],risk_score_table)
		      #risk_score_table<-cbind(risk_score_table,'total_risk_score'=exp(rowSums(risk_score_table))) %>%
		      #  cbind(survival_gene_data[,c('title','OS.time','OS')])
		      risk_score_table<-risk_score_table[,c('sample','OS.time','OS',candidate_genes_for_cox,'total_risk_score')]
		        return(risk_score_table)
	  }

	  risk_score_table_multi_cox_left<-riskscore(gene_survival_data_left,sig_gene_multi_cox_left,multi_variate_cox_left)
	  risk_score_table_multi_cox_left<-risk_score_table_multi_cox_left[risk_score_table_multi_cox_left$OS.time>0,]

	  #save(risk_score_table_multi_cox_left,file='data/risk_score_table_multi_cox_left.Rdata')
	  #draw the ROC curves: define the function
	  multi_ROC<-function(time_vector,risk_score_table){
		    library(survivalROC)
	    single_ROC<-function(single_time){
		        for_ROC<-survivalROC(Stime = risk_score_table$OS.time,
					                              status=risk_score_table$OS,
								                               marker=risk_score_table$total_risk_score,
								                               predict.time = single_time,method = 'KM')
	        data.frame(
			         'True_positive'=for_ROC$TP,'False_positive'=for_ROC$FP,
				       'Cut_values'=for_ROC$cut.values,
				       'Time_point'=rep(single_time,length(for_ROC$TP)),
				             'AUC'=rep(for_ROC$AUC,length(for_ROC$TP))
				           )
		    
		  }
	      multi_ROC_list<-lapply(time_vector, single_ROC)
	      do.call(rbind,multi_ROC_list)
	  }
	  for_multi_ROC_left<-multi_ROC(time_vector = c(365*seq(3,5,1)),risk_score_table = risk_score_table_multi_cox_left)

	  #visualization for the ROC curves of multiple time points
	  #plot ROC
	  year3ROC<-for_multi_ROC_left[for_multi_ROC_left$Time_point==3*365,]
	  year4ROC<-for_multi_ROC_left[for_multi_ROC_left$Time_point==4*365,]
	  year5ROC<-for_multi_ROC_left[for_multi_ROC_left$Time_point==5*365,]
	  TP_3year = year3ROC$True_positive
	  FP_3year = year3ROC$False_positive
	  TP_4year = year4ROC$True_positive
	  FP_4year = year4ROC$False_positive
	  TP_5year = year5ROC$True_positive
	  FP_5year = year5ROC$False_positive
	  time_ROC_df=data.frame(TP_3year,FP_3year,TP_4year,FP_4year,TP_5year,FP_5year)
	  library(ggplot2)
	  tiff(file='Figures/Figure2D.OS_AUC_left.tiff',width=15,height=15,units="cm",compression="lzw",bg="white",res=600)
	  p=ggplot(data = time_ROC_df) +
		    geom_line(aes(x = FP_3year, y = TP_3year), size = 1, color = "#BC3C29FF") +
		      geom_line(aes(x = FP_4year, y = TP_4year), size = 1, color = "#0072B5FF") +
		        geom_line(aes(x = FP_5year, y = TP_5year), size = 1, color = "#E18727FF") +
			  geom_abline(slope = 1, intercept = 0, color = "grey", size = 1, linetype = 2) +
			    theme_bw() +
			      annotate("text",
				                  x = 0.75, y = 0.25, size = 4.5,
						             label = paste0("AUC at 3 years = ", sprintf("%.3f", year3ROC$AUC)), color = "#BC3C29FF"
						    ) +
  annotate("text",
	              x = 0.75, y = 0.15, size = 4.5,
		                 label = paste0("AUC at 4 years = ", sprintf("%.3f", year4ROC$AUC)), color = "#0072B5FF"
		        ) +
  annotate("text",
	              x = 0.75, y = 0.05, size = 4.5,
		                 label = paste0("AUC at 5 years = ", sprintf("%.3f", year5ROC$AUC)), color = "#E18727FF"
		        ) +
  labs(x = "False positive rate", y = "True positive rate") +
    theme(
	      axis.text = element_text(face = "bold", size = 11, color = "black"),
	          axis.title.x = element_text(face = "bold", size = 14, color = "black", margin = margin(c(15, 0, 0, 0))),
	          axis.title.y = element_text(face = "bold", size = 14, color = "black", margin = margin(c(0, 15, 0, 0)))
		    )

  print(p)
  dev.off()

  #K-M analysis
  AUC_max<-max(for_multi_ROC_left$AUC)
  AUC_max_time<-for_multi_ROC_left$Time_point[which(for_multi_ROC_left$AUC==AUC_max)]
  AUC_max_time<-AUC_max_time[!duplicated(AUC_max_time)]
  AUC_max_time<-AUC_max_time[length(AUC_max_time)]
  for_multi_ROC_left$Time_point<-as.factor(for_multi_ROC_left$Time_point)
  optimal_time_ROC_df<-for_multi_ROC_left[which(for_multi_ROC_left$Time_point==AUC_max_time),]
  cut.off<-optimal_time_ROC_df$Cut_values[which.max(optimal_time_ROC_df$True_positive - optimal_time_ROC_df$False_positive)]
  high_low<-risk_score_table_multi_cox_left$total_risk_score>cut.off
  high_low[high_low==TRUE]<-'high'
  high_low[high_low==FALSE]<-'low'
  risk_score_table_multi_forKM<-cbind(risk_score_table_multi_cox_left,high_low)

  risk_score_table_multi_forKM_left<-risk_score_table_multi_forKM
  #KM_plot 
  library(survminer)
  risk_score_table_multi_forKM$OS[which(risk_score_table_multi_forKM$OS.time>AUC_max_time)]<-0
  risk_score_table_multi_forKM$OS.time[which(risk_score_table_multi_forKM$OS.time>AUC_max_time)]<-AUC_max_time
  fit_km<-survfit(Surv(OS.time,OS)~high_low,data = risk_score_table_multi_forKM)
  tiff(file='Figures/Figure2C.OS_riskscore_left.tiff',width=15,height=15,units="cm",compression="lzw",bg="white",res=300)
  p=ggsurvplot(fit_km,conf.int = F,pval = T,legend.title="Total risk score",
	                    legend.labs=c('high','low'),risk.table=T,
			                 #legend.labs=c(paste0('>',as.character(round(cut.off,2))),
			                 # paste0('<=',as.character(round(cut.off,2)))),risk.table=T,
			                 palette=c('red','blue'),surv.median.line='hv')
  print(p)
  dev.off()

  #prognosis for RCC
  TCGA_pheno_survival_right<-TCGA_pheno_survival[TCGA_pheno_survival$site=='Right',]
  group_list<-ifelse(substr(TCGA_pheno_survival_right$sample,14,15)<10,'Tumor','Normal')
  commone_sample_right<-intersect(TCGA_pheno_survival_right$sample,colnames(TCGA_colon_expr_filt))
  tumor_sample_right<-commone_sample_right[as.numeric(substr(commone_sample_right,14,15))<10]

  common_gene_right<-intersect(all_diffgene_u,rownames(TCGA_colon_expr_filt))
  gene_survival_data_right<-cbind(TCGA_pheno_survival[tumor_sample_right,c('sample','OS','OS.time')],t(TCGA_colon_expr_filt[common_gene_right,tumor_sample_right]))

  save(gene_survival_data_right,file = 'data/gene_survival_data_right.Rdata')

  load('data/gene_survival_data_right.Rdata')

  library(survminer)
  library(survival)
  surv_gene_right<-c()

  dir.create('Figures/KM_right/')
  for (i in 4:(ncol(gene_survival_data_right)-1)){
	    gene_survival_data_right$group<-ifelse(gene_survival_data_right[,i]>median(gene_survival_data_right[,i]),'high','low')
    fit <- survfit(Surv(as.numeric(OS.time), as.numeric(OS)) ~gene_survival_data_right$group, data = gene_survival_data_right)

      a<-as.character(surv_pvalue(fit)[4])

      if(as.numeric(gsub('p = ','',a))<0.05){
	         surv_gene_right<-append(surv_gene_right,colnames(gene_survival_data_right)[i])
         # write.table(paste(colnames(gene_survival_data_right)[i],as.character(surv_pvalue(fit)[4]),sep = '\t'),file = survfile,append = T,quote = F,sep = '\t',row.names = F,col.names = F)
          p=ggsurvplot(fit, 
		                        #pval = TRUE, 
		                        pval = toupper(surv_pvalue(fit)[4]),
					                 conf.int = TRUE,
					                 risk.table = TRUE, # Add risk table
							                  risk.table.col = "strata", # Change risk table color by groups
							                  linetype = "strata", # Change line type by groups
									                   surv.median.line = "hv", # Specify median survival
									                   ggtheme = theme_bw(), # Change ggplot2 theme
											                    legend.labs = c("High", "Low"), 
											                    legend.title = paste0("Expression of ",colnames(gene_survival_data_right)[i]), 
													                     palette = c( "#f03b20","#2c7fb8"))
	      
	      filename = paste0('Figures/KM_right/',colnames(gene_survival_data_right)[i],'.tiff')
	      tiff(file=filename,width=15,height=15,units="cm",compression="lzw",bg="white",res=600)
	          print(p)
	          dev.off()
		    }
  }

  surv_gene_right<-gsub(surv_gene_right,pattern = '-',replacement = '_')
  colnames(gene_survival_data_right)<-gsub(colnames(gene_survival_data_right),pattern = '-',replacement = '_')

  uni_cox_df_right<-processUniCOX(surv_gene_right,gene_survival_data_right)

  save(uni_cox_df_right,file = 'data/uni_cox_df_right.Rdata')
  #Table 3
  write.table(uni_cox_df_right,file = 'COX_right/table_uni_cox_right.xls',quote = F,sep='\t')

  #LASSO
  set.seed(1234)
  coefficient_right<-dolasso(uni_cox_df_right,gene_survival_data_right)
  active.Index<-which(as.numeric(coefficient_right)!=0)
  active.coefficients<-as.numeric(coefficient_right)[active.Index]
  sig_gene_multi_cox_right<-rownames(coefficient_right)[active.Index]

  multi_variate_cox_right<-processMultiCOX(sig_gene_multi_cox_right, gene_survival_data_right)

  C_index<-multi_variate_cox_right$concordance['concordance']

  #save(sig_gene_multi_cox_right,multi_variate_cox_right,file='data/model_right_TCGA.Rdata')

  #load('data/model_right_TCGA.Rdata')
  #load('data/gene_survival_data_right.Rdata')

  #forest plot：Figure 2E
  ggforest(model = multi_variate_cox_right,data = gene_survival_data_right,main = 'Hazard ratios of candidate gene',fontsize = 1)
  ggsave(filename = 'Figures/Figure2E.multi_cox_forestplot_right.png',width = 10,height = 4)
  dev.off()

  risk_score_table_multi_cox_right<-riskscore(gene_survival_data_right,sig_gene_multi_cox_right,multi_variate_cox_right)
  risk_score_table_multi_cox_right<-risk_score_table_multi_cox_right[risk_score_table_multi_cox_right$OS.time>0,]

  #evaluate AUCs between 3-5 years.
  for_multi_ROC_right<-multi_ROC(time_vector = c(365*seq(3,5,1)),risk_score_table = risk_score_table_multi_cox_right)

  #visualization for the ROC curves of multiple time points
  #plot ROC
  year3ROC<-for_multi_ROC_right[for_multi_ROC_right$Time_point==3*365,]
  year4ROC<-for_multi_ROC_right[for_multi_ROC_right$Time_point==4*365,]
  year5ROC<-for_multi_ROC_right[for_multi_ROC_right$Time_point==5*365,]
  TP_3year = year3ROC$True_positive
  FP_3year = year3ROC$False_positive
  TP_4year = year4ROC$True_positive
  FP_4year = year4ROC$False_positive
  TP_5year = year5ROC$True_positive
  FP_5year = year5ROC$False_positive
  time_ROC_df=data.frame(TP_3year,FP_3year,TP_4year,FP_4year,TP_5year,FP_5year)
  library(ggplot2)
  tiff(file='Figures/Figure2H.OS_AUC_right.tiff',width=15,height=15,units="cm",compression="lzw",bg="white",res=600)
  p=ggplot(data = time_ROC_df) +
	    geom_line(aes(x = FP_3year, y = TP_3year), size = 1, color = "#BC3C29FF") +
	      geom_line(aes(x = FP_4year, y = TP_4year), size = 1, color = "#0072B5FF") +
	        geom_line(aes(x = FP_5year, y = TP_5year), size = 1, color = "#E18727FF") +
		  geom_abline(slope = 1, intercept = 0, color = "grey", size = 1, linetype = 2) +
		    theme_bw() +
		      annotate("text",
			                  x = 0.75, y = 0.25, size = 4.5,
					             label = paste0("AUC at 3 years = ", sprintf("%.3f", year3ROC$AUC)), color = "#BC3C29FF"
					    ) +
  annotate("text",
	              x = 0.75, y = 0.15, size = 4.5,
		                 label = paste0("AUC at 4 years = ", sprintf("%.3f", year4ROC$AUC)), color = "#0072B5FF"
		        ) +
  annotate("text",
	              x = 0.75, y = 0.05, size = 4.5,
		                 label = paste0("AUC at 5 years = ", sprintf("%.3f", year5ROC$AUC)), color = "#E18727FF"
		        ) +
  labs(x = "False positive rate", y = "True positive rate") +
    theme(
	      axis.text = element_text(face = "bold", size = 11, color = "black"),
	          axis.title.x = element_text(face = "bold", size = 14, color = "black", margin = margin(c(15, 0, 0, 0))),
	          axis.title.y = element_text(face = "bold", size = 14, color = "black", margin = margin(c(0, 15, 0, 0)))
		    )

  print(p)
  dev.off()

  #K-M analysis
  AUC_max<-max(for_multi_ROC_right$AUC)
  AUC_max_time<-for_multi_ROC_right$Time_point[which(for_multi_ROC_right$AUC==AUC_max)]
  AUC_max_time<-AUC_max_time[!duplicated(AUC_max_time)]
  AUC_max_time<-AUC_max_time[length(AUC_max_time)]
  for_multi_ROC_right$Time_point<-as.factor(for_multi_ROC_right$Time_point)
  optimal_time_ROC_df<-for_multi_ROC_right[which(for_multi_ROC_right$Time_point==AUC_max_time),]
  cut.off<-optimal_time_ROC_df$Cut_values[which.max(optimal_time_ROC_df$True_positive - optimal_time_ROC_df$False_positive)]

  high_low<-risk_score_table_multi_cox_right$total_risk_score>cut.off
  high_low[high_low==TRUE]<-'high'
  high_low[high_low==FALSE]<-'low'
  risk_score_table_multi_forKM_right<-cbind(risk_score_table_multi_cox_right,high_low)
  #KM_plot 
  library(survminer)
  risk_score_table_multi_forKM_right$OS[which(risk_score_table_multi_forKM_right$OS.time>AUC_max_time)]<-0
  risk_score_table_multi_forKM_right$OS.time[which(risk_score_table_multi_forKM_right$OS.time>AUC_max_time)]<-AUC_max_time
  fit_km<-survfit(Surv(OS.time,OS)~high_low,data = risk_score_table_multi_forKM_right)
  tiff(file='Figures/Figure2G.OS_riskscore_right.tiff',width=15,height=15,units="cm",compression="lzw",bg="white",res=300)
  p=ggsurvplot(fit_km,conf.int = F,pval = T,legend.title="Total risk score",
	                    legend.labs=c('high','low'),risk.table=T,
			                 #legend.labs=c(paste0('>',as.character(round(cut.off,2))),
			                 # paste0('<=',as.character(round(cut.off,2)))),risk.table=T,
			                 palette=c('red','blue'),surv.median.line='hv')
  print(p)
  dev.off()

  #add microsatellite status
  clinical_data<-read.table('data/COAD_clinicalMatrix.xls',header = T,sep='\t')
  clinical_data$MSI<-ifelse(clinical_data$MSI_updated_Oct6201=="",clinical_data$CDE_ID_3226963,clinical_data$MSI_updated_Oct62011)
  clinical_data$sample<-clinical_data$sampleID
  TCGA_pheno_survival_right<-merge(TCGA_pheno_survival_right,clinical_data,by='sample')
  TCGA_pheno_survival_right<-TCGA_pheno_survival_right[,c('sample','age_at_initial_pathologic_diagnosis.x','gender.demographic','tumor_stage.diagnoses','pathologic_T.x','pathologic_M.x','pathologic_N.x','MSI')] 
  colnames(TCGA_pheno_survival_right)=c('sampleid','Age','Gender','Stage','T','M','N','microsatellite_status')
  rownames(TCGA_pheno_survival_right)<-TCGA_pheno_survival_right$sampleid
  library(stringr)
  TCGA_pheno_survival_right <- TCGA_pheno_survival_right %>% na.omit()
  TCGA_pheno_survival_right$Age<-factor(ifelse(TCGA_pheno_survival_right$Age >60,'>60','<=60'))
  TCGA_pheno_survival_right$Group<-ifelse(substr(rownames(TCGA_pheno_survival_right),14,15)<10,'Tumor','Normal')
  TCGA_pheno_survival_right$Stage<-factor(toupper(str_extract(TCGA_pheno_survival_right$Stage, "i+v*")))
  TCGA_pheno_survival_right$T<-factor(str_extract(TCGA_pheno_survival_right$T, "T[0-9]+"))
  TCGA_pheno_survival_right$M <- factor(str_extract(TCGA_pheno_survival_right$M, "M[0-9]*X*"))
  TCGA_pheno_survival_right$N<-factor(str_extract(TCGA_pheno_survival_right$N, "N[0-9]"))

  TCGA_pheno_survival_right$Group<-factor(TCGA_pheno_survival_right$Group,
					                                          levels=c('Tumor','Normal'),
										                                          labels=c("Tumor", # 第一个作为参考组
																                                                    "Normal"))

  save(TCGA_pheno_survival_right,file = 'data/TCGA_pheno_survival_right.Rdata')
  load('data/TCGA_pheno_survival_right.Rdata')

  library(table1)
  table1(~ factor(Gender) + factor(Age) +factor(microsatellite_status)+ factor(Stage) +
	          factor(T) +  factor(M) +  factor(N) | Group, data=TCGA_pheno_survival_right)

  TCGA_pheno_survival_right$sample<-TCGA_pheno_survival_right$sampleid
  data_right<-merge(TCGA_pheno_survival_right, risk_score_table_multi_cox_right, by = "sample")
  save(risk_score_table_multi_cox_right,file="data/risk_score_table_multi_cox_right.Rdata")

  save(data_right,file='data/data_right.Rdata')
  library(cowplot)
  library(tidyverse)
  library(ggplot2)
  library(ggsci)
  library(ggpubr)

  #load('data/data_right.Rdata')
  #for T stage
  library(dplyr)
  T_data<-data.frame(data_right$sample,data_right$T,data_right$total_risk_score) %>% na.omit()#让分数变为正数
  colnames(T_data)<-c('sample','T stage','level of risk score')
  library(tidyr)
  rownames(T_data)<-T_data$sample
  T_data$sample<-NULL
  T_data$`T stage`<-as.factor(T_data$`T stage`)
  mydata<-T_data %>% gather(key = "T stage",value = "level of risk score") %>% na.omit() 
  my_comparisons <- list( c("T1", "T2"),c("T1", "T3"),c("T1", "T4"),c("T2", "T3"),c("T2", "T4"),c("T3", "T4"))
  tiff(file='Figures/Figure3E.T_stage_score.tiff',width=18,height=15,units="cm",compression="lzw",bg="white",res=300)
  p<-ggboxplot(mydata,x='T stage',y='level of risk score',color='T stage',palette = 'jama',add = 'jitter')
  p+stat_compare_means()
  dev.off()


  #for M stage
  M_data<-data.frame(data_right$sample,data_right$M,data_right$total_risk_score)  %>% na.omit() #让分数变为正数
  colnames(M_data)<-c('sample','M stage','level of risk score')
  M_data<-M_data[M_data$`M stage`!='MX',]
  rownames(M_data)<-M_data$sample
  M_data$sample<-NULL
  M_data$`M stage`<-as.factor(M_data$`M stage`)
  mydata<-M_data %>% gather(key = "M stage",value = "level of risk score") %>% na.omit() 
  my_comparisons <- list( c("M0", "M1"))
  tiff(file='Figures/Figure3A.M_stage_score2.tiff',width=18,height=15,units="cm",compression="lzw",bg="white",res=300)
  p<-ggboxplot(mydata,x='M stage',y='level of risk score',color='M stage',palette = 'jama',add = 'jitter')
  #p+stat_compare_means(comparisons = my_comparisons)
  p+stat_compare_means()
  dev.off()

  #for N stage
  N_data<-data.frame(data_right$sample,data_right$N,data_right$total_risk_score) %>% na.omit()
  colnames(N_data)<-c('sample','N stage','level of risk score')
  rownames(N_data)<-N_data$sample
  N_data$sample<-NULL
  N_data$`N stage`<-as.factor(N_data$`N stage`)
  mydata<-N_data %>% gather(key = "N stage",value = "level of risk score") %>% na.omit() 
  my_comparisons <- list( c("N0", "N1"),c("N0", "N2"),c("N1", "N2"))
  tiff(file='Figures/Figure3B.N_stage_score2.tiff',width=18,height=15,units="cm",compression="lzw",bg="white",res=600)
  p<-ggboxplot(mydata,x='N stage',y='level of risk score',color='N stage',palette = 'jama',add = 'jitter')
  #p+stat_compare_means(comparisons = my_comparisons)
  p+stat_compare_means()
  dev.off()

  #for stage
  stage_data<-data.frame(data_right$sample,data_right$Stage,data_right$total_risk_score) %>% na.omit()#让分数变为正数
  colnames(stage_data)<-c('sample','stage','level of risk score')
  rownames(stage_data)<-stage_data$sample
  stage_data$sample<-NULL
  stage_data$`stage`<-as.factor(stage_data$`stage`)
  mydata<-stage_data %>% gather(key = "stage",value = "level of risk score") %>% na.omit() 
  my_comparisons <- list( c("i", "ii"),c("i", "iii"),c("i", "iv"),c("ii", "iii"),c("ii", "iv"),c("iii", "iv"))
  tiff(file='Figures/Figure3C.stage_score.tiff',width=18,height=15,units="cm",compression="lzw",bg="white",res=300)
  p<-ggboxplot(mydata,x='stage',y='level of risk score',color='stage',palette = 'jama',add = 'jitter')
  #p+stat_compare_means(comparisons = my_comparisons)
  p+stat_compare_means()
  dev.off()

  #for microsatellite status
  data_right$microsatellite_status<-ifelse(substr(data_right$microsatellite_status,1,2)!='MS','NA',data_right$microsatellite_status)
  #过滤掉NA
  data_right<-data_right[data_right$microsatellite_status!='NA',]
  #将MSS和MSI-L合并
  data_right$microsatellite_status<-ifelse(data_right$microsatellite_status=='MSI-H','MSI-H','MSI-L/MSS')

  MS_data<-data.frame(data_right$sample,data_right$microsatellite_status,data_right$total_risk_score) %>% na.omit()#让分数变为正数
  colnames(MS_data)<-c('sample','microsatellite status','level of risk score')
  library(tidyr)
  rownames(MS_data)<-MS_data$sample
  MS_data$sample<-NULL
  MS_data$`microsatellite status`<-as.factor(MS_data$`microsatellite status`)
  mydata<-MS_data %>% gather(key = "microsatellite status",value = "level of risk score") %>% na.omit() 
  my_comparisons <- list( c("MSI-L/MSS", "MSI-H"))
  #my_comparisons<-list(c('MSS','MSI-H'),c('MSI-L','MSI-H'))
  tiff(file='Figures/Figure3D.MSS_score.tiff',width=18,height=15,units="cm",compression="lzw",bg="white",res=300)
  p<-ggboxplot(mydata,x='microsatellite status',y='level of risk score',color='microsatellite status',palette = 'jama',add = 'jitter')
  #p+stat_compare_means(comparisons = my_comparisons)
  p+stat_compare_means()
  dev.off()

  #左肠stage的预测
  TCGA_pheno_survival_left<-merge(TCGA_pheno_survival_left,clinical_data,by='sample')
  TCGA_pheno_survival_left<-TCGA_pheno_survival_left[,c('sample','age_at_initial_pathologic_diagnosis.x','gender.demographic','tumor_stage.diagnoses','pathologic_T.x','pathologic_M.x','pathologic_N.x','MSI')] 
  colnames(TCGA_pheno_survival_left)=c('sampleid','Age','Gender','Stage','T','M','N','microsatellite_status')
  rownames(TCGA_pheno_survival_left)<-TCGA_pheno_survival_left$sampleid
  library(stringr)
  TCGA_pheno_survival_left <- TCGA_pheno_survival_left %>% na.omit()
  TCGA_pheno_survival_left$Age<-factor(ifelse(TCGA_pheno_survival_left$Age >60,'>60','<=60'))
  TCGA_pheno_survival_left$Group<-ifelse(substr(rownames(TCGA_pheno_survival_left),14,15)<10,'Tumor','Normal')
  TCGA_pheno_survival_left$Stage<-factor(toupper(str_extract(TCGA_pheno_survival_left$Stage, "i+v*")))
  TCGA_pheno_survival_left$T<-factor(str_extract(TCGA_pheno_survival_left$T, "T[0-9]+"))
  TCGA_pheno_survival_left$M <- factor(str_extract(TCGA_pheno_survival_left$M, "M[0-9]*X*"))
  TCGA_pheno_survival_left$N<-factor(str_extract(TCGA_pheno_survival_left$N, "N[0-9]"))

  TCGA_pheno_survival_left$Group<-factor(TCGA_pheno_survival_left$Group,
					                                        levels=c('Tumor','Normal'),
										                                       labels=c("Tumor", # 第一个作为参考组
																                                                "Normal"))

  save(TCGA_pheno_survival_left,file = 'data/TCGA_pheno_survival_left.Rdata')
  load('data/TCGA_pheno_survival_left.Rdata')
  load('data/risk_score_table_multi_cox_left.Rdata')
  library(table1)
  table1(~ factor(Gender) + factor(Age) +factor(microsatellite_status)+ factor(Stage) +
	          factor(T) +  factor(M) +  factor(N) | Group, data=TCGA_pheno_survival_left)

  TCGA_pheno_survival_left$sample<-TCGA_pheno_survival_left$sampleid
  data_left<-merge(TCGA_pheno_survival_left, risk_score_table_multi_cox_left, by = "sample")

  save(data_left,file='data/data_left.Rdata')
  library(cowplot)
  library(tidyverse)
  library(ggplot2)
  library(ggsci)
  library(ggpubr)

  load('data/data_left.Rdata')
  #for T stage
  T_data<-data.frame(data_left$sample,data_left$T,data_left$total_risk_score) %>% na.omit()#让分数变为正数
  colnames(T_data)<-c('sample','T stage','level of risk score')
  library(tidyr)
  rownames(T_data)<-T_data$sample
  T_data$sample<-NULL
  T_data$`T stage`<-as.factor(T_data$`T stage`)
  mydata<-T_data %>% gather(key = "T stage",value = "level of risk score") %>% na.omit() 

  tiff(file='Figures/FigureS3E.T_stage_score.tiff',width=18,height=15,units="cm",compression="lzw",bg="white",res=300)
  p<-ggboxplot(mydata,x='T stage',y='level of risk score',color='T stage',palette = 'jama',add = 'jitter')
  p+stat_compare_means()
  dev.off()

  #for M stage
  M_data<-data.frame(data_left$sample,data_left$M,data_left$total_risk_score)  %>% na.omit() #让分数变为正数
  colnames(M_data)<-c('sample','M stage','level of risk score')
  M_data<-M_data[M_data$`M stage`!='MX',]
  rownames(M_data)<-M_data$sample
  M_data$sample<-NULL
  M_data$`M stage`<-as.factor(M_data$`M stage`)
  mydata<-M_data %>% gather(key = "M stage",value = "level of risk score") %>% na.omit() 
  my_comparisons <- list( c("M0", "M1"))
  tiff(file='Figures/Figure S3F.M_stage_score.tiff',width=18,height=15,units="cm",compression="lzw",bg="white",res=300)
  p<-ggboxplot(mydata,x='M stage',y='level of risk score',color='M stage',palette = 'jama',add = 'jitter')
  p+stat_compare_means()
  dev.off()

  #for N stage
  N_data<-data.frame(data_left$sample,data_left$N,data_left$total_risk_score) %>% na.omit()
  colnames(N_data)<-c('sample','N stage','level of risk score')
  rownames(N_data)<-N_data$sample
  N_data$sample<-NULL
  N_data$`N stage`<-as.factor(N_data$`N stage`)
  mydata<-N_data %>% gather(key = "N stage",value = "level of risk score") %>% na.omit() 
  my_comparisons <- list( c("N0", "N1"),c("N0", "N2"),c("N1", "N2"))
  tiff(file='Figures/FigureS3G.N_stage_score.tiff',width=18,height=15,units="cm",compression="lzw",bg="white",res=600)
  p<-ggboxplot(mydata,x='N stage',y='level of risk score',color='N stage',palette = 'jama',add = 'jitter')
  p+stat_compare_means()
  dev.off()

  #for stage
  stage_data<-data.frame(data_left$sample,data_left$Stage,data_left$total_risk_score) %>% na.omit()#让分数变为正数
  colnames(stage_data)<-c('sample','stage','level of risk score')
  rownames(stage_data)<-stage_data$sample
  stage_data$sample<-NULL
  stage_data$`stage`<-as.factor(stage_data$`stage`)
  mydata<-stage_data %>% gather(key = "stage",value = "level of risk score") %>% na.omit() 
  my_comparisons <- list( c("i", "ii"),c("i", "iii"),c("i", "iv"),c("ii", "iii"),c("ii", "iv"),c("iii", "iv"))
  tiff(file='Figures/Figure_S3H.stage_score2.tiff',width=18,height=15,units="cm",compression="lzw",bg="white",res=300)
  p<-ggboxplot(mydata,x='stage',y='level of risk score',color='stage',palette = 'jama',add = 'jitter')
  p+stat_compare_means()
  dev.off()

  #for microsatellite status
  data_left$microsatellite_status<-ifelse(substr(data_left$microsatellite_status,1,2)!='MS','NA',data_left$microsatellite_status)

  data_left<-data_left[data_left$microsatellite_status!='NA',]
  data_left$microsatellite_status<-ifelse(data_left$microsatellite_status=='MSI-H','MSI-H','MSI-L/MSS')
  MS_data<-data.frame(data_left$sample,data_left$microsatellite_status,data_left$total_risk_score) %>% na.omit()#让分数变为正数
  colnames(MS_data)<-c('sample','microsatellite status','level of risk score')
  library(tidyr)
  rownames(MS_data)<-MS_data$sample
  MS_data$sample<-NULL
  MS_data$`microsatellite status`<-as.factor(MS_data$`microsatellite status`)
  mydata<-MS_data %>% gather(key = "microsatellite status",value = "level of risk score") %>% na.omit() 
  my_comparisons <- list( c("MSI-L/MSS", "MSI-H"))
  tiff(file='Figures/Figure_S3I.MSS_score2.tiff',width=18,height=15,units="cm",compression="lzw",bg="white",res=300)
  p<-ggboxplot(mydata,x='microsatellite status',y='level of risk score',color='microsatellite status',palette = 'jama',add = 'jitter')
  p+stat_compare_means(comparisons = my_comparisons)
  dev.off()

  #the expression level of marker genes in TCGA-COAD cohort
  #load('data/gene_survival_data_left.Rdata')
  #load('data/gene_survival_data_right.Rdata')
  #load('data/data_left.Rdata')
  #load('data/data_right.Rdata')
  #load('data/TCGA_colon_expr_filt.Rdata')
  #load('data/TCGA_pheno_survival.Rdata')
  TCGA_markers<-c(colnames(data_left)[13:(ncol(data_left)-1)],colnames(data_right)[13:(ncol(data_right)-1)])

  common_sample<-intersect(colnames(TCGA_colon_expr_filt),substr(TCGA_pheno_survival$sample,1,15))
  tumor_sample<-common_sample[substr(common_sample,14,15)<10]
  TCGA_colon_expr_target<-t(TCGA_colon_expr_filt[TCGA_markers,tumor_sample])
  TCGA_colon_expr_target<- data.frame(rownames(TCGA_colon_expr_target),TCGA_colon_expr_target)
  colnames(TCGA_colon_expr_target)[1]<-'sample'
  TCGA_pheno_survival$sample<-substr(TCGA_pheno_survival$sample,1,15)

  TCGA_colon_expr_target_geno<-merge(TCGA_pheno_survival,TCGA_colon_expr_target,by='sample')
  TCGA_colon_expr_target_geno_left<-TCGA_colon_expr_target_geno[TCGA_colon_expr_target_geno$site=='Left',]
  TCGA_colon_expr_target_geno_right<-TCGA_colon_expr_target_geno[TCGA_colon_expr_target_geno$site=='Right',]

  marker_in=array()
  for (marker in TCGA_markers){
	    
	    uni_marker_left<-data.frame(marker,TCGA_colon_expr_target_geno_left[,marker],TCGA_colon_expr_target_geno_left$sample,'Left')
    uni_marker_right<-data.frame(marker,TCGA_colon_expr_target_geno_right[,marker],TCGA_colon_expr_target_geno_right$sample,'Right')
      colnames(uni_marker_left)<-c('Markers','expression','sample','Site')
      colnames(uni_marker_right)<-c('Markers','expression','sample','Site')
        
        uni_marker_in<-rbind(uni_marker_left,uni_marker_right)
        
        marker_in<-rbind(marker_in,uni_marker_in)
  }

  library(dplyr)
  marker_in<-marker_in %>% na.omit()


  library(cowplot)
  library(tidyverse)
  library(ggplot2)
  library(ggsci)
  library(ggpubr)
  library(dplyr)
  marker_in$Group<-marker_in$Site
  marker_in$Group<-ifelse(marker_in$Group=='Left','LCC','RCC')

  p <- ggboxplot(marker_in, x = "Markers", y = "expression",
		                color = 'Group', palette = "nejm",add = "jitter")+
  stat_compare_means(aes(group = Group),method = "t.test",label = "p.format")+rotate_x_text(angle = 45)
tiff(file='Figures/Figure3H.left_right_marker.tiff',width=25,height=15,units="cm",compression="lzw",bg="white",res=600)
print(p)
dev.off()

#predict the CMS of TCGA-COAD
load('data/data_left.Rdata')
load('data/TCGA_colon_expr_filt.Rdata')
TCGA_COAD<-read.table('data/TCGA-COAD.htseq_counts.tsv',sep='\t',header = T)
TCGA_COAD$Ensembl_ID<-substr(TCGA_COAD$Ensembl_ID,1,15)

library(org.Hs.eg.db)
library(clusterProfiler)
ids<-TCGA_COAD$Ensembl_ID
addID<-bitr(ids, fromType="ENSEMBL", toType=c("SYMBOL","ENTREZID"), OrgDb="org.Hs.eg.db")
colnames(addID)[1]<-colnames(TCGA_COAD)[1]
TCGA_COAD_use<-merge(addID,TCGA_COAD,by='Ensembl_ID')
TCGA_COAD_use<-distinct(TCGA_COAD_use,ENTREZID,.keep_all = TRUE)
rownames(TCGA_COAD_use)<-TCGA_COAD_use$ENTREZID
TCGA_COAD_use$Ensembl_ID<-NULL
save(TCGA_COAD_use,file='data/TCGA_COAD_use.Rdata')
load('TCGA_COAD_use.Rdata')
TCGA_COAD_use$SYMBOL<-NULL
TCGA_COAD_use$ENTREZID<-NULL
TCGA_COAD_use<-TCGA_COAD_use[,substr(colnames(TCGA_COAD_use),14,15)<10]
#TCGA CMS预测
library(Biobase)
library(CMScaller)
emat<-exprs(crcTCGAsubset)
par(mfrow=c(1,3))
res<-CMScaller(TCGA_COAD_use,RNAseq = TRUE, doPlot = FALSE)
#cam <- CMSgsa(TCGA_COAD_use, class=res$prediction, RNAseq=TRUE)
save(res,file = 'data/TCGA_COAD_CMS.Rdata')


left_markers<-c("FOSB","RPL35","REG1A", "TESC","C11orf96")
right_markers<-c("HSPA1A","CD69","GDF15","LGALS2")

load('data/TCGA_COAD_CMS.Rdata')
load('data/TCGA_pheno_survival.Rdata')

res$sample<-gsub(pattern = '\\.',replacement = '-',rownames(res))
res$sample<-substr(res$sample,1,15)
TCGA_pheno_survival_CMS<-merge(res,TCGA_pheno_survival,by='sample')

TCGA_pheno_CMS<-data.frame(TCGA_pheno_survival_CMS$sample,TCGA_pheno_survival_CMS$site,TCGA_pheno_survival_CMS$prediction)

colnames(TCGA_pheno_CMS)<-c('sample','tumor_site','CMS_group')
TCGA_pheno_CMS<-TCGA_pheno_CMS[TCGA_pheno_CMS$tumor_site!='Transverse',]
TCGA_pheno_CMS$tumor_site<-ifelse(TCGA_pheno_CMS$tumor_site=='Left','LCC','RCC')
TCGA_pheno_CMS<-TCGA_pheno_CMS %>% na.omit()

save(TCGA_pheno_CMS,file='data/TCGA_pheno_CMS.Rdata')

load('data/TCGA_pheno_CMS.Rdata')

TCGA_COAD_markers<-TCGA_COAD_use
TCGA_COAD_markers$ENTREZID<-rownames(TCGA_COAD_markers)
TCGA_COAD_markers_ID<-merge(addID,TCGA_COAD_markers,by='ENTREZID')
TCGA_COAD_markers_ID$ENTREZID<-NULL
TCGA_COAD_markers_ID$Ensembl_ID<-NULL

TCGA_pheno_CMS_markers<-TCGA_COAD_markers_ID[TCGA_COAD_markers_ID$SYMBOL %in% c(left_markers,right_markers),]
rownames(TCGA_pheno_CMS_markers)<-TCGA_pheno_CMS_markers$SYMBOL
TCGA_pheno_CMS_markers$SYMBOL<-NULL
TCGA_pheno_CMS_markers_t<-t(TCGA_pheno_CMS_markers)
rownames(TCGA_pheno_CMS_markers_t)<-substr(rownames(TCGA_pheno_CMS_markers_t),1,15)
TCGA_pheno_CMS_markers_t<-as.data.frame(TCGA_pheno_CMS_markers_t)
TCGA_pheno_CMS_markers_t$sample<-gsub(pattern = '\\.',replacement = '-',rownames(TCGA_pheno_CMS_markers_t))
rownames(TCGA_pheno_CMS_markers_t)<-TCGA_pheno_CMS_markers_t$sample

TCGA_pheno_CMS_markers_use<-merge(TCGA_pheno_CMS,TCGA_pheno_CMS_markers_t,by='sample')
TCGA_pheno_CMS_markers_use_OS<-merge(TCGA_pheno_CMS_markers_use,TCGA_pheno_survival[,c('sample','OS','OS.time')],by='sample')
#select the tumor samples
#TCGA_pheno_CMS_markers_use_OS<-TCGA_pheno_CMS_markers_use_OS[as.numeric(substr(TCGA_pheno_CMS_markers_use_OS$sample,14,15))<10,]

save(TCGA_pheno_CMS_markers_use_OS,file='data/TCGA_pheno_CMS_markers_use_OS.Rdata')

load('data/TCGA_pheno_CMS_markers_use_OS.Rdata')

library(cowplot)
library(tidyverse)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(tidyr)

TCGA_pheno_CMS_markers_use_OS_use<-data.frame(TCGA_pheno_CMS_markers_use_OS$sample,TCGA_pheno_CMS_markers_use_OS$OS.time,TCGA_pheno_CMS_markers_use_OS$CMS_group)
colnames(TCGA_pheno_CMS_markers_use_OS_use)<-c('sample','OS.time','CMS_group')

library(tidyr)
TCGA_pheno_CMS_markers_use_OS_use<-distinct(TCGA_pheno_CMS_markers_use_OS_use,sample,.keep_all = TRUE)
rownames(TCGA_pheno_CMS_markers_use_OS_use)<-TCGA_pheno_CMS_markers_use_OS_use$sample

TCGA_pheno_CMS_markers_use_OS_use$sample<-NULL


TCGA_pheno_CMS_markers_use_OS_use$CMS_group<-as.factor(TCGA_pheno_CMS_markers_use_OS_use$CMS_group)

mydata<-TCGA_pheno_CMS_markers_use_OS_use %>% gather(key = "CMS_group",value = "OS.time") %>% na.omit() 
my_comparisons <- list( c("CMS1", "CMS2"),c("CMS1", "CMS3"),c("CMS1", "CMS4"),c("CMS2", "CMS3"),c("CMS2", "CMS4"),c("CMS3", "CMS4"))

tiff(file='Figures/FigureS1C.CMS_OS_TCGA.tiff',width=18,height=15,units="cm",compression="lzw",bg="white",res=600)
p<-ggboxplot(mydata,x='CMS_group',y='OS.time',color='CMS_group',palette = 'jama',add = 'jitter')
p+stat_compare_means(comparisons = my_comparisons)
#p+stat_compare_means()
dev.off()

#the expression differences
TCGA_pheno_CMS_markers_use_OS<-distinct(TCGA_pheno_CMS_markers_use_OS,sample,.keep_all = TRUE)

rownames(TCGA_pheno_CMS_markers_use_OS)<-TCGA_pheno_CMS_markers_use_OS$sample
marker_in=array()
for (marker in c(left_markers,right_markers)){
	  for (sample in rownames(TCGA_pheno_CMS_markers_use_OS)){
		      uni_marker_CMS<-data.frame(marker,TCGA_pheno_CMS_markers_use_OS[sample,marker],sample,TCGA_pheno_CMS_markers_use_OS[sample,'CMS_group'])
    colnames(uni_marker_CMS)<-c('Markers','expression','sample','CMS_group')
        marker_in<-rbind(marker_in,uni_marker_CMS)
      }
}

library(dplyr)
marker_in<-marker_in %>% na.omit()


library(cowplot)
library(tidyverse)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(dplyr)

p <- ggboxplot(marker_in, x = "Markers", y = "expression",
	                      color = 'CMS_group', palette = "nejm",add = "jitter")+
  stat_compare_means(aes(group = CMS_group),label = "p.signif")+rotate_x_text(angle = 45)
tiff(file='Figures/Figure4K.CMS_group_expression.tiff',width=25,height=15,units="cm",compression="lzw",bg="white",res=600)
print(p)
dev.off()

#cibersort  analysis
#load('data/data_left.Rdata')
load('data/TCGA_colon_expr_filt.Rdata')
expr_left<-TCGA_colon_expr_filt[,data_left$sample]
write.table(expr_left,file = 'data/expr_left.file',sep='\t',quote = F)

source('data/CIBERSORT.R')
# Define LM22 file
LM22.file <- "data/LM22.txt"
#for left
left.file<-'expr_left.file'
left.TME.results = CIBERSORT(LM22.file, left.file, perm = 1000, QN = TRUE)  #耗内存3-4G
#save(left.TME.results,file = 'data/cibersort.left.TME.results.Rdata')
load('data/cibersort.left.TME.results.Rdata')
#write.table(left.TME.results, "cibersort.CIBERSORT-Results-left.txt",sep = "\t", row.names = T, col.names = T, quote = F)
library(dplyr) 
library(tidyr) 
library(tibble) 
dd1 <- left.TME.results %>%    as.data.frame() %>%    rownames_to_column("GeneName") %>%    pivot_longer(cols=2:23,names_to= "celltype", values_to = "Proportion")
dd1_filt <- dd1[dd1$`P-value`<=0.05,] #取p-value小于0.05

library(ggplot2)
ggplot(dd1_filt,aes(GeneName,Proportion,fill = celltype)) + geom_bar(position = "stack",stat = "identity")+   theme_bw()+ theme(legend.position = "right")
ggsave("Figures/Figure_S4A.left.CIBERSORT1.tiff",width = 7,height = 4)
dev.off()


left.TME.results<-data.frame(rownames(left.TME.results),left.TME.results)
colnames(left.TME.results)[1]='sample'

left.TME.merge<-merge(left.TME.results,data_left,by='sample')
colnames(left.TME.merge)<-gsub(pattern = "\\.",replacement=' ',colnames(left.TME.merge))
colnames(left.TME.merge)<-gsub(pattern = " Tregs ",replacement='(Tregs)',colnames(left.TME.merge))

immune_cells<-colnames(left.TME.merge)[2:22]
for (cell in immune_cells){
	  x<-left.TME.merge[,cell]
  y<-left.TME.merge[,'total_risk_score']
    corT <- cor.test(x,y)
    z <- lm(y~x)
      estimate <- corT$estimate
      cor <- round(estimate,3)
        pvalue <- corT$p.value   
        pval <- signif(pvalue,4) 
	  pval <- format(pval,scientific=T)
	  picname=paste(cell,"cor.tiff",sep='-')
	    outTiff <- paste('Figures/',picname,sep='/')
	    tiff(file=outTiff,width=15,height=15,units="cm",compression="lzw",bg="white",res=300)
	      plot(x,y,type="p",pch=16,main=paste("Cor=",cor,"(p-value=",pval,")",sep=""),xlab=cell,ylab="Risk Score")
	      lines(x,fitted(z),col=2)
	        dev.off()
	        line=paste(cell,round(corT$estimate,4),corT$p.value,sep='\t')
		  
		  write.table(line, "Figures/cor_result.xls",col.names = F, append = T,row.names = F,quote = F)
}

#for right CRC
load('data/data_right.Rdata')
expr_right<-TCGA_colon_expr_filt[,data_right$sample]
write.table(expr_right,file = 'expr_right.file',sep='\t',quote = F)

source('data/CIBERSORT.R')
# Define LM22 file
LM22.file <- "data/LM22.txt"
#for right
right.file<-'data/expr_right.file'
right.TME.results = CIBERSORT(LM22.file, right.file, perm = 1000, QN = TRUE)  #耗内存3-4G
#right.TME.results<-read.table('data/cibersort.CIBERSORT-Results-right.txt',sep = '\t',header = T)
#save(right.TME.results,file = 'data/cibersort.right.TME.results.Rdata')

#write.table(right.TME.results, "cibersort.CIBERSORT-Results-right.txt",sep = "\t", row.names = T, col.names = T, quote = F)

library(dplyr) 
library(tidyr) 
library(tibble) 

load('cibersort.right.TME.results.Rdata')
dd1 <- right.TME.results %>%    as.data.frame() %>%    rownames_to_column("GeneName") %>%    pivot_longer(cols=2:23,names_to= "celltype", values_to = "Proportion")
dd1_filt <- dd1[dd1$`P-value`<=0.05,] #p-value<=0.05

library(ggplot2)
ggplot(dd1_filt,aes(GeneName,Proportion,fill = celltype)) + geom_bar(position = "stack",stat = "identity")+   theme_bw()+ theme(legend.position = "right")
ggsave("data/Figure_S4B.right.CIBERSORT1.tiff",width = 7,height = 4)
dev.off()

right.TME.results<-data.frame(rownames(right.TME.results),right.TME.results)
colnames(right.TME.results)[1]='sample'

right.TME.merge<-merge(right.TME.results,data_right,by='sample')

colnames(right.TME.merge)<-gsub(pattern = "\\.",replacement=' ',colnames(right.TME.merge))
colnames(right.TME.merge)<-gsub(pattern = " Tregs ",replacement='(Tregs)',colnames(right.TME.merge))

immune_cells<-colnames(right.TME.merge)[2:22]
for (cell in immune_cells){
	  x<-right.TME.merge[,cell]
  y<-right.TME.merge[,'total_risk_score']
    corT <- cor.test(x,y)
    z <- lm(y~x)
      estimate <- corT$estimate
      cor <- round(estimate,3)
        pvalue <- corT$p.value   
        pval <- signif(pvalue,4) 
	  pval <- format(pval,scientific=T)
	  picname=paste(cell,"cor.tiff",sep='-')
	    outTiff <- paste('Figures/',picname,sep='/')
	    tiff(file=outTiff,width=15,height=15,units="cm",compression="lzw",bg="white",res=300)
	      plot(x,y,type="p",pch=16,main=paste("Cor=",cor,"(p-value=",pval,")",sep=""),xlab=cell,ylab="Risk Score")
	      lines(x,fitted(z),col=2)
	        dev.off()
	        line=paste(cell,round(corT$estimate,4),corT$p.value,sep='\t')
		  
		  write.table(line, "Figures/cor_result.xls",col.names = F, append = T,row.names = F,quote = F)
}

load('data/cibersort.right.TME.results.Rdata')
colnames(right.TME.results)<-gsub(pattern = "\\.",replacement=' ',colnames(right.TME.results))
colnames(right.TME.results)<-gsub(pattern = " Tregs ",replacement='(Tregs)',colnames(right.TME.results))

load('data/cibersort.left.TME.results.Rdata')
colnames(left.TME.results)<-gsub(pattern = "\\.",replacement=' ',colnames(left.TME.results))
colnames(left.TME.results)<-gsub(pattern = " Tregs ",replacement='(Tregs)',colnames(left.TME.results))


cellnames=colnames(right.TME.results)[1:22]
cell_in=array()
for (cell in cellnames){
	  uni_cell_left<-data.frame(cell,left.TME.results[,cell],rownames(left.TME.results),'LCC')
  uni_cell_right<-data.frame(cell,right.TME.results[,cell],rownames(right.TME.results),'RCC')
    colnames(uni_cell_left)<-c('celltypes','proportion','sample','Group')
    colnames(uni_cell_right)<-c('celltypes','proportion','sample','Group')
      
      uni_cell_in<-rbind(uni_cell_left,uni_cell_right)
      
      cell_in<-rbind(cell_in,uni_cell_in)
}

library(dplyr)
cell_in<-cell_in %>% na.omit()


library(cowplot)
library(tidyverse)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(dplyr)
p <- ggboxplot(cell_in, x = "celltypes", y = "proportion",
	                      color = 'Group', palette = "nejm")+
  stat_compare_means(aes(group = Group),method = "wilcox.test",label = "p.signif")+rotate_x_text(angle = 45)
tiff(file='Figures/Figure5A.CIBERSORT.right_left.tiff',width=25,height=15,units="cm",compression="lzw",bg="white",res=600)
print(p)
dev.off()

#load('data/data_left.Rdata')
#load('data/data_right.Rdata')

library(maftools)
left_sample<-data.frame(data_left$sample,"left-side")
right_sample<-data.frame(data_right$sample,'right-side')
colnames(left_sample)<-colnames(right_sample)<-c("Tumor_Sample_Barcode",'Location')

TCGA_somatic_location<-rbind(left_sample,right_sample)
#write.table(TCGA_somatic_location,file = 'sample_location.txt',quote = F,sep='\t',row.names = F,col.names = T)

var_maf2 = read.maf(maf ="data/TCGA.COAD.mutect.check.maf",clinicalData = TCGA_somatic_location)

TCGA_somatic_location_right<-TCGA_somatic_location[TCGA_somatic_location$Location=='right',]
TCGA_somatic_location_left<-TCGA_somatic_location[TCGA_somatic_location$Location=='left',]

save(TCGA_somatic_location,TCGA_somatic_location_right,TCGA_somatic_location_left,file='data/TCGA_somatic_location.Rdata')

#load('data/TCGA_somatic_location.Rdata')
#RCC
var_maf_right = read.maf(maf ="data/TCGA.COAD.mutect.check.right.maf",clinicalData = TCGA_somatic_location_right) 
tiff(file='Figures/Figure4G.plotmafSummary_right.tiff',width=15,height=15,units="cm",compression="lzw",bg="white",res=600)
p=plotmafSummary(maf = var_maf_right,rmOutlier = TRUE,addStat = 'median',dashboard = TRUE,titvRaw = FALSE)
print(p)
dev.off()
#LCC
var_maf_left = read.maf(maf ="data/TCGA.COAD.mutect.check.left.maf",clinicalData = TCGA_somatic_location_left)

tiff(file='Figures/Figure4F.plotmafSummary_left.tiff',width=15,height=15,units="cm",compression="lzw",bg="white",res=600)
p=plotmafSummary(maf = var_maf_left,rmOutlier = TRUE,addStat = 'median',dashboard = TRUE,titvRaw = FALSE)
print(p)
dev.off()

#get the P-value of Figure 5A
maf_compare <- mafCompare(m1 = var_maf_right, m2 = var_maf_left,m1Name = 'right-side', m2Name = 'left-side', minMut = 2, useCNV =FALSE)
write.table(maf_compare$results, file="Tables/maf_compare.xls", quote=FALSE, row.names=T,sep="\t")

tiff(file='Figures/Figure5A.coBarplot_right_left.tiff',width=15,height=10,units="cm",compression="lzw",bg="white",res=600)
p=coBarplot(m1 = var_maf_right, m2 = var_maf_left, m1Name = 'right-side', m2Name = 'left-side')
print(p)
dev.off()

left_sample_for_somatic<-data_left[,c('sample','OS','OS.time')]
right_sample_for_somatic<-data_right[,c('sample','OS','OS.time')]
colnames(left_sample_for_somatic)<-colnames(right_sample_for_somatic)<-c('Tumor_Sample_Barcode','OS','OS.time')

#LCC
var_maf_left = read.maf(maf ="data/TCGA.COAD.mutect.check.left.maf",clinicalData = left_sample_for_somatic)

#significant mutated genes: APC/TP53/KRAS/SYNE1/MUC16
for(i in c('APC','TP53','KRAS','SYNE1','MUC16')){
	  tiff(file=paste0("Figures/",i,'_OS_left.tiff'),width=15,height=10,units="cm",compression="lzw",bg="white",res=600)
  p=mafSurvival(maf = var_maf_left, genes = i, time = 'OS.time', Status = 'OS', isTCGA = FALSE)
    print(p)
    dev.off()
}


#RCC
var_maf_right = read.maf(maf ="data/TCGA.COAD.mutect.check.right.maf",clinicalData = right_sample_for_somatic)

# significant mutated genes:APC/TP53/KRAS/SYNE1/MUC16
for(i in c('APC','TP53','KRAS','SYNE1','MUC16')){
	  tiff(file=paste0("Figures/",i,'_OS_right.tiff'),width=15,height=10,units="cm",compression="lzw",bg="white",res=600)
  p=mafSurvival(maf = var_maf_right, genes = i, time = 'OS.time', Status = 'OS', isTCGA = FALSE)
    print(p)
    dev.off()
}

#TMB
tiff(file='Figures/Figure_S4K.tmb_left.tiff',width=10,height=10,units="cm",compression="lzw",bg="white",res=600)
tmb_left= tmb(maf = var_maf_left)
print(tmb_left)
dev.off()

tiff(file='Figures/Figure_S4L.tmb_right.tiff',width=10,height=10,units="cm",compression="lzw",bg="white",res=600)
tmb_right= tmb(maf = var_maf_right)
print(tmb_right)
dev.off()

tmb_left<-data.frame(tmb_left,'left-side')
tmb_right<-data.frame(tmb_right,'right-side')
colnames(tmb_left)[5]<-colnames(tmb_right)[5]<-'Location'

tmb_location<-rbind(tmb_left,tmb_right)
save(tmb_location,file='data/TMB.Rdata')

load('data/TMB.Rdata')

library(cowplot)
library(tidyverse)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(dplyr)

tmb_location$Location<-ifelse(tmb_location$Location=='left-side','LCC','RCC')
tmb_location$Location<-as.factor(tmb_location$Location)

mydata<-data.frame(tmb_location$Tumor_Sample_Barcode,tmb_location$total_perMB_log,tmb_location$Location)
colnames(mydata)<-c('sample','Tumor Mutaion Burden(total perMB)','Group')
rownames(mydata)<-mydata$sample
mydata$sample<-NULL
mydata$Group<-as.factor(mydata$Group)
mydata<-mydata %>% gather(key = "Group",value = "Tumor Mutaion Burden(total perMB)") %>% na.omit() 

my_comparisons<-list(c('RCC','LCC'))
tiff(file='Figures/Figure5F.TMB_location.tiff',width=15,height=15,units="cm",compression="lzw",bg="white",res=600)
p<-ggboxplot(mydata,x='Group',y='Tumor Mutaion Burden(total perMB)',color='Group',palette = 'jama',add = 'jitter')
p+stat_compare_means(comparisons = my_comparisons)
dev.off()

tmb_left$sample<-tmb_left$Tumor_Sample_Barcode
tmb_OS_left<-merge(data_left[,c('sample','OS','OS.time')],tmb_left,by='sample')
tmb_OS_left$`Tumor Mutaion Burden(total perMB)`<-tmb_OS_left$total_perMB

library(survival)
library(survminer)
res.cut<-surv_cutpoint(tmb_OS_left,time='OS.time',event = 'OS',variables = c('Tumor Mutaion Burden(total perMB)'))
res.cat <- surv_categorize(res.cut)


library(cutoff)
mm<-tmb_OS_left
logresult<- logrank(data=mm,
		                        time = 'OS.time',
					                    y='OS', 
					                    x='Tumor Mutaion Burden(total perMB)',
							                        cut.numb=1,
							                        n.per=0.10,
										                    y.per=0.10,
										                    p.cut=0.05,
												                        round=5)
logresult
cutoff<-2.26
tmb_OS_left$TMB<-ifelse(tmb_OS_left$`Tumor Mutaion Burden(total perMB)`>cutoff,'High','Low')

library(survival)
library(survminer)
fit <- survfit(Surv(as.numeric(OS.time), as.numeric(OS)) ~tmb_OS_left$TMB, data = tmb_OS_left)
p=ggsurvplot(fit, 
	                  #pval = TRUE, 
	                  pval = toupper(surv_pvalue(fit)[4]),#p改为大写
			               conf.int = TRUE,
			               #risk.table = TRUE, # Add risk table
			               risk.table.col = "strata", # Change risk table color by groups
				                    linetype = "strata", # Change line type by groups
				                    surv.median.line = "hv", # Specify median survival
						                 ggtheme = theme_bw(), # Change ggplot2 theme
						                 legend.labs = c("High", "Low"), 
								              legend.title = "TMB", 
								              palette = c( "#f03b20","#2c7fb8"))

tiff(file='Figures/Figure5G.TMB_OS_left.tiff',width=15,height=10,units="cm",compression="lzw",bg="white",res=600)
print(p)
dev.off()

#RCC
tmb_right$sample<-tmb_right$Tumor_Sample_Barcode
tmb_OS_right<-merge(data_right[,c('sample','OS','OS.time')],tmb_right,by='sample')
tmb_OS_right$`Tumor Mutaion Burden(total perMB)`<-tmb_OS_right$total_perMB

res.cut<-surv_cutpoint(tmb_OS_right,time='OS.time',event = 'OS',variables = c('Tumor Mutaion Burden(total perMB)'))
res.cat <- surv_categorize(res.cut)

library(cutoff)
mm<-tmb_OS_right
logresult<- logrank(data=mm,
		                        time = 'OS.time',
					                    y='OS', 
					                    x='Tumor Mutaion Burden(total perMB)',
							                        cut.numb=1,
							                        n.per=0.10,
										                    y.per=0.10,
										                    p.cut=0.2,
												                        round=5)
logresult
cutoff<-1.9
tmb_OS_right$TMB<-ifelse(tmb_OS_right$`Tumor Mutaion Burden(total perMB)`>cutoff,'High','Low')

library(survival)
library(survminer)
fit <- survfit(Surv(as.numeric(OS.time), as.numeric(OS)) ~tmb_OS_right$TMB, data = tmb_OS_right)
p=ggsurvplot(fit, 
	                  #pval = TRUE, 
	                  pval = toupper(surv_pvalue(fit)[4]),#p改为大写
			               conf.int = TRUE,
			               #risk.table = TRUE, # Add risk table
			               risk.table.col = "strata", # Change risk table color by groups
				                    linetype = "strata", # Change line type by groups
				                    surv.median.line = "hv", # Specify median survival
						                 ggtheme = theme_bw(), # Change ggplot2 theme
						                 legend.labs = c("High", "Low"), 
								              legend.title = "TMB", 
								              palette = c( "#f03b20","#2c7fb8"))

tiff(file='Figures/Figure5H.TMB_OS_right1.tiff',width=15,height=10,units="cm",compression="lzw",bg="white",res=600)
print(p)
dev.off()

#the expression level of immune checkpoint targets
immune_location<-rbind(left_sample,right_sample)
colnames(immune_location)[1]<-'sample'
imm_RNA<-TCGA_colon_expr_filt[c("CD274","PDCD1","CTLA4",'HAVCR2','LAG3','TIGIT'),] #
imm_RNA<-t(imm_RNA)
imm_RNA<-data.frame(rownames(imm_RNA),imm_RNA)
colnames(imm_RNA)[1]<-'sample'
imm_RNA_location<-merge(imm_RNA,immune_location,by='sample')
rownames(imm_RNA_location)<-imm_RNA_location$sample
imm_RNA_location$sample<-NULL

imm_RNA_location$Group<-imm_RNA_location$Location
imm_RNA_location$Location<-NULL

genenames=colnames(imm_RNA_location)[1:(ncol(imm_RNA_location)-1)]
gene_in=array()
for (gene in genenames){
	  uni_gene_in<-data.frame(gene,imm_RNA_location[,gene],rownames(imm_RNA_location),imm_RNA_location[,'Group'])
 
  colnames(uni_gene_in)<-c('GeneName','Expression','sample','Group')
    gene_in<-rbind(gene_in,uni_gene_in)
}

library(dplyr)
gene_in<-gene_in %>% na.omit()
gene_in$Group<-ifelse(gene_in$Group=='left-side','LCC','RCC')


left_gene_risk<-merge(gene_in,risk_score_table_multi_forKM_left,by='sample')
left_gene_risk$Group<-NULL
left_gene_risk$Group<-left_gene_risk$high_low
left_gene_risk$high_low<-NULL


p <- ggboxplot(left_gene_risk, x = "GeneName", y = "Expression",
	                      color = 'Group', palette = "nejm", add = "jitter")+
  stat_compare_means(aes(group = Group),method = "wilcox.test",label = "p.format")
#+rotate_x_text(angle = 45)
tiff(file='Figures/Figure6C.expression-left.tiff',width=18,height=15,units="cm",compression="lzw",bg="white",res=600)
print(p)
dev.off()

right_gene_risk<-merge(gene_in,risk_score_table_multi_forKM_right,by='sample')
right_gene_risk$Group<-NULL
right_gene_risk$Group<-right_gene_risk$high_low
right_gene_risk$high_low<-NULL


p <- ggboxplot(right_gene_risk, x = "GeneName", y = "Expression",
	                      color = 'Group', palette = "nejm", add = "jitter")+
  stat_compare_means(aes(group = Group),method = "wilcox.test",label = "p.format")
#+rotate_x_text(angle = 45)
tiff(file='Figures/Figure6D.expression-right.tiff',width=18,height=15,units="cm",compression="lzw",bg="white",res=600)
print(p)
dev.off()

library(cowplot)
library(tidyverse)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(dplyr)
p <- ggboxplot(gene_in, x = "GeneName", y = "Expression",
	                      color = 'Group', palette = "nejm", add = "jitter")+
  stat_compare_means(aes(group = Group),method = "wilcox.test",label = "p.format")
#+rotate_x_text(angle = 45)
  tiff(file='Figures/Figure6B.expression-location.tiff',width=18,height=15,units="cm",compression="lzw",bg="white",res=600)
  print(p)
  dev.off()

  #pathways
  tiff(file='Figures/Figure5D.OncogenicPathways_left.tiff',width=15,height=10,units="cm",compression="lzw",bg="white",res=600)
  OncogenicPathways(maf = var_maf_left)
  dev.off()

  tiff(file='maftools/Figure5E.OncogenicPathways_right.tiff',width=15,height=10,units="cm",compression="lzw",bg="white",res=600)
  OncogenicPathways(maf = var_maf_right)
  dev.off()

  #IC50
  library(oncoPredict)
  library(data.table)
  library(gtools)
  library(reshape2)
  library(ggpubr)
  th=theme(axis.text.x = element_text(angle = 45,vjust = 0.5))

  dir='oncopredict/rawData/Training Data/'
  GDSC2_Expr = readRDS(file=file.path(dir,'GDSC2_Expr (RMA Normalized and Log Transformed).rds'))
  GDSC2_Res = readRDS(file = file.path(dir,"GDSC2_Res.rds"))
  GDSC2_Res <- exp(GDSC2_Res) 

  #load('data/data_left.Rdata')
  #load('data/TCGA_colon_expr_filt.Rdata')
  expr_left<-TCGA_colon_expr_filt[,data_left$sample]

  calcPhenotype(trainingExprData = GDSC2_Expr,
		              trainingPtype = GDSC2_Res,
			                    testExprData = as.matrix(expr_left),
			                    batchCorrect = 'eb',  #   "eb" for ComBat  
					                  powerTransformPhenotype = TRUE,
					                  removeLowVaryingGenes = 0.2,
							                minNumSamples = 10, 
							                printOutput = TRUE, 
									              removeLowVaringGenesFrom = 'rawData' )

  library(data.table)
  left_Ptype <- fread('data/calcPhenotype_Output/DrugPredictions_left.csv', data.table = F)

  load('data/data_right.Rdata')
  expr_right<-TCGA_colon_expr_filt[,data_right$sample]

  calcPhenotype(trainingExprData = GDSC2_Expr,
		              trainingPtype = GDSC2_Res,
			                    testExprData = as.matrix(expr_right),
			                    batchCorrect = 'eb',  #   "eb" for ComBat  
					                  powerTransformPhenotype = TRUE,
					                  removeLowVaryingGenes = 0.2,
							                minNumSamples = 10, 
							                printOutput = TRUE, 
									              removeLowVaringGenesFrom = 'rawData' )

  library(data.table)
  right_Ptype <- fread('data/oncopredict/calcPhenotype_Output/DrugPredictions_right.csv', data.table = F)



  #the common drugs for colon cancer
  rownames(left_Ptype)<-left_Ptype$V1
  left_Ptype$V1<-NULL

  drug_left_names<-c()

  drug_left_sen<-array()
  for (i in 1:ncol(left_Ptype)) {
	    if (mean(left_Ptype[,i])<=5) {
		        #drug_left_names<-append(drug_left_names,colnames(left_Ptype)[i])
		        drug_left_sen<-rbind(drug_left_sen,c(colnames(left_Ptype)[i],mean(left_Ptype[,i])))
    }
    
  }

  drug_left_sen<-drug_left_sen %>% na.omit()
  colnames(drug_left_sen)<-c('DrugName','IC50 of LCC')
  drug_right_sen<-array()
  for (i in 1:ncol(right_Ptype)) {
	    if (mean(right_Ptype[,i])<=5) {
		        #drug_right_names<-append(drug_right_names,colnames(right_Ptype)[i])
		        drug_right_sen<-rbind(drug_right_sen,c(colnames(right_Ptype)[i],mean(right_Ptype[,i])))
    }
    
  }

  drug_right_sen<-drug_right_sen %>% na.omit()
  colnames(drug_right_sen)<-c('DrugName','IC50 of RCC')

  drug_left_right_IC50<-merge(drug_left_sen,drug_right_sen,by='DrugName')
  colnames(drug_left_right_IC50)<-c('DrugName','average IC50 of LCC','average IC50 of RCC')

  rownames(right_Ptype)<-right_Ptype$V1
  right_Ptype$V1<-NULL

  drug_right_names<-c()
  for (i in 1:ncol(right_Ptype)) {
	    if (mean(right_Ptype[,i])<=5) {
		        drug_right_names<-append(drug_right_names,colnames(right_Ptype)[i])
    }
    
  }

  common_drug<-intersect(drug_right_names,drug_left_names)

  drug_left<-data.frame(left_Ptype[,c('Camptothecin_1003','Docetaxel_1819','Rapamycin_1084','Paclitaxel_1080',
				      'Teniposide_1809','Mitoxantrone_1810','Dactinomycin_1911','Vinorelbine_2048')],'left-side')

  drug_right<-data.frame(right_Ptype[,c('Camptothecin_1003','Docetaxel_1819','Rapamycin_1084','Paclitaxel_1080',
					                                    'Teniposide_1809','Mitoxantrone_1810','Dactinomycin_1911','Vinorelbine_2048')],'right-side')

  colnames(drug_left)<-colnames(drug_right)<-c('Camptothecin','Docetaxel','Rapamycin','Paclitaxel',
					                                                    'Teniposide','Mitoxantrone','Dactinomycin','Vinorelbine','Location')

  drug_all<-rbind(drug_left,drug_right)

  drugnames=colnames(drug_all)[1:(ncol(drug_all)-1)]
  drug_in=array()
  for (i in 1:length(drugnames)){
	    
	    uni_drug<-data.frame(drugnames[i],drug_all[,drugnames[i]],rownames(drug_all),drug_all[,'Location'])
   
    colnames(uni_drug)<-c('DrugName','IC50','sample','Location')
      
      drug_in<-rbind(drug_in,uni_drug)
  }


  drug_in<-drug_in %>% na.omit()


  library(cowplot)
  library(tidyverse)
  library(ggplot2)
  library(ggsci)
  library(ggpubr)
  library(dplyr)
  p <- ggboxplot(drug_in, x = "DrugName", y = "IC50",
		                color = 'Location', palette = "nejm")+
  stat_compare_means(aes(group = Location),method = "wilcox.test",label = "p.format")+rotate_x_text(angle = 45)
tiff(file='Figures/Figure6E.IC50.right_left.tiff',width=25,height=15,units="cm",compression="lzw",bg="white",res=600)
print(p)
dev.off()

drug_left<-data.frame(rownames(drug_left),drug_left)
colnames(drug_left)[1]<-'sample'
drug_left_risk<-merge(drug_left,risk_score_table_multi_forKM_left[,c('sample','high_low')],by='sample')
drug_left_risk$Risk_Score<-drug_left_risk$high_low
drug_left_risk$high_low<-NULL
rownames(drug_left_risk)<-drug_left_risk$sample
drug_left_risk$sample<-NULL

drug_in_left=array()
for (i in 1:length(drugnames)){
	  
	  uni_drug<-data.frame(drugnames[i],drug_left_risk[,drugnames[i]],rownames(drug_left_risk),drug_left_risk[,'Risk_Score'])
  
  colnames(uni_drug)<-c('DrugName','IC50','sample','Risk_Score')
    
    drug_in_left<-rbind(drug_in_left,uni_drug)
}


drug_in_left<-drug_in_left %>% na.omit()


library(cowplot)
library(tidyverse)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(dplyr)
p <- ggboxplot(drug_in_left, x = "DrugName", y = "IC50",
	                      color = 'Risk_Score', palette = "nejm")+
  stat_compare_means(aes(group = Risk_Score),method = "wilcox.test",label = "p.format")+rotate_x_text(angle = 45)
tiff(file='Figures/Figure6F.IC50.left.tiff',width=25,height=15,units="cm",compression="lzw",bg="white",res=600)
print(p)
dev.off()

#RCC
drug_right<-data.frame(rownames(drug_right),drug_right)
colnames(drug_right)[1]<-'sample'
drug_right_risk<-merge(drug_right,risk_score_table_multi_forKM_right[,c('sample','high_low')],by='sample')
drug_right_risk$Risk_Score<-drug_right_risk$high_low
drug_right_risk$high_low<-NULL
rownames(drug_right_risk)<-drug_right_risk$sample
drug_right_risk$sample<-NULL

drug_in_right=array()
for (i in 1:length(drugnames)){
	  
	  uni_drug<-data.frame(drugnames[i],drug_right_risk[,drugnames[i]],rownames(drug_right_risk),drug_right_risk[,'Risk_Score'])
  
  colnames(uni_drug)<-c('DrugName','IC50','sample','Risk_Score')
    
    drug_in_right<-rbind(drug_in_right,uni_drug)
}


drug_in_right<-drug_in_right %>% na.omit()


library(cowplot)
library(tidyverse)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(dplyr)
 ggboxplot(drug_in_right, x = "DrugName", y = "IC50",
	                  color = 'Risk_Score', palette = "nejm")+
  stat_compare_means(aes(group = Risk_Score),method = "wilcox.test",label = "p.format")+rotate_x_text(angle = 45)

library(cowplot)
library(tidyverse)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(dplyr)

nm<-colnames(drug_right_risk)[1:(ncol(drug_right_risk)-2)]
for(i in c(1:(ncol(drug_right_risk)-2))){
	  p<-ggplot(drug_right_risk[,c(colnames(drug_right_risk)[i],"Risk_Score")], aes(x =Risk_Score,y = drug_right_risk[,i],fill=Risk_Score))+
		      geom_violin(trim = FALSE)+geom_boxplot(width = 0.2)+stat_summary(fun=mean,size = 5, geom = "point",colour = "red")+
		          theme_set(theme_bw())+ylab(paste0('IC50 of ',nm[i]))+xlab("Risk score")+theme(legend.position = 'none')+
			      stat_compare_means(method = "wilcox.test",label.x = 1,label.y = max(drug_right_risk[,i])+min(drug_right_risk[,i]),method.args = list(var.equal = TRUE)) #wilcox.test
		        
		        tiff(file=paste0("Figures/",nm[i],'.tiff'),width=10,height=10,units="cm",compression="lzw",bg="white",res=600)
			  print(p)
			  dev.off()
}


#validation:GSE103479
library(GEOquery)
gset <- getGEO('GSE103479',destdir = ".",
	                      AnnotGPL = F,
			                     getGPL = F)

load('data/GSE103479_gset.Rdata')
GSE103479_exprSet <- data.frame(exprs(gset[[1]])) 

GSE103479_exprSet$ID<-rownames(GSE103479_exprSet)
id2symbol<-read.table(file = 'data/GSE103479.id2symbol.txt',sep='\t',header = T)
GSE103479_exprSet_anno<-merge(id2symbol,GSE103479_exprSet,by='ID')
GSE103479_exprSet_anno$ID<-NULL
library(dplyr)
GSE103479_exprSet_anno<-distinct(GSE103479_exprSet_anno,Gene.Symbol,.keep_all = TRUE)
rownames(GSE103479_exprSet_anno)<-GSE103479_exprSet_anno$Gene.Symbol

GSE103479_pd=pData(gset[[1]])
GSE103479_sam<-intersect(colnames(GSE103479_exprSet_anno),GSE103479_pd$geo_accession)
GSE103479_pd<-GSE103479_pd[GSE103479_pd$geo_accession %in% GSE103479_sam,]
GSE103479_exprSet_anno<-GSE103479_exprSet_anno[,GSE103479_sam]
save(GSE103479_exprSet_anno,file='data/GSE103479_exprSet_anno.Rdata')

GSE103479_pd$tumor_site<-ifelse(
				  toupper(substr(GSE103479_pd$`tumour site:ch1`,1,4))=='LEFT' | 
					      toupper(substr(GSE103479_pd$`tumour site:ch1`,1,4))=='DESC' | 
					          toupper(substr(GSE103479_pd$`tumour site:ch1`,1,5))=='SIGMO' |
						      toupper(substr(GSE103479_pd$`tumour site:ch1`,1,4))=='RECT','left-sided',
					        ifelse(toupper(substr(GSE103479_pd$`tumour site:ch1`,1,4))=='RIGH' |
						                  toupper(substr(GSE103479_pd$`tumour site:ch1`,1,4))=='CAEC' |
								             GSE103479_pd$`tumour site:ch1`=='Hepatic flexure' |
									                GSE103479_pd$`tumour site:ch1`=='Transverse/left colon' |
											           toupper(substr(GSE103479_pd$`tumour site:ch1`,1,3))=='ASC' ,'right-sided','unknown'
											            
											     )
						  
						)

GSE103479_pd$OS.time<-as.numeric(GSE103479_pd$`overall survival time:ch1`)


GSE103479_pd_use<-GSE103479_pd[GSE103479_pd$tumor_site!='unknown',]   
GSE103479_pd_use<-GSE103479_pd_use[GSE103479_pd_use$OS.time >0,]
GSE103479_pd_use<-GSE103479_pd_use %>% na.omit()

GSE103479_pd_use$OS<-ifelse(GSE103479_pd_use$`status alive.dead:ch1`=='Alive',0,1)
GSE103479_pd_compare<-data.frame(GSE103479_pd_use$geo_accession,GSE103479_pd_use$OS.time,GSE103479_pd_use$tumor_site)

colnames(GSE103479_pd_compare)<-c('sample','OS.time','Group')
GSE103479_pd_compare <-GSE103479_pd_compare %>% na.omit()

library(tidyr)
rownames(GSE103479_pd_compare)<-GSE103479_pd_compare$sample 
GSE103479_pd_compare$sample<-NULL
GSE103479_pd_compare$Group<-as.factor(GSE103479_pd_compare$Group)

mydata<-GSE103479_pd_compare %>% gather(key = "Group",value = "OS.time") %>% na.omit() 
my_comparisons <- list( c("left-sided", "right-sided"))

library(cowplot)
library(tidyverse)
library(ggplot2)
library(ggsci)
library(ggpubr)

tiff(file='Figures/Figure_S1B.site_os_GSE103479.tiff',width=15,height=15,units="cm",compression="lzw",bg="white",res=600)
p<-ggboxplot(mydata,x='Group',y='OS.time',color='Group',palette = 'jama',add = 'jitter')
p+stat_compare_means(comparisons = my_comparisons)
dev.off()

#CMS
load('data/GSE103479_pd_use.Rdata')
GSE103479_pd_CMS<-data.frame(GSE103479_pd_use$geo_accession,GSE103479_pd_use$tumor_site,GSE103479_pd_use$`cms subgroup:ch1`)
colnames(GSE103479_pd_CMS)<-c('sample','tumor_site','CMS_group')
GSE103479_pd_CMS<-GSE103479_pd_CMS[GSE103479_pd_CMS$CMS_group!='UNK',]
GSE103479_pd_CMS$tumor_site<-ifelse(GSE103479_pd_CMS$tumor_site=='left-sided','LCC','RCC')
tiff(file='Figures/Figure3J.tumorsite_cmsgroup.tiff',width=10,height=15,units="cm",compression="lzw",bg="white",res=600)
p<-ggplot(data = GSE103479_pd_CMS, aes(x = tumor_site, fill = CMS_group)) + geom_bar( position = "fill")+
	  scale_fill_manual(values = pal_aaas('default')(5))+
	    ylab("proportion")+
	      geom_hline(aes(yintercept = 0.45), colour = "red", size=0.5,linetype="dashed") +
	        geom_hline(aes(yintercept = 0.55), colour = "red", size=0.5,linetype="dashed") +
		  theme(panel.background = element_blank(),panel.grid.major= element_line(color = "white"),panel.grid.minor =element_line(color= "white"),legend.position='right')
	  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
	  print(p)
	  dev.off()

	  #construct the prognostic risk score model for LCC
	  left_markers<-c("FOSB","RPL35","REG1A", "TESC","C11orf96")
	  load('data/GSE103479_pd_use.Rdata')
	  GSE103479_left<-GSE103479_pd_use[GSE103479_pd_use$tumor_site=='left-sided',]

	  GSE103479_gene_survival_data_left<-cbind(GSE103479_left[,c('geo_accession','OS','OS.time')],t(GSE103479_exprSet_anno[left_markers,GSE103479_left$geo_accession]))
	  colnames(GSE103479_gene_survival_data_left)[1]<-'sample'
	  GSE103479_gene_survival_data_left$OS<-ifelse(GSE103479_gene_survival_data_left$OS=='Alive',0,1)

	  GSE103479_gene_survival_data_left[,2:8]=as.numeric(unlist(GSE103479_gene_survival_data_left[,2:8]))

	  save(GSE103479_gene_survival_data_left,file='GSE103479_gene_survival_data_left.Rdata')

	  multi_variate_cox_2<-processMultiCOX(left_markers, GSE103479_gene_survival_data_left)

	  C_index<-multi_variate_cox_2$concordance['concordance']

	  risk_score_table_multi_cox2<-riskscore(GSE103479_gene_survival_data_left,sig_gene_multi_cox_left,multi_variate_cox_2)
	  risk_score_table_multi_cox2<-risk_score_table_multi_cox2[risk_score_table_multi_cox2$OS.time>0,]
	  risk_score_table_multi_cox2_left<-risk_score_table_multi_cox2
	  #save(risk_score_table_multi_cox2_left,file='risk_score_table_multi_cox2_left.Rdata')
	  #evaluate AUCs between 3-5 years.
	  risk_score_table_multi_cox2_left$OS.time<-risk_score_table_multi_cox2_left$OS.time*30
	  for_multi_ROC<-multi_ROC(time_vector = c(365*seq(3,5,1)),risk_score_table = risk_score_table_multi_cox2_left)

	  #visualization for the ROC curves of multiple time points
	  #plot ROC
	  year3ROC<-for_multi_ROC[for_multi_ROC$Time_point==3*365,]
	  year4ROC<-for_multi_ROC[for_multi_ROC$Time_point==4*365,]
	  year5ROC<-for_multi_ROC[for_multi_ROC$Time_point==5*365,]
	  TP_3year = year3ROC$True_positive
	  FP_3year = year3ROC$False_positive
	  TP_4year = year4ROC$True_positive
	  FP_4year = year4ROC$False_positive
	  TP_5year = year5ROC$True_positive
	  FP_5year = year5ROC$False_positive
	  time_ROC_df=data.frame(TP_3year,FP_3year,TP_4year,FP_4year,TP_5year,FP_5year)
	  library(ggplot2)
	  tiff(file='Figures/Figure_S3B.OS_AUC_left.tiff',width=15,height=15,units="cm",compression="lzw",bg="white",res=600)
	  p=ggplot(data = time_ROC_df) +
		    geom_line(aes(x = FP_3year, y = TP_3year), size = 1, color = "#BC3C29FF") +
		      geom_line(aes(x = FP_4year, y = TP_4year), size = 1, color = "#0072B5FF") +
		        geom_line(aes(x = FP_5year, y = TP_5year), size = 1, color = "#E18727FF") +
			  geom_abline(slope = 1, intercept = 0, color = "grey", size = 1, linetype = 2) +
			    theme_bw() +
			      annotate("text",
				                  x = 0.75, y = 0.25, size = 4.5,
						             label = paste0("AUC at 3 years = ", sprintf("%.3f", year3ROC$AUC)), color = "#BC3C29FF"
						    ) +
  annotate("text",
	              x = 0.75, y = 0.15, size = 4.5,
		                 label = paste0("AUC at 4 years = ", sprintf("%.3f", year4ROC$AUC)), color = "#0072B5FF"
		        ) +
  annotate("text",
	              x = 0.75, y = 0.05, size = 4.5,
		                 label = paste0("AUC at 5 years = ", sprintf("%.3f", year5ROC$AUC)), color = "#E18727FF"
		        ) +
  labs(x = "False positive rate", y = "True positive rate") +
    theme(
	      axis.text = element_text(face = "bold", size = 11, color = "black"),
	          #axis.title.x = element_text(face = "bold", size = 14, color = "black", margin = margin(c(15, 0, 0, 0))),
	          #axis.title.y = element_text(face = "bold", size = 14, color = "black", margin = margin(c(0, 15, 0, 0)))
	        )

  print(p)
  dev.off()


  #K-M analysis
  AUC_max<-max(for_multi_ROC$AUC)
  AUC_max_time<-for_multi_ROC$Time_point[which(for_multi_ROC$AUC==AUC_max)]
  AUC_max_time<-AUC_max_time[!duplicated(AUC_max_time)]
  AUC_max_time<-AUC_max_time[length(AUC_max_time)]
  for_multi_ROC$Time_point<-as.factor(for_multi_ROC$Time_point)
  optimal_time_ROC_df<-for_multi_ROC[which(for_multi_ROC$Time_point==AUC_max_time),]
  cut.off<-optimal_time_ROC_df$Cut_values[which.max(optimal_time_ROC_df$True_positive - optimal_time_ROC_df$False_positive)]
  #cut.off<-0.005319551
  #cut.off<-median(risk_score_table_multi_cox2$total_risk_score)
  high_low<-risk_score_table_multi_cox2$total_risk_score>cut.off
  high_low[high_low==TRUE]<-'high'
  high_low[high_low==FALSE]<-'low'
  risk_score_table_multi_forKM<-cbind(risk_score_table_multi_cox2,high_low)

  risk_score_table_multi_forKM_left<-risk_score_table_multi_forKM
  #KM_plot 
  library(survminer)
  library(survival)
  risk_score_table_multi_forKM$OS[which(risk_score_table_multi_forKM$OS.time>AUC_max_time)]<-0
  risk_score_table_multi_forKM$OS.time[which(risk_score_table_multi_forKM$OS.time>AUC_max_time)]<-AUC_max_time
  fit_km<-survfit(Surv(OS.time,OS)~high_low,data = risk_score_table_multi_forKM)
  tiff(file='Figures/Figure_S3B.OS_riskscore_left.tiff',width=15,height=15,units="cm",compression="lzw",bg="white",res=300)
  p=ggsurvplot(fit_km,conf.int = F,pval = T,legend.title="Total risk score",
	                    legend.labs=c('high','low'),risk.table=T,
			                 #legend.labs=c(paste0('>',as.character(round(cut.off,2))),
			                 # paste0('<=',as.character(round(cut.off,2)))),risk.table=T,
			                 palette=c('red','blue'),surv.median.line='hv')
  print(p)
  dev.off()

  #RCC
  right_markers<-c("HSPA1A","CD69","GDF15","LGALS2")
  GSE103479_right<-GSE103479_pd_use[GSE103479_pd_use$tumor_site=='right-sided',]

  GSE103479_gene_survival_data_right<-cbind(GSE103479_right[,c('geo_accession','OS','OS.time')],t(GSE103479_exprSet_anno[right_markers,GSE103479_right$geo_accession]))
  colnames(GSE103479_gene_survival_data_right)[1]<-'sample'
  GSE103479_gene_survival_data_right$OS<-ifelse(GSE103479_gene_survival_data_right$OS=='Alive',0,1)

  GSE103479_gene_survival_data_right[,2:7]=as.numeric(unlist(GSE103479_gene_survival_data_right[,2:7]))

  save(GSE103479_gene_survival_data_right,file='data/GSE103479_gene_survival_data_right.Rdata')


  multi_variate_cox_2<-processMultiCOX(right_markers, GSE103479_gene_survival_data_right)

  C_index<-multi_variate_cox_2$concordance['concordance']

  risk_score_table_multi_cox2<-riskscore(GSE103479_gene_survival_data_right,right_markers,multi_variate_cox_2)
  risk_score_table_multi_cox2<-risk_score_table_multi_cox2[risk_score_table_multi_cox2$OS.time>0,]
  risk_score_table_multi_cox2_right<-risk_score_table_multi_cox2

  #evaluate AUCs between 3-5 years.
  risk_score_table_multi_cox2_right$OS.time<-risk_score_table_multi_cox2_right$OS.time*30
  for_multi_ROC<-multi_ROC(time_vector = c(365*seq(3,5,1)),risk_score_table = risk_score_table_multi_cox2_right)

  #visualization for the ROC curves of multiple time points
  #plot ROC
  year3ROC<-for_multi_ROC[for_multi_ROC$Time_point==3*365,]
  year4ROC<-for_multi_ROC[for_multi_ROC$Time_point==4*365,]
  year5ROC<-for_multi_ROC[for_multi_ROC$Time_point==5*365,]
  TP_3year = year3ROC$True_positive
  FP_3year = year3ROC$False_positive
  TP_4year = year4ROC$True_positive
  FP_4year = year4ROC$False_positive
  TP_5year = year5ROC$True_positive
  FP_5year = year5ROC$False_positive
  time_ROC_df=data.frame(TP_3year,FP_3year,TP_4year,FP_4year,TP_5year,FP_5year)
  library(ggplot2)

  tiff(file='Figures/Figure_S3D.OS_AUC_right.tiff',width=15,height=15,units="cm",compression="lzw",bg="white",res=600)
  p=ggplot(data = time_ROC_df) +
	    geom_line(aes(x = FP_3year, y = TP_3year), size = 1, color = "#BC3C29FF") +
	      geom_line(aes(x = FP_4year, y = TP_4year), size = 1, color = "#0072B5FF") +
	        geom_line(aes(x = FP_5year, y = TP_5year), size = 1, color = "#E18727FF") +
		  geom_abline(slope = 1, intercept = 0, color = "grey", size = 1, linetype = 2) +
		    theme_bw() +
		      annotate("text",
			                  x = 0.75, y = 0.25, size = 4.5,
					             label = paste0("AUC at 3 years = ", sprintf("%.3f", year3ROC$AUC)), color = "#BC3C29FF"
					    ) +
  annotate("text",
	              x = 0.75, y = 0.15, size = 4.5,
		                 label = paste0("AUC at 4 years = ", sprintf("%.3f", year4ROC$AUC)), color = "#0072B5FF"
		        ) +
  annotate("text",
	              x = 0.75, y = 0.05, size = 4.5,
		                 label = paste0("AUC at 5 years = ", sprintf("%.3f", year5ROC$AUC)), color = "#E18727FF"
		        ) +
  labs(x = "False positive rate", y = "True positive rate") +
    theme(
	      axis.text = element_text(face = "bold", size = 11, color = "black"),
	          #axis.title.x = element_text(face = "bold", size = 14, color = "black", margin = margin(c(15, 0, 0, 0))),
	          #axis.title.y = element_text(face = "bold", size = 14, color = "black", margin = margin(c(0, 15, 0, 0)))
	        )

  print(p)
  dev.off()


  #K-M analysis
  AUC_max<-max(for_multi_ROC$AUC)
  AUC_max_time<-for_multi_ROC$Time_point[which(for_multi_ROC$AUC==AUC_max)]
  AUC_max_time<-AUC_max_time[!duplicated(AUC_max_time)]
  AUC_max_time<-AUC_max_time[length(AUC_max_time)]
  for_multi_ROC$Time_point<-as.factor(for_multi_ROC$Time_point)
  optimal_time_ROC_df<-for_multi_ROC[which(for_multi_ROC$Time_point==AUC_max_time),]
  cut.off<-optimal_time_ROC_df$Cut_values[which.max(optimal_time_ROC_df$True_positive - optimal_time_ROC_df$False_positive)]
  #cut.off<-0.005319551
  #cut.off<-median(risk_score_table_multi_cox2$total_risk_score)
  high_low<-risk_score_table_multi_cox2$total_risk_score>cut.off
  high_low[high_low==TRUE]<-'high'
  high_low[high_low==FALSE]<-'low'
  risk_score_table_multi_forKM<-cbind(risk_score_table_multi_cox2,high_low)

  risk_score_table_multi_forKM_right<-risk_score_table_multi_forKM
  #KM_plot 
  library(survminer)
  library(survival)
  risk_score_table_multi_forKM$OS[which(risk_score_table_multi_forKM$OS.time>AUC_max_time)]<-0
  risk_score_table_multi_forKM$OS.time[which(risk_score_table_multi_forKM$OS.time>AUC_max_time)]<-AUC_max_time
  fit_km<-survfit(Surv(OS.time,OS)~high_low,data = risk_score_table_multi_forKM)
  tiff(file='Figures/Figure_S3A.OS_riskscore_right.tiff',width=15,height=15,units="cm",compression="lzw",bg="white",res=300)
  p=ggsurvplot(fit_km,conf.int = F,pval = T,legend.title="Total risk score",
	                    legend.labs=c('high','low'),risk.table=T,
			                 #legend.labs=c(paste0('>',as.character(round(cut.off,2))),
			                 # paste0('<=',as.character(round(cut.off,2)))),risk.table=T,
			                 palette=c('red','blue'),surv.median.line='hv')
  print(p)
  dev.off()

  #table
  View(GSE103479_pd_use)
  GSE103479_pd_use$Age<-ifelse(as.numeric(GSE103479_pd_use$`age diagnosis:ch1`) >60,'>60','<=60')
  GSE103479_pd_use$Gender<-GSE103479_pd_use$`gender:ch1`
  GSE103479_pd_use$stage<-GSE103479_pd_use$`stage ii\\iii:ch1`
  GSE103479_pd_use$T<-GSE103479_pd_use$`tstage:ch1`
  GSE103479_pd_use$M<-GSE103479_pd_use$`mstage:ch1`
  GSE103479_pd_use$N<-GSE103479_pd_use$`nstage:ch1`
  GSE103479_pd_use$CMS<-GSE103479_pd_use$`cms subgroup:ch1`

  library(stringr)
  GSE103479_pd_use <- GSE103479_pd_use %>% na.omit()

  GSE103479_pd_use$stage<-factor(toupper(str_extract(GSE103479_pd_use$stage, "I+V*")))
  GSE103479_pd_use$T<-factor(str_extract(GSE103479_pd_use$T, "T[0-9]+"))
  GSE103479_pd_use$M <- factor(str_extract(GSE103479_pd_use$M, "M[0-9]*X*"))
  GSE103479_pd_use$N<-factor(str_extract(GSE103479_pd_use$N, "N[0-9]"))
  GSE103479_pd_use$CMS<-factor(GSE103479_pd_use$CMS)

  GSE103479_pd_use$OS<-factor(GSE103479_pd_use$OS,
			                                  levels=c(0,1),
							                              labels=c("Alive",  
											                                            "Death"))

  #save(GSE103479_pd_use,file = 'data/GSE103479_pd_use.Rdata')

  library(table1)
  table1(~ Gender + Age + CMS + stage + CMS  + OS + OS.time+
	          T  +   M  + N  | tumor_site, data=GSE103479_pd_use)

  GSE103479_pheno_CMS_markers_use_OS_use<-data.frame(GSE103479_pd_use$geo_accession,GSE103479_pd_use$OS.time,GSE103479_pd_use$CMS)
  colnames(GSE103479_pheno_CMS_markers_use_OS_use)<-c('sample','OS.time','CMS_group')
  rownames(GSE103479_pheno_CMS_markers_use_OS_use)<-GSE103479_pheno_CMS_markers_use_OS_use$sample

  GSE103479_pheno_CMS_markers_use_OS_use$sample<-NULL
  GSE103479_pheno_CMS_markers_use_OS_use<-GSE103479_pheno_CMS_markers_use_OS_use[GSE103479_pheno_CMS_markers_use_OS_use$CMS_group!='UNK',]

  GSE103479_pheno_CMS_markers_use_OS_use$CMS_group<-as.factor(GSE103479_pheno_CMS_markers_use_OS_use$CMS_group)

  mydata<-GSE103479_pheno_CMS_markers_use_OS_use %>% gather(key = "CMS_group",value = "OS.time") %>% na.omit() 
  my_comparisons <- list( c("CMS1", "CMS2"),c("CMS1", "CMS3"),c("CMS1", "CMS4"),c("CMS2", "CMS3"),c("CMS2", "CMS4"),c("CMS3", "CMS4"))

  library(cowplot)
  library(tidyverse)
  library(ggplot2)
  library(ggsci)
  library(ggpubr)
  library(tidyr)
  tiff(file='Figures/Figure_S1D.CMS_OS_GSE103479.tiff',width=18,height=15,units="cm",compression="lzw",bg="white",res=600)
  p<-ggboxplot(mydata,x='CMS_group',y='OS.time',color='CMS_group',palette = 'jama',add = 'jitter')
  p+stat_compare_means(comparisons = my_comparisons)
  dev.off()


  #reply to reviewers
  #1. validation the expression of markers
  #add GSE103479
  load('data/GSE103479_gene_survival_data_left.Rdata')
  load('data/GSE103479_gene_survival_data_right.Rdata')
  load('data/GSE103479_exprSet_anno.Rdata')

  GSE103479_sample<-c(GSE103479_gene_survival_data_left$sample,GSE103479_gene_survival_data_right$sample)
  GSE103479_expr_target<-GSE103479_exprSet_anno[c(left_markers,right_markers),GSE103479_sample]

  marker_in=array()
  for (marker in c(left_markers,right_markers)){
	    for (sam in GSE103479_sample){
		        
		        if(sam %in% GSE103479_gene_survival_data_left$sample){
				      left_marker_in<-data.frame(marker,GSE103479_expr_target[marker,sam],sam,'LCC')
        colnames(left_marker_in)<-c('Markers','expression','sample','Group')
	      marker_in<-rbind(marker_in,left_marker_in)
	    }
     
      if(sam %in% GSE103479_gene_survival_data_right$sample){
	            right_marker_in<-data.frame(marker,GSE103479_expr_target[marker,sam],sam,'RCC')
            colnames(right_marker_in)<-c('Markers','expression','sample','Group')
	          marker_in<-rbind(marker_in,right_marker_in)
	        }
          

        }
  }


  library(dplyr)
  marker_in<-marker_in %>% na.omit()


  library(cowplot)
  library(tidyverse)
  library(ggplot2)
  library(ggsci)
  library(ggpubr)
  library(dplyr)

  p <- ggboxplot(marker_in, x = "Markers", y = "expression",
		                color = 'Group', palette = "nejm",add = "jitter")+
  stat_compare_means(aes(group = Group),method = "t.test",label = "p.format")+rotate_x_text(angle = 45)
tiff(file='Figures/GSE103479.left_right_marker.tiff',width=25,height=15,units="cm",compression="lzw",bg="white",res=600)
print(p)
dev.off()

load('data/GSE103479_pd_use.Rdata')
common_sample<-intersect(colnames(GSE103479_exprSet_anno),rownames(GSE103479_pd_use))
GSE103479_marker_expr<-t(GSE103479_exprSet_anno[c(left_markers,right_markers),common_sample])
GSE103479_marker_expr<-as.data.frame(GSE103479_marker_expr)
GSE103479_marker_expr$sample<-rownames(GSE103479_marker_expr)

GSE103479_pd_use$sample<-GSE103479_pd_use$geo_accession
GSE103479_marker_expr_CMS<-merge(GSE103479_pd_use[,c('sample','CMS')],GSE103479_marker_expr,by='sample')
rownames(GSE103479_marker_expr_CMS)<-GSE103479_marker_expr_CMS$sample

save(GSE103479_marker_expr_CMS,file='data/GSE103479_marker_expr_CMS.Rdata')

GSE103479_marker_expr_CMS<-GSE103479_marker_expr_CMS[GSE103479_marker_expr_CMS$CMS!='UNK',]
marker_in=array()
for (marker in c(left_markers,right_markers)){
	  for (sample in rownames(GSE103479_marker_expr_CMS)){
		      uni_marker_CMS<-data.frame(marker,GSE103479_marker_expr_CMS[sample,marker],sample,GSE103479_marker_expr_CMS[sample,'CMS'])
    colnames(uni_marker_CMS)<-c('Markers','expression','sample','CMS_group')
        marker_in<-rbind(marker_in,uni_marker_CMS)
      }
}

library(dplyr)
marker_in<-marker_in %>% na.omit()


library(cowplot)
library(tidyverse)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(dplyr)

p <- ggboxplot(marker_in, x = "Markers", y = "expression",
	                      color = 'CMS_group', palette = "nejm",add = "jitter")+
  stat_compare_means(aes(group = CMS_group),label = "p.signif")+rotate_x_text(angle = 45)
tiff(file='Figures/Figure_add2.CMS_group_expression.tiff',width=25,height=15,units="cm",compression="lzw",bg="white",res=600)
print(p)
dev.off()

#validation 2:https://ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144735
#https://pubmed.ncbi.nlm.nih.gov/32451460/
#https://sci-hub.st/10.1038/s41588-020-0636-z
library(data.table)
GSE144735_anno<-read.table('data/GSE144735_processed_KUL3_CRC_10X_annotation.txt',sep = '\t',header = T)
rownames(GSE144735_anno)<-GSE144735_anno$Index
GSE144735_info<-read.table('data/E-MTAB-8410.sdrf.txt',sep='\t',header = T)
#GSE144735_expr<-fread('data/GSE144735_processed_KUL3_CRC_10X_raw_UMI_count_matrix.txt.gz',sep='\t',header = T)
GSE144735_expr_Tumor<-as.data.frame(GSE144735_expr)[,c('Index',GSE144735_expr_Tumor_cells)]
#select the tumor sample
GSE144735_tumor_sample<-unique(GSE144735_anno[GSE144735_anno$Class=='Tumor',]$Sample)
rownames(GSE144735_expr_Tumor)<-GSE144735_expr_Tumor$Index
GSE144735_expr_Tumor$Index<-NULL
save(GSE144735_expr_Tumor,file='GSE144735_expr_Tumor.Rdata')

load('data/GSE144735_expr_Tumor.Rdata')
GSE144735_expr_Tumor_cells<-colnames(GSE144735_expr)[substr(colnames(GSE144735_expr),7,7)=='T']
GSE144735_expr_Tumor<-as.data.frame(GSE144735_expr)[,GSE144735_expr_Tumor_cells]

GSE144735_tumor<-CreateSeuratObject(counts = GSE144735_expr_Tumor, min.cells = 3, project = "GSE144735")
GSE144735_tumor[["percent.mt"]] <- PercentageFeatureSet(GSE144735_tumor, pattern = "^MT-")

GSE144735_T2L<-read.table('data/KUL_info.txt',sep='\t',header = T)
GSE144735_anno_tumor<-GSE144735_anno[GSE144735_anno$Class=='Tumor',]
GSE144735_anno_tumor_Location<-merge(GSE144735_anno_tumor,GSE144735_T2L,by='Sample')
rownames(GSE144735_anno_tumor_Location)<-GSE144735_anno_tumor_Location$Index
GSE144735_anno_tumor_Location<-GSE144735_anno_tumor_Location[rownames(GSE144735_tumor@meta.data),]
GSE144735_tumor@meta.data$cell_type<-GSE144735_anno_tumor_Location$Cell_type
GSE144735_tumor@meta.data$cell_subtype<-GSE144735_anno_tumor_Location$Cell_subtype

#right: ascending colon  caecum
#left:splenic descending sigmoid junction rectum
GSE144735_tumor@meta.data$Group<-ifelse(GSE144735_anno_tumor_Location$Location=='ascending colon' |
					                                          GSE144735_anno_tumor_Location$Location=='caecum','RCC','LCC')

head(GSE144735_tumor@meta.data)
table(GSE144735_tumor@meta.data$Group)

save(GSE144735_tumor,file='data/GSE144735_tumor.Rdata')

load('data/GSE144735_tumor.Rdata')

Idents(object = GSE144735_tumor)<-"cell_type"
T_cell.markers <- FindMarkers(GSE144735_tumor, ident.1="LCC",ident.2="RCC", subset.ident = c("T cells"),
			                                      #logfc.threshold = 0.585, 
			                                      logfc.threshold = 0,
							                                      group.by = 'Group',
							                                      test.use = "roc")

Epithelial_cells.markers <- FindMarkers(GSE144735_tumor, ident.1="LCC",ident.2="RCC", subset.ident = c("Epithelial cells"),
					                              #logfc.threshold = 0.585, 
					                              logfc.threshold = 0,
								                                    group.by = 'Group',
								                                    test.use = "roc")

B_cells.markers <- FindMarkers(GSE144735_tumor, ident.1="LCC",ident.2="RCC", subset.ident = c("B cells"),
			                                      #logfc.threshold = 0.585,
			                                      logfc.threshold = 0,
							                                     group.by = 'Group',
							                                     test.use = "roc")

Stromal_cells.markers <- FindMarkers(GSE144735_tumor, ident.1="LCC",ident.2="RCC", subset.ident = c("Stromal cells"),
				                                    #logfc.threshold = 0.585,
				                                    logfc.threshold = 0,
								                                   group.by = 'Group',
								                                   test.use = "roc")

T_cell.markers$celltype='T Cells'
T_cell.markers$markers<-rownames(T_cell.markers)
Epithelial_cells.markers$celltype='Epithelial cells'
Epithelial_cells.markers$markers<-rownames(Epithelial_cells.markers)

B_cells.markers$celltype<-'B cells'
B_cells.markers$markers<-rownames(B_cells.markers)
Stromal_cells.markers$celltype<-'Stromal cells'
Stromal_cells.markers$markers<-rownames(Stromal_cells.markers)

left_markers<-c("FOSB","RPL35","REG1A", "TESC","C11orf96")
right_markers<-c("HSPA1A","CD69","GDF15","LGALS2")

markers <- c(left_markers,right_markers)
#the expression differences of the candidate markers
GSE144735_marker_expression<-rbind(T_cell.markers[markers,],B_cells.markers[markers,],
				                                      Epithelial_cells.markers[markers,],Stromal_cells.markers[markers,])

save(T_cell.markers,B_cells.markers,Epithelial_cells.markers,Stromal_cells.markers,GSE144735_marker_expression,file='data/GSE144735_marker_expression.Rdata')





save(T_cell.markers,file='data/T_cell.markers.Rdata')
save(B_cells.markers,file='data/B_cells.markers.Rdata')
save(Epithelial_cells.markers,file='data/Epithelial_cells.markers.Rdata')

load('data/T_cell.markers.Rdata')
load('data/B_cells.markers.Rdata')
load('data/Epithelial_cells.markers.Rdata')



dir.create('review/')
dir.create('review/Tables')
write.table(T_cell.markers,file='review/Tables/Table_AS1.tcell_DEG.xls',quote = F,sep = '\t')
write.table(B_cells.markers,file='review/Tables/Table_AS1.bcell_DEG.xls',quote = F,sep = '\t')
write.table(Epithelial_cells.markers,file='review/Tables/Table_AS1.Epithelial_cells_DEG.xls',quote = F,sep = '\t')

#the distribution
dir.create('review/Figures')
GSE144735_tumor@meta.data$active.ident<-GSE144735_tumor@active.ident
p=VlnPlot(object = GSE144735_tumor, features = left_markers,group.by = "active.ident",ncol=3,pt.size = 0.5)
tiff(file='review/Figures/markers_expression_left.tiff',width=20,height=15,units="cm",compression="lzw",bg="white",res=600)
print(p)
dev.off()

p=VlnPlot(object = GSE144735_tumor, features = right_markers,group.by = "active.ident",ncol=3,pt.size = 0.5)
tiff(file='review/Figures/markers_expression_right.tiff',width=20,height=15,units="cm",compression="lzw",bg="white",res=600)
print(p)
dev.off()


load('data/GSE144735_marker_expression.Rdata')
GSE144735_marker_expression <- GSE144735_marker_expression %>% na.omit()
write.table(GSE144735_marker_expression,file = 'review/Tables/GSE144735_marker_expression.xls',sep = '\t',quote = F)

