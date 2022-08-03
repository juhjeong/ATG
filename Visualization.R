library(dplyr)
library(patchwork)
library(Seurat)
library(RColorBrewer)
library(colorRamps)
library(ggplot2)
library(Nebulosa)
library(viridis)

DefaultAssay(Pan02_integrated)<-"RNA"
DefaultAssay(Pan02_CAF)<-"RNA"
DefaultAssay(Pan02_TNK)<-"RNA"
DefaultAssay(Pan02_CD8)<-"RNA"
DefaultAssay(Pan02_Neut)<-"RNA"
DefaultAssay(Elyada_CAF)<-"RNA"
DefaultAssay(Dominguez_CAF)<-"RNA"

#Feature plot
FeaturePlot(Pan02_integrated, features = "Nrp1",pt.size =0.8, cols = c("ivory2","red"))
FeaturePlot(Pan02_integrated, features = "Nrp2",pt.size =0.8, cols = c("ivory2","red"))
FeaturePlot(Pan02_CAF, features = "Inhba",pt.size =0.8, cols =brewer.pal(9,"Reds"))
FeaturePlot(Pan02_TNK, features = "Cd8b1",pt.size =0.8, cols =brewer.pal(9,"Reds"), split.by="orig.ident")
plot_density(Pan02_TNK, "Ccr7")
plot_density(Pan02_TNK, "Sell")
plot_density(Pan02_TNK, "Cd44")
plot_density(Pan02_TNK, "Prf1")
plot_density(Pan02_TNK, "Ifng")
plot_density(Pan02_TNK, "Gzmb")

#Dot plot
DotPlot(Pan02_integrated, features = c("Pdgfra","Pdgfrb","Acta2","Col1a1","Cd3d","Trac","Cd8a","Nkg7","Ms4a1","Cd79a","Igkc","Cd19","C1qa","Adgre1","Csf1r","Csf3r","Cxcr2","S100a8","Flt3","Xcr1","Irf8","Siglech","Fcer1a","Tpsb2","Ms4a2","Krt19","Krt18","Krt8"),cols = c("ivory2","red"))
DotPlot(Pan02_CAF, features = c("Flt1","Kdr","Flt4","Nrp1","Nrp2"),cols = c("ivory2","red"))
DotPlot(Pan02_CAF, features = c("Thbd","Cd74"),cols = c("ivory2","red"))
DotPlot(Pan02_CAF, features = c("Spp1","Acta2","Ecscr","Tnc","Rcn3","Col12a1","Thbs2","H19","Cilp","Sdc1","Col1a1","Tgfb1","Sfrp1","Cthrc1","Thy1","Tagln","Crlf1","Cxcl14","Serpine2","Igfbp3"),cols = c("blue","yellow"))
DotPlot(Pan02_CAF, features = c("Ogn","Mfap5","Rarres2","Islr","Fxyd1","Thbd","Prss23","Sema3c","Clec3b","Tnxb","Ly6a","Plpp3","Efemp1","Ly6c1","Ifi27l2a","Dpt","Svep1","Pcolce2","Ccl11","C3"),cols = c("blue","yellow"))
DotPlot(Pan02_CAF, features = c("Cd74","H2-DMb1","H2-Ab1","Stmn1","Msln","Slpi","Saa3","Serpinb2","Dmkn","Lgals7","Upk3b","Nkain4","Calca","Slc9a3r1","Clu","Ptgis","Cxadr","Ezr","Cav1","Angptl7"),cols = c("blue","yellow"))
DotPlot(Pan02_CAF, features = c("Gstm1","Sema3c","Prss23","Fxyd1","Mfap5","Ogn","Ly6c1","Efemp1","Cd248","Ly6a","Pi16","Pcolce2","Ackr3","Anxa3","Scara3","Col14a1","Plpp3","Lrrn4cl","Clec3b","Adgrd1"),cols = c("blue","yellow"))
DotPlot(Pan02_CAF, features = c("Col12a1","Ecscr","Acta2","Ctxn1","Col8a1","Rbp1","Olfml3","F2r","Actn1","Cdh11","Mdk","Fzd1","Palld","Pmepa1","Csrp2","Marcksl1","Tgfb3","Csrp1","Timp3","Dkk3"),cols = c("blue","yellow"))
DotPlot(Pan02_CAF, features = c("Nfkbia","Cxcl1","Fosb","Ifrd1","Hk2","Ier3","Csrnp1","Gm12840","Nr4a1","Nfkbiz","Tnfaip6","Zfp36","Serpine1","Btg2","Il6","Klf4","Ccnl1","Ccl2","Has2","Atf3"),cols = c("blue","yellow"))
DotPlot(Pan02_Neut, features = "Cd274",cols = c("ivory2","red"))
DotPlot(Pan02_TNK, features = c("Cd3e","Cd8b1","Cd4","Nkg7","Foxp3","Cd7","Rorc","Mki67","Mcm2","Csf3r","mt-Co2"),cols = c("ivory2","red"))

#Violin plot
VlnPlot(Pan02_TNK, features = "Gzmb", group.by = "orig.ident")
VlnPlot(Pan02_TNK, features = "Ifng", group.by = "orig.ident")
VlnPlot(Pan02_TNK, features = "Cd28", group.by = "orig.ident")

#Heatmap(Average expression)
orig.levels <- levels(Pan02_CAF)
Idents(Pan02_CAF) <- gsub(pattern = " ", replacement = " ", x = Idents(Pan02_CAF))
orig.levels <- gsub(pattern = " ", replacement = " ", x = orig.levels)
levels(Pan02_CAF) <- orig.levels
cluster.averages_Pan02_CAF <- AverageExpression(Pan02_CAF, return.seurat = TRUE)
my_cols <- c("CAF-1"="#04B4EC","CAF-2"="#FF514C","CAF-3"="#F48FEB","CAF-4"="#FF9F25","CAF-5"="#a433f5","CAF-6"="#6ACB61")
my_cols2 <- my_cols[order(as.integer(names(my_cols)))]
mapal <- colorRampPalette(RColorBrewer::brewer.pal(9,"Reds"))(256)
DoHeatmap( cluster.averages_Pan02_CAF,features = c("Plac8","Bst2","Apol9a","Gbp3","Gbp2","Apol9b","Ly6c1","Psmb8","Stat1","Igtp","Ogn","Prelp","Wnt4","Mfap5","Rarres2","Islr","Angpt4","Fzd2","Ecrg4","Spry1","Col6a3","Scara5","Olfml3","Dclk1","Lcn2","Gda","Fmo1","Lbp","Orm1","Kng2","Egr1","Btg2","Fos","Nr4a1","Jun","Dusp1","Zfp36","Trib1","Ppp1r15a","Ier3","Wfdc2","Col12a1","Rcn3","Cpe","Aqp1","Cd24a","Onecut2","Ecscr","Cxcl2","Cdh13","Basp1","Cyba","Pltp","Tnc","Atp1b1","Stxbp6","Mtap","Arhgef28","Cops5","Mrpl15"
), size = 3, draw.lines = FALSE,  label = F, group.colors  = my_cols2)+ scale_fill_gradientn(colours = mapal) + 
  theme(text = element_text(size = 15))

orig.levels <- levels(Elyada_CAF)
Idents(Elyada_CAF) <- gsub(pattern = " ", replacement = " ", x = Idents(Elyada_CAF))
orig.levels <- gsub(pattern = " ", replacement = " ", x = orig.levels)
levels(Elyada_CAF) <- orig.levels
cluster.averages_Elyada_CAF <- AverageExpression(Elyada_CAF, return.seurat = TRUE)
DoHeatmap(cluster.averages_Elyada_CAF,features = c("Col8a1","Tagln","Clec3b","Ccl7","H2-Ab1","Cd74"), size = 3, draw.lines = FALSE,  label = F)+ scale_fill_viridis()

orig.levels <- levels(Dominguez_CAF)
Idents(Dominguez_CAF) <- gsub(pattern = " ", replacement = " ", x = Idents(Dominguez_CAF))
orig.levels <- gsub(pattern = " ", replacement = " ", x = orig.levels)
levels(Dominguez_CAF) <- orig.levels
cluster.averages_Dominguez_CAF <- AverageExpression(Dominguez_CAF, return.seurat = TRUE)
DoHeatmap(cluster.averages_Dominguez_CAF,features = c("Ly6c1","Pi16","Serpine2","Gpx3","Tagln","Cthrc1","Cxcl1","Has1"), size = 3, draw.lines = FALSE,  label = F)+ scale_fill_viridis()

orig.levels <- levels(Pan02_CD8)
Idents(Pan02_CD8) <- gsub(pattern = " ", replacement = " ", x = Idents(Pan02_CD8))
orig.levels <- gsub(pattern = " ", replacement = " ", x = orig.levels)
levels(Pan02_CD8) <- orig.levels
cluster.averages_Pan02_CD8 <- AverageExpression(Pan02_CD8, return.seurat = TRUE)
DoHeatmap(cluster.averages_Pan02_CD8,features = c("Tnfrsf9","Prf1","Gzmb","Ifng","Cd160","Xcl1"), size = 3, draw.lines = FALSE,  label = F)+ scale_fill_viridis()

#Feature plot (Average expression)
CAF2DEG_list <- list(c("Ogn","Prelp","Wnt4","Mfap5","Rarres2","Islr","Angpt4","Fzd2","Ecrg4","Spry1","Fmod","Inmt","Plk2","Inhba","Omd","Phyhd1","Ltbp4","Fxyd1","Ckb","Thbd","Rdh10","Svil","Lbp","Prss23","Sema3c","Dhrs3","Galk2","Rabgap1l","Sult1a1","Ptprf"
))
CAF_scored <- AddModuleScore(object=Pan02_integrated, features = CAF2DEG_list, name = "CAF2_score", replace=TRUE)
FeaturePlot(CAF_scored, features = "CAF2_score1",pt.size =0.8, cols = matlab.like(8))


#Dimplot
DimPlot(Pan02_integrated, reduction = "umap", label = TRUE, pt.size = 0.4)+NoLegend()
my_cols <- c("CAF-1"="#04B4EC","CAF-2"="#FF514C","CAF-3"="#F48FEB","CAF-4"="#FF9F25","CAF-5"="#a433f5","CAF-6"="#6ACB61")
my_cols2 <- my_cols[order(as.integer(names(my_cols)))]
DimPlot(Pan02_CAF, reduction = "umap", label = TRUE, pt.size = 0.5, cols = my_cols2, split.by = "orig.ident")+NoLegend()
DimPlot(Pan02_Neut, reduction = "umap", label = TRUE, pt.size = 0.5, split.by="orig.ident")+NoLegend()
DimPlot(Pan02_TNK, reduction = "umap", label = TRUE, pt.size = 0.5)+NoLegend()

