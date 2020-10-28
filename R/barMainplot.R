#
# j = "Phylum"
# j = "Class"
# j = "Order"
# j =  "Family"
# j = "Genus"
#
#
# ##k 是否过滤或者设置过滤值为多少
# k= 0.01
#
# ##重复数量
# rep = 6
#

#----因为微生物群落物种数量太多，所以我在这里要设置将低丰度的定义为other，为了突出主要的功能
# ps = readRDS("./ps_liu.rds")
# ps
# otu = NA
# tax = NA
# map = NA
# ps = ps
# j = "Phylum"
#
# rep = 6
# axis_ord = NA
# label = TRUE
# sd = FALSE


# ps = readRDS("./ps_liu.rds")
# result = barMainplot(ps = ps,j = "Phylum",rep = 6,axis_ord = NA,label = FALSE ,sd = FALSE,Top = 10)
# result[[1]]
# result[[2]]
#
# result[[3]]
# # 提取丰度最高的前十个物种作展示：
#
# otu =NULL
# tax = NULL
# map = NULL
#
#
# ps = ps
# j = "Phylum"
# rep = 6
# axis_ord = NA
# group = "Group"
# label = FALSE
# sd = FALSE
# Top = 10

# # ps = readRDS("./ps16S.rds")
# otu =NULL
# tax = NULL
# map = NULL
# ps = ps_rela
# ps
# j = "group"
# rep = 6
# group = "Group"
# axis_ord = NULL
# label = FALSE
# sd = FALSE
# Top = 10
#
# tran = TRUE

library(phyloseq)
library(tidyverse)
library(vegan)
library(reshape2)
library("plyr")
library(ggalluvial)
library(ggplot2)


barMainplot = function(otu = NULL,tax = NULL,map = NULL,tree = NULL ,ps = NULL,group  = "Group",
                       j = "Phylum",rep = 6,axis_ord = NULL,label = TRUE ,sd = FALSE,Top = 10,tran = TRUE){


  if (is.null(axis_ord)) {
    axis_order = NA
  }else{
    axis_order = strsplit(basename(axis_ord), "-")[[1]]
  }

  ps = inputMicro(otu,tax,map,tree,ps,group  = group)

  psdata = ps %>%
    tax_glom(taxrank = j)

  # transform to relative abundance
  if (tran == TRUE) {
    psdata = psdata%>%
      transform_sample_counts(function(x) {x/sum(x)} )
  }




  otu = otu_table(psdata)
  tax = tax_table(psdata)



for (i in 1:dim(tax)[1]) {
  if (row.names(tax)[i] %in% names(sort(rowSums(otu), decreasing = TRUE)[1:Top])) {

    tax[i,j] =tax[i,j]
  } else {
    tax[i,j]= "others"
  }
}
  tax_table(psdata)= tax

  Taxonomies <- psdata %>% # Transform to rel. abundance
    psmelt()


  # 这里我们看到有很过属，因此颜色上就会出现不能很好区分的现象
  colbar <- dim(unique(dplyr::select(Taxonomies, one_of(j))))[1]


  colors = colorRampPalette(c( "#CBD588", "#599861", "orange","#DA5724", "#508578", "#CD9BCD",
                               "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
                               "#8569D5", "#5E738F","#D1A33D", "#8A7C64","black"))(colbar)



  Taxonomies$Abundance = Taxonomies$Abundance * 100
  # Taxonomies$Abundance = Taxonomies$Abundance/sum(Taxonomies$Abundance)

  Taxonomies$Abundance = Taxonomies$Abundance/rep


  #按照分组求均值
  colnames(Taxonomies) <- gsub(j,"aa",colnames(Taxonomies))
  by_cyl <- dplyr::group_by(Taxonomies, Group,aa)
  zhnagxu2 = dplyr::summarise(by_cyl, sum(Abundance), sd(Abundance))


  iris_groups<- dplyr::group_by(Taxonomies, aa)
  cc<- dplyr::summarise(iris_groups, sum(Abundance))
  head(cc)
  colnames(cc)= c("aa","allsum")
  cc<- dplyr::arrange(cc, desc(allsum))

  ##使用这个属的因子对下面数据进行排序

  colnames(zhnagxu2) <- c("group","aa","Abundance","sd")
  zhnagxu2$aa = factor(zhnagxu2$aa,order = T,levels = cc$aa)

  zhnagxu3 = zhnagxu2

  ##制作标签坐标，标签位于顶端
  # Taxonomies_x = ddply(zhnagxu3,"group", summarize, label_y = cumsum(Abundance))
  # head(Taxonomies_x )
  #标签位于中部
  # Taxonomies_x1 = ddply(zhnagxu3,"group", transform, label_y = cumsum(Abundance) - 0.5*Abundance)
  Taxonomies_x = plyr::ddply(zhnagxu3,"group", summarize,label_sd = cumsum(Abundance),label_y = cumsum(Abundance) - 0.5*Abundance)
  # Taxonomies_x$label_y =
  Taxonomies_x = cbind(as.data.frame(zhnagxu3),as.data.frame(Taxonomies_x)[,-1])

  head(Taxonomies_x,6 )
  Taxonomies_x$label = Taxonomies_x$aa
  # #使用循环将堆叠柱状图柱子比较窄的别写标签，仅仅宽柱子写上标签
  # for(i in 1:nrow(Taxonomies_x)){
  #   if(Taxonomies_x[i,3] > 3){
  #     Taxonomies_x[i,5] = Taxonomies_x[i,5]
  #   }else{
  #     Taxonomies_x[i,5] = NA
  #   }
  # }



  Taxonomies_x$aa = factor(Taxonomies_x$aa,order = T,levels = c(as.character(cc$aa)))



##普通柱状图


  p4 <- ggplot(Taxonomies_x , aes(x =  group, y = Abundance, fill = aa, order = aa)) +
    geom_bar(stat = "identity",width = 0.5,color = "black") +
    scale_fill_manual(values = colors) +
    theme(axis.title.x = element_blank()) +
    theme(legend.text=element_text(size=6)) +
    scale_y_continuous(name = "Relative abundance (%)")

  p4
  if (is.na(axis_order)) {
    p4 = p4
  }else{
    p4 = p4 +scale_x_discrete(limits = axis_order)
  }


  if (sd == TRUE) {
    p4 =  p4 +
      geom_errorbar(aes(ymin=label_sd-sd, ymax=label_sd +sd), width=.2)
  }

  if (label == TRUE) {
    p4 = p4 +
      geom_text(aes(y = label_y, label = label ),size = 4,fontface = "bold.italic")
  }

  # print(p4)

  # install.packages("ggalluvial")
  p4 = p4+theme_bw()+
    scale_y_continuous(expand = c(0,0))+

    theme(

      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      text=element_text(face = "bold"),
      plot.title = element_text(vjust = -8.5,hjust = 0.1),
      axis.title.y =element_text(size = 20,face = "bold",colour = "black"),
      axis.title.x =element_text(size = 24,face = "bold",colour = "black"),
      axis.text = element_text(size = 20,face = "bold"),
      axis.text.x = element_text(colour = "black",size = 14,),
      axis.text.y = element_text(colour = "black",size = 14),

      legend.text = element_text(size = 15)
      #legend.position = "none"#是否删除图例

    )
  p4
  map = as.data.frame(sample_data(ps))
  if (length(unique(map$Group))>3){	p4=p4+theme(axis.text.x=element_text(angle=45,vjust=1, hjust=1))}

  #-------------冲击图
  head(Taxonomies_x )
  cs = Taxonomies_x $aa

  # head(cs)
  # as.factor(Taxonomies_x $Genus)
  # cs = as.character(Taxonomies_x $Genus)
  # cs1 = as.factor(cs)
  cs1 = cs
  #提取真正的因子的数量
  lengthfactor = length(levels(cs1))
  #提取每个因子对应的数量
  cs3 = summary(as.factor(cs1))
  cs4 = as.data.frame(cs3)
  cs4$id = row.names(cs4)
  #对因子进行排序
  df_arrange<- arrange(cs4, id)
  #对Taxonomies_x 对应的列进行排序
  Taxonomies_x1<- arrange(Taxonomies_x , aa)
  head(Taxonomies_x1)
  #构建flow的映射列Taxonomies_x
  Taxonomies_x1$ID = factor(rep(c(1:lengthfactor), cs4$cs3))

  #colour = "black",size = 2,,aes(color = "black",size = 0.8)
  head(Taxonomies_x1)
  p3 = ggplot(Taxonomies_x1,
              aes(x = group, stratum = aa, alluvium = ID,
                  weight = Abundance,
                  fill = aa, label = aa)) +
    geom_flow(stat = "alluvium", lode.guidance = "rightleft",
              color = "black",size = 0.2,width = 0.3,alpha = .2) +
    geom_bar(width = 0.45)+
    geom_stratum(width = 0.45,size = 0.2) +
    #geom_text(stat = "stratum", size = 3,family="Times New Roman",fontface = "bold.italic") +
    #theme(legend.position = "none") +
    scale_fill_manual(values = colors)+
    #ggtitle("fow_plot")+
    # scale_x_discrete(limits = axis_order)+
    # geom_text(aes(y = label_y, label = label ),size = 4,fontface = "bold.italic")+
    labs(x="group",
         y="Relative abundance (%)",
         title="")
  p3


  if (is.na(axis_order)) {
    p3 = p3
  }else{
    p3 = p3 +scale_x_discrete(limits = axis_order)
  }
  # p3
  if (label == TRUE) {
    p3 = p3 +
      geom_text(aes(y = label_y, label = label ),size = 4,fontface = "bold.italic")
  }

  if (sd == TRUE) {
    p3 =  p3 +
      geom_errorbar(aes(ymin=label_sd-sd, ymax=label_sd +sd), width=.2)
  }
  p3 =p3+theme_bw()+
    scale_y_continuous(expand = c(0,0))+
    #geom_hline(aes(yintercept=0), colour="black", linetype=2) +
    #geom_vline(aes(xintercept=0), colour="black", linetype="dashed") +
    #scale_fill_manual(values = mi, guide = guide_legend(title = NULL))+
    theme(

      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      text=element_text(face = "bold"),
      plot.title = element_text(vjust = -8.5,hjust = 0.1),
      axis.title.y =element_text(size = 20,face = "bold",colour = "black"),
      axis.title.x =element_text(size = 24,face = "bold",colour = "black"),
      axis.text = element_text(size = 20,face = "bold"),
      axis.text.x = element_text(colour = "black",size = 14),
      axis.text.y = element_text(colour = "black",size = 14),

      legend.text = element_text(size = 15,face = "bold.italic")
      #legend.position = "none"#是否删除图例

    )
  p3
  if (length(unique(map$Group))>3){	p3=p3+theme(axis.text.x=element_text(angle=45,vjust=1, hjust=1))}

  return(list(p4,Taxonomies_x,p3))
}





