library(ggsignif)
library(ggpubr)
library(ggpmisc)
library(Hmisc)
library(clusterProfiler)
library(enrichplot)
library(ggraph)
library(igraph)
library(tidyverse)


sum.id.nkt<-table(paste(sc.data.nkt$sampleID,sc.data.nkt$cluster_annot_nkt,sep = "/"))

sc.data.nkt$cluster_annot_nkt<-factor(sc.data.nkt$cluster_annot_nkt)
sc.data.nkt$sampleID<-factor(sc.data.nkt$sampleID)

sample.sum<-data.frame(matrix(0,
                              nrow = length(levels(sc.data.nkt$cluster_annot_nkt)),
                              ncol = length(levels(sc.data$sampleID))
                              )
                       )
rownames(sample.sum)<-levels(sc.data.nkt$cluster_annot_nkt)
colnames(sample.sum)<-levels(sc.data$sampleID)

sum.id.nkt.split<-tstrsplit(names(sum.id.nkt),split = "/")
for(i in 1:length(sum.id.nkt))
{
  sample.sum[sum.id.nkt.split[[2]][i],sum.id.nkt.split[[1]][i]]<-sum.id.nkt[i]
}

sum.id.cell<-table(paste(sc.data$sampleID,sc.data$cluster_tree_annot_merge,sep = "/"))

sc.data$sampleID<-factor(sc.data$sampleID)

sample.sum2<-data.frame(matrix(0,
                              nrow = length(levels(sc.data$cluster_tree_annot_merge)),
                              ncol = length(levels(sc.data$sampleID))
                              )
                       )

rownames(sample.sum2)<-levels(sc.data$cluster_tree_annot_merge)
colnames(sample.sum2)<-levels(sc.data$sampleID)

sum.id.cell.split<-tstrsplit(names(sum.id.cell),split="/")
for(i in 1:length(sum.id.cell))
{
  sample.sum2[sum.id.cell.split[[2]][i],sum.id.cell.split[[1]][i]]<-sum.id.cell[i]
}

sample.sum<-rbind(sample.sum,sample.sum2)
a<-apply(sample.sum,2,sum)

sample.sum<-t(sample.sum)
sample.sum<-as.data.frame(sample.sum)
sample.sum$Sample.name<-rownames(sample.sum)
sample.sum$cell_count<-a[match(sample.sum$Sample.name,names(a))]

rm(sample.sum2)
rm(sum.id.cell)
rm(sum.id.cell.split)
rm(sum.id.nkt)
rm(sum.id.nkt.split)


meltP<-data.frame(a)
meltP$sample<-rownames(meltP)

ggplot(meltP)+geom_col(aes(x=sample,y=a))

sc.meta.cell<-merge(sc.meta.pbmc,sample.sum)
sc.meta.cell$X<-NULL
sc.meta.cell<-sc.meta.cell[sc.meta.cell$cell_count>500,]


colnames(sc.meta.cell)<-gsub("characteristics..","",colnames(sc.meta.cell),fixed = T)

sc.meta.cell$`Sampling.day..Days.after.symptom.onset.`<-as.numeric(sc.meta.cell$`Sampling.day..Days.after.symptom.onset.`)

sc.meta.cell$Lym<-sc.meta.cell$Lymphocytes..G.L.
sc.meta.cell$Lym[grep("[^0-9.]+",sc.meta.cell$Lym)]<-NA
sc.meta.cell$Lym<-as.numeric(sc.meta.cell$Lym)

sc.meta.cell$Leu<-sc.meta.cell$Leukocytes..G.L.
sc.meta.cell$Leu[grep("[^0-9.]+",sc.meta.cell$Leu)]<-NA
sc.meta.cell$Leu<-as.numeric(sc.meta.cell$Leu)

sc.meta.cell$Neu<-sc.meta.cell$Neutrophils..G.L.
sc.meta.cell$Neu[grep("[^0-9.]+",sc.meta.cell$Neu)]<-NA
sc.meta.cell$Neu<-as.numeric(sc.meta.cell$Neu)

sc.meta.cell$Mon<-sc.meta.cell$Leu-sc.meta.cell$Lym-sc.meta.cell$Neu

sc.meta.cell$Lym.est<-apply(sc.meta.cell[,grep("^T_cell|^B_cell|^NK_cell|^Plasma_cell",colnames(sc.meta.cell))],1,sum)
sc.meta.cell$Mon.est<-apply(sc.meta.cell[,grep("^Monocyte",colnames(sc.meta.cell)),drop=F],1,sum)
sc.meta.cell$Tcell.est<-apply(sc.meta.cell[,grep("^T_cell|NK_cell",colnames(sc.meta.cell))],1,sum)

sc.meta.cell$Outcome<-factor(sc.meta.cell$Outcome,levels = c("control","discharged","deceased"))

sc.meta.cell$NKT.perc<-sc.meta.cell$T_cell_NKT/sc.meta.cell$Tcell.est
sc.meta.cell<-sc.meta.cell[sc.meta.cell$T_cell_NKT>=20|sc.meta.cell$NKT.perc>0.01,]

sc.meta$Sev_Outcome<-paste(sc.meta$characteristics..CoVID.19.severity,sc.meta$characteristics..Outcome,sep="_")
sc.meta$Sev_Outcome[grep("control.*",sc.meta$Sev_Outcome)]<-"control"
sc.meta$Sev_Outcome[grep("mild*",sc.meta$Sev_Outcome)]<-"mild/moderate"
sc.meta$Sev_Outcome<-factor(sc.meta$Sev_Outcome,levels = c("control","mild/moderate","severe/critical_discharged","severe/critical_deceased"))

sc.meta.cell$Sev_Outcome<-paste(sc.meta.cell$CoVID.19.severity,sc.meta.cell$Outcome,sep="_")
sc.meta.cell$Sev_Outcome[grep("control.*",sc.meta.cell$Sev_Outcome)]<-"control"
sc.meta.cell$Sev_Outcome[grep("mild*",sc.meta.cell$Sev_Outcome)]<-"mild/moderate"
sc.meta.cell$Sev_Outcome<-factor(sc.meta.cell$Sev_Outcome,levels = c("control","mild/moderate","severe/critical_discharged","severe/critical_deceased"))

sc.data$Sev_Outcome<-sc.meta$Sev_Outcome[match(sc.data$sampleID,sc.meta$Sample.name)]
sc.data.nkt$Sev_Outcome<-sc.meta$Sev_Outcome[match(sc.data.nkt$sampleID,sc.meta$Sample.name)]

#==============================================================show.cell funciton====================
show.cell<-function(data,x,y,...){
  eval(parse(text=sprintf("ggplot(data,...)+
                          geom_boxplot(aes(x=%s,y=%s))+
                          geom_jitter(aes(x=%s,y=%s))",x,y,x,y)
             )
       )
}

show.cell.col<-function(data,x,y,comparison=NULL,comparison1=NULL,comparison2=NULL,signif_level=T,
                        show.point=F,show.pvalue=F,show.allpvalue=F,decreasing=NULL,
                        signif_ensemble=F,pvalue_cut=c(0.05,0.01,1e-3,1e-4,1e-10,1e-20,1e-50,1e-90),
                        test.args=NULL,vjust=0.6,height.step=0.5,text_size=1.5,tip_length=0.2,width=0.1){
  eval(parse(text=sprintf("data.a<-data%%>%%
    dplyr::group_by(%s)%%>%%
    dplyr::summarise(mean=mean(%s,na.rm=T),
              sd=sd(%s,na.rm = T),
              N=length(%s),
              se=sd/sqrt(N))",
                          x,y,y,y)))
  data.a<-as.data.frame(data.a)
  if(!is.null(decreasing))
  {
    data.a<-data.a[order(data.a$mean,decreasing = decreasing),]
    data.a[,x]<-factor(data.a[,x],levels = data.a[,x])
  }
  p<-ggplot(data.a)+
    geom_col(aes(x=eval(parse(text=x)),
                 y=mean,
                 fill=eval(parse(text=x))),
             width = 0.5)
  if(show.point)
  {
    p<-p+geom_jitter(data=data,
                    mapping=aes(x=eval(parse(text=x)),
                                y=eval(parse(text=y))
                                )
                    )
  }
  p<-p+
    geom_errorbar(aes(x=eval(parse(text=x)),
                      ymin=mean-se,
                      ymax=mean+se,
                      group=eval(parse(text=x)),
                      width=0.2),
                  position = position_dodge(width=0.8),
                  size=0.2)
  if(class(data[,x])=="factor")
  {
    comp<-levels(data[,x])
    comp<-comp[comp%in%data[,x]]
  } else{
    if(is.null(decreasing))
    {
      comp<-levels(factor(data[,x]))
    } else{
      comp<-data.a[,x]
    }
  }
  data.a<-na.omit(data.a)
  height<-max(data.a$mean+data.a$se)
  comp1<-NULL
  comp2<-NULL
  if(!is.null(comparison1)&!is.null(comparison2))
  {
    comp1<-comparison1
    comp2<-comparison2
  }
  if(!is.null(comparison))
  {
    a<-tstrsplit(comparison,split="-")
    comp1<-a[[1]]
    comp2<-a[[2]]
  }
  if(is.null(comp1))
  {
    comp1<-comp
    comp2<-comp
  }
  

  if(signif_level)
  {
    p.annot<-c()
    if(!is.null(comparison))
    {
      for(i in 1:length(comparison))
      {
        if(!paste(comp1[i],comp2[i],sep="-")%in%names(p.annot)&!paste(comp2[i],comp1[i],sep="-")%in%names(p.annot))
        {
          p.annot[paste(comp1[i],comp2[i],sep = "-")]<-do.call(wilcox.test,c(list(formula=eval(parse(text=sprintf("%s~%s",y,x))),
                                                                                  data=data[data[,x]%in%c(comp1[i],comp2[i]),]),
                                                                             test.args))$p.value
        }
      }
    } else{
      for(j in comp1)
      {
        for(k in comp2)
        {
          if(j!=k)
          {
            if(!paste(j,k,sep="-")%in%names(p.annot)&!paste(k,j,sep="-")%in%names(p.annot))
            {
              p.annot[paste(j,k,sep = "-")]<-do.call(wilcox.test,c(list(formula=eval(parse(text=sprintf("%s~%s",y,x))),
                                                                        data=data[data[,x]%in%c(j,k),]),
                                                                   test.args))$p.value
            }
          }
        }
      }
    }
    p.annot.text<-sapply(p.annot,function(p){
      if(show.allpvalue) return(sprintf("%.3g",p))
      if(p>0.05) return(NA)
      if(p<=0.05&p>0.01) return(sprintf("*"))
      if(p<=0.01&p>0.001) return(sprintf("**"))
      if(p<=0.001&p>0.0001) return(sprintf("***"))
      if(show.pvalue) return(sprintf("%.3g",p))
      return("****")
    })
    if(signif_ensemble&!all(is.na(p.annot.text)))
    {
      data.signif<-data.frame(From=tstrsplit(names(p.annot),split="-")[[1]],
                              To=tstrsplit(names(p.annot),split="-")[[2]],
                              signif=p.annot)
      data.signif$text<-p.annot.text
      data.signif$From.x<-as.numeric(factor(data.signif$From,levels = factor(comp)))
      data.signif$To.x<-as.numeric(factor(data.signif$To,levels = factor(comp)))
      data.signif<-data.signif[order(data.signif$From.x,data.signif$To.x),]
      if(is.null(pvalue_cut)) pvalue_cut<-c(0.05,0.01,1e-3,1e-4)
      data.signif$pvalue_cut<-sapply(data.signif$signif,function(p){g<-0;for(i in 1:length(pvalue_cut)){if(p<pvalue_cut[i]){g<-g+1}};return(g)})
      data.signif$group<-1
      for(i in 1:(nrow(data.signif)-1))
      {
        if(data.signif$From[i]==data.signif$From[i+1]&data.signif$To.x[i]+1==data.signif$To.x[i+1])
        {
          if(data.signif$pvalue_cut[i]==data.signif$pvalue_cut[i+1])
          {
            data.signif$group[i+1]<-data.signif$group[i]
          } else{
            data.signif$group[i+1]<-data.signif$group[i]+1
          }
        } else{
          data.signif$group[i+1]<-data.signif$group[i]+1
        }
      }
      data.signif$length<-sapply(data.signif$group,function(p){return(length(which(data.signif$group==p)))})
      data.signif$To.x.mod<-sapply(data.signif$group,function(p){return(mean(data.signif$To.x[data.signif$group==p]))})
      data.signif<-data.signif[!duplicated(data.signif$group),]
      data.signif<-na.omit(data.signif)
      y_position = height+(1:nrow(data.signif))*height.step*height/nrow(data.signif)
      data.signif$y_postion<-y_position
      data.signif$length<-data.signif$length-1
      data.signif$is.length<-as.numeric(data.signif$length>0)
      tip_length2<-tip_length*height/5
      vjust2<-vjust*height/5
      text_size2<-text_size*height/5
      if(nrow(data.signif)>0)
      {
        p<-p+geom_segment(data=data.signif,mapping = aes(x=From,xend=To.x.mod,
                                                         y=y_postion,yend=y_postion),size=width)+
          geom_segment(data=data.signif,mapping = aes(x=From,xend=From,y=y_postion,yend=y_position-tip_length2),size=width)+
          geom_segment(data=data.signif,mapping = aes(x=To.x.mod,xend=To.x.mod,y=y_postion,yend=y_position-tip_length2),size=width)+
          geom_segment(data=data.signif,mapping = aes(x=To.x.mod-length/2,xend=To.x.mod+length/2,
                                                      y=y_position-tip_length2,yend=y_position-tip_length2),size=width)+
          geom_segment(data=data.signif,mapping = aes(x=To.x.mod-length/2,xend=To.x.mod-length/2,
                                                      y=y_position-2*tip_length2*is.length,yend=y_position-tip_length2*is.length),size=width)+
          geom_segment(data=data.signif,mapping = aes(x=To.x.mod+length/2,xend=To.x.mod+length/2,
                                                      y=y_position-2*tip_length2*is.length,yend=y_position-tip_length2*is.length),size=width)+
          annotate("text",x=(data.signif$From.x+data.signif$To.x.mod)/2,
                   y=data.signif$y_postion+vjust*height/10,
                   label = paste("<",pvalue_cut[data.signif$pvalue_cut],sep=""),
                   size=text_size)
      }
    } else{
      n<-length(comp)
      n<-n*(n-1)/2
      xycomp<-strsplit(names(p.annot),split="-")
      suppressWarnings(p<-p+geom_signif(data=data,
                                        aes(x=eval(parse(text=x)),
                                            y=eval(parse(text=y))
                                        ),
                                        comparisons = xycomp,
                                        # y_position = height+(1:n)*0.5*height/n,
                                        test = "wilcox.test",
                                        test.args = test.args,
                                        tip_length = 0,
                                        map_signif_level = function(p){
                                          if(show.allpvalue) return(sprintf("%.3g",p))
                                          if(p>0.05) return(NA)
                                          if(p<=0.05&p>0.01) return(sprintf("*"))
                                          if(p<=0.01&p>0.001) return(sprintf("**"))
                                          if(p<=0.001&p>0.0001) return(sprintf("***"))
                                          if(show.pvalue) return(sprintf("%.3g",p))
                                          return("****")
                                        },
                                        y_position = {which=which(!is.na(p.annot.text));a<-1:n*0;a[which]<-(1:length(which))*height.step*height/length(which);height+a},
                                        vjust = vjust,
                                        size = 0.2,textsize = text_size
      )
      )
    }
  }
  # rm(data.a,comp,xycomp,height,n,p.annot)
  p
}

show.cell.gene.outcome.col<-function(data,target.gene,cell.subgroup,outcome,
                                     show.pvalue=F,show.allpvalue=F,show.point=F,
                                     signif_ensemble=F,pvalue_cut=c(0.05,0.01,1e-3,1e-4,1e-10,1e-20,1e-50,1e-90),
                                     vjust=0.6,text_size=1.5,height.step=0.05,tip_length=0.2,width=0.1,
                                     inner.p=T,outer.p=T,sp.p=c())
{
  if(class(data)%in%c("matrix","data.frame"))
  {
    data.gene<-data.frame(gene=data[target.gene,],
                          NKT_TIM3=cell.subgroup,
                          Outcome=outcome)
  } else{
    data.gene<-data.frame(gene=data,
                          NKT_TIM3=cell.subgroup,
                          Outcome=outcome)
  }
  
  
  data.gene$NKT_TIM3<-factor(data.gene$NKT_TIM3)
  data.gene$Outcome<-factor(data.gene$Outcome)
  
  data.a<-data.gene%>%
    group_by(NKT_TIM3,Outcome)%>%
    dplyr::summarise(mean=mean(gene,na.rm=T),
              sd=sd(gene,na.rm = T),
              N=length(gene),
              se=sd/sqrt(N))
  
  p.ans<-ggplot()+
    geom_col(data=data.a,
             mapping = aes(x=NKT_TIM3,
                           y=mean,
                           fill=Outcome),
             width = 0.5,
             position = position_dodge(0.5))+
    geom_errorbar(mapping = aes(x=NKT_TIM3,
                                ymin=mean-se,
                                ymax=mean+se,
                                group=Outcome,
                                width=0.2),
                  data = data.a,
                  position = position_dodge(0.5),
                  size=0.2)
  if(show.point)
  {
    p.ans<-p.ans+geom_jitter(data=data.gene,
                             mapping=aes(x=NKT_TIM3,y=gene))
  }
  
  p.annot.NKT<-c()
  
  if(inner.p)
  {
    for(i in levels(data.gene$NKT_TIM3))
    {
      for(j in levels(data.gene$Outcome))
      {
        for(k in levels(data.gene$Outcome))
        {
          if(j!=k)
          {
            if(!paste(i,j,k,sep="-")%in%names(p.annot.NKT)&!paste(i,k,j,sep="-")%in%names(p.annot.NKT))
            {
              p.annot.NKT[paste(i,j,k,sep = "-")]<-wilcox.test(gene~Outcome,
                                                               data.gene[data.gene$NKT_TIM3==i&data.gene$Outcome%in%c(j,k),])$p.value
            }
          }
        }
      }
    }
  }
  
  p.annot.Outcome<-c()
  
  if(outer.p)
  {
    for(i in levels(data.gene$Outcome))
    {
      for(j in levels(data.gene$NKT_TIM3))
      {
        for(k in levels(data.gene$NKT_TIM3))
        {
          if(j!=k)
          {
            if(!paste(i,j,k,sep="-")%in%names(p.annot.Outcome)&!paste(i,k,j,sep="-")%in%names(p.annot.Outcome))
            {
              p.annot.Outcome[paste(i,j,k,sep = "-")]<-wilcox.test(gene~NKT_TIM3,
                                                                   data.gene[data.gene$Outcome==i&data.gene$NKT_TIM3%in%c(j,k),])$p.value
            }
          }
        }
      }
    }
  }
  p.annot.NKT.text<-sapply(p.annot.NKT,function(p){
    if(show.allpvalue) return(sprintf("%.3g",p))
    if(is.nan(p)) return(NA)
    if(is.na(p)) return(NA)
    if(is.null(p)) return(NA)
    if(p>0.05) return(NA)
    if(p<=0.05&p>0.01) return(sprintf("*"))
    if(p<=0.01&p>0.001) return(sprintf("**"))
    if(p<=0.001&p>0.0001) return(sprintf("***"))
    if(show.pvalue) return(sprintf("%.3g",p))
    return(sprintf("****"))
  })
  if(length(p.annot.NKT.text)==0) p.annot.NKT.text<-c()
  p.annot.Outcome.text<-sapply(p.annot.Outcome,function(p){
    if(show.allpvalue) return(sprintf("%.3g",p))
    if(is.nan(p)) return(NA)
    if(is.na(p)) return(NA)
    if(is.null(p)) return(NA)
    if(p>0.05) return(NA)
    if(p<=0.05&p>0.01) return(sprintf("*"))
    if(p<=0.01&p>0.001) return(sprintf("**"))
    if(p<=0.001&p>0.0001) return(sprintf("***"))
    if(show.pvalue) return(sprintf("%.3g",p))
    return(sprintf("****"))
  })
  if(length(p.annot.Outcome.text)==0) p.annot.Outcome.text<-c()
  height=max(data.a$mean+data.a$se)
  n.in<-length(levels(data.gene$Outcome))
  n.out<-length(levels(data.gene$NKT_TIM3))
  if(signif_ensemble&!all(is.na(c(p.annot.NKT.text,p.annot.Outcome.text))))
  {
    if(!is.null(p.annot.NKT))
    {
      data.signif<-data.frame(From1=tstrsplit(names(p.annot.NKT),split="-")[[1]],
                              From2=tstrsplit(names(p.annot.NKT),split="-")[[2]],
                              To1=tstrsplit(names(p.annot.NKT),split="-")[[1]],
                              To2=tstrsplit(names(p.annot.NKT),split="-")[[3]],
                              signif=p.annot.NKT)
    } else{
      data.signif<-data.frame(From1=NA,
                              From2=NA,
                              To1=NA,
                              To2=NA,
                              signif=NA)
      data.signif<-data.signif[-1,]
    }
    
    data.signif<-rbind(data.signif,data.frame(From1=tstrsplit(names(p.annot.Outcome),split="-")[[2]],
                                              From2=tstrsplit(names(p.annot.Outcome),split="-")[[1]],
                                              To1=tstrsplit(names(p.annot.Outcome),split="-")[[3]],
                                              To2=tstrsplit(names(p.annot.Outcome),split="-")[[1]],
                                              signif=p.annot.Outcome))
    data.signif$text<-c(p.annot.NKT.text,p.annot.Outcome.text)
    data.signif$From1.x<-as.numeric(factor(data.signif$From1,levels = levels(data.gene$NKT_TIM3)))
    data.signif$From2.x<-as.numeric(factor(data.signif$From2,levels = levels(data.gene$Outcome)))
    data.signif$To1.x<-as.numeric(factor(data.signif$To1,levels = levels(data.gene$NKT_TIM3)))
    data.signif$To2.x<-as.numeric(factor(data.signif$To2,levels = levels(data.gene$Outcome)))
    if(is.null(pvalue_cut)) pvalue_cut<-c(0.05,0.01,1e-3,1e-4)
    data.signif$pvalue_cut<-sapply(data.signif$signif,function(p){g<-0;for(i in 1:length(pvalue_cut)){if(p<pvalue_cut[i]){g<-g+1}};return(g)})
    data.signif$From<-data.signif$From1.x+((data.signif$From2.x-mean(unique(data.signif$From2.x)))/max(data.signif$From2.x))/2
    data.signif$To<-data.signif$To1.x+((data.signif$To2.x-mean(unique(data.signif$To2.x)))/max(data.signif$To2.x))/2
    data.signif$group<-1
    for(i in 1:(nrow(data.signif)-1))
    {
      if(all(data.signif[i,c("From1.x","From2.x")]==data.signif[i+1,c("From1.x","From2.x")])&data.signif$To1.x[i]==data.signif$To1.x[i]&data.signif$To2.x[i]+1==data.signif$To2.x[i+1])
      {
        if(data.signif$pvalue_cut[i]==data.signif$pvalue_cut[i+1])
        {
          data.signif$group[i+1]<-data.signif$group[i]
        } else{
          data.signif$group[i+1]<-data.signif$group[i]+1
        }
      } else{
        data.signif$group[i+1]<-data.signif$group[i]+1
      }
    }
    data.signif$length<-sapply(data.signif$group,function(p){return(length(which(data.signif$group==p)))})
    data.signif$To.x.mod<-sapply(data.signif$group,function(p){return(mean(data.signif$To[data.signif$group==p]))})
    y_position<-height*(1+c(rep(1:(n.in*(n.in-1)/2),time=n.out))*height.step)
    if(!inner.p) y_position<-c()
    y_position2<-height*(1+c(1:(n.in*n.out*(n.out-1)/2)+n.in*(n.in-1)/2)*height.step)
    if(!outer.p) y_position2<-c()
    data.signif$y<-c(y_position,y_position2)
    data.signif<-data.signif[!duplicated(data.signif$group),]
    data.signif<-na.omit(data.signif)
    data.signif<-data.signif[data.signif$signif<max(pvalue_cut),]
    data.signif$length<-data.signif$length-1
    data.signif$is.length<-as.numeric(data.signif$length>0)
    data.signif$y<-data.signif$y-min(data.signif$y-height*(1+height.step))
    if(!is.null(sp.p)) data.signif<-data.signif[data.signif$From1.x==sp.p|data.signif$To1.x==sp.p,]
    tip_length2<-tip_length*height/5
    vjust2<-vjust*height/20
    text_size2<-text_size*height/5
    if(nrow(data.signif)>0)
    {
      p.ans<-p.ans+geom_segment(data=data.signif,mapping = aes(x=From,xend=To.x.mod,
                                                               y=y,yend=y),size=width)+
        geom_segment(data=data.signif,mapping = aes(x=From,xend=From,y=y,yend=y-tip_length2),size=width)+
        geom_segment(data=data.signif,mapping = aes(x=To.x.mod,xend=To.x.mod,y=y,yend=y-tip_length2),size=width)+
        geom_segment(data=data.signif,mapping = aes(x=To.x.mod-length/2,xend=To.x.mod+length/2,
                                                    y=y-tip_length2,yend=y-tip_length2),size=width)+
        geom_segment(data=data.signif,mapping = aes(x=To.x.mod-length/2,xend=To.x.mod-length/2,
                                                    y=y-2*tip_length2*is.length,yend=y-tip_length2*is.length),size=width)+
        geom_segment(data=data.signif,mapping = aes(x=To.x.mod+length/2,xend=To.x.mod+length/2,
                                                    y=y-2*tip_length2*is.length,yend=y-tip_length2*is.length),size=width)+
        annotate("text",x=(data.signif$From+data.signif$To.x.mod)/2,
                 y=data.signif$y+vjust2,
                 label = paste("<",pvalue_cut[data.signif$pvalue_cut],sep=""),
                 size=text_size)
      p.ans<-p.ans+ylab(target.gene)
    }
  } else{
    p.annot<-c(p.annot.NKT.text,p.annot.Outcome.text)
    p.ans.xmin<-c()
    p.ans.xmax<-c()
    y_position<-c()
    if(inner.p)
    {
      p.ans.xmin<-rep(as.numeric(tstrsplit(group.combine(1:n.in),"-")[[1]])/n.in*0.5,time=n.out)+rep(1:n.out,each=n.in*(n.in-1)/2)
      p.ans.xmax<-rep(as.numeric(tstrsplit(group.combine(1:n.in),"-")[[2]])/n.in*0.5,time=n.out)+rep(1:n.out,each=n.in*(n.in-1)/2)
      y_position<-height*(1+c(rep(1:(n.in*(n.in-1)/2),time=n.out))*height.step)
    }
    if(outer.p)
    {
      p.ans.xmin2<-rep(as.numeric(tstrsplit(group.combine(1:n.out),"-")[[1]]),each=n.in)+rep(1:n.in/n.in,time=n.out*(n.out-1)/2)*0.5
      p.ans.xmax2<-rep(as.numeric(tstrsplit(group.combine(1:n.out),"-")[[2]]),each=n.in)+rep(1:n.in/n.in,time=n.out*(n.out-1)/2)*0.5
      y_position2<-height*(1+c(1:(n.in*n.out*(n.out-1)/2)+n.in*(n.in-1)/2)*height.step)
      p.ans.xmin<-c(p.ans.xmin,p.ans.xmin2)
      p.ans.xmax<-c(p.ans.xmax,p.ans.xmax2)
      y_position<-c(y_position,y_position2)
    }
    p.ans.x.just<--0.25*(1+1/n.in)
    
    p.ans.xy<-data.frame(xmin=p.ans.xmin,
                         xmax=p.ans.xmax,
                         y_position=y_position,
                         annot=p.annot)
    p.ans.xy$xmin.f<-floor(p.ans.xy$xmin)
    p.ans.xy$xmax.f<-floor(p.ans.xy$xmax)
    if(!is.null(sp.p))
    {
      p.ans.xy<-p.ans.xy[p.ans.xy$xmin.f%in%sp.p|p.ans.xy$xmax.f%in%sp.p,]
    }
    p.ans.xy$y_position<-p.ans.xy$y_position-(min(p.ans.xy$y_position)-height*(1+height.step))
    
    p.ans<-p.ans+geom_signif(
      data=data.a,
      mapping = aes(x=NKT_TIM3,
                    y=mean),
      y_position = p.ans.xy$y_position,
      xmin = p.ans.xy$xmin+p.ans.x.just,
      xmax = p.ans.xy$xmax+p.ans.x.just,
      tip_length = 0,
      annotations = p.ans.xy$annot,
      size = 0.2,
      vjust=vjust,
      textsize = text_size)
    
    p.ans<-p.ans+
      ylab(target.gene)
  }
  p.ans
  # rm(data.a,data.gene,p.ans,p.annot,p.annot.NKT,p.annot.Outcome,
  # p.ans.x.just,p.ans.xmax,p.ans.xmax2,p.ans.xmin,p.ans.xmin2,y_position,y_position2,p.ans.xy,n.in,n.out)
}

show.cell.gene.col<-function(data,cell.subgroup,target.gene)
{
  
  for(i in 1:length(target.gene))
  {
    eval(parse(text=sprintf("p.ans%d<-show.cell.col(data.frame(gene=data[target.gene[%d],],
                                NKT_subgroup=cell.subgroup),
                     x=\"NKT_subgroup\",
                     y=\"gene\",
                     show.point = F)",
                            i,i)
    )
    )
  }
  
  
  p.ans<-ggplot(cbind(eval(parse(text=sprintf("rbind(%s)",paste("p.ans",1:length(target.gene),"$data",sep = "",collapse = ",")))),data.frame(Gene=rep(target.gene,each=2))),
                aes(x=Gene,y=mean,fill=NKT_subgroup))+
    geom_col(position = "dodge")+
    geom_errorbar(aes(x=Gene,
                      ymin=mean-se,
                      ymax=mean+se,
                      width=0.2),
                  position = position_dodge(0.8))
  
  p.ans.annot<-eval(parse(text=sprintf("c(%s)",paste("wilcox.test(gene~NKT_subgroup,p.ans",
                                                     1:length(target.gene),
                                                     "$layers[[3]]$data)$p.value",
                                                     sep="",
                                                     collapse = ","
  )
  )))
  
  p.ans.x<-1:length(target.gene)-0.2
  
  p.ans<-p.ans+geom_signif(y_position = max(p.ans$data$mean+p.ans$data$se)*1.05,
                           xmin = p.ans.x,
                           xmax = p.ans.x+0.4,
                           annotations = signif(p.ans.annot,3))
  
  p.ans<-p.ans+
    scale_fill_manual(values = RColorBrewer::brewer.pal(3,"Set1")[2:1])+
    xlab("Target Gene")+
    ylab("Target Gene Expression")+
    theme(panel.background = element_blank(),
          axis.line = element_line(),
          plot.title = element_text(hjust = 0.5,vjust = 0.5))
  
  eval(parse(text=sprintf("rm(p.ans.x,p.ans.annot,%s)",
                          paste("p.ans",1:length(target.gene),sep = "",collapse = ",")
  )
  ))
  p.ans
}

show.cell.gene.cor<-function(data,gene,target.gene,...)
{
  for(i in 1:length(target.gene))
  {
    eval(parse(text=sprintf("p.ans%d<-show.gene.cor(data,\"%s\",\"%s\",...)",
                            i,gene,target.gene[i]
    )
    ))
    eval(parse(text=sprintf("p.ans%d<-p.ans%d+xlab(label=\"TIM3\")+
                            labs(title=\"TIM3-%s Correlation\")",
                            i,i,target.gene[i]
    )
    ))
  }
  p.ans<-eval(parse(text=sprintf("list(%s)",
                                 paste("p.ans",1:length(target.gene),sep = "",collapse = ",")
  )))
  p.ans
}



#==============================================================Figure 1=========================
p.1a1<-DimPlot(sc.data,reduction = "umap",label = T,group.by = "cluster_tree_annot_merge",
              raster = T,pt.size = 0.01,repel = T,
              label.size = 2.5)+
  labs(title = "")+
  theme(panel.background = element_blank(),
        axis.line = element_line(),
        plot.title = element_text(hjust = 0.5,vjust = 0.5,size=10,face = "bold"),
        legend.position = "top",
        legend.justification = "center",
        legend.text = element_text(size=6),
        legend.key.size = unit(0.1,"lines"),
        legend.spacing = unit(0,"cm"),
        legend.margin = margin(0,0,0,0,"pt"),
        legend.box.margin = margin(0,0,0,0,"pt"),
        axis.title = element_text(size=8),
        axis.text = element_text(size=7))
p.1a1$layers[[1]]$geom$default_aes$alpha<-0.3
p.1a1$guides$colour$override.aes<-list(alpha=1)

p.1a2<-DimPlot(sc.data,reduction = "umap",label = T,group.by = "majortype",
        raster = T,pt.size = 0.01,repel = T,
        label.size = 2.5)+
  guides(color=guide_legend(nrow=1,byrow = T))+
  labs(title = "")+
  theme(panel.background = element_blank(),
        axis.line = element_line(),
        plot.title = element_text(hjust = 0.5,vjust = 0.5,size=10,face = "bold"),
        plot.margin = unit(c(0,0,0,0),"cm"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "top",
        legend.justification = "center",
        legend.text = element_text(size=6),
        legend.key.size = unit(0.1,"lines"),
        legend.spacing = unit(0,"cm"),
        legend.margin = margin(0,0,0,0,"pt"),
        legend.box.margin = margin(0,0,0,0,"pt"),
        axis.title = element_text(size=8),
        axis.text = element_text(size=7))
p.1a2$layers[[1]]$geom$default_aes$alpha<-0.3
p.1a2$guides$colour$override.aes<-list(alpha=1)


p.1a3<-DimPlot(sc.data,reduction = "umap",label = F,group.by = "Sev_Outcome",
               raster = T,pt.size = 0.01)+
  scale_color_manual(values = c("#74C476","#FEB24C","#E31A1C","#800026"),
                     labels=function(x){gsub("_","\n",x)})+
  labs(title = "")+
  theme(panel.background = element_blank(),
        axis.line = element_line(),
        plot.title = element_text(hjust = 0.5,vjust = 0.5,size=10,face = "bold"),
        plot.margin = unit(c(0,0,0,0),"cm"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "top",
        legend.justification = "center",
        legend.text = element_text(size=6),
        legend.key.size = unit(0.1,"lines"),
        legend.spacing = unit(0,"cm"),
        legend.margin = margin(0,0,0,0,"pt"),
        legend.box.margin = margin(0,0,0,0,"pt"),
        axis.title = element_text(size=8),
        axis.text = element_text(size=7))
p.1a3$layers[[1]]$geom$default_aes$alpha<-0.3
p.1a3$guides$colour$override.aes<-list(alpha=1)

sc.data$sampleID_plot<-factor(sc.data$sampleID)
sc.data$sampleID_plot<-as.numeric(sc.data$sampleID_plot)
# sc.data$sampleID_plot<-paste("S",sprintf("%03d",sc.data$sampleID_plot),sep = "")

meltP<-data.frame(UMAP1=sc.data@reductions$umap@cell.embeddings[,1],
                  UMAP2=sc.data@reductions$umap@cell.embeddings[,2],
                  sample=sc.data$sampleID_plot)

nSample=length(unique(sc.data$sampleID_plot))
p.1a4<-ggplot(meltP,aes(x=UMAP1,y=UMAP2,color=sample))+
  scattermore::geom_scattermore(alpha=0.3)+
  scale_color_stepsn(guide=guide_colorbar(title = "Sample",title.position = "left",title.hjust = 0.5,title.vjust = 1,
                                          nbin=nSample,ticks.colour = "black",ticks.linewidth = 0.2,
                                          barwidth = unit(0.4,"npc"),barheight = unit(0.8/nSample,"npc"),label = F),
                     colors = colorRampPalette(brewer.pal(9, "Set1"))(nSample)[sample(1:nSample)],
                     n.breaks=nSample,
                     show.limits=T)+
  theme(panel.background = element_blank(),
        axis.line = element_line(),
        plot.title = element_text(hjust = 0.5,vjust = 0.5,size=10,face = "bold"),
        plot.margin = unit(c(0,0,0,0),"cm"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        legend.direction = "horizontal",
        legend.position = c(0.5,1.1),
        legend.justification = "center",
        legend.title = element_text(hjust = 0.5,vjust = 0.5,size = 7),
        legend.text = element_text(size=6),
        legend.key.size = unit(0.01,"lines"),
        legend.box.spacing = unit(0.01,"lines"),
        legend.spacing = unit(0,"cm"),
        legend.margin = margin(0,0,0,0,"pt"),
        legend.box.margin = margin(0,0,0,0,"pt"),
        axis.title = element_text(size=8),
        axis.text = element_text(size=7))


sc.data$annot_plot<-sc.data$annot
sc.data$annot_plot[!sc.data$annot_plot%in%c("NKT","NK_cell","T_cell:MAI")]<-NA

p.1a5<-DimPlot(sc.data,reduction = "umap",label = T,group.by = "annot_plot",
               raster = T,pt.size = 0.01,repel = T,
               label.size = 2.5)+
  guides(color=guide_legend(nrow=1,byrow = T))+
  labs(title = "")+
  theme(panel.background = element_blank(),
        axis.line = element_line(),
        plot.title = element_text(hjust = 0.5,vjust = 0.5,size=10,face = "bold"),
        plot.margin = unit(c(0,0,0,0),"cm"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "top",
        legend.justification = "center",
        legend.text = element_text(size=6),
        legend.key.size = unit(0.1,"lines"),
        legend.spacing = unit(0,"cm"),
        legend.margin = margin(0,0,0,0,"pt"),
        legend.box.margin = margin(0,0,0,0,"pt"),
        axis.title = element_text(size=8),
        axis.text = element_text(size=7))
p.1a5$layers[[1]]$geom$default_aes$alpha<-0.3
p.1a5$guides$colour$override.aes<-list(alpha=1)

sc_edge<-as.data.frame(sc.data@tools$BuildClusterTree$edge)
sc_node<-data.frame(node=unique(c(sc_edge[,2],sc_edge[,1])))
node_name<-sc.data@tools$BuildClusterTree$tip.label
sc_node$name<-"NA"
sc_node$name[match(1:length(node_name),sc_node$node)]<-node_name
sc_node$cell_type<-gsub("^[0-9]+-(.*)","\\1",sc_node$name)
sc_node$id=NA
myleaves=which(is.na(match(sc_node$node,sc_edge$V1)))
sc_node$id[myleaves]<-1:length(myleaves)
sc_node$angle<-90-360*sc_node$id/length(myleaves)
sc_node$hjust<-ifelse(sc_node$angle< -90,0,1)
sc_node$angle<-ifelse(sc_node$angle< -90,sc_node$angle+180,sc_node$angle)
sc_node$size<-as.data.frame(table(sc.data$cluster_tree_annot)[match(sc_node$name,levels(sc.data$cluster_tree_annot))])$Freq

scgraph<-graph_from_data_frame(sc_edge,vertices = sc_node)

p.1a6<-ggraph(scgraph,layout = "dendrogram",circular=T)+
  geom_edge_elbow(color="grey",lineend = "square",linejoin = "mitre")+
  geom_node_text(aes(x=x*1.15,y=y*1.15,label=name,filter=leaf,color=cell_type,
                     angle=angle,hjust=hjust),size=2)+
  geom_node_point(aes(x=x*1.07,y=y*1.07,filter=leaf,color=cell_type,size=size),alpha=0.6)+
  # scale_colour_manual(values= rep( brewer.pal(9,"Paired") , 30)) +
  scale_size_continuous( range = c(0.5,4) ) +
  theme_void() +
  theme(
    legend.position="none"
    # plot.margin=unit(c(1,1,1,1),"cm"),
  ) +
  expand_limits(x = c(-2.3, 2.3), y = c(-2.3, 2.3))


aligned1<-cowplot::align_plots(p.1a3,
                               p.1a4,
                               align="hv",axis="tblr")
aligned2<-cowplot::align_plots(p.1a2,
                               p.1a5,
                               align="hv",axis="tblr")


p.1a<-ggarrange(ggarrange(p.1a1+theme(plot.margin = unit(c(0,0,0,0),"cm")),
                          ggarrange(cowplot::ggdraw(aligned1[[1]])+theme(plot.margin = unit(c(0,0,0,0),"cm")),
                                          cowplot::ggdraw(aligned1[[2]])+theme(plot.margin = unit(c(0,0,0,0),"cm")),
                                          ncol = 2,align = "hv"),
                          nrow = 2,heights = c(2,1)),
                ggarrange(cowplot::ggdraw(aligned2[[1]])+theme(plot.margin = unit(c(0,0,0,0),"cm")),
                          cowplot::ggdraw(aligned2[[2]])+theme(plot.margin = unit(c(0,0,0,0),"cm")),
                          p.1a6+theme(plot.margin = unit(c(0,0,0,0),"cm")),
                          nrow = 3),
                ncol = 2,widths = c(2,1))

nkr_gene<-c("FCGR3A","NCR3","NCR1","CD244","KLRK1","IL2RB","NCAM1","CD160")
p.1b.list<-list()
for(i in 1:length(nkr_gene))
{
  p.1b.list[[i]]<-show.cell.col(data.frame(cell_cluster=sc.data$cluster_tree_annot_merge,
                                            gene=sc.data@assays$RNA@data[nkr_gene[i],]),
                                 "cell_cluster",
                                 "gene",
                                 signif_level = F,
                                 decreasing = NULL
  )+
    scale_y_continuous(expand=expansion(c(0,0.05)))+
    labs(title = "")+
    xlab("")+
    ylab(nkr_gene[i])+
    theme(panel.background = element_blank(),
          legend.position = "none",
          legend.box.just = "center",
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y.left = element_text(angle=0),
          axis.ticks = element_blank(),
          plot.title = element_text(hjust = 0.5,vjust = 0.5,size=10,face = "bold"),
          plot.margin = unit(c(0, 0, 0, 0), "cm"))
}
for(i in 1:length(p.1b.list))
{
  p.1b.list[[i]]<-p.1b.list[[i]]+theme(axis.line = element_line(size=0.2))
}
p.1b.list[[length(p.1b.list)]]<-p.1b.list[[length(p.1b.list)]]+theme(axis.text.x = element_text(angle = 45,hjust=1,vjust=1))

p.1b<-patchwork::wrap_plots(plotlist = p.1b.list, ncol = 1,heights = 1)+
  theme(panel.spacing = unit(0,"lines"))

p.1c.list<-list()
p.1c.gene<-c("CD3D","CD4","CD8A","TRGV9","TRDV2",
             "KLRF1","MS4A1","MZB1",
             "CD14","FCGR3A","LYZ","PPBP")

DefaultAssay(sc.data)<-"RNA"
p.1c.list<-FeaturePlot(sc.data,p.1c.gene,combine = F,slot = "data")
DefaultAssay(sc.data)<-"integrated"
for(i in 1:length(p.1c.gene))
{
  p.1c.list[[i]]<-p.1c.list[[i]]+scale_color_gradientn(colors=c("grey","#0000CC","#0000FF"))
    theme_blank()+NoLegend()+theme(plot.title = element_text(size=6,hjust = 0.5),
                                                                plot.margin = unit(c(0,0,0,0),"cm"))
}
p.1c.list[[2]]<-p.1c.list[[2]]+scale_color_gradientn(colors=c("grey","blue","blue","darkblue"))


p.1c<-do.call(ggarrange,c(p.1c.list,
                           list(
                             nrow=2,ncol=6,
                             hjust=0,vjust=0,
                             align="hv")))


p.1d<-show.cell.col(sc.meta.cell[sc.meta.cell$Sev_Outcome!="control",],
                    "Sev_Outcome",
                    "Lym/Leu",
                    test.args = list(alternative="greater"),
                    signif_ensemble = T,text_size = 2,vjust = 0.3,tip_length = 0.15,width = 0.1,height.step = 0.2)+
  scale_fill_manual(values = c("#FEB24C","#E31A1C","#800026"),
                    labels=function(x){gsub("_","\n",x)},
                    guide=guide_legend(title="COVID-19 Severity"))+
  scale_y_continuous(expand=expansion(c(0,0.05)))+
  xlab("")+
  ylab("Lymphocyte Percentage")+
  theme(panel.background = element_blank(),
        axis.line = element_line(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5,vjust = 0.5),
        legend.position = "top")

p.1e<-show.cell.col(sc.meta.cell,
                    "Sev_Outcome",show.allpvalue = F,
                    "(T_cell_NKT+T_cell_gdT+T_cell_MAIT)/cell_count",
                    signif_ensemble = T,text_size = 2,vjust = 0.3,tip_length = 0.15,width = 0.1,height.step = 0.2)+
  scale_fill_manual(values = c("#74C476","#FEB24C","#E31A1C","#800026"),
                    labels=function(x){gsub("_","\n",x)},
                    guide=guide_legend(title="COVID-19 Severity"))+
  scale_y_continuous(expand=expansion(c(0,0.05)))+
  xlab("")+
  ylab("NKT cell Percentage")+
  theme(panel.background = element_blank(),
        axis.line = element_line(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5,vjust = 0.5,size=10,face = "bold"),
        legend.position = "top")


myget_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

library(patchwork)
p.legend<-myget_legend(p.1e+
                         theme(legend.key.height = unit(c(0.5),"cm"),
                               legend.key.width = unit(c(0.5),"cm"),
                               legend.text = element_text(size=9,vjust=0.5,hjust=0.5),
                               legend.title = element_text(size=9,hjust=0.5,vjust=0.5))
)

p<-ggarrange(ggarrange(ggarrange(p.1a,p.1b,ncol=2,widths = c(3,1),labels = c("Fig.1A","Fig.1B"),hjust = 0,vjust = 0,align = "v",
                                 font.label = list(size=10)),
                       p.1c,
                       nrow = 2,heights = c(5,2),labels = c("","Fig.1C"),hjust = 0,vjust = 0,
                       font.label = list(size=10)),
             ggarrange(p.legend,
                       ggarrange(p.1d,p.1e,ncol = 2,labels = c("Fig.1D","Fig.1E"),hjust = 0,vjust = 0,align = "h",common.legend = T,legend="none",
                                 font.label = list(size=10)),
                       nrow=2,
                       heights = c(2,11)),
             nrow=2,
             heights = c(3,1),
             hjust = 0)+
  theme(plot.margin = unit(c(1,1,1,1),units = "cm"))

pdf("Figure mod1.pdf",width=14,height=20)
p
dev.off()


p.1s.list<-list()

p.theme<-theme(panel.background = element_blank(),
               axis.line = element_line(),
               axis.text.x = element_blank(),
               axis.ticks.x = element_blank(),
               axis.title.y = element_text(size=8),
               plot.title = element_text(hjust = 0.5,vjust = 0.5,size=10,face = "bold"),)

p.1s.y<-c("(T_cell_CD4)/cell_count","(T_cell_CD8)/cell_count","(T_cell_MAIT)/cell_count","(T_cell_gdT)/cell_count",
          "(NK_cell)/cell_count","(B_cell)/cell_count","(Plasma_cell)/cell_count",
          "(`Monocyte_CD14`)/cell_count","(`Monocyte_CD16`)/cell_count","(`mono_DC`)/cell_count","(pDC)/cell_count",
          "(Platelet)/cell_count")
p.1s.ylab<-c("CD4+ T cell Percentage","CD8+ T cell Percentage","MAIT cell Percentage","gdT cell Percentage",
             "NK cell Percentage","B cell Percentage","Plasma cell Percentage",
             "CD14+ Monocyte Percentage","CD16+ Monocyte Percentage","mono-DC Percentage","pDC Percentage",
             "Platelet Percentage")

for(i in 1:length(p.1s.y))
{
  p.1s.list[[i]]<-show.cell.col(sc.meta.cell,
                                "Sev_Outcome",
                                p.1s.y[i],
                                signif_ensemble = T,text_size = 1.5,vjust = 0.5,tip_length = 0.15,width = 0.1)+
    scale_y_continuous(expand=expansion(c(0,0.05)))+
    scale_fill_manual(values = c("#74C476","#FEB24C","#E31A1C","#800026"),
                      labels=function(x){gsub("_","\n",x)},
                      guide=guide_legend(title="COVID-19 Severity"))+
    xlab("")+
    ylab(p.1s.ylab[i])+
    p.theme
}

p.legend<-myget_legend(p.1s.list[[1]]+
                         theme(legend.key.height = unit(c(0.5),"cm"),
                               legend.key.width = unit(c(0.5),"cm"),
                               legend.text = element_text(size=9,vjust=0.5,hjust=0.5),
                               legend.title = element_text(size=9,hjust=0.5,vjust=0.5))
)
p.s<-do.call(ggarrange,c(p.1s.list,
                         list(
                           ncol=3,nrow=4,
                           hjust=0,vjust=0,
                           common.legend=T)))
p.s<-p.s+theme(plot.margin = unit(c(1,1,1,1),units = "cm"))

pdf("Figure S1.pdf",width=8,height=12)
p.s
dev.off()

rm(p.1a1,p.1a2,p.1a3,p.1a4,p.1a5,p.1a6,p.1a,p.1b,p.1b.list,p.1c,p.1c.list,p.1d,p.1e,p,p.theme,p.legend,
   aligned1,aligned2,sc_edge,sc_node,node_name,nSample,meltP,nkr_gene,scgraph,myleaves,
   p.1s.list,p.s,p.1s.ylab,p.1s.gene,p.1s.y,p.1c.gene)




#==============================================================Figure 2=========================



p.2a2<-show.cell.col(data.frame(NKT_subgroup=sc.data.nkt$cluster_annot_nkt,
                               TIM3=sc.data.nkt@assays$RNA@data["HAVCR2",]),
                    "NKT_subgroup",
                    "TIM3",
                    signif_level = T,
                    decreasing = T,
                    signif_ensemble = T,text_size = 2,vjust = 0.5,tip_length = 0.2,width = 0.1,height.step = 1
)+
  scale_y_continuous(expand=expansion(c(0,0.05)))+
  labs(title = "")+
  xlab("")+
  ylab("TIM3")+
  theme(panel.background = element_blank(),
        plot.margin = unit(c(0,0,0,0),"cm"),
        legend.position = "top",
        legend.box.just = "center",
        axis.line = element_line(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5,vjust = 0.5,size=10,face = "bold"))

sc.data.nkt$cluster_annot_nkt_plot<-factor(as.character(sc.data.nkt$cluster_annot_nkt),
                                           levels = p.2b$data$NKT_subgroup)

sc.data.nkt<-sc.data.nkt[,colnames(sc.data.nkt@assays$RNA@counts)!="GTGCATATCTTGTTTG-283"]

DefaultAssay(sc.data.nkt)<-"RNA"

p.2a1<-DimPlot(sc.data.nkt,group.by = "cluster_annot_nkt_plot",
              raster = T,
              pt.size = 0.1,
              repel = T)+
  scale_color_manual(labels=function(x){gsub("TIM3.*","TIM3",x)},
                    guide=guide_legend(title = "NKT subtype"),
                    values = hue_pal()(6))+
  labs(title = "")+
  theme(panel.background = element_blank(),
        plot.margin = unit(c(0,0,0,0),"cm"),
        axis.line = element_line(),
        plot.title = element_text(hjust = 0.5,vjust = 0.5,size=10,face = "bold"),
        legend.position = "top",
        legend.justification = "center",
        legend.text.align = 0,
        legend.text = element_text(size=6),
        legend.key.size = unit(0.1,"lines"),
        legend.spacing = unit(0,"cm"),
        legend.margin = margin(0,0,0,0,"pt"),
        legend.box.margin = margin(0,0,0,0,"pt"),
        axis.title = element_text(size=8),
        axis.text = element_text(size=7))
p.2a1$layers[[1]]$geom$default_aes$alpha<-0.3
p.2a1$guides$colour$override.aes$alpha<-1


sc_edge<-as.data.frame(sc.data.nkt@tools$BuildClusterTree$edge)
sc_node<-data.frame(node=unique(c(sc_edge[,2],sc_edge[,1])))
node_name<-sc.data.nkt@tools$BuildClusterTree$tip.label
sc_node$name<-"NA"
sc_node$name[match(1:length(node_name),sc_node$node)]<-node_name
sc_node$cell_type<-gsub("^[0-9]+-(.*)","\\1",sc_node$name)
sc_node$label<-gsub("^[0-9]+-","",sc_node$name)
sc_node$label<-gsub("TIM3_CD62L","TIM3",sc_node$label)
sc_node$id=NA
myleaves=which(is.na(match(sc_node$node,sc_edge$V1)))
sc_node$id[myleaves]<-1:length(myleaves)
sc_node$angle<-90-360*sc_node$id/length(myleaves)+51
sc_node$hjust<-ifelse(sc_node$angle>90|sc_node$angle< -90,0,1)
sc_node$angle<-ifelse(sc_node$angle< -90,sc_node$angle+180,sc_node$angle)
sc_node$angle<-ifelse(sc_node$angle> 90,sc_node$angle-180,sc_node$angle)
sc_node$size<-as.data.frame(table(sc.data.nkt$cluster_annot_nkt_label)[match(sc_node$node,levels(sc.data.nkt$seurat_clusters))])$Freq

scgraph<-graph_from_data_frame(sc_edge,vertices = sc_node)

p.2a3<-ggraph(scgraph,layout = "dendrogram",circular=T)+
  geom_edge_elbow(color="grey",lineend = "square",linejoin = "mitre")+
  geom_node_text(aes(x=x*1.15,y=y*1.15,label=label,filter=leaf,color=cell_type,
                     angle=angle,hjust=hjust),size=2)+
  geom_node_point(aes(x=x*1.07,y=y*1.07,filter=leaf,color=cell_type,size=size),alpha=0.6)+
  # scale_colour_manual(values= rep( brewer.pal(9,"Paired") , 30)) +
  scale_size_continuous( range = c(0.5,4) ) +
  theme_void() +
  theme(
    legend.position="none",
    plot.margin = unit(c(0,0,0,0),"cm")
    # plot.margin=unit(c(1,1,1,1),"cm"),
  ) +
  expand_limits(x = c(-2, 2), y = c(-2, 2))


p.2b1<-DimPlot(sc.data.nkt,group.by = "Sev_Outcome",
              pt.size = 0.1,
              raster = T,
              repel = T)+
  scale_color_manual(labels=function(x){gsub("_","\n",x)},
                     values = c("#74C476","#FEB24C","#E31A1C","#800026"))+
  labs(title = "")+
  theme(panel.background = element_blank(),
        plot.margin = unit(c(0,0,0,0),"cm"),
        axis.line = element_line(),
        plot.title = element_text(hjust = 0.5,vjust = 0.5,size=10,face = "bold"),
        legend.position = "top",
        legend.justification = "center",
        legend.text.align = 0,
        legend.text = element_text(size=6),
        legend.key.size = unit(0.1,"lines"),
        legend.spacing = unit(0,"cm"),
        legend.margin = margin(0,0,0,0,"pt"),
        legend.box.margin = margin(0,0,0,0,"pt"),
        axis.title = element_text(size=8),
        axis.text = element_text(size=7))
p.2b1$layers[[1]]$geom$default_aes$alpha<-0.3
p.2b1$guides$colour$override.aes$alpha<-1

p.2b2<-show.cell.col(sc.meta.cell,
                    "Sev_Outcome",
                    "(NKT_CD4_TIM3_CD62L+NKT_CD8_TIM3_CD62L)/(T_cell_NKT)",
                    test.args = list(alternative="less"),
                    signif_ensemble = T,text_size = 2,vjust = 0.5,tip_length = 0.2,width = 0.1)+
  scale_y_continuous(expand=expansion(c(0,0.05)))+
  scale_fill_manual(labels=function(x){gsub("_","\n",x)},
                    guide=guide_legend(title = "COVID-19 Severity"),
                    values = c("#74C476","#FEB24C","#E31A1C","#800026"))+
  xlab("")+
  ylab("Tim-3+ NKT cell Percentage")+
  theme(panel.background = element_blank(),
        plot.margin = unit(c(0,0,0,0),"cm"),
        legend.position = "top",
        legend.box.just = "center",
        axis.line = element_line(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y  = element_text(size = 8),
        plot.title = element_text(hjust = 0.5,vjust = 0.5,size=10,face = "bold"))


p.2b3<-show.cell.col(sc.meta.cell,
                    "Sev_Outcome",
                    "NKT_CD8_TIM3_CD62L/(T_cell_NKT)",
                    test.args = list(alternative="less"),
                    signif_ensemble = T,text_size = 2,vjust = 0.5,tip_length = 0.2,width = 0.1)+
  scale_y_continuous(expand=expansion(c(0,0.05)))+
  scale_fill_manual(values = c("#74C476","#FEB24C","#E31A1C","#800026"))+
  xlab("")+
  ylab("Tim-3+ CD8+ NKT cell Percentage")+
  NoLegend()+
  theme(panel.background = element_blank(),
        plot.margin = unit(c(0,0,0,0),"cm"),
        axis.line = element_line(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y  = element_text(size = 8),
        plot.title = element_text(hjust = 0.5,vjust = 0.5,size=10,face = "bold"))

p.2b4<-show.cell.col(sc.meta.cell,
                    "Sev_Outcome",
                    "NKT_CD4_TIM3_CD62L/(T_cell_NKT)",
                    test.args = list(alternative="less"),
                    signif_ensemble = T,text_size = 2,vjust = 0.5,tip_length = 0.2,width = 0.1)+
  scale_y_continuous(expand=expansion(c(0,0.05)))+
  scale_fill_manual(values = c("#74C476","#FEB24C","#E31A1C","#800026"))+
  xlab("")+
  ylab("Tim-3+ CD4+ NKT cell Percentage")+
  NoLegend()+
  theme(panel.background = element_blank(),
        plot.margin = unit(c(0,0,0,0),"cm"),
        axis.line = element_line(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y  = element_text(size = 8),
        plot.title = element_text(hjust = 0.5,vjust = 0.5,size=10,face = "bold"))

myget_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

p.legend1<-myget_legend(p.2a1+scale_fill_manual(labels=function(x){gsub("TIM3.*","TIM3",x)},
                                               guide=guide_legend(title = "NKT subtype"),
                                               values = hue_pal()(6))+
                         theme(legend.key.height = unit(c(0.5),"cm"),
                               legend.key.width = unit(c(0.5),"cm"),
                               legend.text = element_text(size=9,vjust=0.5,hjust=0.5),
                               legend.title = element_text(size=9,hjust=0.5,vjust=0.5))
)

p.legend2<-myget_legend(p.2b1+
                        theme(legend.key.height = unit(c(0.5),"cm"),
                              legend.key.width = unit(c(0.5),"cm"),
                              legend.text = element_text(size=9,vjust=0.5,hjust=0.5),
                              legend.title = element_text(size=9,hjust=0.5,vjust=0.5))
)

aligned1<-cowplot::align_plots(p.2a1+NoLegend(),
                               p.2a2+NoLegend(),
                               p.2a3,
                               align = "hv",axis = "tblr")
aligned2<-cowplot::align_plots(p.2b1+NoLegend(),
                               p.2b2+NoLegend(),
                               p.2b3,
                               p.2b4,align="hv",axis="tblr")

p<-ggarrange(ggarrange(p.legend1,
                       ggarrange(cowplot::ggdraw(aligned1[[1]]),
                                 cowplot::ggdraw(aligned1[[2]]),
                                 cowplot::ggdraw(aligned1[[3]]),
                                 heights = c(3,2,4),
                                 nrow=3,hjust = 0,vjust = 1,align = "v",
                                 labels = c("Fig.2A","Fig.2B","Fig.2C"),
                                 common.legend = T,legend="none",font.label = list(size=10)),
                       nrow = 2,
                       heights = c(1,11)),
             ggarrange(p.legend2,
                       ggarrange(cowplot::ggdraw(aligned2[[1]]),
                                 cowplot::ggdraw(aligned2[[2]]),
                                 cowplot::ggdraw(aligned2[[3]]),
                                 cowplot::ggdraw(aligned2[[4]]),
                                 heights = c(3,2,2,2),
                                 ncol = 1,nrow=4,
                                 labels = c("Fig.2D","Fig.2E","Fig.2F","Fig.2G"),
                                 hjust = 0,vjust = 1,align = "hv",
                                 font.label = list(size=10)),
                       nrow=2,
                       heights = c(1,11)),
             ncol=2,
             vjust = 0,
             hjust=0
             )+
  theme(plot.margin = unit(c(1,1,1,1),units = "cm"))

pdf("Figure 2.pdf",width=10,height=13)
p
dev.off()

p.2s.list<-list()

p.theme<-theme(panel.background = element_blank(),
               axis.line = element_line(),
               axis.text.x = element_blank(),
               axis.ticks.x = element_blank(),
               axis.title.y  = element_text(size = 8),
               plot.title = element_text(hjust = 0.5,vjust = 0.5,size=10,face = "bold"))

p.2s.y<-c("NKT_DN_ITGAX/(T_cell_NKT+T_cell_gdT+T_cell_MAIT)","NKT_CD8/(T_cell_NKT+T_cell_gdT+T_cell_MAIT)",
          "NKT_CD4_CD40LG/(T_cell_NKT+T_cell_gdT+T_cell_MAIT)","NKT_CD8_CD62L/(T_cell_NKT+T_cell_gdT+T_cell_MAIT)")
p.2s.ylab<-c("NKT_DN_ITGAX NKT cell Percentage","NKT_CD8 NKT cell Percentage",
             "NKT_CD4_CD40LG NKT cell Percentage","NKT_CD8_CD62L NKT cell Percentage")

for(i in 1:length(p.2s.y))
{
  p.2s.list[[i]]<-show.cell.col(sc.meta.cell,
                                "Sev_Outcome",
                                p.2s.y[i],
                                signif_ensemble = T,text_size = 2,vjust = 0.6,tip_length = 0.15,width = 0.1)+
    scale_y_continuous(expand=expansion(c(0,0.05)))+
    scale_fill_manual(labels=function(x){gsub("_","\n",x)},
                      guide=guide_legend(title = "COVID-19 Severity"),
                      values = c("#74C476","#FEB24C","#E31A1C","#800026"))+
    xlab("")+
    ylab(p.2s.ylab[i])+
    NoLegend()+p.theme
}

p.s<-do.call(ggarrange,c(p.2s.list,
                         list(
                           ncol=2,nrow=2,
                           hjust=0,vjust=0,
                           common.legend=T)))

p.s<-p.s+theme(plot.margin = unit(c(1,1,1,1),units = "cm"))

pdf("Figure S2.pdf",width=8,height=8)
p.s
dev.off()

rm(p.2a1,p.2a2,p.2a3,p.2b1,p.2b2,p.2b3,p.2b4,aligned1,aligned2,p,p.legend1,p.legend2,sc_edge,sc_node,myleaves,scgraph,node_name,
   p.2s.list,p.s)



#==============================================================gene.cor function==========================
gene.cor<-function(data,gene1,gene2=NULL,method="pearson")
{
  if(!gene1%in%rownames(data))
  {
    stop("No gene name found\n")
  }
  if(is.null(gene2))
  {
    cor.a<-cor(t(data[gene1,]),t(data),method = "pearson")
  } else{
    if(!gene2%in%rownames(data))
    {
      stop("No gene name found\n")
    }
    cor.a<-cor(t(data[gene1,]),t(data[gene2,]),method = "pearson")
  }
  cor.a<-t(cor.a)
  cor.a<-cor.a[!is.na(cor.a),]
  cor.a<-cor.a[order(abs(cor.a),decreasing = T)]
  return(cor.a)
}


#==============================================================show.gene.cor function=======================
show.gene.cor<-function(data,gene1,gene2,method="pearson",filter=F,meta=NULL,show.formula=F,show.pvalue=T)
{
  meltP<-data[,c(gene1,gene2)]
  meltP<-na.omit(meltP)
  if(filter)
  {
    meltP<-meltP[meltP[,gene1]>0&meltP[,gene2]>0,]
  }
  if(!is.null(meta))
  {
    if(length(meta)!=nrow(meltP))
    {
      if(!is.null(meta)&!is.null(rownames(meltP)))
      {
        meltP<-meltP[rownames(meltP)%in%names(meta),]
        meta<-meta[names(meta)%in%rownames(meltP)]
      } else{
        meta<-NULL
      }
    }
  }
  p<-ggplot(meltP,
            aes(x=eval(parse(text=gene1)),
                y=eval(parse(text=gene2)),
                color=meta,
                group=meta
            )
  )+
    geom_point()+
    xlab(gene1)+
    ylab(gene2)+
    labs(title =sprintf("%s-%s Correlation",gene1,gene2))
  if(is.null(meta))
  {
    corr<-rcorr(as.matrix(meltP[,c(gene1,gene2)]),type = method)
    cor.r<-corr$r[1,2]
    cor.p<-corr$P[1,2]
  } else{
    meltP2<-meltP
    meltP2$meta<-meta
    cor.r<-c()
    cor.p<-c()
    for(i in levels(factor(meta)))
    {
      corr<-rcorr(as.matrix(meltP2[meltP2$meta==i,c(gene1,gene2)]),type = method)
      cor.r<-c(cor.r,corr$r[1,2])
      cor.p<-c(cor.p,corr$P[1,2])
    }
  }
  
  
  if(show.formula){
    if(is.null(meta))
    {
      p<-p+stat_poly_eq(aes(label = paste(..eq.label.., 
                                          paste("R==",signif(cor.r,3),sep = ""), 
                                          paste("p==",signif(cor.p,3),sep = ""), 
                                          sep = '~~~~')),
                        color="red",
                        formula = y ~ x, parse = T)
    } else{
      p<-p+stat_poly_eq(aes(label = paste(..eq.label.., 
                                          paste("R==",signif(cor.r,3),sep = ""), 
                                          paste("p==",signif(cor.p,3),sep = ""), 
                                          sep = '~~~~')), 
                        formula = y ~ x, parse = T)
      
    }
    
  } else{
    if(show.pvalue){
      if(is.null(meta))
      {
        p<-p+stat_poly_eq(aes(label = paste(..grp.label..,
                                            paste("R==",signif(cor.r,3),sep = ""), 
                                            paste("p==",signif(cor.p,3),sep = ""), 
                                            sep = '~~~~')), 
                          color="red",
                          formula = y ~ x, parse = T)

        
      } else{
        p<-p+stat_poly_eq(aes(label = paste(..grp.label..,
                                            paste("R==",signif(cor.r,3),sep = ""), 
                                            paste("p==",signif(cor.p,3),sep = ""), 
                                            sep = '~~~~')), 
                          formula = y ~ x, parse = T)
      }
      
    } else{
      if(is.null(meta))
      {
        p<-p+stat_poly_eq(aes(label = paste(..grp.label..,
                                            paste("R==",signif(cor.r,3),sep = ""), 
                                            sep = '~~~~')), 
                          color="red",
                          formula = y ~ x, parse = T)
        
      } else{
        p<-p+stat_poly_eq(aes(label = paste(..grp.label..,
                                            paste("R==",signif(cor.r,3),sep = ""), 
                                            sep = '~~~~')), 
                          formula = y ~ x, parse = T)
      }
    }
  }
    
  if(is.null(meta))
  {
    p<-p+geom_smooth(method = "lm",
                     se=F,
                     color="red")
  } else{
    p<-p+geom_smooth(method = "lm",
                     se=F)
  }
  p<-p+theme(
    panel.background = element_blank(),
    axis.line = element_line(),
    plot.title = element_text(hjust = 0.5,vjust = 0.5),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.key = element_blank()
  )
  return(p)
}

show.gene.cor2<-function(data1,data2,gene1,gene2,method="pearson",filter=F,meta=NULL)
{
  data<-data.frame(data1[,gene1],row.names = rownames(data1))
  colnames(data)<-gene1
  data<-data[rownames(data)%in%rownames(data2),,drop=FALSE]
  data<-data.frame(data[,gene1],data2[rownames(data2)%in%rownames(data),gene2])
  rownames(data)<-rownames(data1)
  meltP<-data
  if(filter)
  {
    meltP<-meltP[meltP[,1]>0&meltP[,2]>0,]
  }
  if(is.null(meta))
  {
    col<-"black"
  }else{
    col<-meta
  }
  p<-ggplot(meltP,
            aes(x=eval(parse(text=colnames(meltP)[1])),
                y=eval(parse(text=colnames(meltP)[2])),
                col=meta,
                group=meta
            )
  )+
    geom_point()+
    annotate("text",
             x=max(meltP[,1])*0.8,
             y=max(meltP[,2])*0.9,
             label=sprintf("r=%.3f",
                           cor(meltP[,1],meltP[,2,],method = method)
             )
    )+
    stat_poly_eq(aes(label = paste(..eq.label.., 
                                   sprintf("R==%.3g",sign(..r.squared..)*sqrt(abs(..r.squared..))), 
                                   ..p.value.label..,
                                   sep = '~~~~')), 
                 formula = y ~ x, parse = T)+
    xlab(gene1)+
    ylab(gene2)+
    labs(title =sprintf("%s-%s Correlation",gene1,gene2))
  if(is.null(meta))
  {
    p<-p+geom_smooth(method = "lm",
                     se=F,
                     color="red")
  } else{
    p<-p+geom_smooth(method = "lm",
                     se=F)
  }
  p<-p+theme(
    panel.background = element_blank(),
    axis.line = element_line(),
    plot.title = element_text(hjust = 0.5,vjust = 0.5),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.key = element_blank()
  )
  return(p)
}


#==============================================================Figure 3=========================


tim3.marker<-FindMarkers(sc.data.nkt[,sc.data.nkt$sampleID%in%sc.meta.cell$Sample.name[sc.meta.cell$.Datasets%in%c("Batch02","Batch07")]],
                         ident.1 = levels(sc.data.nkt$cluster_annot_nkt_plot)[grep("TIM3",levels(sc.data.nkt$cluster_annot_nkt_plot))],
                         ident.2 = levels(sc.data.nkt$cluster_annot_nkt_plot)[-grep("TIM3",levels(sc.data.nkt$cluster_annot_nkt_plot))],
                         group.by = "cluster_annot_nkt_plot",
                         slot = "counts",
                         assay = "RNA",
                         logfc.threshold=0,
                         min.pct=0.02,
                         test.use = "wilcox")

gene_fc<-tim3.marker$avg_log2FC
names(gene_fc)<-rownames(tim3.marker)
gene_fc<-gene_fc[order(gene_fc,decreasing = T)]
gene_fc<-na.omit(gene_fc)

gsea_result<-gseGO(gene_fc,
                   ont="BP",
                   OrgDb = org.Hs.eg.db,
                   keyType = "SYMBOL",
                   minGSSize = 4,
                   maxGSSize = 2000,
                   eps = 0,
                   pvalueCutoff = 1,
                   nPermSimple=100000)

p.3a<-gseaplot2(gsea_result,"GO:0006915")+
  annotate("text",x=0.85,y=0.85,
           label = sprintf("p=%.3g\np.adj=%.3g",
                           gsea_result@result$pvalue[gsea_result@result$ID=="GO:0006915"],
                           gsea_result@result$p.adjust[gsea_result@result$ID=="GO:0006915"]),
           hjust=0,vjust=0.6)+
  annotate("text",x=0.65,y=0.85,
           label = "apoptotic process",
           hjust=0,vjust=0,
           color="green",
           size=5)+theme(plot.margin = unit(c(1,1,1,1),units = "cm"))
  

p.3b<-show.cell.col(data.frame(NKT_subgroup=sc.data.nkt$cluster_annot_nkt,
                               MT=sc.data.nkt$percent.mt
),
"NKT_subgroup",
"MT",
# comparison = group.combine(levels(sc.data.nkt$cluster_annot_nkt_plot))[grep("TIM3",group.combine(levels(sc.data.nkt$cluster_annot_nkt_plot)))],
comparison1 = c("NKT_CD4_TIM3_CD62L","NKT_CD8_TIM3_CD62L"),
comparison2 = levels(sc.data.nkt$cluster_annot_nkt_plot),
signif_level = T,
decreasing = T,
signif_ensemble = T,
height.step = 0.5,text_size = 2,vjust=0.5,tip_length = 0.2
)+
  scale_y_continuous(expand=expansion(c(0,0.05)))+
  labs(title = "")+
  xlab("")+
  ylab("MT gene percentage")+
  theme(panel.background = element_blank(),
        axis.line = element_line(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5,vjust = 0.5,size=10,face = "bold"))


a<-data.frame(percent.mt=sc.data$percent.mt,
              sample.id=sc.data$sampleID) %>%
  group_by(sample.id) %>%
  dplyr::summarise(mt_median=median(percent.mt),
            mt_total=length(percent.mt))

mt.ratio<-data.frame(percent.mt=sc.data.nkt$percent.mt,
                     sample.id=sc.data.nkt$sampleID,
                     nkt=sc.data.nkt$cluster_annot_nkt) %>%
  group_by(sample.id,nkt) %>%
  dplyr::summarise(mt_total=length(percent.mt),
            mt_high=length(percent.mt[percent.mt>2*a$mt_median[match(sample.id,a$sample.id)]]),
            # mt_high=length(percent.mt[percent.mt>7]),
            mt_ratio=mt_high/mt_total)


mt.ratio<-as.data.frame(mt.ratio)


p.3c<-show.cell.col(mt.ratio,
                    "nkt",
                    "mt_ratio",
                    signif_level = T,
                    decreasing = T,
                    signif_ensemble = T,
                    height.step = 1,text_size = 2,vjust=0.5,tip_length = 0.2
)+
  scale_y_continuous(expand=expansion(c(0,0.05)))+
  NoLegend()+
  labs(title = "")+
  xlab("")+
  ylab("MT gene > 2Med percentage")+
  theme(panel.background = element_blank(),
        axis.line = element_line(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5,vjust = 0.5,size=10,face = "bold"))

meta.batch<-c("Batch02","Batch07","Batch15")

p.theme<-theme(panel.background = element_blank(),
               axis.line = element_line(),
               legend.position = "top",
               plot.title = element_text(hjust = 0.5,vjust = 0.5),
               legend.key.size = unit(c(0.2,0.2),"cm"),
               legend.background = element_blank(),
               legend.key = element_blank(),
               legend.spacing = unit(c(0.1),"cm"),
               legend.text = element_text(size=5,vjust = 0.5),
               legend.title = element_text(size=6,vjust = 0.85),
               axis.text.x = element_text(size = 8))

p.3d<-plot_cell_trajectory(cds,color_by="Pseudotime",cell_size = 0.1, size=1,show_backbone=TRUE)+
  theme(legend.position = "top",
        legend.title = element_text(size = 8))

p.3e<-plot_cell_trajectory(cds,color_by="cluster_annot_nkt_plot",cell_size = 0.1, size=1,show_backbone=TRUE)+
  scale_color_manual(labels=function(x){gsub("TIM3.*","TIM3",x)},
                    guide=guide_legend(title = "",override.aes = list(size=2)),
                    values = hue_pal()(6))
  # guides(color=guide_legend(title = "NKT cell subset",override.aes = list(size=2)))+
  theme(legend.position = "top",
        legend.title = element_text(size=8),
        legend.text = element_text(size=6),
        legend.key.size = unit(0.1,"lines"),
        legend.spacing = unit(0,"cm"),
        legend.margin = margin(0,0,0,0,"pt"),
        legend.box.margin = margin(0,0,0,0,"pt"))

p.3f.list<-list()

gene_int<-c("CASP2","CASP3","CASP6","CASP9")
gene_int_name<-gene_int

meta.batch<-c("Batch02","Batch07","Batch15","Batch01","Batch03","Batch07")

for(i in 1:length(gene_int))
{
  p.3f.list[[i]]<-show.cell.gene.outcome.col(data=rep(sc.data.nkt@assays$RNA@data[gene_int[i],sc.data.nkt$batch%in%meta.batch],time=2),
                                             target.gene = gene_int_name[i],
                                             cell.subgroup = c(sc.data.nkt$cluster_nkt_tim3[sc.data.nkt$batch%in%meta.batch],
                                                               rep("NKT",time=length(sc.data.nkt$cluster_nkt_tim3[sc.data.nkt$batch%in%meta.batch]))),
                                             outcome = rep(sc.meta$Sev_Outcome[match(sc.data.nkt$sampleID[sc.data.nkt$batch%in%meta.batch],sc.meta$Sample.name)],time=2),
                                             inner.p = F,show.allpvalue = F,text_size = 2,show.point = F,
                                             signif_ensemble = T,height.step = 0.08,tip_length = 0,
                                             sp.p=c(3)
  )+
    scale_y_continuous(expand=expansion(c(0,0.05)))+
    scale_fill_manual(values = c("#74C476","#FEB24C","#E31A1C","#800026"))+
    xlab("")+
    p.theme
  
}


myget_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

p.legend1<-myget_legend(p.3b+scale_fill_manual(labels=function(x){gsub("TIM3.*","TIM3",x)},
                                               guide=guide_legend(title = "NKT subtype"),
                                               values = hue_pal()(6))+
                          theme(legend.key.height = unit(c(0.5),"cm"),
                                legend.key.width = unit(c(0.5),"cm"),
                                legend.position = "top",
                                legend.justification = "center",
                                legend.text = element_text(size=9,vjust=0.5,hjust=0.5),
                                legend.title = element_text(size=9,hjust=0.5,vjust=0.5))
)

p.legend2<-myget_legend(p.3f.list[[1]]+
                          scale_fill_manual(values = c("#74C476","#FEB24C","#E31A1C","#800026"),
                                            labels=function(x){gsub("_","\n",x)},
                                            guide=guide_legend(title="COVID-19 Severity"))+
                          theme(legend.key.height = unit(c(0.5),"cm"),
                                legend.key.width = unit(c(0.5),"cm"),
                                legend.text = element_text(size=9,vjust=0.5,hjust=0.5),
                                legend.title = element_text(size=9,hjust=0.5,vjust=0.5))
)

p<-ggarrange(ggarrange(p.3a,
                       ggarrange(p.legend1,
                                 ggarrange(p.3b+NoLegend(),p.3c,nrow=2,hjust=0,vjust=0,
                                           labels = c("Fig.3B","Fig.3C"),align = "hv",
                                           font.label = list(size=10)),
                                 nrow=2,
                                 heights = c(2,11)),
                       ncol = 2,
                       widths = c(2,1),hjust = 0,vjust = 0,labels = c("Fig.3A",""),
                       font.label = list(size=10)),
             ggarrange(ggarrange(p.3d,p.3e,nrow = 2,hjust=0,vjust=0,labels = c("Fig.3D","Fig.3E",align="hv"),
                                 font.label = list(size=10)),
                       ggarrange(p.legend2,
                                 do.call(ggarrange,c(p.3f.list,list(ncol=2,
                                                                    nrow=2,
                                                                    hjust=0,
                                                                    vjust=0,
                                                                    align="hv",
                                                                    common.legend=T,
                                                                    legend="none"))),
                                 nrow = 2,
                                 heights = c(1,11),
                                 hjust = 0,vjust = 0,labels = c("Fig.3F",""),
                                 font.label = list(size=10)),
                       widths = c(1,2)),
             nrow = 2,
             hjust = 0,vjust = 0)+
  theme(plot.margin = unit(c(1,1,1,1),units = "cm"))
  


pdf("Figure 3.pdf",width=16,height=16)
p
dev.off()

rm(p.3a,p.3b,p.3c,p.3d,p.3e,p.3f.list,p.legend1,p.legend2,p.theme,p,a,mt.ratio,meta.batch)



#==============================================================Figure 4=========================
p.theme<-theme(panel.background = element_blank(),
        axis.line = element_line(),
        legend.position = "top",
        plot.title = element_text(hjust = 0.5,vjust = 0.5),
        legend.key.size = unit(c(0.2,0.2),"cm"),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.spacing = unit(c(0.1),"cm"),
        legend.text = element_text(size=5,vjust = 0.5),
        legend.title = element_text(size=6,vjust = 0.85),
        axis.text.x = element_text(size = 8))

p.4.list<-list()

gene_int<-c("PDCD1","CTLA4","LAG3","IFNG","IL10","IL4","KLRK1","PRF1","GZMB")
gene_int_name<-c("PD-1","CTLA4","LAG3","IFNG","IL10","IL4","NKG2D","PRF1","GZMB")

meta.batch<-c("Batch02","Batch07","Batch15")

for(i in 1:length(gene_int))
{
  p.4.list[[i]]<-show.cell.gene.outcome.col(data=rep(sc.data.nkt@assays$RNA@data[gene_int[i],sc.data.nkt$batch%in%meta.batch],time=2),
                                            target.gene = gene_int_name[i],
                                            cell.subgroup = c(sc.data.nkt$cluster_nkt_tim3[sc.data.nkt$batch%in%meta.batch],
                                                              rep("NKT",time=length(sc.data.nkt$cluster_nkt_tim3[sc.data.nkt$batch%in%meta.batch]))),
                                            outcome = rep(sc.meta$Sev_Outcome[match(sc.data.nkt$sampleID[sc.data.nkt$batch%in%meta.batch],sc.meta$Sample.name)],time=2),
                                            inner.p = F,show.allpvalue = F,text_size = 2,show.point = F,
                                            height.step = 0.08,tip_length = 0,
                                            signif_ensemble = T,
                                            pvalue_cut = c(0.01,1e-3,1e-4,1e-10,1e-20,1e-50,1e-90),
                                            sp.p=c(3)
  )+
    scale_y_continuous(expand=expansion(c(0,0.05)))+
    scale_fill_manual(values = c("#74C476","#FEB24C","#E31A1C","#800026"))+
    xlab("")+
    p.theme
}

p.legend<-myget_legend(p.4.list[[1]]+
                         scale_fill_manual(values = c("#74C476","#FEB24C","#E31A1C","#800026"),
                                           guide=guide_legend(title="COVID-19 Severity"),
                                           labels=function(x){gsub("_","\n",x)})+
                         theme(legend.key.width = unit(0.8,"cm"),
                               legend.key.height = unit(0.8,"cm"),
                               legend.text = element_text(size=10,vjust=0.5,hjust=0.5),
                               legend.title = element_text(size=12,hjust=0.5,vjust=0.5))
                       )

p<-do.call(ggarrange,c(p.4.list,list(ncol=3,
                                     nrow=3,
                                     hjust=0,
                                     vjust=0,
                                     align="hv",
                                     common.legend=T,
                                     legend="none",
                                     labels=paste("Fig.4",LETTERS[1:length(p.4.list)],sep = ""),
                                     font.label = list(size=10))))


pdf("Figure 4.pdf",width=12,height=12)
ggarrange(p.legend,
          p,
          nrow = 2,
          heights = c(1,11))+
  theme(plot.margin = unit(c(1,1,1,1),"cm"))
dev.off()

rm(p.4.list,p.legend,p.theme)





#==============================================================Figure 5=========================
p.5a<-show.cell.col(data.frame(Outcome=sc.data@assays$RNA@data["IL12A",sc.data$cluster_tree_annot_merge%in%c("Monocyte_CD14","Monocyte_CD16","mono_DC","pDC")],
                               IL12A=sc.data$Sev_Outcome[sc.data$cluster_tree_annot_merge%in%c("Monocyte_CD14","Monocyte_CD16","mono_DC","pDC")]),
                    "IL12A",
                    "Outcome",
                    signif_level = T,
                    signif_ensemble = T,text_size = 3,vjust = 0.5,tip_length = 0.1,width = 0.1,height.step = 0.4
)+
  scale_y_continuous(expand=expansion(c(0,0.05)))+
  scale_fill_manual(values = c("#74C476","#FEB24C","#E31A1C","#800026"),
                    guide=guide_legend(title = "COVID-19 Severity"))+
  xlab("")+
  ylab("IL12 in DC/Monocyte")+
  theme(panel.background = element_blank(),
        axis.line = element_line(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5,vjust = 0.5,size=10,face = "bold"),
        axis.text.x.bottom = element_blank())
  

data1<-data.frame(gene=(sc.data@assays$RNA@data["IL12B",grep("DC|Monocyte",sc.data$cluster_tree_annot_merge)]+
                          sc.data@assays$RNA@data["IL12A",grep("DC|Monocyte",sc.data$cluster_tree_annot_merge)])>0,
                  sampleID=sc.data$sampleID[grep("DC|Monocyte",sc.data$cluster_tree_annot_merge)])
data1<-table(data1)
data1.ratio<-data1["TRUE",]/(data1["FALSE",]+data1["TRUE",])

data2<-data.frame(gene=(sc.data@assays$RNA@data["HAVCR2",grep("NKT",sc.data$cluster_tree_annot_merge)])>0,
                  sampleID=sc.data$sampleID[grep("NKT",sc.data$cluster_tree_annot_merge)])
data2<-table(data2)
data2.ratio<-data2["TRUE",]/(data2["FALSE",]+data2["TRUE",])


meta<-sc.meta.cell$Sev_Outcome
names(meta)<-sc.meta.cell$Sample.name
meta<-meta[meta!="severe/critical_deceased"]

data1.ratio<-data1.ratio[names(data1.ratio)%in%names(meta)]
data2.ratio<-data2.ratio[names(data2.ratio)%in%names(meta)]

p.5b<-show.gene.cor(data.frame("IL12"=data1.ratio,"TIM3"=data2.ratio),
                    "IL12","TIM3",meta=meta,method = "spearman")+
  scale_color_manual(values = c("#74C476","#FEB24C","#E31A1C","#800026"))
p.5b<-p.5b+xlab("IL12 highly expressed DC/Monocyte percentage")+
  ylab("TIM3 highly expressed NKT cell percentage")+
  labs(title = "")+
  guides(color=guide_legend(title="COVID-19 Severity"))+
  theme(axis.title.x.bottom  = element_text(vjust=0,hjust=0.5),
        legend.position = "top"
        )



p<-ggarrange(p.5a,p.5b,ncol=2,
             common.legend = T,
             align = "h",
             labels = c("Fig.5A","Fig.5B"),
             hjust = 0,vjust = 1,
             font.label = list(size=10))+
  theme(plot.margin = unit(c(1,1,1,1),"cm"))


pdf("Figure 5.pdf",width=12,height=8)
p
dev.off()

rm(data1,data2,data1.ratio,data2.ratio,p.5a,p.5b)



