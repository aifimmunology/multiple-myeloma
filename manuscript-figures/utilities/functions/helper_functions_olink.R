plot_longitudinal_protein <- function(assay, dat, alpha=0.3, size=0.5){
  df = dat[Assay==assay]
  df[, NPX_combined := mean(NPX), by=list(Subject,visitDetails)]
  
  ggplot(df, aes(
    x=visitDetails,
    y=NPX_combined, 
    group = Subject,
    col=Subject)
  )+
   
    geom_boxplot(data=df,
                 aes(x=visitDetails,
                     y=NPX_combined,
                     group=visitDetails),
                 alpha=0.4)+
   geom_line(alpha=alpha)+ geom_point(size=size)+
    ggtitle(assay)+
    theme_minimal()+
    xlab('Visit')+
    ylab('NPX')
  
}

### write function to test proteins 
test_protein <- function(assay, Plot=F, mat){
  
  ### subset matrix to 
  ### include only protein of interest
  tmp = mat[Assay ==assay]
  
  ### filter to the columns of interest
  dat = tmp[,c('NPX','Subject','visitDetails')]
  
  ### it looks like there are 
  ### duplicate runs 
  ### per subject + visit combo
  ### so the choice is to take mean
  dat =reshape2::dcast(dat, Subject ~ visitDetails, fun.aggregate = mean, value.var='NPX')
  
  ### remove subjects that have no paired data
  dat = dat[!is.na(dat[,2]) & !is.na(dat[,3]),]
  
  ### test for differences 
  ### Using a MWU paired test 
  res = data.frame(
    Assay = assay,
    Pvalue = wilcox.test(dat[,3], dat[,2], paired=T)$p.value,
    Log2FC = median(dat[,3] - dat[,2]),
    N = nrow(dat)
  )
  
  ### add (optional) plotting function 
  if(Plot){
    df= melt(dat)
    p=ggplot(df,
           aes(x=variable,
               y=value,))+geom_boxplot()+geom_point()+
      geom_line(aes(x=variable,
                    y=value,
                    group=Subject))+
      ggpubr::stat_compare_means(method='wilcox.test', paired=T)+
      ggtitle(paste(res[1,]))
    #print(p)
    return(p)
    } else{ 
    return(res)
    }
}

#### summarize contrast 
summarize_contrast <- function(contrast_matrix, tName ='Transplant'){
  
  ### Run differential test 
  results <- lapply(proteins,
                              function(x)
                                try(test_protein(x,Plot = F, mat=contrast_matrix)))
  ### Remove proteins 
  ### that failed differential testing
  results = rbindlist(results[sapply(results, class) !='try-error'])
  
  ### transform p-values into q-values 
  results$qvalue = qvalue::qvalue(results$Pvalue)$qvalues
  results$contrast =tName
  
  ### Calculate # DEPs
  nhits = sum(results$qvalue < 0.1)
  
  ### plot histogram
  ### to show global behavior 
  p1= ggplot(results,
             aes(x=Pvalue))+geom_histogram()+
    ggtitle('Pvalue Distribution')
  
  ### show q-value results 
  p2= ggplot(results,
             aes(x=Log2FC,
                 y= -log10(qvalue),
                 col = qvalue < 0.1,
                 label=Assay)) +geom_point()+
    
    ggrepel::geom_text_repel(size=4)+
    ggtitle(paste(tName,':', nhits, 'hits'))+theme_classic()+
    theme(legend.position = 'none', text=element_text(size=25))
  
  ### combine into 1
  #ggpubr::ggarrange(p1,p2, ncol=2)
  print(p2)
  
  return(results)
}
