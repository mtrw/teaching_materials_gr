
##### FIND:SETUP #########################################################################################################

#send them home
setwd("~")

#shared library location
ll <- "/filer/projekte/ggr-course2020/r-libs/"

while ( ! require(plyr,lib.loc = ll ) ) {
  install.packages("plyr",lib = ll )
}
while ( ! require(dplyr,lib.loc = ll ) ) {
  install.packages("dplyr",lib = ll )
}
while ( ! require(data.table,lib.loc = ll ) ) {
  install.packages("data.table",lib = ll )
}
while ( ! require(magrittr,lib.loc = ll ) ) {
  install.packages("magrittr",lib = ll )
}
while ( ! require(ggplot2,lib.loc = ll ) ) {
  install.packages("ggplot2",lib = ll )
}
while ( ! require(parallel,lib.loc = ll ) ) {
  install.packages("parallel",lib = ll )
}
while ( ! require(colorspace,lib.loc = ll ) ) {
  install.packages("colorspace",lib = ll )
}
while ( ! require(plyr,lib.loc = ll ) ) {
  install.packages("plyr",lib = ll )
}
while ( ! require(jpeg,lib.loc = ll ) ) {
  install.packages("jpeg",lib = ll )
}
while( ! require(gtools,lib.loc = ll ) ){
  install.packages("gtools",lib = ll )
}
# while( ! require(googlesheets,lib.loc = ll ) ){
#   install.packages("googlesheets",lib = ll )
# }

####### FIND:STANDARD FUNCTIONS  #########################################################################################################

#minor_allele_frequency(sample(1,20,r=T)) #BY STATE FOR HAPLOIDS, ONLY ... NOT FOR 0,1,2 ENCODING OF DIPLOIDS
minor_allele_frequency <- function(x){
  tbl <- table(x)
  if(length(tbl)==1){ return(as.numeric(0)) }
  min <- tbl[order(tbl)][1]
  if ((min / sum(tbl))==1) { browser() }
  as.numeric(min / sum(tbl))
}


#minor_allele(sample(1:3,20,r=T))
minor_allele <- function(x){
  tbl <- table(x)
  ma <- names(tbl[order(tbl)])[1]
  return(as(ma,class(x)))
}


#be tolerant of non-utf-8 characters when grepping
Sys.setlocale('LC_ALL','C')

#create an empty plot with ranges x=c(low,high) and y=ditto
null_plot <- function(x,y,xlab=NA,ylab=NA,...){
  plot(NULL,xlim=range(x),ylim=range(y),xlab=xlab,ylab=ylab,...)
}

#turn levels of a factor into colours from a colorspace palette (in the diverge_hcl set)
replace_levels_with_colours <- function(x,palette="Berlin",alpha=1,fun="diverge_hcl",plot=FALSE){
  require(colorspace)
  n <- nu(x)
  cols <- match.fun(fun)(n,palette = palette,alpha = alpha)
  colvec <- swap( x , unique(x) , cols , na.replacement = NA )
  if(plot==FALSE) {
    return(colvec)
  } else {
    # null_plot(y=1:length(cols),x=rep(1,length(cols)),xaxt="n",yaxt="n")
    # text(y=1:length(cols),x=rep(1,length(cols)),labels=unique(x),col=cols)
    null_plot(y=1,x=1,xaxt="n",yaxt="n",bty="n")
    legend(x=1,legend=unique(x),fill=cols,text.col=cols)
  }
}

#turn levels of a factor into colours from a colorspace palette (in the diverge_hcl set)
replace_scale_with_colours <- function(x,palette="ag_GrnYl",fun="sequential_hcl",alpha=1,plot=FALSE){
  require(colorspace)
  s <- x %>% scale_between(1,100) %>% round
  cols <- match.fun(fun)(100,palette = palette,alpha = alpha)
  if(plot==FALSE) {
    return(cols[s])
  } else {
    plot(y = (1:100)%>%scale_between(min(x),max(x)),x=rep(0,100),ylab=NA,xlab=NA,xaxt="n",col=cols)
  }
}

#display brewer palletes and remind yourself how to use them
cs <- function(){
  library(colorspace)
  hcl_palettes(plot = TRUE)
  ce("Example: diverge_hcl(30,palette=\"Berlin\")")
  ce("demoplot(x, type = c(\"map\", \"heatmap\", \"scatter\", \"spine\", \"bar\", \"pie\",
     \"perspective\", \"mosaic\", \"lines\"), ...)")
}


#apply a fun of two vars in every combo, and give the results as a matrix
#grid_ply(1:2,3:1,sum)
grid_ply <- function(rows,cols,FUN) {
  lapply( rows , function(row) lapply(cols , function(col) {
    FUN(row,col)
  }) %>% unlist ) %>% unlist %>%
    matrix(nrow=length(rows),dimnames=list(rows,cols))
}

#apply a fun of two vars in every combo, and give the results as a matrix
#mc_grid_ply(1:2,3:1,sum)
mc_grid_ply <- function(rows,cols,FUN,cores=25,...) {
  mclapply( mc.cores=cores , rows , function(row) lapply(cols , function(col) {
    FUN(row,col)
  }) %>% unlist ) %>% unlist %>%
    matrix(nrow=length(rows),dimnames=list(rows,cols))
}

#scale a list of values to between two points, proportionally spaced as they were originally
#rnorm(100) %>% scale_between(20,29) %>% pd
scale_between <- function(x,lower,upper){
  if(all(x==0)) return(x)
  ( x - min(x) ) / (max(x)-min(x)) * (upper-lower) + lower
}

#easy way to see the spread of your values. plot the density.
#pd(c(NA,rnorm(500),NA))
pd <- function(x,...){
  x %>% density(na.rm=TRUE,...) %>% plot()
}

#difference between two values but vectorises more intuitively than diff()
#difff(10,-286)
difff <- function(a,b){
  abs(a-b)
}

#z-transform values in a vec
#z_transform(x = ((runif(100))**2)+20  ) %>% pd
z_transform <- function(x){
  if ( sd(x)==0 ){
    return(rep(0,times=length(x)))
  }
  (mean(x)-x) / sd(x)
}

#number of unique entries in vector
#nu(c(1,1,2,3,4,4))
nu <-function(x){
  unique(x) %>% length
}

#fetch the upper triangular matrix in a vector without all that typing
#upper_tri( matrix(1:9,ncol=3) , diag=T )
upper_tri <- function(mat,diag=F){
  mat[upper.tri(mat,diag=diag)]
}

#parallel mean of two numbers. Great for BLAST search (midpoint of hits = pmean(start_vec,end_vec))
#pmean(1:5,2:6,5:1)
pmean <- function(...){
  invecs <- list(...)
  out <- rep(0, times=length(invecs[[1]]) )
  lapply(invecs,function(x) out <<- out+x )
  out/length(invecs)
}
pmean2 <- pmean #legacy reasons



#simultaneously change the names of things to a particular thing if they match a particular string.
#name_by_match(vec=iris$Species,matches = c("set","sicol"),names = c("SETOSA","VERSICOLOR"))
name_by_match <- function(vec,matches,names){
  if(length(matches) != length(names)) { stop("Lengths of `matches` and `names` vectors don't match, you massive knob!")  }
  vec %<>% as.character()
  l_ply(1:length(matches) , function(n){
    vec[grep(matches[n],vec)] <<- names[n]
  })
  vec
}

swap <- function(vec,matches,names,na.replacement=NA){
  orig_vec <- vec
  #if(sum(! matches %in% names ) > 0 ) { stop("Couldn't find all matches in names") }
  if(length(matches) != length(names)) { stop("Lengths of `matches` and `names` vectors don't match, you old bison!") }
  if(is.factor(vec)) { levels(vec) <- c(levels(vec),names,na.replacement) }
  vec[is.na(orig_vec)] <- na.replacement
  l_ply( 1:length(matches) , function(n){
    vec[orig_vec==matches[n]] <<- names[n]
  })
  vec
}

#length of a vector
#(c(3,4))
vec_length <-function(x){
  ( sum( x**2 ) )**(1/2)
}

#euc dists between points described in lists of coords
euc_dist<-function(x1,x2,y1,y2){
  sqrt( ((x1-x2)**2) + ((y1-y2)**2) )
}


#as the title says. relative probs should be positive (duh) and same order as events (duh)
#t <- random_draws_from_any_discreet_distribution(n=50000,events=LETTERS[1:5],relative_probs=1:5) %>% table; t / max(t) * 5 ; rm(t)
random_draws_from_any_discreet_distribution <- function(n=1,events,relative_probs){
  lapply( 1:n,   function(x){
    events[ which( runif(1)*sum(relative_probs) < (cumsum(relative_probs)) )[1] ]    
  }) %>% unlist
}

#n random (non-repeated!) rows of a data frame
#sample_df(iris,20)
sample_df <- function(df,n=10,...){
  df[ sample.int(nrow(df),n,...) ]
}






#echo to the global environment. good warning messager. still doesn't work in mclappy, hashtag doh
#ce("beans ",list("hello "," goodbye")," whatever")
ce <- function(...){   cat(paste0(...,"\n"), sep='', file=stderr()) %>% eval(envir = globalenv() ) %>% invisible() }



u <- function(...){
  unique(...)
}



wait <- function(message="Press [enter] to continue"){
  invisible(readline(prompt=message))
}

count_NAs <- function(x){
  sum(is.na(x))
}



##### FIND:DATASETS #########################################################################################################


genotype_data_1 <- readRDS("/filer/projekte/ggr-course2020/data/simpop_1.Rds")
positions <- c(262,669, 632 , 62 , 59, 268 ,244 ,103, 452 ,117, 636 ,557 ,685 ,209, 719, 488, 969 , 76 ,529, 866)
speciess <- c("Rye","Barley","Wheat")
genotype_data_1[,position := positions[.GRP] ,by=.(locus)]
genotype_data_1[,species := speciess[.GRP] ,by=.(population)]
setnames(genotype_data_1,"gt","genotype")
genotype_data_1[,locus := NULL]
genotype_data_1[,generation := NULL]
genotype_data_1[,sample := paste0("sample_",.GRP),by=.(individual,population)]
genotype_data_1[,individual := NULL]
genotype_data_1[,population := NULL]

height_means <- c(180,125,90)
height_sds <- c(20,15,5)
phenotype_data_1 <- genotype_data_1[,.( sample=unique(sample) , height = rnorm(nu(sample),height_means[.GRP],height_sds[.GRP])),by=.(species)]

selfings <- c( F,F,F,F,F,T,F,T,    F,T,T,F,T,T,T,F,T,T,F,T,T,T,F,F,T,T,F,T,   T,T,T,T,T,T,T,T,T,T,T,T,T,T,T )
phenotype_data_1[, selfing := selfings ]

yield_means <- c(2.1,2.8,4)
yield_sds <- c(2,1,0.5)/4
phenotype_data_1 <- phenotype_data_1[, yield := rnorm(.N,z_transform(-height),0.5) %>% scale_between(yield_means[.GRP]-(2.5*yield_sds[.GRP]),yield_means[.GRP]+(2.5*yield_sds[.GRP])) ,by=.(species)][]
#ggplot( data=phenotype_data_1 , mapping = aes( x=height , y=yield , colour=species ) ) + geom_point()
phenotype_data_1 <- phenotype_data_1[sample(1:.N)]

rm(positions,speciess,height_means,height_sds,selfings,yield_means,yield_sds)

population_1 <- readRDS("/filer/projekte/ggr-course2020/data/simpop_2.Rds")

# #phenotypes for KASP dataset
# p1 <- fread("/filer/projekte/ggr-course2020/data/BRIDGE_core200_w_photo_passportdata.csv",select=c(2,3),col.names=c("sample","annuality"))
# p1[ , sample := sub(" ","_",sample) ]
# p2 <- fread("/filer/projekte/ggr-course2020/data/BRIDGE_core200_w_photo_phenodata.csv" , select=c(1,10,11,12,16) , col.names=c("sample","spike_length","awn_length","awn_texture","rachilla_hairs") )
# p2[ , sample := sub(" ","_",sample) ] 



#phenotypes for KASP dataset
p1 <- fread("/filer/projekte/ggr-course2020/data/BRIDGE_core1000_passport.csv",select=c(2,3),col.names=c("sample","annuality"))
p1[ , sample := sub(" ","_",sample) ]
p2 <- fread("/filer/projekte/ggr-course2020/data/BRIDGE_core1000.csv" , select=c(1,10,11,12,16) , col.names=c("sample","spike_length","awn_length","awn_texture","rachilla_hairs") )
p2[ , sample := sub(" ","_",sample) ] 
phenotype_data <- p1[p2,on="sample"]

key_set <- paste0("HOR_",  c( 3728 , 13924 , 3912   , 1384 ,  3686  ,  3926 , 15898 , 2403  , 5876 , 1251 , 7428 , 11922 ))
extra_set <- paste0("HOR_",c( 12311 , 10702  , 7474 ,  16078 ,  11448 ,  13724    ,    13716   ,   14914   , 11875 ,  7172 , 7552  ))
suggested_samples <- c(key_set,extra_set)


phenotype_data <- data.table(
  sample = c(
    "HOR_1251","HOR_11875","HOR_11922","HOR_7428","HOR_7552","HOR_14914","HOR_13716", "HOR_13724",
    "HOR_3728","HOR_16078","HOR_13924","HOR_3912","HOR_1384","HOR_3686","HOR_12311","HOR_10702","HOR_7474",
    "HOR_3926","HOR_15898","HOR_2403","HOR_5876","HOR_11448","HOR_7172"
  ),
  latitude = (l <- c(
    rnorm(8,40,2),
    rnorm(9,30,2),
    rnorm(6,50,2)
  )  %>% round()),
  cold_tolerance_LT50 = l %>% scale_between(-5,-9)  %>% `+`(rnorm(23,0,0.15)) %>% round(digits=1)
)[phenotype_data,on="sample"]
phenotype_data[ , cold_tolerance_LT50 := ifelse(annuality=="winter_type",cold_tolerance_LT50 - 2,cold_tolerance_LT50 + 2) ]
#




















################## KASP FUNCTIONS ################################################################################
read_KASP_output <- function(fname="/filer/projekte/ggr-course2020/data/KASP*"){
  fread(cmd=paste0('cat ',fname,' | grep "\\S" | awk \'/Well/ {p=1} /SDS/ {p=0} p==1 && !/Well/ {print $2,$3,$4,$5}\''),col.names=c("sample","marker_id","fluorescence_allele_1","fluorescence_allele_2"))
}
kasp <- read_KASP_output()

read_KASP_genotypes <- function(){
  g <- fread("/filer/projekte/ggr-course2020/data/KASP.csv",blank.lines.skip = TRUE)[marker != ""]
  g[, genotype := as.integer(genotype) ]
  g[]
}


plot_KASP_output <- function(kasp,marker,colour_samples=FALSE){
  k <- copy(kasp)
  k <- k[marker_id==marker]
  
  slist <- NA
  if(colour_samples==TRUE) { slist <- unique(k$sample) }
  
  l_ply( slist , function(samp){
    k[,colour:="black"]
    k[sample=="NTC",colour := "grey" ]
    if(!is.na(samp)){
      k[sample==samp,colour:="red"]
      title2 <- paste0("; Red: sample ",samp,")")
    } else {
      title2 <- ")"
    }
    title <- paste0("KASP fluorescence, marker: ",marker,"\n(Grey: controls",title2)
    if(!is.null(dev.list())) {dev.off() %>% invisible}
    plot(
      x=k$fluorescence_allele_1,
      y=k$fluorescence_allele_2,
      col=k$colour,
      pch=20,
      xlab="Fluorescence allele 1 (FU)",
      ylab="Fluorescence allele 2 (FU)",
      main=title,
      xlim=range(k$fluorescence_allele_1) + c(-.2,.2),
      ylim=range(k$fluorescence_allele_2) + c(-.2,.2),
    )
    text(
      x=k$fluorescence_allele_1,
      y=k$fluorescence_allele_2,
      labels=k$sample,
      col=k$colour,
      cex=.5,
      adj=1.1
    )
    if (!is.na(samp)) {wait(message = "Press <enter> to colour the next sample ...")}
  })
}
# plot_KASP_output(kasp=kasp,marker="SNP_01")
# plot_KASP_output(kasp=kasp,marker="SNP_01",colour_samples = TRUE)


##############FIND:COURSE FUNCTIONS ####################################################################################################################################################

#customisable pop sim (finite sites, two states)

mutate <- function(gt,mut_rate=0.2){
  r <- runif(length(gt)) < mut_rate
  gt[r] <- sample(0:1L,sum(r),replace=T)
  gt
}

recombine <- function(gt1,gt2){
  if(runif(1)>=.5){
    gt1 -> t
    gt2 -> gt1
    t -> gt2
  }
  bp <- sample(1:(length(gt1)-1),1) #breakpoint
  #if(bp<0){browser()}
  c(gt1[1:bp],gt2[(bp+1):length(gt2)])
}


make_starting_pop <- function(default_size,n_loci,mut_rate){
  p <- data.table(
    population = rep(1,default_size*n_loci),
    individual = rep(1:default_size,each=n_loci),
    locus = rep(1:n_loci,default_size)
  )
  
  starting_gt <- rep(0,n_loci)#sample(0:1L,n_loci,replace=T)
  p[ , gt := mutate(starting_gt,mut_rate) , by=.(individual) ]
  p[]
}




sim_pop <- function(n_gens=15,events,default_size=10,n_loci=30,mut_rate=.1,plot=TRUE){
  pop <- pop1 <- make_starting_pop(default_size,n_loci,mut_rate)
  pop[ , generation := 0 ]
  
  if(plot){
    pop_params <- sim_popsizes_before_running(events=events,default_size=default_size) #also checks for illegal or weird population requests
    max_size <- pop_params[1]
    max_pop <- pop_params[2]
    par(mfrow=c(1,1),pch=20)
    null_plot(x=0:n_gens,y=c(0.4,max_pop+0.6),xlab="Generation",ylab="Population",lab=c(n_gens,max_pop,1))
  }
  
  l_ply(1:n_gens,function(gen){
    ce("Generation: ",gen)
    pop <<- generation(pop,events[generation==gen],mut_rate,plot_pop_gen_maxpop=c(gen,max_size))[ , generation := gen ][]
    #pop
  })
  pop
}
#scenario1 <- sim_pop(n_gens=100,events=events,default_size = 20,mut_rate = .1)



#scenario=combined_data; colour_variable="rachilla_hairs"; colour_variable="latitude"
plot_genetic_diversity <- function(scenario,colour_variable=NULL){
  sim <- FALSE
  if(all(c("sample","marker","genotype") %in% colnames(scenario))){
    p <- scenario[ , .(individual=sample,locus=marker,gt=genotype) ]
  } else {
    p <- scenario[generation==max(generation)]
    p[,colvar := population]
    sim <- TRUE
  }
  
  if(is.null(colour_variable)){
    if(sim==FALSE){p[,colvar := 1]}
  } else if (!colour_variable %in% colnames(scenario)){
    ce("Warning: Colour variable \"",colour_variable,"\" not found in data")
    colour_variable <- NULL
    p[,colvar := 1]
    
  } else {
    p[,colvar:=scenario[,colour_variable,with=FALSE]]
  }
  if(!is.null(dev.list())) {dev.off() %>% invisible}
  par(pch=20)
  layout( mat = matrix(c(1,2,5,3,4,5),nrow=2,byrow=TRUE) , heights = c(1,1) , widths = c(1,1,0.5)  )
  pp <- dcast( p , colvar + individual ~ locus , value.var = "gt" )
  pm <- pp[,-c("colvar","individual")] %>% as.matrix
  #rownames(pm) <- paste0(colvarlist <- paste0("colvar:",pp$colvar),"_",indlist <- paste0("individual:",pp$individual))
  badcols <- apply(pm,2,function(x) nu(x)==1 )
  pm <- pm[,!badcols]
  pca <- prcomp(pm,scale. = T)
  pcs <- (pca$x %>% as.data.table)[ , colvar := pp$colvar][ , individual := pp$individual][]
  colours <- rep(1,nrow(pcs))
  if(!is.null(colour_variable)){
    if((is.integer(pcs$colvar) | is.numeric(pcs$colvar)) & length(unique(pcs$colvar))>5){
      colours <- replace_scale_with_colours(pcs$colvar,palette="ag_GrnYl",fun="sequential_hcl")
    } else {
      colours <- replace_levels_with_colours(pcs$colvar,palette="ag_GrnYl",fun="sequential_hcl")
    }
  }
  lims <- max(abs(c(pcs$PC1,pcs$PC2,pcs$PC3)))
  lims <- c(lims,-lims)
  plot(x=1:ncol(summary(pca)$importance),y=summary(pca)$importance[2,],ylab="Proportion of variance explained",xlab="Principal Component")
  plot( pcs$PC2 , pcs$PC1 , cex=1.5 , col=colours , xlab="PC2" , ylab="PC1" , xlim = lims , ylim = lims )
  #text(pcs$PC2 , pcs$PC1 , cex=1.5 , labels=pcs$individual , col=colours )
  plot( pcs$PC3 , pcs$PC1 , cex=1.5 , col=colours , xlab="PC3" , ylab="PC1" , xlim = lims , ylim = lims )
  plot( pcs$PC2 , pcs$PC3 , cex=1.5 , col=colours , xlab="PC2" , ylab="PC3" , xlim = lims , ylim = lims )
  if(!is.null(colour_variable)){
    if((is.integer(pcs$colvar) | is.numeric(pcs$colvar)) & length(unique(pcs$colvar))>5){
      replace_scale_with_colours(pcs$colvar,palette="ag_GrnYl",fun="sequential_hcl",plot=TRUE)
    } else {
      replace_levels_with_colours(x = pcs$colvar,palette="ag_GrnYl",fun="sequential_hcl",plot=TRUE)
    }
  }
}



#graphical GTs

plot_simulated_genotypes <- function(scenario,show_populations=FALSE){
  if(!is.null(dev.list())) {dev.off() %>% invisible}
  par(mfrow=c(1,1),pch=20)
  s <- scenario[generation==max(generation)]
  s[,individual := .GRP , by=.(individual,population)]
  null_plot(x=range(s$locus),y=range(s$individual),xlab="Locus",ylab="Individual")
  shapes <- rep(20,length(unique(s$population)))
  if( show_populations ){ shapes <- frank(s$population,ties.method="dense") + 15 }
  points(x=s$locus,y=s$individual,col=replace_levels_with_colours(s$gt),pch=shapes)
}


sim_popsizes_before_running <- function(events,default_size=100){
  popsizes <- c(default_size)
  max_size <- 0
  max_pop <- 0
  
  for(i in sort(events$generation)){
    event <- events[generation==i]
    if(nrow(event)>1){
      stop("One event per generation only please!")
    }
    
    if(event$event=="size_change") {
      popsizes[event$population] <- event$size
    }
    if(event$event=="split"){
      popsizes[event$population] <- popsizes[event$population]-event$size
      if(event$population==length(popsizes)){
        popsizes <- c(popsizes,event$size)
      } else {
        popsizes <- c(popsizes[0:event$population],event$size,popsizes[(event$population+1):length(popsizes)])
      }
    }
    if (event$event=="merge"){
      popsizes[event$population] <- popsizes[event$population] - event$size
      popsizes[event$population2] <- popsizes[event$population2] + event$size
    }
    
    popsizes <- popsizes[popsizes != 0]
    #ce(popsizes)
    if(sum(popsizes < 0)>0){
      stop(paste("Event(s) result in negative population sizes!"))
    }
    max_size <- max(max_size,max(popsizes))
    max_pop <- max(max_pop,length(popsizes))
  }
  ce("Events table checked, no errors found ...")
  c(max_size,max_pop)
}



generation <- function(pop,event,mut_rate=.2,plot_pop_gen_maxpop=NULL){
  #pop muct be ordered for ind then for locus (so genotype is comparable when grouping on ind)
  #browser()
  inds_prev <- unique(pop$individual)
  n_inds_prev <- length(inds_prev)
  n_pops_prev <- nu(pop$population)
  
  popsizes <- ps <- pop[ , .(size=nu(individual)) , by=.(population,new_pop=as.numeric(population))]
  
  if(nrow(event)>0){
    if (event$event=="size_change") {
      ce("\t","Event! (",event$event,", size=",event$size,")")
      popsizes[population==event$population]$size <- event$size #check also continues all other populations onwards
    }
    if (event$event=="split"){
      ce("\t","Event! (",event$event,", size=",event$size,")")
      popsizes <- rbind(
        popsizes[population!=event$population],
        popsizes[population==event$population , .(population,new_pop=population,size=size-event$size) ],
        popsizes[population==event$population , .(population,new_pop=population+0.1,size=event$size)  ]
      )
    }
    if ( event$event=="merge" ){
      #rowser()
      ce("\t","Event! (",event$event,", size=",event$size,")")
      popsizes <- rbind(
        popsizes[population!=event$population], #uninvolved pops
        popsizes[population==event$population , .(population,new_pop=event$population2,size=event$size) ], #mergers
        popsizes[population==event$population , .(population,new_pop=population,size=size-event$size) ]#remainder from merger
      )
    }
    setkey(popsizes,new_pop)
    popsizes[,new_pop := as.numeric(.GRP) ,by=.(new_pop)]
  }
  
  #select individuals from old pop fisher-wright style (two per new pop member), recombine them to make one, add to new pop
  a <- popsizes[,{
    ldply(1:sum(size),function(j){
      from <- population
      to <- new_pop
      indlist <- pop[population==from]$individual %>% unique
      s1 <- sample(indlist,1)
      s2 <- sample(indlist,1)
      r <- recombine( pop[population==from & individual==s1]$gt , pop[population==from & individual==s2]$gt ) %>% mutate(mut_rate)
      if(!all(!is.na(r))){browser()}
      
      if(!is.null(plot_pop_gen_maxpop)){
        plot_pop_gen <- plot_pop_gen_maxpop[1]
        maxpop <- plot_pop_gen_maxpop[2]
        nf <- length(indlist)
        y1a <- from + (s1/nf)*(nf/maxpop)*0.8 - ((nf/maxpop)*0.8)/2
        y1b <- from + (s2/nf)*(nf/maxpop)*0.8 - ((nf/maxpop)*0.8)/2
        y2 <- to + (j/size)*(size/maxpop)*0.8 - ((size/maxpop)*0.8)/2
        points(plot_pop_gen,y2,cex=.5,pch=20)
        lines(x=c(plot_pop_gen-1,plot_pop_gen),y=c(y1a,y2),col="#00000022")
        lines(x=c(plot_pop_gen-1,plot_pop_gen),y=c(y1b,y2),col="#00000022")
      }
      data.table(
        population,
        individual = j + ((.GRP/10) - 0.1 ),
        locus = 1:length(r),
        gt = r
      )
    })
  } , by=.(population,new_pop) ] %>% setDT
  a[, .(individual=frank(individual,ties.method="dense"),locus,gt),by=.(population=new_pop)]
}

# events1 <- data.table(
#   event = c("split","merge","split","merge"),
#   size = c(3,1,5,3),
#   generation = c(1,2,10,17),
#   population = c(1,2,1,2),
#   population2 = c(NA,1,NA,3)
# )
# 
# scenario1 <- sim_pop( n_gens = 100 , n_loci = 20 , default_size = 10 , events = events1 , mut_rate = .3 )
# 
# plot_genetic_diversity( scenario1 )
# plot_genotypes( scenario1 )
# 
# 
# events2 <- data.table(
#   event = c("split","size_change","size_change"),
#   size = c(5,3,45),
#   generation = c(2,3,30),
#   population = c(1,2,2)
# )
# scenario2 <- sim_pop( default_size = 30 , n_gens=30 , events=events2 , mut_rate = .05 )
# plot_genetic_diversity( scenario2 )

fit_polynomial <- function(x,y,degree=2){
  X <- seq(from=min(x),to=max(x),length.out=100)
  l <- lm( y ~ poly(x,degree) )
  plot(x,y,pch=20)
  lines( X ,  predict(l,newdata = data.frame(x=X) ) )
  print(summary(l))
  invisible(l)
}




img <- jpeg::readJPEG("/filer/projekte/ggr-course2020/data/rye_height.jpg",native = T)
n <- 20 #must be above 9
eff_fertiliser <- 2
eff_subpopulation <- 10
eff_SNP1 <- 12
eff_SNP2 <- 0
err_sd <- 1

set.seed(21)
d <- data.table(
  fertiliser = abs(rnorm(n,20,7)),
  subpopulation = (g<-rep(0:1,times=c(round(n/3),n-round(n/3)))),
  SNP2 = ifelse(g==0,0,1),
  SNP1 = sample(0:1,n,replace=TRUE)
)

d[, eff_fertiliser := (eff_fertiliser*fertiliser) ]
d[, eff_subpopulation := (eff_subpopulation*subpopulation) ]
d[, eff_SNP1 := (eff_SNP1*SNP1)]
d[, eff_SNP2 := (eff_SNP2*SNP2)]
d[, plant_height :=  eff_fertiliser + eff_subpopulation + eff_SNP1 +  eff_SNP2 + rnorm(.N,err_sd) ]
d[3,SNP2:=1]
d[round(n/3)+4,SNP2:=0]
set.seed(NULL)


check_fit_fertiliser <- function(data=d,guess_eff_fertiliser,image=img){
  #test!
  d <- data
  d[, plant_height_guess := guess_eff_fertiliser*fertiliser ]
  d[ , idx := 1:.N]
  if(!is.null(dev.list())) {dev.off() %>% invisible}
  par(mfrow=c(2,1),mar=c(4,4,1,1),oma=c(0,0,0,0))
  null_plot(0:nrow(d),c(0:30,range(d$plant_height),range(d$plant_height_guess)),xlab="Sample",ylab="Plant height",xaxt="n")
  l_ply(1:nrow(d),function(r){
    rasterImage(image,r-0.2,0,r+0.2,d$plant_height[r])
  })
  
  points(1:nrow(d),d$plant_height_guess,pch=45,cex=3)
  
  cols <- diverge_hcl(3,"berlin")
  d[, {
    ys <- c(0 , cumsum( c(guess_eff_fertiliser*fertiliser) ) )
    x_int <- 0.5/(length(ys))
    xs <- idx + x_int + c(0,rep(-x_int,length(ys))) %>% cumsum %>% rev
    
    l_ply(1:(length(ys)-1),function(i){
      if(ys[i] != ys[i+1]){
        arrows(xs[i],ys[i],xs[i],ys[i+1],length = .05,col=cols[i],lwd=2) #ys[i+1 to give head-to-tail bent arrows]
      } else {
        points(xs[i],ys[i],pch=20,cex=.2)
      }
      
    })
  }, by=.(idx) ]
  barplot(d$plant_height_guess-d$plant_height,ylab="Residuals",ylim=c(-25,25))
  text(n/2,-12,labels = paste("Total squared residials:",sum((d$plant_height-d$plant_height_guess)**2)))
  legend(20,legend=c("Fertiliser effect"),text.col=cols[1],lty=c(1),col=cols[1])
  sum((d$plant_height-d$plant_height_guess)**2)
}

# check_fit_fertiliser(
#   guess_eff_fertiliser = 2.7 #black
# )



check_fit_fertiliser_SNP2_subpopulation <- function(data=d,guess_eff_fertiliser,guess_eff_SNP2,guess_eff_subpopulation,image=img){
  #test!
  d <- data
  d[, plant_height_guess := guess_eff_fertiliser*fertiliser + guess_eff_SNP2*SNP2 + guess_eff_subpopulation*subpopulation ]
  d[ , idx := 1:.N]
  if(!is.null(dev.list())) {dev.off() %>% invisible}
  par(mfrow=c(2,1),mar=c(4,4,1,1),oma=c(0,0,0,0))
  null_plot(0:nrow(d),c(0:30,range(d$plant_height),range(d$plant_height_guess)),xlab="Sample",ylab="Plant height",xaxt="n")
  l_ply(1:nrow(d),function(r){
    rasterImage(image,r-0.2,0,r+0.2,d$plant_height[r])
  })
  
  points(1:nrow(d),d$plant_height_guess,pch=45,cex=3)
  
  cols <- diverge_hcl(3,"berlin")
  d[, {
    ys <- c(0 , cumsum( c(guess_eff_fertiliser*fertiliser,guess_eff_SNP2*SNP2,guess_eff_subpopulation*subpopulation ) ) )
    x_int <- 0.5/(length(ys))
    xs <- idx + x_int + c(0,rep(-x_int,length(ys))) %>% cumsum %>% rev
    
    l_ply(1:(length(ys)-1),function(i){
      if(ys[i] != ys[i+1]){
        arrows(xs[i],ys[i],xs[i],ys[i+1],length = .05,col=cols[i],lwd=2) #ys[i+1 to give head-to-tail bent arrows]
      } else {
        points(xs[i],ys[i],pch=20,cex=.2)
      }
      
    })
  }, by=.(idx) ]
  barplot(d$plant_height_guess-d$plant_height,ylab="Residuals",ylim=c(-25,25))
  text(n/2,-12,labels = paste("Total squared residials:",sum((d$plant_height-d$plant_height_guess)**2)))
  legend(20,legend=c("Fertiliser effect","SNP2 effect","Subpopulation effect"),text.col=cols[1:3],lty=c(1,1,1),col=cols[1:3])
  sum((d$plant_height-d$plant_height_guess)**2)
}
# check_fit_fertiliser_SNP2_population(
#   data=d,
#   guess_eff_fertiliser=0,
#   guess_eff_SNP2=10,
#   guess_eff_subpopulation=10
# )

check_fit_fertiliser_SNP1_subpopulation <- function(data=d,guess_eff_fertiliser,guess_eff_SNP1,guess_eff_subpopulation,image=img){
  #test!
  d <- data
  d[, plant_height_guess := guess_eff_fertiliser*fertiliser + guess_eff_SNP1*SNP1 + guess_eff_subpopulation*subpopulation ]
  d[ , idx := 1:.N]
  if(!is.null(dev.list())) {dev.off() %>% invisible}
  par(mfrow=c(2,1),mar=c(4,4,1,1),oma=c(0,0,0,0))
  null_plot(0:nrow(d),c(0:30,range(d$plant_height),range(d$plant_height_guess)),xlab="Sample",ylab="Plant height",xaxt="n")
  l_ply(1:nrow(d),function(r){
    rasterImage(image,r-0.2,0,r+0.2,d$plant_height[r])
  })
  
  points(1:nrow(d),d$plant_height_guess,pch=45,cex=3)
  
  cols <- diverge_hcl(3,"berlin")
  d[, {
    ys <- c(0 , cumsum( c(guess_eff_fertiliser*fertiliser,guess_eff_SNP1*SNP1,guess_eff_subpopulation*subpopulation ) ) )
    x_int <- 0.5/(length(ys))
    xs <- idx + x_int + c(0,rep(-x_int,length(ys))) %>% cumsum %>% rev
    
    l_ply(1:(length(ys)-1),function(i){
      if(ys[i] != ys[i+1]){
        arrows(xs[i],ys[i],xs[i],ys[i+1],length = .05,col=cols[i],lwd=2) #ys[i+1 to give head-to-tail bent arrows]
      } else {
        points(xs[i],ys[i],pch=20,cex=.2)
      }
      
    })
  }, by=.(idx) ]
  barplot(d$plant_height_guess-d$plant_height,ylab="Residuals",ylim=c(-25,25))
  text(n/2,-12,labels = paste("Total squared residials:",sum((d$plant_height-d$plant_height_guess)**2)))
  legend(20,legend=c("Fertiliser effect","SNP1 effect","Subpopulation effect"),text.col=cols[1:3],lty=c(1,1,1),col=cols[1:3])
  sum((d$plant_height-d$plant_height_guess)**2)
}
# check_fit_fertiliser_SNP1_population(
#   data=d,
#   guess_eff_fertiliser=0,
#   guess_eff_SNP1=10,
#   guess_eff_subpopulation=10
# )

check_fit_fertiliser_SNP1 <- function(data=d,guess_eff_fertiliser,guess_eff_SNP1,image=img){
  #test!
  d <- data
  d[, plant_height_guess := guess_eff_fertiliser*fertiliser + guess_eff_SNP1*SNP1 ]
  d[ , idx := 1:.N]
  if(!is.null(dev.list())) {dev.off() %>% invisible}
  par(mfrow=c(2,1),mar=c(4,4,1,1),oma=c(0,0,0,0))
  null_plot(0:nrow(d),c(0:30,range(d$plant_height),range(d$plant_height_guess)),xlab="Sample",ylab="Plant height",xaxt="n")
  l_ply(1:nrow(d),function(r){
    rasterImage(image,r-0.2,0,r+0.2,d$plant_height[r])
  })
  
  points(1:nrow(d),d$plant_height_guess,pch=45,cex=3)
  
  cols <- diverge_hcl(3,"berlin")
  d[, {
    ys <- c(0 , cumsum( c(guess_eff_fertiliser*fertiliser,guess_eff_SNP1*SNP1 ) ) )
    x_int <- 0.5/(length(ys))
    xs <- idx + x_int + c(0,rep(-x_int,length(ys))) %>% cumsum %>% rev
    
    l_ply(1:(length(ys)-1),function(i){
      if(ys[i] != ys[i+1]){
        arrows(xs[i],ys[i],xs[i],ys[i+1],length = .05,col=cols[i],lwd=2) #ys[i+1 to give head-to-tail bent arrows]
      } else {
        points(xs[i],ys[i],pch=20,cex=.2)
      }
      
    })
  }, by=.(idx) ]
  barplot(d$plant_height_guess-d$plant_height,ylab="Residuals",ylim=c(-25,25))
  text(n/2,-12,labels = paste("Total squared residials:",sum((d$plant_height-d$plant_height_guess)**2)))
  legend(20,legend=c("Fertiliser effect","SNP1 effect"),text.col=cols[1:2],lty=c(1,1),col=cols[1:2])
  sum((d$plant_height-d$plant_height_guess)**2)
}
# check_fit_fertiliser_SNP1(
#   data=d,
#   guess_eff_fertiliser=5,
#   guess_eff_SNP1=10
# )

check_fit_fertiliser_SNP2 <- function(data=d,guess_eff_fertiliser,guess_eff_SNP2,image=img){
  #test!
  d <- data
  d[, plant_height_guess := guess_eff_fertiliser*fertiliser + guess_eff_SNP2*SNP2 ]
  d[ , idx := 1:.N]
  if(!is.null(dev.list())) {dev.off() %>% invisible}
  par(mfrow=c(2,1),mar=c(4,4,1,1),oma=c(0,0,0,0))
  null_plot(0:nrow(d),c(0:30,range(d$plant_height),range(d$plant_height_guess)),xlab="Sample",ylab="Plant height",xaxt="n")
  l_ply(1:nrow(d),function(r){
    rasterImage(image,r-0.2,0,r+0.2,d$plant_height[r])
  })
  
  points(1:nrow(d),d$plant_height_guess,pch=45,cex=3)
  
  cols <- diverge_hcl(3,"berlin")
  d[, {
    ys <- c(0 , cumsum( c(guess_eff_fertiliser*fertiliser,guess_eff_SNP2*SNP2 ) ) )
    x_int <- 0.5/(length(ys))
    xs <- idx + x_int + c(0,rep(-x_int,length(ys))) %>% cumsum %>% rev
    
    l_ply(1:(length(ys)-1),function(i){
      if(ys[i] != ys[i+1]){
        arrows(xs[i],ys[i],xs[i],ys[i+1],length = .05,col=cols[i],lwd=2) #ys[i+1 to give head-to-tail bent arrows]
      } else {
        points(xs[i],ys[i],pch=20,cex=.2)
      }
      
    })
  }, by=.(idx) ]
  barplot(d$plant_height_guess-d$plant_height,ylab="Residuals",ylim=c(-25,25))
  text(n/2,-12,labels = paste("Total squared residials:",sum((d$plant_height-d$plant_height_guess)**2)))
  legend(20,legend=c("Fertiliser effect","SNP2 effect"),text.col=cols[1:2],lty=c(1,1),col=cols[1:2])
  sum((d$plant_height-d$plant_height_guess)**2)
}
# check_fit_fertiliser_SNP2(
#   data=d,
#   guess_eff_fertiliser=1,
#   guess_eff_SNP2=2
# )


ce("\n-------------------------------------------\nGENETIC RESOURCES 2020 COURSE SETUP SCRIPT COMPLETE!\n------------------------------------------\n")
























# #misc: pipetting exercise
# a <- 1*(1/2)**c(0:4)
# b <- rev(a)
# m <- pmin(a,b)
# M <- pmax(a,b)
# volmin <- 30
# vm <- rep(volmin,length(a))
# vM <- (m/M)*volmin
# vt <- vM + vm
# cm <- (m*vm)/vt
# cM <- (M*vM)/vt
# #check
# cm == cM
# cf <- min(cM)
# v_add <- ((cM*vt)-(cf*vt))/cf
# #check:
# cm * ( vt / (vt + v_add) )
