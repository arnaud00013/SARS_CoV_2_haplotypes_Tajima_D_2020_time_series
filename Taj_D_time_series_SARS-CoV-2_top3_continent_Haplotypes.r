#Time code execution
start_time <- Sys.time()

#import libraries
library("ggplot2")
library("seqinr")
library("grid")
library("RColorBrewer")
library("randomcoloR")
library("gplots")
library("lmPerm")
library("ggpubr")
library("gridExtra")
library("RColorBrewer")
library("tidyr")
library("Cairo")
library("parallel")
library("foreach")
library("doParallel")
library("FD")
library("vegan")
library("plyr")
library("lme4")
library("lmerTest")
library("MuMIn")
library("AICcmodavg")
library("EnvStats")
library("session")
#import script arguments
#ABSOLUTE Path of the folder containing the data 
output_workspace <- as.character(commandArgs(TRUE)[1]) 

v_lst_filenames_table <- c("Table_time_series_Taj_D_with_resamplings_first_wave_Africa_top3_haplotypes.csv","Table_time_series_Taj_D_with_resamplings_first_wave_Asia_top3_haplotypes.csv","Table_time_series_Taj_D_with_resamplings_first_wave_Europe_top3_haplotypes.csv","Table_time_series_Taj_D_with_resamplings_first_wave_North_America_top3_haplotypes.csv","Table_time_series_Taj_D_with_resamplings_first_wave_South_America_top3_haplotypes.csv","Table_time_series_Taj_D_with_resamplings_first_wave_Oceania_top3_haplotypes.csv","Table_time_series_Taj_D_with_resamplings_second_wave_Africa_top3_haplotypes.csv","Table_time_series_Taj_D_with_resamplings_second_wave_Asia_top3_haplotypes.csv","Table_time_series_Taj_D_with_resamplings_second_wave_Europe_top3_haplotypes.csv","Table_time_series_Taj_D_with_resamplings_second_wave_North_America_top3_haplotypes.csv","Table_time_series_Taj_D_with_resamplings_second_wave_South_America_top3_haplotypes.csv","Table_time_series_Taj_D_with_resamplings_second_wave_Oceania_top3_haplotypes.csv")
palette_haplotypes <-c("#003FFF", "#28573F","#FF0000","#FF964F", "#FFCAE4", "#CB5590", "#A6D4FF", "#AD26FF", "#A9A9A9", "#B3F396", "#40E0D0", "#8F1C55","#03C0A4", "#FFD700","#000000","#009900","#E3B778")
names(palette_haplotypes) <- c("I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV","XV","XVI","XVII")
#Set language as English for date formatting
Sys.setlocale("LC_ALL","English")

#import metadata 
df_metadata_samples <- read.csv2(file = paste0(output_workspace,"HussinGroup_Interdb_metadata.csv"),sep = ",",header = TRUE,stringsAsFactors = FALSE)
df_metadata_samples$X <- NULL
rownames(df_metadata_samples) <- df_metadata_samples$GISAID_id
df_metadata_samples <- subset(df_metadata_samples,(!is.na(date))&(!is.na(seq_country)))
#focus on human hosted viral sequences without too much ambiguous calls
df_metadata_samples <- subset(df_metadata_samples,(host_species=="Human")&(seqdata=="yes"))

v_seq_country_to_seq_continent <- unique(df_metadata_samples[,c("seq_country","seq_continent")])[,"seq_continent"]
names(v_seq_country_to_seq_continent) <- unique(df_metadata_samples[,c("seq_country","seq_continent")])[,"seq_country"]

#Bootstrap for estimating Tajima's D
df_time_period_Taj_D_with_resamplings <- NULL
for (current_fn in v_lst_filenames_table){
  df_time_period_Taj_D_with_resamplings <- rbind(df_time_period_Taj_D_with_resamplings,read.csv2(file = paste0(output_workspace,current_fn),sep = ",",header = TRUE,stringsAsFactors = FALSE))
}
df_time_period_Taj_D_with_resamplings$Tajima_D <- as.numeric(df_time_period_Taj_D_with_resamplings$Tajima_D)

time_period_length <- 1 #number of month(s)

#function for plotting linear model
ggplotRegression <- function (fit,ggsave_path,the_filename,xlabl=NA,ylabl=NA) {
  library(ggplot2)
  bool_gg_save <- TRUE
  if(is.na(xlabl)){
    xlabl <- names(fit$model)[2]
  }
  if(is.na(ylabl)){
    ylabl <- names(fit$model)[1]
  }
  adj_r_sq <- formatC(summary(fit)$adj.r.squared, format = "e", digits = 3)
  slope <-formatC(summary(fit)$coefficients[,1][2], format = "e", digits = 3)
  p_val <- ifelse(test = "Iter"%in%colnames(summary(fit)$coefficients),yes=ifelse(test = unname(summary(fit)$coefficients[,3][2])<(2E-4),yes = "<2E-4",no = formatC(unname(summary(fit)$coefficients[,3][2]), format = "e", digits = 3)),no=formatC(broom::glance(fit)$p.value, format = "e", digits = 3))
  tryCatch(expr = {ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) +
      geom_point() +
      stat_smooth(method = "lm", col = "red") +
      xlab(xlabl)+
      ylab(ylabl)+
      labs(title = paste("Adj R2 = ",adj_r_sq,
                         " Slope =",slope,
                         " P =",p_val))+ theme(plot.title=element_text(hjust=0,size=12))},error=function(e) bool_gg_save <- FALSE)
  
  if (bool_gg_save){
    ggsave(filename = the_filename, path=ggsave_path, width = 15, height = 10, units = "cm")
  }else{
    print(paste0(the_filename, "won't be created because of it is irrelevant for gene in path ", ggsave_path))
  }
  #return result as the real float numbers
  adj_r_sq <- unname(summary(fit)$adj.r.squared)
  slope <-unname(summary(fit)$coefficients[,1][2])
  p_val <- ifelse(test = "Iter"%in%colnames(summary(fit)$coefficients),yes=ifelse(test = unname(summary(fit)$coefficients[,3][2])<(2E-4),yes = "<2E-4",no = unname(summary(fit)$coefficients[,3][2])),no=broom::glance(fit)$p.value)
  return(list(adj_r_sq_current_lm = adj_r_sq,slope_current_lm = slope,p_val_current_lm=p_val))
}

ggplotRegression_export_eps <- function (fit,ggsave_path,the_filename,xlabl=NA,ylabl=NA) {
  library(ggplot2)
  bool_gg_save <- TRUE
  if(is.na(xlabl)){
    xlabl <- names(fit$model)[2]
  }
  if(is.na(ylabl)){
    ylabl <- names(fit$model)[1]
  }
  adj_r_sq <- formatC(summary(fit)$adj.r.squared, format = "e", digits = 3)
  slope <-formatC(summary(fit)$coefficients[,1][2], format = "e", digits = 3)
  p_val <- ifelse(test = "Iter"%in%colnames(summary(fit)$coefficients),yes=ifelse(test = unname(summary(fit)$coefficients[,3][2])<(2E-4),yes = "<2E-4",no = formatC(unname(summary(fit)$coefficients[,3][2]), format = "e", digits = 3)),no=formatC(broom::glance(fit)$p.value, format = "e", digits = 3))
  tryCatch(expr = {ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) +
      geom_point() +
      stat_smooth(method = "lm", col = "red") +
      xlab(xlabl)+
      ylab(ylabl)+
      labs(title = paste("Adj R2 = ",adj_r_sq,
                         " Slope =",slope,
                         " P =",p_val))+ theme(plot.title=element_text(hjust=0,size=12))},error=function(e) bool_gg_save <- FALSE)
  
  if (bool_gg_save){
    ggsave(filename = the_filename, path=ggsave_path, width = 15, height = 10, units = "cm", device = cairo_ps)
  }else{
    print(paste0(the_filename, "won't be created because of it is irrelevant for gene in path ", ggsave_path))
  }
  #return result as the real float numbers
  adj_r_sq <- unname(summary(fit)$adj.r.squared)
  slope <-unname(summary(fit)$coefficients[,1][2])
  p_val <- ifelse(test = "Iter"%in%colnames(summary(fit)$coefficients),yes=ifelse(test = unname(summary(fit)$coefficients[,3][2])<(2E-4),yes = "<2E-4",no = unname(summary(fit)$coefficients[,3][2])),no=ifelse(test = unname(summary(fit)$coefficients[,4][2])<(2e-16),yes = "<2e-16",no = unname(summary(fit)$coefficients[,4][2])))
  return(list(adj_r_sq_current_lm = adj_r_sq,slope_current_lm = slope,p_val_current_lm=p_val))
}

##################
df_time_period_Taj_D_with_resamplings$start_time_period_formatted <- format(as.Date(df_time_period_Taj_D_with_resamplings$start_time_period,"%Y-%m-%d"),"%B %d, %Y")
df_time_period_Taj_D_with_resamplings_North_America_top3 <- subset(df_time_period_Taj_D_with_resamplings, ((seq_continent=="North America")))
df_time_period_Taj_D_with_resamplings_Europe_top3 <- subset(df_time_period_Taj_D_with_resamplings, ((seq_continent=="Europe")))
df_time_period_Taj_D_with_resamplings <- subset(df_time_period_Taj_D_with_resamplings, ((seq_continent=="North America")&(Hussin_name%in%c("VIII","II","IX"))) | ((seq_continent=="Europe")&(Hussin_name%in%c("VIII","XV","XIV"))) | (!seq_continent%in%c("Europe","North America")))

#plot Taj D time series
ggplot(data = df_time_period_Taj_D_with_resamplings,mapping = aes(x=factor(start_time_period_formatted,levels=as.character(format(sort(unique(as.Date(start_time_period))),"%B %d, %Y"))),y=Tajima_D)) + geom_boxplot(aes(col=Hussin_name)) + ylab("Tajima's D")+ xlab(paste0("start of the time period\n(1 month)")) + theme_bw() + theme(title =  element_text(size=12),axis.text.x = element_text(angle = 60,hjust = 1,size=6),axis.text.y = element_text(size=10),axis.title.x = element_text(size=12),axis.title.y = element_text(size=12),legend.position="right") + labs(col="Haplotype") + scale_color_manual(values = palette_haplotypes) + facet_wrap(~seq_continent,ncol=1)
ggsave(filename = paste0(time_period_length,"month_Taj_D_with_resampling_across_time_HAPLOTYPES_horizontal_comparison.png"), path=output_workspace, width = 9, height = 23, units = "cm",dpi = 1200)
ggsave(filename = paste0(time_period_length,"month_Taj_D_with_resampling_across_time_HAPLOTYPES_horizontal_comparison.svg"), path=output_workspace, width = 9, height = 23, units = "cm",dpi = 1200,device=svg)

ggplot(data = df_time_period_Taj_D_with_resamplings,mapping = aes(x=factor(start_time_period_formatted,levels=as.character(format(sort(unique(as.Date(start_time_period))),"%B %d, %Y"))),y=Tajima_D)) + geom_boxplot(aes(col=Hussin_name)) + ylab("Tajima's D")+ xlab(paste0("start of the time period\n(1 month)")) + theme_bw() + theme(title =  element_text(size=12),axis.text.x = element_text(angle = 60,hjust = 1,size=6),axis.text.y = element_text(size=10),axis.title.x = element_text(size=12),axis.title.y = element_text(size=12),legend.position="right") + labs(col="Haplotype") + scale_color_manual(values = palette_haplotypes) + facet_wrap(~seq_continent,nrow=1)
ggsave(filename = paste0(time_period_length,"month_Taj_D_with_resampling_across_time_HAPLOTYPES_vertical_comparison.png"), path=output_workspace, width = 35, height = 15, units = "cm",dpi = 1200)
ggsave(filename = paste0(time_period_length,"month_Taj_D_with_resampling_across_time_HAPLOTYPES_vertical_comparison.svg"), path=output_workspace, width = 35, height = 15, units = "cm",dpi = 1200,device=svg)

ggplot(data = df_time_period_Taj_D_with_resamplings,mapping = aes(x=factor(start_time_period_formatted,levels=as.character(format(sort(unique(as.Date(start_time_period))),"%B %d, %Y"))),y=Tajima_D)) + geom_boxplot(aes(col=Hussin_name)) + ylab("Tajima's D")+ xlab(paste0("start of the time period\n(1 month)")) + theme_bw() + theme(text= element_text(size=14),title =  element_text(size=14),axis.text.x = element_text(angle = 60,hjust = 1,size=10),axis.text.y = element_text(size=12),axis.title.x = element_text(size=14),axis.title.y = element_text(size=14),legend.position="right") + labs(col="Haplotype") + scale_color_manual(values = palette_haplotypes) + facet_grid(Hussin_name~seq_continent)
ggsave(filename = paste0("Grid",time_period_length,"month_Taj_D_with_resampling_across_time_HAPLOTYPES_vertical_comparison.png"), path=output_workspace, width = 35, height = 30, units = "cm",dpi = 1200)
ggsave(filename = paste0("Grid",time_period_length,"month_Taj_D_with_resampling_across_time_HAPLOTYPES_vertical_comparison.svg"), path=output_workspace, width = 35, height = 30, units = "cm",dpi = 1200,device=svg)

#plot North America top3 only
ggplot(data = df_time_period_Taj_D_with_resamplings_North_America_top3,mapping = aes(x=factor(start_time_period_formatted,levels=as.character(format(sort(unique(as.Date(start_time_period))),"%B %d, %Y"))),y=Tajima_D)) + geom_boxplot(aes(col=Hussin_name)) + ylab("Tajima's D")+ xlab(paste0("start of the time period\n(1 month)")) + theme_bw() + theme(title =  element_text(size=12),axis.text.x = element_text(angle = 60,hjust = 1,size=6),axis.text.y = element_text(size=10),axis.title.x = element_text(size=12),axis.title.y = element_text(size=12),legend.position="right") + labs(col="Haplotype") + scale_color_manual(values = palette_haplotypes) + facet_wrap(~seq_continent,nrow=1)
ggsave(filename = paste0(time_period_length,"month_Taj_D_with_resampling_across_time_top3_North_America_HAPLOTYPES_vertical_comparison.png"), path=output_workspace, width = 17.5, height = 15, units = "cm",dpi = 1200)
ggsave(filename = paste0(time_period_length,"month_Taj_D_with_resampling_across_time_top3_North_America_HAPLOTYPES_vertical_comparison.svg"), path=output_workspace, width = 17.5, height = 15, units = "cm",dpi = 1200,device=svg)

#plot Europe top3 only
ggplot(data = df_time_period_Taj_D_with_resamplings_Europe_top3,mapping = aes(x=factor(start_time_period_formatted,levels=as.character(format(sort(unique(as.Date(start_time_period))),"%B %d, %Y"))),y=Tajima_D)) + geom_boxplot(aes(col=Hussin_name)) + ylab("Tajima's D")+ xlab(paste0("start of the time period\n(1 month)")) + theme_bw() + theme(title =  element_text(size=12),axis.text.x = element_text(angle = 60,hjust = 1,size=6),axis.text.y = element_text(size=10),axis.title.x = element_text(size=12),axis.title.y = element_text(size=12),legend.position="right") + labs(col="Haplotype") + scale_color_manual(values = palette_haplotypes) + facet_wrap(~seq_continent,nrow=1)
ggsave(filename = paste0(time_period_length,"month_Taj_D_with_resampling_across_time_top3_Europe_HAPLOTYPES_vertical_comparison.png"), path=output_workspace, width = 17.5, height = 15, units = "cm",dpi = 1200)
ggsave(filename = paste0(time_period_length,"month_Taj_D_with_resampling_across_time_top3_Europe_HAPLOTYPES_vertical_comparison.svg"), path=output_workspace, width = 17.5, height = 15, units = "cm",dpi = 1200,device=svg)


library("lubridate")
df_time_period_Taj_D_with_resamplings$stop_time_period <- as.character(as.Date(df_time_period_Taj_D_with_resamplings$start_time_period) %m+% months(time_period_length) - 1)
detach("package:lubridate", unload=TRUE)

v_unique_haplotype_keys <- sort(unique(paste0(df_time_period_Taj_D_with_resamplings$start_time_period,"/",df_time_period_Taj_D_with_resamplings$stop_time_period,"/",df_time_period_Taj_D_with_resamplings$seq_continent,"/",df_time_period_Taj_D_with_resamplings$Hussin_name)))
v_unique_haplotype_keys_nb_seqs <- NULL
i=1
for (current_hap_key in v_unique_haplotype_keys){
  pos_sep <- gregexpr(pattern = "/",text = current_hap_key,fixed = T)[[1]]
  current_start <- substr(current_hap_key,1,pos_sep[1]-1)
  current_stop <- substr(current_hap_key,pos_sep[1]+1,pos_sep[2]-1)
  current_continent <- substr(current_hap_key,pos_sep[2]+1,pos_sep[3]-1)
  current_hap <- substr(current_hap_key,pos_sep[3]+1,nchar(current_hap_key))
  hap_nb_seq_to_add <- nrow(subset(df_metadata_samples,((as.Date(date)>=as.Date(current_start))&(as.Date(date)<=as.Date(current_stop))&(seq_continent==current_continent)&(Hussin_name==current_hap))))
  names(hap_nb_seq_to_add) <- current_hap_key
  v_unique_haplotype_keys_nb_seqs <- c(v_unique_haplotype_keys_nb_seqs,hap_nb_seq_to_add)
  print(paste0(i, " iterations done out of ",length(v_unique_haplotype_keys),"!"))
  i=i+1
}
df_time_period_Taj_D_with_resamplings$label_hap_key <- paste0(df_time_period_Taj_D_with_resamplings$start_time_period,"/",df_time_period_Taj_D_with_resamplings$stop_time_period,"/",df_time_period_Taj_D_with_resamplings$seq_continent,"/",df_time_period_Taj_D_with_resamplings$Hussin_name)
df_time_period_Taj_D_with_resamplings$nb_sequences <- v_unique_haplotype_keys_nb_seqs[df_time_period_Taj_D_with_resamplings$label_hap_key]#unname(vapply(X = 1:nrow(df_time_period_Taj_D_with_resamplings),FUN = function(i) nrow(subset(df_metadata_samples,((as.Date(date)>=df_time_period_Taj_D_with_resamplings$start_time_period[i])&(as.Date(date)<=df_time_period_Taj_D_with_resamplings$stop_time_period[i])&(seq_continent==df_time_period_Taj_D_with_resamplings$seq_continent[i])&(Hussin_name==df_time_period_Taj_D_with_resamplings$Hussin_name[i])))),FUN.VALUE=c(0)))
df_time_period_Taj_D_with_resamplings$log10_nb_seqs <- log10(df_time_period_Taj_D_with_resamplings$nb_sequences+1)
ggplotRegression(fit = lmp(formula = Tajima_D~nb_sequences,data = df_time_period_Taj_D_with_resamplings,Iter=99999,center = FALSE),ggsave_path = output_workspace,the_filename = paste0("Correlation_",time_period_length,"month_Taj_D_with_resampling_vs_nb_sequences_across_continents.png"),xlabl = paste0("Number of sequenced samples per time period (",time_period_length," month)"),ylabl = "Tajima's D")
ggplotRegression(fit = lmp(formula = Tajima_D~log10_nb_seqs,data = df_time_period_Taj_D_with_resamplings,Iter=99999,center = FALSE),ggsave_path = output_workspace,the_filename = paste0("Correlation_",time_period_length,"month_Taj_D_with_resampling_vs_log10_nb_sequences_across_continents.png"),xlabl = paste0("log10(Number of sequenced samples + 1)\n(time period =",time_period_length," month)"),ylabl = "Tajima's D")

ggplot(data = df_time_period_Taj_D_with_resamplings,mapping = aes(x=factor(start_time_period_formatted,levels=as.character(format(sort(unique(as.Date(start_time_period))),"%B %d, %Y"))),y=log10_nb_seqs,col=as.factor(Hussin_name),group=as.factor(Hussin_name))) + geom_point(stat='summary') + stat_summary(geom="line") + ylab("log10(Number of monthly sequences+1)")+ xlab(paste0("start of the time period\n(1 month)")) + theme_bw() + theme(title =  element_text(size=12),axis.text.x = element_text(angle = 60,hjust = 1,size=6),axis.text.y = element_text(size=10),axis.title.x = element_text(size=12),axis.title.y = element_text(size=12),legend.position="right") + labs(col="Haplotype") + scale_color_manual(values = palette_haplotypes) + facet_wrap(~seq_continent,ncol=1)
ggsave(filename = paste0("Time_series_log10_Nb_new_sequences_HAPLOTYPES_across_country_and_continent_horizontal_comparison.png"), path=output_workspace, width = 9.3, height = 23, units = "cm",dpi = 1200)
ggsave(filename = paste0("Time_series_log10_Nb_new_sequences_across_country_and_continent_horizontal_comparison.svg"), path=output_workspace, width = 9.3, height = 23, units = "cm",dpi = 1200,device=svg)

ggplot(data = df_time_period_Taj_D_with_resamplings,mapping = aes(x=factor(start_time_period_formatted,levels=as.character(format(sort(unique(as.Date(start_time_period))),"%B %d, %Y"))),y=log10_nb_seqs,col=as.factor(Hussin_name),group=as.factor(Hussin_name))) + geom_point(stat='summary') + stat_summary(geom="line") + ylab("log10(Number of monthly sequences+1)")+ xlab(paste0("start of the time period\n(1 month)")) + theme_bw() + theme(title =  element_text(size=12),axis.text.x = element_text(angle = 60,hjust = 1,size=6),axis.text.y = element_text(size=10),axis.title.x = element_text(size=12),axis.title.y = element_text(size=12),legend.position="right") + labs(col="Haplotype") + scale_color_manual(values = palette_haplotypes) + facet_wrap(~seq_continent,nrow=1)
ggsave(filename = paste0("Time_series_log10_Nb_new_sequences_HAPLOTYPES_across_country_and_continent_vertical_comparison.png"), path=output_workspace, width = 35, height = 15, units = "cm",dpi = 1200)
ggsave(filename = paste0("Time_series_log10_Nb_new_sequences_across_country_and_continent_vertical_comparison.svg"), path=output_workspace, width = 35, height = 15, units = "cm",dpi = 1200,device=svg)

ggplot(data = df_time_period_Taj_D_with_resamplings,mapping = aes(x=factor(start_time_period_formatted,levels=as.character(format(sort(unique(as.Date(start_time_period))),"%B %d, %Y"))),y=log10_nb_seqs,col=as.factor(Hussin_name),group=as.factor(Hussin_name))) + geom_point(stat='summary') + stat_summary(geom="line") + ylab("log10(Number of monthly sequences+1)")+ xlab(paste0("start of the time period\n(1 month)")) + theme_bw() + theme(text= element_text(size=14),title =  element_text(size=14),axis.text.x = element_text(angle = 60,hjust = 1,size=10),axis.text.y = element_text(size=12),axis.title.x = element_text(size=14),axis.title.y = element_text(size=14),legend.position="right") + labs(col="Haplotype") + scale_color_manual(values = palette_haplotypes) + facet_grid(Hussin_name~seq_continent)
ggsave(filename = paste0("Grid_Time_series_log10_Nb_new_sequences_HAPLOTYPES_across_country_and_continent_horizontal_comparison.png"), path=output_workspace, width = 35, height = 30, units = "cm",dpi = 1200)
ggsave(filename = paste0("Grid_Time_series_log10_Nb_new_sequences_across_country_and_continent_horizontal_comparison.svg"), path=output_workspace, width = 35, height = 30, units = "cm",dpi = 1200,device=svg)

#function for getting linear regression metrics
regression=function(df){
  #setting the regression function. 
  reg_fun<-lmp(formula=Tajima_D~log10_nb_seqs,data=df,Iter=5000,center = FALSE) #regression function
  #getting the slope, intercept, R square and adjusted R squared of 
  #the regression function (with 3 decimals).
  slope<- unname(round(unname(summary(reg_fun)$coefficients[,1][2]),3))
  R2<- round(as.numeric(broom::glance(reg_fun)$r.squared),3)
  R2.Adj<- round(as.numeric(broom::glance(reg_fun)$adj.r.squared),3)
  pval <- round(as.numeric(broom::glance(reg_fun)$p.value),3)
  c(slope,R2,R2.Adj,pval)
}
#Regressions per haplotype
regressions_data <- ddply(df_time_period_Taj_D_with_resamplings,"Hussin_name",regression)
colnames(regressions_data)<-c ("Hussin_name","slope","R2","R2.Adj","p.value")

qplot(log10_nb_seqs, Tajima_D, data = df_time_period_Taj_D_with_resamplings, size=I(2))+geom_smooth(method="lm")+
  geom_label(data=regressions_data, inherit.aes=FALSE, aes(x = 2.2, y = 2,
                                                           label=paste("slope=",slope,";","R^2=",R2,";","R^2.Adj=",R2.Adj,";","p:",ifelse(p.value==0,"<2e-16",p.value))),label.size = 0.3 ) + theme_bw() + ylab("Tajima's D")+ xlab(paste0("log10(Number of sequences per haplotype per continent per month + 1)")) + theme(plot.title = element_text(size=8))+
  facet_wrap(~Hussin_name,ncol=3)
ggsave(filename = paste0(paste0("Grid_correlation_",time_period_length,"month_Taj_D_with_resampling_vs_log10_nb_sequences_across_haplotypes.png")), path=output_workspace, width = 35, height = 30, units = "cm",dpi = 1200)
ggsave(filename = paste0(paste0("Grid_correlation_",time_period_length,"month_Taj_D_with_resampling_vs_log10_nb_sequences_across_haplotypes.svg")), path=output_workspace, width = 35, height = 30, units = "cm",dpi = 1200,device=svg)

#Regressions per continent
regressions_data <- ddply(df_time_period_Taj_D_with_resamplings,"seq_continent",regression)
colnames(regressions_data)<-c ("seq_continent","slope","R2","R2.Adj","p.value")

qplot(log10_nb_seqs, Tajima_D, data = df_time_period_Taj_D_with_resamplings, size=I(2))+geom_smooth(method="lm")+
  geom_label(data=regressions_data, inherit.aes=FALSE, aes(x = 2.2, y = 2,
                                                           label=paste("slope=",slope,";","R^2=",R2,";","R^2.Adj=",R2.Adj,";","p:",ifelse(p.value==0,"<2e-16",p.value))),label.size = 0.3 ) + theme_bw() + ylab("log10(Number of monthly sequences+1)")+ xlab(paste0("log10(Number of sequenced samples + 1)\n(time period =",time_period_length," month)")) + theme(plot.title = element_text(size=8))+
  facet_wrap(~seq_continent,ncol=2)

ggsave(filename = paste0(paste0("Grid_correlation_",time_period_length,"month_Taj_D_with_resampling_vs_log10_nb_sequences_across_continents.png")), path=output_workspace, width = 25, height = 25, units = "cm",dpi = 1200)
ggsave(filename = paste0(paste0("Grid_correlation_",time_period_length,"month_Taj_D_with_resampling_vs_log10_nb_sequences_across_continents.svg")), path=output_workspace, width = 25, height = 25, units = "cm",dpi = 1200,device=svg)

#Linear mixed models Tajima_D ~ log10_nb_seqs + Hussin_name + seq_continent
df_data_lmm_normalized <- df_time_period_Taj_D_with_resamplings
#boxcox_lamba_log10_nb_seqs <- boxcox(df_time_period_Taj_D_with_resamplings$log10_nb_seqs,optimize = TRUE)$lambda
#df_data_lmm_normalized$log10_nb_seqs <- scale(df_time_period_Taj_D_with_resamplings$log10_nb_seqs,T,T)#boxcoxTransform(df_time_period_Taj_D_with_resamplings$log10_nb_seqs,lambda = boxcox_lamba_log10_nb_seqs)
#boxcox_lamba_Tajima_D<- boxcox(df_time_period_Taj_D_with_resamplings$value,optimize = TRUE)$lambda
df_data_lmm_normalized$Tajima_D <- scale(df_time_period_Taj_D_with_resamplings$Tajima_D,T,T)#boxcoxTransform(df_time_period_Taj_D_with_resamplings$value,lambda = boxcox_lamba_Tajima_D)
M1 <- lmerTest::lmer(Tajima_D~log10_nb_seqs + (1+log10_nb_seqs|Hussin_name) + (1+log10_nb_seqs|seq_continent), data=df_data_lmm_normalized, REML=FALSE)

E1 <- resid(M1)
F1<-fitted(M1)
#Check non-violation of lmm model assumption (residuals normality)
EnvStats::qqPlot(resid(M1), add.line = TRUE,ylab= "Quantiles of regression residuals",main="Normal Q-Q plot for regression residuals")
hist(E1,main = "Regression residuals histogram")
#coefficients
summary(M1)
coef(M1)

#Nested models comparison (detect significant random effects)
#without seq_continent
M2 <- lmerTest::lmer(Tajima_D~log10_nb_seqs + (1+log10_nb_seqs|Hussin_name), data=df_data_lmm_normalized, REML=FALSE)
#without Hussin_name
M3 <- lmerTest::lmer(Tajima_D~log10_nb_seqs + (1+log10_nb_seqs|seq_continent), data=df_data_lmm_normalized, REML=FALSE)

AICc<-c(AICc(M1), AICc(M2), AICc(M3))
# Put values into one table for easy comparision
Model<-c("M1", "M2", "M3")
AICtable<-data.frame(Model=Model, AICc=AICc)
View(AICtable)

anova(M1,M2)
anova(M1,M3)

#Individual explained variance
MuMIn::r.squaredGLMM(M1)
MuMIn::r.squaredGLMM(M2)
MuMIn::r.squaredGLMM(M3)

#Time code execution
end_time <- Sys.time()
paste0("Execution time: ",end_time - start_time,"!")


############################
library("session")
save.session(file = paste0(output_workspace,"Taj_D_top3_Haplotypes_per_continent_",gsub(pattern = ":",replacement = "_",x = gsub(pattern = " ",replacement = "_",x = date())),"_RSession.Rda"))
