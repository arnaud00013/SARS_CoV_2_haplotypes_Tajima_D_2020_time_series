#@Author=Arnaud NG
#This script compute Tajima's D using SARS-CoV-2 consensus sequences of the 3 most widespread haplotypes in a continent (NCBI)

#Time code execution
start_time <- Sys.time()

#import libraries
library("ggplot2")
library("seqinr")
library("genbankr")
library("grid")
library("RColorBrewer")
library("randomcoloR")
library("gplots")
library("lmPerm")
library("gridExtra")
library("RColorBrewer")
library("parallel")
library("foreach")
library("doParallel")
library("session")
#import script arguments
#ABSOLUTE Path of the folder containing the data 
output_workspace <- as.character(commandArgs(TRUE)[1]) 
#Number of cpus for Tajima's D analysis. It corresponds to the number of month that will be analyzed in parallel. Thus, it should be <=5 for Wave 2.
nb_cores <- as.integer(commandArgs(TRUE)[2]) 
nb_cores <- ifelse(test = nb_cores>5,yes=5,no=nb_cores)
#The continent analyzed
the_continent <- as.character(commandArgs(TRUE)[3]) 
the_continent <- ifelse(test = grepl(pattern = "North",x = the_continent,fixed = T),yes="North America",no=the_continent)
the_continent <- ifelse(test = grepl(pattern = "South",x = the_continent,fixed = T),yes="South America",no=the_continent)

#Set language as English for date formatting
Sys.setlocale("LC_ALL","English")

#create dataframe with samples consensus seq 
df_fasta_consensus_seq <- read.csv2(file = paste0(output_workspace,"Hussingroup_Inter_db_consensus_sequences.fasta"),sep = ",",header = F,stringsAsFactors = FALSE)
df_consensus_seq <- data.frame(Sample=unname(vapply(X = df_fasta_consensus_seq[which((1:nrow(df_fasta_consensus_seq))%%2==1),1],FUN= function(x) substr(x,2,nchar(x)),FUN.VALUE = c(""))),Consensus_seq=df_fasta_consensus_seq[which((1:nrow(df_fasta_consensus_seq))%%2==0),1],stringsAsFactors = F)
rownames(df_consensus_seq) <- df_consensus_seq$Sample

#import metadata 
df_metadata_samples <- read.csv2(file = paste0(output_workspace,"HussinGroup_Interdb_metadata.csv"),sep = ",",header = TRUE,stringsAsFactors = FALSE)
df_metadata_samples$X <- NULL
rownames(df_metadata_samples) <- df_metadata_samples$GISAID_id
df_metadata_samples <- subset(df_metadata_samples,(!is.na(date))&(!is.na(seq_country)))
#focus on human hosted viral sequences without too much ambiguous calls
df_metadata_samples <- subset(df_metadata_samples,(host_species=="Human")&(seqdata=="yes"))
#find top 3 most widespread haplotype in current continent
if ((length(unique(subset(df_metadata_samples,seq_continent==the_continent)$Hussin_name))>=3)){
  the_top3_haplotypes <- names(sort(table(subset(df_metadata_samples,seq_continent==the_continent)$Hussin_name),dec=T))[1:3]
}else{
  the_top3_haplotypes <- names(sort(table(subset(df_metadata_samples,seq_continent==the_continent)$Hussin_name),dec=T))
}
#focus on second wave
df_metadata_samples <- subset(df_metadata_samples,(as.Date(date)>=as.Date("2020-08-01"))&((as.Date(date)<=as.Date("2020-12-31"))))

df_consensus_seq <- df_consensus_seq[rownames(df_metadata_samples),]

#focus on a certain continent
df_metadata_samples <- subset(df_metadata_samples,seq_continent==the_continent)


v_seq_country_to_seq_continent <- unique(df_metadata_samples[,c("seq_country","seq_continent")])[,"seq_continent"]
names(v_seq_country_to_seq_continent) <- unique(df_metadata_samples[,c("seq_country","seq_continent")])[,"seq_country"]

v_continent <- unique(df_metadata_samples$seq_continent)
v_continent <- v_continent[!is.na(v_continent)]
v_continent <- v_continent[v_continent!="NA"]

library(lubridate)
time_period_length <- 1 #number of month
nb_time_periods <- ceiling(as.numeric(as.Date("2020-12-31")-as.Date("2020-08-01"))/31)
df_time_periods <- data.frame(time_period=1:nb_time_periods,start=as.Date("2020-08-01") %m+% months(time_period_length*(0:(nb_time_periods-1))),stop=(as.Date("2020-08-01")+(time_period_length-1))%m+% months(time_period_length*(1:(nb_time_periods)))-1,Tajima_D=NA,stringsAsFactors = FALSE)
df_time_periods$middle_day <- as.Date(df_time_periods$start)+((as.Date(df_time_periods$stop)-as.Date(df_time_periods$start))/2)
detach("package:lubridate", unload=TRUE)
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

#create a function that does all possible pariwise COMPARISONS between a set of unique sequences 
#THE FUNCTION TAKES INTO ACCOUNT THE PRESENCE OF AMBIGUOUS CONSENSUS CALLS ("N")
#Execution time is O(n^2)
calculate_nb_pwd <- function(v_seqs){
  output <- 0 #initialize the variable that will contain the number of pairwise differences in single nucleotide sites
  pwdiff_with_others<- function(the_sequence,lst_unique_sequences){
    if (any(duplicated(lst_unique_sequences))){
      stop("There are doublons in the list of unique sequences.  Retry!")
    }
    result<-rep(NA,length(lst_unique_sequences))
    ind_res_list<-1
    for (current_diff_sequence in lst_unique_sequences){
      current_result<-0
      for (i in 1:nchar(v_seqs[1])){
        if ((substr(x = the_sequence,i,i)!=substr(x = current_diff_sequence,i,i))&(substr(x = the_sequence,i,i)!="N")&("N"!=substr(x = current_diff_sequence,i,i))){
          current_result <- current_result +1
        }
      }
      result[ind_res_list] <- current_result
      ind_res_list <- ind_res_list + 1
    }
    return(result)
  }
  v_unique_freqs <- list(table(v_seqs))[[1]]
  
  #comparisons are not duplicated here! 
  for (current_unique_sequence in names(v_unique_freqs)){
    output <- output + sum(v_unique_freqs[current_unique_sequence]*pwdiff_with_others(the_sequence = current_unique_sequence,lst_unique_sequences = names(v_unique_freqs)[which(names(v_unique_freqs)==current_unique_sequence):length(names(v_unique_freqs))]))
  }
  return(output)
}  
#Tajima's D with resampling of samples from a certain time period
get_time_period_Tajima_D_with_resampling <- function(the_time_period,n,k,the_seq_continent,the_current_lineage){
  v_consensus_seq_samples_sequenced_at_this_time_period <- df_consensus_seq[intersect(df_metadata_samples[(as.Date(df_metadata_samples$date)>=df_time_periods$start[the_time_period])&(as.Date(df_metadata_samples$date)<=df_time_periods$stop[the_time_period])&(df_metadata_samples$seq_continent==the_seq_continent)&(df_metadata_samples$Hussin_name==the_current_lineage),"GISAID_id"],rownames(df_consensus_seq)),"Consensus_seq"]
  if ((length(v_consensus_seq_samples_sequenced_at_this_time_period)==0)){
    return(data.frame(time_period=the_time_period,Tajima_D=NA,start_time_period=as.Date(df_time_periods[df_time_periods$time_period==the_time_period,"start"]),num_resampling=NA,seq_continent=the_seq_continent,Hussin_name=the_current_lineage,stringsAsFactors=FALSE))
  }else if (n>=length(v_consensus_seq_samples_sequenced_at_this_time_period)){
    #print(length(v_consensus_seq_samples_sequenced_at_this_time_period))
    warning("Sample size for the period ",the_time_period," is too small. Thus, we assigned it a missing Tajima's D value!")
    return(data.frame(time_period=the_time_period,Tajima_D=NA,start_time_period=as.Date(df_time_periods[df_time_periods$time_period==the_time_period,"start"]),num_resampling=NA,seq_continent=the_seq_continent,Hussin_name=the_current_lineage,stringsAsFactors=FALSE))
  }else if ((n<=0)||(k<=0)){
    stop("Re-sampling size and number of resamplings should be positive integers!")
  }else{
    df_downsamples_Taj_D <-  NULL
    for (i in 1:k){
      v_current_resampling_consensus_seq <- sample(x = v_consensus_seq_samples_sequenced_at_this_time_period,size = n,replace = FALSE)
      if (is.na(v_current_resampling_consensus_seq[1])){
        print(the_date)
      }
      nb_copy_genome <- n
      #calculate number of pairwise differences without considering "N" (ambiguous consensus call)
      nb_pwdiffs_current_resampling <- calculate_nb_pwd(v_current_resampling_consensus_seq)
      #calculate number of segregating sites without considering "N" (ambiguous consensus call)
      nb_segreg_sites_current_resampling <- sum(vapply(X = 1:(max(nchar(df_consensus_seq$Consensus_seq),na.rm=T)),FUN = function(the_indx) return(length(table(substr(x = v_current_resampling_consensus_seq,start = the_indx,stop = the_indx))[names(table(substr(x = v_current_resampling_consensus_seq,start = the_indx,stop = the_indx)))!="N"])>1),FUN.VALUE = TRUE))
      #parameters to calculate expected sqrt_variance
      a1_current <- sum((1:(nb_copy_genome-1))^-1)
      a2_current <- sum((1:(nb_copy_genome-1))^-2)
      b1_current <- (nb_copy_genome+1)/(3*(nb_copy_genome-1))
      b2_current <- (2*((nb_copy_genome^2)+nb_copy_genome+3))/((9*nb_copy_genome)*(nb_copy_genome-1))
      c1_current <- b1_current - (1/a1_current)
      c2_current <- b2_current - ((nb_copy_genome+2)/(a1_current*nb_copy_genome)) + (a2_current/(a1_current^2))
      e1_current <- c1_current/a1_current
      e2_current <- c2_current/((a1_current^2)+a2_current)
      
      #Find S , find a1_current, calculate expected sqrt_variance with formula form Tajima (1989) and calculate Tajima's D and Ne_S_Taj according to it
      sqrt_expected_variance_current <- sqrt((e1_current*nb_segreg_sites_current_resampling)+((e2_current*nb_segreg_sites_current_resampling)*(nb_segreg_sites_current_resampling-1)))
      
      #Thetas and Tajima's D
      Theta_pi_current <- nb_pwdiffs_current_resampling/(choose(k = 2,n = nb_copy_genome))
      Theta_w_current <- nb_segreg_sites_current_resampling/(a1_current)
      df_downsamples_Taj_D <- rbind(df_downsamples_Taj_D,data.frame(time_period=the_time_period,Tajima_D=(((Theta_pi_current - Theta_w_current))/sqrt_expected_variance_current),start_time_period=as.Date(df_time_periods[df_time_periods$time_period==the_time_period,"start"]),num_resampling=i,seq_continent=the_seq_continent,Hussin_name=the_current_lineage,stringsAsFactors=FALSE))
      print(paste0(i," iterations out of ",k," done for period ",the_time_period,"!"))
    }
  }
  print(paste0("Tajima's D analysis with resampling is done for period #",the_time_period))
  print("*****************************************************************************")
  return(df_downsamples_Taj_D)
}

#Resamplings' Tajima's D
v_time_periods <- sort(unique(df_time_periods$time_period))
nb_time_periods <- length(v_time_periods)
lst_splits <- split(1:nb_time_periods, ceiling(seq_along(1:nb_time_periods)/(nb_time_periods/nb_cores)))
the_f_parallel <- function(i_cl){
  the_vec<- lst_splits[[i_cl]]
  i_time_period <- 1
  current_time_period_df_Taj_D_from_all_resamplings <- NULL
  for (current_time_period in v_time_periods[the_vec]){
    for (current_haplotype in the_top3_haplotypes){
      current_time_period_df_Taj_D_from_all_resamplings <- rbind(current_time_period_df_Taj_D_from_all_resamplings,get_time_period_Tajima_D_with_resampling(the_time_period = current_time_period, n = 20,k = 200,the_seq_continent=the_continent,the_current_lineage = current_haplotype))
      print(paste0("Core ",i_cl,";ITERATION ",i_time_period,"; Analysis done for haplotype ",current_haplotype," in ",the_continent,"!"))
    }
    print(paste0("Core ",i_cl,": ",i_time_period," iterations done out of ",length(v_time_periods[the_vec]),"!"))
    i_time_period <- i_time_period + 1
  }
  return(current_time_period_df_Taj_D_from_all_resamplings)
}
cl <- makeCluster(nb_cores,outfile=paste0(output_workspace,"LOG_second_wave_",the_continent,"_top3_haplotypes_resamplings_Tajima_D.txt"))
registerDoParallel(cl)
df_time_period_Taj_D_with_resamplings <- foreach(i_cl = 1:nb_cores, .combine = rbind, .packages=c("ggplot2","seqinr","grid","RColorBrewer","randomcoloR","gplots","RColorBrewer","tidyr","infotheo","parallel","foreach","doParallel"))  %dopar% the_f_parallel(i_cl)
stopCluster(cl)

df_time_period_Taj_D_with_resamplings$start_time_period_formatted <- format(df_time_period_Taj_D_with_resamplings$start_time_period,"%B %d")
#save table 
write.table(x = df_time_period_Taj_D_with_resamplings,file = paste0(output_workspace,"Table_time_series_Taj_D_with_resamplings_second_wave_",the_continent,"_top3_haplotypes.csv"),row.names = F,col.names = T,sep=",")

#Time code execution
end_time <- Sys.time()
paste0("Execution time: ",end_time - start_time,"!")
############################
library("session")
save.session(file = paste0(output_workspace,"Taj_D_second_wave_",the_continent,"_top3_haplotypes_RSession.Rda"))
