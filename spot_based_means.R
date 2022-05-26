#Create an R script to investigate and plot the mean expression with different number of samples (12, 4, 2). This stratification is always
#spot based bcs there is no coherence between samples and between layers. NA, 0 and relative positions

#Load the library and the spe object  
#library(spatialLIBD)
#spe <- fetch_data('spe')

#check the repetition of spots names 
table(colnames(spe))
#Check for unique spots and how many they are 
unique(colnames(spe))
lenght(unique(colnames(spe)))

#Functions ####
##Prepare a common set --> definition of the coordinates of the spots ####

#Access the coordinates of the spots. Since each spot is repeated n times we must choose from which list of values take the coordinates --> arbitrarely set to 1st entry
#It is not a problem if a spot in not present in all samples since the value of expression would be NA 
#The coords function takes in input the spot data frame and gives in output the xx dataframe with the x and y coordinates from the first entry of the indeces of each spot (since each spot is repeated n times. 1 < n < 12)
#The entries of the spot are all accessed and stored as indeces in list_index. This vector is filtered for values in a certain range (i_start and i_stop, defined by the spot data frame passed and coherent with the number and IDs of samples to be used for mean computation), 
#So that we can access the first of the list with no problems of being outside the range of values of the current studied samples 
#If the value is missing bcs no spot (expecialy for the 12 samples --> here we should be looking at whole set so smth present)
coords <- function(spot_frame, i_start, i_stop){
  xx <- data.frame()
  for (i in spot_frame[, 1]){
    name <- i
    list_index_c <- which(rownames(colData(spe)) == name)
    list_index_c <- list_index_c[list_index_c >= i_start & list_index_c <= i_stop][1] #subset the index list that should always have smth since everything derive form the unique subset --> at least 1 should be guaranteed 
    here <- spatialCoords(spe)[list_index_c,]
    xx <- rbind(xx, here)
  }
  names(xx) <- c('x', 'y')
  return(xx)
}

##Prepare a common set --> definition of layers ####

#Access the Layer_Guess_reordered_short for each spot and compute the majority of votes for that spot 
#The function layer_guess takes in input the dataframe of spots and returns in output a list of consensus layer for each spot. 
#If the spot has multiple max voted layers, the first one is picked. If the layer is not assigned (NA) it does not appear in the table  
#The rownames list is subset to have just the values corresponding to the sample/subject/position spots of interest
layer_guess <- function(spot_frame, i_start, i_stop){
  list_consensus <- NULL #Initialize the list of consensus layers found as the value for which votes are max
  for (i in spot_frame[, 1]){
    list_index_l <- NULL #list of indeces of the spots the layer of which is searched
    list_index_l <- which(rownames(colData(spe)) == i)
    list_index_l <- list_index_l[list_index_l >= i_start & list_index_l <= i_stop] #remove the elements out of the boundaries of my searching area and use just those ones 
    values <- colData(spe)$layer_guess_reordered_short[list_index_l] 
    values_data <- as.data.frame(table(values)) #some spots have been assigned as NA so do not appear 
    if (max(values_data$Freq) != 0){
      putative <- which(values_data$Freq == max(values_data$Freq)) #Some spots can have two equally max values --> pick the first one 
      list_consensus <- append(list_consensus, values_data$values[putative[1]] ) 
    } else{ list_consensus <- append(list_consensus, '0' ) #for NA dealing 
    }
  }
  return(list_consensus)
}


##Prepare a common set --> expression computation ####

#compute_expression is a function that takes in input the subset of colnames (range) and the indeces of it (i_start, i_stop) used to look for expression,
#the index of the gene of interest (row_index) used to look in the assays and the number value that indicates how many spots should be present in the mean expression computation (number)
#Access the assays of that gene name indexand compute the mean expression for the spots. The mean expression is retained only if the spots 
#used are n (12, 4, 2) according to the sample of interest. 
#It returns the list of mean expressions for the spots. If < n sptots are used for mean computation the mean is set to NA.

compute_expression <- function(range, i_start, i_stop, row_index, number){
  list_expr <- NULL #define the list of mean expression in every spot that is repeated number (12, 4, 2) times (so present in each sample, subject, position)
  uniques <- unique(range) #take the set of unique spots 
  area <- assays(spe)$logcounts[row_index, i_start:i_stop] #store this to have faster acess 
  length(area)
  for (i in uniques){
    name <- i #define the spot to be searched
    list_index <- which(range == name) #define the indeces at which search for it --> you will consider just the full spots so no NA 
    #list where to store the single expression in single samples --> if this is shorter than number it is transformed in NA. Thrn included in list_expr that will form the column of expression 
    local_expr <- area[list_index]  
    if (length(local_expr) == number){
      local_expr <- mean(local_expr)
    } else { local_expr <- NA
    }
    list_expr <- append(list_expr, local_expr)
  }
  return(list_expr) 
}

##Plot the dataframe####
#plot_it is a function that takes in input the spot_frame of interest and displays the plot.  
#Data are plotted according to the layer and the expression --> color for the single layers (concordant with the paper) an diameter concordant with expression (20*logcounts).
#Dataframe passed in input is filtered according to remove NAs rows (expression) and '0'-layer (NA layer guess)
plot_it <- function(spot_frame, name_matrix, gene){
  spot_frame <- spot_frame[complete.cases(spot_frame), ] #remove the NA rows 
  spot_frame <- spot_frame[spot_frame$layer != '0', ] #remove the rows that have no layer guess 
  colors <- c("#F0027F","#377EB8","#4DAF4A","#984EA3","#FFD700","#FF7F00","#666666")
  symbols(spot_frame$x, spot_frame$y, circles = 20*(spot_frame$expr), fg = colors[as.numeric(spot_frame$layer)], inches = F, xlab = 'x', ylab = 'y', main = paste0(gene , ' expression'), sub = name_matrix)
  #different indexing of legend to deal with no L1 and L2 in Br5595
  if (grepl('Br5595', name_matrix, fixed = T)){
    legend('topright', legend = levels(as.factor(spot_frame$layer)), col = colors[-c(1,2)], cex = 0.5, pch = 0.5) 
  } else {
    legend('topright', legend = levels(as.factor(spot_frame$layer)), col = colors, cex = 0.5, pch = 0.5)
  }
}



#Data ####
#The indexes used are taken from the matrix_dim object of expr_no_layer.R script 
##Prepare a common set --> Define the matrices ####
## Study based - Define the rows of the dataframe to be used as a map #####
i_start_full <- 1
i_stop_full <- ncol(spe)
range <- colnames(spe)
only_spots <- data.frame(unique(range)) #--> 4941 unique spots 
only_spots <- cbind(only_spots, layer_guess(only_spots, i_start_full, i_stop_full), coords(only_spots, i_start_full, i_stop_full))
names(only_spots) <- c('spot', 'layer', 'x', 'y')

## Subject based - Define the matrix to be used as a map ####
### Subject 1 ####
i_start_1_4 <- 1
i_stop_1_4 <- 18033
range_1 <- colnames(spe[, 1:18033])
only_spots_1 <- data.frame(unique(colnames(spe[, 1:18033]))) #--> 4902 unique spots 
only_spots_1 <- cbind(only_spots_1, layer_guess(only_spots_1, i_start_1_4, i_stop_1_4), coords(only_spots_1, i_start_1_4, i_stop_1_4))
names(only_spots_1) <- c('spot', 'layer', 'x', 'y')
### Subject 2 ####
i_start_2_4 <- 18034
i_stop_2_4 <- 18033 +15284
range_2 <- colnames(spe[, 18034:(18033 +15284)])
only_spots_2 <- data.frame(unique(colnames(spe[, 18034:(18033 +15284)]))) #--> 4481 unique spots 
only_spots_2 <- cbind(only_spots_2, layer_guess(only_spots_2, i_start_2_4, i_stop_2_4), coords(only_spots_2, i_start_2_4, i_stop_2_4))
names(only_spots_2) <- c('spot', 'layer', 'x', 'y')
### Subject 3 ####
i_start_3_4 <- 18033 +15284 + 1
i_stop_3_4 <- 18033 + 15284 + 14364
range_3 <- colnames(spe[ , (18033 + 15285):(18033 + 15284 + 14364)])
only_spots_3 <- data.frame(unique(colnames(spe[, (18033 + 15285):(18033 + 15284 + 14364)]))) #--> 4186 unique spots 
only_spots_3 <- cbind(only_spots_3, layer_guess(only_spots_3, i_start_3_4, i_stop_3_4 ), coords(only_spots_3, i_start_3_4, i_stop_3_4))
names(only_spots_3) <- c('spot', 'layer', 'x', 'y')

##Position based - Define the matrix to be used as a map ####
### Subject 1 ####
#### Position 0 ####
i_start_1_0 <- 1
i_stop_1_0 <- (4226 + 4384)
range_1_0 <- colnames(spe[, 1:(4226 + 4384)])
only_spots_1_0 <- data.frame(unique(colnames(spe[, 1: (4226 + 4384)]))) #--> 4428 unique spots 
only_spots_1_0 <- cbind(only_spots_1_0, layer_guess(only_spots_1_0, i_start_1_0, i_stop_1_0), coords(only_spots_1_0, i_start_1_0, i_stop_1_0))
names(only_spots_1_0) <- c('spot', 'layer', 'x', 'y')
#### Position 300 ####
i_start_1_300 <- 4226 + 4384 + 1
i_stop_1_300 <- 4226 + 4384 + 4789 + 4634
range_1_300 <- colnames(spe[, (4226 + 4384 + 1):(4226 + 4384 + 4789 + 4634)])
only_spots_1_300 <- data.frame(unique(colnames(spe[, (4226 + 4384 + 1):(4226 + 4384 + 4789 + 4634)]))) #--> 4803 unique spots
only_spots_1_300 <- cbind(only_spots_1_300, layer_guess(only_spots_1_300, i_start_1_300, i_stop_1_300), coords(only_spots_1_300, i_start_1_300, i_stop_1_300))
names(only_spots_1_300) <- c('spot', 'layer', 'x', 'y')
### Subject 2 ####
####Position 0 ####
i_start_2_0 <- 18033 + 1
i_stop_2_0 <- 18033 + 3661 + 3498
range_2_0 <- colnames(spe[, (18033 + 1) : (18033 + 3661 + 3498)])
only_spots_2_0 <- data.frame(unique(colnames(spe[, (18033 + 1) : (18033 + 3661 + 3498)]))) #--> 3796 unique spots 
only_spots_2_0 <- cbind(only_spots_2_0, layer_guess(only_spots_2_0, i_start_2_0, i_stop_2_0), coords(only_spots_2_0, i_start_2_0, i_stop_2_0) )
names(only_spots_2_0) <- c('spot', 'layer', 'x', 'y')
#### Position 300 ####
i_start_2_300 <- 18033 + 1 + 3661 + 3498
i_stop_2_300 <- 18033 + 3661 + 3498 + 4110 + 4015
range_2_300 <- colnames(spe[, (18033 + 1 + 3661 + 3498 ) : (18033 + 3661 + 3498 + 4110 + 4015)])
only_spots_2_300 <- data.frame(unique(colnames(spe[, (18033 + 1 + 3661 + 3498 ) : (18033 + 3661 + 3498 + 4110 + 4015)]))) #--> 4427 unique spots 
only_spots_2_300 <- cbind(only_spots_2_300, layer_guess(only_spots_2_300, i_start_2_300, i_stop_2_300), coords(only_spots_2_300, i_start_2_300, i_stop_2_300))
names(only_spots_2_300) <- c('spot', 'layer', 'x', 'y')
### Subject 3 ####
#### Position 0 ####
i_start_3_0 <- 33317 + 1
i_stop_3_0 <- 33317 + 3639 + 3673
range_3_0 <- colnames(spe[, (33317 + 1) : (33317 + 3639 + 3673)])
only_spots_3_0 <- data.frame(unique(colnames(spe[, (33317 + 1) : (33317 + 3639 + 3673)]))) #--> 3810 unique spots 
only_spots_3_0 <- cbind(only_spots_3_0, layer_guess(only_spots_3_0, i_start_3_0, i_stop_3_0), coords(only_spots_3_0,i_start_3_0, i_stop_3_0))
names(only_spots_3_0) <- c('spot', 'layer', 'x', 'y')
#### Position 300 ####
i_start_3_300 <- 33317 + 1 + 3639 + 3673
i_stop_3_300 <- 33317 + 3639 + 3673 + 3592 + 3460
range_3_300 <- colnames(spe[, (33317 + 1 + 3639 + 3673) : (33317 + 3639 + 3673 + 3592 + 3460)])
only_spots_3_300 <- data.frame(unique(colnames(spe[, (33317 + 1 + 3639 + 3673) : (33317 + 3639 + 3673 + 3592 + 3460)]))) #--> 3994 unique spots 
only_spots_3_300 <- cbind(only_spots_3_300, layer_guess(only_spots_3_300, i_start_3_300, i_stop_3_300), coords(only_spots_3_300, i_start_3_300, i_stop_3_300))
names(only_spots_3_300) <- c('spot', 'layer', 'x', 'y')


#Gene expression ####
gene <- 'SNAP25'
row_index <- which(rowData(spe)$'gene_name' == gene)

###Study-based - Fill the matrix cells with the mean expression of the gene for that spot across all samples --> mean of 12 points  ####
samples_12 <- cbind(only_spots, compute_expression(range,i_start_full, i_stop_full,row_index, 12))
names(samples_12)[5] <- 'expr'
plot_it(samples_12, 'Mean_all',gene)

### Subject-based - Fill the matrix cells with the mean expression of the gene for that spot across all samples --> mean of 4 points  ####
samples_1_4 <- cbind(only_spots_1, compute_expression(range_1,i_start_1_4,i_stop_1_4,row_index, 4))
names(samples_1_4)[5] <- 'expr'
plot_it(samples_1_4, 'Mean_subj_Br5292', gene)
###Subject-based - Fill the matrix cells with the mean expression of the gene for that spot across all samples --> mean of 4 points  ####
samples_2_4 <- cbind(only_spots_2, compute_expression(range_2,i_start_2_4,i_stop_2_4, row_index, 4))
names(samples_2_4)[5] <- 'expr'
plot_it(samples_2_4, 'Mean_subj_Br5595', gene)
###Subject-based - Fill the matrix cells with the mean expression of the gene for that spot across all samples --> mean of 4 points  ####
samples_3_4 <- cbind(only_spots_3, compute_expression(range_3,i_start_3_4,i_stop_3_4,row_index, 4))
names(samples_3_4)[5] <- 'expr'
plot_it(samples_3_4, 'Mean_subj_Br8100', gene)
#At the end you should have three matrices unique_rows x 2, one for subject. 
#Each cell should be the mean of 4 observation of expression 

###Position-based - Fill the matrix cells with the mean expression of the gene for that spot across all samples --> mean of 2 points  ####
samples_1_0 <- cbind(only_spots_1_0, compute_expression(range_1_0,i_start_1_0,i_stop_1_0,row_index, 2))
names(samples_1_0)[5] <- 'expr'
plot_it(samples_1_0, 'Mean_subj_Br5292_0', gene)
###Position-based -Fill the matrix cells with the mean expression of the gene for that spot across all samples --> mean of 2 points  ####
samples_1_300 <- cbind(only_spots_1_300, compute_expression(range_1_300, i_start_1_300, i_stop_1_300,row_index, 2))
names(samples_1_300)[5] <- 'expr'
plot_it(samples_1_300, 'Mean_subj_Br5292_300', gene)
###Position-based - Fill the matrix cells with the mean expression of the gene for that spot across all samples --> mean of 2 points  ####
samples_2_0 <- cbind(only_spots_2_0, compute_expression(range_2_0,i_start_2_0,i_stop_2_0,row_index, 2))
names(samples_2_0)[5] <- 'expr'
plot_it(samples_2_0, 'Mean_subj_Br5595_0', gene)
###Position-based - Fill the matrix cells with the mean expression of the gene for that spot across all samples --> mean of 4 points  ####
samples_2_300 <- cbind(only_spots_2_300, compute_expression(range_2_300,i_start_2_300,i_stop_2_300,row_index, 2))
names(samples_2_300)[5] <- 'expr'
plot_it(samples_2_300, 'Mean_subj_Br5595_300', gene)
###Position-based - Fill the matrix cells with the mean expression of the gene for that spot across all samples --> mean of 2 points  ####
samples_3_0 <- cbind(only_spots_3_0, compute_expression(range_3_0,i_start_3_0,i_stop_3_0,row_index, 2))
names(samples_3_0)[5] <- 'expr'
plot_it(samples_3_0, 'Mean_subj_Br100_0', gene)
###Position-based - Fill the matrix cells with the mean expression of the gene for that spot across all samples --> mean of 2 points  ####
samples_3_300 <- cbind(only_spots_3_300, compute_expression(range_3_300,i_start_3_300,i_stop_3_300,row_index, 2))
names(samples_3_300)[5] <- 'expr'
plot_it(samples_3_300, 'Mean_subj_Br8100_300', gene)


























##Prepare a common frame --> extract layer-spots ####
#layer_spot function works on modified colData(spe) and extract all the spots at layer level 
#modified --> remove of NA spots --> spot_frame dataset
#Function --> from the entire data set it returns a smaller with the spot name, the index in spe, the spot key, the subject, the position and the expression for the given gene  

remove_spots <- which(is.na(colData(spe)$layer_guess_reordered_short)) #remove NA 
spot_frame <- colData(spe)[-remove_spots, ] 

layer_spot <- function(spot_frame, layer_name, row_index){
  area <- assays(spe)$logcounts[row_index, ] #define the assay to access fastly 
  kept_spots <- spot_frame[spot_frame$layer_guess_reordered_short == layer_name, ] #keep the spots assigned to layer 2 
  #table(kept_spots) #all unique for spot name
  #length(unique(rownames(kept_spots))) #the number of unique spots 
  expression_list <- currents <- NULL #list where to store the expressions 
  for(i in 1:nrow(kept_spots)){
    current_index <- which(colData(spe)$key == kept_spots$key[i]) #find the index of the spot in the general --> so find the expression of the spot given the patient and the layer 
    currents <- append(currents, current_index)
    expression_list <-append(expression_list, area[current_index])
  }
  dataframe_kept <- data.frame('spots' = rownames(kept_spots), 'index' = currents, 'key' = kept_spots$key, 'sub' = kept_spots$subject, 'sub_pos' = kept_spots$position,  'expr' = expression_list)
  return(dataframe_kept)
}

#Why not L1 in subj 2? --> probelm with NA factors --> they are not reported in table so the first max value is corresponds to L1

#Compute the means at all levels you want --> these are means of all spots so spatial/coordinate is lost 
compute_means <-- function(spot_frame){
  means_i_2_small <- NULL #store here the pos means --> vector long 6 (2 x 3)
  subject <- c("Br5292","Br5595","Br8100", '0', '300')
  names <- c('L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'WM')
  means <- data.frame('Layer' = c('L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'WM'))
  means_i_4 <- matrix(NA, 7, 3)
  means_i_2 <- matrix(NA, 7, 6)
  means_12 <- NULL
  for (i in 1:length(names)){
    spots_L <- layer_spot(spot_frame, layer_name = names[i], row_index = row_index)
    #length(unique(dataframe_kept$spots)) --> 1740 
    means_12 <- append(means_12, mean(spots_L$expr))
    for (subj in 1:3){
      if (subject[subj] %in% unique(spots_L$sub)){
        means_i_4[i, subj] <- mean(subset(spots_L$expr, spots_L$sub == subject[subj] ))
        for (pos in 1:2){
          means_i_2_small<-append(means_i_2_small, mean(subset(spots_L$expr, spots_L$sub == subject[subj] & spots_L$sub_pos == subject[(3+pos)])))
        }
      }else { means_i_4[i, subj] <- NA
      means_i_2_small <- append(means_i_2_small, c(NA, NA))
      }
    }
    means_i_2[i, ] <- means_i_2_snall
  }
  means <- cbind(means, means_12, means_i_4, means_i_2)
  names(means) <- c('all_samples', '1_4', '2_4', '3_4', '1_2_0', '1_2_300', '2_2_0', '2_2_300', '3_2_0', '3_2_300')
  return(means)
}




##STuuf to be discussed ####
#Coordinates to be used --> of which subject?
#The problem with assignment of layers --> what to do at the margins? --> 
#Example with 12 samples mean computation --> report the max_voted spot layer 
#[1] 4 6 
#[1] 3
#[1] 6 7
#[1] 3 4 5 6 7
#[1] 6
#So the computation is precise at spot level --> that spot experienced that measure 
#But the layer is not always the same --> maybe the tissue positioned differently or just different 
#How to deal? remove the spots that have not the concordance? Order by layer and see if for that spot at tat layer u have n obs? --> then how many obs u have?
#sparse matrices 