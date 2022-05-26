#Create an R script to create a set of matrices that store the information regarding clustering per rows (same y) to be used for computations based 
#on spots distances. The script for each sample can plot the spots colored according to 'columns' of belonging (plot) and 
#the variation of the expression according to the flowing along rows and layers (interaction.plot). 
#The script has a chapter for k-mean approach test for the clustering of layers based on the gene expression  

#Load the library
library(spatialLIBD)
spe <- fetch_data('spe')

#Create the matrix at spot level retrieving information about each single spot. 
#Data frame with 4781 rows and 8 columns corresponding to 
#the spot name, the layer of belonging, the position inside the patient, the patient name, the key for the univocity of the spot
#the x coord and the y coord (Ho preso x come row pixel e y come column pixel)


spots_frame <- data.frame('Spot' = rownames(colData(spe)), 'Layer' = spe$layer_guess_reordered_short, 'Pos' = spe$position,
                          'j_subj' = spe$subject, 'Cells' = spe$cell_count, 'ID' = spe$key,
                          'x' = spatialCoords(spe)[, 1], 'y' = spatialCoords(spe)[,2]) # --> this dataframe has (352,9) sub-dataframe of NA in Layer column --> will be removed after the assignment of single columns for single spots  

#If we get for each spot the coordinates and then we order them according to the y coord so that we have spots 
#separated per columns  --> correspondance : 
#table(spatialCoords(spe)[as.integer(rownames(spots_frame)[spots_frame$j_subj == 'Br5292']), ] == spatialCoords(spe)[1:18033, ]) --> TRUE 36066 



#apply_class is a function that allows for classify the spots from the single replicates. It takes the general spots_frame dataframe and for each 
#sample_id (12 replicates) it creates a sub_datafarme over which work. This is subset by the creat_columns function such that we can have the labelling of each spot 
#according to the replicate of origin --> in the general spots_frame we discriminate class 1 of 151507 from class 1 of 151508 by means of the key of the spot
#The function returns the spots_frame with additional column. Since the spots_fram already has the gene expression we can associate the gene expression to the single spots and spots' column 
apply_class <- function(spots_frame, gene){
  colors <- c("#F0027F","#377EB8","#4DAF4A","#984EA3","#FFD700","#FF7F00","#666666")
  #Creat_columns is a function that takes the dataframe and the number parameters. The dataframe is divided according to the number parameter and 
  #a number of classes for the set of spots is created in order to have the way to discriminate between spots in different columns (y-based)
  #The function returns the list of single class-assignements of the spots. The 'number' parameter is assigned by personal decision based on 
  #the table of distances. The striking point of distance has been taken but some data are 'dirty'. For each case this is specified in alongside comments 
  creat_columns <- function(frame, number){
    #Access the spots starting from the first spot and checking if the following spot wrt the current is further than number. 
    #If yes report the following spot that is the changing point. If the first point is self standing than we see it bcs there is 2 as the first change.
    #The cycle goes until the second to last so the last, if different, is reported otherwise is included in the previous column 
    list_i <- NULL 
    for(i in 1:(nrow(frame)-1)){
      if ((frame$y[(i+1)] - frame$y[i]) > number){
        list_i <- append(list_i, i)} 
    } #create the indeces for the columns 
    #Create the list for the storage of the length of each class --> list_index
    list_index <- append((list_index = NULL), list_i[1])
    for(i in 2:length(list_i)){
      list_index[i] <- (list_i[i] - list_i[(i-1)])
    }
    #Create the sequence of labels for each column and the list where to store the label for each spot --> class_found and class_assign 
    class_found <- seq(1, length(list_i), length = length(list_i))
    class_assign <- NULL
    for(i in 1 : length(list_index)){
      class_assign <- append(class_assign, rep(class_found[i], list_index[i])) 
    }
    #You need to add the final set of spots for the final column --> you have gone through all spaces until the last changing point
    #But you are not sure that this last point is the last spot too --> it might be a last space to fill --> check for the difference in spots with the original frame and assign last observations 
    last_class <- rep((tail(class_assign, n = 1) + 1), (nrow(frame) - length(class_assign)))
    class_assign <- append(class_assign, last_class)
    return(class_assign)
  }
  
  
  #remeber the IDs and the number of spots for each
  #151507 151508 151509 151510 151669 151670 151671 151672 151673 151674 151675 151676 
  #4226   4384   4789   4634   3661   3498   4110   4015   3639   3673   3592   3460 
  
  ## Subject 1 ####
  ### Pos 0 ####
  Br_1_1_0 <- spots_frame[grepl('151507', as.character(spots_frame$ID), fixed = TRUE), ] #4226 8 
  Br_1_1_0 <- Br_1_1_0[order(Br_1_1_0$y, Br_1_1_0$x), ] #ordering the dataframe for the value of y and then x 
  #table(Br_1_1_0[2: nrow(Br_1_1_0), ]$x - Br_1_1_0$x) #check the difference btw the columns/y values to check for the columns
  #It seems that for this sample we have 93.99 as the striking point of change --> divide for this value 
  #basically you consider a pot to be in the same column if it is less than 93 pxl from the current spots 
  #We will divide the dataframe every time we observe a difference of 93.
  Br_1_1_0$'Column' <- creat_columns(Br_1_1_0, 93) # 4226 9 
  #plot(Br_1_1_0$x, Br_1_1_0$y, col = rainbow(78)[Br_1_1_0$Column], xlab = 'x', ylab = 'y', main = 'Br5292 pos 0 rep 1') #plot the points colored for the column
  #interaction.plot(Br_1_1_0$Column, Br_1_1_0$Layer, response = Br_1_1_0$Expr, fun = mean, main = paste(gene, 'mean Expression going through the x-axis', sep = '_'),
  #                 col = colors, lty = 1, lwd = 1, legend = T, xlab = 'X-column', ylab = 'Mean Expression', xpd = F,
  #                 trace.label = 'Br5292_Pos_0_Rep_1', fixed = TRUE)
  
  Br_1_2_0 <- spots_frame[grepl('151508', as.character(spots_frame$ID), fixed = TRUE), ] #4384 8
  Br_1_2_0 <- Br_1_2_0[order(Br_1_2_0$y, Br_1_2_0$x), ]
  #table(Br_1_2_0[2: nrow(Br_1_2_0), ]$x - Br_1_2_0$x) # 92 
  Br_1_2_0$'Column' <- creat_columns(Br_1_2_0, 92)
  #plot(Br_1_2_0$x, Br_1_2_0$y, col = rainbow(78)[Br_1_2_0$Column], xlab = 'x', ylab = 'y', main = 'Br5292 pos 0 rep 2')
  #interaction.plot(Br_1_2_0$Column, Br_1_2_0$Layer, response = Br_1_2_0$Expr, fun = mean, main = paste(gene, 'mean Expression going through the x-axis', sep = '_'),
  #                col = colors, lty = 1, lwd = 1, legend = T, xlab = 'X-column', ylab = 'Mean Expression', xpd = F,
  #                 trace.label = 'Br5292_Pos_0_Rep_2', fixed = TRUE)
 
  
  ### Pos 300 ####
  Br_1_1_300 <- spots_frame[grepl('151509', as.character(spots_frame$ID), fixed = TRUE), ] # 4789 8
  Br_1_1_300 <- Br_1_1_300[order(Br_1_1_300$y, Br_1_1_300$x), ] #order
  #table(Br_1_1_300[2:nrow(Br_1_1_300), ]$x - Br_1_1_300$x) # 115
  Br_1_1_300$'Column' <- creat_columns(Br_1_1_300, 115)
  #plot(Br_1_1_300$x, Br_1_1_300$y, col = rainbow(78)[Br_1_1_300$Column], xlab = 'x', ylab = 'y', main = 'Br5292 pos 300 rep 1')
  #interaction.plot(Br_1_1_300$Column, Br_1_1_300$Layer, response = Br_1_1_300$Expr, fun = mean, main = paste(gene, 'mean Expression going through the x-axis', sep = '_'),
  #                 col = colors, lty = 1, lwd = 1, legend = T, xlab = 'X-column', ylab = 'Mean Expression', xpd = F,
  #                 trace.label = 'Br5292_Pos_300_Rep_1', fixed = TRUE)
  
  Br_1_2_300 <- spots_frame[grepl('151510', as.character(spots_frame$ID), fixed = TRUE), ] #4634 8
  Br_1_2_300 <- Br_1_2_300[order(Br_1_2_300$y, Br_1_2_300$x), ]
  #table(Br_1_2_300[2:nrow(Br_1_2_300), ]$x - Br_1_2_300$x)#53 but there are low values that stain this value
  Br_1_2_300$'Column' <- creat_columns(Br_1_2_300, 53)
  #plot(Br_1_2_300$x, Br_1_2_300$y, col = rainbow(78)[Br_1_2_300$Column], xlab = 'x', ylab = 'y', main = 'Br5292 pos 300 rep 2')
  #interaction.plot(Br_1_2_300$Column, Br_1_2_300$Layer, response = Br_1_2_300$Expr, fun = mean, main = paste(gene, 'mean Expression going through the x-axis', sep = '_'),
  #                 col = colors, lty = 1, lwd = 1, legend = T, xlab = 'X-column', ylab = 'Mean Expression', xpd = F,
  #                 trace.label = 'Br5292_Pos_300_Rep_2', fixed = TRUE)
  
  
  #-------------------------------------------------------------------------------------------------------
  ## Subject 2 ####
  ### Pos 0 ####
  
  Br_2_1_0 <- spots_frame[grepl('151669', as.character(spots_frame$ID), fixed = TRUE), ] #3661  8
  Br_2_1_0 <- Br_2_1_0[order(Br_2_1_0$y, Br_2_1_0$x), ] #ordering the dataframe for the value of x and then y 
  #table(Br_2_1_0[2: nrow(Br_2_1_0), ]$x - Br_2_1_0$x) #101
  Br_2_1_0$'Column' <- creat_columns(Br_2_1_0, 101) 
  #plot(Br_2_1_0$x, Br_2_1_0$y, col = rainbow(78)[Br_2_1_0$Column], xlab = 'x', ylab = 'y', main = 'Br5595 pos 0 rep 1') #plot the points colored for the column
  #interaction.plot(Br_2_1_0$Column, Br_2_1_0$Layer, response = Br_2_1_0$Expr, fun = mean, main = paste(gene, 'mean Expression going through the x-axis', sep = '_'),
  #                 col = colors, lty = 1, lwd = 1, legend = T, xlab = 'X-column', ylab = 'Mean Expression', xpd = F,
  #                 trace.label = 'Br5595_Pos_0_Rep_1', fixed = TRUE)
  
  Br_2_2_0 <- spots_frame[grepl('151670', as.character(spots_frame$ID), fixed = TRUE), ] #3498 8
  Br_2_2_0 <- Br_2_2_0[order(Br_2_2_0$y, Br_2_2_0$x), ]
  #table(Br_2_2_0[2: nrow(Br_2_2_0), ]$x - Br_2_2_0$x) # 98
  Br_2_2_0$'Column' <- creat_columns(Br_2_2_0, 98)
  #plot(Br_2_2_0$x, Br_2_2_0$y, col = rainbow(78)[Br_2_2_0$Column], xlab = 'x', ylab = 'y', main = 'Br5595 pos 0 rep 2')
  #interaction.plot(Br_2_2_0$Column, Br_2_2_0$Layer, response = Br_2_2_0$Expr, fun = mean, main = paste(gene, 'mean Expression going through the x-axis', sep = '_'),
  #                 col = colors, lty = 1, lwd = 1, legend = T, xlab = 'X-column', ylab = 'Mean Expression', xpd = F,
  #                 trace.label = 'Br5595_Pos_0_Rep_2', fixed = TRUE)
  
  ### Pos 300 ####
  Br_2_1_300 <- spots_frame[grepl('151671', as.character(spots_frame$ID), fixed = TRUE), ] # 4789 8
  Br_2_1_300 <- Br_2_1_300[order(Br_2_1_300$y, Br_2_1_300$x), ] #order
  #table(Br_2_1_300[2:nrow(Br_2_1_300), ]$x - Br_2_1_300$x) # 101
  Br_2_1_300$'Column' <- creat_columns(Br_2_1_300, 101)
  #plot(Br_2_1_300$x, Br_2_1_300$y, col = rainbow(78)[Br_2_1_300$Column], xlab = 'x', ylab = 'y', main = 'Br5595 pos 300 rep 1')
  #interaction.plot(Br_2_1_300$Column, Br_2_1_300$Layer, response = Br_2_1_300$Expr, fun = mean, main = paste(gene, 'mean Expression going through the x-axis', sep = '_'),
  #                 col = colors, lty = 1, lwd = 1, legend = T, xlab = 'X-column', ylab = 'Mean Expression', xpd = F,
  #                 trace.label = 'Br5595_Pos_300_Rep_1', fixed = TRUE)
  
  Br_2_2_300 <- spots_frame[grepl('151672', as.character(spots_frame$ID), fixed = TRUE), ] #4634 8
  Br_2_2_300 <- Br_2_2_300[order(Br_2_2_300$y, Br_2_2_300$x), ]
  #table(Br_2_2_300[2:nrow(Br_2_2_300), ]$x - Br_2_2_300$x)#98 but there are low values that stain this value
  Br_2_2_300$'Column' <- creat_columns(Br_2_2_300, 98)
  #plot(Br_2_2_300$x, Br_2_2_300$y, col = rainbow(78)[Br_2_2_300$Column], xlab = 'x', ylab = 'y', main = 'Br5595 pos 300 rep 2')
  #interaction.plot(Br_2_2_300$Column, Br_2_2_300$Layer, response = Br_2_2_300$Expr, fun = mean, main = paste(gene, 'mean Expression going through the x-axis', sep = '_'),
  #                col = colors, lty = 1, lwd = 1, legend = T, xlab = 'X-column', ylab = 'Mean Expression', xpd = F,
  #                 trace.label = 'Br5595_Pos_300_Rep_2', fixed = TRUE)
  
  #------------------------------------------------------------------------------------------------------- 
  ## Subject 3 ####
  ### Pos 0 ####
  Br_3_1_0 <- spots_frame[grepl('151673', as.character(spots_frame$ID), fixed = TRUE), ] #3661  8
  Br_3_1_0 <- Br_3_1_0[order(Br_3_1_0$y, Br_3_1_0$x), ] #ordering the dataframe for the value of x and then y 
  #table(Br_3_1_0[2: nrow(Br_3_1_0), ]$x - Br_3_1_0$x) #63 but there are lots of low values that stain this value
  Br_3_1_0$'Column' <- creat_columns(Br_3_1_0, 63) 
  #plot(Br_3_1_0$x, Br_3_1_0$y, col = rainbow(78)[Br_3_1_0$Column], xlab = 'x', ylab = 'y', main = 'Br8100 pos 0 rep 1') #plot the points colored for the column
  #interaction.plot(Br_3_1_0$Column, Br_3_1_0$Layer, response = Br_3_1_0$Expr, fun = mean, main = paste(gene, 'mean Expression going through the x-axis', sep = '_'),
  #                 col = colors, lty = 1, lwd = 1, legend = T, xlab = 'X-column', ylab = 'Mean Expression', xpd = F,
  #                 trace.label = 'Br8100_Pos_0_Rep_1', fixed = TRUE)
  
  Br_3_2_0 <- spots_frame[grepl('151674', as.character(spots_frame$ID), fixed = TRUE), ] #3498 8
  Br_3_2_0 <- Br_3_2_0[order(Br_3_2_0$y, Br_3_2_0$x), ]
  #table(Br_3_2_0[2: nrow(Br_3_2_0), ]$x - Br_3_2_0$x) # 62 but there are lots of low values that stain this value
  Br_3_2_0$'Column' <- creat_columns(Br_3_2_0, 62)
  #plot(Br_3_2_0$x, Br_3_2_0$y, col = rainbow(78)[Br_3_2_0$Column], xlab = 'x', ylab = 'y', main = 'Br8100 pos 0 rep 2')
  #interaction.plot(Br_3_2_0$Column, Br_3_2_0$Layer, response = Br_3_2_0$Expr, fun = mean, main = paste(gene, 'mean Expression going through the x-axis', sep = '_'),
  #                 col = colors, lty = 1, lwd = 1, legend = T, xlab = 'X-column', ylab = 'Mean Expression', xpd = F,
  #                 trace.label = 'Br8100_Pos_0_Rep_2', fixed = TRUE)
  
  ### Pos 300 ####
  Br_3_1_300 <- spots_frame[grepl('151675', as.character(spots_frame$ID), fixed = TRUE), ] # 4789 8
  Br_3_1_300 <- Br_3_1_300[order(Br_3_1_300$y, Br_3_1_300$x), ] #order
  #table(Br_3_1_300[2:nrow(Br_3_1_300), ]$x - Br_3_1_300$x) # 62 but there are lots of low values that stain this value
  Br_3_1_300$'Column' <- creat_columns(Br_3_1_300, 62)
  #plot(Br_3_1_300$x, Br_3_1_300$y, col = rainbow(78)[Br_3_1_300$Column], xlab = 'x', ylab = 'y', main = 'Br8100 pos 300 rep 1')
  #interaction.plot(Br_3_1_300$Column, Br_3_1_300$Layer, response = Br_3_1_300$Expr, fun = mean, main = paste(gene, 'mean Expression going through the x-axis', sep = '_'),
  #                 col = colors, lty = 1, lwd = 1, legend = T, xlab = 'X-column', ylab = 'Mean Expression', xpd = F,
  #                 trace.label = 'Br8100_Pos_300_Rep_1', fixed = TRUE)
  
  
  Br_3_2_300 <- spots_frame[grepl('151676', as.character(spots_frame$ID), fixed = TRUE), ] #4634 8
  Br_3_2_300 <- Br_3_2_300[order(Br_3_2_300$y, Br_3_2_300$x), ]
  #table(Br_3_2_300[2:nrow(Br_3_2_300), ]$x - Br_3_2_300$x)#59 but there are low values that stain this value 
  Br_3_2_300$'Column' <- creat_columns(Br_3_2_300, 59)
  #plot(Br_3_2_300$x, Br_3_2_300$y, col = rainbow(78)[Br_3_2_300$Column], xlab = 'x', ylab = 'y', main = 'Br8100 pos 300 rep 2')
  #interaction.plot(Br_3_2_300$Column, Br_3_2_300$Layer, response = Br_3_2_300$Expr, fun = mean, main = paste(gene, 'mean Expression going through the x-axis', sep = '_'),
  #                 col = colors, lty = 1, lwd = 1, legend = T, xlab = 'X-column', ylab = 'Mean Expression', xpd = F,
  #                 trace.label = 'Br8100_Pos_300_Rep_2', fixed = TRUE)
  
  #-------------------------------------------------------------------------------------------------------
  ## General spots_frame ####
  #Apply the new computation results to spots_frame 
  spots_frame_col <- rbind(Br_1_1_0, Br_1_2_0, Br_1_1_300, Br_1_2_300, Br_2_1_0, 
                           Br_2_2_0, Br_2_1_300, Br_2_2_300, Br_3_1_0, Br_3_2_0, Br_3_1_300, Br_3_2_300)  # 47681     10
  #But in this way you still need the expression of the gene --> the spots' disposition is different --> do it before  
  return(spots_frame_col)
}

#Function to get the expression of the desired gene --> DO NOT REMOVE THE NA AND DO NOT ORDER FOR THE SPOT NAME 
get_expression <- function(gene, frame){
  gene_id <-  which(rowData(spe)$'gene_name' == gene) #find the index corresponding to the gene of interest --> row 
  expr <- assays(spe)$counts[gene_id, ] #retrieve the expression information of that row --> 47681 --> the row is the first column of colData so coherence
  frame <- cbind(frame, 'Expr' = expr) #from 47681 6 to 47681 7 
  return(frame)
}


#---------------------------------------------------------------------------------------------------------

#SNAP25 --> GM marker
frame_SNAP25 <- get_expression('SNAP25', spots_frame) # 47329     6
frame_SNAP25_all <- apply_class(frame_SNAP25, 'SNAP25')


#MOBP --> Wm/ligodendrocytes marker 
frame_MOBP <- get_expression('MOBP', spots_frame) # 47329     6
frame_MOBP_all <- apply_class(frame_MOBP, 'MOBP')

#HPCAL1 --> L2 and striking expression decrease through L4, L5, L6 
frame_HPCAL1 <- get_expression('HPCAL1', spots_frame)# 47329     6
frame_HPCAL1_all <- apply_class(frame_HPCAL1, 'HPCAL1')

#KRT17 --> L6 and progressive decrease in expression from L6 through L5, L4, L3, L2
frame_KRT17 <- get_expression('KRT17', spots_frame)# 47329     6
frame_KRT17_all <- apply_class(frame_KRT17, 'KRT17')

#HBB, IGKC, NPY --> SVG (subunit of hemoglobin, inhibitory interneurons, constant region of light chian in immunoglobins)
frame_IGKC <- get_expression('IGKC', spots_frame)# 47329     6
frame_IGKC_all <- apply_class(frame_IGKC, 'IGKC')



##K-means trial ####

library(FCPS)

set.seed(1234)
km_SNAP25_1 <- kmeans(frame_SNAP25[frame_SNAP25$j_subj == 'Br5292', ]$Expr, 7, nstart = 250)
table(frame_SNAP25[frame_SNAP25$j_subj == 'Br5292', ]$Layer, km_SNAP25_1$cluster)
km_SNAP25_1 <- kmeans(frame_SNAP25[frame_SNAP25$j_subj == 'Br5595', ]$Expr, 7, nstart = 250)
table(frame_SNAP25[frame_SNAP25$j_subj == 'Br5595', ]$Layer, km_SNAP25_1$cluster)
km_SNAP25_1 <- kmeans(frame_SNAP25[frame_SNAP25$j_subj == 'Br8100', ]$Expr, 7, nstart = 250)
table(frame_SNAP25[frame_SNAP25$j_subj == 'Br8100', ]$Layer, km_SNAP25_1$cluster)

km_SNAP25_1 <- kmeans(frame_SNAP25[frame_SNAP25$j_subj == 'Br5292' & frame_SNAP25$Pos == '0', ]$Expr, 7, nstart = 250)
table(frame_SNAP25[frame_SNAP25$j_subj == 'Br5292'  & frame_SNAP25$Pos == '0', ]$Layer, km_SNAP25_1$cluster)
km_SNAP25_1 <- kmeans(frame_SNAP25[frame_SNAP25$j_subj == 'Br5292' & frame_SNAP25$Pos == '300', ]$Expr, 7, nstart = 250)
table(frame_SNAP25[frame_SNAP25$j_subj == 'Br5292'  & frame_SNAP25$Pos == '300', ]$Layer, km_SNAP25_1$cluster)



