#Create an R script to compute different statistics and plot the results for genes of interest ####
#Load the library 
library(spatialLIBD)
spe <- fetch_data('spe')

#Create the matrix at spot level retrieving information about each single spot. 
#Data frame with 4781 rows and 6 columns corresponding to the spot name, the layer of belonging, the position inside the patient, the patient name, the cell counts and the key for the univocity of the spot 

spots_frame <- data.frame('Spot' = rownames(colData(spe)), 'Layer' = spe$layer_guess_reordered_short, 'Pos' = spe$position,
                          'j_subj' = spe$subject, 'Cells' = spe$cell_count, 'ID' = spe$key)



#get_expression is a function that takes in input the gene name of interest and the general spots_frame data structure. 
#It returns a data frame that contains an additional column for the logcounts expression of the gene of interest
#The dataframe is filtered for any NA (Layer, expression or whatever) and it is ordered according to the spot name and the layer
get_expression <- function(gene, frame){
  gene_id <-  which(rowData(spe)$'gene_name' == gene) #find the index corresponding to the gene of interest --> row 
  expr <- assays(spe)$counts[gene_id, ] #retrieve the expression information of that row --> 47681 --> the row is the first column of colData so coherence
  
  frame <- cbind(frame, 'Expr' = expr) #from 47681 6 to 47681 7 
  frame <- frame[complete.cases(frame), ] #Keep the spots that are not NA in every column --> 47329 7. Here the rownames are kept as in the original --> check in appendix for more 
  frame <- frame[order(frame$Spot, frame$Layer), ] #order the dataframe according to 1) spots name, 2) layer
  return(frame)
}



#remove_0 is a function that takes a dataframe and remove the values of Expression != 0
remove_0 <- function(frame){
  return(frame[frame$Expr !=0, ])
}

#get_numbers is a function that takes in input a dataframe and prints same numerical considerations about the expression levels of the frame
get_numbers <- function(frame){
  #Remove the Expr = 0
  frame_no0 <- remove_0(frame)
  
  print('Range of expression: ')
  print(range(frame$Expr))
  
  #Tables 
  table(frame$Expr, frame$Layer)
  table(frame$Expr, frame$j_subj)
  table(frame$Expr, frame$Pos)
  table(frame$Expr, frame$Layer, frame$j_subj)
  table(frame$Expr, frame$Layer, frame$Pos)
  table(frame$Expr, frame$Pos, frame$j_subj)
  
  #Tables Expr > 0
  table(frame_no0$Expr, frame_no0$Layer)
  table(frame_no0$Expr, frame_no0$j_subj)
  table(frame_no0$Expr, frame_no0$Pos)
  table(frame_no0$Expr, frame_no0$Layer, frame_no0$j_subj)
  table(frame_no0$Expr, frame_no0$Layer, frame_no0$Pos)
  table(frame_no0$Expr, frame_no0$Pos, frame_no0$j_subj)
  
}


#plot_box is a function that takes a dataframe and performs some box-plots over it and the filtered version of the dataframe (frame_no0). Input: dataframe and gene name
plot_box <- function(frame, gene){
  #Remove the Expr = 0 elements
  frame_no0 <- remove_0(frame)
  colors <- c("#F0027F","#377EB8","#4DAF4A","#984EA3","#FFD700","#FF7F00","#666666")
  
  #Plot the expression given the Subject with all spots and with the filtered ones 
  par(mfrow = c(1,2))
  boxplot(frame$Expr ~ frame$j_subj, xlab = 'Subject', ylab = 'Expr', main = paste(gene, 'expression in each Subject', sep = '_'))
  abline(h = 0, col = 'red', lty = 2, lwd = 0.5 )
  boxplot(frame_no0$Expr ~ frame_no0$j_subj, xlab = 'Subject', ylab = 'Expr', main = paste(gene, 'expression(> 0) in each Subject', sep = '_'))
  abline(h = 0, col = 'red', lty = 2, lwd = 0.5 )
  
  #Plot the expression given the Layer with all spots and with the filtered ones 
  par(mfrow = c(1,2))
  boxplot(frame$Expr ~ frame$Layer, xlab = 'Layer', ylab = 'Expr', main = paste(gene, 'expression in each Layer', sep = '_'),  ylim = c(0, (max(frame$Expr) + 1)))
  abline(h = 0, col = 'red', lty = 2, lwd = 0.5 )
  boxplot(frame_no0$Expr ~ frame_no0$Layer, xlab = 'Layer', ylab = 'Expr', main = paste(gene, 'expression (> 0) in each Layer', sep = '_'),  ylim = c(0, (max(frame_no0$Expr) + 1)))
  abline(h = 0, col = 'red', lty = 2, lwd = 0.5 )
  
  #Plot the expression given the Position with all spots and with the filtered ones 
  par(mfrow = c(1,2))
  boxplot(frame$Expr ~ frame$Pos, xlab = 'Position', ylab = 'Expr',main = paste(gene, 'expression in each Posotion', sep = '_'),  ylim = c(0, (max(frame$Expr) + 1)) )
  abline(h = 0, col = 'red', lty = 2, lwd = 0.5 )
  boxplot(frame_no0$Expr ~ frame_no0$Pos, xlab = 'Position', ylab = 'Expr', main = paste(gene, 'expression (> 0) in each Position', sep = '_'),  ylim = c(0, (max(frame_no0$Expr) + 1)))
  abline(h = 0, col = 'red', lty = 2, lwd = 0.5 )
  
  #PLot the expression given the cell counts with all spots and the filtered ones 
  par(mfrow = c(1,2))
  boxplot(frame$Expr ~ frame$Cells, xlab = 'Cells count', ylab = 'Expr',main = paste(gene, 'expression for each cells count', sep = '_'),  ylim = c(0, (max(frame$Expr) + 1)))
  abline(h = 0, col = 'red', lty = 2, lwd = 0.5 )
  boxplot(frame_no0$Expr ~ frame_no0$Cells, xlab = 'Cells count', ylab = 'Expr', main = paste(gene, 'expression (> 0) for each cells count', sep = '_'),  ylim = c(0, (max(frame_no0$Expr) + 1)))
  abline(h = 0, col = 'red', lty = 2, lwd = 0.5 )
  
  #Plot the boxplots for Expression given the position, position and layer, position and layer and subject
  par(mfrow = c(1,3))
  
  boxplot(frame$Expr ~ frame$Pos, xlab = 'Position', ylab = 'Expression', ylim = c(0, max(frame$Expr)), main = paste(gene, 'expression', sep = '_'))
  boxplot(frame$Expr ~ frame$Pos + frame$Layer, xlab = 'Position:Layer', ylab = 'Expression', ylim = c(0, max(frame$Expr)), main = paste(gene, 'expression', sep = '_'))
  abline(v = c(2.5, 4.5, 6.5, 8.5, 10.5, 12.5), col = colors[2:length(colors)], lty = 2)
  boxplot(frame$Expr ~ frame$Pos + frame$Layer + frame$j_subj, xlab = 'Position:Layer:Subject', ylab = 'Expression', ylim = c(0, max(frame$Expr)), main = paste(gene, 'expression', sep = '_'))
  abline( v = c(14.5, 28.5), col = c('red', 'blue'), lty = 2)
  
  
  boxplot(frame_no0$Expr ~ frame_no0$Pos, xlab = 'Position', ylab = 'Expression', ylim = c(0, max(frame_no0$Expr)), main = paste(gene, 'expression (> 0)', sep = '_'))
  boxplot(frame_no0$Expr ~ frame_no0$Pos + frame_no0$Layer, xlab = 'Position:Layer', ylab = 'Expression', ylim = c(0, max(frame_no0$Expr)), main = paste(gene, 'expression (> 0)', sep = '_'))
  abline(v = c(2.5, 4.5, 6.5, 8.5, 10.5, 12.5), col = colors[2:length(colors)], lty = 2)
  boxplot(frame_no0$Expr ~ frame_no0$Pos + frame_no0$Layer + frame_no0$j_subj, xlab = 'Position:Leyer:Subject', ylab = 'Expression', ylim = c(0, max(frame_no0$Expr)), main = paste(gene, 'expression (> 0)', sep = '_'))
  abline( v = c(14.5, 28.5), col = c('red', 'blue'), lty = 2)
  
}


#plot_hist is a function that takes in input a dataframe and gives some histogram plots of it and of the filtered version of it. Input is the dataframe and the gene name  
plot_hist <- function(frame, gene){ 
  frame_no0 <- remove_0(frame)
  
  #Plot the hist of the expression, expression > 0 and the cells count 
  par(mfrow = c(1,3))
  hist(frame$Expr, xlab = 'Counts', ylab = 'Frequency',main = paste(gene, 'expression', sep = '_'), xlim = c(min(frame$Expr), (max(frame$Expr) + 1)), freq = FALSE)
  hist(frame_no0$Expr, xlab = 'Counts', ylab = 'Frequency',main = paste(gene, 'expression (> 0)', sep = '_'), xlim = c(min(frame_no0$Expr), (max(frame_no0$Expr) + 1)), freq = FALSE)
  hist(frame$Cells , xlab = 'Cells',  ylab = 'Frequency',main = 'Cells count', xlim = c(0, (max(frame$Cells) + 1)), freq = FALSE)
}  


#plot_interaction is a function that takes a dataframe and gives some interaction plots of it and of the filtered version of it. Input is dataframe and the gene name 
plot_interaction <- function(frame, gene){ 
  frame_no0 <- remove_0(frame)
  
  par(mfrow = c(1, 1))
  #Plot the mean expression given the layer and the subject in absolute values y
  interaction.plot(frame$Layer, frame$j_subj, response = frame$Expr, fun = mean, lty = 1, lwd = 1, xlab = 'Layer', 
                   ylab = 'Mean expression', main = paste(gene, 'expession in each Layer and Subject', sep = '_'), ylim = c(0.0, 30), xpd = F, trace.label = 'Subject', col = c('red', 'cornflowerblue', 'pink'), fixed = TRUE)
  #interaction.plot(frame_no0$Layer, frame_no0$j_subj, response = frame_no0$Expr, fun = mean, lty = 1, lwd = 1, xlab = 'Layer', 
  #                 ylab = 'Mean expression > 0', main = paste(gene, 'expession (>0) in each Layer and Subject', sep = '_'), ylim = c(0.0, 30), xpd = F, trace.label = 'Subject', col = c('red', 'cornflowerblue', 'pink'), fixed = TRUE)
  
  #Plot the mean expression given the layer and the subject in tailored y values
  interaction.plot(frame$Layer, frame$j_subj, response = frame$Expr, fun = mean, lty = 1, lwd = 1, xlab = 'Layer', 
                   ylab = 'Mean expression', main = paste(gene, 'expession in each Layer and Subject', sep = '_'), xpd = F, trace.label = 'Subject', col = c('red', 'cornflowerblue', 'pink'), fixed = TRUE)
  #interaction.plot(frame_no0$Layer, frame_no0$j_subj, response = frame_no0$Expr, fun = mean, lty = 1, lwd = 1 ,xlab = 'Layer', 
  #                 ylab = 'Mean expression > 0', main = paste(gene, 'expession (>0) in each Layer and Subject', sep = '_'), xpd = F, trace.label = 'Subject', col = c('red', 'cornflowerblue', 'pink'), fixed = TRUE) 
  
  
  
  #Plot the mean expression given the layer and position in absolute values of y
  interaction.plot(frame$Layer, frame$Pos, response = frame$Expr, fun = mean, lty = 1, lwd = 1,xlab = 'Layer', 
                   ylab = 'Mean expression', main = paste(gene, 'expession in each Layer and Position', sep = '_'), ylim = c(0.0, 30), xpd = F, trace.label = 'Position', col = c('orange', 'forestgreen'), fixed = TRUE)
  #interaction.plot(frame_no0$Layer, frame_no0$Pos, response = frame_no0$Expr, fun = mean, lty = 1, lwd = 1,xlab = 'Layer', 
  #                 ylab = 'Mean expression > 0', main = paste(gene, 'expession (>0) in each Layer and Position', sep = '_'), ylim = c(0.0, 30), xpd = F, trace.label = 'Position', col = c('orange', 'forestgreen'), fixed = TRUE)
  
  #Plot the mean expression given the layer and the subject in tailored y values
  interaction.plot(frame$Layer, frame$Pos, response = frame$Expr, fun = mean,lty = 1, lwd = 1 ,xlab = 'Layer', 
                   ylab = 'Mean expression', main = paste(gene, 'expession in each Layer and Position', sep = '_'), xpd = F, trace.label = 'Position', col = c('orange', 'forestgreen'), fixed = TRUE)
  #interaction.plot(frame_no0$Layer, frame_no0$Pos, response = frame_no0$Expr, fun = mean,lty = 1, lwd = 1 ,xlab = 'Layer', 
  #                 ylab = 'Mean expression > 0', main = paste(gene, 'expession (>0) in each Layer and Position', sep = '_'), xpd = F, trace.label = 'Position', col = c('orange', 'forestgreen'), fixed = TRUE) 
  
  
  
  #Plot the mean expression given the position and subject in absolute values of y
  interaction.plot(frame$Pos, frame$j_subj, response = frame$Expr, fun = mean,lty = 1, lwd = 1, xlab = 'Position', 
                   ylab = 'Mean expression', main = paste(gene, 'expession in each Position and Subject', sep = '_'), ylim = c(0.0, 30), xpd = F, trace.label = 'Subject', col = c('red', 'cornflowerblue', 'pink'), fixed = TRUE)
  #interaction.plot(frame_no0$Pos, frame_no0$j_subj, response = frame_no0$Expr, fun = mean, lty = 1, lwd = 1, xlab = 'Position', 
  #                 ylab = 'Mean expression > 0', main = paste(gene, 'expession (>0) in each Position and Subject', sep = '_'), ylim = c(0.0, 30), xpd = F, trace.label = 'Subject', col = c('red', 'cornflowerblue', 'pink'), fixed = TRUE)
  
  #Plot the mean expression given the position and the subject in tailored y values
  interaction.plot(frame$Pos, frame$j_subj, response = frame$Expr, fun = mean, lty = 1, lwd = 1, xlab = 'Position', 
                   ylab = 'Mean expression', main = paste(gene, 'expession in each Position and Subject', sep = '_'), xpd = F, trace.label = 'Subject', col = c('red', 'cornflowerblue', 'pink'), fixed = TRUE)
  #interaction.plot(frame_no0$Pos, frame_no0$j_subj, response = frame_no0$Expr, fun = mean, lty = 1, lwd = 1, xlab = 'Position', 
  #                 ylab = 'Mean expression > 0', main = paste(gene, 'expession (>0) in each Position and Subject', sep = '_'), xpd = F, trace.label = 'Subject', col = c('red', 'cornflowerblue', 'pink'), fixed = TRUE) 
  
  
  
  #Plot the mean and sd expression given the Layer and subject in absolute values of y
  interaction.plot(frame$Layer, frame$j_subj, response = frame$Expr, fun = mean, lty = 1, lwd = 1, xlab = 'Layer', 
                   ylab = 'Mean expression', main = paste(gene, 'expession in each Layer and Subject', sep = '_'), xpd = F, trace.label = 'Subject', col = c('red', 'cornflowerblue', 'pink'), fixed = TRUE)
  interaction.plot(frame$Layer, frame$j_subj, response = frame$Expr, fun = sd, lty = 1, lwd = 1, xlab = 'Layer', 
                   ylab = 'SD expression', main = paste(gene, 'expession in each Layer and Subject', sep = '_'), xpd = F, trace.label = 'Subject', col = c('red', 'cornflowerblue', 'pink'), fixed = TRUE)
  
  #interaction.plot(frame_no0$Layer, frame_no0$j_subj, response = frame_no0$Expr, fun = mean, lty = 1, lwd = 1, xlab = 'Layer', 
  #                 ylab = 'Mean expression (> 0)', main = paste(gene, 'expession (>0) in each Layer and Subject', sep = '_'), xpd = F, trace.label = 'Subject', col = c('red', 'cornflowerblue', 'pink'), fixed = TRUE) 
  #interaction.plot(frame_no0$Layer, frame_no0$j_subj, response = frame_no0$Expr, fun = sd, lty = 1, lwd = 1,  xlab = 'Layer', 
  #                 ylab = 'SD expression (> 0)', main = paste(gene, 'expession (>0) in each Layer and Subject', sep = '_'), xpd = F, trace.label = 'Subject', col = c('red', 'cornflowerblue', 'pink'), fixed = TRUE)
  
  
  #Plot the mean expression given the cell counts in a given layer
  par(mfrow = c(1, 1))
  interaction.plot(frame$Cells, frame$Layer, response = frame$Expr, fun = mean, 
                   lty = 1, lwd = 1, legend = T, xlab = 'Cells count', ylab = 'Mean expression', 
                   main = paste(gene, 'Mean expression given cells count in each Layer', sep = '_'), xpd = F, trace.label = 'Layer', col = colors, fixed = TRUE)
  #interaction.plot(frame_no0$Cells, frame_no0$Layer, response = frame_no0$Expr, fun = mean, 
  #                 lty = 1, lwd = 1,legend = T, xlab = 'Cells count', ylab = 'Mean expression', 
  #                 main = paste(gene, 'Mean expression given cells count in each Layer', sep = '_'),xpd = F, trace.label = 'Layer', col = colors, fixed = TRUE)
  
  #then for subject
  interaction.plot(frame[frame$j_subj == 'Br5292', ]$Cells, frame[frame$j_subj == 'Br5292', ]$Layer, response = frame[frame$j_subj == 'Br5292', ]$Expr, fun = mean, 
                   lty = 1, lwd = 1, legend = T, xlab = 'Cells count', ylab = 'Mean expression', main = paste(gene, 'Mean expression given cells count in each Layer (Br5292)', sep = '_'),
                   xpd = F, trace.label = 'Layer subject Br5292', col = colors, fixed = TRUE)
  interaction.plot(frame[frame$j_subj == 'Br5595', ]$Cells, frame[frame$j_subj == 'Br5595', ]$Layer, response = frame[frame$j_subj == 'Br5595', ]$Expr, fun = mean, 
                   lty = 1, lwd = 1, legend = T, xlab = 'Cells count', ylab = 'Mean expression', main = paste(gene, 'Mean expression given cells count in each Layer (Br5595)', sep = '_'), 
                   xpd = F, trace.label = 'Layer subject Br5595', col = colors, fixed = TRUE)
  interaction.plot(frame[frame$j_subj == 'Br8100', ]$Cells, frame[frame$j_subj == 'Br8100', ]$Layer, response = frame[frame$j_subj == 'Br8100', ]$Expr, fun = mean, 
                   lty = 1, lwd = 1, legend = T, xlab = 'Cells count', ylab = 'Mean expression', main = paste(gene, 'Mean expression given cells count in each Layer (Br8100)', sep = '_'), 
                   xpd = F, trace.label = 'Layer subject Br8100', col = colors, fixed = TRUE)
  
} 



#As interaction plot but at single level --> visualize one at time
plot_int_zoom <- function(frame, gene){ 
  #Plot the mean Expression for each patient against the Layer --> Not same number of observations per layer per subject and per spots --> comparison is not really fair
  par(mfrow = c(3, 2))
  x_three <- table(frame$Expr, frame$Layer, frame$j_subj)
  for(j in 1:3){
    for(i in 1:ncol(x_three[, , j])){x_three[, i, j] <- sum(as.numeric(rownames(x_three[, , j])) * x_three[, i, j]) / sum(x_three[, i, j])}
    plot(x_three[1, , j], xlab = 'Layers', ylab = 'Mean_expression', main = paste('Subject', j, sep = '_'), type = 'b', ylim = c(0.0, 20.0))
    plot(x_three[1, , j], xlab = 'Layers', ylab = 'Mean_expression', main = paste('Subject', j, sep = '_'), type = 'b')
  }
  par(mfrow = c(1, 1))
  
  
  #Plot the mean expression for each patient against the Layer given the position 
  par(mfrow = c(2, 2))
  x_three <- table(frame$Expr, frame$Layer, frame$Pos, frame$j_subj)
  for(j in 1:3){
    for(pos in 1:2){
      for(i in 1:ncol(x_three[, , pos,j])){x_three[, i,pos, j] <- sum(as.numeric(rownames(x_three[, , pos, j])) * x_three[, i, pos, j]) / sum(x_three[, i,pos,j])}
      plot(x_three[1, ,pos, j], xlab = 'Layers', ylab = 'Mean_expression', main = paste(paste('Subject', j, sep = '_'), paste('Pos', pos, sep = '_'), sep = '_'), type = 'b', ylim = c(0.0,25.0))
      plot(x_three[1, ,pos, j], xlab = 'Layers', ylab = 'Mean_expression', main = paste(paste('Subject', j, sep = '_'), paste('Pos', pos, sep = '_'), sep = '_'), type = 'b')
    }
  }  
  par(mfrow = c(1, 1))
  
}

plot_cc <- function(frame){
  #Plot the cells count for each layer given the position 
  interaction.plot(frame$Layer, frame$Pos, response = frame$Cell, fun = mean, 
                   lty = 1, lwd = 1, legend = T, xlab = 'Layer', ylab = 'Mean Cells count', 
                   main = 'Mean cells count in each Layer given the position', xpd = F, trace.label = 'Position', col = c('orange', 'forestgreen'), fixed = TRUE)
  
  #Plot the cells count for each Layer given the position for each subject
  interaction.plot(frame[frame$j_subj == 'Br5292', ]$Layer, frame[frame$j_subj == 'Br5292', ]$Pos, response = frame[frame$j_subj == 'Br5292', ]$Cells, fun = mean, 
                   lty = 1, lwd = 1, legend = T, xlab = 'Layer', ylab = 'Mean Cells count', main = 'Mean cells count in each Layer given the position (Br5292)',
                   xpd = F, trace.label = 'Br5292', col = c('orange', 'forestgreen'), fixed = TRUE)
  interaction.plot(frame[frame$j_subj == 'Br5595', ]$Layer, frame[frame$j_subj == 'Br5595', ]$Pos, response = frame[frame$j_subj == 'Br5595', ]$Cells, fun = mean, 
                   lty = 1, lwd = 1, legend = T, xlab = 'Layer', ylab = 'Mean Cells count', main ='Mean cells count in each Layer given the position (Br5595)', 
                   xpd = F, trace.label = 'Br5595', col = c('orange', 'forestgreen'), fixed = TRUE)
  interaction.plot(frame[frame$j_subj == 'Br8100', ]$Layer, frame[frame$j_subj == 'Br8100', ]$Pos, response = frame[frame$j_subj == 'Br8100', ]$Cells, fun = mean, 
                   lty = 1, lwd = 1, legend = T, xlab = 'Layer', ylab = 'Mean Cells count', main = 'Mean cells count in each Layer given the position (Br8100)', 
                   xpd = F, trace.label = 'Br8100', col = c('orange', 'forestgreen'), fixed = TRUE)
  
}

#Create Dataframes for single genes of interest --> 47329, 6 are the dimentions all coherent between the data frames bcs the spots are those --> NA possible in Layer asignment only in this case
#SNAP25 --> GM marker
frame_SNAP25 <- get_expression('SNAP25', spots_frame) # 47329     6
plot_interaction(frame_SNAP25, 'SNAP25')

#MOBP --> Wm/ligodendrocytes marker 
frame_MOBP <- get_expression('MOBP', spots_frame) # 47329     6
plot_interaction(frame_MOBP, 'MOBP')

#AQP4--> cell-type marker in cortical layers and middle temporal gyrus and L1 marker 
frame_AQP4 <- get_expression('AQP4', spots_frame)# 47329     6

#HPCAL1 --> L2 and striking expression decrease through L4, L5, L6 
frame_HPCAL1 <- get_expression('HPCAL1', spots_frame)# 47329     6
plot_interaction(frame_HPCAL1, 'HPCAL1')

#FREM3 --> L3
frame_FREM3 <- get_expression('FREM3', spots_frame)#  47329     6

#TRABD2A --> L5
frame_TRABD2A <- get_expression('TRABD2A',spots_frame)#  47329     6

#KRT17 --> L6 and progressive decrease in expression from L6 through L5, L4, L3, L2
frame_KRT17 <- get_expression('KRT17', spots_frame)# 47329     6
plot_interaction(frame_KRT17, 'KRT17')

#FABP7 --> new markers found 
frame_FABP7 <- get_expression('FABP7', spots_frame)# 47329     6

#BCL11B --> rodent marker found now weak 
frame_BCL11B <- get_expression('BCL11B', spots_frame)# 47329, 6; Expr> 0 --> 1801, 6

#MKX--> ASD-dominant traits enriched in L5
frame_MKX <- get_expression('MKX', spots_frame)# 47329     6

#TCF4 --> neurodevelopmental delay enricehd in L2
frame_TCF4 <- get_expression('TCF4', spots_frame) # --> 47329     6

#HBB, IGKC, NPY --> SVG (subunit of hemoglobin, inhibitory interneurons, constant region of light chian in immunoglobins)
frame_IGKC <- get_expression('IGKC', spots_frame)# 47329     6
plot_interaction(frame_IGKC, 'IGKC')


##Appendix ####
#frame <- spots_frame
#frame <- cbind(frame, 'Expr' = expr) #from 47681 5 to 47681 6
#this <- rownames(frame)
#table(as.numeric(this[2: length(this)]) - as.numeric(this[1: (length(this) -1 )])) --> 1 : 47680 --> onne is missing bcs we start from the second postion 
#frame <- frame[complete.cases(frame), ] #Keep the spots that are not NA in every column --> 47329 6
#this <- rownames(frame)
#table(as.numeric(this[2: length(this)]) - as.numeric(this[1: (length(this) -1 )])) --> 1: 46983   2: 339    3: 5    4: 1  --> toto 47328. Since the frame is not yet ordered we are in the same order of before --> we keep the original indeces/rownmaes 
#dim(frame) --> 47329 6 

