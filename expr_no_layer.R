#Create a script to access infos and display a simple symbol plot to show the expression of the gene of interest in samples

#load the library and create a background set of objects useful to study the structures
library(spatialLIBD)
spe <- fetch_data('spe')
indeces <- c("151507","151508", "151509", "151510", "151669", "151670", "151671", "151672", "151673", "151674", "151675", "151676")
sub_pos <- c("Br5292_pos0","Br5292_pos300", "Br5595_pos0", "Br5595_pos300", "Br8100_pos0","Br8100_pos300")
matrix_dim <- matrix(data = c(0, 0), nrow = 12, ncol = 2, byrow = T)
rownames(matrix_dim) <- indeces
colnames(matrix_dim) <- c('genes', 'spots')
for(i in 1:12){
  matrix_dim[i, ] <- c(as.integer(dim(spe[, spe$sample_id == indeces[i]])))
}
matrix_dim <- cbind(matrix_dim, c(rep(unique(spe$subject_position), each = '2'))) #--> trasforma tutto in char 
#use as.integer(matrix_dim[1, 1:2]) to retrieve 
colnames(matrix_dim)[3] <- 'subj_pos'



##Plots of gene expression for sub_pos ####
#List of putative genes 
list_of_interest <- c('SNAP25', 'MOBP', 'IGKC')
list_genes_ASD_L5 <- c('TBR1', 'SATB1', 'ANK2', 'RORB', 'MKX', 'CELF4', 'PPP5C', 'AP2S1')
list_genes_ASD_L2 <- c('TCF4', 'CACNA1E', 'MYT1L', 'SCN2A', 'TBL1XR1', 'NR3C2', 'SYNGAP1', 'GRIN2B', 'IRF2BPL', 'GABRB3', 'RAI1', 'ADNP')

#gene_expr is a function to retrieve the expression of a gene in single pairs and plot the mean expression at sub_pos level. It gets in input the gene of interest's name
#and parses the spe object to find the coordinates of points in samples pairs of replicates (ex. Br5292 pair at postion 0 and pair at position 300)
#and the expression of the gene (counts and logcounts). The parsing is performed by slinding through the spe object by indeces stored in the matrix_dim matrix created before. 
#The concordance between source and data is kept for the ordered structure of the spe object (matrices of matrices). 
#The function gives in output the pair of plots for the position 0 and position 300 for each subject. 
#The mean expression is computed considering same named spots and if one has NA the mean is set to NA. All NAs are removed. No layer structuring is presented.
#The Expression is represented in terms of logcounts (20*logcounts for the diameter of the circles). Blue stands for expression > 0, red for exprssion <= 0. 
gene_expr <- function(the_gene){
  
  spce <- as.data.frame(spatialCoords(spe))
  basic <- data.frame(colnames(spe), spce$pxl_col_in_fullres, spce$pxl_row_in_fullres)
  names(basic) <- c('spot', 'x', 'y') 
  
  #Find the gene index
  index_gene <- which(rowData(spe)$gene_name == the_gene)
  spatial_counts_gene <- data.frame(assays(spe)$counts[index_gene, ], assays(spe)$logcounts[index_gene, ])
  spatial_counts_gene <- cbind(basic, spatial_counts_gene) #there is correspondance of indeces 
  names(spatial_counts_gene) <- c('spot','x', 'y', 'counts', 'logcounts')  
  
  #Retrieve the expressions in terms of logcounts and counts. Compute the mean expression for replicates at same position (by = 2). If a replicate has NA in that spot --> NA the mean. All NAs are then removed for plotting 
  par(mfrow = c(1,2))
  last_index <- 0 #first time is 0. Then it is used to characterize the start of the following couple of spots sets
  for (j in seq(1, 12, by = 2)){
    log_M <- count_M <- NULL
    index_first <- as.integer(matrix_dim[j, 2]) #access the first replicate's number of spots
    index_second <- as.integer(matrix_dim[j+1, 2]) #acess the second replicates's number of spots
    subset <- spatial_counts_gene[(last_index+1):(last_index + index_first + index_second), ]
    index_list <- NULL
    #Find the spots that have the same name 
    for (i in 1: index_first){
      index_list <- append(index_list, which(subset$spot == subset$spot[i])[2])
    }
    for (i in 1:length(index_list)){
      if (!is.na(index_list[i])){
        log_M[i] <- mean(c(subset$logcounts[i], subset$logcounts[index_list[i]]))
        count_M[i] <- mean(c(subset$counts[i], subset$counts[index_list[i]]))
      } else { log_M[i] <- NA
      count_M[i] <- NA
      }
      
    }
    last_index <- last_index + nrow(subset) #so to know how long it was --> next starting point 
    subset <- subset[1: index_first, ] #keep the coordinates of the first replicate only --> there is not a division in layers so no crucial relative positions alterations
    subset <- cbind(subset, count_M, log_M)
    names(subset)[6:7] <- c('count_M', 'log_M')
    subset <- subset[complete.cases(subset), ]
    #PLot this subset --> symbols to plot the expression as the diameter of circles. Color discrimination of the expression --> blue > 0, red <= 0. Logcounts used. 
    symbols(subset$x, subset$y, circles = 20*(subset$logcounts), fg = ifelse(subset$logcounts <= 0, 'red', 'blue'), inches = F, xlab = 'x', ylab = 'y', main = paste(the_gene,sub_pos[floor(j/2)+1], sep = '_'))
    remove(subset)
  }
}

gene_expr('SNAP25')

