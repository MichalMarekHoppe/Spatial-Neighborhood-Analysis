# https://github.com/MichalMarekHoppe/Spatial-analysis---calculating-delta-between-cellular-phenotypes.git
# author: Michal Marek Hoppe (mmlhoppe@gmail.com)
# source: https://aacrjournals.org/cancerdiscovery/article/13/5/1144/726201/Patterns-of-Oncogene-Coexpression-at-Single-Cell
# 
# Spatial analysis - calculating delta % between cellular phenotypes among k nearest neighbours
# spatial_sample.txt -  sample data from three high-power microscopic images of DLBCL cases. input data consists of
# a) id - sample id
# b) label - different phenotypes / grouping
# d) x and y - spatial coordiantes for each cell.
# 
# Number of k neighbours can be defined at the beginning of the process.
# 
# Summary results sheet provides results for interactions among all phenotypes for each sample. 

## START

# manually define number of k neighbours to be interrogated
k <- 20

library(RANN) 
# import data
data <- read.delim("spatial_sample.txt",
                   sep = "\t")

# define images
id <- sort(unique(data[,"id"]))

# define cell populations
labels <- sort(unique(data[,"label"]),
               decreasing =  TRUE)

# create summary results sheet
results <- 
  as.data.frame(matrix(ncol = (length(labels)^2)*2 + length(labels),
                       nrow = length(id),
                       dimnames = list(id,
                                       c(paste0("perc_", labels),
                                         paste0("subp_",
                                                rep(labels, each = 9),
                                                "_",
                                                rep(labels, 9)),
                                         paste0("delt_",
                                                rep(labels, each = 9),
                                                "_",
                                                rep(labels, 9))))))
for (i in seq_along(id)){
  #populate summary "perc" results
  results[rownames(results) == id[i], grep("perc_",
                                           colnames(results))] <- 
    table(data[data[,"id"] == id[i],
               "label"])[match(names(table(data[data[,"id"] == id[i],
                                                "label"])),
                               gsub("perc_",
                                    "",
                                    colnames(results[,grep("perc_",
                                                           colnames(results))])))] / 
      nrow(data[data[,"id"] == id[i],]) * 100
  
  # create individual sample result data
  result_id <- data[data[,"id"] == id[i],]
  rownames(result_id) <- 1:nrow(result_id)
  
  # identify nearest neighbours
  nearest_id <- RANN::nn2(result_id[,c("x", "y")],
                    result_id[,c("x", "y")],
                    k = k + 1)[["nn.idx"]] # +1 as always to account for the centre point
  nearest_id <- nearest_id[,-1] #remove the centre point
  colnames(nearest_id) <- paste0("Neighbour_", seq(1:k), "_id")
  
  # identify cell phenotype
  nearest_label <- nearest_id
  colnames(nearest_label) <- gsub("id", "label", colnames(nearest_label))
  for (n in 1:ncol(nearest_label)){
    nearest_label[,n] <- result_id[match(nearest_label[,n],
                                         as.integer(rownames(result_id))),"label"]
  }; rm(n)
  
  # count local populations
  local_label_summaries <- 
    as.data.frame(matrix(ncol = length(labels) * 2,
                                       nrow = nrow(result_id),
                                       dimnames = list(rownames(result_id),
                                                       c(paste0(labels,
                                                                "_k",
                                                                k,
                                                                "_perc"),
                                                         paste0(labels,
                                                                "_k",
                                                                k,
                                                                "_delta")))))
  local_label_summaries[,grep("perc",
                              colnames(local_label_summaries))] <- 0
  for (n in 1:nrow(result_id)){
     # calculate local percentages
     local_label_summaries[n,
                           match(names(table(nearest_label[n,])),
                                 gsub(paste0("_k",
                                          k,
                                          "_perc"),
                                      "",
                                      colnames(local_label_summaries[
                                        n,grep("perc",
                                               colnames(local_label_summaries))])))] <- 
       table(nearest_label[n,]) / k * 100
     
     # calculate local delta
     local_label_summaries[n, grep("delta",
                                   colnames(local_label_summaries))] <-
       local_label_summaries[n, grep("perc",
                                     colnames(local_label_summaries))] - 
         results[rownames(results) == id[i], grep("perc", colnames(results))]
     cat(paste0("Calculating local percentages for ", id[i], " - ", 
                format(round(n/nrow(result_id)*100,1), nsmall = 1), "% \t"), "\r")
  }; rm(n)
  # collate per-cell results
  result_id <- cbind(result_id, local_label_summaries, nearest_id, nearest_label)
  write.csv(result_id, paste0("results_", id[i], "_k", k, ".csv"))
  
  # fill-in "mean" and "delta" summary results
  for (l in labels) {
  results[rownames(results) == id[i], 
          intersect(grep("subp_", 
                         colnames(results)),
                       seq(1:ncol(results))[substr(colnames(results),
                                                6,
                                                11) == l])] <- 
    colMeans(result_id[result_id[,"label"] == l,
                       grep("perc",
                            colnames(result_id))])
  
  results[rownames(results) == id[i], 
          intersect(grep("delt_", 
                         colnames(results)),
                    seq(1:ncol(results))[substr(colnames(results),
                                                6,
                                                11) == l])] <- 
    colMeans(result_id[result_id[,"label"] == l,
                       grep("delta",
                            colnames(result_id))])
  }; rm(l)
  cat(paste0("Done ", id[i], "! ", i, "/", length(id), " ",
             format(round(i/length(id)*100,0), nsmall = 0),
             "%                                                 "), "\r")
  cat("\n")
}; rm(i, data, nearest_id, nearest_label, local_label_summaries, result_id)
write.csv(results, paste0("results_summary", "_k", k, ".csv"))
rm(k, id, labels, results)

  
