library(motifStack)
library(jsonlite)

aa <- c("G", "A", "L", "M", "F", "W", "K", "Q", "E", "S", "P", "V", "I", "C", "Y", "H", "R", "N", "D", "T")
colset <- rep("gray", 20)
total_entries <- 79248

analyse_mut <- function(data, ref) {
  len <- nchar(ref)
  pfm <- matrix(0, nrow = len, ncol = length(aa))
  colnames(pfm) <- aa
  rownames(pfm) <- 1:len
  rates <- rep(0, len)
  vars <- rep(1, len)
  
  for (i in 1:len) {
    for_cur_res <- data[data[,"position"] == i, c("entry", "count")]
    pfm[i, substr(ref, i, i)] <- total_entries - sum(for_cur_res[, "count"])
    
    if (nrow(for_cur_res) == 0) {
      next
    }
    
    for (j in 1:nrow(for_cur_res)) {
      pfm[i, for_cur_res[j, "entry"]] <- for_cur_res[j, "count"]
    }
    
    rates[i] <- sum(for_cur_res[, "count"])
    vars[i] <- nrow(for_cur_res)
  }
  
  return(list(pfm=t(pfm/total_entries), rates=rates, vars=vars))
}

files <- list.files(pattern = "\\.tsv")

ref <- read.table("ref", header = F, sep = "\t", quote = "")
seqs <- ref[,2]
names(seqs) <- ref[,1]
color_data <- ref[,2]
color_data <- lapply(ref[,3], fromJSON)
names(color_data) <- ref[,1]

for (file in files) {
  data <- read.table(file, header = T)
  name <- unlist(strsplit(file, "\\.")[[1]])[1]
  seq <- unlist(seqs[name])
  res <- analyse_mut(data, seq)
  
  svg(paste0(name, "_pfm.svg"), height = 3, width = len/2)
  plotMotifLogo(res$pfm, p=rep(0.05, 20), colset = colset, ic.scale = T)
  dev.off()
  
  colors <- rep("black", nchar(seq))
  col_spec <- color_data[[name]]
  for (col in names(col_spec)) {
    colors[unlist(col_spec[col])] <- col
  }
  
  svg(paste0(name, "_rates.svg"), height = 7, width = 0.075*len)
  plot(1:nchar(seq), col=colors, res$rates, type = "h", xlab = "position", ylab = "mutation rate")
  dev.off()
  
  svg(paste0(name, "_vars.svg"), height = 7, width = 0.075*len)
  plot(1:nchar(seq), col=colors, res$vars, type = "h", xlab = "position", ylab = "mutation rate")
  dev.off()
}
