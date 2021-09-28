load.DIA.data <- function(syn){
  data <- querySynapseTable(syn)
  columns <- c("Sample", "PatientID", "CellType", 
               "TrametinibTreatment", "CellCount", "Date")
  meta <- data[, columns] %>%
    unique()
  meta <- meta[order(meta$CellType, meta$Trametinib), ]
  
  
  return(list("Long form data" = data, "Metadata" = meta))
}
