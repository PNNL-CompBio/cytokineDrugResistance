load.DIA.data <- function(syn){
  data <- querySynapseTable(syn)
  columns <- c("Sample", "PatientID", "CellType", 
               "TrametinibTreatment", "CellCount", "Date")
  meta <- data[, columns] %>%
    unique()
  
  return(list("Long form data" = data, "Metadata" = meta))
}