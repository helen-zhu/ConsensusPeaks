
# Simple function for loading text input files
import.files = function(
  files,
  type = c("bed", "bam", "tsv"),
  sample.names = NULL
){
  
  # Error checking
  if(!all(file.exists)){stop("Please check that all files exist!")}
  type = match.arg(type)
  if(!is.null(sample.names) & length(sample.names) != length(files)){stop("The number of sample names does not match the number of files!")}
  
  # Importing txt files & bed files using readr
  if(type == "tsv"){
    elements = readr::read_delim(files)
  } else if (type == "bed"){
    elements = readr::read_delim(files)
  
  # Importing bam files using rsamtools
  } else if (type == "bam"){
    
  }
  
  # Adding Sample Names
  
  # Formatting to the Correct Format
  
  # Returns a data frame
  return(elements)
}