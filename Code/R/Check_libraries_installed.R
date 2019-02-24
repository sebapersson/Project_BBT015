# This program checks that all required packages are installed. If a package isn't installed the program will 
# install it locally. Since dependencies are tricky among the bioconductor packages the program will exit if 
# bioconductor package isn't installed (and won't try to install it).

if(!require(Rsamtools)){
  message("Rsamtool isn't installed, program will exit")
  quit(status=1)
}

if(!require(GenomicAlignments)){
  message("GenomicAlignments isn't installed, program will exit")
  quit(status=1)
}

if(!require(GenomicFeatures)){
  write("GenomicFeatures isn't installed, program will exit")
  quit(status=1)
}

if(!require(HDF5Array)){
  message("HDF5Array isn't installed, program will exit")
  quit(status=1)
}
    
if(!require(DESeq2)){
  message("Rsamtool isn't installed, program will exit")
  quit(status=1)
}

if(!require(pheatmap)){
  print("pheatmap isn't installed, will install:")
  install.packages("pheatmap")
}

if(!require(RColorBrewer)){
  print("RColorBrewer isn't installed, will install:")
  install.packages("RColorBrewer")
}

if(!require(ggplot2)){
  print("ggplot2 isn't installed")
  install.packages("ggplot2")
}

if(!require(PoiClaClu)){
  print("PoiClaClu isn't installed, will install:")
  install.packages("PoiClaClu")
}

message("All packages are installed")

quit(status = 0)






