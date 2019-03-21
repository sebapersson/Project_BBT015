# This program checks that all required packages are installed. If a package isn't installed the program will 
# install it locally. Since dependencies are tricky among the bioconductor packages the program will exit if 
# bioconductor package isn't installed (and won't try to install it).

if(!require(Rsamtools)){
  message("Rsamtool isn't installed, program will exit")
  quit(status=1)
}

if(!require(GenomicAlignments)){
  message("GenomicAlignments isn't installed, program will exit")
  quit(status=2)
}

if(!require(GenomicFeatures)){
  write("GenomicFeatures isn't installed, program will exit")
  quit(status=3)
}

if(!require(DESeq2)){
  message("Rsamtool isn't installed, program will exit")
  quit(status=4)
}

if(!require(pheatmap)){
  print("pheatmap isn't installed, program will exit:")
  quit(status=5)	
}

if(!require(xtable)){
  print("pheatmap isn't installed, program will exit:")
  quit(status=6)	
}

if(!require(RColorBrewer)){
  quit(status=7)	
  print("RColorBrewer isn't installed, program will exit:")
}

if(!require(ggplot2)){
  print("ggplot2 isn't installed, program will exit")
  quit(status=8)	
}

if(!require(PoiClaClu)){
  print("PoiClaClu isn't installed, program will exit")
  quit(status=9)	
}

message("All packages are installed")

quit(status = 0)

