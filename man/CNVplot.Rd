\name{CNVplot}
\alias{CNVplot}
\title{CNVplot}
\usage{
CNVplot(df,Start,End,copyNumber,genome,title,yLabel)
}
\description{
Plots CNV data (eg: LogRRatio, Coverage) with a genomic viewer.
}
\arguments{
  \item{df}{Imput data frame: columns must be ordered <chromosome, position, proband, mother, father>. For plots with just one sample only pass <chromosome, position, sample_data>}

  \item{Start}{Start position of the CNV}
  \item{End}{End position of the CNV}
  \item{copyNumber}{Copy number state of the variant. Deletion arguments include :DEL,del,<2. Duplication arguments inlcude: DUP,dup,>2}
  \item{genome}{Reference genome build either hg19 or hg38}
  \item{title}{Title of the plot}
  \item{yLabel}{Y label}
}
\examples{
#
# each df contains one locus
#df <- read.table("myTrio_16p11.2_deletion.txt")
# ncol(df)
# 5
# colnames(df)
# 'chr', 'position', 'proband', 'mother', 'father'
# MAINTAIN THIS ORDER
# deletion in a trio,hg19
#CNVplot(df, "16p11.2 duplication example",29673954,30198600,1,"hg19","Coverage")


#df <- read.table("NA12878_chr1_dup.txt")
#ncol(df)
#> 3
# colnames(df)
# 'chr', 'position', 'sample'
# MAINTAIN THIS ORDER
# duplication in one sample, hg38
#CNVplot(df, "chr1:148895234-148895234 1|1 DUP",148895234,148895234,"DUP","hg38","Coverage")
}
