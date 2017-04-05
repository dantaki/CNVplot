# CNVplot

Plot CNV data with a genome viewer in R.

![alt text](CNVplot.jpg "de novo deletion in an ASD proband")

>[<sub>Brandler, W. M., Antaki, D., & Gujral, M. (2015). Frequency and complexity of de novo structural mutation in autism. Am J Hum Genet. doi:10.1101/030270</sub>](https://www.ncbi.nlm.nih.gov/pubmed/27018473)

## Installation

Install in R with [devtools](https://github.com/hadley/devtools)

```
install.packages("devtools")
library(devtools)
install_github("dantaki/CNVplot")
```

## Usage

```
CNVplot(df,Start,End,copyNumber,genome,title,yLabel)
```
### Arguments
  
ARG | Description 
--- | ---- 
df | Input data frame: **See Inputs for formatting requirements**
Start | Start position of the CNV
End | End position of the CNV
copyNumber | Copy number state. Deletion arguments: DEL,del,<2. Duplication arguments: DUP,dup,>2
genome | Reference genome build either hg19 or hg38
title | Title of the plot
yLabel | Y label

## Inputs

* Single individual

CHROM | POSITION | SAMPLE_DATA
--- | --- | --- 
chr8 | 1000 | -0.56
chr8 | 1023 | 0.12

* Trio

CHROM | POSITION | PROBAND | MOTHER | FATHER
------| -------- | ------- | ------ | ------
chr8 | 1000 | -0.56 | 0.13 | -0.43 
chr8 | 1023 | 0.12 | -0.63 | 0.26

## Author:

* Danny Antaki
  * dantaki@ucsd.edu

### Acknowledgements:

* William Brandler
* Jonathan Sebat
   * Sebat Lab http://sebatlab.ucsd.edu

## License 
MIT License

Copyright (c) 2016 Danny Antaki

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
