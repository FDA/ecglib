# File:     generatereport.R
# Author:   Jose Vicente <Jose.VicenteRuiz@fda.hhs.gov
# Date:     May 2017
# Version:  1.0.0
#
# **Disclaimer**
#   
# This code does not necessarily reflect any position of the Government or the Food and Drug Administration.
# 
# This software and documentation (the "Software") were developed at the Food and Drug Administration (FDA) by employees of the Federal Government in the course of their official duties. Pursuant to Title 17, Section 105 of the United States Code, this work is not subject to copyright protection and is in the public domain. Permission is hereby granted, free of charge, to any person obtaining a copy of the Software, to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, or sell copies of the Software or derivatives, and to permit persons to whom the Software is furnished to do so. FDA assumes no responsibility whatsoever for use by other parties of the Software, its source code, documentation or compiled executables, and makes no guarantees, expressed or implied, about its quality, reliability, or any other characteristic. Further, use of this code in no way implies endorsement by the FDA or confers any advantage in regulatory decisions. Although this software can be redistributed and/or modified freely, we ask that any derivative works bear some notice that they are derived from it, and any modified versions bear some notice that they have been modified.
#


# The line below was included to avoid the error message "error pandoc version 1.12.3 or higher is required" when running knitr outside rstudio in linux. You may need to comment out the ine below depending on your system's settings.
Sys.setenv(RSTUDIO_PANDOC="/usr/lib/rstudio/bin/pandoc")

#Set refence and annotations to compare from:
#  * reference: annotations from physionet atr files
#  * algorithm: annotations from running the code in github
theref <- 'reference'
themethod <- 'algorithm'
# Generate html report from Rmd file for output from algorithm using filter
filterwasenabled <- 1
rmarkdown::render("FDAStudy1Comparison.Rmd",output_file=paste0('FDAStudy1Comparison-',theref,'-vs-',themethod,'-Filtered.html'))
# Store measures in ISCE's 2017 J-Tpeak initiative csv format for ER analysis
erdataset <- allintervalswithtime %>% filter(annotator=='algorithm' & param=='JTPEAK') %>% 
  mutate(PARTICIPANT='Filtered',BS_FLAG=2) %>% rename(ECGID=EGREFID,JT_VM=value) %>% select(PARTICIPANT,BS_FLAG,ECGID,RR,JT_VM)
dir.create('submissions',showWarnings = F)
write.table(erdataset,'submissions/filtered.csv',row.names = F,sep=',')

# Generate html report from Rmd file for output from algorithm with no filter
filterwasenabled <- 0
rmarkdown::render("FDAStudy1Comparison.Rmd",output_file=paste0('FDAStudy1Comparison-',theref,'-vs-',themethod,'-NO-Filtered.html'))
# Store measures in ISCE's 2017 J-Tpeak initiative csv format for ER analysis
erdataset <- allintervalswithtime %>% filter(annotator=='algorithm' & param=='JTPEAK') %>% 
  mutate(PARTICIPANT='NOfiltered',BS_FLAG=2) %>% rename(ECGID=EGREFID,JT_VM=value) %>% select(PARTICIPANT,BS_FLAG,ECGID,RR,JT_VM)
dir.create('submissions',showWarnings = F)
write.table(erdataset,'submissions/nofiltered.csv',row.names = F,sep=',')

# Generate ER analysis reports
participant <- 'Filtered'
biomarker <- 'JT_VMc'
rmarkdown::render("ISCE-JTpeak-Time.and.ER.Rmd",output_file='ISCE-JTpeak-Time.and.ER-Filtered.html')
participant <- 'NOfiltered'
biomarker <- 'JT_VMc'
rmarkdown::render("ISCE-JTpeak-Time.and.ER.Rmd",output_file='ISCE-JTpeak-Time.and.ER-NO-Filtered.html')

#Generate Filter effects report
rmarkdown::render("validate.Rmd",output_file='validate.html')

#Generate Filter effects report
examplerecord <- 'ecgrdvq/medians/1005/643ee40b-9432-4dbc-9e67-a191c14c4843'
rmarkdown::render("FilterEffect.Rmd",output_file='FilterEffect.html')
