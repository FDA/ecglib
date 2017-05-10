# EXAMPLES

We have updated the example and the contents of this folder as follows:

* **filter**: this new folder contains several files implementating a Butterworth filter. We are providing this implementation so those who do not have a C++ signal processing library available can use this filter to pre-filter the ECG before calling the T-wave delineator. See line 60 in [twaveDelineator.cpp](../ecglib/src/delineators/twave/ecglib/delineator/twave/twaveDelineator.cpp) for more information about filtering requirements of the T-wave delineator.
* **twavedelineatorphysionet.cpp**: we have updated the example so it can use the *filter* provided above in case that a signal processing library is not installed. See twavedelineatorphysionet.cpp documentation for more information.
* **getdbannotations.cpp**: this program retrieves the annotations from a list of physionet files and write them in the standard output. See getdbannotations.cpp documentation for more information
* **delineateall.sh**: this script parses the annotations list output produced by *getdbannotations* together with the study clinical data file (SCR-002.Clinical.Data.csv) and calls *twavedelineatorphysionet* with and without the filtering enabled for each median ECG record.
* **FDAStudy1Comparison.Rmd**: R script that compares two annotations datasets and produces a report using Markdown syntax.
* **ISCE-JTpeak-Time.and.ER.Rmd**: R script that computes time-profiles and exposure-response plots for the J-Tpeak interval for each treatment in FDA Study 1 using the specified annotations dataset. Exprosure-response plot in the Markdown report is similar to Figure 3.K in [Vicente et. al JAHA 2015](https://doi.org/10.1161/JAHA.114.001615).
* **FilterEffect.Rmd**: R script that run the example with and without filter enabled and produces plots comparing signal sent to the T-wave delineator and delineation results.
* **validation.csv**:  File containing output annotations from T-wave delineator. This file can be used to validate the results obtained when running the T-wave delineator example for all ECGs.
* **validate.Rmd**: R script that compartes results obtained when delineating all ECGs vs. the validation.csv file above.
* **generatereport.R**: R script that generates Markdown reports for filtered and non-filtered annotation sets produced by *delineate.sh* script above. This script calls all .Rmd scripts above.

## Delineating all ECGS

Section *T-wave delineator example* below describes how to build and test the T-wave delineator example with one individual ECG record. File [HOWTO.md](HOWTO.md) describes how to delineate all ECG records in [FDA Study 1 database](http://physionet.org/physiobank/database/ecgrdvq) as well as how to produce statistical analysis reports using the source code and scripts enumerated above.
## T-wave delineator example
This program delineates Tpeak and Tend in a physionet record (.hea, .dat) using ecglib T-wave delineator.

### Requirements
- Median/representative beat ECG signal including vector magnitude lead sampled at 1000 Hz.
- QRS onset, R peak and QRS offset annotations in the median beat.
- Average RR interval in ms for the 10 second ECG strip from which the median beat was derived.

### Compile
To compile the T-wave delineator example:

```
g++ -std=c++11 -o twavedelineatorphysionet twavedelineatorphysionet.cpp filters/butterworth/butter.cpp -lecglib-core -lecglib-delineators-twave -lboost_program_options -lboost_filesystem -lboost_system -lboost_date_time -lwfdb
```

### Command line options

To see command line options and their defaults type:

```
./twavedelineatorphysionet --help
```

Current options and defaults

```
Usage: twavedelineatorphysionet [options]

T-wave delineator command line options:
  -h [ --help ]                         print this help
  --record arg (=ecgrdvq/medians/1001/00ed2097-cd14-4f03-ab33-853da5be5550)
                                        Physionet record
  --qon arg (=297)                      QRS onset in ms
  --rpeak arg (=350)                    R peak in ms
  --qoff arg (=392)                     QRS offset in ms
  --rr arg (=808)                       mean RR interval in ms
  --printheader arg (=1)                Flag to print header to standard output
  --filterecg arg (=1)                  If true then filter the ECG with a 5th 
                                        order butterworth filter before calling
                                        the T-wave delineator
  --vcgmag2file arg (=0)                Export vector magnitude lead delineated
                                        by the T-wave delineator to text file
```

### Output

Running with default parameters in default record should produce the following output:

```
RECORD,ERROR,FILTER,RR,QON,RPEAK,QOFF,TPEAK,TPPEAK,TEND
ecgrdvq/medians/1001/00ed2097-cd14-4f03-ab33-853da5be5550,,1,808,297,350,392,600,-1,670
```

### LICENSE

twavedelineatorphysionet example code is in the public domain within the United States, and copyright and related rights in the work worldwide are waived through the CC0 1.0 Universal Public Domain Dedication. This example is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See DISCLAIMER section below, https://github.com/FDA/ecglib/, https://creativecommons.org/publicdomain/zero/1.0/ and https://www.gnu.org/licenses/gpl-faq.html for more details.

### DISCLAIMER

This software and documentation were developed by the authors in their capacities as  Oak Ridge Institute for Science and Education (ORISE) research fellows at the U.S. Food and Drug Administration (FDA).

FDA assumes no responsibility whatsoever for use by other parties of the Software, its source code, documentation or compiled executables, and makes no guarantees, expressed or implied, about its quality, reliability, or any other characteristic.  Further, FDA makes no representations that the use of the Software will not infringe any patent or proprietary rights of third parties.   The use of this code in no way implies endorsement by the FDA or confers any advantage in regulatory decisions.
 
