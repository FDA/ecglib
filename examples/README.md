# T-wave delineator example
This program delineates Tpeak and Tend in a physionet record (.hea, .dat) using ecglib t-wave delineator.

## Requirements
- Median/representative beat ECG signal including vector magnitude lead sampled at 1000 Hz.
- QRS onset, R peak and QRS offset annotations in the median beat.
- Average RR interval in ms for the 10 second ECG strip from which the median beat was derived.

## Compile
g++ -std=c++11 -o twavedelineatorphysionet twavedelineatorphysionet.cpp -lecglib-core -lecglib-delineators-twave -lboost_program_options -lboost_filesystem -lboost_system -lwfdb

## Command line options

- record    : ECG record in Physionet or Physionet format ECG file
- qon       : QRS onset in ms
- rpeak     : R peak in ms
- qoff      : QRS offset in ms
- rr        : Mean RR in ms


To see command line options and their defaults type:

```
$ ./twavedelineatorphysionet --help
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
  --vcgmag2file arg (=0)                Export vector magnitude lead to text 
                                        file
```

## Output

Running with default parameters in default record should produce the following output:

```
RR,QON,RPEAK,QOFF,TPEAK,TPPEAK,TEND
808,297,350,392,601,-1,670
```

## LICENSE

twavedelineatorphysionet example code is in the public domain within the United States, and copyright and related rights in the work worldwide are waived through the CC0 1.0 Universal Public Domain Dedication. This example is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See DISCLAIMER section and https://www.gnu.org/licenses/gpl-faq.html for more details.
 
## DISCLAIMER

This software and documentation were developed by the authors in their capacities as  Oak Ridge Institute for Science and Education (ORISE) research fellows at the U.S. Food and Drug Administration (FDA).

FDA assumes no responsibility whatsoever for use by other parties of the Software, its source code, documentation or compiled executables, and makes no guarantees, expressed or implied, about its quality, reliability, or any other characteristic.  Further, FDA makes no representations that the use of the Software will not infringe any patent or proprietary rights of third parties.   The use of this code in no way implies endorsement by the FDA or confers any advantage in regulatory decisions.
 
