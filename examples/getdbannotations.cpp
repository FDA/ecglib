/**
 * @file delineateallrecords.cpp
 * @author Jose Vicente <jose.vicenteruiz@fda.hhs.gov>
 * @author Dustin C McAfee <dustin.mcafee@fda.hhs.gov
 * @author Meisam Hosseini <meisam.hosseini@fda.hhs.gov>
 * @author Lars Johannesen <lars.johannesen@fda.hhs.gov>
 *
 * @version 1.0.0
 *
 * @section LICENSE
 * getdbannotations example code is in the public domain within the United States, and copyright and related rights in the work worldwide are waived through the CC0 1.0 Universal Public Domain Dedication. This example is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See DISCLAIMER section below, https://github.com/FDA/ecglib/, https://creativecommons.org/publicdomain/zero/1.0/ and https://www.gnu.org/licenses/gpl-faq.html for more details.
 *
 * @section DISCLAIMER
 * This software and documentation were developed by the authors in their capacities as  Oak Ridge Institute for Science and Education (ORISE) research fellows at the U.S. Food and Drug Administration (FDA).
 * FDA assumes no responsibility whatsoever for use by other parties of the Software, its source code, documentation or compiled executables, and makes no guarantees, expressed or implied, about its quality, reliability, or any other characteristic.  Further, FDA makes no representations that the use of the Software will not infringe any patent or proprietary rights of third parties.   The use of this code in no way implies endorsement by the FDA or confers any advantage in regulatory decisions.
 *
 * @section DESCRIPTION
 * This program retrieves the annotations from the list of records specified in index and writes them to the standard output.
 *
 * Arguments:
 *	- index		: file containing list of ecg RECORDS 
 */


//Include WFDB Library
#include <wfdb/wfdb.h>
#include <wfdb/ecgcodes.h>

#include <stdexcept>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <tuple>

#include <cmath>
#include <sstream>
#include <cstdlib>
#include <map>

//Include Boost
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>
#include <boost/bimap.hpp>
#include <boost/assign.hpp>


using namespace std;

/**
 * @brief Load physionet annotations
 *
 * @param pon		P wave onset location
 * @param qon		QRS complex onset location
 * @param qoff		QRS complex offset (J-point) location
 * @param tpeak		Location of the first peak of the T-wave
 * @param tpeak2	Location of second peak of the T-wave (set to -1 if not present)
 * @param toff		T-wave end location
 * @param rec 		Physionet's record to delineate
 */
void load_physionet_ann(int &pon, int &qon, int &qoff, int &tpeak, int &tpeak2, int &toff, const char *rec);

int main(int argc, char **argv) {
    try{
        boost::program_options::options_description generic("getdbannotations command line options");
        generic.add_options()
                ("help,h", "print this help")
                ("index",boost::program_options::value<std::string>()->default_value("allmedians.csv"),"List of records to get annotations from ecgrdvq Physionet RECORDS files");

        boost::program_options::variables_map vm;
        boost::program_options::store(boost::program_options::parse_command_line(argc, argv, generic), vm);
        boost::program_options::notify(vm);

        if(vm.count("help")) {
            std::cout << std::endl << "Usage: getdbannotations [options]" << std::endl;
            std::cout << std::endl << generic << std::endl;
            return EXIT_SUCCESS;
        }

        // Get input parameters from command line options
        std::string index =  vm["index"].as<std::string>();

	// Load index file
	std::ifstream ifindex(index);
	// Loop through all ECG records
	std::string record;
	int pon = -1;
	int qon = -1;
	int qoff = -1;
	int tpeak = -1;
	int tpeak2 = -1;
	int toff = -1;
	std::cout << "RECORD,EGREFID,PON,QON,QOFF,TPEAK,TPPEAK,TOFF" << std::endl;
	while (std::getline(ifindex, record)){
		std::string fullrecord = "ecgrdvq/" + record;
		// Load ECG record
        	load_physionet_ann(pon,qon,qoff,tpeak,tpeak2,toff,fullrecord.c_str());
		// Parse record name to get EGREFID string
		std::string delim = "/";
		size_t pos = 0;
		std::string egrefid = record;
		while ((pos = egrefid.find(delim)) != std::string::npos) {
			egrefid.erase(0, pos + delim.length());
		}
		// Call twavedelineator with parameters
		std::cout << fullrecord << "," << egrefid << "," << pon << "," << qon << "," << qoff << "," << tpeak << "," << tpeak2 << "," << toff << std::endl;
	}
        return EXIT_SUCCESS;
    }catch(const std::exception &e){
        std::string line = std::string("Exception caught: ") + e.what();
        std::cerr << std::endl << line << std::endl;
    }catch(const char* e){
        std::string line = std::string("Exception caught: ") + e;
        std::cerr << std::endl << line << std::endl;
    }catch(...){
	    std::string line = std::string("Unknown exception caught");
	    std::cerr << std::endl << line << std::endl;
    }
    return EXIT_FAILURE;
}

// ----------------------------
// Functions implementations
// ----------------------------

void load_physionet_ann(int &pon, int &qon, int &qoff, int &tpeak, int &tpeak2, int &toff, const char *rec) {
	char *rec2 = const_cast<char*>(rec);

	WFDB_Anninfo a;
	WFDB_Annotation annot;
	a.name = const_cast<char*>("atr"); a.stat = WFDB_READ;
	pon = -1;
	qon = -1;
	qoff = -1;
	tpeak = -1;
	tpeak2 = -1;
	toff = -1;
	int res;
	if (annopen(rec2, &a, 1) >= 0){
		// P onset
		res = getann(0, &annot);
		if (res == 0){
			if (annot.anntyp == WFON){
				if (annot.time < 300){
					// P onset
					pon = annot.time;
				} else {
					//These are two records with no P onset
					// QRS onset
					qon = annot.time;
				}
			}
		}
		if (qon < 0) {
			// QRS onset
			res = getann(0, &annot);
			if (res == 0){
				if (annot.anntyp == WFON){
					qon = annot.time;
				}
			}
		}
		// QRS offset (J-point)
		res = getann(0, &annot);
		if (res == 0){
			if (annot.anntyp == WFOFF){
				qoff = annot.time;
			}
		}
		// T wave
		res = getann(0, &annot);
		if (res == 0){
			// T peak
			if (annot.anntyp == TWAVE){
				tpeak = annot.time;
			}
			res = getann(0, &annot);
			if (res == 0){
				// Second T peak
				if (annot.anntyp == TWAVE){
					tpeak2 = annot.time;
					// T offset (tend)
					res = getann(0, &annot);
					if (res == 0){
						if (annot.anntyp == WFOFF) {
							// T offset (Tend)
							toff = annot.time;
						}
					}
				} else if (annot.anntyp == WFOFF) {
					// T offset (Tend)
					toff = annot.time;
				}
			}
		}
	}
}
