/**
 * @file twavedelineatorphysionet.cpp
 * @author Jose Vicente <jose.vicenteruiz@fda.hhs.gov>
 * @author Dustin C McAfee <dustin.mcafee@fda.hhs.gov
 * @author Meisam Hosseini <meisam.hosseini@fda.hhs.gov>
 * @author Lars Johannesen <lars.johannesen@fda.hhs.gov>
 *
 * @version 1.0.0
 *
 * @section LICENSE
 * twavedelineatorphysionet example code is in the public domain within the United States, and copyright and related rights in the work worldwide are waived through the CC0 1.0 Universal Public Domain Dedication. This example is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See DISCLAIMER section and https://www.gnu.org/licenses/gpl-faq.html for more details.
 *
 * @section DISCLAIMER
 * This software and documentation were developed by the authors in their capacities as  Oak Ridge Institute for Science and Education (ORISE) research fellows at the U.S. Food and Drug Administration (FDA).
 * FDA assumes no responsibility whatsoever for use by other parties of the Software, its source code, documentation or compiled executables, and makes no guarantees, expressed or implied, about its quality, reliability, or any other characteristic.  Further, FDA makes no representations that the use of the Software will not infringe any patent or proprietary rights of third parties.   The use of this code in no way implies endorsement by the FDA or confers any advantage in regulatory decisions.
 *
 * @section DESCRIPTION
 * Example of applying ecglib's T-wave delineator.
 * This Program delineates Tpeak and Tend in a physionet record (.hea, .dat) using ecglib t-wave delineator.
 * Requirements:
 *      - Median/representative beat ECG signal including vector magnitude lead sampled at 1000 Hz.
 *      - QRS onset, R peak and QRS offset annotations in the median beat.
 *      - Average RR interval in ms for the 10 second ECG strip from which the median beat was derived.
 * Arguments:
 *      - record    : ECG record
 *      - qon       : QRS onset in ms
 *      - rpeak     : R peak in ms
 *      - qoff      : QRS offset in ms
 *      - rr        : Mean RR in ms
 *
 *
 */

#include <ecglib.hpp>
#include <ecglib/delineator/twave.hpp>

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
#include <stdexcept>
#include <vector>
#include <map>

//Include Boost
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>
#include <boost/bimap.hpp>
#include <boost/assign.hpp>



using namespace ecglib;
using namespace arma;
using namespace std;


/**
 * @brief leadmap for WFDB library
 *
 * @param siarray array of WFDB signal info
 * @param nsig number of signals
 */
ecgdata::leadmap getleadnames(WFDB_Siginfo *siarray, const int nsig);

/**
 * @brief Load physionet
 *
 * @param e ECGdata record
 * @param rec Record
 */
void load_physionet(ecgdata &e, const char *rec);

int main(int argc, char **argv) {
    try{
        boost::program_options::options_description generic("T-wave delineator command line options");
        generic.add_options()
                ("help,h", "print this help")
                ("record",boost::program_options::value<std::string>()->default_value("ecgrdvq/medians/1001/00ed2097-cd14-4f03-ab33-853da5be5550"),"Physionet record")
                ("qon",boost::program_options::value<int>()->default_value(297),"QRS onset in ms")
                ("rpeak",boost::program_options::value<int>()->default_value(350),"R peak in ms")
                ("qoff",boost::program_options::value<int>()->default_value(392),"QRS offset in ms")
                ("rr",boost::program_options::value<double>()->default_value(808),"mean RR interval in ms")
                ("vcgmag2file",boost::program_options::value<bool>()->default_value(false),"Export vector magnitude lead to text file");

        boost::program_options::variables_map vm;
        boost::program_options::store(boost::program_options::parse_command_line(argc, argv, generic), vm);
        boost::program_options::notify(vm);

        if(vm.count("help")) {
            std::cout << std::endl << "Usage: twavedelineatorphysionet [options]" << std::endl;
            std::cout << std::endl << generic << std::endl;
            return EXIT_SUCCESS;
        }

        // Get input parameters from command line options
        std::string record =  vm["record"].as<std::string>();
        int qon =  vm["qon"].as<int>();
        int rpeak =  vm["rpeak"].as<int>();
        int qoff =  vm["qoff"].as<int>();
        double rr =  vm["rr"].as<double>();
        bool vcgmag2file =  vm["vcgmag2file"].as<bool>();

        //ecg variable
        ecgdata ecg;

        // Load ECG record
        load_physionet(ecg, record.c_str());

        // Rescale the ECG signal from mV to uV: Physionet files are stored in mV, but t-wave delineator expects the ECG to be in uV
        arma::mat data = ecg.data()*1000.0;
        ecg.data() = data;

        int vcgidx = ecg.leadnum(ecglead::VCGMAG);

        if (vcgmag2file){
            arma::mat ecgsignal = ecg.lead(vcgidx);
            ecgsignal.save("vcgmag.txt", arma_ascii);
        }

        // Check ECG sampling frequency is 1000 Hz
        if(ecg.fs() == 1000){
            //Create qon, rpeak and qoff annotations in the ECG
            pointmap pm = ecg.pointsmap();
            pm[GLOBAL_LEAD][qon] = annotation(qon, annotation_type::QON, GLOBAL_LEAD);
            pm[GLOBAL_LEAD][rpeak] = annotation(rpeak, annotation_type::RPEAK, GLOBAL_LEAD);
            pm[GLOBAL_LEAD][qoff] = annotation(qoff, annotation_type::QOFF, GLOBAL_LEAD);

            //Assign mean RR as a property of ecgdata
            const Type typ = ecglib::Type::Double;
            const ecglib::property prop(typ, rr);
            ecg.setproperty("meanrr", prop);

            int outtpeak = -1;
            int outtpeakp = -1;
            int outtend = -1;

            // T-wave delineation
            try{
                // Use T-wave delineator default settings
                twaveDelineator_config tcfg;
                // Delineate T-wave in the ecg
                std::tuple<pointmap, ecglib::twaveDelineate::annotation> res = twaveDelineators(ecg, pm, tcfg);

                // Retrieve post-delineation annotation locations
                // Get T-wave annotations from the vector magnitude lead (VCGMAG)
                vector<annotation> locs;
                pointmap pmout = get<0>(res);
                // Get Tpeak
                get_annotations(pmout, vcgidx, annotation_type::TPEAK, locs);
                if(locs.size() == 1) {
                    outtpeak = locs[0].location();
                }
                locs.clear();

                // Get Tend
                get_annotations(pmout, vcgidx, annotation_type::TOFF, locs);
                if(locs.size() == 1) {
                    outtend = locs[0].location();
                }
                locs.clear();

                // Get secondary Tpeak
                get_annotations(pmout, vcgidx, annotation_type::TPPEAK, locs);
                if(locs.size() == 1) {
                    outtpeakp = locs[0].location();
                }
                locs.clear();

                std::cout << "RR,QON,RPEAK,QOFF,TPEAK,TPPEAK,TEND" << std::endl;
                std::cout << rr << "," << qon << "," << rpeak << "," << qoff << "," << outtpeak << "," << outtpeakp << "," << outtend << std::endl;
            }catch(const std::exception &e){
                std::string line = std::string("Exception caught in t-wave Delineator: ") + e.what();
                std::cerr << std::endl << line << std::endl;
            }catch(const char* e){
                std::string line = std::string("Exception caught in t-wave Delineator: ") + e;
                std::cerr << std::endl << line << std::endl;
            }catch(...){
                std::string line = std::string("Unknown exception caught in t-wave Delineator");
                std::cerr << std::endl << line << std::endl;
            }

        } else {
            //else (frequency != 1kHz) report error
            std::cerr << "** ERROR: " << "Sampling frequency for " << record << " is at " << ecg.fs() << " but 1000 Hz is required." << std::endl;
        } // end if(ecg.fs() == 1000)

        std::cout << std::endl;
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

ecgdata::leadmap getleadnames(WFDB_Siginfo *siarray, const int nsig) {
        long nsamp = siarray[0].nsamp;
        ecgdata::leadmap lm;
        std::string nam2;

        for(int i = 0; i < nsig; ++i) {
                if(nsamp != siarray[i].nsamp){
                        std::string line = std::string("Mismatch nsamp");
                        std::cerr << line;
                        throw std::runtime_error(line);
                }

                std::string nam(siarray[i].desc);
                if(nam != std::string("ECG")) {
                        nam2 = nam;
                        std::transform(nam2.begin(), nam2.end(), nam2.begin(), ::tolower);
                        if(nam2 == std::string("vx") || nam2 == std::string("vy") || nam2 == std::string("vz")){
                                nam.erase(0,1);
                        }
                        boost::optional<ecglead> ll = ecglead::get_by_name(nam.c_str());
                        if(ll) {
                                ecgdata::leadmap::value_type lval(i,*ll);
                                lm.insert(lval);
                        }
                }
        }

        return lm;
}

void load_physionet(ecgdata &e, const char *rec) {
        char *rec2 = const_cast<char*>(rec);

        double fs = sampfreq(rec2);

        WFDB_Siginfo *siarray;

        int nsig = isigopen(rec2,NULL,0);

        siarray = (WFDB_Siginfo*)malloc(nsig*sizeof(WFDB_Siginfo));
        nsig = isigopen(rec2,siarray,nsig);

        if(nsig == -1) {
                std::string line = std::string("Could not open file at ") + std::string(rec);
                std::cerr << line << std::endl;
                throw std::logic_error(line);
        }
        long nsamp = siarray[0].nsamp;
        ecgdata::leadmap leadnames = getleadnames(siarray, nsig);
        ecglib::ecgdata er(nsamp, nsig);
        er.leadnames(leadnames);

        WFDB_Sample *vin;
        vin = (WFDB_Sample*)malloc(nsig*sizeof(WFDB_Sample));
        for (int i = 0; i < nsamp; i++) {
                getvec(vin);

                for(int j = 0; j < nsig; ++j) er(j,i) = aduphys(j,vin[j]);
        }

        er.fs(fs);

        free(siarray);

        e = er;
}
