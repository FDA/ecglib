/**
 * @file delineators/twave/ecglib/delineator/twave/twaveDelineator.cpp
 * @author Meisam Hosseini <meisam.hosseini@fda.hhs.gov>
 * @author Jose Vicente <jose.vicenteruiz@fda.hhs.gov>
 * @author Lars Johannesen <lars.johannesen@fda.hhs.gov>
 * @author Dustin C McAfee <dustin.mcafee@fda.hhs.gov
 *
 * @version 1.0.0
 *
 * @section LICENSE
 * ecglib is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License (GNU LGPL3), or (at your option) any later version.
 * ecglib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
 * You should have received a copy of the GNU Lesser General Public License along with this program.  If not, see http://www.gnu.org/licenses/.
 *
 * @section DISCLAIMER
 * ecglib software and documentation were developed by the authors in their capacities as  Oak Ridge Institute for Science and Education (ORISE) research fellows at the U.S. Food and Drug Administration (FDA). .
 * FDA assumes no responsibility whatsoever for use by other parties of the Software, its source code, documentation or compiled executables, and makes no guarantees, expressed or implied, about its quality, reliability, or any other characteristic.  Further, FDA makes no representations that the use of the Software will not infringe any patent or proprietary rights of third parties.   The use of this code in no way implies endorsement by the FDA or confers any advantage in regulatory decisions.
 *
 * @section DESCRIPTION
 * General input/output data to/form twaveDelineator
 */

#include <functional> // make_tuple

#include "delineate.hpp"

#ifdef ECGLIB_PREPROCESSORS
#include <ecglib/preprocessors/filters.hpp>
#endif

#include <ecglib/delineator/twave/twaveDelineator.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/any.hpp>

namespace ecglib {
    // prepartion of thresholds of classification rules based on decision tree
	// Input is an '_' separated string of thresholds and output is vector of thresholds
	std::vector<std::vector<double> > featursThresholdPreparation(std::string thresholds); // function prototype

	// main entrance into twaveDelineator for calculating twave annotations
	std::tuple<pointmap, ecglib::twaveDelineate::annotation> twaveDelineators(const ecglib::ecgdata &e, const ecglib::pointmap &pmin, const ecglib::twaveDelineator_config &cfg) {
		if(e.fs() != 1000){ // check the valid frequency
			std::string line = std::string("frequency should be 1000Hz");
			std::cerr << line;
			throw std::logic_error(line);
		}

		ecglib::ecgdata ecg(e);		// internal ecg variable
		ecglib::pointmap pm(pmin);	// internal annotation variable

		int vcgIndex (0);		// index of VCG inside ecgdata
		arma::vec filt = zeros<vec>(1);	// filter instantiation


		vcgIndex = ecg.leadnum(ecglead::VCGMAG); // index of VCG

		/* step 01: filter input ecg */
		filt(0) = cfg.get<double>("filterHighCutoff"); // high cutoff 25 Hz
#ifdef ECGLIB_PREPROCESSORS
        // Preprocessing and filtering methods are not released in version 1.0.0 of ecglib, but ecg should be filter as follows
		ecglib::filter filterData(cfg.get<int>("filterOrder"), filt, false); // 5th order, with cutoff 25 HZ and not a stop filter, i.e. a lowpass filter
		filterData(ecg);
#endif
		/* step 02: preparation of Twave range */
			// Determine seed points:
			// Strategy 1: Grab globals
			int seedoff = -1;
			std::vector<annotation> locs;
			get_annotations(pmin, GLOBAL_LEAD, annotation_type::QOFF, locs);

			if(locs.size() == 1) {	
				seedoff = locs[0].location();
			}
			locs.clear();
			// Strategy 2: Average across leads
			if(seedoff == -1) {
				get_annotations(pmin, annotation_type::QOFF, locs);
				vec offs = zeros<vec>(locs.size());
				std::copy(locs.begin(),locs.end(),offs.begin());
				seedoff = mean(offs);
				locs.clear();
			}

			// Determine rpeak for re-adjusting toff
			// Strategy 1: Grab globals
			int rpeak = -1;
			get_annotations(pmin, GLOBAL_LEAD, annotation_type::RPEAK, locs);

			if(locs.size() == 1) {	
				rpeak = locs[0].location();
			}
			locs.clear();
			// Strategy 2: Average across leads
			if(rpeak == -1) {
				get_annotations(pmin, annotation_type::RPEAK, locs);
				vec peaks = zeros<vec>(locs.size());
				std::copy(locs.begin(),locs.end(),peaks.begin());
				rpeak = mean(peaks);
				locs.clear();
				if(e.hasproperty("precut")) { // if there is not a rpeak, the rpeak will be precut
					rpeak = ecg.nsamples() - boost::any_cast<double>((e.getproperty("precut")).value);
				}else{ // if there is not a rpeak or precut, the rpeak will be 250 which is just an arbitary value
					rpeak = 250;
				}
			}

		/* step 03: calculates twave boundries [Qoff+a	Qoff+b] */
		double rr = 0;
		if(e.hasproperty("meanrr")) {
			rr = boost::any_cast<double>((e.getproperty("meanrr")).value); // mean value of rr based on property
		}else if(e.hasproperty("precut")) { // if there is not meanrr, the length of rr will be length of ecg - precut
			rr = ecg.nsamples() - boost::any_cast<double>((e.getproperty("precut")).value);
		}else{ // if there is not meanrr or precut, the approximate length of rr will be 80/100 of length of ecg
			rr = (80./100.*ecg.nsamples());
		}

		int pointStart = seedoff + 25; // 25 uses for avoiding j-point in calculations
		int pointEnd  = pointStart + (rr*cfg.get<double>("approximateRangeOfTsegment")); // for testing purpose 'pointEnd = pointStart + 300' got used
		if (pointEnd >= static_cast<int>(ecg.nsamples())) pointEnd = ecg.nsamples()-1;
		arma::vec twave = ecg.lead(vcgIndex, pointStart, pointEnd);

		/* step 04: call twave annotators functions*/
		ecglib::twaveDelineate::delineate deli;
		std::vector<std::vector<double> > featursThreshold = featursThresholdPreparation(cfg.get<std::string>("featursThreshold")); // threshoulds of classification rules based on decision tree
		ecglib::twaveDelineate::annotation anns = deli.delineator(twave.t(), pointStart, featursThreshold, cfg.get<int>("candidateFinder") , cfg.get<double>("deltaStepSlope"), cfg.get<int>("looseWindow"), cfg.get<int>("minPoints"), cfg.get<double>("deltaAmplitude"), cfg.get<double>("minVoltageMainPeak"), cfg.get<double>("percentMainePeak"), cfg.get<double>("minVoltage"), cfg.get<double>("percentPeak"), cfg.get<double>("maxDelatAplitudeNotches"), cfg.get<double>("minAmplitudeFlatness"), cfg.get<double>("minValidAmplitudePeak"), cfg.get<double>("measurable"));
		twave.clear(); // parameters of twave annotators

		/* step 05: re-adjusts toff place */
		arma::vec orignwave = ecg.lead(vcgIndex, 0, ecg.nsamples()-1);
		double toff_new = deli.readjustToff(orignwave.t(), anns, rr, rpeak);
		orignwave.clear();

		/* step 06: propagate the output delineators */

		// clear old annotations of twave
		locs.clear();
		get_annotations(pmin, vcgIndex, annotation_type::TON, locs);
		get_annotations(pmin, vcgIndex, annotation_type::TOFF, locs);
		get_annotations(pmin, vcgIndex, annotation_type::TPEAK, locs);
		get_annotations(pmin, vcgIndex, annotation_type::TPPEAK, locs);
		for(std::size_t i = 0; i < locs.size(); ++i) {
			pm[vcgIndex].erase(locs[i].location());
		}

		if (anns.peak.size() > 0) {
			// new toff
			double toff = anns.off;
			anns.rulesHit["toff_maxslope"] = toff;

			// re-assign toff based on value of re-adjused toff
			if (toff_new == -1) toff = 0; // calculation of toff had problem
			else toff = toff_new;

			// makes sure that qoff is not overwritten by ton
			locs.clear();
			get_annotations(pmin, vcgIndex, annotation_type::QOFF, locs);
			int qofftmp = -1;
			if(locs.size() == 1) {
				qofftmp = locs[0].location();
			}
			if (anns.on == qofftmp && anns.on > 0)
				pm[vcgIndex][(anns.on)+1] = annotation((anns.on)+1, annotation_type::TON, vcgIndex);
			else if (anns.on > qofftmp && anns.on > 0) // add ton when it is greater than qoff
				pm[vcgIndex][anns.on] = annotation(anns.on, annotation_type::TON, vcgIndex);

			if (toff < (pointStart+rr*cfg.get<double>("approximateBoundaryOfToff")) && toff != 0) { // dont expose Toff if it's too far or shorter than tpeak
				pm[vcgIndex][toff] = annotation(toff, annotation_type::TOFF, vcgIndex);
			}

			pm[vcgIndex][anns.peak[0]] = annotation(anns.peak[0], annotation_type::TPEAK, vcgIndex);				
			if (anns.peak.size() > 1) {
				pm[vcgIndex][anns.peak[1]] = annotation(anns.peak[1], annotation_type::TPPEAK, vcgIndex);
			}
			anns.rulesHit["hasDelineators"] = 1; // twaveDelineator has anns
		}
		else {
			anns.rulesHit["hasDelineators"] = 0; // twaveDelineator has not any anns
		}

		return std::make_tuple(pm,anns);
	}

	// Default config value of twaveDelineator's parameters
	void twaveDelineator_config::defaults() {
		add("filterHighCutoff",property(Type::Double,25.,"high cutoff a butterworth filter in Hz for filtering input ecg"));
		add("filterOrder",property(Type::Int,5,"order of butterworth filter for filtering input ecg"));
		add("candidateFinder",property(Type::Int,1,"finds candidates based on moving zero crossing line (1) or first/second derivative (2) functions"));
		add("featursThreshold",property(Type::String,std::string{"20,0,0_10,10,1.5_0,0,1.7"},"thresholds of extracted rules by Decision-Tree for slur classifier"));
		add("deltaStepSlope",property(Type::Double,10.,"interval value for calculating the number of moving zero crossing lines"));
		add("looseWindow",property(Type::Int,10,"min points of a valid candidate"));
		add("minPoints",property(Type::Int,10,"min points that make a candidate"));
		add("deltaAmplitude",property(Type::Double,5.,"delta amplitude of points that make a peak of candidate"));
		add("minVoltageMainPeak",property(Type::Double,150.,"minimum acceptable voltage of main peak"));
		add("percentMainePeak",property(Type::Double,80./100.,"percentage of main peak for evaluating the other candidates"));
		add("minVoltage",property(Type::Double,100.,"minimum acceptable voltage"));
		add("percentPeak",property(Type::Double,30./100.,"percentage of main peak for evaluating the other peaks"));
		add("maxDelatAplitudeNotches",property(Type::Double,50.," acceptable delta amplitudes of two peaks"));
		add("minAmplitudeFlatness",property(Type::Double,7.,"delta amplitudes of flatness "));
		add("minValidAmplitudePeak",property(Type::Double,7.,"threshold of peak candidate that declares small angle"));
		add("approximateRangeOfTsegment",property(Type::Double,40./100.,"aproximate range of Tsegment based on RR percentage"));
		add("approximateBoundaryOfToff",property(Type::Double,75./100.,"aproximate boundry of Toff based on RR percentage"));
		add("measurable",property(Type::Double,100.,"min threshsold of Tpeak amplitude as a measurable ecg"));
	}

    // prepartion of thresholds of classification rules based on decision tree
	// Input is an '_' separated string of thresholds and output is vector of thresholds
	std::vector<std::vector<double> > featursThresholdPreparation(std::string thresholds) {
		std::vector<std::vector<double> > featursThreshold; 	// each row contains related rules
		std::vector<std::string> rules;
		boost::split(rules, thresholds, boost::is_any_of("_")); // converto to related rules
		for (std::size_t i = 0; i < rules.size(); ++i) {
			std::vector<std::string> rule;
			std::vector<double> ruleValue;
			boost::split(rule, rules[i], boost::is_any_of(",")); // retrieve rules
		    	std::transform(rule.begin(), rule.end(), std::back_inserter(ruleValue), [](const std::string& str) { return std::stod(str); });
			featursThreshold.push_back(ruleValue);
		}
		return featursThreshold;
	}
}

