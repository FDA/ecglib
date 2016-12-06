/**
 * @file delineators/twave/ecglib/delineator/twave/delineate.cpp
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
 * Twave Delineator functionality
 */

#include "math.h" // for calculating sigmoid function
#include "delineate.hpp"

using namespace ecglib::twaveDelineate;


/// Main function for finding twave delineation based on twave's candidates containing peaks & slurs
annotation delineate::delineator(const arma::rowvec& twave, int pointStart, const std::vector<std::vector<double> >& featursThreshold,
				int candidateFinderFlag, double deltaStepSlope, int looseWindow, int minPoints, double deltaAmplitude, double minVoltageMainPeak, double percentMainePeak,
				double minVoltage, double percentPeak, double maxDelatAplitudeNotches, double minAmplitudeFlatness, double minValidAmplitudePeak, double measurableVoltage) {
	/*
	// This function takes a twave as an input and returns twave's annotators ton/ tpeak/ tppeak/ toff/ flatness/ skewness/ rotation
	// twave should be a filtered wave otherwise this function could not find proper annotators
	*/

	annotation anns;				// main output structure of annotators
	std::vector<candidate> candids; 		// main internal structure for keeping info of all candidates
	std::vector<int> badCandidates; 		// temporary container for removing candidates in different steps
	try {	
		arma::rowvec derivative = twave(arma::span(1,twave.n_elem-1)) - twave(arma::span(0,twave.n_elem-2)); // first derivative of twave

		/* step01: finds candidates based on 1) moving the zero crossing line on the first derivative 2) moving on the first derivative */
		arma::rowvec movedZeroCrossingAmplitutedPage;	// matrix of [different slopes, candidate Amplitude]
		arma::uvec candidatePeaksPosition;		// place of peaks (candidate has intersection with slope = 0)
		//double deltaStepSlope = 10; // should be greater than 1. This values is used as inverted value (1/deltaStepSlope).
				    // bigger number increases the time complexity and implies more accurate calculation of slopes based on moving zero crossing line.
				    // smaller number can reduce the number of candidates.
				    // 10 is about 5.7' degree in which removes slurs less than this range of slope.	
		if (candidateFinderFlag != 2)
			delineate::candidateFinder(twave, derivative, movedZeroCrossingAmplitutedPage, candidatePeaksPosition, deltaStepSlope); // finds the candidates
		else
			delineate::candidateFinder2DerivativeBased(twave, derivative, movedZeroCrossingAmplitutedPage, candidatePeaksPosition); // finds the candidates	

		/* step02: finds a range of each candidate for later processing */
		//int looseWindow = 10; // defines min points of a candidate: bigger size can marge candidates and smaller size can generate more candidates
		delineate::candidateRangeInfoFinder(derivative, movedZeroCrossingAmplitutedPage, candids, looseWindow); // find the candidate ranges containing start of rising slope of candidate & end of falling slope of candidates
		movedZeroCrossingAmplitutedPage.clear();

		/* step03 (rule pre_process): removes small candidates in term of number of points */
		//int minPoints = 10; // min points that make a candidate
		badCandidates = preProcessingRules::fewPointsCandidates(candids, minPoints);
		anns.rulesHit["fewPointsCandidates"] = badCandidates.size();
		delineate::cleanUpCandidates(candids, badCandidates);

		/* step04: labels candidates by slur(0) & peak(2) */
		delineate::labellingPeaks(candids, candidatePeaksPosition);
		candidatePeaksPosition.clear();

		/* step05: finds annotation of each candidate */
		//double deltaAmplitude = 5; 	/* pre defined threshold for calculating peak of each candidate
		//				   The candidat's peak can contain couple of points with highest amplitude <= deltaAmplitude */
		for (candidate& candid : candids) {
			delineateFinder::delineatorsInfo(twave, derivative, candid, deltaAmplitude); // finds rising slope/ peak/ falling slope/ skewness/ distortion & flatness of candidate
		}
		derivative.clear();

		/* *********		   ********* */
		/* * * * * 		    * * * *  */
		/*  * * * 		     * * * * */
		/* ***** post processing steps ***** */

		/* step06 (rule post_process): removes candidates when main peak has a low amplitude*/
		badCandidates = postProcessingRules::lowAmplitudeMainPeak(candids, minVoltageMainPeak, percentMainePeak);
		anns.rulesHit["lowAmplitudeMainPeak"] = badCandidates.size();
		delineate::cleanUpCandidates(candids, badCandidates);

		/* step07 (rule post_process): removes peaks with low amplitude based on max peak */
		badCandidates = postProcessingRules::lowAmplitudePeaks(candids, minVoltage, percentPeak);
		anns.rulesHit["lowAmplitudePeaks"] = badCandidates.size();
		delineate::cleanUpCandidates(candids, badCandidates);

		/* step08 (rule post_process): remove a peak candidate, if the amplitude differences between max peak and this peak is considerable */
		badCandidates = postProcessingRules::inconsistentPeaks(candids, maxDelatAplitudeNotches);
		anns.rulesHit["inconsistentPeaks"] = badCandidates.size();
		delineate::cleanUpCandidates(candids, badCandidates);

		/* step09: re-labeles the slur candidates
		 	   slur(0)[contains slur before or after slur & slur with inconsistent slope],
		           rising slur(1) [contains slur before a peak with rising slope], and
		           falling slur(-1) [contains slur after a peak with falling slope] */
		delineate::reLabellingSlurs(candids);

		/* step10 (rule post_process): removes unrelated slur(0)
			   caution: this step should not be commented or ignored. It has an effects on output preparation steps */
		badCandidates = postProcessingRules::unrelatedSlure(candids);
		anns.rulesHit["unrelatedSlure"] = badCandidates.size();
		delineate::cleanUpCandidates(candids, badCandidates);

		/* step11 (rule post_process): merges two consecutive candidates  for shaping a flat candidate based on amplitude */
		badCandidates = postProcessingRules::meargingCandidates(twave, candids, minAmplitudeFlatness);
		anns.rulesHit["meargingCandidates"] = badCandidates.size() / 2;
		delineate::cleanUpCandidates(candids, badCandidates);

		/* step12 (rule post_process): distinguishes between good slur and bad slur to clean up the slur's candidates (using extracted rules by decision-tree) */
		badCandidates = postProcessingRules::slurClassifier(candids, featursThreshold);
		anns.rulesHit["slurClassifier"] = badCandidates.size();
		cleanUpCandidates(candids, badCandidates);

		/* step13 (rule post_process): keeps max two peaks based on max amplitude
		          this step changes the label of peaks after second peak */
		badCandidates = postProcessingRules::keepJustTwoPeaks(candids);
		anns.rulesHit["keepJustTwoPeaks"] = badCandidates.size();
		for (std::size_t i = 0; i < badCandidates.size() && candids.size() > 0; ++i)
			candids[badCandidates[i]].set_label(candidateLabel::peakUnrelated);	// converts to unrelated peak
		eraseList(badCandidates);	/* do not clear peaks after second peak from list just simply do not use them for final output.
						   Still the slope of these peaks can change the place of on/off set
						   if do not remove them, they will assume as slur. if want to remove them, should remove dependent slurs with those peaks as well */

		/* step14 (rule post_process): converts a peak to slure when one of its angle is really small based on delta ampiltude of that angle with local minima */
		badCandidates = postProcessingRules::convertPeakToSlur(twave, candids, minValidAmplitudePeak);	
		anns.rulesHit["convertPeakToSlur"] = badCandidates.size();
		for (std::size_t i = 0; i < badCandidates.size(); ++i)
			candids[badCandidates[i]].set_label(candidateLabel::sluredPeak);	// converts to slured-peak
		eraseList(badCandidates);	/* do not clear peaks after second peak from list just simply do not use them for final output.
						   Still the slope of these peaks can change the place of on/off set
						   if do not remove them, they will assume as slur. if want to remove them, should remove dependent slurs with those peaks as well */

		/* step15 (rule post_process): just an informative flag to label a non-measurable signal (low main peak amplitude) */
		bool nonMeasurable = postProcessingRules::nonMeasurableSignal(candids, measurableVoltage);
		anns.rulesHit["non-measurable"] = (nonMeasurable? 1: 0);

		/* step16: output preparation */
		std::vector<int> indexPeakCandidates = postProcessingRules::peakCandidates(candids);
		if (indexPeakCandidates.size() > 0) {   // if this condition hits false: main peak has got removed based on current rules
			anns.on = std::round(pointStart - (candids[0].get_b0() / candids[0].get_a0())); // intersection between max slope of first left candidate with line amplitude  = 0
			anns.off = std::round(pointStart - (candids[candids.size()-1].get_b1() / candids[candids.size()-1].get_a1())); // intersection between min slope of last right candidate with line amplitude  = 0
			anns.lastCandidate = pointStart + candids[candids.size()-1].get_x();

			for (std::size_t i = 0; i < indexPeakCandidates.size(); ++i) {
				anns.peak.push_back(pointStart + candids[indexPeakCandidates[i]].get_x());
				anns.flatness.push_back(candids[indexPeakCandidates[i]].get_flatnessSamples());
				anns.distortion.push_back(candids[indexPeakCandidates[i]].get_distortion());
				anns.skewness.push_back(candids[indexPeakCandidates[i]].get_skewness());
			}
		}

	} catch(const std::exception &e){
		std::string line = std::string("Caught exception (Could not delineate ECG, error in t-wave delineation): ") + e.what();
		std::cerr << line;
		throw;
	} catch(const char* e){
		std::string line = std::string("Caught exception (Could not delineate ECG, error in t-wave delineation): ") + e;
		std::cerr << line;
		throw;
	} catch(...){
		std::string line = std::string("Unknown caught exception (Could not delineate ECG, error in t-wave delineation)");
		std::cerr << line;
		throw;
	}

	return anns;
}

/// finds all candidates of twave based on moving zero crossing line on a clean first derivative vectore
void delineate::candidateFinder(const arma::rowvec& wave, const arma::rowvec& derivative, arma::rowvec& movedZeroCrossingAmplitutedMax, arma::uvec& candidatePeaksPosition, double deltaStepSlope) {
	// moving the slope origin from max slop of first derivative to min slope of first derivative to find all candidates (peak/slur) and range of each
	// this function takes care of falt slopes (derivative = zero) as well

	double startingSlope = std::trunc(derivative.max() * deltaStepSlope) / deltaStepSlope;
	double endingSlope   = std::trunc(derivative.min() * deltaStepSlope) / deltaStepSlope;
	int numberSlope      = ((startingSlope - endingSlope ) * deltaStepSlope) + 1; // number of different slopes when moving the slope origin line

	int zeroSlopeIndex = static_cast<int>(startingSlope * deltaStepSlope); // index slope when it's equal to zero: for sanity check 'zeroSlopeIndex - 1' & 'zeroSlopeIndex + 1' should be checked too
	arma::mat movedZeroCrossingPage(numberSlope, derivative.n_elem -1, arma::fill::zeros); // keeps the intersection of moving slope line with derivative
	arma::mat peakCandidates(3, derivative.n_elem -1, arma::fill::zeros); // 3 cells: for 'zeroSlopeIndex' & 'zeroSlopeIndex - 1' & 'zeroSlopeIndex + 1'

	for (int j = 0; j < numberSlope; ++j) {
		double movingOriginZeroSlope = startingSlope + (j * (-1 / deltaStepSlope));
		bool signFlag = true;
        	arma::uvec indexFlatSlope;

		for (int i = 1; i < static_cast<int>(derivative.n_elem); ++i) {
			double deltaSlope = std::trunc((derivative(i) - movingOriginZeroSlope) * deltaStepSlope);
			int signDeltaSlope = (deltaSlope > 0) ? 1 : ((deltaSlope < 0) ? -1 : 0);
			double deltaSlopePrevious = std::trunc((derivative(i-1) - movingOriginZeroSlope) * deltaStepSlope);
			int signDeltaSlopePrevious = (deltaSlopePrevious > 0) ? 1 : ((deltaSlopePrevious < 0) ? -1 : 0);

			if (signDeltaSlope == 1) {
				signFlag = true;
				indexFlatSlope.clear();
				if ((i-1)-1 >= 0) indexFlatSlope = arma::find(movedZeroCrossingPage(j,arma::span(0,(i-1)-1)) == 2);  // index of flat slopes on derivative (temporary index)
				arma::urowvec jj(1); jj = j;
				if (!indexFlatSlope.is_empty()) movedZeroCrossingPage(jj,indexFlatSlope) = arma::zeros<arma::mat>(1,indexFlatSlope.n_elem); // reset index of flat slopes on derivative (temporary index)
			}
			else if (signDeltaSlope == 0) {
				if (signDeltaSlopePrevious == 1)
					movedZeroCrossingPage(j,i-1) = 2; //set index of flat slope on derivative (temporary index)
				else if (signDeltaSlopePrevious == 0) {
					if (signFlag == true)
						movedZeroCrossingPage(j,i-1) = 2; //set index of flat slopes on derivative (temporary index)
					else // if (signFlag == false)
						signFlag = false;
					}
				else if (signDeltaSlopePrevious == -1)
					signFlag = false;
			}
 			else if(signDeltaSlope == -1) {
				if (signDeltaSlopePrevious == 1) {
					movedZeroCrossingPage(j,i-1) = 1; // moved zero crossing line has intersection with derivative (one point of a candidate )
					int jz = (j == zeroSlopeIndex -1) ? 0 : (j == zeroSlopeIndex ? 1 : (j == zeroSlopeIndex +1 ? 2 : -1));
					if (jz != -1) peakCandidates(jz,i-1) = 1; // derivative has intersection with moved zero crossing line = 0 (slope == 0)
				}
				else if (signDeltaSlopePrevious == 0) {
					if (signFlag == true) {
						movedZeroCrossingPage(j,i-1) = 1;
						indexFlatSlope.clear();
						if ((i-1)-1 >= 0) indexFlatSlope = arma::find(movedZeroCrossingPage(j,arma::span(0,(i-1)-1)) == 2); // index of flat slopes on derivative (temporary index)
						arma::urowvec jj(1); jj = j;
						if (!indexFlatSlope.is_empty()) movedZeroCrossingPage(jj,indexFlatSlope) = arma::ones<arma::mat>(1,indexFlatSlope.n_elem);// moved zero crossing line has intersection with derivative (based on temporary index of flat slopes)

						int jz = (j == zeroSlopeIndex -1) ? 0 : (j == zeroSlopeIndex ? 1 : (j == zeroSlopeIndex +1 ? 2 : -1));
						if (jz != -1) {
							peakCandidates(jz,i-1) = 1; // moved zero crossing line has intersection with derivative (one point of a candidate )
							arma::urowvec jj(1); jj = jz;
							if (!indexFlatSlope.is_empty()) peakCandidates(jj,indexFlatSlope) = arma::ones<arma::mat>(1,indexFlatSlope.n_elem); // moved zero crossing line has intersection with derivative (based on temporary index of flat slopes)
						}
					}
					else { //if (signFlag == false)
						indexFlatSlope.clear();
						if ((i-1)-1 >= 0) indexFlatSlope = arma::find(movedZeroCrossingPage(j,arma::span(0,(i-1)-1)) == 2); // index of flat slopes on derivative (temporary index)
						arma::urowvec jj(1); jj = j;
						if (!indexFlatSlope.is_empty()) movedZeroCrossingPage(jj,indexFlatSlope) = arma::zeros<arma::mat>(1,indexFlatSlope.n_elem); // reset index of flat slopes on derivative (temporary index)
					}
				} // end of: if(signDeltaSlopePrevious == 0)
				else if (signDeltaSlopePrevious == -1) {
					// do nothing
				}
				signFlag = false;
				indexFlatSlope.clear();
				if ((i-1)-1 >= 0) indexFlatSlope = arma::find(movedZeroCrossingPage(j,arma::span(0,(i-1)-1)) == 2); // index of flat slopes on derivative (temporary index)
				arma::urowvec jj(1); jj = j;
				if (!indexFlatSlope.is_empty()) movedZeroCrossingPage(jj,indexFlatSlope) = arma::zeros<arma::mat>(1,indexFlatSlope.n_elem); // reset index of flat slopes on derivative (temporary index)
			} // end of: if(signDeltaSlope == -1)
		} // end of for: i

		indexFlatSlope = arma::find(movedZeroCrossingPage(j,arma::span::all) == 2);
		arma::urowvec jj(1); jj = j;
		if (!indexFlatSlope.is_empty()) movedZeroCrossingPage(jj,indexFlatSlope) = arma::zeros<arma::mat>(1,indexFlatSlope.n_elem); // reset index of flat slopes on derivative (temporary index)
	} // end of for: j

	movedZeroCrossingAmplitutedMax = arma::max(movedZeroCrossingPage,0) % wave(arma::span(1, wave.n_elem -2)); // vectore of intersection of moving zero crossing line * amplitude of each point
	arma::rowvec peakCandidatesMax = arma::max(peakCandidates,0);
	candidatePeaksPosition = arma::find(peakCandidatesMax > 0); // place of peaks (intersection with slope = 0 ) based on first derivative
}

/// finds all candidates of twave based on a clean derivative (first/second) vectore
void delineate::candidateFinder2DerivativeBased(const arma::rowvec& wave, const arma::rowvec& derivative, arma::rowvec& candidateRangeBasedOnDerivative, arma::uvec& candidatePeaksPosition) {
	// walking on a first derivative to find all candidates (peak/slur) and range of each due to second derivative
	// this function takes care of falt slopes (derivative = zero) as well
	double deltaStepSlope = 100;
	std::vector<double> precision {0, 0.1, 0.2};

	arma::mat walkOnDerivative(precision.size(), derivative.n_elem -1, arma::fill::zeros); 	// keeps the falling curve of derivative
	arma::mat peakCandidates(precision.size(), derivative.n_elem -1, arma::fill::zeros); 	// peak candidates based on derivative

	for (int j = 0; j < static_cast<int>(precision.size()); ++j) {
		bool signFlag = true;
		bool signFlagZero = true;
        	arma::uvec indexFlatSlope;

		for (int i = 1; i < static_cast<int>(derivative.n_elem); ++i) {
			double deltaSlopeCurrent = std::trunc((derivative(i) - precision[j]) * deltaStepSlope);
			int signDeltaSlopeCurrent = (deltaSlopeCurrent > 0) ? 1 : ((deltaSlopeCurrent < 0) ? -1 : 0);
			double deltaSlopePrevious = std::trunc((derivative(i-1) - precision[j]) * deltaStepSlope);
			int signDeltaSlopePrevious = (deltaSlopePrevious > 0) ? 1 : ((deltaSlopePrevious < 0) ? -1 : 0);
			double deltaSlope = deltaSlopeCurrent - deltaSlopePrevious;
			int signDeltaSlope = (deltaSlope > 0) ? 1 : ((deltaSlope < 0) ? -1 : 0);

			// candidates with local maxima (crossing the zero line of second derivative)
			if (signDeltaSlopeCurrent == 1) {
				signFlagZero = true;
				indexFlatSlope.clear();
				indexFlatSlope = arma::find(peakCandidates(j,arma::span::all) == 2);
				arma::urowvec jj(1); jj = j;
				if (!indexFlatSlope.is_empty()) peakCandidates(jj,indexFlatSlope) = arma::zeros<arma::mat>(0,indexFlatSlope.n_elem); // reset index of local maxima on derivative (temporary index)
			}
			else if (signDeltaSlopeCurrent == 0) {
				if (signFlagZero == true) {
					peakCandidates(j,i-1) = 2; // set temporary index for flat slope
				}
			}
			else if (signDeltaSlopeCurrent == -1) {
				if (signDeltaSlopePrevious == 1) {
					peakCandidates(j,i-1) = 1; // intersection with derivative
					indexFlatSlope.clear();
					indexFlatSlope = arma::find(peakCandidates(j,arma::span::all) == 2);
					arma::urowvec jj(1); jj = j;
					if (!indexFlatSlope.is_empty()) peakCandidates(jj,indexFlatSlope) = arma::ones<arma::mat>(0,indexFlatSlope.n_elem); // set index of local maxima on derivative
				}
				else if (signDeltaSlopePrevious == 0) {
				    if (signFlagZero == 1) {
					peakCandidates(j,i-1) = 1; // set index of local maxima on derivative
					indexFlatSlope.clear();
					indexFlatSlope = arma::find(peakCandidates(j,arma::span::all) == 2);
					arma::urowvec jj(1); jj = j;
					if (!indexFlatSlope.is_empty()) peakCandidates(jj,indexFlatSlope) = arma::ones<arma::mat>(1,indexFlatSlope.n_elem); // set index of local maxima on derivative derivative			       
				    }
				}
				else if (signDeltaSlopePrevious == -1) {
				    // do nothing
				}
				signFlagZero = false;
				indexFlatSlope.clear();
				indexFlatSlope = arma::find(peakCandidates(j,arma::span::all) == 2);
				arma::urowvec jj(1); jj = j;
				if (!indexFlatSlope.is_empty()) peakCandidates(jj,indexFlatSlope) = arma::zeros<arma::mat>(1,indexFlatSlope.n_elem); // reset index of flat slopes on derivative
			}

			// range of candidates based on first derivative
			if (signDeltaSlope == 1) {
				signFlag = true;
			}
			else if (signDeltaSlope == 0) {
				if (signFlag == false) {
					walkOnDerivative(j,i-1) = 1; //set index of flat slope on derivative (temporary index)
				}
			}
 			else if(signDeltaSlope == -1) {
				if (signFlag == false) {
					walkOnDerivative(j,i-1) = 1; //set index of flat slope on derivative (temporary index)
				}
				signFlag = false;
			}
		} // end of for: i
		arma::urowvec jj(1); jj = j;

		indexFlatSlope.clear();
		indexFlatSlope = arma::find(walkOnDerivative(j,arma::span::all) == 2);
		if (!indexFlatSlope.is_empty()) walkOnDerivative(jj,indexFlatSlope) = arma::zeros<arma::mat>(1,indexFlatSlope.n_elem); // reset index of flat slopes on derivative (temporary index)

		indexFlatSlope.clear();
		indexFlatSlope = arma::find(peakCandidates(j,arma::span::all) == 2);
		if (!indexFlatSlope.is_empty()) peakCandidates(jj,indexFlatSlope) = arma::zeros<arma::mat>(1,indexFlatSlope.n_elem); // reset index of local maxima on derivative (temporary index)

	} // end of for: j

    candidateRangeBasedOnDerivative = arma::max(walkOnDerivative,0) % wave(arma::span(1, wave.n_elem -2)); // vector of intersection of moving zero crossing line * amplitude of each point
	arma::rowvec peakCandidatesMax = arma::max(peakCandidates,0);
	candidatePeaksPosition = arma::find(peakCandidatesMax > 0); // place of peaks (intersection with slope = 0 ) based on first derivative
}

/// finds the range of each candidate based on derivative and previous/next candidate info
void delineate::candidateRangeInfoFinder(const arma::rowvec& derivative, const arma::rowvec& movedZeroCrossingAmplitutedMax, std::vector<candidate>& candids, const int looseWindow) {

	candidate candidTmp; // temporary candidate
    	arma::uvec candidatesPoints = arma::find( movedZeroCrossingAmplitutedMax != 0); // points of signal that have intersection with zero crossing line

	if (candidatesPoints.n_elem == 0)
		return; // no candidate gets funded

    	arma::uvec index = arma::find(derivative(arma::span(0,candidatesPoints(0))) == arma::min(derivative(arma::span(0,candidatesPoints(0)))));

	candidTmp.set_candidateRangeInfo(index(index.n_elem -1), 0);
	candidTmp.set_candidateRangeInfo(candidatesPoints(0), 1);
	candidTmp.set_risingRangeInfo(candidatesPoints(0), 0);
	candids.push_back(candidTmp);

	int i = 1;
	while (i < static_cast<int>(candidatesPoints.n_elem)) {
		int wndowUpperBound = std::min(static_cast<int>(candidatesPoints.n_elem-1), i+looseWindow);
		arma::uvec candidateIndexesRange = candidatesPoints(arma::span(i,wndowUpperBound)) -  arma::as_scalar(candidatesPoints(i));
		if (static_cast<int>(accu(candidateIndexesRange)) > ((looseWindow * (looseWindow+1) / 2))) {
			int fraction = 0;
			for (int j = 0; j < looseWindow; ++j) {
				if (candidateIndexesRange(j) + 1 != candidateIndexesRange(j+1)) {
					fraction = j;
					break;
				}
			}

			i += fraction;
			arma::uvec indexMinSlope = arma::find(derivative(arma::span(candidatesPoints(i) +1, candidatesPoints(i+1) -1)) == arma::min(derivative(arma::span(candidatesPoints(i) +1, candidatesPoints(i+1) -1))));
			candids[candids.size()-1].set_candidateRangeInfo(indexMinSlope(indexMinSlope.n_elem-1) + candidatesPoints(i) +1, 1) ;
			candids[candids.size()-1].set_risingRangeInfo(candidatesPoints(i), 1);

			{ // new candidate
				candidTmp.set_candidateRangeInfo(candids[candids.size()-1].get_candidateRangeInfo(1), 0);
				candidTmp.set_risingRangeInfo(candidatesPoints(i+1), 0);
				candids.push_back(candidTmp);
			}
		}
		++i;
	}

	index = arma::find(derivative(arma::span(candidatesPoints(i-1), derivative.n_elem-1)) == arma::min(derivative(arma::span(candidatesPoints(i-1), derivative.n_elem-1))));
	candids[candids.size()-1].set_candidateRangeInfo(index(index.n_elem-1) + candidatesPoints(i-1), 1);
	candids[candids.size()-1].set_risingRangeInfo(candidatesPoints(i-1), 1);
}

/// labelling the candidates into: slur/ peak
void delineate::labellingPeaks(std::vector<candidate> &candids, const arma::uvec& candidatePeaksPosition) {
 
	 int deltaSample = 3;  // tolerance sample point number of each candidate (loose filter calculation)
	 for (auto& candid : candids) {
		 for(int position: candidatePeaksPosition) {
			 if ((position >= (candid.get_risingRangeInfo(0) - deltaSample)) && (position <= (candid.get_risingRangeInfo(1) + deltaSample))) {
				candid.set_label(candidateLabel::peak); // peak
				break;
			 }
		 }
	 }
}

/// re_labelling the slur candidates to: slur/ rising slur/ falling slur
void delineate::reLabellingSlurs(std::vector<candidate>& candids) {
	for (int i = 0; i < static_cast<int>(candids.size()); ++i) {
		if (candids[i].get_label() == 0 && candids[std::max(0,i-1)].get_label() == 2) { // falling slur
			if (candids[i].get_a0() <= 0 && candids[i].get_a1() <= 0) //  slur after a peak with negative slopes
			candids[i].set_label(candidateLabel::slurFalling);
		}
		else if (candids[i].get_label() == 2 && candids[std::max(0,i-1)].get_label() == 0) { // rising slur
			if (candids[i-1].get_a0() >= 0 && candids[i-1].get_a1() >= 0) // slur before a peak with positive slopes
				candids[i-1].set_label(candidateLabel::slurRising);
			}
		}
}

/// re-adjusts the end of Twave based on Tpeak and current Toff
double delineate::readjustToff(const arma::rowvec wave, annotation &anns, double rr, double rpeak) {
	if (anns.peak.size() < 1) return -1;    // no annotation
	if (anns.peak[0] > anns.off) return -1; // incorrect toff
	if (anns.peak[0] < rpeak) return -1;    // incorrect tpeak

	double toff = anns.off > arma::as_scalar(wave.n_elem) -1 ? arma::as_scalar(wave.n_elem) -1 : anns.off; // re-adjust toff: if toff > length(wave)
	double tpeakAmp = arma::as_scalar(wave(anns.peak[0]));
	double lastCandidate = anns.lastCandidate;

	if (anns.peak.size() > 1) { // notch
		if (anns.peak[1] > anns.off) return -1; // incorrect toff
		if (wave(anns.peak[1]) > tpeakAmp) {
			tpeakAmp = arma::as_scalar(wave(anns.peak[1]));
		}
	}

	double adjusted_toff = toff;
	double lowerBound = lastCandidate+(std::trunc((toff-lastCandidate)/2));
	arma::uvec highAmpIndex = arma::find(wave(arma::span(lowerBound,toff)) > tpeakAmp,1);
	if (highAmpIndex.n_elem > 0) {
		adjusted_toff = lowerBound + arma::as_scalar(highAmpIndex(0));
	}

	arma::rowvec lastCandidToff_segment = wave(arma::span(lastCandidate,adjusted_toff));

	double toff_new = delineate::newToffFunc(lastCandidToff_segment, rr, rpeak, lastCandidate);

	return (toff_new > 0) ? toff_new : adjusted_toff;
}

/// finds a new Toff based on energy/cost function
double delineate::newToffFunc(const arma::rowvec& lastCandidToff_segment, double rr, double rpeak, double lastCandid) {

	// this function finds new toff candidates based on energy of lastCandidToff_segment and chooces one of them based on cost function (F(distance,energy))
	double toff_new = 0;
	try {
		if (lastCandidToff_segment.n_elem < 2) {
			return toff_new;
		}

		/* step01 : calculates smooth derivative */
		arma::rowvec derivative = lastCandidToff_segment(arma::span(1,lastCandidToff_segment.n_elem-1)) - lastCandidToff_segment(arma::span(0,lastCandidToff_segment.n_elem-2)); // derivative of lastCandidToff_segment
		delineate::smoothWaveFunc(derivative, std::min(5.0,arma::as_scalar(derivative.n_elem)+0.0)); // smooth filter of derivative

		/* step02 : calculates energy of lastCandidToff_segment. The peak has higher energy and toff should have least amount of energy if toff is global minima in lastCandidToff_segment.
				also finds the new toff candidates (indexlocalMinimaEnergy) */
		arma::rowvec energy = arma::zeros<arma::rowvec>(derivative.n_elem); // energy of lastCandidToff_segment: higher amplitude has less energy
		std::vector<int> indexlocalMinimaEnergy; // toff candidates: index of local minima of energy signal
		energy(0) = lastCandidToff_segment(0);
		bool flagmaxima = true;
		for (std::size_t i = 1; i < derivative.n_elem; ++i) {
			if (derivative(i) < 0) {     // negative derivative: decreases energy
				energy(i) = energy(i-1) - lastCandidToff_segment(i);
				if (flagmaxima == true) indexlocalMinimaEnergy.push_back(i);
				indexlocalMinimaEnergy[indexlocalMinimaEnergy.size()-1] = i;
				flagmaxima = false;
			}
			else if (derivative(i) > 0) { // positive derivative: increases energy
				energy(i) = energy(i-1) + lastCandidToff_segment(i);
				flagmaxima = true;
			}
			else {  // zero derivative: preserves energy
				energy(i) = energy(i-1);
			}
		}

		arma::rowvec energyNormal = 100.0*(energy-energy.min())/(energy.max()-energy.min()); // normalize energy between [0-100] due to current lastCandidToff. Higher amount of energy has higher amplitude

		/* step03 : finds first faling point on the energy. energy(0) should have highest value if lastCandid (tpeak) is chosen appropriately */
		int startEnergy = 0;
		for (std::size_t i = 0; i < energyNormal.n_elem-1; ++i) {
			if (energyNormal(i) > energyNormal(i+1)) {
				startEnergy = i;
				break;
			}
		}
		/* step04 : re-adjust the index of toff candidates (indexlocalMinimaEnergy) based on derivative of energy. The greatest localmaxima will be a new toff of each candidate */
		arma::rowvec derivativeEnergy = energyNormal(arma::span(1,energyNormal.n_elem-1)) - energyNormal(arma::span(0,energyNormal.n_elem-2));
		std::vector<int> indexNewToffCandidates(indexlocalMinimaEnergy.size(),0); //= arma::zeros<arma::uvec>(indexlocalMinimaEnergy.size());
		int j = startEnergy;
		for (int i = 0; i < static_cast<int>(indexlocalMinimaEnergy.size()); ++i) {

			int toff_index2 = indexlocalMinimaEnergy[i];

			arma::uvec indexMin = arma::find(derivative(arma::span(j,toff_index2)) == min(derivative(arma::span(j,toff_index2))),1);
			int toff_index1 = j + arma::as_scalar(indexMin(0));

			// adjust toff_index1 and toff_index2
			if (toff_index1 > toff_index2) {
				int tmp = toff_index1;
				toff_index1 = toff_index2;
				toff_index2 = tmp;
			}
			if(toff_index1 == toff_index2){toff_index1--;} //Indices out of bounds check.

			// derivative of current toff candidate
			arma::rowvec derivativeEnergyCandidate = derivativeEnergy(arma::span(toff_index1,toff_index2-1));

			// cautious!!!: 3 steps of smooth and precisions can affect the new toff index
			//(1.) smooth the derivative of current toff candidate
			delineate::smoothWaveFunc(derivativeEnergyCandidate, std::trunc(arma::as_scalar(derivativeEnergyCandidate.n_elem)/10) + 1);
			//(2.) fixed the precision of derivative of current toff candidate
			derivativeEnergyCandidate = arma::floor(derivativeEnergyCandidate*10000.0)/10000.0;

			// new place of toff based on current candidate is the greatest local maxima of derivativeEnergy of this candidate
			bool falingFlag = false;
			int indexToffCandidate = 1;
			for (std::size_t k = 1; k < derivativeEnergyCandidate.n_elem; ++k) {

				if (arma::as_scalar(derivativeEnergyCandidate(k-1) - derivativeEnergyCandidate(k)) < -0.001) { // (3.) -0.001 is for using precision
					falingFlag = false;
				}

				// rising of derivative
				if (arma::as_scalar(derivativeEnergyCandidate(k-1) - derivativeEnergyCandidate(k)) > -0.001) { // (3.) -0.001 is for using precision
					if (falingFlag == false) {
						falingFlag = true; // local maxima
						if (derivativeEnergyCandidate(k) > derivativeEnergyCandidate(indexToffCandidate)) {
							indexToffCandidate = k;// first point of local maxima
						}
					}
				}
			}

			// if derivativeEnergy has not any local maxima, use the maximum size of it as a new place of toff
			indexToffCandidate = (indexToffCandidate == 1) ? derivativeEnergyCandidate.n_elem-1 : indexToffCandidate;

			indexNewToffCandidates[i] = toff_index1 + indexToffCandidate;
			j = indexlocalMinimaEnergy[i] + 1;
		}

		/* step05 : chooses one of the new toff candidates for re-adjusted toff based on cost function (F(distance,energy))*/
		// distance parameter of cost function
		double a = lastCandid - rpeak;
		arma::rowvec b = (lastCandid+arma::conv_to<arma::rowvec>::from(indexNewToffCandidates)) - rpeak;
		arma::rowvec x = b - a;
		double y = rr - a;
		arma::rowvec D = x/y;

		// energy parameter of cost function
		arma::uvec indextmp = arma::conv_to<arma::uvec>::from(indexNewToffCandidates);
		arma::rowvec E = energyNormal.cols(indextmp)/100.0;

		// cost function
		arma::rowvec costFunc = E + D;

		// new toff obtains by minimizing cost function
		arma::uvec indexCandid = find(costFunc == costFunc.min(),1);
		toff_new = lastCandid + indexNewToffCandidates[indexCandid(0)]; // final toff candidate

	} catch(const std::exception &e){
		std::string line = std::string("Caught exception (Could not delineate ECG, error in TOFF function): ") + e.what();
		std::cerr << line;
		throw;
	} catch(const char* e){
		std::string line = std::string("Caught exception (Could not delineate ECG, error in TOFF function): ") + e;
		std::cerr << line;
		throw;
	} catch(...){
		std::string line = std::string("Unknown caught exception (Could not delineate ECG, error in TOFF function)");
		std::cerr << line;
		throw;
	}

	return toff_new;
}

/// makes a smooth filtering by taking an average in a sliding window on a signal
void delineate::smoothWaveFunc(arma::rowvec& wave, int windowSize){
	arma::rowvec waveTmp = wave;
	int sampleSize = arma::as_scalar(wave.n_elem);

	for ( int i = 0; i < sampleSize; ++i) {
		int windowLowerBound = std::max(0,i-windowSize);
		int windowUpperBound = std::min(sampleSize-1,i+windowSize);
		waveTmp(i) = arma::mean(wave(arma::span(windowLowerBound,windowUpperBound)));
	}

	wave = waveTmp; // output
}

/// cleans unwanted candidates
void delineate::cleanUpCandidates(std::vector<candidate>& candids, std::vector<int>& indexCandidates) {
	std::size_t i = 0;	
	while(indexCandidates.size()) {
		candids.erase(candids.begin() + indexCandidates.at(0) - i);
		indexCandidates.erase(indexCandidates.begin());
		++i;
	}
	eraseList(indexCandidates); // clean up the container
}

/// cleans a list
template <class T>
void delineate::eraseList(std::vector<T>& list) {
	list.erase(list.begin(),list.end());

}

