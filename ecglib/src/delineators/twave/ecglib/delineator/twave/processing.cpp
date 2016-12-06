/**
 * @file delineators/twave/ecglib/delineator/twave/processing.cpp
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
* Processing rules of twave delineation and candidates.
*/

#include "processing.hpp"
#include <algorithm>    // std::min_element
#include <iterator>     // std::distance
#include <cmath>     	// std::trunc
#include "delineateFinder.hpp" // delineateFinder::peakOriginFinder()

using namespace ecglib::twaveDelineate;

/*                          */
/*-----   PreProcesing   ---*/
/*                          */

/// finds candidates with few number of points in X axes (time)
std::vector<int> preProcessingRules::fewPointsCandidates(const std::vector<candidate>& candids, double minPoints) {

	std::vector<int> indexCandidates;
	for(std::size_t i = 0; i < candids.size(); ++i) {
		if ((candids[i].get_risingRangeInfo(1) - candids[i].get_risingRangeInfo(0)) < minPoints) {
			indexCandidates.push_back(i);
		}
	}

	return indexCandidates;
}

/*                          */
/*-----   PostProcesing  ---*/
/*                          */

///  finds candidates in a certain percentage lower than main peak, if its amplitued is too low
std::vector<int> postProcessingRules::lowAmplitudeMainPeak(const std::vector<candidate>& candids, double minValidAmplitudeMainPeak, double percentMainePeak) {

	std::vector<int> allCandidates; // index of candidates

	std::vector<int> indexPeakCandidates = postProcessingRules::peakCandidates(candids); // index of peak candidates

	int mainPeakIndex = -1;
	double mainPeakAmplitude = 0;
	postProcessingRules::mainPeak(candids, indexPeakCandidates, mainPeakIndex, mainPeakAmplitude); // main peak

	if (mainPeakAmplitude < minValidAmplitudeMainPeak) {
		for (std::size_t i = 0; i < candids.size(); ++i) {
			if (candids[i].get_y() < mainPeakAmplitude*percentMainePeak) {
				allCandidates.push_back(i);
			}
		}
	}

	return allCandidates;
}

/// finds peak candidates with lower amplitude in compare with a percentage of highest peak and minimum valid amplitude
std::vector<int> postProcessingRules::lowAmplitudePeaks(const std::vector<candidate>& candids, double minValidAmplitude, double percentPeak) {

	std::vector<int> indexPeakCandidates = postProcessingRules::peakCandidates(candids); // index of peak candidates

	int mainPeakIndex = -1;
	double mainPeakAmplitude = 0;
	postProcessingRules::mainPeak(candids, indexPeakCandidates, mainPeakIndex, mainPeakAmplitude); // main peak

	std::vector<int> indexLowAmplitude; // index of peak candidates	with low amplitude
	for(std::size_t i = 0; i < indexPeakCandidates.size(); ++i) {
		double peakMostlySlur = (mainPeakAmplitude > minValidAmplitude) ? (mainPeakAmplitude - minValidAmplitude) * percentPeak - (candids[indexPeakCandidates[i]].get_y() - minValidAmplitude) : mainPeakAmplitude * (1-percentPeak) - candids[indexPeakCandidates[i]].get_y();
		if (peakMostlySlur > 5e-3) { // 5e-3 == 0.004 just for sanity check
			indexLowAmplitude.push_back(indexPeakCandidates[i]);
		}
	}

	return indexLowAmplitude;
}

/// finds peaks are not close to main peak in terms of amplitude
std::vector<int> postProcessingRules::inconsistentPeaks(const std::vector<candidate>& candids, double maxDelatAplitudeNotches) {

	std::vector<int> indexInconsistentPeaks;
	std::vector<int> indexPeakCandidates = postProcessingRules::peakCandidates(candids); // index of peak candidates

	int mainPeakIndex = -1;
	double mainPeakAmplitude = 0;
	postProcessingRules::mainPeak(candids, indexPeakCandidates, mainPeakIndex, mainPeakAmplitude); // main peak

	for  (std::size_t i = 0; i < indexPeakCandidates.size(); ++i) {
		double deltaAmplitudePeaks = mainPeakAmplitude - candids[indexPeakCandidates[i]].get_y();
		if (deltaAmplitudePeaks > maxDelatAplitudeNotches) {
			indexInconsistentPeaks.push_back(indexPeakCandidates[i]);
		}
	}
	return indexInconsistentPeaks;
}

/// finds slur candidates
std::vector<int> postProcessingRules::unrelatedSlure(const std::vector<candidate>& candids) {
	std::vector<int> slurIndex = postProcessingRules::slurCandidates(candids); // index of slur candidates
	return slurIndex;
}

/// merges two consequent condidates if they are close in terms of amplitude
std::vector<int> postProcessingRules::meargingCandidates(const arma::rowvec& wave, std::vector<candidate>& candids, double minAmplitudeFlatness){

	std::vector<int> indexMargedCandidates;
	std::size_t i = 1;
	while(i < candids.size()) {
		if (candids[i].get_label() != 0 && candids[i-1].get_label() != 0) {
			if (abs(candids[i].get_label()) != 1 || abs(candids[i-1].get_label()) != 1) {
				int xIntersect = 0;
				postProcessingRules::intersectionTwoCandidates(candids[i-1], candids[i], xIntersect);

                		if ((std::max(candids[i-1].get_y(), candids[i].get_y()) - wave[xIntersect] < minAmplitudeFlatness) && (std::abs(candids[i-1].get_y() - candids[i].get_y()) < minAmplitudeFlatness)) {

					candidate mergedCandidate;
					mergedCandidate.set_a0(candids[i-1].get_a0());
					mergedCandidate.set_b0(candids[i-1].get_b0());

					mergedCandidate.set_a1(candids[i].get_a1());
					mergedCandidate.set_b1(candids[i].get_b1());

					mergedCandidate.set_x(std::trunc((candids[i-1].get_x() + candids[i].get_x())/2));
					mergedCandidate.set_y(wave[mergedCandidate.get_x()]);

					mergedCandidate.set_flatnessSamples(candids[i].get_x() - candids[i-1].get_x() +1 + std::round((candids[i].get_flatnessSamples() + candids[i-1].get_flatnessSamples())/2));

					mergedCandidate.set_candidateRangeInfo(candids[i-1].get_candidateRangeInfo(0), 0);
					mergedCandidate.set_candidateRangeInfo(candids[i].get_candidateRangeInfo(1), 1);

					mergedCandidate.set_risingRangeInfo(candids[i-1].get_risingRangeInfo(0), 0);
					mergedCandidate.set_risingRangeInfo(candids[i].get_risingRangeInfo(1), 1);

					// re-calculates of peak origin based on new merged slopes line
					int xIntersect(0);
					double yIntersect(0), angleIntersect(0);
					delineateFinder delineat;
					delineat.peakOriginFinder(wave(arma::span(mergedCandidate.get_candidateRangeInfo(0),mergedCandidate.get_candidateRangeInfo(1))),
									mergedCandidate.get_candidateRangeInfo(0), mergedCandidate.get_a0(), mergedCandidate.get_b0(), mergedCandidate.get_a1(), 
									mergedCandidate.get_b1(), xIntersect, yIntersect, angleIntersect);			

					mergedCandidate.set_xOrigin(xIntersect);
					mergedCandidate.set_yOrigin(yIntersect);

					mergedCandidate.set_skewness(angleIntersect);

					double peakDistX = std::pow(mergedCandidate.get_x() - mergedCandidate.get_xOrigin(), 2.);
					double peakDistY = std::pow(mergedCandidate.get_y() - mergedCandidate.get_yOrigin(), 2.);
					mergedCandidate.set_distortion(std::sqrt(peakDistX + peakDistY));

					mergedCandidate.set_label(candidateLabel::peak);

					// these two indexes should be removed
					indexMargedCandidates.push_back(i-1);
					indexMargedCandidates.push_back(i);

					// the new merged candidate should be added to rest of candidates
					candids.insert(candids.begin() + i+1, mergedCandidate);
					++i;
				}
			}
		}
		++i;
	} // end of while
	return indexMargedCandidates;
}

/// distinguishes between good slur and bad slur based on classification rules
std::vector<int> postProcessingRules::slurClassifier(const std::vector<candidate>& candids, const std::vector<std::vector<double> >& featursThreshold) {		

	// making feature set for classifying slur_peak to distinguish between slur and non-slur (bad slur)
	// if the output of classifier is 1: slur should be removed.
	// if the output of classifier is 2: slur should be preserved.
	// for creating model of classification, 5 features were used after feature extraction and feature reduction.
	// the results is based on 630 input samples contains 160 train samples and 470 validation samples(470: 370 separated samples & 100 samples from same study of train samples)
	// in total contains 40 positive sample from 630.
	// after training and enhancing the model by SVM and decision-tree, few rules get extracted by decision-tree that are used here.

	// feature0 : abs(atand(a0(slur))-atand(a1(slur)));						// angle between two slopes of slur 
	// feature1 : abs(atand(a0(peak))-atand(a1(slur))) || abs(atand(a1(peak))-atand(a0(slur)));	// angle between slopes of slur and peak
	// feature2 : _y(peak)/_yOrigin(slur));								// ratio between amplitude of peak and slur
    	// feature3 : _yOrigin(peak));									// amplitude of slur
	// feature4 : abs(_xOrigin(slur)-intersection(peak,slur));					// distance between slur and interconnection

	std::vector<int> notValidSlur;

	std::vector<int> slurIndex = postProcessingRules::risingFallingslurCandidates(candids);
	const double dpi = 180 / 3.1415; // for converting radian to degree
	for (std::size_t i = 0; i < slurIndex.size(); ++i) {
		std::vector<double> featureSet;

		int s = slurIndex[i]; // s: slur index
		int p = s; 	      // p: peak index
		if (candids[s].get_label() == 1)
			p = p + 1;
		else if (candids[s].get_label() == -1)
		    	p = p - 1;

		double at2_1 = std::atan(candids[s].get_a0())*dpi;   // slur
		double at2_2 = std::atan(candids[s].get_a1())*dpi;   // slur

		double at1_1 = std::atan(candids[p].get_a0())*dpi;   // peak
		double at1_2 = std::atan(candids[p].get_a1())*dpi;   // peak

		featureSet.push_back(std::abs(at2_1 - at2_2));				   // F0: angle between two slopes of slur 

		if (candids[s].get_label() == 1)					   // F1: angle between slopes of slur and peak
		    featureSet.push_back(std::abs(at1_1 - at2_1));
		else if (candids[s].get_label() == -1)
		    featureSet.push_back(std::abs(at1_2 - at2_2));

		featureSet.push_back(candids[p].get_y()/candids[s].get_yOrigin());	   // F2: ratio between amplitude of peak and slur
		/*
		featureSet.push_back(candids[s].get_yOrigin());				   // F3: amplitude of slur

		int xIntersect = 0;
		postProcessingRules::intersectionTwoCandidates(candids(p), candids(s), xIntersect);
		featureSet.push_back(std::abs(candids[s].get_xOrigin() - xIntersect));  // F4: distance between slur and interconnection

		// classifier
		int newLable = predict(MLmodel, featureSet);
		if (newLable == 1) //  candidate is not slur then we should remove it from candidate list
		    notValidSlur.push_back(j);
		*/

		// extracted rules
		if (featureSet[0] < featursThreshold[0][0]/*20*/) {
		    if (!(featureSet[0] > featursThreshold[1][0]/*10*/ && featureSet[1] > featursThreshold[1][1]/*10*/ && featureSet[2] <featursThreshold[1][2]/*1.5*/)) // candidate is a bad slur
			notValidSlur.push_back(s);
		}
		else {
			if (featureSet[2] > featursThreshold[2][2]/*1.7*/) // candidate is a bad slur
				notValidSlur.push_back(s);
		}

	} // end of for

	return notValidSlur;
}


/// keeps main and second peaks
std::vector<int> postProcessingRules::keepJustTwoPeaks(const std::vector<candidate>& candids) {

	std::vector<int> indexPeakCandidates = postProcessingRules::peakCandidates(candids);

	std::vector<candidate> candidsTMP;
	for (std::size_t i = 0; i < indexPeakCandidates.size(); ++i) {
		candidsTMP.push_back(candids[indexPeakCandidates[i]]);
	}

	int mainPeakIndex = -1;
	double mainPeakAmplitude = 0;
	postProcessingRules:: candidateMax(candidsTMP, mainPeakAmplitude, mainPeakIndex);

	int secondPeakIndex = -1;
	double secondPeakAmplitude = 0;
	postProcessingRules:: candidateSecondMax(candidsTMP, secondPeakAmplitude, secondPeakIndex);
	candidsTMP.clear();


    	indexPeakCandidates.erase(indexPeakCandidates.begin() + mainPeakIndex); // remove mainPeakIndex

        if (secondPeakIndex != -1){ // remove secondPeakIndex
		if (secondPeakIndex > mainPeakIndex){
			indexPeakCandidates.erase(indexPeakCandidates.begin() + secondPeakIndex - 1);
		}else{
			indexPeakCandidates.erase(indexPeakCandidates.begin() + secondPeakIndex);
		}
	}
	return indexPeakCandidates; // index peak candidates without index of main(first) & second peaks
}

///  finds all peaks that have one flat side and converts them to slur
std::vector<int> postProcessingRules::convertPeakToSlur(const arma::rowvec& wave, const std::vector<candidate>& candids, double minValidAmplitudePeak) {

	std::vector<int> newSlurCondidates;
	std::vector<int> indexPeakCandidates = postProcessingRules::peakCandidates(candids);

	int mainPeakIndex = -1;
	double mainPeakAmplitude = 0;
	postProcessingRules::mainPeak(candids, indexPeakCandidates, mainPeakIndex, mainPeakAmplitude);

	if(indexPeakCandidates.size() > 0) indexPeakCandidates.erase(indexPeakCandidates.begin() + mainPeakIndex); // remove mainPeakIndex

	for(std::size_t i = 0; i < indexPeakCandidates.size() && indexPeakCandidates.size() > 0; ++i) {
		if(indexPeakCandidates[i] > 0)  { // there is a local minima before peak
			int xIntersect = 0;
			postProcessingRules::intersectionTwoCandidates(candids[indexPeakCandidates[i]-1], candids[indexPeakCandidates[i]], xIntersect);
			if ((candids[indexPeakCandidates[i]].get_y() - wave(xIntersect)) < minValidAmplitudePeak) {
				newSlurCondidates.push_back(indexPeakCandidates[i]);
				continue; // dont look at the other side of peak
		    	}
		}
		 if(indexPeakCandidates[i] < static_cast<int>(candids.size()) -1)  { // there is a local minima after peak
			int xIntersect = 0;
			postProcessingRules::intersectionTwoCandidates(candids[indexPeakCandidates[i]], candids[indexPeakCandidates[i]+1], xIntersect);
			if((candids[indexPeakCandidates[i]].get_y() - wave(xIntersect)) < minValidAmplitudePeak) {
				newSlurCondidates.push_back(indexPeakCandidates[i]);
		 	}
		}
	}

	return newSlurCondidates;
}

///  labels a non-measurable signal
bool postProcessingRules::nonMeasurableSignal(std::vector<candidate>& candids, double nonMeasurableVoltage) {
	std::vector<int> indexPeakCandidates = postProcessingRules::peakCandidates(candids);

	int mainPeakIndex = -1;
	double mainPeakAmplitude = 0;
	postProcessingRules::mainPeak(candids, indexPeakCandidates, mainPeakIndex, mainPeakAmplitude);
	return ( mainPeakAmplitude < nonMeasurableVoltage); // returns true if it is non-measurable
}

/// finds an intersection between slops of two candidates
void postProcessingRules::intersectionTwoCandidates(const candidate& candids1, const candidate& candids2, int& xIntersect) {

	// y = ax + b
	// y = cx + d
	// => ax + b = cx + d => x = (d-b)/(a-c)
	double db(0), ac(0);
	if (candids2.get_x() > candids1.get_x()) {
		db = candids1.get_b1() - candids2.get_b0();
		ac = candids2.get_a0() - candids1.get_a1();
	}
	else {
		db = candids2.get_b1() - candids1.get_b0();
		ac = candids1.get_a0() - candids2.get_a1();
	}
	if (ac > 1e-10) // checks for dividing by zero
		xIntersect = std::ceil(std::abs(db/ac));
	else
		xIntersect = -1;
	if (xIntersect < candids1.get_risingRangeInfo(2) || xIntersect > candids2.get_risingRangeInfo(1)) { // two candidates have a same slope
		xIntersect = std::ceil((candids1.get_risingRangeInfo(2) + xIntersect > candids2.get_risingRangeInfo(1))/2);
	}
}

std::vector<int> postProcessingRules::typeCandidates(const std::vector<candidate>& candids, int type) { // index of candidates specified by type
	std::vector<int> indexCandidates;
	for(std::size_t i = 0; i < candids.size(); ++i) {
		if (candids[i].get_label() == type)
			indexCandidates.push_back(i);
	}
	return indexCandidates;
}

std::vector<int> postProcessingRules::peakCandidates(const std::vector<candidate>& candids) { // index of peak candidates
	return postProcessingRules::typeCandidates(candids, 2);
}

std::vector<int> postProcessingRules::risingSlurCandidates(const std::vector<candidate>& candids) { // index of rising slur candidates
	return postProcessingRules::typeCandidates(candids, 1);
}

std::vector<int> postProcessingRules::fallingSlurCandidates(const std::vector<candidate>& candids) { // index of falling slur candidates
	return postProcessingRules::typeCandidates(candids, -1);
}

std::vector<int> postProcessingRules::slurCandidates(const std::vector<candidate>& candids) { // index of slur candidates
	return postProcessingRules::typeCandidates(candids, 0);
}

std::vector<int> postProcessingRules::risingFallingslurCandidates(const std::vector<candidate>& candids) { // index of rising & falling slur candidates
	std::vector<int> slurIndex = postProcessingRules::risingSlurCandidates(candids);
	std::vector<int> fslurIndex = postProcessingRules::fallingSlurCandidates(candids);
	slurIndex.insert(slurIndex.end(), fslurIndex.begin(), fslurIndex.end());
	return slurIndex;
}

/// returns index and amplitude of main peak
void postProcessingRules::mainPeak(const std::vector<candidate>& candids, const std::vector<int>& indexPeakCandidates, int& maxPeakIndex, double& maxPeakAmplitude) {

	std::vector<candidate> candidsTMP;
	for (std::size_t i = 0; i < indexPeakCandidates.size(); ++i) {
		candidsTMP.push_back(candids[indexPeakCandidates[i]]);
	}

	postProcessingRules::candidateMax(candidsTMP, maxPeakAmplitude, maxPeakIndex);
	candidsTMP.clear();

	return;
}

/// internal comparision function
bool cmp(const candidate& lhs, const candidate& rhs) {
  return lhs.get_y() < rhs.get_y();
}

/// returns index and amplitude of main peak in a peak candidats' list
void postProcessingRules::candidateMax(const std::vector<candidate>& candids, double& amp, int& index) {
	if (candids.size() > 0) {
		auto max_it = std::max_element(candids.begin(), candids.end(), cmp);
		amp = max_it->get_y();
		index = std::distance(candids.begin(), max_it);
	}
}

/// returns index and amplitude of second peak
void postProcessingRules::candidateSecondMax(const std::vector<candidate>& candids, double& amp, int& index) {

	if (candids.size() > 0) {
		auto max_it = std::max_element(candids.begin(), candids.end(), cmp);
		int max_index = std::distance(candids.begin(), max_it);
		auto max_it_part1 = std::max_element(candids.begin(), candids.begin() + max_index, cmp);
		int max_index_part1 = std::distance(candids.begin(), max_it_part1);
		auto max_it_part2 = std::max_element(candids.begin() + max_index + 1, candids.end(), cmp);
		int max_index_part2 = std::distance(candids.begin(), max_it_part2);

		index = max_index != 0 ? (max_index < static_cast<int>(candids.size()) -1 ? (max_it_part1->get_y() > max_it_part2->get_y() ? max_index_part1: max_index_part2) : max_index_part1) : (max_index < static_cast<int>(candids.size()) -1 ? max_index_part2 : -1);
		amp = index != -1 ? (index != max_index_part1 ? max_it_part2->get_y() : max_it_part1->get_y()) : 0;
	}
}

