/**
* @file delineators/twave/ecglib/delineator/twave/processing.hpp
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

#ifndef ECGLIB_DELINEATORS_TWAVE_CANDIDATES_PROCESSING_MH_2015_12_09
#define ECGLIB_DELINEATORS_TWAVE_CANDIDATES_PROCESSING_MH_2015_12_09 1

#include <armadillo>
#include "generalstructure.hpp"

namespace ecglib {

	/*! \addtogroup delineator-twave
	 * T-wave delineator classes and functions
	 * @{
	 */


namespace twaveDelineate {

struct candidate; // forward declaration

/**
 * @brief pre-processing rules of twave delineation
 */
class preProcessingRules {
	public:
		/**
		 * @brief constructor of class
		 */
		preProcessingRules() {};

		/**
		 * @brief destructor of class
		 */
		~preProcessingRules() {};

	protected:
		/**
		 * @brief Finds candidates with few number of points
		 *
	 	 * @param candids candidate list
		 * @param minPoints min point that each candidate should have to pass a rule
		 *
		 * @return index of candidates that should be removed based on this rule
		 */
		std::vector<int> fewPointsCandidates(const std::vector<candidate>& candids, double minPoints);
};

/**
 * @brief post-processing rules of twave delineation
 */
class postProcessingRules {
	public:
		/**
		 * @brief constructor of class
		 */
		postProcessingRules() {}

		/**
		 * @brief destructor of class
		 */
		~postProcessingRules() {}

	protected:
		/* 					*/
		/* 		rules			*/
		/* 					*/


		/**
		 * @brief Finds candidates in a certain percentage lower than main peak, if its amplitued is too low
		 *
	 	 * @param candids candidate list
		 * @param minValidAmplitudeMainPeak the minimum valid amplitude of a main peak
		 * @param percentMainePeak the valid percentage of a main peak to keep candidates
		 *
		 * @return index of candidates that should be removed based on this rule
		 */
		std::vector<int> lowAmplitudeMainPeak(const std::vector<candidate>& candids, double minValidAmplitudeMainPeak, double percentMainePeak);

		/**
		 * @brief finds peak candidates with lower amplitude in compare with a percentage of highest peak and minimum valid amplitude
		 *
	 	 * @param candids candidate list
		 * @param minValidAmplitude the minimum valid amplitude of candidates
		 * @param percent the percentage of highest peak amplitude 
		 *
		 * @return index of candidates that should be removed based on this rule
		 */
		std::vector<int> lowAmplitudePeaks(const std::vector<candidate>& candids, double minValidAmplitude, double percentPeak);

		/**
		 * @brief Finds peaks which are not close to the main peak in term of amplitude
		 *
		 * @param candids candidate list
		 * @param maxDelatAplitudeNotches the valid distance of having notch in terms of amplitude
		 *
		 * @return index of candidates that should be removed based on this rule
		 */
		std::vector<int> inconsistentPeaks(const std::vector<candidate>& candids, double maxDelatAplitudeNotches);

		/**
		 * @brief Finds slur candidates
		 *
		 * @param candids candidate list
		 *
		 * @return index of candidates that should be removed based on this rule
		 */
		std::vector<int> unrelatedSlure(const std::vector<candidate>& candids);

		/**
		 * @brief Merges two consequent condidates if they are close in terms of amplitude
		 *
		 * @param wave twave signal
		 * @param candids candidate list
		 * @param minAmplitudeFlatness merging criterion threshold
		 *
		 * @return index of candidates that should be removed based on this rule
		 */
		std::vector<int> meargingCandidates(const arma::rowvec& wave, std::vector<candidate>& candids, double minAmplitudeFlatness);

		/**
		 * @brief Distinguishes good slur from bad slur based on classification rules (bad slur will be removed)
		 *
		 * @param candids candidate list
		 * @param featursThreshold thresholds of extracted rules (these values are based on extracted rules of Decision-Tree)
		 *
		 * @return index of candidates that should be removed based on this rule
		 */
		std::vector<int> slurClassifier(const std::vector<candidate>& candids, const std::vector<std::vector<double> >& featursThreshold);

		/**
		 * @brief Keeps main and second peaks
		 *
		 * @param candids candidate list
		 *
		 * @return index of candidates that should be removed based on this rule
		 */
		std::vector<int> keepJustTwoPeaks(const std::vector<candidate>& candids);

		/**
		 * @brief Finds all peaks that have one flat side and converts them to slur
		 *
		 * @param wave twave signal
		 * @param candids candidate list
		 * @param minValidAmplitudePeak the minimum ampitude that uses in calculating a flat side of a candidate
		 *
		 * @return index of candidates that should be removed based on this rule
		 */
		std::vector<int> convertPeakToSlur(const arma::rowvec& wave, const std::vector<candidate>& candids, double minValidAmplitudePeak);

		/**
		 * @brief Just an informative falg to label a non-measurable signal (low main peak amplitude)
		 *
		 * @param candids candidate list
		 * @param nonMeasurableVoltage the thresold amplitude to label non-measurable
		 *
		 * @return true if signal is non-measurable 
		 */
		bool nonMeasurableSignal(std::vector<candidate>& candids, double nonMeasurableVoltage);

        /*                          */
	 	/* 		end of rules		*/
        /*                          */


		/**
		 * @brief Finds one specific type of candidates in a candidate list
		 *
		 * @param candids candidate list
		 * @param type type of candidate slure, peak, rising slur, falling slur
		 *
		 * @return index of candidates
		 */
		std::vector<int> typeCandidates(const std::vector<candidate>& candids, int type);

		/**
		 * @brief Finds peak candidates in a candidate list
		 *
		 * @param candids candidate list
		 *
		 * @return index of candidates
		 */
		std::vector<int> peakCandidates(const std::vector<candidate>& candids);

		/**
		 * @brief Finds slur candidates in a candidate list
		 *
		 * @param candids candidate list
		 *
		 * @return index of candidates
		 */
		std::vector<int> slurCandidates(const std::vector<candidate>& candids);

		/**
		 * @brief Finds rising slur candidates in a candidate list
		 *
		 * @param candids candidate list
		 *
		 * @return index of candidates
		 */
		std::vector<int> risingSlurCandidates(const std::vector<candidate>& candids);

		/**
		 * @brief Finds falling slur candidates in a candidate list
		 *
		 * @param candids candidate list
		 *
		 * @return index of candidates
		 */
		std::vector<int> fallingSlurCandidates(const std::vector<candidate>& candids);

		/**
		 * @brief Finds rising and falling slur candidates in a candidate list
		 *
		 * @param candids candidate list
		 *
		 * @return index of candidates
		 */
		std::vector<int> risingFallingslurCandidates(const std::vector<candidate>& candids);

		/**
		 * @brief Finds main peak between all candidates
		 *
		 * @param candids candidate list
		 * @param indexPeakCandidates index of all peaks in candidate list
		 * @param maxPeakIndex (out var) index of main peak
		 * @param maxPeakAmplitude (out var) amplitude of main peak
		 */
		void mainPeak(const std::vector<candidate>& candids, const std::vector<int>& indexPeakCandidates, int& maxPeakIndex, double& maxPeakAmplitude);

		/**
		 * @brief Finds main peak in a peak list
		 *
		 * @param candids candidate list of all peaks
		 * @param amp (out var) amplitude of main peak
		 * @param index (out var) index of main peak
		 */
		void candidateMax(const std::vector<candidate>& candids, double& amp, int& index);

		/**
		 * @brief Finds second peak in a peak list
		 *
		 * @param candids candidate list of all peaks
		 * @param amp (out var) amplitude of second peak
		 * @param index (out var) index of second peak
		 */
		void candidateSecondMax(const std::vector<candidate>& candids, double& amp, int& index);

		/**
		 * @brief Finds intersection between two candidates based on intersection of their slopes line
		 *
		 * @param candids1 first candidate
		 * @param candids2 second candidate
		 * @param xIntersect (out var) index of intersection
		 */
		void intersectionTwoCandidates(const candidate& candids1, const candidate& candids2, int& xIntersect);
};
}
/*!
 * @}
 */
}
#endif
