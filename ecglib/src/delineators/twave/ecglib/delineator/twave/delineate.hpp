/**
 * @file delineators/twave/ecglib/delineator/twave/delineate.hpp
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

#ifndef ECGLIB_DELINEATORS_TWAVE_DELINEATE_CANDIDATES_MH_2015_12_09
#define ECGLIB_DELINEATORS_TWAVE_DELINEATE_CANDIDATES_MH_2015_12_09 1

#include <vector>
#include <armadillo>

#include "generalstructure.hpp"
#include "processing.hpp"
#include "delineateFinder.hpp"

namespace ecglib {
	/*! \addtogroup delineator-twave
	 * T-wave delineator classes and functions
	 * @{
	 */

/**
 * @brief Twave delineate namespace
 */
namespace twaveDelineate {

	/**
	 * @brief Finds the twave delineations
	 */
	class delineate: // inherited from
			protected preProcessingRules,  // pre-processing rules class
			protected postProcessingRules, // post-processing rules class
			protected delineateFinder {    // delineator finder functionalities

		public:
			delineate() {};		/**< @brief constructor for delineate */
			~delineate() {};	/**< @brief destructor for delineate */

				/**
				 * @brief Main function for finding twave delineation based on twave's candidates containing peaks & slurs
				 *
			 	 * @param twave Input filtered twave
				 * @param pointStart Start point of a twave on a ecg signal
				 * @param featursThreshold Threshold vectors of extracted rules by Decision-Tree for slur classifier
				 * @param candidateFinderFlag Option for finding candidates based on moving zero crossing line (1) or first/second derivative (2) functions
				 * @param deltaStepSlope Interval value for calculating the number of moving zero crossing lines
				 * @param looseWindow Min points of a valid candidate
				 * @param minPoints Min points which make a candidate
				 * @param deltaAmplitude Delta amplitude of points that make a peak of candidate
				 * @param minVoltageMainPeak Minimum acceptable voltage of main peak
				 * @param percentMainePeak Percentage of main peak for evaluating the other candidates
				 * @param minVoltage Minimum acceptable voltage
				 * @param percentPeak Percentage of main peak for evaluating the other peaks
				 * @param maxDelatAplitudeNotches Acceptable delta amplitude differences of two peaks
				 * @param minAmplitudeFlatness Delta amplitudes of flatness
				 * @param minValidAmplitudePeak Threshsold of peak candidate that declares small angle
				 * @param measurableVoltage Min threshsold of peak amplitude as a measurable signal
				 *
				 * @return annotator Annotations of input twave
				 */
			annotation delineator(const arma::rowvec& twave, int pointStart, const std::vector<std::vector<double> >& featursThreshold,
					int candidateFinderFlag = 1, double deltaStepSlope = 10, int looseWindow = 10, int minPoints = 10, double deltaAmplitude = 5, double minVoltageMainPeak = 150,
	 				double percentMainePeak = 0.8, double minVoltage = 100, double percentPeak = 0.3, double maxDelatAplitudeNotches  = 50, 
					double minAmplitudeFlatness = 7, double minValidAmplitudePeak = 7, double measurableVoltage = 100);

				/**
				 * @brief Re-adjusts the end of twave based on tpeak and current toff
				 *
			 	 * @param wave Input filtered wave (vcg)
				 * @param annotation Twave annotation contains toff
			 	 * @param tpeaktoff_segment Input filtered wave (tpeak to toff segment)
				 * @param rr Value of rr interval
				 * @param rpeak Place of rpeak in an input wave (vcg)
				 *
				 * @return new toff index or -1 if toff has problem
				 */
			double readjustToff(const arma::rowvec wave, annotation &anns, double rr, double rpeak);

		private:
				/**
				 * @brief Finds all candidates of twave based on a clean (without major noise) moving zero crossing line on a first derivative vectore
				 *
			 	 * @param wave An input filtered wave for calculation of first derivative
				 * @param derivative Clean first derivative
				 * @param movedZeroCrossingAmplitutedMax (out var)Map of derivative to a vector. this vector shows the amplitude of signal when derivative has intersection with moving zero crossing line
				 * @param candidatePeaksPosition (out var)Shows the position of real peaks based on intersection of slope = 0 with first derivative
				 * @param deltaStepSlope Interval value of derivativ steps for calculating the number of moving zero crossing line
				 */
			void candidateFinder(const arma::rowvec& wave, const arma::rowvec& derivative, arma::rowvec& movedZeroCrossingAmplitutedMax, arma::uvec& candidatePeaksPosition, double deltaStepSlope);

				/**
				 * @brief Finds all candidates of twave based on a clean (first/second) derivative vectore
				 *
			 	 * @param wave An input filtered wave for calculation of first derivative
				 * @param derivative Clean first derivative
				 * @param candidateRangeBasedOnDerivative (out var)Shows the amplitude of wave signal when first derivative has falling curve
				 * @param candidatePeaksPosition (out var)Shows the position of real peaks based on intersection of slope = 0 with first derivative
				 */
			void candidateFinder2DerivativeBased(const arma::rowvec& wave, const arma::rowvec& derivative, arma::rowvec& candidateRangeBasedOnDerivative, arma::uvec& candidatePeaksPosition);

				/**
				 * @brief Finds the range of each candidate for finding candidate's delineation 
				 *
			 	 * @param derivative Clean first derivative
				 * @param movedZeroCrossingAmplitutedMax Intersection derivative with moving zero crossing line.
				 * @param candids (in/out var)Vector of candidates
				 * @param looseWindow Defines minimum acceptale points of a candidate
				 */
			void candidateRangeInfoFinder(const arma::rowvec& derivative, const arma::rowvec& movedZeroCrossingAmplitutedMax, std::vector<candidate>& candids, const int looseWindow);

				/**
				 * @brief Adds labels (slur & peak) into candidates 
				 *
				 * @param candids (in/out var)Vector of candidates
				 * @param candidatePeaksPosition The position of real peaks on the first derivative
				 */
			void labellingPeaks(std::vector<candidate>& candids, const arma::uvec& candidatePeaksPosition);

				/**
				 * @brief Relabels slur candidates by slur, rising-slur & falling-slur
				 *
				 * @param candids (in/out var)Vector of candidates
				 */
			void reLabellingSlurs(std::vector<candidate>& candids);

				/**
				 * @brief Find new Toff based on energy/cost function
				 *
			 	 * @param lastCandidToff_segment An input filtered wave (lastCandid(/tpeak) to toff segment)
				 * @param rr Value of rr interval
				 * @param rpeak Place of rpeak in an input wave (vcg)
				 * @param lastCandid Place of tpeak in an input wave (vcg)
				 */
			double newToffFunc(const arma::rowvec& lastCandidToff_segment, double rr, double rpeak, double lastCandid);

				/**
				 * @brief smooth Filters an input wave based on taking avarage by windowing
				 *
			 	 * @param wave (in/out var)An input wave (vcg)
				 * @param windowSize Size of smooth filtering
				 */
			void smoothWaveFunc(arma::rowvec& wave, int windowSize);

				/**
				 * @brief Remove part of candidates' vector
				 *
				 * @param candids (in/out var)Vector of candidates
				 * @param indexCandidates (in/out var)Candidates that should be removed
				 */
			void cleanUpCandidates(std::vector<candidate>& candids, std::vector<int>& indexCandidates);

				/**
				 * @brief Cleans a list
				 *
				 * @param list (in/out var)Input list that should get cleaned
				 */
			template <class T>
			void eraseList(std::vector<T>& list);
	}; // end of delineate class
} // end of namespace twaveDelineate

	/*!
	 * @}
	 */
} // end of namespace ecglib
#endif

