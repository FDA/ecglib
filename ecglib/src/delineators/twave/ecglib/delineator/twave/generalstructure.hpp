/**
* @file delineators/twave/ecglib/delineator/twave/generalstructure.hpp
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
* General structure of twave delineation processes.
*/

#ifndef ECGLIB_DELINEATORS_TWAVE_CANDIDATES_STRUCTURES_MH_2015_12_09
#define ECGLIB_DELINEATORS_TWAVE_CANDIDATES_STRUCTURES_MH_2015_12_09 1

#include <string>
#include <vector>
#include <unordered_map>
namespace ecglib {
namespace twaveDelineate {

	/*! \addtogroup delineator-twave
	 * T-wave delineator classes and functions
	 * @{
	 */

/**
 * @brief annotators' exposure of Twave to out side
 */
struct annotation {
	public:
		/**
		 * @brief index (time) of start point of signal
		 */
		double on;
		/**
		 * @brief index (time) of peaks of signal (_x)
		 */
		std::vector<double> peak;
		/**
		 * @brief index (time) of end point of signal
		 */
		double off;
		/**
		 * @brief index (time) of last candidate(peak or slur) of signal for readjusting Toff
		 */
		double lastCandidate;
		/**
		 * @brief flatness number of samples in which creates peak of a candidate
		 */
		std::vector<double> flatness;
		/**
		 * @brief distortion of each peak; euclidean distance of peak[real value] & origin peak: dist((_x,_y) , (_xOrigin,_yOrigin))
		 */
		std::vector<double> distortion;
		/**
		 * @brief skewness of each peak; rotation of origin peak from vertical line
		 */
		std::vector<double> skewness;
		/**
		 * @brief map of <rule, number of hit>
		 */
		std::unordered_map<std::string,unsigned int> rulesHit;
	public:
		/**
		 * @brief constructor of class
		 */
		annotation(): on(-1), off(-1) { clear();}

		/**
		 * @brief destructor of class
		 */
		~annotation() {peak.clear(); flatness.clear(); distortion.clear(); skewness.clear(); rulesHit.clear();}
		/**
		 * @brief sets all parameters into default value
		 */
		void clear() {on = -1; off = -1; peak.clear(); flatness.clear(); distortion.clear(); skewness.clear(); rulesHit.clear();}
};

/**
 * @brief enum of all valid label value of a candidate
 */
enum candidateLabel {slurUnrelated = 0, slurRising = 1, slurFalling = -1, peak = 2, peakUnrelated = -2, sluredPeak = 3, sluredPeakUnrelated = -3}; // different label value of a candidate

/**
 * @brief details information of a candidate
 */
struct candidate {
	private:
		/**
		 * @brief lable of candidate contains: (0) unrelated slur(default val), (1) rising slur, (-1) falling slur, (2) peak, (-2) unrelated peak
		 */
		candidateLabel _label;

		// Y = aX + b: is a line that can be fitted to the candidate based on a pre-defined criteria of 'max slope'
		/**
		 * @brief rising slope
		 */
		double _a0;
		/**
		 * @brief rising intercept
		 */
		double _b0;
		/**
		 * @brief falling slopes
		 */
		double _a1;
		/**
		 * @brief falling intercept
		 */
		double _b1;

		/**
		 * @brief max amplitude of candidate
		 */
		double _y;
		/**
		 * @brief index (time) of max amplitude of candidate
		 */
		int _x;
		/**
		 * @brief max amplitude of imaginary candidate based on fitted line
		 */
		double _yOrigin;
		/**
		 * @brief index (time) of max amplitude of imaginary candidate based on fitted line
		 */
		int _xOrigin;
		/**
		 * @brief number of samples which shape '_y'
		 */
		double _flatnessSamples;
		/**
		 * @brief skewness of candidate; rotation of fitted lines from vertical line
		 */
		double _skewness;
		/**
		 * @brief distortion of candidate; euclidean distance of candidate from candidate's origin
		 */
		double _distortion;
		/**
		 * @brief start point and end point of a candidate
		 */
		int _risingRangeInfo[2];// [0]: start point of rising slope of candidate based on first derivative && [1]: end point of rising slope of candidate based on derivative
					//	It shows the start point of rising edge of candidate till peak of the candidate 		
		/**
		 * @brief start poin and end point of a candidate based on respective candidates
		 */
		int _candidateRangeInfo[2];// [0]: start point of candidate range based on previous candidate on a real signal && [1]: end point of candidate range based on next candidate on a real signal 
					   // 	It covers _risingRangeInfo containing rising edge/peak/falling edge. The start points gets started from last sample of end point of previous candidate
		/**
		 * @brief feature set of candidate for classifier
		 */
		std::vector<double> _featureSet;

	public:
		/**
		 * @brief constructor of class
		 */
		candidate(): _label(candidateLabel::slurUnrelated), _a0(0), _b0(0), _a1(0), _b1(0), _y(0), _x(0), _yOrigin(0), _xOrigin(0), _flatnessSamples(0), _skewness(0), _distortion(0) {} 

		// set functions
		void set_label(candidateLabel val) { _label = val;}
		void set_a0(double val) { _a0 = val;}
		void set_b0(double val) { _b0 = val;}
		void set_a1(double val) { _a1 = val;}
		void set_b1(double val) { _b1 = val;}
		void set_y(double val)  { _y = val;}
		void set_x(int val)     { _x = val;}
		void set_yOrigin(double val) { _yOrigin = val;}
		void set_xOrigin(int val)    { _xOrigin = val;}
		void set_flatnessSamples(int val) { _flatnessSamples = val;}
		void set_skewness(double val)     { _skewness = val;}
		void set_distortion(double val)   { _distortion = val;}
		void set_risingRangeInfo(int val, int index)    {index = (index >= 1) ? 1 : 0; _risingRangeInfo[index] = val;}
		void set_candidateRangeInfo(int val, int index) {index = (index >= 1) ? 1 : 0; _candidateRangeInfo[index] = val;}
		void set_featureSet(std::vector<double> val)	{_featureSet = val;}

		// get functions
		candidateLabel    get_label() const	{ return _label;}
		double get_a0() const		{ return _a0;}
		double get_b0() const		{ return _b0;}
		double get_a1() const		{ return _a1;}
		double get_b1()	const		{ return _b1;}
		double get_y()	const		{ return _y;}
		int    get_x()  const		{ return _x;}
		double get_yOrigin() const	{ return _yOrigin;}
		int    get_xOrigin() const	{ return _xOrigin;}
		double get_flatnessSamples() const { return _flatnessSamples;}
		double get_skewness()	const	{ return _skewness;}
		double get_distortion()	const	{ return _distortion;}
		int		get_risingRangeInfo(const int index) const {return _risingRangeInfo[(index >= 1) ? 1 : 0];}
		int		get_candidateRangeInfo(const int index)	const {return _candidateRangeInfo[(index >= 1) ? 1 : 0];}
		std::vector<double> get_featureSet()	const {return _featureSet;}
	};

}

/*!
 * @}
 */

}
#endif

