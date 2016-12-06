/**
 * @file delineators/twave/ecglib/delineator/twave/delineateFinder.hpp
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
 * Methods for finding delineators and peaks of candidates
 */
#ifndef ECGLIB_DELINEATORS_TWAVE_DELINEATEFINDER_MH_2015_12_09
#define ECGLIB_DELINEATORS_TWAVE_DELINEATEFINDER_MH_2015_12_09 1

#include <armadillo>

#include "generalstructure.hpp"

namespace ecglib {
	/*! \addtogroup delineator-twave
	 * T-wave delineator classes and functions
	 * @{
	 */

namespace twaveDelineate {

/**
 * @brief Finds the candidate delineations
 */
class delineateFinder {
	public:
		/**
		 * @brief constructor of class
		 */
		delineateFinder() {}

		/**
		 * @brief desstructor of class
		 */
		~delineateFinder() {}

		/**
		 * @brief Finds a peak of candidate based on intersection between 'bisector of slopes' and 'signal' (not based on max amplitude)
		 *
	 	 * @param wave part of candidate wave
		 * @param shift index of first point of partial wave on a candidate wave
		 * @param a0 rising slope of candidate
		 * @param b0 rising intercept of candidate
		 * @param a1 falling slope of candidate
		 * @param b1 falling intercept of candidate
		 * @param x (out var) x of peak
		 * @param y (out var) y of peak
		 * @param bisectAngle (out var) angle of bisector beetween a0 & a1
		 */
		void peakOriginFinder(const arma::rowvec& wave, int shift, double a0, double b0, double a1, double b1, int& x, double& y, double& bisectAngle);
	protected:

		/**
		 * @brief Main function for finding delineators of a candidate
		 *
	 	 * @param wave input twave 
		 * @param derivative first derivative of twave
		 * @param candid (out var) current candidate for finding delineators
		 * @param deltaAmplitude threshold for finding real peak of candidate
		 */
		void delineatorsInfo(const arma::rowvec& wave, const arma::rowvec& derivative, candidate& candid, double deltaAmplitude);

        /**
		 * @brief Linear regression [x,y]
		 *
	 	 * @param y value of points 
		 * @param leftBoundX start point of x
		 * @param rightBoundX end point of x
		 * @param a (out var) slope of linear regression
		 * @param b (out var) intercept of linear regression
		 */
		void linearRegression(const arma::rowvec& y, int leftBoundX, int rightBoundX, double& a, double& b);

		/**
		 * @brief Finds an intersection point of two lines
		 *
		 * @param a0 slope of line0
		 * @param b0 intercept of line0
		 * @param a1 slope of line1
		 * @param b1 intercept of line1
		 * @param x (out var) x of intersection
		 * @param y (out var) y of intersection
		 * @param angle (out var) bisector angle between two lines
		 */
		void intersectionLine(double a0, double b0, double a1, double b1, double& a, double& y, double& angle);

		/**
		 * @brief Finds a peak of candidate based on max amplitude
		 *
		 * @param wave partial twave
		 * @param deltaAmplitude delta value for fining a peak based on max amplitude
		 * @param x (out var) x of peak
		 * @param y (out var) y of peak
		 * @param flatness (out var) number of points that have an amplitude <= |max amplitude - deltaAmplitude| and the middle point of them gets shown by x (peak of candidate) and median of them shown by y 
		 */
		void peakFinder(const arma::rowvec& wave, double deltaAmplitude, int& x, double& y, int& flatness);
	};
}

/*!
 * @}
 */

}
#endif

