/**
 * @file delineators/twave/ecglib/delineator/twave/delineateFinder.cpp
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

#include <cmath>

#include "delineateFinder.hpp"

using namespace ecglib::twaveDelineate;

/// main function for finding delineators of a candidate
void delineateFinder::delineatorsInfo(const arma::rowvec& wave, const arma::rowvec& derivative, candidate& candid, double deltaAmplitude) {

	int waveStart = candid.get_candidateRangeInfo(0);
	int waveEnd = candid.get_candidateRangeInfo(1);

	int risingRangeStart = candid.get_risingRangeInfo(0) - waveStart;
	int risingRangeEnd = candid.get_risingRangeInfo(1) - waveStart;
	
	arma::rowvec waveCandidate = wave(arma::span(waveStart, waveEnd));
	arma::rowvec derivativeCandidate = derivative(arma::span(waveStart, waveEnd-1));

	// for fining slope/intercept is better to do a light filtering such as loose median window filtering. 
	// pre_assumption is, the signal has got filterd before this step and consequently the filtering has not implied for improving the time complexity

	/* step 01: finds rising & falling slopes/intercepts of a candidate */ 
	arma::uvec maxSlopeIndex = arma::find( derivativeCandidate(arma::span(0, risingRangeEnd)) == arma::max(derivativeCandidate(arma::span(0, risingRangeEnd))), 1);
	arma::uvec minSlopeIndex = risingRangeEnd + arma::find( derivativeCandidate(arma::span(risingRangeEnd, derivativeCandidate.n_elem -1)) == arma::min(derivativeCandidate(arma::span(risingRangeEnd, derivativeCandidate.n_elem -1))), 1);
	int window = 5; // number of window samples for finding intercept

		/* step 01-1: rising line */
		int windowLowerBound0 = std::max(0, static_cast<int>(maxSlopeIndex(0)) - window);
		int windowUpperBound0 = std::min(static_cast<int>(waveCandidate.n_elem) -1, static_cast<int>(maxSlopeIndex(0)) + window);
		double a0(0), b0(0);
		delineateFinder::linearRegression(waveCandidate(arma::span(windowLowerBound0,windowUpperBound0)), windowLowerBound0, windowUpperBound0, a0, b0);

		/* step 01-2: falling line */
		int windowLowerBound1 = std::max(0, static_cast<int>(minSlopeIndex(0)) - window);
		int windowUpperBound1 = std::min(static_cast<int>(waveCandidate.n_elem) -1, static_cast<int>(minSlopeIndex(0)) + window);
		double a1(0), b1(0);
		delineateFinder::linearRegression(waveCandidate(arma::span(windowLowerBound1,windowUpperBound1)), windowLowerBound1, windowUpperBound1, a1, b1);

	/* step 02: finds peak of candidate based on intersection between signal and bisector of slopes.
		    this peak is an imaginery peak that shows the place of original peak before some distortion based on regression lines*/
	int xIntersect(0);
	double yIntersect(0), angleIntersect(0);
	delineateFinder::peakOriginFinder(waveCandidate(arma::span(arma::as_scalar(maxSlopeIndex(0)), arma::as_scalar(minSlopeIndex(0)))), maxSlopeIndex(0), a0, b0, a1, b1, xIntersect, yIntersect, angleIntersect);

	/* step 03: finds a peak of candidate based on max amplitude.
		    this peak is a median of points with high amplitude. 
		    deltaAmplitude uses for a range that median should be calculated by.
		    flatness shows the number of points with high amplitude which are finded as peak */
	int x(0), flatness(0);
	double y(0);
	delineateFinder::peakFinder(waveCandidate(arma::span(risingRangeStart, risingRangeEnd)), deltaAmplitude, x, y, flatness);
	x += risingRangeStart; // update x based on partial wave

	/* step 04: preparation of candidate delineators*/
	candid.set_a0(a0);
	candid.set_b0(b0 - (a0*waveStart)); // recalculated based on waveStart

	candid.set_a1(a1);
	candid.set_b1(b1 - (a1*waveStart)); // recalculated based on waveStart

	candid.set_x(x + waveStart); 	    // recalculated based on waveStart
	candid.set_y(y);

	candid.set_xOrigin(xIntersect + waveStart); // recalculated based on waveStart
	candid.set_yOrigin(yIntersect);

	candid.set_flatnessSamples(flatness);

	candid.set_skewness(angleIntersect);

	double peakDistX = std::pow(candid.get_x() - candid.get_xOrigin(), 2.);
	double peakDistY = std::pow(candid.get_y() - candid.get_yOrigin(), 2.);
	candid.set_distortion(std::sqrt(peakDistX + peakDistY));
}

/// linear regression
void delineateFinder::linearRegression(const arma::rowvec& y, int leftBoundX, int rightBoundX, double& a, double& b){
    	// y = ax + b
    
	int n = y.n_elem;
    	arma::rowvec x = arma::linspace<arma::rowvec>(leftBoundX, rightBoundX, rightBoundX - leftBoundX +1);
	double divisor = (n * arma::accu(arma::pow(x,2.))) - std::pow(arma::accu(x),2.);
    	a = 0; b = 0;
	if (std::abs(divisor) > 1e-4) {
		a = ( (n* arma::accu(x%y)) - (arma::accu(x) * arma::accu(y)) )  / divisor;
		b = ( (arma::accu(arma::pow(x,2))*arma::accu(y)) - (arma::accu(x)*arma::accu(x%y)) )  / divisor;
	}
}

/// find a peak based on intersection of bisector between 'rising & falling slopes' and signal
void delineateFinder::peakOriginFinder(const arma::rowvec& wave, int shift, double a0, double b0, double a1, double b1, int& x, double& y, double& bisectAngle) {

	double aBisect(0), bBisect(0);
	delineateFinder::intersectionLine(a0, b0, a1, b1, aBisect, bBisect, bisectAngle);

	arma::rowvec rangeX = arma::linspace<arma::rowvec>(1, wave.n_elem, wave.n_elem);
	arma::rowvec rangeY = aBisect * (rangeX + shift) + bBisect; // calculate y line based on shifted X
	arma::rowvec rangeWave = wave;
	arma::rowvec deltaError = arma::pow(rangeY - rangeWave, 2.);
    	arma::uvec xIntersectionIndex = arma::find(deltaError == deltaError.min(), 1);

	// x is the nearest point of intersection bisector and signal to the signal
	x = std::min(arma::as_scalar(rangeX(xIntersectionIndex)), static_cast<double>(wave.n_elem -1));
	y = wave(x);

	x += shift; // update x based on start point of partial wave 
}

/// finds an intersection point of two lines
void delineateFinder::intersectionLine(double a0, double b0, double a1, double b1, double& a, double& b, double& angle) {

	// Intersection of two lines
	// y = ax + b
	// y = cx + d
	// => ax + b = cx + d => x = (d-b)/(a-c)

	double d_b = b1 - b0;
	double a_c = a0 - a1;

	double x = (a_c != 0) ? d_b/a_c : 0;
	double y = a0*x + b0;

	const double pi = 3.1416;
	double at0 = atan(a0)*180/pi;
	double at1 = atan(a1)*180/pi;

	double bisector = 90 + (at0+at1)/2;
	if (bisector == 90) bisector = atan(100)*180/pi; // ~ 89.5 degree
	
	angle = 89.5 - bisector;
	a = tan(bisector*pi/180);
	b = y - a*x;
}

/// finds a peak of candidate based on max amplitude
void delineateFinder::peakFinder(const arma::rowvec& wave, double deltaAmplitude, int& x, double& y, int& flatness) {

	// deltaAmplitude: is used for finding a peak. peak is not just a point it is median of couple of points with highest amplitude

    	arma::uvec indexPeaks = arma::find( wave >= (wave.max() - deltaAmplitude)); // find all points with high amplitude: max amplitude - deltaAmplitude
	int indexPeakMid = std::round(static_cast<double>(indexPeaks.n_elem)/2.);
	x = indexPeaks(std::max(0,indexPeakMid-1));
	y = wave(x);
	flatness = indexPeaks.n_elem;
}
