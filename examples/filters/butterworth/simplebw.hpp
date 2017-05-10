/**
 * @file simplebw.hpp
 * @author Jose Vicente <jose.vicenteruiz@fda.hhs.gov>
 * @author Dustin C McAfee <dustin.mcafee@fda.hhs.gov
 * @author Meisam Hosseini <meisam.hosseini@fda.hhs.gov>
 * @author Lars Johannesen <lars.johannesen@fda.hhs.gov>
 *
 * @version 1.0.0
 *
 * @section LICENSE
 * this example code is in the public domain within the United States, and copyright and related rights in the work worldwide are waived through the CC0 1.0 Universal Public Domain Dedication. This example is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See DISCLAIMER section and https://www.gnu.org/licenses/gpl-faq.html for more details.
 *
 * @section DISCLAIMER
 * This software and documentation were developed by the authors in their capacities as  Oak Ridge Institute for Science and Education (ORISE) research fellows at the U.S. Food and Drug Administration (FDA).
 * FDA assumes no responsibility whatsoever for use by other parties of the Software, its source code, documentation or compiled executables, and makes no guarantees, expressed or implied, about its quality, reliability, or any other characteristic.  Further, FDA makes no representations that the use of the Software will not infringe any patent or proprietary rights of third parties.   The use of this code in no way implies endorsement by the FDA or confers any advantage in regulatory decisions.
 *
 * @section DESCRIPTION
 * The functions in this file are used to apply a butterworth filter with a fixed set of coefficients
 */

#ifndef ECGLIB_SIMPLEBW
#define ECGLIB_SIMPLEBW 1

#include <armadillo>

namespace ecglib {
	using namespace arma;
	/**
         * @brief filters using coefficients a/b
         *
         * @param signalStart Start of signal
	 * @param signalStop Stop of signal
	 * @param filteredSignal Output
         * @param a a coefficients
         * @param b b coefficients
         * @param order order of filter
         */
	template <class T,class O,class B>
	void filter(T signalStart, T signalStop, O filteredSignal, B a, B b, unsigned int order) {
		filteredSignal[0]=b[0]*signalStart[0];

		for(unsigned int i=1; i<order+1;++i) {
			filteredSignal[i]=0.0f;
			for (unsigned int j=0;j<i+1;++j) {
				filteredSignal[i]=filteredSignal[i]+b[j]*signalStart[i-j];
			}

			for (unsigned int j=0;j<i;j++) {
				filteredSignal[i]=filteredSignal[i]-a[j+1]*filteredSignal[i-j-1];
			}
		}

		unsigned int numberOfPoints = signalStop-signalStart;
		for (unsigned int i=order+1;i<numberOfPoints;++i) {
			filteredSignal[i]=0.0f;
			for (unsigned int j=0;j<order+1;j++){
				filteredSignal[i]=filteredSignal[i]+b[j]*signalStart[i-j];
			}
			for (unsigned int j=0;j<order;j++){
				filteredSignal[i]=filteredSignal[i]-a[j+1]*filteredSignal[i-j-1];
			}
		}
	}

	/*
         * @brief forward-backward filter - zerophase
         *
         * @param signalStart Start of signal
	 * @param signalStop Stop of signal
	 * @param output Output
         * @param a a coefficients
         * @param b b coefficients
         * @param order order of filter
         */
	template <class T,class O,class B>
	void filtfilt(T signalStart, T signalStop, O output, B a, B b, unsigned int order){
		int numberOfPoints = signalStop - signalStart;

		typedef typename std::iterator_traits<T>::value_type V;

		std::vector<V> tmp(numberOfPoints);

		filter(signalStart,signalStop,output,a,b,order);

		std::reverse_copy(output,output+numberOfPoints,tmp.begin());

		filter(tmp.begin(),tmp.end(),output,a,b,order);

		std::reverse(output,output+numberOfPoints);

	}
}

#endif
