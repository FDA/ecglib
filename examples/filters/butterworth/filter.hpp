/**
 * @file filter.hpp
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
 * The class in this file is used for filtering a digital ECG signal applying a butterworth filter with a fixed set of coefficients
 */


#ifndef ECGLIB_FILTERS_HPP
#define ECGLIB_FILTERS_HPP 1

#include "simplebw.hpp"
#include <ecglib/ecgdata.hpp>
#include <ecglib/beat.hpp>
#include <armadillo>

#include <cmath>

#include <map>
#include <vector>
#include <algorithm>

using namespace arma;
using namespace ecglib;

namespace ecglibfilter{
	/**
	 * @brief Butter worth IIR filter
	 */
	class filter {
		public:

			/**
			 * @brief copy constructor for class
			 *
			 * @param e ecgdata
			 */
			void operator()(ecgdata &e) {
				const std::size_t nleads = e.nleads();
				vec tmp(e.nsamples());

				/**
				 * Filter coefficients for a 5th order, with cutoff 25 HZ and not a stop filter, i.e. a lowpass Butterworth filter
				 * These coefficients were obtained in Octave 4.0.1 with signal-1.3.2 package (e.g.: [b, a] = butter(5, 2*25/1000);)
				 */

				vec a = {1 , -4.491830965077046 , 8.0940554178266471 , -7.3120812801503829 , 3.3110475619883983 , -0.60111582285983844};
				vec b = {2.3409914930614002e-006 , 1.1704957465307002e-005 , 2.3409914930614004e-005 , 2.3409914930614004e-005 , 1.1704957465307002e-005 , 2.3409914930614002e-006};

				for(std::size_t i = 0; i < nleads; ++i) {
					ecglib::filtfilt(e.begin_lead(i),e.end_lead(i),tmp.begin(),a.begin(),b.begin(),a.n_elem-1);

					noisemap[i] = e.lead(i) - tmp;
					e.lead(i) = tmp;
				}
			}

			/**
			 * @brief get noise vec from leadnumber
			 *
			 * @param lnum leadnumber
			 *
			 * @return noise vector
			 */
			vec noise(leadnumber lnum) const {
				std::map<leadnumber,vec>::const_iterator i = noisemap.find(lnum);

				if(i == noisemap.end()){
					std::string line = std::string("No such lead: ") + std::to_string(static_cast<int>(lnum));
					std::cerr << line << std::endl;
					throw std::runtime_error(line);
				}else{
					return i->second;
				}
			}

		private:
			std::map<leadnumber, vec> noisemap;	/**< @brief noisemap for class */
	};

}

#endif
