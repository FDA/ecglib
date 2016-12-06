/**
 * @file core/ecglib/util/util.hpp
 * @author Lars Johannesen <lars.johannesen@fda.hhs.gov>
 * @author Jose Vicente <jose.vicenteruiz@fda.hhs.gov>
 * @author Meisam Hosseini <meisam.hosseini@fda.hhs.gov>
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
 * This file defines ecglib samples to/from assuming 1000 Hz sampling rate
 **/

#ifndef ECGLIB_CORE_UTIL_LJ_2015_12_09
#define ECGLIB_CORE_UTIL_LJ_2015_12_09 1

#include <ecglib/ecglib.hpp>

namespace ecglib { 
	/*! \addtogroup core
	 * Core ECGlib classes and functions
	 * @{
	 */

	/**
	 * @brief Converts sample to time
	 *
	 * @param sample Sample number
	 * @param fs Sampling frequency
	 *
	 * @return Time in ms
	 */
	inline timems sample_to_time(const double sample, const double fs) {
		return static_cast<unsigned int>(round(1000./fs * sample));
	}

	/**
	 * @brief Converts time to sample
	 *
	 * @param time Time in ms
	 * @param fs Sampling rate
	 *
	 * @return Sample
	 */
	inline sample time_to_sample(const timems time, const double fs) {
		return static_cast<unsigned int>(round(fs/1000. * time));
	}

	/*!
	 *@}
	 */
}

#endif
