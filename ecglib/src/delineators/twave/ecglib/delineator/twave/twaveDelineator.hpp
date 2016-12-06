/**
 * @file delineators/twave/ecglib/delineator/twave/twaveDelineator.hpp
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
 * General input/output data to/form twaveDelineator
 */

#ifndef ECGLIB_DELINEATORS_TWAVE_TWAVEDELINEATOR_MH_2015_12_09
#define ECGLIB_DELINEATORS_TWAVE_TWAVEDELINEATOR_MH_2015_12_09 1

#include <tuple>
#include <ecglib/delineator/twave/delineate.hpp>
#include <ecglib/ecglib.hpp>
#include <ecglib/ecgdata.hpp>
#include <ecglib/annotation.hpp>
#include <ecglib/util/config.hpp>

namespace ecglib {
	/*! \addtogroup delineator-twave
	 * T-wave delineator classes and functions
	 * @{
	 */

	/**
	 * @brief Configuration class of twave delineator which inherited from config class
	 */
	class twaveDelineator_config : public config {
		public:
			/**
			 * @brief Construction of twaveDelineator_config
			 */
			twaveDelineator_config(): config() {
				defaults();
			}

		protected:
			/**
			 * @brief Default configuration of twaveDelineator
			 */
			void defaults();
	};

	/**
	 * @brief Main entrance into twaveDelineator
	 *
 	 * @param e Input ecg data
 	 * @param pmin Pointmap to use as source
	 * @param cfg
	 *
	 * @return tuple<pointmap, annotation>: pointmap as output and annotation contains twaveDelineator
	 */
	std::tuple<pointmap, ecglib::twaveDelineate::annotation> twaveDelineators(const ecglib::ecgdata &e, const ecglib::pointmap &pmin, const ecglib::twaveDelineator_config &cfg);

	/*! 
	 * @}
	 */
}

#endif
