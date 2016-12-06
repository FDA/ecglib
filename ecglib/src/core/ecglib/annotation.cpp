/**
 * @file core/ecglib/annotation.cpp
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
 * This file has the annotation class definition
 */
#include "ecgdata.hpp"

#include "annotation.hpp"

namespace ecglib {
	void get_annotations(const ecgdata::annotationset &pm, const ecglib::annotation_type &typ, std::vector<annotation> &out) {
		typedef ecglib::ecgdata::const_pointiterator piter;
		piter pend = pm.end();
		for(piter pi = pm.begin(); pi != pend; ++pi) {
			if(pi->second.type() == typ) {
				out.push_back(pi->second);
			}
		}
	}

	void get_annotations(const ecgdata::pointmap &pm, const ecglib::annotation_type &typ, std::vector<annotation> &out) {
		typedef ecglib::ecgdata::const_pointmapiterator pmiter;

		pmiter pmend = pm.end();

		for(pmiter pmi = pm.begin(); pmi != pmend; ++pmi) {
			typedef ecglib::ecgdata::const_pointiterator piter;

			piter pend = pmi->second.end();
			for(piter pi = pmi->second.begin(); pi != pend; ++pi) {
				if(pi->second.type() == typ) {
					out.push_back(pi->second);
				}
			}
		}
	}

	ecgdata::annotationset rebase(const ecgdata::annotationset &ann, const int newstart) {
		typedef ecgdata::const_pointiterator piter;

		piter pend = ann.end();

		ecgdata::annotationset rebased;

		for(piter pi = ann.begin(); pi != pend; ++pi) {
			ecglib::annotation newann = pi->second;
			newann.location(newann.location()-newstart);

			rebased[newann.location()] = newann;
		}

		return rebased;
	}

	void rebase(ecgdata::pointmap &pm, const int newstart) {
		typedef ecglib::ecgdata::pointmapiterator pmiter;
	
		pmiter pmend = pm.end();
		for(pmiter pmi = pm.begin(); pmi != pmend; ++pmi) {
			pmi->second = rebase(pmi->second, newstart);
		}
	}

	void get_annotations(const pointmap &pm, const leadnumber l, const annotation_type &typ, std::vector<annotation> &anns) {
		auto pmi = pm.find(l);
		if (pmi == pm.end()) return;

		for(auto &anniter : pmi->second) {
			if(anniter.second.type() == typ) {
				anns.push_back(anniter.second);
			}
		}
	}
}
