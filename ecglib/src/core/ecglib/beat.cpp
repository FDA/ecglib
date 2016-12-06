/**
 * @file core/ecglib/beat.cpp
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
 * The classes in this file are used for describing beats
 */

#include <ecglib/annotation.hpp>
#include <ecglib/beat.hpp>
#include <ecglib/detail/detail.hpp>

#include <cassert>

namespace ecglib { 

	//Definitions for annotation_type, annotation_subtype, ecglead, and position classes.
	constexpr const char* annotation_type::thevalues[16];
	constexpr const char* annotation_subtype::thevalues[19];
	constexpr const char* ecglead::thevalues[19];
	constexpr const char* position::thevalues[3];

        const enum_map annotation_type::annotation_typemap = boost::assign::list_of<enum_map::relation>
	        (1,"PON") (2,"PPEAK") (3,"POFF") (4,"QON")
                (5,"QPEAK") (6,"RPEAK") (7,"RPPEAK") (8,"SPEAK")
                (9,"QOFF") (10,"TON") (11,"TPEAK") (12,"TPPEAK")
                (13,"TOFF") (14,"UON") (15,"UPEAK") (16,"UOFF");

        const enum_map annotation_subtype::annotation_subtypemap = boost::assign::list_of<enum_map::relation>
	        (1,"NONE") (2,"NORMAL") (3,"UNKNOWN") (4,"VPC")
                (5,"APB") (6,"LBBB") (7,"RBBB") (8,"UBBB")
                (9,"AAPB") (10,"NJPB") (11,"SPB") (12,"RONT")
                (13,"FVNB") (14,"AESC") (15,"NJESC") (16,"SVESC")
                (17,"VESC") (18,"PACED") (19,"FPN");

        const enum_map ecglead::ecgleadmap = boost::assign::list_of<enum_map::relation>
	        (-1,"GLOBAL") (0,"UNKNOWN1") (1,"I") (2,"II")
                (3,"III") (4,"AVR") (5,"AVL") (6,"AVF")
                (7,"V1") (8,"V2") (9,"V3") (10,"V4")
                (11,"V5") (12,"V6") (13,"VCGMAG") (14,"X")
                (15,"Y") (16,"Z") (17,"UNKNOWN2");

        const enum_map position::positionmap = boost::assign::list_of<enum_map::relation>
		(0,"BEFORE") (1,"AFTER") (2,"NEUTRAL");

	using namespace arma;

	std::vector<beat> create_all_beats(ecglib::ecgdata::pointmap points, const std::vector<ecglib::annotation> &locs, int nsamples, double fs, bool keepall, int precutwin) {
		std::vector<beat> beats;

		const unsigned int precut = static_cast<int>(round(precutwin*(fs/1000.0)));

		typedef ecglib::ecgdata::pointmapiterator citer;
		typedef ecglib::ecgdata::pointiterator piter;

		if(locs.size() == 1) {
			beat b;

			b.beatlabel = locs[0].subtype();
			b.start = 0;
			b.stop = nsamples;
			b.rpeak = locs[0].location();
			b.rr = -1;
			b.length = b.stop - b.start;

			beats.push_back(b);
		} else {
			std::size_t start = 0, stop = locs.size();
			if(!keepall) {
				// Find first beat that has a start within the signal
				for(std::size_t i = 1; i < locs.size()-1; ++i) {
					if (locs[i] > precut) {
						start = i;
						break;
					}
				}

				// Exclude the last beat as the stop is not defined
				stop = locs.size()-1;

				if(stop == start) {
					std::cerr << "No extractable beats";
					throw ecglib_exception("No extractable beats");
				}
			}

			for(std::size_t i = start; i < stop; ++i) {
				beat b;
				b.beatlabel = locs[i].subtype();

				b.rpeak = locs[i];	
				b.start = locs[i] - precut;
				if(locs[i] < precut) {
					b.start = 0;
				}

				if(i == locs.size()-1) {
					b.stop = nsamples;
				} else {
					b.stop = locs[i+1];
				}

				b.length = b.stop - b.start;

				if(i == 0) {
					b.rr = -1;
				} else {
					b.rr = locs[i] - locs[i-1];
				}

				beats.push_back(b);
			}
		}

		// Foreach lead with annotations
		citer endp = points.end();
		for(citer iterp = points.begin(); iterp != endp; ++iterp) {
			// Foreach annotation in a given lead
			piter pi = iterp->second.begin();
			piter pend = iterp->second.end();

			while(pi != pend) {
				int bnum = 0;
				bool found = false;
				unsigned int lowdiff = 0;

				// Foreach beat - find best suiter
				for(std::size_t i = 0; i < beats.size(); ++i) {
					// Is annotation inside a given beat ?
					if(pi->first > beats[i].start && pi->first < beats[i].stop) {
						std::size_t rpeak = beats[i].rpeak;

						typedef std::map<annotation_type, detail::position>::const_iterator pmap_iter;
						pmap_iter pmi = detail::position_map.find(pi->second.type());
						// Make sure annotation is right (ie. p-waves  before rpeak, t/u-waves after rpeak)
						bool ok = false;
						unsigned int curdiff;
						if(pmi != detail::position_map.end()) {
							if(pmi->second == detail::position::BEFORE) {
								if(pi->first < rpeak) {
									ok = true;
									curdiff = beats[i].rpeak - pi->first;
								}
							} else if(pmi->second == detail::position::AFTER) {
								if(pi->first > rpeak) {
									ok = true;
									curdiff = pi->first - beats[i].rpeak;
								}
							} else {
								ok = true;
								curdiff = abs(pi->first - beats[i].rpeak);
							}
						}

						if (ok) {
							if(!found) {
								found = true;
								lowdiff = curdiff;
								bnum = i;
							} else if(curdiff < lowdiff) {
								lowdiff = curdiff;
								bnum = i;
							}
						}
					}
				} // End foreach beat
				
				if(found) {
					beats[bnum].points[iterp->first][pi->first] = pi->second;
					piter pio = pi;
					++pi;
					iterp->second.erase(pio->first);
				} else {
					++pi;
				}
			} // End foreach annotation for a given lead
		} // End foreach lead with annotations

		return beats;
	}

	pointmap beats_to_pointmap(const std::vector<beat> &beats) {
		pointmap pm;

		for(std::size_t i = 0; i < beats.size(); ++i) {
			for(auto &annsetiter : beats[i].points) {
				for(auto &anniter : annsetiter.second) {
					pm[annsetiter.first][anniter.first] = anniter.second;

					if(annsetiter.first == GLOBAL_LEAD && anniter.second.type() == annotation_type::RPEAK) {
						pm[annsetiter.first][anniter.first].subtype(beats[i].beatlabel);
					}
				}
			}
		}

		return pm;
	}
}
