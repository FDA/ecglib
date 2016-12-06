/**
 * @file core/ecglib/beat.hpp
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

#ifndef ECGLIB_CORE_BEAT_LJ_2015_12_09
#define ECGLIB_CORE_BEAT_LJ_2015_12_09 1

#include <ecglib/ecgdata.hpp>

#include <armadillo>

#include <vector>
#include <map>
#include <stdexcept>

namespace ecglib {
	/*! \addtogroup core
	 * Core ECGlib classes and functions
	 * @{
	 */

	using namespace arma;

        /**
        * @brief Describes position 
        */
	class position
	{
	  public:
                /**
                * @brief enum for class
                */
		enum domain {BEFORE=0,AFTER,NEUTRAL};

                /**
                * @brief Descriptions of enum
                */
		static constexpr const char* thevalues[3] = {"Before","After","Neutral"};

                /**
                * @brief Class index value
                */
                domain index;

                /**
                * @brief enum_map position & string representation bimap
                */
                static const enum_map positionmap;

                /**
                * @brief Empty Constructor
                */
                position()
		{
			this->index = NEUTRAL;
		}

                /**
                * @brief Copy constructor for position
                *
                * @param i Other position
                */
                position(const position& i)
                {
                        this->index = i.index;
                }

                /**
                * @brief Copy constructor for position
                *
                * @param i other position
                */
                position(int i)
                {
                        this->index = static_cast<domain>(i);
                }

                /**
                * @brief Equal constructor for position
                *
                * @param p other position
                */
                position operator = (position p)
                {
                        this->index = p.index;
			return *this;
                }

                /**
                * @brief Equal constructor for position
                *
                * @param i other position
                */
                position operator = (int i)
                {
                        this->index = static_cast<domain>(i);
			return *this;
                }

                /**
                * @brief Boost optional type for position
                */
		typedef boost::optional<position> position_optional;

                /**
                * @brief String to position conversion
                *
                * @param str string representation of position value
                *
                * @return optional position value corresponding to the string
                */
		static position_optional get_by_name(const char* str)
		{
			return position_optional(static_cast<domain>(positionmap.right.at(str)));
		}

	private:
                /**
                * @brief Position to string conversion
                *
                * @param index index of position
                *
                * @return const char* representation of the Position
                */
		static const char* names(domain index)
                {
			int tmp = static_cast<int>(index);
                        if(tmp > 2 || tmp < 0){return NULL;}
			return positionmap.left.at(tmp).c_str();
		}

                /**
                * @brief Position to description conversion
                *
                * @param index index of position
                *
                * @return Optional string description of the position
                */
		static optional_value values(domain index)
		{
			int tmp = static_cast<int>(index);
                        if(tmp > 2 || tmp < 0){return optional_value();}
                        return optional_value(thevalues[tmp]);
		}

	public:

                /**
                * @brief Get const char* representation of position
                *
                * @return const char* representation
                */
                const char* str() const
                {
                        const char* ret = names(this->index);
                        return ret;
                }

                /**
                * @brief Get description of position
                *
                * @return const char* description
                */
                const char* value() const
                {
                        const char* ret = values(this->index).get().c_str();
                        return ret;
                }

                /**
                * @brief Equality operator for positions
                *
                * @param a other position
                *
                * @return True if positions are identical
                */
                bool operator == (const position& a) const
                {
                        return static_cast<int>(this->index) == static_cast<int>(a.index);
                }

                /**
                * @brief Equality operator for positions
                *
                * @param a other position
                *
                * @return True if positions are identical
                */
                bool operator == (const int a) const
                {
                        return static_cast<int>(this->index) == a;
                }

                /**
                * @brief Inequality operator for positions
                *
                * @param a other position
                *
                * @return True if positions are different
                */
                bool operator != (const position& a) const
                {
                        return static_cast<int>(this->index) != static_cast<int>(a.index);
                }

                /**
                * @brief Inequality operator for positions
                *
                * @param a Other position
                *
                * @return True if positions are different
                */
                bool operator != (const int a) const
                {
                        return static_cast<int>(this->index) != a;
                }

                /**
                * @brief Less-than operator for positions
                *
                * @param a Other position
                *
                * @return True if this position is less than a
                */
                bool operator < (const position& a) const
                {
                        return static_cast<int>(this->index) < static_cast<int>(a.index);
                }

                /**
                * @brief Less-than operator for positions
                *
                * @param a Other position
                *
                * @return True if this position is less than a
                */
                bool operator < (const int a) const
                {
                        return static_cast<int>(this->index) < a;
                }
                /**
                * @brief greater-than operator for positions
                *
                * @param a Other position
                *
                * @return True if this position is greater than a
                */
                bool operator > (const position& a) const
                {
                        return static_cast<int>(this->index) > static_cast<int>(a.index);
                }

                /**
                * @brief greater-than operator for positions
                *
                * @param a Other position
                *
                * @return True if this position is greater than a
                */
                bool operator > (const int a) const
                {
                        return static_cast<int>(this->index) > a;
                }
	};

	/**
	* @brief Describes cardiac beats
	*/
	struct beat {
		beat() : beatlabel(ecglib::annotation_subtype::UNKNOWN) {
		}

		/**
		* @brief Start of cardiac beat
		*/
		std::size_t start;

		/**
		* @brief Stop of cardiac beat with allowed overlap to the following cardiac beat
		*/
		std::size_t stop;

		/**
		* @brief The location of the Rpeak used for synchronization
		*/
		std::size_t rpeak;

		/**
		* @brief Length of cardiac beat(with overlap)
		*/
		std::size_t length;

		/**
		* @brief RR interval of a given cardiac beat or -1 if it is the first
		*/
		int rr;

		/**
		* @brief Class id of beat
		*/
		int classid;

		/**
		* @brief Points within start and stop by lead
		*/
		ecglib::ecgdata::pointmap points;
		/**
		 * @brief Subtype used commonly to store beat labels
		 */
		ecglib::annotation_subtype beatlabel;

		/**
		* @brief Two beats are identical if they have same start/stop
		*
		* @param other Other beat
		*
		* @return True if equal
		*/
		bool operator==(const beat &other) {
			return start == other.start && stop == other.stop;
		}

		/**
		* @brief Inequality operator for beats
		*
		* @param other Other beat
		*
		* @return True if beats are different (different start/stop)
		*/
		bool operator!=(const beat &other) {
			return !(*this==other);
		}
	};

	/**
	 * @brief Map annotations to beats
	 *
	 * @param points Pointmap to use as source
	 * @param locs Beat locations
	 * @param nsamples Number of samples
	 * @param fs Sampling frequench (Hz)
	 * @param keepall Keep all beats, i.e. keep very first and very last beat
	 * @param precutwin Precut window, i.e. the left of the beat is cut using a prefix window, as the PR interval exhibits little rate dependency
	 *
	 * @return Vector of beats
	 */
	std::vector<beat> create_all_beats(ecglib::ecgdata::pointmap points, const std::vector<ecglib::annotation> &locs, int nsamples, double fs, bool keepall, int precutwin);

	/**
	 * @brief Converts beats to a pointmap
	 *
	 * @param beats Beats to use
	 *
	 * @return Pointmap
	 */
	pointmap beats_to_pointmap(const std::vector<beat> &beats);

	/**
	* @brief Compute global annotations, i.e. average Tpeak is global Tpeak
	*
	* @tparam F Globalizer
	* @param pm Input pointmap
	* @param nleads Number of leads
	* @param anntyp Annotation type to globalize
	* @param f Globalizer
	*
	* @return input pointmap that contains global annotations
	*/
	template<class F>
	pointmap make_global(const pointmap &pm, const int nleads, const ecglib::annotation_type &anntyp, const F &f) {
		uvec vals = zeros<uvec>(nleads);
		int cnt = 0;

		pointmap pmo;

		// Loop over all annotationsets (across leads)
		auto pmend = pm.end();
		for(auto pmi = pm.begin(); pmi != pmend; ++pmi) {
			if(pmi->first == ecglib::GLOBAL_LEAD){
				continue;
			}

			// Loop across all annotations from an annotationset
			auto piend = pmi->second.end();
			for(auto pi = pmi->second.begin(); pi != piend; ++pi) {
				// Store annotations of the type we are globalizing
				if(pi->second.type() == anntyp) {
					// Should not have more than nleads!
					if(cnt >= nleads){
						std::cerr << "Too many points";
						throw std::logic_error("Too many points");
					}

					vals[cnt++] = pi->first;
				}

				// Always keep the annotations
				pmo[pmi->first][pi->first] = pi->second;
			}
		}

		// If we had more than 1 annotation, create global!
		if(cnt > 1) {
			vals = vals.subvec(0,cnt-2); // subvec is to [0,n] so need to move two down as cnt is one above what it should be
			unsigned int newloc = f(vals);

			// Create and save global annotation
			ecglib::annotation newann(newloc, anntyp, ecglib::GLOBAL_LEAD);

			pmo[ecglib::GLOBAL_LEAD][newloc] = newann;
		}

		return pmo;
	}


	/**
	 * @brief Compute global annotations
	 *
	 * @tparam F Globalizer type
	 * @param beat Beat
	 * @param nleads Number of leads
	 * @param anntyp Annotation type
	 * @param f Globalizer functor
	 */
	template<class F>
	void make_global(ecglib::beat &beat, const int nleads, const ecglib::annotation_type &anntyp, const F &f) {
		beat.points = make_global(beat.points, nleads, anntyp, f);
	}

	/**
	 * @brief Compute global annotations
	 *
	 * @tparam F Globalizer type
	 * @param beats Beat vector
	 * @param nleads Number of leads
	 * @param anntyp Annotation type
	 * @param f Globalizer functor
	 */
	template<class F>
	void make_global(std::vector<ecglib::beat> &beats, const int nleads, const ecglib::annotation_type &anntyp, const F &f) {
		for(std::size_t i = 0; i < beats.size(); ++i) {
			make_global(beats[i].points, nleads, anntyp, f);
		}
	}

	/*!
	 *@}
	 */
}

#endif
