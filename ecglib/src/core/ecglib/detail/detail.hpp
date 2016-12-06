/**
 * @file core/ecglib/detail/detail.hpp
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
 * The classes in this file are used for describing beats (detail)
 */

#ifndef ECGLIB_CORE_DETAIL_LJ_2015_12_09
#define ECGLIB_CORE_DETAIL_LJ_2015_12_09 1

#include <ecglib/ecgdata.hpp>
#include <ecglib/beat.hpp>

#include <armadillo>

#include <vector>
#include <map>
#include <stdexcept>

namespace ecglib { namespace detail {
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
	 * @brief Rules for matching annotations to beats, which annotations should occur when relative to R-peak
	 *
	 * @return Rulemap
	 */
	const std::map<annotation_type, position> position_map = {
		{annotation_type::PON, position::BEFORE},
		{annotation_type::PPEAK, position::BEFORE},
		{annotation_type::POFF, position::BEFORE},
		{annotation_type::QON, position::NEUTRAL},
		{annotation_type::QPEAK, position::NEUTRAL},
		{annotation_type::RPEAK, position::NEUTRAL},
		{annotation_type::RPPEAK, position::NEUTRAL},
		{annotation_type::SPEAK, position::NEUTRAL},
		{annotation_type::QOFF, position::NEUTRAL},
		{annotation_type::TON, position::AFTER},
		{annotation_type::TPEAK, position::AFTER},
		{annotation_type::TPPEAK, position::AFTER},
		{annotation_type::TOFF, position::AFTER},
		{annotation_type::UON, position::AFTER},
		{annotation_type::UPEAK, position::AFTER},
		{annotation_type::UOFF, position::AFTER}
	};
}}

#endif
