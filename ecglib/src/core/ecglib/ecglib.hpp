/**
 * @file core/ecglib/ecglib.hpp
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
 * The enums in this file describes leads, annotations, etc
 */

#ifndef ECGLIB_CORE_ECGLIB_LJ_2015_12_09
#define ECGLIB_CORE_ECGLIB_LJ_2015_12_09 1

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/assert.hpp>
#include <boost/any.hpp>
#include <boost/bimap.hpp>
#include <boost/assign.hpp>
#include <boost/optional.hpp>

#include <armadillo>

#include <algorithm>

#include <vector>
#include <string>

#include <map>
#include <stdexcept>

/**
* @brief Core ecglib namespace
*/
namespace ecglib {
	/*! \addtogroup core
	 * Core ECGlib classes and functions
	 * @{
	 */

	using namespace arma;

	/**
	* @brief Constant defining the number in the pointmap of annotations expected to be for global annotations
	*/
	const int GLOBAL_LEAD=-1;

	typedef int leadnumber;						/**< @brief lead number */

	/**
	* @brief Exception used in ecglib
	*/
	struct ecglib_exception : std::logic_error {
		explicit ecglib_exception(const std::string &what) : std::logic_error(what) {
		}
	};

	typedef unsigned int sample;					/**< @brief sample number */
	typedef unsigned int timems;					/**< @brief time in ms */
        typedef boost::optional<std::string> optional_value;		/**< @brief boost optional string value */

	/**
	* @brief boost bimap type for int to string to static map enum classes to strings.
	*/
	typedef boost::bimap<int, std::string> enum_map;

        /**
        * @brief Describes annotation_type
        */
	class annotation_type
	{
	  public:

                /**
                * @brief enum for class
                */
		enum domain {PON=1,PPEAK,POFF,QON,QPEAK,RPEAK,RPPEAK,SPEAK,QOFF,TON,TPEAK,TPPEAK,TOFF,UON,UPEAK,UOFF,UNKNOWN};

		/**
		* @brief enum_map annotation_type & string representation bimap
		*/
		static const enum_map annotation_typemap;

                /**
                * @brief Class index value
                */
		domain index;

                /**
                * @brief Descriptions of enum
                */
		static constexpr const char* thevalues [16] ={"P-wave onset","P-wave peak","P-wave offset","QRS onset","Q-peak","R-wave peak","R'-wave peak","S-wave peak",
					"QRS offset","T-wave onset","T-wave peak","T-wave second peak","T-wave offset","U-wave onset","U-wave peak","U-wave offset"};

                /**
                * @brief Empty Constructor
                */
                annotation_type()
		{
			this->index = UNKNOWN;
		}

                /**
                * @brief Copy constructor for annotation_type
                *
                * @param i Other annotation_type
                */
                annotation_type(const annotation_type& i)
		{
			this->index = i.index;
		}

                /**
                * @brief Copy constructor for annotation_type
                *
                * @param i other annotation_type
                */
		annotation_type(int i)
		{
			this->index = static_cast<domain>(i);
		}

                /**
                * @brief Equal constructor for annotation_type
                *
                * @param p other annotation_type
                */
		annotation_type operator = (annotation_type ann)
		{
			this->index = ann.index;
			return *this;
		}

                /**
                * @brief Equal constructor for annotation_type
                *
                * @param i other annotation_type
                */
		annotation_type operator = (int i)
		{
			this->index = static_cast<domain>(i);
			return *this;
		}

                /**
                * @brief Boost optional type for annotation_type
                */
		typedef boost::optional<annotation_type> ann_type_optional;

                /**
                * @brief String to annotation_type conversion
                *
                * @param str string representation of annotation_type value
                *
                * @return optional annotation_type value corresponding to the string
                */
		static ann_type_optional get_by_name(const char* str)
		{
			std::string name = str;
			std::transform(name.begin(), name.end(), name.begin(), ::toupper);
			return ann_type_optional(static_cast<domain>(annotation_typemap.right.at(name.c_str())));
		}

	private:

                /**
                * @brief annotation_type to string conversion
                *
                * @param index index annotation_type
                *
                * @return const char* representation of the annotation_type
                */
		static const char* names(domain index)
		{
			int tmp = static_cast<int>(index);
			if(tmp > 16 || tmp < 1){return NULL;}
			return annotation_typemap.left.at(tmp).c_str();
		}

                /**
                * @brief annotation_type to description conversion
                *
                * @param index index annotation_type
                *
                * @return Optional string description of the annotation_type
                */
		static optional_value values(domain index)
		{
			int tmp = static_cast<int>(index);
			if(tmp > 16 || tmp < 1){return optional_value();}
			return optional_value(thevalues[tmp-1]);
		}

	public:

                /**
                * @brief Get const char* representation of annotation_type
                *
                * @return const char* representation
                */
                const char* str() const
                {
                        const char* ret = names(this->index);
                        return ret;
                }

                /**
                * @brief Get description of annotation_type
                *
                * @return const char* description
                */
		const char* value() const
		{
			const char* ret = values(this->index).get().c_str();
			return ret;
		}

                /**
                * @brief Equality operator for annotation_type
                *
                * @param a other annotation_type
                *
                * @return True if annotation_types are identical
                */
                bool operator == (const annotation_type& a) const
                {
                        return static_cast<int>(this->index) == static_cast<int>(a.index);
                }

                /**
                * @brief Equality operator for annotation_type
                *
                * @param a other annotation_type
                *
                * @return True if annotation_types are identical
                */
		bool operator == (const int a) const
		{
			return static_cast<int>(this->index) == a;
		}

                /**
                * @brief Inequality operator for annotation_type
                *
                * @param a other annotation_type
                *
                * @return True if annotation_types are different
                */
                bool operator != (const annotation_type& a) const
                {
                        return static_cast<int>(this->index) != static_cast<int>(a.index);
                }

                /**
                * @brief Inequality operator for annotation_type
                *
                * @param a other annotation_type
                *
                * @return True if annotation_types are different
                */
                bool operator != (const int a) const
                {
                        return static_cast<int>(this->index) != a;
                }

                /**
                * @brief Less-than operator for annotation_type
                *
                * @param a Other annotation_type
                *
                * @return True if this annotation_type is less than a
                */
                bool operator < (const annotation_type& a) const
                {
                        return static_cast<int>(this->index) < static_cast<int>(a.index);
                }

                /**
                * @brief Less-than operator for annotation_type
                *
                * @param a Other annotation_type
                *
                * @return True if this annotation_type is less than a
                */
                bool operator < (const int a) const
                {
                        return static_cast<int>(this->index) < a;
                }

                /**
                * @brief Greater-than operator for annotation_type
                *
                * @param a Other annotation_type
                *
                * @return True if this annotation_type is greater than a
                */
                bool operator > (const annotation_type& a) const
                {
                        return static_cast<int>(this->index) > static_cast<int>(a.index);
                }

                /**
                * @brief Greater-than operator for annotation_type
                *
                * @param a Other annotation_type
                *
                * @return True if this annotation_type is greater than a
                */
                bool operator > (const int a) const
                {
                        return static_cast<int>(this->index) > a;
                }
	};



        /**
        * @brief Describes annotation_subtype 
        */
	class annotation_subtype
	{
	  public:

                /**
                * @brief enum for class
                */
		enum domain {NONE=1,NORMAL,UNKNOWN,VPC,APB,LBBB,RBBB,UBBB,AAPB,NJPB,SPB,RONT,FVNB,AESC,NJESC,SVESC,VESC,PACED,FPN};

                /**
                * @brief Class index value
                */
                domain index;

                /**
                * @brief enum_map annotation_subtype & string representation bimap
                */
                static const enum_map annotation_subtypemap;

                /**
                * @brief Descriptions of enum
                */
		static constexpr const char* thevalues [19] = {"None","Normal","Unknown","Ventricular premature contraction","Atrial premature beat","Left bundle branch block beat",
					"Right bundle branch block beat","(unspecified) bundle branch block beat","Aberrated atrial premature beat",
					"Nodal (junctional) premature beat","Supraventricular premature or ectopic beat (atrial/nodal)",
					"R-on-T premature ventricular contraction","Fusion of ventricular and normal beat","Atrial escape beat",
					"Nodal (junctional) premature beat","Supraventricular escape beat (atrial or nodal)","Ventricular escape beat",
					"Paced beat","Fusion of paced and normal beat"};

                /**
                * @brief Empty Constructor
                */
                annotation_subtype()
		{
			this->index = UNKNOWN;
		}

                /**
                * @brief Copy constructor for annotation_subtype
                *
                * @param i Other annotation_subtype
                */
                annotation_subtype(const annotation_subtype& i)
                {
                        this->index = i.index;
                }

                /**
                * @brief Copy constructor for annotation_subtype
                *
                * @param i other annotation_subtype
                */
                annotation_subtype(int i)
                {
                        this->index = static_cast<domain>(i);
                }

                /**
                * @brief Equal constructor for annotation_subtype
                *
                * @param ann other annotation_subtype
                */
                annotation_subtype operator = (annotation_subtype ann)
                {
                        this->index = ann.index;
			return *this;
                }

                /**
                * @brief Equal constructor for annotation_subtype
                *
                * @param i other annotation_subtype
                */
                annotation_subtype operator = (int i)
                {
                        this->index = static_cast<domain>(i);
			return *this;
                }

                /**
                * @brief Boost optional type for annotation_subtype
                */
		typedef boost::optional<annotation_subtype> ann_subtype_optional;

                /**
                * @brief String to annotation_subtype conversion
                *
                * @param str string representation of annotation_subtype value
                *
                * @return optional annotation_subtype value corresponding to the string
                */
		static ann_subtype_optional get_by_name(const char* str)
		{
			std::string name = str;
			std::transform(name.begin(), name.end(), name.begin(), ::toupper);
			return ann_subtype_optional(static_cast<domain>(annotation_subtypemap.right.at(name.c_str())));
		}
	private:

                /**
                * @brief annotation_subtype to string conversion
                *
                * @param index index of annotation_subtype
                *
                * @return const char* representation of the annotation_subtype
                */
		static const char* names(domain index)
		{
			int tmp = static_cast<int>(index);
                        if(tmp > 19 || tmp < 1){return NULL;}
			return annotation_subtypemap.left.at(tmp).c_str();
		}

                /**
                * @brief annotation_subtype to description conversion
                *
                * @param index index of annotation_subtype
                *
                * @return Optional string description of the annotation_subtype
                */
		static optional_value values(domain index)
		{
			int tmp = static_cast<int>(index);
                        if(tmp > 19 || tmp < 1){return optional_value();}
                        return optional_value(thevalues[tmp-1]);
		}

	public:

                /**
                * @brief Get const char* representation of annotation_subtype
                *
                * @return const char* representation
                */
                const char* str() const
                {
                        const char* ret = names(this->index);
                        return ret;
                }

                /**
                * @brief Get description of annotation_subtype
                *
                * @return const char* description
                */
                const char* value() const
                {
                        const char* ret = values(this->index).get().c_str();
                        return ret;
                }

                /**
                * @brief Equality operator for annotation_subtypes
                *
                * @param a other annotation_subtype
                *
                * @return True if annotation_subtypes are identical
                */
		bool operator == (const annotation_subtype& a) const
                {
                        return static_cast<int>(this->index) == static_cast<int>(a.index);
                }

                /**
                * @brief Equality operator for annotation_subtypes
                *
                * @param a other annotation_subtype
                *
                * @return True if annotation_subtypes are identical
                */
                bool operator == (const int a) const
                {
                        return static_cast<int>(this->index) == a;
                }

                /**
                * @brief Inequality operator for annotation_subtypes
                *
                * @param a other annotation_subtype
                *
                * @return True if annotation_subtypes are different
                */
                bool operator != (const annotation_subtype& a) const
                {
                        return static_cast<int>(this->index) != static_cast<int>(a.index);
                }

                /**
                * @brief Inequality operator for annotation_subtypes
                *
                * @param a other annotation_subtype
                *
                * @return True if annotation_subtypes are different
                */
                bool operator != (const int a) const
                {
                        return static_cast<int>(this->index) != a;
                }

                /**
                * @brief Less-than operator for annotation_subtypes
                *
                * @param a Other annotation_subtype
                *
                * @return True if this annotation_subtype is less than a
                */
                bool operator < (const annotation_subtype& a) const
                {
                        return static_cast<int>(this->index) < static_cast<int>(a.index);
                }

                /**
                * @brief Less-than operator for annotation_subtypes
                *
                * @param a Other annotation_subtype
                *
                * @return True if this annotation_subtype is less than a
                */
                bool operator < (const int a) const
                {
                        return static_cast<int>(this->index) < a;
                }

                /**
                * @brief Greater-than operator for annotation_subtype
                *
                * @param a Other annotation_subtype
                *
                * @return True if this annotation_subtype is greater than a
                */
                bool operator > (const annotation_subtype& a) const
                {
                        return static_cast<int>(this->index) > static_cast<int>(a.index);
                }

                /**
                * @brief Greater-than operator for annotation_subtype
                *
                * @param a Other annotation_subtype
                *
                * @return True if this annotation_subtype is greater than a
                */
                bool operator > (const int a) const
                {
                        return static_cast<int>(this->index) > a;
                }
	};

        /**
        * @brief Describes ecglead
	*/
	class ecglead
	{
	  public:

                /**
                * @brief enum for class
                */
		enum domain {GLOBAL=-1,UNKNOWN1,I,II,III,AVR,AVL,AVF,V1,V2,V3,V4,V5,V6,VCGMAG,X,Y,Z,UNKNOWN2};

                /**
                * @brief Class index value
                */
                domain index;

                /**
                * @brief enum_map ecglead & string representation bimap
                */
                static const enum_map ecgleadmap;

                /**
                * @brief Descriptions of enum
                */
		static constexpr const char* thevalues[19] = {"Global","Unknown1","I","II","III","avR","avL","avF","V1","V2","V3","V4","V5","V6","VCGMAG",
								"X","Y","Z","Unknown2"};
                /**
                * @brief Empty constructor
                */
                ecglead()
		{
			this->index = UNKNOWN2;
		}

                /**
                * @brief Copy constructor for ecglead
                *
                * @param i Other ecglead
                */
                ecglead(const ecglead& i)
                {
                        this->index = i.index;
                }

                /**
                * @brief Copy constructor for ecglead
                *
                * @param i other ecglead
                */
                ecglead(int i)
                {
                        this->index = static_cast<domain>(i);
                }

                /**
                * @brief Equal constructor for ecglead
                *
                * @param lead other ecglead
                */
                ecglead operator = (ecglead lead)
                {
                        this->index = lead.index;
			return *this;
                }

                /**
                * @brief Equal constructor for ecglead
                *
                * @param i other ecglead
                */
                ecglead operator = (int i)
                {
                        this->index = static_cast<domain>(i);
			return *this;
                }

                /**
                * @brief Boost optional type for ecglead
                */
		typedef boost::optional<ecglead> ecglead_optional;

                /**
                * @brief String to ecglead conversion
                *
                * @param str string representation of ecglead value
                *
                * @return optional ecglead value corresponding to the string
                */
		static ecglead_optional get_by_name(const char* str)
		{
			std::string name = str;
			std::transform(name.begin(), name.end(), name.begin(), ::toupper);
			return ecglead_optional(static_cast<domain>(ecgleadmap.right.at(name.c_str())));
		}

	private:

                /**
                * @brief ecglead to string conversion
                *
                * @param index index of ecglead
                *
                * @return const char* representation of the ecglead
                */
		static const char* names(domain index)
                {
			int tmp = static_cast<int>(index);
                        if(tmp > 17 || tmp < -1){return NULL;}
			return ecgleadmap.left.at(tmp).c_str();
		}

                /**
                * @brief ecglead to description conversion
                *
                * @param index index of ecglead
                *
                * @return Optional string description of the ecglead
                */
		static optional_value values(domain index)
		{
			int tmp = static_cast<int>(index);
                        if(tmp > 17 || tmp < -1){return optional_value();}
                        return optional_value(thevalues[tmp+1]);
		}

	public:

                /**
                * @brief Get const char* representation of ecglead
                *
                * @return const char* representation
                */
                const char* str() const
                {
                        const char* ret = names(this->index);
                        return ret;
                }

                /**
                * @brief Get description of ecglead
                *
                * @return const char* description
                */
                const char* value() const
                {
                        const char* ret = values(this->index).get().c_str();
                        return ret;
                }

                /**
                * @brief Equality operator for ecgleads
                *
                * @param a other ecglead
                *
                * @return True if ecgleads are identical
                */
                bool operator == (const ecglead& a) const
                {
                        return static_cast<int>(this->index) == static_cast<int>(a.index);
                }

                /**
                * @brief Equality operator for ecgleads
                *
                * @param a other ecglead
                *
                * @return True if ecgleads are identical
                */
                bool operator == (const int a) const
                {
                        return static_cast<int>(this->index) == a;
                }

                /**
                * @brief Inequality operator for ecgleads
                *
                * @param a other ecglead
                *
                * @return True if ecgleads are different
                */
                bool operator != (const ecglead& a) const
                {
                        return static_cast<int>(this->index) != static_cast<int>(a.index);
                }

                /**
                * @brief Inequality operator for ecgleads
                *
                * @param a Other ecglead
                *
                * @return True if ecgleads are different
                */
                bool operator != (const int a) const
                {
                        return static_cast<int>(this->index) != a;
                }

                /**
                * @brief Less-than operator for ecgleads
                *
                * @param a Other ecglead
                *
                * @return True if this ecglead is less than a
                */
                bool operator < (const ecglead& a) const
                {
                        return static_cast<int>(this->index) < static_cast<int>(a.index);
                }

                /**
                * @brief Less-than operator for ecgleads
                *
                * @param a Other ecglead
                *
                * @return True if this ecglead is less than a
                */
                bool operator < (const int a) const
                {
                        return static_cast<int>(this->index) < a;
                }

                /**
                * @brief Greater-than operator for ecglead
                *
                * @param a Other ecglead
                *
                * @return True if this ecglead is greater than a
                */
                bool operator > (const ecglead& a) const
                {
                        return static_cast<int>(this->index) > static_cast<int>(a.index);
                }

                /**
                * @brief Greater-than operator for ecglead
                *
                * @param a Other ecglead
                * 
                * @return True if this ecglead is greater than a
                */
                bool operator > (const int a) const
                {
                        return static_cast<int>(this->index) > a;
                }

	};

	/**
	* @brief Generic ECG header struct for continuous formats, largely used by bigecg formats
	*/
	struct ecgheader{
		/**
		* @brief Filename from where the header was loaded 
		*/
		std::string filename;

		/**
		* @brief Assuming same sampling frequency for all leads (loaders will need to interpolate up/down as needed)
		*/
		double fs;

		/**
		* @brief Assuming same number of samples for all leads (loaders will need to interpolate up/down as needed)
		*/
		int nsamples;

		/**
		* @brief Number of leads
		*/
		int nleads;

		/**
		* @brief Provided for backward compatility with ENUM leads
		*/
		std::vector<ecglead> leads;

		/**
		* @brief Resolution by lead
		*/
		std::vector<int> resolution;

		/**
		* @brief offset in bytes where the ecg data starts in the file. To be use by big files (bigfoot) access 
		*/
		unsigned int startoffset;
	};

	enum class Type {String, Double, Int, Uint};		/**< @brief enum of {String, Double, Int, Uint} */

	/**
	* @brief Property class for ECG properties, ECG index properties or program_options
	*/
	struct property {
		/**
		* @brief Type
		*/
		Type type;

		/**
		* @brief Value, public as it can have any value
		*/
		boost::any value;

		/**
		* @brief Description, public as it can have any value
		*/
		std::string desc;

		/**
		* @brief Create property based on type
		*
		* @param typ Type
		*/
		property(const Type typ) : type(typ) {
		}

		/**
		* @brief Create property based on type / value
		*
		* @param typ Type
		* @param val Value
		*/
		property(const Type typ, const boost::any &val) : type(typ), value(val) {
		}

		/**
		* @brief Create property based on type / value with a description
		*
		* @param typ Type
		* @param val Value
		* @param description Description
		*/
		property(const Type typ, const boost::any &val, const std::string &description) : type(typ), value(val), desc(description) {
		}

		/**
		* @brief Operator== for property
		*
		* @param other Comparison property
		*
		* @return True if equal
		*/
		bool operator==(const property &other) const {
			// Different types are always different
			if(type != other.type) return false;

			// Perform comparison based on type
			if(type == Type::String) {
				std::string s1 = boost::any_cast<std::string>(value);
				std::string s2 = boost::any_cast<std::string>(other.value);
				return s1==s2;
			} else if(type == Type::Double) {
				double s1 = boost::any_cast<double>(value);
				double s2 = boost::any_cast<double>(other.value);
				return s1==s2;
			} else if(type == Type::Int) {
				int s1 = boost::any_cast<int>(value);
				int s2 = boost::any_cast<int>(other.value);
				return s1==s2;
			} else if(type == Type::Uint) {
				unsigned int s1 = boost::any_cast<unsigned int>(value);
				unsigned int s2 = boost::any_cast<unsigned int>(other.value);
				return s1==s2;
			}

                        // If type if something else not defined, they are not different.
			return false;
		}

		/**
		* @brief Operator!= for property
		*
		* @param other Comparison property
		*
		* @return True if not equal
		*/
		bool operator!=(const property &other) const {
			return !(*this==other);
		}
	};

	typedef std::map<std::string, property> propertymap;			/**< @brief Convenience declaration of a property map type */

	/*!
	 *@}
	 */
}

#endif
