/**
 * @file core/ecglib/annotation.hpp
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

#ifndef ECGLIB_CORE_ANNOTATION_LJ_2015_12_09
#define ECGLIB_CORE_ANNOTATION_LJ_2015_12_09 1

#include <ecglib/util/util.hpp>
#include <ecglib/ecglib.hpp>

#include <algorithm>
#include <vector>
#include <string>
#include <map>

#include <armadillo>

namespace ecglib {
	/*! \addtogroup core
	 * Core ECGlib classes and functions
	 * @{
	 */

	using namespace arma;

	/**
	* @brief ECG annotation
	*/
	class annotation {
		// Constructors
		public:

			/**
			* @brief Empty annotation
			*/
			annotation() {}

			/**
			* @brief Create annotation in global lead
			*
			* @param location Location of annotation
			* @param type Annotation type
			*/
			annotation(timems location, const annotation_type type) : _location(location), _type(type), _subtype(annotation_subtype::NONE), _lead(GLOBAL_LEAD) {}

			/**
			* @brief Create an annotation in a specific lead
			*
			* @param location Location of annotation
			* @param type Annotation type
			* @param l Annotation lead
			*/
			annotation(timems location, const annotation_type type, const leadnumber l) : _location(location), _type(type), _subtype(annotation_subtype::NONE), _lead(l) {}

			/**
			* @brief Create an annotation in a specific lead
			*
			* @param location Location of annotation
			* @param type Annotation type
			* @param l Annotation lead
			* @param sub Subannotation type
			*/
			annotation(timems location, const annotation_type type, const leadnumber l, const annotation_subtype sub) : _location(location), _type(type), _subtype(sub), _lead(l) {}

		// Methods:
		public:
			/**
			* @brief comparison to other annotation (less than)
			*
			* @param other Other annotation
			*
			* @return True if less than other annotation
			*/
			bool operator<(const annotation &other) const {
				return _location < other._location;
			}

			/**
			* @brief comparison to other annotation (greater than)
			*
			* @param other Other annotation
			*
			* @return True if greater than other annotation
			*/
			bool operator>(const annotation &other) const {
				return _location > other._location;
			}

			/**
			* @brief comparison to other time point (less than)
			*
			* @param other Other time point
			*
			* @return True if less than other time point
			*/
			bool operator<(const timems &other) const {
				return _location < other;
			}

			/**
			* @brief comparison to other time point (greater than)
			*
			* @param other Other time point
			*
			* @return True if greater than other time point
			*/
			bool operator>(const timems &other) const {
				return _location > other;
			}

			/**
			* @brief Decrement operation
			*
			* @param subt Amount to decrement by
			*
			* @return Updated annotation
			*/
			annotation operator-(const timems subt) const {
				return annotation(_location-subt,_type,_lead);
			}

			/**
			* @brief Cast to timems
			*
			* @return Location
			*/
			operator timems() const {
				return _location;
			}

		// Setters/getters:
		public:

			/**
			* @brief Get ecglead of annotation
			*
			* @return ecglead
			*/
			leadnumber lead() const {
				return _lead;
			}

			/**
			* @brief Set ecglead for annotation
			*
			* @param lead ecglead
			*/
			void lead(leadnumber lead) {
				_lead = lead;
			}

			/**
			* @brief Get location of annotation
			*
			* @return location location of annotation
			*/
			timems location() const {
				return _location;
			}

			/**
			* @brief Set location of annotation
			*
			* @param location location to set
			*/
			void location(timems location) {
				_location = location;
			}

			/**
			* @brief Get type of annotation
			*
			* @return Annotation type
			*/
			annotation_type type() const {
				return _type;
			}

			/**
			* @brief Get annotation type
			*
			* @param type Annotation type
			*/
			void type(annotation_type type) {
				_type = type;
			}

			/**
			* @brief Get annotation subtype (used to represent beat labels eg.)
			*
			* @return annotation subtype
			*/
			annotation_subtype subtype() const {
				return _subtype;
			}

			/**
			* @brief Set annotation subtype
			*
			* @param subtype Annotation subtype
			*/
			void subtype(annotation_subtype subtype) {
				_subtype = subtype;
			}

		// Helpers:
		private:

		// Attributes:
		protected:
			/**
			 * @brief Location in ms
			 */
			timems _location;

			/**
			* @brief Annotation type
			*/
			annotation_type _type;

			/**
			* @brief Annotation subtype
			*/
			annotation_subtype _subtype;

			/**
			* @brief Lead number
			*/
			leadnumber _lead;
	};


	/**
	 * @brief Annotationset class - contains annotations for a lead
	 */
	class annotationset {
		public:
			typedef std::map<timems, annotation>::iterator iterator;			/** @brief iterator for annotationset */
			typedef std::map<timems, annotation>::const_iterator const_iterator;		/** @brief constant iterator for annotationset */

		public:
			/**
			* @brief empty constructor for annotationset
			*/
			annotationset() {
			}

		public:
			/**
			* @brief Retrieve annotation at location idx
			*
			* @param idx Time location in ms
			*
			* @return Annotation copy
			*/
			annotation& operator[](const timems idx) {
				return _annset[idx];
			}

			/**
			* @brief Retrieve reference to location
			*
			* @param idx Time location in ms
			*
			* @return Constant reference to annotation 
			*/
			const annotation& operator[](const timems idx) const {
				const_iterator iter = _annset.find(idx);

				if(iter == _annset.end()) {
					std::stringstream ss;
					ss << "No annotation at: " << idx;
					std::cerr << ss.str();
					throw std::logic_error(ss.str());
				}

				return iter->second;
			}

			/**
			* @brief Iterator to start of map
			*
			* @return Iterator to start of map
			*/
			iterator begin() {
				return _annset.begin();
			}

			/**
			* @brief Iterator to end of map
			*
			* @return Iterator to end of map
			*/
			iterator end() {
				return _annset.end();
			}

			/**
			* @brief Return iterator pointing to first element containing idx
			*
			* @param idx time location for lower bound
			*
			* @return iterator pointing to first element including idx
			*/
			iterator lower_bound(const timems idx) {
				return _annset.lower_bound(idx);
			}

			/**
			* @brief Return iterator pointing to one element outside idx, so [start,end)
			*
			* @param idx time location for upper bound
			*
			* @return iterator pointing to second element including idx
			*/
			iterator upper_bound(const timems idx) {
				return _annset.upper_bound(idx);
			}

			/**
			* @brief Const correct begin
			*
			* @return Start iterator
			*/
			const_iterator begin() const {
				return _annset.begin();
			}

			/**
			* @brief Const correct end
			*
			* @return Stop iterator
			*/
			const_iterator end() const {
				return _annset.end();
			}

			/**
			* @brief Const correct lower_bound
			*
			* @param idx time location for lower bound
			*
			* @return Lower bound const iterator
			*/
			const_iterator lower_bound(const timems idx) const {
				return _annset.lower_bound(idx);
			}

			/**
			* @brief Const correct upper_bound
			*
			* @param idx time location for upper bound
			*
			* @return Upper bound const iterator
			*/
			const_iterator upper_bound(const timems idx) const {
				return _annset.upper_bound(idx);
			}

			/**
			* @brief Find annotation based on location
			*
			* @param idx Location in ms
			*
			* @return Iterator pointing to element
			*/
			iterator find(timems idx) {
				return _annset.find(idx);
			}

			/**
			* @brief Find annotation based on location
			*
			* @param idx Location in ms
			*
			* @return Const iterator to annotation
			*/
			const_iterator find(timems idx) const {
				return _annset.find(idx);
			}

			/**
			* @brief Size of map
			*
			* @return Size of map
			*/
			std::size_t size() const {
				return _annset.size();
			}

			/**
			* @brief Clear map
			*/
			void clear() {
				_annset.clear();
			}

			/**
			* @brief Erase annotation at location idx
			*
			* @param idx Location idx
			*
			* @return 1 if successful, otherwise 0
			*/
			std::size_t erase(timems idx) {
				iterator it = _annset.find(idx);

				if(it == _annset.end()) {
					return 0;
				}

				erase(it);

				return 1;
			}

			/**
			* @brief Erase annotation pointed to by iterator it
			*
			* @param it Iterator to be deleted
			*
			* @return True if successful
			*/
			bool erase(iterator it) {
				_annset.erase(it);

				return 1;
			}

			/**
			* @brief Remove annotation if PRED is true
			*
			* @tparam PRED Predicate function
			* @param pred Predicate function
			*
			* @return Iterator pointing to new end
			*/
			template<typename PRED>
			iterator remove_if(PRED &pred) {
				iterator it = begin();
				iterator endit = end();

				while(it != endit) {
					if(pred(*it)) {
						erase(it++);
					} else {
						++it;
					}
				}

				return it;
			}

			/**
			* @brief Insert from another map container of annotations
			*
			* @tparam INPUT_ITER Iterator type for other container
			* @param a1 Start iterator
			* @param a2 Stop iterator
			*
			* @return Iterator pointing at the inserted item
			*/
			template<typename INPUT_ITER>
			iterator insert(INPUT_ITER a1, INPUT_ITER a2) {
				iterator ret;

				while(a1 != a2) {
					ret = _annset.insert(*a1).first;

					++a1;
				}

				return ret;
			}

		private:
			/**
			* @brief Annotation map
			*/
			std::map<timems, annotation> _annset;
	};

	/**
	 * @brief Pointmap, links annotationsets with ecgleads
	 */
	class pointmap {
		public:
			typedef std::map<leadnumber, annotationset>::iterator iterator;			/**< @brief iterator for pointmap */
			typedef std::map<leadnumber, annotationset>::const_iterator const_iterator;	/**< @brief const iterator for pointmap */

		public:
			/**
			* @brief empty constructor for pointmap
			*/
			pointmap() {
			}

		public:
			/**
			* @brief Retrieve annotation set for lead l
			*
			* @param l Leadnumber l
			*
			* @return Annotationset reference
			*/
			annotationset& operator[](const leadnumber l) {
				return _pm[l];
			}

			/**
			* @brief Retrieve annotation set for lead l (const correct)
			*
			* @param l Leadnumber l
			*
			* @return Annotationset reference
			*/
			const annotationset& operator[](const leadnumber l) const {
				const_iterator iter = _pm.find(l);

				if(iter == _pm.end()) {
					std::stringstream ss;
					ss << "No such lead number: " << l;
					std::cerr << ss.str();
					throw std::logic_error(ss.str());
				}

				return iter->second;
			}

			/**
			* @brief Beginning of pointmap
			*
			* @return Iterator for beginning
			*/
			iterator begin() {
				return _pm.begin();
			}

			/**
			* @brief End of pointmap
			*
			* @return Iterator for end
			*/
			iterator end() {
				return _pm.end();
			}

			/**
			* @brief Beginning of pointmap (const correct)
			*
			* @return Iterator for beginning
			*/
			const_iterator begin() const {
				return _pm.begin();
			}

			/**
			* @brief End of pointmap (const correct)
			*
			* @return Iterator for end
			*/
			const_iterator end() const {
				return _pm.end();
			}

			/**
			* @brief Number of leads in pointmap
			*
			* @return Number of leads
			*/
			std::size_t size() const {
				return _pm.size();
			}

			/**
			* @brief Find annotations for lead l
			*
			* @param l Lead l
			*
			* @return Iterator for lead l
			*/
			iterator find(ecglib::leadnumber l) {
				return _pm.find(l);
			}

			/**
			 * @brief Find annotations for lead l (const correct)
			 *
			 * @param l Lead l
			 *
			 * @return Iterator for lead l
			 */
			const_iterator find(ecglib::leadnumber l) const {
				return _pm.find(l);
			}

			/**
			* @brief Clear pointmap
			*/
			void clear() {
				_pm.clear();
			}

			/**
			* @brief Erase annotationset for leads pointed to by iterator it
			*
			* @param it Iterator pointing to leads to remove annotations from
			*/
			void erase(iterator it) {
				_pm.erase(it);
			}

			/**
			* @brief Insert annotationsets from other leads into this pointmap
			*
			* @tparam INPUT_ITER Iterator type
			* @param a1 Start iterator
			* @param a2 Stop iterator
			*
			* @return Iterator to last inserted element
			*/
			template<typename INPUT_ITER>
			iterator insert(INPUT_ITER a1, INPUT_ITER a2) {
				iterator ret;

				while(a1 != a2) {
					ret = _pm.insert(*a1);

					++a1;
				}

				return ret;
			}

			/**
			* @brief Return number of annotations in a pointmap
			*
			* @return Number of annotations in pointmap
			*/
			std::size_t nanns() const {
				std::size_t n = 0;

				for(auto &r : _pm) {
					n += r.second.size();
				}

				return n;
			}

		private:
			std::map<leadnumber, annotationset> _pm;
	};

	/**
	* @brief Predicate function for if annotation in
	*/
	class annotationtype_in : public std::unary_function<annotation, bool> {
		public:
			annotationtype_in(const std::vector<annotation> &list) : _alist(list) {
			}

		public:
			result_type operator()(argument_type in) {
				for(std::size_t i = 0; i < _alist.size(); ++i) {
					if(_alist[i] == in) return true;
				}

				return false;
			}

		private:
			std::vector<annotation> _alist;
	};

	/**
	 * @brief Get annotations by type
	 *
	 * @param pm Annotation set
	 * @param typ Annotation type
	 * @param[out] out Vector of annotations
	 */
	void get_annotations(const annotationset &pm, const ecglib::annotation_type &typ, std::vector<annotation> &out);

	/**
	 * @brief Get annotations across leads by type
	 *
	 * @param pm Pointmap
	 * @param typ Annotation type
	 * @param[out] out Vector of annotations
	 */
	void get_annotations(const pointmap &pm, const ecglib::annotation_type &typ, std::vector<annotation> &out);

	/**
	* @brief Get annotations across leads by type / type
	*
	* @param pm Pointmap
	* @param l Leadnumber
	* @param typ Type of annotation
	* @param[out] out Output
	*/
	void get_annotations(const pointmap &pm, const leadnumber l, const ecglib::annotation_type &typ, std::vector<annotation> &out);

	/**
	 * @brief Shift all annotations to a new start, used for chopping an ECG
	 *
	 * @param ann Annotationset
	 * @param newstart New start
	 *
	 * @return New annotationset
	 */
	annotationset rebase(const annotationset &ann, const int newstart);

	/**
	 * @brief Shifts all annotations to a new start, used for chopping an ECG
	 *
	 * @param pm Pointmap
	 * @param newstart New start
	 */
	void rebase(pointmap &pm, const int newstart);

	/**
	 * @brief Used to generate global annotations, i.e. across leads using mean.
	 *
     * This globalizer works by taking the number (passed on constructor) from a direction left(-1)/right(1)
	 * and computes the mean. Example:
	 *  Following Ponset locations are passed through: 1, 2, 3, 4 and 10
	 * A globalizer with direction -1 is used (as we are looking at the onset) with nnumber=3, this will result
	 * in the mean being: mean([1,2,3]) had direction 1 been used it would have been mean([3,4,10]).
	 */
	class globalizer {
		public:
			/**
			 * @brief Globalize constructor
			 *
			 * @param direction Direction (-1: from the left, 1: from the right)
			 * @param nnum Number of annotations to use
			 */
			globalizer(const int direction, const int nnum) : _direction(static_cast<direc>(direction)), _nnum(nnum) {
			}

		public:
			/**
			* @brief compute mean time in ms
			*
			* @param vals values to compute
			*/
			timems operator()(const uvec &vals) const {
				uvec svals = sort(vals);
				int cnt = svals.n_elem;

				if(cnt <= _nnum) {
					return static_cast<timems>(round(as_scalar(mean(svals.subvec(0,cnt-1)))));
				}

				if(_direction == -1) {
					return static_cast<timems>(round(as_scalar(mean(svals.subvec(0,_nnum-1)))));
				} else {
					return static_cast<timems>(round(as_scalar(mean(svals.subvec(cnt-_nnum,cnt-1)))));
				}
			}
		private:
			enum direc {left = -1, right = 1};		/**< @brief enumeration for direction */
			direc _direction;				/**< @brief -1: from the left, 1: from the right */
			const int _nnum;				/**< @brief number of elements */
	};


	/**
	 * @brief Median globalizer, similar to the mean globalizer (globalizer)
	 *
	 * But, instead of computing the mean a percentile is used, which when set to 0.5, corresponds to the median.
	 * So, if direction left(-1) and perc 0.25 then the 25th percentile from the left would be the onset. If perc
	 * was set to either 0 (for left) or 1 (for right) it corresponds to earliest and latest.
	 */
	class mglobalizer {
		public:
			/**
			 * @brief Constructor for median globalizer
			 *
			 * @param direction Left(-1) or right(1)
			 * @param perc Percentile to use
			 */
			mglobalizer(const int direction, const double perc) : _direction(static_cast<direc>(direction)), _perc(perc) {
			}

		public:
			/**
			* @brief compute the median time in ms
			*
			* @parm vals values to compute
			*/
			timems operator()(const uvec &vals) const {
				uvec svals = sort(vals);

				double perc = _perc;
				if(_direction == -1) perc = 1.0 - perc;

				int idx = static_cast<int>(round(perc * svals.n_elem));
				idx = idx == static_cast<int>(svals.n_elem) ? idx-1 : idx;
				idx = idx < 0 ? 0 : idx;

				return svals[idx];
			}

		private:
			enum direc {left = -1, right = 1};	/**< @brief enumeration for direction */
			const int _direction;			/**< @brief -1: from the left, 1: from the right */
			const double _perc;			/**< @brief percentile */
	};

	/*!
	 *@}
	 */
}

#endif
