/**
 * @file core/ecglib/ecgdata.hpp
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
 * The classes in this file are used for describing characteristics about ecg records
 */

#ifndef ECGLIB_CORE_ECGDATA_LJ_2015_12_09
#define ECGLIB_CORE_ECGDATA_LJ_2015_12_09 1

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/assert.hpp>
#include <boost/bimap.hpp>

#include <armadillo>

#include <ecglib/util/util.hpp>
#include <ecglib/annotation.hpp>
#include <ecglib/ecglib.hpp>


#include <algorithm>
#include <vector>
#include <string>
#include <map>

namespace ecglib {
	/*! \addtogroup core
	 * Core ECGlib classes and functions
	 * @{
	 */

	using namespace arma;

	/**
	* @brief ECGdata class, contains a matrix (nsamples x nleads), vector of lead names and annotations by lead
	*
	* ECGdata consists of:
	*  - Armadillo of type T to define samples per samples in units of uV (assumed for now) resolution defined in fs
	*    Samples x Leads
	*  - Properties Auxilary properties, e.g. configuration used for earlier processing step, heart rate for a median beat, start of recording etc
	*  - Annotations, Annotations for ECG
	*  - Lead map, defines lead types for leads with known type
	*
	* Annotations are stored in a pointmap which links between columns (leads) and annotations.
	*/
	template<class T>
	class Ecgdata {
		public:
		// NOTE: Loops all samples for a leads - lj
		typedef mat::col_iterator leaditerator;				/**< @brief iterator for ecgleads */
		typedef mat::const_col_iterator const_leaditerator;		/**< @brief constant iterator for ecgleads */

		// NOTE: Loops all samples for all leads - lj
		typedef mat::row_iterator sampleiterator;			/**< @brief iterator for samples across leads */
		typedef mat::const_row_iterator const_sampleiterator;		/**< @brief constant iterator for samples across leads */

		typedef ecglib::annotationset annotationset;			/**< @brief annotation set */
		typedef ecglib::pointmap pointmap;				/**< @brief pointmap */
		typedef pointmap::iterator pointmapiterator;			/**< @brief iterator for pointmap */
		typedef pointmap::const_iterator const_pointmapiterator;	/**< @brief constant iterator for pointmap */

		typedef annotationset::iterator pointiterator;			/**< @brief iterator for annotation set */
		typedef annotationset::const_iterator const_pointiterator;	/**< @brief constant iterator for annotation set */

		typedef boost::bimap<leadnumber, ecglead> leadmap;		/**< @brief bimap of leadnumber and ecglead is a leadmap */

		// Constructors
		public:

			/**
			* @brief Creates empty ecgdata
			*/
			Ecgdata() : _nsamples(0), _res(1), _nleads(0) {}

			/**
			* @brief Creates an ecgdata object of length nsamples, but number of leads based on ecgheader information
			*
			* @param nsamples Number of samples (could be equal of bigger than nsample  from ecgheader
			* @param ecgheader
			*/
            		Ecgdata(unsigned int nsamples, ecgheader eh) : _data(nsamples,eh.nleads), _nsamples(nsamples),_res(1) {

			}

			/**
			* @brief Creates an ecgdata object of length nsamples, number of columns/leads corresponding to leadnames size
			*
			* @param nsamples Number of samples
			* @param leadnames Lead names of the leads
			*/
			Ecgdata(unsigned int nsamples, const std::vector<ecglead> &leadnames) : _data(nsamples,leadnames.size()), _nsamples(nsamples), _res(1), _nleads(leadnames.size()) { 
				typedef leadmap::value_type lval;

				for(std::size_t i = 0; i < leadnames.size(); ++i) {
					_leadmap.insert(lval(i, leadnames[i]));
				}
			}

			/**
			* @brief Prepare empty ECG data class
			*
			* @param nsamples Number of samples
			* @param nleads Number of leads
			*/
			Ecgdata(unsigned int nsamples, unsigned int nleads) : _data(nsamples, nleads), _nsamples(nsamples), _res(1), _nleads(nleads) {
			}

			/**
			* @brief Create ECG data from armadillo matrix with a lead map
			*
			* @param indata Input data
			* @param lm Leadmap
			*/
			Ecgdata(const mat &indata, const leadmap &lm) : _data(indata), _nsamples(indata.n_rows), _leadmap(lm), _res(1), _nleads(indata.n_cols) {
			}

			/**
			* @brief Create ECG data from armadillo matrix with a fs/res
			*
			* @param indata Input data
			* @param fs Sampling frequency
			* @param res Resolution
			*/
			Ecgdata(const mat &indata, const double fs, const double res=1) : _data(indata), _nsamples(indata.n_rows), _fs(fs), _res(res), _nleads(indata.n_cols) {
			}

			/**
			* @brief Creates an ecgdata class from a matrix of data and a set of lead names (first column is first element in leadnames)
			*
			* @param indata Matrix of data
			* @param leadnames Lead names of the data
			*/
			Ecgdata(const mat &indata, const std::vector<ecglead> &leadnames) : _data(indata), _nsamples(indata.n_rows), _nleads(leadnames.size()) {
				if(indata.n_cols != _nleads) {
					std::cerr << "ecglib::constructor::Length of leadnames does not match column count";
					throw ecglib::ecglib_exception("ecglib::constructor::Length of leadnames does not match column count");
				}

				typedef leadmap::value_type lval;

				for(std::size_t i = 0; i < leadnames.size(); ++i) {
					_leadmap.insert(lval(i, leadnames[i]));
				}
			}

		// Public methods
		public:

			/**
			* @brief Gets the column number of the lead. Returns false if it cannot find that lead
			*
			* @param lead Lead to get the column number for
			*
			* @return bool true if has lead number
			*/
			bool hasleadnum(const ecglead lead) const {
				leadmap::right_const_iterator rit = _leadmap.right.find(lead.index);

				if(rit == _leadmap.right.end()) {
					return false;
				}

				return true;
			}


			/**
			* @brief Gets the column number of the lead. Throws an exception if the lead is not in the data
			*
			* @param lead Lead to get the column number for
			*
			* @return Column number
			*/
			int leadnum(const ecglead lead) const {
				leadmap::right_const_iterator rit = _leadmap.right.find(lead.index);

				if(rit == _leadmap.right.end()) {
					std::cerr << "No such lead";
					throw ecglib::ecglib_exception("No such lead");
				}

				return rit->second;
			}

			/**
			* @brief Get leadname for a column number
			*
			* @param leadnum Column number
			*
			* @return lead name
			*/
			ecglead leadname(int leadnum) const {
				leadmap::left_const_iterator lit = _leadmap.left.find(leadnum);

				if(lit == _leadmap.left.end()) {
					std::cerr << "Lead has no lead name";
					throw ecglib::ecglib_exception("Lead has no lead name");
				}

				return lit->second;
			}

			/**
			* @brief Adds a lead to the container. Causes container to resize internal storage (not efficient)
			*
			* @tparam ITER An iterator, must be usable by std::copy
			* @param lead Lead name for the lead to add
			* @param start Start of what to copy
			* @param stop Stop of what to copy
			*/
			template<class ITER>
			void add_lead(const ecglead lead, ITER start, ITER stop) {
				std::size_t length = stop-start;
				int newlead = _nleads;

				if(_nsamples == 0) {
					_nsamples = length;
				} else if(_nsamples != length) {
					std::cerr << "ecglib::add_lead: _nsamples != length";
					throw ecglib::ecglib_exception("ecglib::add_lead: _nsamples != length");
				} else {
					leadmap::right_iterator rit = _leadmap.right.find(lead.index);

					if(rit != _leadmap.right.end()) {
						std::cerr << "ecglib::add_lead:lead is not new";
						throw ecglib::ecglib_exception("ecglib::add_lead:lead is not new");
					}
				}

				_data.reshape(_nsamples,_nleads+1);

				std::copy(start,stop,_data.begin_col(_nleads));

				leadmap::value_type lval(newlead, lead.index);
				_leadmap.insert(lval);
				++_nleads;
			}

		// Setters / Getters
		public:
			/**
			* @brief Get sampling frequency
			*
			* @return Sampling frequency
			*/
			double fs() const {
				return _fs;
			}

			/**
			* @brief Set sampling frequency
			*
			* @param fs2 Sampling frequency
			*/
			void fs(const double fs2) {
				_fs = fs2;
			}

			/**
			* @brief Get resolution
			*
			* @return Resolution how many units per uV, i.e. if data is in uV it is 1
			*/
			double resolution() const {
				return _res;
			}

			/**
			* @brief Set resolution
			*
			* @param res2 resolution
			*/
			void resolution(const double res2) {
				_res = res2;
			}

			/**
			* @brief Get number of samples
			*
			* @return Nsamples
			*/
			std::size_t nsamples() const {
				return _nsamples;
			}

			/**
			* @brief Get number of leads
			*
			* @return Nleads
			*/
			std::size_t nleads() const {
				return _nleads;
			}

			/**
			* @brief Get annotations for data
			*
			* @return annotations
			*/
			pointmap pointsmap() const {
				return _points;
			}

			/**
			* @brief Set annotaitons for data
			*
			* @param pts annotations
			*/
			void pointsmap(const pointmap &pts) {
				_points = pts;
			}

			/**
			* @brief Get vector of leads stored
			*
			* @return Vector leads
			*/
			leadmap leadnames() const {
				return _leadmap;
			}

			/**
			* @brief Get vector of leads stored
			*
			* @return Vector leads
			*/
			void leadnames(const leadmap &lm) {
				_leadmap = lm;
			}

			/**
			* @brief Set all properties
			*
			* @param props Props
			*/
			void setproperties(const propertymap &props) {
				_props = props;
			}

			/**
			* @brief Insert property map into internal map
			*
			* @param props
			*/
			void insertproperties(const propertymap &props) {
				for(auto &r : props) {
					_props.insert(std::pair<std::string, property>(r.first,r.second));
				}
			}

			/**
			* @brief Get all properties
			*
			* @return Property map
			*/
			propertymap getproperties() const {
				return _props;
			}

			/**
			* @brief Set a property
			*
			* @param prop Property name
			* @param value Value
			*/
			void setproperty(const std::string &nam, const property &prop) {
				auto propit = _props.find(nam);
				if(propit != _props.end()) {
					_props.erase(propit);
				}

				_props.insert(std::pair<std::string, property>(nam,prop));
			}

			/**
			* @brief Determine if a property is set
			*
			* @param prop Property name
			*
			* @return True/false
			*/
			bool hasproperty(const std::string &prop) const {
				auto it = _props.find(prop);

				return !(it == _props.end());
			}

			/**
			* @brief Number of properties for an ECG
			*
			* @return Number of properties
			*/
			std::size_t nproperties() const {
				return _props.size();
			}

			/**
			* @brief Get a property
			*
			* @param prop Property
			*
			* @return Value of property
			*/
			property getproperty(const std::string &prop) const {
				if(!hasproperty(prop)) {
					std::string line = std::string("No such property: ") + prop;
					std::cerr << line;
					throw std::logic_error(line);
				}

				auto it = _props.find(prop);

				return it->second;
			}

			/**
			* @brief Beginning of properties
			*
			* @return Const iterator to start
			*/
			propertymap::const_iterator begin_properties() const {
				return _props.begin();
			}

			/**
			* @brief End of properties
			*
			* @return Const iterator to end
			*/
			propertymap::const_iterator end_properties() const {
				return _props.end();
			}

		// Container methods
		public:
			/**
			* @brief Get sample reference by lead,sample
			*
			* @param lead lead
			* @param sample Sample number
			*
			* @return sample reference
			*/
			double& operator()(const ecglead &lead, const int sample) {
				int leadn = leadnum(lead);		

				return _data(sample,leadnum(leadn));
			}

			/**
			* @brief Get sample value by lead,sample
			*
			* @param lead lead
			* @param sample Sample number
			*
			* @return sample value
			*/
			double operator()(const ecglead &lead, const int sample) const {
				int leadn = leadnum(lead);		

				return _data(sample,leadnum(lead));
			}

			/**
			* @brief Get subview of ECG for a lead
			*
			* @param rowspan Row span for samples
			* @param colnum Column number
			*
			* @return col vec of subview of lead of ECG
			*/
			const subview_col<T> operator()(const span& rowspan, const uword colnum) const {
				return _data(rowspan, colnum);
			}

			/**
			* @brief Get sample reference by lead,sample
			*
			* @param lead lead number
			* @param sample sample
			*
			* @return sample reference
			*/
			double& operator()(const int lead, const int sample) {
				return _data(sample,lead);
			}

			/**
			* @brief Get sample value by lead,sample
			*
			* @param lead lead number
			* @param sample sample
			*
			* @return sample value
			*/
			double operator()(const int lead, const int sample) const {
				return _data(sample,lead);
			}

			/**
			* @brief Get lead by number
			*
			* @param lead lead number
			*
			* @return lead (subview)
			*/
			subview_col<double> lead(const ecglead lead) {
				int leadn = leadnum(lead);		

				return _data.col(leadn);
			}

			/**
			* @brief Const-correct get lead by number
			*
			* @param lead lead number
			*
			* @return lead (subview) constant
			*/
			const subview_col<double> lead(const ecglead lead) const {
				int leadn = leadnum(lead);		

				return _data.col(leadn);
			}

			/**
			* @brief Get lead by number
			*
			* @param lead lead number
			*
			* @return lead (subview)
			*/
			subview_col<double> lead(const int lead) {
				return _data.col(lead);
			}

			/**
			* @brief Const-correct get lead by number
			*
			* @param lead lead number
			*
			* @return lead (subview) constant
			*/
			const subview_col<double> lead(const int lead) const {
				return _data.col(lead);
			}

			/**
			* @brief Get part of a lead
			*
			* @param lead lead number
			* @param start start sample
			* @param stop stop sample
			*
			* @return part of lead 
			*/
			subview_col<double> lead(const int lead, const int start, const int stop) {
				return _data(span(start,stop),lead);
			}

			/**
			* @brief Get part of a lead, const-correct
			*
			* @param lead lead number
			* @param start start sample
			* @param stop stop sample
			*
			* @return part of lead constant
			*/
			const subview_col<double> lead(const int lead, const int start, const int stop) const {
				return _data(span(start,stop),lead);
			}

			/**
			* @brief Get part of a lead
			*
			* @param lead lead number
			* @param start start sample
			* @param stop stop sample
			*
			* @return part of lead 
			*/
			subview_col<double> lead(const ecglead lead, const int start, const int stop) {
				int leadn = leadnum(lead);		

				return _data(span(start,stop),leadn);
			}

			/**
			* @brief Get part of a lead, const-correct
			*
			* @param lead lead number
			* @param start start sample
			* @param stop stop sample
			*
			* @return part of lead constant
			*/
			const subview_col<double> lead(const ecglead lead, const int start, const int stop) const {
				int leadn = leadnum(lead);

				return _data(span(start,stop),leadn);
			}

			/**
			* @brief get data or matrix of all data
			*
			* @return data
			*/
			subview<double> data() {
				return _data(span::all,span::all);
			}

			/**
			* @brief get data or matrix of all data, const-correct
			*
			* @return data const
			*/
			subview<double> data() const {
				return _data(span::all,span::all);
			}

			/**
			* @brief Begin lead iterator
			*
			* @param lead lead
			*
			* @return lead iterator
			*/
			leaditerator begin_lead(const ecglead &lead) {
				return _data.begin_col(leadnum(lead));
			}

			/**
			* @brief End lead iterator
			*
			* @param lead lead
			*
			* @return lead iterator
			*/
			leaditerator end_lead(const ecglead &lead) {
				return _data.end_col(leadnum(lead));
			}

			/**
			* @brief Begin lead iterator
			*
			* @param lead lead
			*
			* @return lead iterator const
			*/
			const_leaditerator begin_lead(const ecglead &lead) const {
				return _data.begin_col(leadnum(lead));
			}

			/**
			* @brief End lead iterator
			*
			* @param lead lead
			*
			* @return lead iterator const
			*/
			const_leaditerator end_lead(const ecglead &lead) const {
				return _data.end_col(leadnum(lead));
			}

			/**
			* @brief Begin lead iterator
			*
			* @param lead lead number
			*
			* @return lead iterator
			*/
			leaditerator begin_lead(int num) {
				return _data.begin_col(num);
			}

			/**
			* @brief Begin lead iterator
			*
			* @param lead lead number
			*
			* @return lead iterator
			*/
			leaditerator end_lead(int num) {
				return _data.end_col(num);
			}

			/**
			* @brief Begin lead iterator
			*
			* @param lead lead number
			*
			* @return lead iterator const
			*/
			const_leaditerator begin_lead(int num) const {
				return _data.begin_col(num);
			}

			/**
			* @brief End lead iterator
			*
			* @param lead lead number
			*
			* @return lead iterator const
			*/
			const_leaditerator end_lead(int num) const {
				return _data.end_col(num);
			}

			/**
			* @brief Begin sample (across leads)
			*
			* @param sample Sample number
			*
			* @return sample iterator
			*/
			sampleiterator begin_sample(int sample) {
				return _data.begin_row(sample);
			}

			/**
			* @brief End sample (across leads)
			*
			* @param sample Sample number
			*
			* @return sample iterator
			*/
			sampleiterator end_sample(int sample) {
				return _data.end_row(sample);
			}

			/**
			* @brief Begin sample (across leads)
			*
			* @param sample Sample number
			*
			* @return sample iterator const
			*/
			const_sampleiterator begin_sample(int sample) const {
				return _data.begin_row(sample);
			}

			/**
			* @brief End sample (across leads)
			*
			* @param sample Sample number
			*
			* @return sample iterator const
			*/
			const_sampleiterator end_sample(int sample) const {
				return _data.end_row(sample);
			}

			/**
			* @brief Get begin pointmap const iterator
			*
			* @param lead lead
			*
			* @return iterator to end
			*/
			pointmapiterator pointmap_begin() {
				return _points.begin();
			}

			/**
			* @brief Get end pointmap iterator
			*
			* @param lead lead
			*
			* @return iterator to end
			*/
			pointmapiterator pointmap_end() {
				return _points.end();
			}

			/**
			* @brief Get begin pointmap const iterator
			*
			* @param lead lead
			*
			* @return iterator to begin
			*/
			const_pointmapiterator pointmap_begin() const {
				return _points.begin();
			}

			/**
			* @brief Get end pointmap const iterator
			*
			* @param lead lead
			*
			* @return iterator to end
			*/
			const_pointmapiterator pointmap_end() const {
				return _points.end();
			}

			/**
			* @brief Get begin-point iterator by lead
			*
			* @param lead lead
			*
			* @return iterator to begin
			*/
			pointiterator points_begin(const ecglead &lead) {
				int leadn = leadnum(lead);

				pointmapiterator iter = _points.find(leadn);

				if(iter == _points.end()){
					std::cerr << "Could not find annotations for lead";
					throw std::runtime_error("Could not find annotations for lead");
				}

				return iter->second.begin();
			}

			/**
			* @brief Get begin-point iterator by lead
			*
			* @param lead lead
			*
			* @return iterator to begin
			*/
			pointiterator points_begin(const leadnumber leadn) {
				pointmapiterator iter = _points.find(leadn);

				if(iter == _points.end()){
					std::cerr << "Could not find annotations for lead";
					throw std::runtime_error("Could not find annotations for lead");
				}

				return iter->second.begin();
			}

			/**
			* @brief Get end-point iterator by lead
			*
			* @param lead lead
			*
			* @return iterator to begin
			*/
			pointiterator points_end(const ecglead &lead) {
				int leadn = leadnum(lead);

				pointmapiterator iter = _points.find(leadn);

				if(iter == _points.end()){
					std::cerr << "Could not find annotations for lead";
					throw std::runtime_error("Could not find annotations for lead");
				}

				return iter->second.end();
			}

			/**
			* @brief Get end-point iterator by lead
			*
			* @param lead lead
			*
			* @return iterator to begin
			*/
			pointiterator points_end(const leadnumber leadn) {
				pointmapiterator iter = _points.find(leadn);

				if(iter == _points.end()){
					std::cerr << "Could not find annotations for lead";
					throw std::runtime_error("Could not find annotations for lead");
				}

				return iter->second.end();
			}

			/**
			* @brief Get begin-point iterator by lead
			*
			* @param lead lead
			*
			* @return const iterator to begin
			*/
			const_pointiterator points_begin(const leadnumber leadn) const {
				const_pointmapiterator iter = _points.find(leadn);

				if(iter == _points.end()){
					std::cerr << "Could not find annotations for lead";
					throw std::runtime_error("Could not find annotations for lead");
				}

				return iter->second.begin();
			}

			/**
			* @brief Get end-point iterator by lead
			*
			* @param lead lead
			*
			* @return const iterator to end
			*/
			const_pointiterator points_end(const leadnumber leadn) const {
				const_pointmapiterator iter = _points.find(leadn);

				if(iter == _points.end()){
					std::cerr << "Could not find annotations for lead";
					throw std::runtime_error("Could not find annotations for lead");
				}

				return iter->second.end();
			}

			/**
			* @brief Get begin-point iterator by lead
			*
			* @param lead lead
			*
			* @return const iterator to begin
			*/
			const_pointiterator points_begin(const ecglead &lead) const {
				int leadn = leadnum(lead);

				const_pointmapiterator iter = _points.find(leadn);

				if(iter == _points.end()){
					std::cerr << "Could not find annotations for lead";
					throw std::runtime_error("Could not find annotations for lead");
				}
				return iter->second.begin();
			}

			/**
			* @brief Get end-point iterator by lead
			*
			* @param lead lead
			*
			* @return const iterator to end
			*/
			const_pointiterator points_end(const ecglead &lead) const {
				int leadn = leadnum(lead);

				const_pointmapiterator iter = _points.find(leadn);

				if(iter == _points.end()){
					std::cerr << "Could not find annotations for lead";
					throw std::runtime_error("Could not find annotations for lead");
				}

				return iter->second.end();
			}

			/**
			* @brief Extracts a sub-part of an ecgdata class
			*
			* @param start First sample
			* @param stop Last sample
			*
			* @return ECGData of the sub part (both data/annotations)
			*/
			Ecgdata<T> subpart(int starttime, int stoptime) const {
				if (stoptime <= 0) {
					unsigned int length = static_cast<unsigned int>(round( (_nsamples/_fs)*1000. ));

					stoptime += length-1;
				}
				unsigned int start = static_cast<unsigned int>(round((starttime/1000.) * _fs)); 
				unsigned int stop = static_cast<unsigned int>(round((stoptime/1000.) * _fs)); 

				if(start < 0 || stop > _nsamples || start > stop) {
					std::string line = std::string("Start_or_stop_incorrect");
					std::cerr << line;
					throw ecglib::ecglib_exception(line);
				}

				Ecgdata<T> e(_data(span(start,stop),span::all), _leadmap);

				e.fs( _fs);
				e.resolution(_res);

				// Copy annotations
				pointmap pm;

				const_pointmapiterator pmend = _points.end();
				for(const_pointmapiterator pmi = _points.begin(); pmi != pmend; ++pmi) {
					const_pointiterator pistart = pmi->second.lower_bound(starttime);
					const_pointiterator pistop = pmi->second.upper_bound(stoptime);

					annotationset as;

					while(pistart != pistop) {
						as[pistart->first-starttime] = pistart->second-static_cast<unsigned int>(starttime);

						++pistart;
					}

					pm[pmi->first] = as;
				}

				e.pointsmap(pm);

				return e; 
			}

			/**
			 * @brief Extracts a sub-part of an ecgdata class
			 *
			 * @param leads ECG leads
			 *
			 * @return ECGdata
			 */
			Ecgdata<T> subpart(const std::vector<ecglib::ecglead> &leads) const {
				mat ecg = zeros<mat>(_nsamples+1, leads.size());

				for(std::size_t i = 0; i < leads.size(); ++i) {
					std::copy(begin_lead(leads[i]),end_lead(leads[i]), ecg.begin_col(i));
				}

				Ecgdata<T> e(ecg, leads);
				e.fs( _fs);
				e.resolution(_res);
				Ecgdata<T>::pointmap pm;
				typedef Ecgdata<T>::const_pointmapiterator pmiter;
				pmiter pmend = _points.end();

				typedef std::vector<ecglead>::const_iterator liter;
				liter lend = leads.end();

				for(pmiter pmi = _points.begin(); pmi != pmend; ++pmi) {
					if(pmi->first == ecglib::GLOBAL_LEAD) {
						pm[pmi->first] = pmi->second;
					} else {
						for(liter li = leads.begin(); li != lend; ++li) {
							if(*li == pmi->first) {
								pm[pmi->first] = pmi->second;
							}
						}
					}
				}

				e.pointsmap(pm);

				return e;
			}

		// Helpers
		private:

		// Attributes
		protected:
			/**
			* @brief Data, stored as samples x leadsc
			*/
			Mat<T> _data;

			/**
			* @brief Number of samples
			*/
			std::size_t _nsamples;

			/**
			* @brief Lead map, maps column numbers with known lead names
			*/
			leadmap _leadmap;

			/**
			* @brief Sampling frequency
			*/
			double _fs;

			/**
			* @brief Resolution..
			*/
			double _res;

			/**
			* @brief Number of leads
			*/
			std::size_t _nleads;

			/**
			* @brief Pointmap
			*/
			pointmap _points;

			/**
			* @brief Property map
			*/
			propertymap _props;
	};

	/**
	* @brief ECGdata of type double is the currently accepted input for all functions.
	*/
	typedef Ecgdata<double> ecgdata;

	/*!
	 *@}
	 */
}


#endif
