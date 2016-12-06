/**
 * @file core/ecglib/util/config.hpp
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
 * The classes in this file are used for converting to/from command line options and properties for ecg records and analysis methods
 */

#ifndef ECGLIB_CORE_UTIL_CONFIG_LJ_2015_12_09
#define ECGLIB_CORE_UTIL_CONFIG_LJ_2015_12_09 1

#include <boost/any.hpp>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include <map>
#include <ecglib/ecglib.hpp>

namespace ecglib { 
	/*! \addtogroup core
	 * Core ECGlib classes and functions
	 * @{
	 */

	/**
	* @brief Converts from boost::program_options to propertymap
	*/
	class from_programoptions {
		public:
			/**
			* @brief Constructor for class
			*
			* @param prefix Prefix of the program options name
			* @param vm Name of program options to convert
			*/
			from_programoptions(const std::string &prefix, const po::variables_map &vm) : _prefix(prefix), _vm(vm) {
			}
		public:
			/**
                        * @brief get option from map
                        *
                        * @param nam name of option to get
                        * @param typ type to cast it to (String, Double, Int, Uint)
                        */
			boost::any get(const std::string &nam, const Type typ) const {
				std::string fullnam = _prefix + nam;

				if(typ == Type::String) {
					return _vm[fullnam].as<std::string>();
				} else if(typ == Type::Double) {
					return _vm[fullnam].as<double>();
				} else if(typ == Type::Int) {
					return _vm[fullnam].as<int>();
				} else if(typ == Type::Uint) {
					return _vm[fullnam].as<unsigned int>();
				} else {
					std::cerr << "Unsupported type";
					throw std::logic_error("Unsupported type");
				}
			}

		private:
			std::string _prefix;	/**< @brief prefix of options/properties */
			po::variables_map _vm;	/**< @brief program options to map */
	};

	/**
	* @brief Converts from propertymap to boost::program_options
	*/
	class to_programoptions {
		public:
			/**
                        * @brief Constructor for conversion class
                        *
                        * @param prefix Prefix of properties/ program options name
                        * @param name Name of properties/ program options
                        */
			to_programoptions(const std::string &prefix, const std::string &name) : _opts(name), _prefix(prefix) {

			}
		public:
                        /**
                        * @brief set option from property
                        *
                        * @param nam Name of option to set
                        * @param prop The property to set
                        */
			void set(const std::string &nam, const property &prop) {
				std::string fullnam = _prefix + nam;

				if(prop.type == Type::String) {
					boost::shared_ptr<po::option_description> desc(new po::option_description(fullnam.c_str(), po::value<std::string>()->default_value(boost::any_cast<std::string>(prop.value)),prop.desc.c_str()));
					_opts.add(desc);
				} else if(prop.type == Type::Double) {
					boost::shared_ptr<po::option_description> desc(new po::option_description(fullnam.c_str(), po::value<double>()->default_value(boost::any_cast<double>(prop.value)),prop.desc.c_str()));
					_opts.add(desc);
				} else if(prop.type == Type::Int) {
					boost::shared_ptr<po::option_description> desc(new po::option_description(fullnam.c_str(), po::value<int>()->default_value(boost::any_cast<int>(prop.value)),prop.desc.c_str()));
					_opts.add(desc);
				} else if(prop.type == Type::Uint) {
					boost::shared_ptr<po::option_description> desc(new po::option_description(fullnam.c_str(), po::value<unsigned int>()->default_value(boost::any_cast<unsigned int>(prop.value)),prop.desc.c_str()));
					_opts.add(desc);
				} else {
					std::cerr << "Unsupported type";
					throw std::logic_error("Unsupported type");
				}
			}

			po::options_description opts() const {
				return _opts;
			}

		private:
			po::options_description _opts;	/**< @brief the properties as program options */
			std::string _prefix;		/**< @brief the prefix of the name of the properties */
	};

	/**
	* @brief ECGlib config class
	*/
	class config {
		public:
			typedef property propertytype;

		public:
                        /**
                        * @brief Empty constructor
                        */
			config() {
			}

		public:
			/**
			* @brief Load configuration (e.g. program_options)
			*
			* @tparam T Configuration loader
			* @param t Configuration loader
			*/
			template<class T>
			void load(const T &t) {
				for(auto &r : _props) {
					r.second.value = t.get(r.first, r.second.type);
				}
			}

			/**
			* @brief Save configuration (e.g. program_options)
			*
			* @tparam T Configuration saver
			* @param t Configurationsaver
			*/
			template<class T>
			void save(T &t) const {
				for(auto const &r : _props) {
					t.set(r.first, r.second);
				}
			}

			/**
			* @brief Set configuration
			*
			* @param nam Name of property
			* @param typ Type of property
			* @param val Value
			*/
			void set(const std::string &nam, const Type typ, const boost::any &val) {
				auto it = _props.find(nam);

				if(it == _props.end()) {
					std::cerr << "No such key";
					throw std::logic_error("No such key");
				}

				if(it->second.type != typ) {
					std::cerr << "Wrong type";
					throw std::logic_error("Wrong type");
				}

				it->second.value = val;
			}

			/**
			* @brief Get configuration
			*
			* @tparam T Convert to type T
			* @param nam Name of configuration
			*
			* @return Value of configuration cast to T
			*/
			template<class T>
			T get(const std::string &nam) const {
				auto it = _props.find(nam);

				if(it == _props.end()) {
					std::string line = std::string("No such key: ") + nam;
					std::cerr << line;
					throw std::logic_error(line);
				}

				return boost::any_cast<T>(it->second.value);
			}

			/**
			* @brief Build propertymap based on config and add prefix
			*
			* @param prefix Prefix
			*
			* @return Propertymap
			*/
			propertymap pm(const std::string &prefix) const {
				propertymap pm;

				for(auto &r : _props) {
					std::string onam = prefix + r.first;
					pm.insert(std::pair<std::string, property>(onam, r.second));
				}

				return pm;
			}

			/**
			* @brief Get a copy of propertymap
			*
			* @return Propertymap
			*/
			propertymap pm() const {
				return _props;
			}

		protected:
			/**
			* @brief Protected function inserting properties
			*
			* @param nam Name string
			* @param prop Rvalue of property
			*/
			void add(const char *nam, const property &&prop) {
				_props.insert(std::pair<std::string, property>(nam, prop));
			}

		protected:
			/**
			* @brief Internal map of string, property
			*/
			std::map<std::string, property> _props;
	};

	/*!
	 *@}
	 */
}

#endif
