#include "../include/config_parser.h"

ConfigParser::ConfigParser() {
    //do nothing
}
ConfigParser::ConfigParser(const std::string &file){
    // Extract the file extension from the filename
    std::string ext = boost::filesystem::extension(file);
    if (ext == ".json") {
        // parse json and write into 'tree'
        boost::property_tree::read_json(file, tree);
    } else if (ext == ".info"){
        // parse info and write into 'tree'
        boost::property_tree::read_info(file, tree);
    } else {
        //BOOST_THROW_EXCEPTION(std::logic_error("Unsupported file extension"));
        std::string errMsg = ext + ": unsupported file extension (supported: .json, .info)";
        BOOST_THROW_EXCEPTION(std::invalid_argument(errMsg));
    }
}

ConfigParser::ConfigParser(const boost::property_tree::ptree &subtree) : tree { subtree } {}

ConfigParser ConfigParser::getObj(const std::string &key){
    return ConfigParser(tree.get_child(key));
}

std::list<ConfigParser> ConfigParser::getObjList(const std::string &key) {
    std::list <ConfigParser> lst_;
    BOOST_FOREACH(const boost::property_tree::ptree::value_type &val,tree.get_child(key))
    {
        if(val.second.empty()){
            std::cerr << "List does not contain objects. Please use 'getList<T>(const std::string &key)'instead."
                                  << " - Returning empty list." << std::endl;
        } else {
            // storing instances of ConfigParsers each containing a subtree instead
            lst_.push_back(ConfigParser(val.second.get_child("")));
        }
    }
    return lst_;
}