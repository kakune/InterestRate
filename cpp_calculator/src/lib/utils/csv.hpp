#ifndef UTILS_CSV_HPP
#define UTILS_CSV_HPP

#include <cmath>
#include <map>
#include <string>
#include <vector>

#include "utils/parameters.hpp"

namespace Utils::CSV
{

std::map<std::string, std::vector<double>> readFile(
    const std::string& inPathCSV );

void prepareZCBColumn( std::map<std::string, std::vector<double>>& inMapMarket,
                       const Parameters& inParams );

}  // namespace Utils::CSV

#endif