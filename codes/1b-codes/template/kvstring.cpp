#include <cassert>
#include <sstream>
#include <iomanip>
#include <stdexcept>

#include "kvstring.h"

namespace {

static const std::string whitespaces(" \t\f\v\n\r");

void strip_ws(std::string& str)
{
    size_t found = str.find_first_not_of(whitespaces);

    if (found != std::string::npos) {
        str.erase(0, found);

        found = str.find_last_not_of(whitespaces);
        if (found != std::string::npos)
            str.erase(found + 1);
    } else
        str.clear();
}

} // namespace

namespace kit {

void kv_parse(const std::string& str, std::string& k, double& v)
{
    k.clear();

    const size_t eq = str.find('=');
    if (eq == std::string::npos || str.find('=', eq + 1) != std::string::npos) {
        std::ostringstream oss;
        oss << "malformed KV pair '" << str << "'";
        throw std::runtime_error(oss.str());
    }

    std::istringstream iss(str.substr(eq + 1));
    iss >> v >> std::ws;

    if (iss.fail() || !iss.eof()) {
        std::ostringstream oss;
        oss << "KV pair '" << str << "' : V is not a real number";
        throw std::runtime_error(oss.str());
    }

    k = str.substr(0, eq);
    strip_ws(k);

    if (k.empty()) {
        std::ostringstream oss;
        oss << "KV pair '" << str << "' : K is empty";
        throw std::runtime_error(oss.str());
    }
}

//
// k1=v1:k2=v2:...:kN=vN
//

void kvstring_parse(const std::string& str, std::map<std::string, double>& kv)
{
    size_t b = 0;
    size_t e = str.find(':');

    while (b != std::string::npos) {

        const std::string kvstr = str.substr(b, e - b);

        if (!kvstr.empty()) {
            std::pair<std::string, double> kvp;

            kv_parse(kvstr, kvp.first, kvp.second);

            if (!kvp.first.empty()) {
                std::pair<std::map<std::string, double>::iterator, bool>
                    retval = kv.insert(kvp);

                if (!retval.second)
                    retval.first->second = kvp.second;
            }
        }

        if (e == std::string::npos)
            b = e;
        else
            b = e + 1;

        e = str.find(':', b);
    }
}

} // namespace kit
