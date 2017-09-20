#ifndef KVSTRING_H
#define KVSTRING_H

#include <map>
#include <string>

namespace kit {

void kv_parse(const std::string&, std::string&, double&);
void kvstring_parse(const std::string&, std::map<std::string, double>&);

} // namespace kit

#endif // KVSTRING_H
