#ifndef FILL_VALUE_H
#define FILL_VALUE_H

// If we were using C++14 we could use variable templates
// but this work around keeps things from getting too verbose
template<typename T> static constexpr T fill_value();

template<> constexpr float fill_value()              { return -999999.0F; }
template<> constexpr double fill_value()             { return -999999.0; }
template<> constexpr signed char fill_value()        { return -127; }
template<> constexpr short fill_value()              { return -32767; }
template<> constexpr int fill_value()                { return -2147483647; }
template<> constexpr long long fill_value()          { return -9223372036854775807LL; }
template<> constexpr unsigned char fill_value()      { return 254; }
template<> constexpr unsigned short fill_value()     { return 65534; }      
template<> constexpr unsigned int fill_value()       { return 4294967294U; }
template<> constexpr unsigned long long fill_value() { return 18446744073709551614ULL; }

#endif
