#include <tuple>
#include <boost/serialization/serialization.hpp>
namespace boost {
namespace serialization {
template<class Archive, typename... Args>
void serialize(Archive & ar, std::tuple<Args...> & t, const unsigned int version)
{
    ar & std::get<0>(t);
    ar & std::get<1>(t);
    ar & std::get<2>(t);
}
}
}