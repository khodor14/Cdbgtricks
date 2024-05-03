#include <tuple>

#include <unordered_map>

#include "unitig.h"

#ifndef FileSerializer_H
#define FileSerializer_H
namespace boost {
  namespace serialization {
    template < class Archive, typename...Args >
      void serialize(Archive & ar, std::tuple < Args... > & t,
        const unsigned int version) {
        ar & std::get < 0 > (t);
        ar & std::get < 1 > (t);
        ar & std::get < 2 > (t);
      }
  }
}
#endif