#include <cstdio>
class FileSerializer 
{
public:
    // serialize basic types to FILE
    // -----------------------------
    
    template <class T>
    bool operator()(FILE *fp, const T& value) 
    {
        return fwrite((const void *)&value, sizeof(value), 1, fp) == 1;
    }

    template <class T>
    bool operator()(FILE *fp, T* value) 
    {
        return fread((void *)value, sizeof(*value), 1, fp) == 1;
    }
    // serialize std::pair<const A, B> to FILE - needed for maps
    // ---------------------------------------------------------
    template <class A, class B>
    bool operator()(FILE *fp, const std::pair<const A, B>& value)
    {
        return (*this)(fp, value.first) && (*this)(fp, value.second);
    }

    template <class A, class B>
    bool operator()(FILE *fp, std::pair<const A, B> *value) 
    {
        return (*this)(fp, (A *)&value->first) && (*this)(fp, &value->second);
    }
};
