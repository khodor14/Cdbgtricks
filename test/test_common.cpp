#include <gtest/gtest.h>
#include "../src/CommonUtils.h"
TEST(common,hash){
    std::string s="ACCTGGT";
    uint64_t h=hash(s);
    EXPECT_EQ(h,1515);
}
TEST(common,revcomp){
    uint64_t km=0;
    uint64_t revcmp=reverse_complement(km,32);
    EXPECT_EQ(revcmp,0xffffffffffffffff);
}