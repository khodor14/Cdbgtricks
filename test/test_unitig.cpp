#include <gtest/gtest.h>
#include "../src/unitig.h"
#include <vector>
TEST(unitig,unused){
    std::string s="ACCCTGGACC";
    Unitig u=Unitig(s);
    EXPECT_EQ(u.get_left_unused(),0);
    EXPECT_EQ(u.get_right_unused(),2);
}
TEST(unitig,encoding){
    std::string s="ACCCTGGACC";
    Unitig u=Unitig(s);
    std::vector<uint8_t> vec={21,232,80};
    std::vector<uint8_t> enc=u.get_encoding();
    EXPECT_EQ(enc,vec);
}
TEST(unitig,length){
    std::string s="ACCCTGGACC";
    Unitig u=Unitig(s);
    EXPECT_EQ(u.unitig_length(),s.length());
    std::string seq="AACCTTAATGTAGGTGAATGGTA";
    Unitig u_Seq=Unitig(seq);
    EXPECT_EQ(u_Seq.unitig_length(),seq.length());
}
TEST(unitig,ith_mer){
    std::string s="ACCCTGGACC";
    Unitig u=Unitig(s);
    uint64_t i_thmer=u.get_ith_mer(0,5);//ACCCT=0001010111=87
    uint64_t h_ithmer=87;
    EXPECT_EQ(i_thmer,h_ithmer);

    i_thmer=u.get_ith_mer(1,5);//CCCTG=0101011110=350
    EXPECT_EQ(i_thmer,350);
}
TEST(unitig,next_mer){
    std::string s="ACCCTGGACC";
    Unitig u=Unitig(s);
    uint64_t i_thmer=u.get_ith_mer(0,5);//ACCCT=0001010111=87
    i_thmer=u.get_next_mer(i_thmer,1,5);//CCCTG=0101011110=350
    EXPECT_EQ(i_thmer,350);
    i_thmer=u.get_next_mer(i_thmer,2,5);//CCTGG=0101111010=378
    EXPECT_EQ(i_thmer,378);
}