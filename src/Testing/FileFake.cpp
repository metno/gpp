#include "../File/Fake.h"
#include "../Util.h"
#include "../Downscaler/Downscaler.h"
#include <gtest/gtest.h>

namespace {
   class FileFakeTest : public ::testing::Test {
   };

   TEST_F(FileFakeTest, getLatLons) {
      FileFake file(Options("nLat=3 nLon=3 nEns=1 nTime=1"));
      vec2 lats = file.getLats();
      vec2 lons = file.getLons();
      ASSERT_EQ(3, lats.size());
      ASSERT_EQ(3, lats[0].size());
      ASSERT_EQ(3, lons.size());
      ASSERT_EQ(3, lons[0].size());
   }
   TEST_F(FileFakeTest, setLatLons) {
      FileFake file(Options("nLat=3 nLon=3 nEns=1 nTime=1"));
      vec2 origLats = file.getLats();
      vec2 origLons = file.getLons();
      vec2 origElevs = file.getElevs();
      vec2 lats(3);
      vec2 lons(3);
      vec2 elevs(3);
      for(int i = 0; i < 3; i++) {
         lats[i].resize(3, 1.9);
         lons[i].resize(3, 3.1);
         elevs[i].resize(3, 4.1);
      }

      EXPECT_TRUE(file.setLats(lats));
      EXPECT_TRUE(file.setLons(lons));
      EXPECT_TRUE(file.setElevs(elevs));

      vec2 newLats = file.getLats();
      vec2 newLons = file.getLons();
      vec2 newElevs = file.getElevs();

      ASSERT_EQ(3,  newLats.size());
      ASSERT_EQ(3,  newLons.size());
      ASSERT_EQ(3,  newElevs.size());
      ASSERT_EQ(3,  newLats[0].size());
      ASSERT_EQ(3,  newLons[0].size());
      ASSERT_EQ(3,  newElevs[0].size());
      for(int i = 0; i < origLats.size(); i++) {
         for(int j = 0; j < origLons.size(); j++) {
            EXPECT_EQ(lats[i][j], newLats[i][j]);
            EXPECT_EQ(lons[i][j], newLons[i][j]);
            EXPECT_EQ(elevs[i][j], newElevs[i][j]);
         }
      }
   }
   TEST_F(FileFakeTest, setLatLonsInvalid) {
      FileFake file(Options("nLat=3 nLon=3 nEns=1 nTime=1"));
      vec2 origLats = file.getLats();
      vec2 origLons = file.getLons();
      vec2 origElevs = file.getElevs();
      vec2 lats(2);
      vec2 lons(2);
      vec2 elevs(2);
      for(int i = 0; i < 2; i++) {
         lats[i].resize(4, 3);
         lons[i].resize(4, 3);
         elevs[i].resize(4, 3);
      }

      EXPECT_FALSE(file.setLats(lats));
      EXPECT_FALSE(file.setLons(lons));
      EXPECT_FALSE(file.setElevs(elevs));

      vec2 newLats = file.getLats();
      vec2 newLons = file.getLons();
      vec2 newElevs = file.getElevs();
      ASSERT_EQ(origLats.size(),  newLats.size());
      ASSERT_EQ(origLons.size(),  newLons.size());
      ASSERT_EQ(origElevs.size(), newElevs.size());
      for(int i = 0; i < origLats.size(); i++) {
         for(int j = 0; j < origLons.size(); j++) {
            EXPECT_EQ(origLats[i][j], newLats[i][j]);
            EXPECT_EQ(origLons[i][j], newLons[i][j]);
            EXPECT_EQ(origElevs[i][j], newElevs[i][j]);
         }
      }
   }
   TEST_F(FileFakeTest, getDimensionSizes) {
      FileFake file(Options("nLat=4 nLon=2 nEns=5 nTime=3"));
      EXPECT_EQ(4, file.getNumY());
      EXPECT_EQ(2, file.getNumX());
      EXPECT_EQ(5, file.getNumEns());
      EXPECT_EQ(3, file.getNumTime());
   }
   TEST_F(FileFakeTest, invalid) {
      ::testing::FLAGS_gtest_death_test_style = "threadsafe";
      Util::setShowError(false);
      EXPECT_DEATH(FileFake file(Options("nLat=0 nLon=1 nEns=2 nTime=3")), ".*");
      EXPECT_DEATH(FileFake file(Options("nLat=1 nLon=0 nEns=2 nTime=3")), ".*");
      EXPECT_DEATH(FileFake file(Options("nLat=1 nLon=1 nEns=0 nTime=3")), ".*");
      EXPECT_DEATH(FileFake file(Options("nLat=1 nLon=1 nEns=1 nTime=0")), ".*");
      EXPECT_DEATH(FileFake file(Options("nLat=-999 nLon=1 nEns=1 nTime=1")), ".*");
      EXPECT_DEATH(FileFake file(Options("nLat=1 nLon=-999 nEns=1 nTime=1")), ".*");
      EXPECT_DEATH(FileFake file(Options("nLat=1 nLon=1 nEns=-999 nTime=1")), ".*");
      EXPECT_DEATH(FileFake file(Options("nLat=1 nLon=1 nEns=1 nTime=-999")), ".*");
      EXPECT_DEATH(FileFake file(Options("nLat=1 nLon=0 nEns=1 nTime=-999")), ".*");
   }
}
int main(int argc, char **argv) {
     ::testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
