#ifndef FIELD_H
#define FIELD_H
#include <boost/shared_ptr.hpp>
#include <vector>
#include "Util.h"

//! Encapsulates gridded data in 3 dimensions: latitude, longitude, ensemble member.
//! Latitude generally represents the north-south direction and longitude the east-west, but the
//! grid does not necessarily need to follow a lat/lon grid. Any 2D grid will do.
// TODO: Rename latitude to x and longitude to y, as this is more generally correct.
class Field {
   public:
      //! Initialize 3D field
      //! @param nLat number of latitudes
      //! @param nLon number of longitudes
      //! @param nLon number of ensemble members
      //! @param iFillValue initialize all values in field with this
      Field(int nY, int nX, int nEns, float iFillValue=Util::MV);

      //! Access to data. Inlined for improved performance on some systems
      //! @param i latitude index
      //! @param j longitude index
      //! @param k ensemble index
      //! @return data at specified coordinate
      float      & operator()(unsigned int y, unsigned int x, unsigned int e) {
         int index = getIndex(y, x, e);
         return mValues[index];
      };
      float const& operator()(unsigned int y, unsigned int x, unsigned int e) const {
         int index = getIndex(y, x, e);
         return mValues[index];
      };

      //! Access to an ensemble for a specific grid point
      //! @param y y-axis index
      //! @param x x-axis index
      //! @return ensemble of values
      std::vector<float> operator()(unsigned int y, unsigned int x) const;

      //! Are all values (for all lat/lon/ens) in fields identical?
      bool operator==(const Field& iField) const;
      bool operator!=(const Field& iField) const;

      //! Number of gridpoints in the y direction
      int getNumY() const;

      //! Number of gridpoints in the x direction
      int getNumX() const;

      //! Number of ensemble members
      int getNumEns() const;
   private:
      //! Data values stored in a flat array. Index for ensemble changes fastest.
      std::vector<float> mValues;
      int mNY;
      int mNX;
      int mNEns;
      //! Index into flat array that corresponds to coordinate. Inlined for performance reasons
      int getIndex(unsigned int y, unsigned int x, unsigned int e) const {
         // Don't use an if statement, since this seems to be slow
         assert(y < mNY && x < mNX && e < mNEns);
         int index = e + x*mNEns + y*mNX*mNEns;
         return index;
      };
};
typedef boost::shared_ptr<Field> FieldPtr;
#endif
