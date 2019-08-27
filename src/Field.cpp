#include <sstream>
#include "Field.h"
Field::Field(int nY, int nX, int nEns, float iFillValue) :
      mNY(nY), mNX(nX), mNEns(nEns) {
   if(Util::isValid(nY) && Util::isValid(nX) && Util::isValid(nEns)
         && nY >= 0 && nX >= 0 && nEns >= 0) {
      mValues.resize(nY*nX*nEns, iFillValue);
   }
   else {
      std::stringstream ss;
      ss << "Cannot create field of size [" << nY << "," << nX << "," << nEns << "]";
      Util::error(ss.str());
   }
}

std::vector<float> Field::operator()(unsigned int y, unsigned int x) const {
   std::vector<float> values(mValues.begin()+getIndex(y, x, 0), mValues.begin()+getIndex(y, x, mNEns-1)+1);
   return values;
}

int Field::getNumY() const {
   return mNY;
}
int Field::getNumX() const {
   return mNX;
}
int Field::getNumEns() const {
   return mNEns;
}

bool Field::operator==(const Field& iField) const {
   return mValues == iField.mValues;
}
bool Field::operator!=(const Field& iField) const {
   return mValues != iField.mValues;
}
