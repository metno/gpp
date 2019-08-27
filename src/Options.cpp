#include "Options.h"
#include "Util.h"

Options::Options() {

}
void Options::addOption(std::string iKey, std::string iValue) {
   if(iKey == "" || iValue == "")
      return;

   // Overwrite existing value if key already exists
   for(int i = 0; i < mPairs.size(); i++) {
      std::string key = mPairs[i].first;
      if(key == iKey) {
         mPairs[i].second = iValue;
         return;
      }
   }
   mPairs.push_back(KeyValue(iKey, iValue));
}

void Options::addOptions(std::string iOptionString) {
   std::vector<std::string> options = Util::split(iOptionString);
   for(int i = 0; i < options.size(); i++) {
      addOption(options[i]);
   }
}
void Options::clear() {
   mPairs.clear();
}

std::string Options::toString() const {
   std::stringstream ss;
   for(int i = 0; i < mPairs.size(); i++) {
      if(i != 0)
         ss << " ";
      ss << mPairs[i].first << "=" << mPairs[i].second;
   }
   return ss.str();
}

Options::Options(std::string iOptionString) {
   addOptions(iOptionString);
}

void Options::addOption(std::string iOptionString) {
   int nextEqual = iOptionString.find('=');
   if(nextEqual < 0) {
      std::string key = iOptionString;
      std::string value = "1";
      addOption(key, value);
   }
   else {
      std::string key   = iOptionString.substr(0, nextEqual);
      std::string value = iOptionString.substr(nextEqual+1);
      addOption(key, value);
   }
}
bool Options::hasValue(const std::string& iKey) const {
   for(int i = 0; i < mPairs.size(); i++) {
      if(mPairs[i].first == iKey)
         return true;
   }
   return false;
}

bool Options::check() const {
   std::vector<std::string> unChecked;
   for(int i = 0; i < mPairs.size(); i++) {
      std::string key = mPairs[i].first;
      if(mHasBeenAccessed.find(key) == mHasBeenAccessed.end()) {
         unChecked.push_back(key);
      }
   }
   if(unChecked.size() > 0) {
      std::stringstream ss;
      if(unChecked.size() == 1) {
         ss << "The key '" << unChecked[0] << "' has";
      }
      else {
         ss << "The keys ";
         for(int i = 0; i < unChecked.size()-1; i++) {
            ss << "'" << unChecked[i] << "', ";
         }
         ss << "and '" << unChecked[unChecked.size()-1] << "' have";
      }
      ss <<  " never been used in '" << toString() << "'" << std::endl;
      Util::warning(ss.str());
      return false;
   }
   return true;
}
bool Options::operator==(const Options &right) const {
   for(int i = 0; i < mPairs.size(); i++) {
      std::string key = mPairs[i].first;
      std::string value = mPairs[i].second;
      if(!right.hasValue(key))
         return false;
      std::string s = "";
      right.getValue(key, s);
      if(value != s)
         return false;
   }
   for(int i = 0; i < right.mPairs.size(); i++) {
      std::string key = right.mPairs[i].first;
      std::string value = right.mPairs[i].second;
      if(!hasValue(key))
         return false;
      std::string s = "";
      getValue(key, s);
      if(value != s)
         return false;
   }
   return true;
}
bool Options::operator!=(const Options &right) const {
   return !(*this == right);
}
