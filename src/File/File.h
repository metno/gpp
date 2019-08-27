#ifndef FILE_H
#define FILE_H
#include <vector>
#include <map>
#include <boost/shared_ptr.hpp>
#include "../Variable.h"
#include "../Uuid.h"
#include "../Field.h"

class Options;

// 3D array of data: [y][x][ensemble_member]
typedef std::vector<std::vector<float> > vec2; // Y, X

//! Represents a data file containing spatial and temporal data initialized at one particular time
class File {
   public:
      File(std::string iFilename, const Options& iOptions);
      virtual ~File();

      void setVariables(std::vector<Variable> iVariables);

      //! Insantiates a file. Returns null if file does not exist or cannot be parsed.
      static File* getScheme(std::string iFilename, const Options& iOptions, bool iReadOnly=false);

      //! Get a handle to a field
      //! @param iSkipRead If 'true', don't read the field, just initialize a field with missing values
      FieldPtr getField(const Variable& iVariable, int iTime, bool iSkipRead=false) const;
      FieldPtr getField(std::string iVariable, int iTime) const;

      //! Get a new field initialized with missing values. Don't register the field to any variable.
      FieldPtr getEmptyField(float iFillValue=Util::MV) const;

      // Register a field to the file, overwriting existing ones (if necessary)
      void addField(FieldPtr iField, const Variable& iVariable, int iTime) const;

      // Write these variables to file
      void write(std::vector<Variable> iVariables, std::string iMessage="");

      // Dimension sizes

      // Latitudes and longitudes must be set early on in the code, as this sets
      // the X and Y dimension
      bool setLats(vec2 iLats);
      bool setLons(vec2 iLons);
      int getNumY() const;
      int getNumX() const;
      int getNumEns() const;
      int getNumTime() const;
      vec2 getLats() const;
      vec2 getLons() const;
      vec2 getElevs() const;
      vec2 getLandFractions() const;
      bool setElevs(vec2 iElevs);
      bool setLandFractions(vec2 iLandFractions);
      bool setNumEns(int iNum);

      //! Does this file provide the variable (deriving it if necessary)?
      bool hasVariable(const Variable& iVariable) const;

      std::string getFilename() const;
      bool hasDefinedVariable(Variable iVariable) const;

      bool hasSameDimensions(const File& iOther) const;
      std::string getDimenionString() const;
      void initNewVariable(const Variable& iVariable);
      virtual std::string name() const = 0;
      //! Clear the retrieved/computed fields stored in cache
      void clear();
      //! How many bytes of retrieved/computed  data are stored in cache?
      //! @return Number of bytes
      long getCacheSize() const;

      //! Returns a tag that uniquely identifies the latitude/longitude grid
      //! If the grid changes, a new tag is issued. Two files with the same grid
      //! will not have the same unique tag.
      Uuid getUniqueTag() const;

      //! Set the time that the file is issued
      //! @ param iTime The number of seconds since 1970-01-01 00:00:00 +00:00
      void setReferenceTime(double iTime);
      double getReferenceTime() const;
      //! Set the times of the timesteps in the file
      //! @ param iTimes vector of number of seconds since 1970-01-01 00:00:00 +00:00
      void setTimes(std::vector<double> iTimes);
      std::vector<double> getTimes() const;
      bool getVariable(std::string iVariableName, Variable& iVariable) const;
      static std::string getDescriptions();
      void addVariableAlias(std::string iAlias, Variable iVariable);
      bool hasElevs() const;
   protected:
      virtual FieldPtr getFieldCore(const Variable& iVariable, int iTime) const = 0;
      // File must save variables, but also altitudes, in case they got changed
      virtual void writeCore(std::vector<Variable> iVariables, std::string iMessage="") = 0;
      //! Does the subclass provide this variable without deriving it?
      virtual bool hasVariableCore(const Variable& iVariable) const = 0;

      // Subclasses must fill these fields in the constructor:
      vec2 mLandFractions;
      int mNEns;

      //! These must be populated on initialization by subclass
      mutable std::vector<Variable> mVariables;
      std::map<std::string, Variable> mVariableAliases;
   private:
      std::string mFilename;
      mutable std::map<Variable, std::vector<FieldPtr> > mFields;  // Variable, offset
      mutable Uuid mTag;
      void createNewTag() const;
      FieldPtr getEmptyField(int nY, int nX, int nEns, float iFillValue=Util::MV) const;
      double mReferenceTime;
      std::vector<double> mTimes;
      static Uuid mNextTag;
      vec2 mElevs;
      bool mHasElevs;
      vec2 mLats;
      vec2 mLons;
};
#include "Netcdf.h"
#include "Fake.h"
#include "Point.h"
#include "NorcomQnh.h"
#include "Text.h"
#endif
