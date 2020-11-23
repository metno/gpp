#include "gridpp.h"

using namespace gridpp;

float gridpp::dewpoint(float temperature, float relative_humidity) {
    if(gridpp::is_valid(temperature) && gridpp::is_valid(relative_humidity)) {
        // Taken from https://github.com/metno/wdb2ts
        float tempC = temperature - 273.15;
        float e = (relative_humidity)*0.611*exp( (17.63 * tempC) / (tempC + 243.04) );
        float tdC = (116.9 + 243.04 * log( e ))/(16.78-log( e ));
        float td = tdC + 273.15;
        return (td<=temperature ? td : temperature);
        // Taken from https://github.com/WFRT/comps
        // float es = 0.611*exp(5423.0*(1/273.15 - 1/(temperature)));
        // float e  = relative_humidity*(es);
        // float td = 1/(1/273.15 - 1.844e-4*log(e/0.611));
        // return td;
    }
    else
        return gridpp::MV;
}
vec gridpp::dewpoint(const vec& temperature, const vec& relative_humidity) {
    if(temperature.size() != relative_humidity.size())
        throw std::invalid_argument("Temperature and relative_humidity vectors are not the same size");

    vec ret(temperature.size());
    #pragma omp parallel for
    for(int y = 0; y < temperature.size(); y++) {
        ret[y] = gridpp::dewpoint(temperature[y], relative_humidity[y]);
    }
    return ret;
}
float gridpp::relative_humidity(float temperature, float dewpoint) {
    float mEwt[41] = { .000034,.000089,.000220,.000517,.001155,.002472,
                                .005080,.01005, .01921, .03553, .06356, .1111,
                                .1891,  .3139,  .5088,  .8070,  1.2540, 1.9118,
                                2.8627, 4.2148, 6.1078, 8.7192, 12.272, 17.044,
                                23.373, 31.671, 42.430, 56.236, 73.777, 95.855,
                                123.40, 157.46, 199.26, 250.16, 311.69, 385.56,
                                473.67, 578.09, 701.13, 845.28, 1013.25 };
   // Taken from https://github.com/metno/wdb2ts
   if(gridpp::is_valid(temperature) && gridpp::is_valid(dewpoint)) {
      if(temperature <= dewpoint)
         return 1;

      float x, et, etd, rh;
      int l;

      x = (temperature - 173.16) * 0.2;
      if(x < 0)
         x = 0;
      else if(x > 39)
         x = 39;
      l = int(x);
      assert(l >= 0);
      assert(l +1 < 41);
      et = mEwt[l] + (mEwt[l + 1] - mEwt[l]) * (x - float(l));

      x = (dewpoint - 173.16) * 0.2;
      if(x < 0)
         x = 0;
      else if(x > 39)
         x = 39;
      l = int(x);
      assert(l >= 0);
      assert(l +1 < 41);
      etd = mEwt[l] + (mEwt[l + 1] - mEwt[l]) * (x - float(l));
      rh = etd / et;

      if(rh < 0)
         rh = 0;
      if(rh > 1)
         rh = 1;

      return rh;
   }
   else
      return gridpp::MV;
}
vec gridpp::relative_humidity(const vec& temperature, const vec& dewpoint) {
    if(temperature.size() != dewpoint.size())
        throw std::invalid_argument("Temperature and dewpoint vectors are not the same size");

    vec ret(temperature.size());
    #pragma omp parallel for
    for(int y = 0; y < temperature.size(); y++) {
        ret[y] = gridpp::relative_humidity(temperature[y], dewpoint[y]);
    }
    return ret;
}
float gridpp::wetbulb(float temperature, float pressure, float relative_humidity) {
   float temperatureC = temperature - 273.15;
   if(temperatureC <= -243.04 || relative_humidity <= 0)
      return gridpp::MV;
   if(gridpp::is_valid(temperatureC) && gridpp::is_valid(pressure) && gridpp::is_valid(relative_humidity)) {
      float e  = (relative_humidity)*0.611*exp((17.63*temperatureC)/(temperatureC+243.04));
      float Td = (116.9 + 243.04*log(e))/(16.78-log(e));
      float gamma = 0.00066 * pressure/1000;
      float delta = (4098*e)/pow(Td+243.04,2);
      if(gamma + delta == 0)
         return gridpp::MV;
      float wetbulbTemperature   = (gamma * temperatureC + delta * Td)/(gamma + delta);
      float wetbulbTemperatureK  = wetbulbTemperature + 273.15;
      return wetbulbTemperatureK;
   }
   else {
      return gridpp::MV;
   }
}
vec gridpp::wetbulb(const vec& temperature, const vec& pressure, const vec& relative_humidity){
    if(temperature.size() != pressure.size())
        throw std::invalid_argument("Temperature and pressure vectors are not the same size");
    if(temperature.size() != relative_humidity.size())
        throw std::invalid_argument("Temperature and relative_humidity vectors are not the same size");

    vec ret(temperature.size());
    #pragma omp parallel for
    for(int y = 0; y < temperature.size(); y++) {
        ret[y] = gridpp::wetbulb(temperature[y], pressure[y], relative_humidity[y]);
    }
    return ret;
}
