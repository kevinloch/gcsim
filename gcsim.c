/*
 * Copyright (c) 2019, Kevin M. Loch
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the <organization> nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "dalsgaard.h"

//#define WSPHERE // 1kg tungsten sphere
#define EARTH
//#define SUN
//#define MOON  // do not yet have density model for moon

#ifdef WSPHERE
const char body[]=           "1kg tungsten sphere";
const char densitymodel[]=   "Uniform density";
const double r_body=         0.02312694;  // meters, radius
const double p_body=         19300.0;     // kg/m^3, density
const double m_body=         1.0;         // kg,     mass
const double samplelength=   0.000005;    // meters, sampling scale of body
const double altitude=       0.9768731;   // meters, 1m from center
const double sensor_angle=   0.001000;    // degrees
#endif

#ifdef EARTH
const char body[]=           "Earth";
const char densitymodel[]=   "PREM";
const double r_body=         6371.0E3;    // meters, earth mean radius
const double p_body=         5514.0;      // kg/m^3, earth mean density
const double m_body=         6.97237E24;  // kg,     earth mass
const double samplelength=   1000.0;      // meters, sampling scale of body

const double altitude=       1.0E0;       // meters, surface
const double sensor_angle=   0.1000;      // degrees
//const double altitude=     362600.0E3;  // meters, lunar distance
//const double sensor_angle= 0.0010;      // degrees
#endif

#ifdef SUN
const char body[]=           "Sun";
const char densitymodel[]=   "Christensen-Dalsgaard";
const double r_body=         695700E3;    // meters, sun mean radius
const double p_body=         1408;        // kg/m^3, sun mean density
const double m_body=         7.342E22;    // kg,     sun mass
const double samplelength=   250000.0;    // meters, sampling scale of body

const double altitude=       1.496E11;    // meters, 1 AU
const double sensor_angle=   0.0005;      // degrees
//const double altitude=     1.0E0;       // surface
//const double sensor_angle= 0.1;         // degrees
#endif

#ifdef MOON
const char body[]=           "Moon - internal density model needed";
const char densitymodel[]=   "*Model Needed* using uniform density";
const double r_body=         1737.1E3;    // meters, moon mean radius
const double p_body=         3344.0;      // kg/m^3, moon mean density
const double m_body=         7.342E22;    // kg,     moon mass
const double samplelength=   250.0;       // meters, sampling scale of body

const double altitude=       362600.0E3;  // meters, lunar distance
const double sensor_angle=   0.0005;      // degrees
//const double altitude=     1.0E0;       // surface
//const double sensor_angle= 0.1;         // degrees
#endif

// reference values
const double G_ref=          6.67408E-11; // codata
const double testmass=       1.000; // kg  - keep this at 1.0 for standard intensity or plots will be wrong

typedef struct {
  double theta;
  long double force;
  long double forcez;
  double spotratio;
} sensorring;

double getBodyDensity(double r) {
  double p; // density
  double x; // normalized radius;
  int i;

#ifdef EARTH
  // PREM model
  x=r/6371.0E3;
  if (r < 1221.5E3) {
    p=1.0E3 * (13.0885 - (8.8381 * pow(x, 2.0)));
    return(p);
  } else if (r < 3480.0E3) {
    p=1.0E3 * (12.5815 - (1.2638 * x) - (3.6426 * pow(x, 2.0)) - (5.5281 * pow(x, 3.0)));
    return(p);
  } else if (r < 5701.0E3) {
    p=1.0E3 * (7.9565 - (6.4761 * x) + (5.5283 * pow(x, 2.0)) - (3.0807 * pow(x, 3.0))); 
    return(p);
  } else if (r < 5771.0E3) {
    p=1.0E3 * (5.3197 - (1.4836 * x));
    return(p);
  } else if (r < 5971.0E3) {
    p=1.0E3 * (11.2494 - (8.0298 * x)); 
    return(p);
  } else if (r < 6151.0E3) {
    p=1.0E3 * (7.1089 - (3.8045 * x));
    return(p);
  } else if (r < 6346.6E3) {
    p=1.0E3 * (2.691 + (0.6924 * x));
    return(p);
  } else if (r < 6356.0E3) {
    return(2900.0);
  } else if (r < 6368.0E3) {
    return(2600.0); 
  } else {
    return(1020.0);
  }
#endif

#ifdef SUN
  // Dalsgaard Model 'S' http://users-phys.au.dk/jcd/solar_models/
  x=r/r_body;
  for (i=0; i < 2482; i++) {
    if (x >= solardensity[i].radius) {
      p=1000.0 * solardensity[i].density;
      //printf("radius: %.6e, density: %6e\n", x, p);
      return(p);
    }
  }
#endif

  // default return mean density
  return(p_body);
};

int main() {
  long int i;
  int x,y,z;
  int xybegin, xyend;
  int zend;
  int ecenterz;
  double rcenter, rexp;
  double theta;
  double density;
  double force;
  struct timespec starttime;
  struct timespec endtime;
  struct timespec totstarttime;
  struct timespec totendtime;
  double elapsedtime;
  double xm, ym, zm, ecenterzm;
  double rxy;
  long double totalforce;
  double forcez;
  long double totalforcez;
  double expz; // height of text mass above body surface
  int stopxy;
  double sl2;
  double ringvolume, ringmass;
  long double exactforcez;
  double angular_size;
  sensorring *sensorrings;
  int outputslots;
  double xmi, xmo;
  double rexpi, rexpo;
  double samplethetai, samplethetao;
  double outtheta;
  int outindex;
  int j;
  long double intensity_factor;
  double sa, sa2;
  double sensorthetai, sensorthetao;
  int overlap;
  int outindexi, outindexo;
  double sampleringarea;
  double sensorringarea;

  sl2 = samplelength / 2.0;
  sa = M_PI * sensor_angle / 180.0; // sensor angle in radians
  sa2 = sa / 2.0;
  intensity_factor = 1.0 / pow((sensor_angle / 2.0), 2.0); // to normalize to deg^2 // only valid if sensor angle is smaller than angular size

  xyend=(int)((r_body - ((double)samplelength / 2.0)) / (double)samplelength);
  ecenterz=-xyend; 
  ecenterzm=(double)(ecenterz * samplelength);
  zend=-2 * xyend;
  expz=altitude + ((double)samplelength / 2.0);


  // determine angular size from sensor
  angular_size=2.0 *  asin(r_body / (r_body + altitude));
  // allocate output buffer
  outputslots = (int)((angular_size / (2.0 * sa)) + 2.0);
  sensorrings = malloc(outputslots * sizeof(sensorring));  
  if (sensorrings == NULL) {
    printf("error, could not allocate output buffer\n");
    exit(1);
  }
  for (j=0; j<outputslots; j++) {
    sensorrings[j].theta=0.0;
    sensorrings[j].force=0.0;
    sensorrings[j].forcez=0.0;
    if (j == 0) {
      sensorrings[j].spotratio=1.0;
    } else {
      sensorrings[j].spotratio=(1.0 / (pow(((2.0 * j) + 1.0), 2.0) - pow(((2.0 * j) - 1.0), 2.0))); 
    }
  }

  printf("sttaus, Observing: %s\n", body);
  printf("status, Density model: %s\n", densitymodel);
  printf("status, Sampling scale: %.3e m\n", samplelength);
  printf("status, Altitude above surface: %.3e m\n", altitude);
  printf("status, Test mass: %.3e kg\n", testmass);
  printf("status, Angular resolution: %.3e deg\n", sensor_angle);
  printf("status, xyend: %d, ecenterz: %d, zend: %d, angular size: %.3e, output slots: %d\n", xyend, ecenterz, zend, (180.0 * angular_size / M_PI), outputslots);

/*
 * Quantize body by sampling nested/stacked cylinders and rings and calculate force from each to sensor test mass
 */

  i=0;
  totalforce=0.0;
  totalforcez=0.0;
  clock_gettime(CLOCK_REALTIME, &totstarttime);
  clock_gettime(CLOCK_REALTIME, &starttime);
  for (z=0; z >= zend; z--) {
    zm=((double)z * samplelength);
    stopxy=xyend;
    for (x=0; x <= stopxy; x++) {
      xm=((double)x * samplelength);
      
      // sample distance from center of body
      rcenter=hypot(xm, (ecenterzm-zm));

      if (rcenter <= r_body) {
        i++;
        // sample distance from experiment
        rexp=hypot(xm, (zm-expz));

        // sample angle from experiment (vs center of body)
        theta=acos((zm-expz) / -rexp);

        // cylinder or ring volume
        if (x == 0) { // cylinder
          ringvolume=M_PI * samplelength * pow(sl2, 2.0);
        } else { // ring
          ringvolume=M_PI * samplelength * (pow((xm + sl2), 2.0) - pow((xm - sl2), 2.0));
        }
        
        // mass of cylinder or ring
        density=getBodyDensity(rcenter);
        ringmass=ringvolume * density;
 
        // force from cylinder or ring
        force=G_ref * testmass * ringmass / pow(rexp, 2.0);
        totalforce = totalforce + (long double)force;
        forcez=cos(theta) * force;
        totalforcez=totalforcez + (long double)forcez;

/*
 *  now map sample ring force to sensor "rings"
 *  these ficticious sensor rings will be converted to normalized sensor spot apertures during output
 */
        //  find inner and outer sample ring theta
        xmi=xm - sl2;
        if (xmi < 0) {  // special cyliner case
          xmi=0;
        }
        rexpi=hypot(xmi, (zm-expz));
        samplethetai=acos((zm-expz) / -rexpi);
        xmo=xm + sl2;
        rexpo=hypot(xmo, (zm-expz));
        samplethetao=acos((zm-expz) / -rexpo);

        // find start and end sendor ring indexes
        outindexi=(int)((samplethetai / sa) + 0.5);
        outindexo=(int)((samplethetao / sa) + 0.5);

        if (outindexi == outindexo) {  // sample ring is entirely within one sensor ring
          outindex=outindexo;
          sensorrings[outindex].theta=(outindex * sa) + sa2;
          sensorrings[outindex].force=sensorrings[outindex].force + (long double)force;
          sensorrings[outindex].forcez=sensorrings[outindex].forcez + (long double)forcez;
        } else { 
          // assign sample ring force proportionally to sensor rings by area
          sampleringarea=pow(samplethetao, 2.0) - pow(samplethetai, 2.0); // we can ignore pi as this is only used in ratios with other rings
          for (j=outindexi; j <=outindexo; j++) {
            if (j == outindexi) { // first ring
              sensorringarea=pow(((j * sa) + sa2), 2.0) - pow(samplethetai, 2.0);
              sensorrings[j].theta=(j * sa) + sa2;
              sensorrings[j].force=sensorrings[j].force + ((long double)force * (long double)sensorringarea / (long double)sampleringarea);
              sensorrings[j].forcez=sensorrings[j].forcez + ((long double)forcez * (long double)sensorringarea / (long double)sampleringarea);
            } else if (j == outindexo) { // last ring
              sensorringarea=pow(samplethetao, 2.0) - pow(((j * sa) - sa2), 2.0);
              sensorrings[j].theta=(j * sa) + sa2;
              sensorrings[j].force=sensorrings[j].force + ((long double)force * (long double)sensorringarea / (long double)sampleringarea);
              sensorrings[j].forcez=sensorrings[j].forcez + ((long double)forcez * (long double)sensorringarea / (long double)sampleringarea);
            } else { // middle
              sensorringarea=pow(((j * sa) + sa2), 2.0) - pow(((j * sa) - sa2), 2.0);
              sensorrings[j].theta=(j * sa) + sa2;
              sensorrings[j].force=sensorrings[j].force + ((long double)force * (long double)sensorringarea / (long double)sampleringarea);
              sensorrings[j].forcez=sensorrings[j].forcez + ((long double)forcez * (long double)sensorringarea / (long double)sampleringarea);
            }
          }
        }

        // find indexed sensor theta and inner/outer theta limis for that index
        outindex=(int)((theta / sa) + 0.5);
        sensorrings[outindex].theta=outindex * sa + sa2;
        sensorthetai = sensorrings[outindex].theta - sa2;
        if (sensorthetai < 0) {  // special cyliner case
          sensorthetai = 0;
        }
        sensorthetao = sensorrings[outindex].theta + sa2;
      } else {
        if (stopxy == xyend) {
          stopxy=x;
        }
      }
    }
  }
  clock_gettime(CLOCK_REALTIME, &totendtime);
  elapsedtime=((double)(totendtime.tv_sec - 1500000000) + ((double)totendtime.tv_nsec / 1.0E9)) - ((double)(totstarttime.tv_sec - 1500000000) + ((double)totstarttime.tv_nsec) / 1.0E9);

  exactforcez=G_ref * (4.0 * M_PI * pow(r_body, 3.0) / 3.0) * p_body * testmass / pow((r_body + altitude), 2.0);
  printf("status, samples: %ld,  totalforce: %.6Le, totalforcez: %.6Le, exactforcez: %.6Le (%6.6fs)\n", i, totalforce, totalforcez, exactforcez, elapsedtime);

  printf("status, index, angle, force, forcez, spotratio, spotforce, spotforcez, normalized spotforce, normalized spotforcez\n");
  for (j=0; j<outputslots; j++) {
    if (sensorrings[j].force != 0.0) {
      printf("result, %d, %.6e, %.6Le, %.6Le, %.6e, %.6Le, %.6Le, %.6Le, %.6Le\n", j, (180.0 * sensorrings[j].theta / M_PI), sensorrings[j].force, sensorrings[j].forcez, sensorrings[j].spotratio\
, (sensorrings[j].force * sensorrings[j].spotratio), (sensorrings[j].forcez * sensorrings[j].spotratio), (sensorrings[j].force * sensorrings[j].spotratio * intensity_factor), (sensorrings[j].forcez * sensorrings[j].spotratio) * intensity_factor);
    }
  }

  return(0); 
}
