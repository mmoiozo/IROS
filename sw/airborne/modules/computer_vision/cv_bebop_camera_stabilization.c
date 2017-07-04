/*
 * Copyright (C) Wilco Vlenterie (wv-tud)
 *
 * This file is part of paparazzi
 *
 * paparazzi is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2, or (at your option)
 * any later version.
 *
 * paparazzi is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with paparazzi; see the file COPYING.  If not, see
 * <http://www.gnu.org/licenses/>.
 */
/**
 * @file "modules/computer_vision//cv_bebop_camera_stabilization.c"
 * @author Wilco Vlenterie (wv-tud)
 * Active random sampling colour filter
 */

#include "modules/computer_vision/cv.h"
#include "modules/computer_vision/cv_bebop_camera_stabilization.h"
#include "modules/computer_vision/cv_image_pose.h"
#include "boards/bebop.h"                       ///< C header used for bebop specific settings
#include  <state.h>                              ///< C header used for state functions and data
#include "math/pprz_algebra_float.h"

#include BOARD_CONFIG

#ifndef CAM_STAB_CAMERA
#define CAM_STAB_CAMERA front_camera
#endif
PRINT_CONFIG_VAR(CAM_STAB_CAMERA);

#ifndef CAM_STAB_FOVY
#define CAM_STAB_FOVY 45.0
#endif
PRINT_CONFIG_VAR(CAM_STAB_FOVY);

#define PRINT(string,...) fprintf(stderr, "[CAM-STAB->%s()] " string,__FUNCTION__ , ##__VA_ARGS__)

#define CAM_STAB_VERBOSE FALSE
#if CAM_STAB_VERBOSE
#define VERBOSE_PRINT PRINT
#else
#define VERBOSE_PRINT(...)
#endif

struct image_t*   cv_cam_stab_func    ( struct image_t* img );
static uint16_t   horizontalLinePixel ( double x_min, double x_max, double y_angle, int8_t dir );
/** Fisheye correction **/
static double     correctRadius       ( double r, double f, double k );
static double     invertRadius        ( double r, double f, double k );
/** Perspective correction **/
static void       inputCoord          ( double x_out, double y_out, double *x_in, double *y_in );
static void       outputCoord         ( double x_in, double y_in, double *x_out, double *y_out );
/** Stabilization functions **/
static void       getMVP              ( struct FloatEulers*  eulerAngles, double MVP[16] );
static void       setPerspectiveMat   ( double m[16] );
static void       view_set_lookat     ( double result[16], double eye[4], double center[4], double up[4] );
/** Helper functions **/
static void       setRotationMat      ( float, float, float, double m[16] );
static void       setIdentityMat      ( double m[16] );
static void       translate_xyz       ( double result[16], const float translatex, const float translatey, const float translatez );
static void       float_mat4_mul      ( double m1[16], double m2[16], double result[16] );
static void       float_mat4vec_mul   ( double m1[16], double v1[4], double result[4] );
static void       rotateVector        ( float x, float y, float z, double vector[4] );


/** Set up perspective and correction parameters **/
double                      viewR               =  0.00105;                         ///< Perspective viewing distance
double                      default_k           =  1.2247445;                       ///< Fisheye correction factor (1.22474604174 max)
double                      firstOrder_comp     = -0.255;                           ///< First order remaing correction factor
double                      secondOrder_comp    =  0.155;                           ///< Second order counter remaing correction factor
double                      default_orbDiag     =  2 * CFG_MT9F002_FISHEYE_RADIUS;  ///< Diagonal size of the fisheye picture in full sensor resolution
float                       angleOfView         =  179.872;                         ///< Perspective angle of view ( < 180 )
static float                near                =  0.075;                           ///< Perspective near clipping plane
static float                far                 =  1.5;                             ///< Perspective far clipping plane
static double               crop_fovY           =  CAM_STAB_FOVY * M_PI / 180.0;    ///< (in Radians) FOV centered around the horizon to search for contours
/** Initialize parameters to be assigned during runtime **/
uint16_t                    ispWidth            = 0;                                ///< Maximum width of ISP after applied scaling
uint16_t                    ispHeight           = 0;                                ///< Maximum height of ISP after applied scaling
uint16_t                    initialWidth        = 0;                                ///< Initial width of ISP after applied scaling
uint16_t                    initialHeight       = 0;                                ///< Initial height of ISP after applied scaling
uint16_t                    cropCol             = 0;                                ///< Column from which the ISP is cropped relative to MIN
int16_t                     fillHeight          = 0;                                ///< Extra height used to make sure the horizon is in the centre
double                      ispScalar           = 1.0;                              ///< Applied scalar by the ISP
static uint16_t             runCount            = 0;                                ///< Total number of frames processed
static double               MVP[16];                                                ///< The model view projection matrix of the current frame

void cv_cam_stab_init() {
  ispScalar       = mt9f002.output_scaler * 2.0 / ( (double) mt9f002.y_odd_inc + 1.0 );
  ispHeight       = round( (CFG_MT9F002_X_ADDR_MAX - CFG_MT9F002_X_ADDR_MIN) * ispScalar );
  ispWidth        = round( (CFG_MT9F002_Y_ADDR_MAX - CFG_MT9F002_Y_ADDR_MIN) * ispScalar );
  initialWidth    = mt9f002.output_width;
  initialHeight   = mt9f002.output_height;
  cv_add_to_device(&CAM_STAB_CAMERA, cv_cam_stab_func);
}

struct image_t* cv_cam_stab_func(struct image_t* img)
{
  getMVP( &cv_image_pose.eulers, MVP );
  double   xAngle     = 75.0 / 180.0 * M_PI;
  uint16_t horizPos   = horizontalLinePixel(-xAngle, xAngle, 0.0, 0);
  uint16_t top        = horizontalLinePixel(-xAngle, xAngle,  0.5 * crop_fovY, +1);
  int16_t  bottom     = horizontalLinePixel(-xAngle, xAngle, -0.5 * crop_fovY, -1);
  uint16_t desHeight  = (top - bottom);
  // Center horizon by adding or removing extra space
  fillHeight          = (int16_t) round( 0.5*initialWidth - (horizPos - bottom) );
  bottom             -= fillHeight;
  desHeight          += fillHeight;
  // Check if we want a frame larger than initial output size
  if(desHeight > initialWidth){
    bottom           += (uint16_t) round((desHeight - initialWidth) / 2.0);
    desHeight         = initialWidth;
  }
  // Check if we want to start at negative pixels
  if(bottom < -CFG_MT9F002_X_ADDR_MIN * ispScalar){
    bottom              = -CFG_MT9F002_X_ADDR_MIN * ispScalar;
  }
  // Check if the desired height is larger than the sensor
  if(desHeight > ispHeight){
    bottom              = 0;
    desHeight           = ispHeight;
  }
  // Check if top of the requested frame is on the sensor
  if((desHeight + bottom) > ispHeight){
    desHeight           = ispHeight - bottom;
  }
  // See that fillheight is an even number
    if((fillHeight & 1) != 0){
      fillHeight--;
      bottom++;
      desHeight--;
    }
  // See that the bottom is an even number
  if((bottom & 1) != 0){
    bottom--;
    desHeight++;
  }
  // Set the chosen frame
  cropCol              = bottom + fillHeight;
  mt9f002.offset_x     = CFG_MT9F002_X_ADDR_MIN + ((uint16_t) bottom / ispScalar);
  mt9f002.output_width = desHeight;
  mt9f002.sensor_width = (uint16_t) (desHeight / ispScalar);
  if(desHeight > img->w || (desHeight - fillHeight) <= 0){
    fillHeight              = 0;
  }
  mt9f002_update_resolution(&mt9f002);
  runCount++;
  return NULL;
}

uint16_t horizontalLinePixel(double x_min, double x_max, double y_angle, int8_t dir){
  double x_in, y_in, center[2];
  angles2point(0.0, y_angle, &x_in, &y_in);
  point2pixel(x_in,y_in,&center[0],&center[1]);
  if(dir == 0){
    y_in = center[1];
  }
  else{
    double left[2], right[2];
    angles2point(x_min, y_angle, &x_in, &y_in);
    point2pixel(x_in,y_in,&left[0],&left[1]);

    angles2point(x_max, y_angle, &x_in, &y_in);
    point2pixel(x_in,y_in,&right[0],&right[1]);
    if(dir > 0){
      // Top of the line
      y_in = fmax( right[1], fmax( left[1], center[1]));
    }
    else{
      // Bottom of the line
      y_in = fmin( right[1], fmin( left[1], center[1]));
    }
  }
  return (uint16_t) round(y_in);
}

void pixel2point(double x_in, double y_in, double *x_out, double *y_out){
  double f, r, theta, corR, maxR;
  f       = default_orbDiag * ispScalar / (4 * sin(M_PI / 4));
  x_in   -= ispWidth * 0.5;
  y_in   -= ispHeight * 0.5;
  r       = sqrt( pow(x_in, 2.0) + pow(y_in, 2.0) );
  theta   = atan2( y_in, x_in );
  // Approximation of 2nd order inversion
  double r_start, r_end;
  r_start = r;
  r_end   = r_start * 1.0 / (1.0 + pow(r    , 1.0) * firstOrder_comp / pow(CFG_MT9F002_FISHEYE_RADIUS * ispScalar, 1.0) + pow(r    , 2.0) * secondOrder_comp / pow(CFG_MT9F002_FISHEYE_RADIUS * ispScalar, 2.0));
  r_end   = r_start * 1.0 / (1.0 + pow(r_end, 1.0) * firstOrder_comp / pow(CFG_MT9F002_FISHEYE_RADIUS * ispScalar, 1.0) + pow(r_end, 2.0) * secondOrder_comp / pow(CFG_MT9F002_FISHEYE_RADIUS * ispScalar, 2.0));
  r_end   = r_start * 1.0 / (1.0 + pow(r_end, 1.0) * firstOrder_comp / pow(CFG_MT9F002_FISHEYE_RADIUS * ispScalar, 1.0) + pow(r_end, 2.0) * secondOrder_comp / pow(CFG_MT9F002_FISHEYE_RADIUS * ispScalar, 2.0));
  corR    = correctRadius(r_end, f, default_k);
  maxR    = correctRadius(CFG_MT9F002_FISHEYE_RADIUS * ispScalar, f, default_k);
  x_in    = (corR * cos(theta)) / maxR;
  y_in    = (corR * sin(theta)) / maxR;
  outputCoord(x_in, y_in, y_out, x_out);
}

void point2pixel(double x_out, double y_out, double *x_in, double *y_in){
  double f, r, theta,corR, maxR;
  f           = default_orbDiag * ispScalar / (4 * sin(M_PI / 4));
  inputCoord(x_out, y_out, x_in, y_in);
  maxR        = correctRadius(CFG_MT9F002_FISHEYE_RADIUS * ispScalar, f, default_k);
  corR        = maxR * sqrt( pow( *x_in, 2.0 ) + pow( *y_in, 2.0 ) );
  theta       = atan2( *y_in, *x_in);
  r           = invertRadius(corR, f, default_k);
  r           = r * (1.0 + pow(r, 1.0) * firstOrder_comp / pow(CFG_MT9F002_FISHEYE_RADIUS * ispScalar, 1.0) + pow(r, 2.0) * secondOrder_comp / pow(CFG_MT9F002_FISHEYE_RADIUS * ispScalar, 2.0));
  *x_in       = ispWidth * 0.5 + cos(theta) * (r);
  *y_in       = ispHeight * 0.5 + sin(theta) * (r);
}

void angles2point(double xAngle, double yAngle, double *x_out, double * y_out){
  *x_out      = viewR * tan(xAngle);
  *y_out      = viewR * tan(yAngle);
}

void point2angles(double x_out, double y_out, double *xAngle, double *yAngle){

  *xAngle     = atan(x_out / viewR);
  *yAngle     = atan(y_out / viewR);
}

double correctRadius(double r, double f, double k){
  // This function calculates the corrected radius for radial distortion
  // According to the article "A Generic Non-Linear Method for Fisheye Correction" by Dhane, Kutty and Bangadkar
  double frac = sin( atan( r / f ) ) * k;
  if(frac < 1.0){
    return f * tan( asin( frac ) );
  }
  else{
    return f * tan( asin( sin( atan( CFG_MT9F002_FISHEYE_RADIUS * ispScalar / f ) ) * k ) );
  }
}

double invertRadius(double r, double f, double k){
  // This function calculates the inverted radius for corrected radial distortion
  // According to the article "A Generic Non-Linear Method for Fisheye Correction" by Dhane, Kutty and Bangadkar
  double frac = sin( atan( r / f ) ) / k;
  if(frac < 1.0){
    return f * tan( asin( frac ) );
  }
  else{
    return f * tan( asin( sin( atan( CFG_MT9F002_FISHEYE_RADIUS * ispScalar / f ) ) / k ) );
  }
}

void outputCoord(double x_in, double y_in, double *x_out, double *y_out){
  *x_out = (MVP[0] * x_in + MVP[4] * y_in + MVP[12]) / (MVP[3] * x_in + MVP[7] * y_in + MVP[15]);
  *y_out = (MVP[1] * x_in + MVP[5] * y_in + MVP[13]) / (MVP[3] * x_in + MVP[7] * y_in + MVP[15]);
}

void inputCoord(double x_out, double y_out, double *x_in, double *y_in){
  double den     = ( MVP[0]*MVP[5]  - MVP[0]*MVP[7]*y_out  - MVP[4]*MVP[1]+ MVP[4]*MVP[3]*y_out  + MVP[1]*MVP[7]*x_out  - MVP[5]*MVP[3]*x_out);

  double x_num   = (-MVP[4]*MVP[13] + MVP[4]*MVP[15]*y_out + MVP[12]*MVP[5]                  - MVP[12]*MVP[7]*y_out - MVP[5]*MVP[15]*x_out + MVP[13]*MVP[7]*x_out);
  double y_num   = (-MVP[0]*MVP[13] + MVP[0]*MVP[15]*y_out + MVP[12]*(MVP[1] - MVP[3]*y_out) - MVP[1]*MVP[15]*x_out + MVP[13]*MVP[3]*x_out);

  *x_in          =  x_num / (-den);
  *y_in          =  y_num / den;
}

void getMVP( struct FloatEulers*  eulerAngles, double MVP[16] ){
  double modelMat[16], viewMat[16], modelviewMat[16], projectionMat[16];
  double eye[4]       = {0.0, 0.0, -viewR, 1.0};
  double forward[4]   = {0.0, 0.0, 1.0, 1.0};
  double up[4]        = {0.0, 1.0, 0.0, 1.0};
  // Create The model matrix
  setRotationMat(0.0, M_PI, 0.0, modelMat);
  // Create the view matrix
  rotateVector(-MT9F002_THETA_OFFSET - eulerAngles->theta, 0.0, -eulerAngles->phi, forward);
  rotateVector(-MT9F002_THETA_OFFSET - eulerAngles->theta, 0.0, -eulerAngles->phi, up);
  double center[4] = {eye[0] + forward[0], eye[1] + forward[1], eye[2] + forward[2], 1.0};
  view_set_lookat(viewMat, eye,  center,  up);
  // Create the projection matrix
  setPerspectiveMat(projectionMat);
  float_mat4_mul(viewMat, modelMat, modelviewMat);
  float_mat4_mul(projectionMat, modelviewMat, MVP);
}

void setPerspectiveMat(double perspectiveMat[16])
{
  // Create perspective matrix
  float size      =  near * tanf((angleOfView / 180.0 * M_PI) / 2.0);
  // First Column
  perspectiveMat[0]  = 2 * near / (2 * size);
  perspectiveMat[1]  = 0.0;
  perspectiveMat[2]  = 0.0;
  perspectiveMat[3]  = 0.0;
  // Second Column
  perspectiveMat[4]  = 0.0;
  perspectiveMat[5]  = 2 * near / (2 * size);
  perspectiveMat[6]  = 0.0;
  perspectiveMat[7]  = 0.0;
  // Third Column
  perspectiveMat[8]  = 0.0;
  perspectiveMat[9]  = 0.0;
  perspectiveMat[10] = -(far + near) / (far - near);
  perspectiveMat[11] = -1;
  // Fourth Column
  perspectiveMat[12] = 0.0;
  perspectiveMat[13] = 0.0;
  perspectiveMat[14] = -(2 * far * near) / (far - near);
  perspectiveMat[15] = 0.0;
}

void view_set_lookat(double result[16], double eye[4], double center[4], double up[4]) {
  struct FloatVect3 f, s, u;
  f.x = center[0] - eye[0];
  f.y = center[1] - eye[1];
  f.z = center[2] - eye[2];
  // normalize f
  float_vect3_normalize(&f);
  // compute s = f x up (x means "cross product")
  s.x = f.y * up[2] - f.z * up[1];
  s.y = f.z * up[0] - f.x * up[2];
  s.z = f.x * up[1] - f.y * up[0];
  // and normalize s
  float_vect3_normalize(&s);
  //The up vector must not be parallel to the line of sight from the eye point to the reference point.
  if((0 == s.x)&&(0 == s.y)&&(0 == s.z))
    return;
  // compute u = s x f
  u.x = s.y * f.z - s.z * f.y;
  u.y = s.z * f.x - s.x * f.z;
  u.z = s.x * f.y - s.y * f.x;
  result[0] =  s.x;   result[4] =  s.y;   result[8]  =  s.z;  result[12] =  0.0f;
  result[1] =  u.x;   result[5] =  u.y;   result[9]  =  u.z;  result[13] =  0.0f;
  result[2] = -f.x;   result[6] = -f.y;   result[10] = -f.z;  result[14] =  0.0f;
  result[3] =  0.0f;  result[7] =  0.0f;  result[11] =  0.0f; result[15] =  1.0f;
  translate_xyz(result, -eye[0], -eye[1], -eye[2]);
}

void translate_xyz(double result[16], const float translatex,
    const float translatey, const float translatez) {
  result[12] += result[0] * translatex + result[4] * translatey
      + result[8] * translatez;
  result[13] += result[1] * translatex + result[5] * translatey
      + result[9] * translatez;
  result[14] += result[2] * translatex + result[6] * translatey
      + result[10] * translatez;
  result[15] += result[3] * translatex + result[7] * translatey
      + result[11] * translatez;
}

void rotateVector(float x, float y, float z, double vector[4])
{
  double rotMat[16], output[4];
  setRotationMat(x, y, z, rotMat);
  float_mat4vec_mul(rotMat, vector, output);
  vector[0] = output[0];
  vector[1] = output[1];
  vector[2] = output[2];
  vector[3] = output[3];
}


void float_mat4_mul(double m1[16], double m2[16], double result[16])
{
  // Fisrt Column
  result[0]  = m1[0]*m2[0] + m1[4]*m2[1] + m1[8]*m2[2] + m1[12]*m2[3];
  result[1]  = m1[1]*m2[0] + m1[5]*m2[1] + m1[9]*m2[2] + m1[13]*m2[3];
  result[2]  = m1[2]*m2[0] + m1[6]*m2[1] + m1[10]*m2[2] + m1[14]*m2[3];
  result[3]  = m1[3]*m2[0] + m1[7]*m2[1] + m1[11]*m2[2] + m1[15]*m2[3];
  // Second Column
  result[4]  = m1[0]*m2[4] + m1[4]*m2[5] + m1[8]*m2[6] + m1[12]*m2[7];
  result[5]  = m1[1]*m2[4] + m1[5]*m2[5] + m1[9]*m2[6] + m1[13]*m2[7];
  result[6]  = m1[2]*m2[4] + m1[6]*m2[5] + m1[10]*m2[6] + m1[14]*m2[7];
  result[7]  = m1[3]*m2[4] + m1[7]*m2[5] + m1[11]*m2[6] + m1[15]*m2[7];
  // Third Column
  result[8]  = m1[0]*m2[8] + m1[4]*m2[9] + m1[8]*m2[10] + m1[12]*m2[11];
  result[9]  = m1[1]*m2[8] + m1[5]*m2[9] + m1[9]*m2[10] + m1[13]*m2[11];
  result[10] = m1[2]*m2[8] + m1[6]*m2[9] + m1[10]*m2[10] + m1[14]*m2[11];
  result[11] = m1[3]*m2[8] + m1[7]*m2[9] + m1[11]*m2[10] + m1[15]*m2[11];
  // Fourth Column
  result[12] = m1[0]*m2[12] + m1[4]*m2[13] + m1[8]*m2[14] + m1[12]*m2[15];
  result[13] = m1[1]*m2[12] + m1[5]*m2[13] + m1[9]*m2[14] + m1[13]*m2[15];
  result[14] = m1[2]*m2[12] + m1[6]*m2[13] + m1[10]*m2[14] + m1[14]*m2[15];
  result[15] = m1[3]*m2[12] + m1[7]*m2[13] + m1[11]*m2[14] + m1[15]*m2[15];
}

void float_mat4vec_mul(double m1[16], double v1[4], double result[4])
{
  // Fisrt Column
  result[0]  = m1[0]*v1[0] + m1[4]*v1[1] + m1[8]*v1[2] + m1[12]*v1[3];
  result[1]  = m1[1]*v1[0] + m1[5]*v1[1] + m1[9]*v1[2] + m1[13]*v1[3];
  result[2]  = m1[2]*v1[0] + m1[6]*v1[1] + m1[10]*v1[2] + m1[14]*v1[3];
  result[3]  = m1[3]*v1[0] + m1[7]*v1[1] + m1[11]*v1[2] + m1[15]*v1[3];
}

void setIdentityMat(double m[16])
{
  m[0]  = m[5]  = m[10] = m[15] = 1.0;
  m[1]  = m[2]  = m[3]  = m[4]  = 0.0;
  m[6]  = m[7]  = m[8]  = m[9]  = 0.0;
  m[11] = m[12] = m[13] = m[14] = 0.0;
}

void setRotationMat(float x, float y, float z, double outputMat[16])
{
  double rotX[16];
  double rotY[16];
  double rotZ[16];
  setIdentityMat(outputMat);
  setIdentityMat(rotX);
  setIdentityMat(rotY);
  setIdentityMat(rotZ);
  // Set z rotation
  rotZ[0]  =  cosf(z);
  rotZ[1]  =  sinf(z);
  rotZ[4]  = -rotZ[1];
  rotZ[5]  =  rotZ[0];
  // Set y rotation
  rotY[0]  =  cosf(y);
  rotY[2]  =  sinf(y);
  rotY[8]  = -rotY[2];
  rotY[10] =  rotY[0];
  // Set x rotation
  rotX[5] = cosf(x);
  rotX[6] = -sinf(x);
  rotX[9] = -rotX[6];
  rotX[10] = rotX[5];
  // Multiply
  double rotYZ[16];
  float_mat4_mul(rotY, rotZ, rotYZ);
  float_mat4_mul(rotX, rotYZ, outputMat);
}
