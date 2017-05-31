/*
 * Copyright (C) 2014 Hann Woei Ho
 *               2015 Freek van Tienen <freek.v.tienen@gmail.com>
 *               2016 Kimberly McGuire <k.n.mcguire@tudelft.nl
 *
 * This file is part of Paparazzi.
 *
 * Paparazzi is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2, or (at your option)
 * any later version.
 *
 * Paparazzi is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Paparazzi; see the file COPYING.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

/**
 * @file modules/computer_vision/opticflow/opticflow_calculator.c
 * @brief Estimate velocity from optic flow.
 *
 * Using images from a vertical camera and IMU sensor data.
 */

#include "std.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

// Own Header
#include "flow_speed_calculator.h"

// Computer Vision
#include "lib/vision/image.h"
#include "lib/vision/lucas_kanade.h"
#include "lib/vision/fast_rosten.h"
#include "lib/vision/edge_flow.h"
#include "size_divergence.h"
#include "linear_flow_fit.h"

#include "math/pprz_algebra.h"
#include "math/pprz_algebra_float.h"

#define PI 3.1415926

// whether to show the flow and corners:
#define OPTICFLOW_SHOW_FLOW 1
#define OPTICFLOW_SHOW_CORNERS 0

// What methods are run to determine divergence, lateral flow, etc.
// SIZE_DIV looks at line sizes and only calculates divergence
#define SIZE_DIV 1
// LINEAR_FIT makes a linear optical flow field fit and extracts a lot of information:
// relative velocities in x, y, z (divergence / time to contact), the slope of the surface, and the surface roughness.
#define LINEAR_FIT 1

// Camera parameters (defaults are from an ARDrone 2)
#ifndef OPTICFLOW_FOV_W
#define OPTICFLOW_FOV_W 0.89360857702
#endif
PRINT_CONFIG_VAR(OPTICFLOW_FOV_W)

#ifndef OPTICFLOW_FOV_H
#define OPTICFLOW_FOV_H 0.67020643276
#endif
PRINT_CONFIG_VAR(OPTICFLOW_FOV_H)

#ifndef OPTICFLOW_FX
#define OPTICFLOW_FX 343.1211
#endif
PRINT_CONFIG_VAR(OPTICFLOW_FX)

#ifndef OPTICFLOW_FY
#define OPTICFLOW_FY 348.5053
#endif
PRINT_CONFIG_VAR(OPTICFLOW_FY)

/* Set the default values */
#ifndef OPTICFLOW_MAX_TRACK_CORNERS
#define OPTICFLOW_MAX_TRACK_CORNERS 50//25
#endif
PRINT_CONFIG_VAR(OPTICFLOW_MAX_TRACK_CORNERS)

#ifndef OPTICFLOW_WINDOW_SIZE
#define OPTICFLOW_WINDOW_SIZE 10
#endif
PRINT_CONFIG_VAR(OPTICFLOW_WINDOW_SIZE)

#ifndef OPTICFLOW_SEARCH_DISTANCE
#define OPTICFLOW_SEARCH_DISTANCE 20
#endif
PRINT_CONFIG_VAR(OPTICFLOW_MAX_SEARCH_DISTANCE)

#ifndef OPTICFLOW_SUBPIXEL_FACTOR
#define OPTICFLOW_SUBPIXEL_FACTOR 10
#endif
PRINT_CONFIG_VAR(OPTICFLOW_SUBPIXEL_FACTOR)

#ifndef OPTICFLOW_MAX_ITERATIONS
#define OPTICFLOW_MAX_ITERATIONS 10
#endif
PRINT_CONFIG_VAR(OPTICFLOW_MAX_ITERATIONS)

#ifndef OPTICFLOW_THRESHOLD_VEC
#define OPTICFLOW_THRESHOLD_VEC 2
#endif
PRINT_CONFIG_VAR(OPTICFLOW_THRESHOLD_VEC)

#ifndef OPTICFLOW_PYRAMID_LEVEL
#define OPTICFLOW_PYRAMID_LEVEL 2//0
#endif
PRINT_CONFIG_VAR(OPTICFLOW_PYRAMID_LEVEL)

#ifndef OPTICFLOW_FAST9_ADAPTIVE
#define OPTICFLOW_FAST9_ADAPTIVE TRUE
#endif
PRINT_CONFIG_VAR(OPTICFLOW_FAST9_ADAPTIVE)

#ifndef OPTICFLOW_FAST9_THRESHOLD
#define OPTICFLOW_FAST9_THRESHOLD 20
#endif
PRINT_CONFIG_VAR(OPTICFLOW_FAST9_THRESHOLD)

#ifndef OPTICFLOW_FAST9_MIN_DISTANCE
#define OPTICFLOW_FAST9_MIN_DISTANCE 10
#endif
PRINT_CONFIG_VAR(OPTICFLOW_FAST9_MIN_DISTANCE)

#ifndef OPTICFLOW_FAST9_PADDING
#define OPTICFLOW_FAST9_PADDING 20
#endif
PRINT_CONFIG_VAR(OPTICFLOW_FAST9_PADDING)

// thresholds FAST9 that are currently not set from the GCS:
#define FAST9_LOW_THRESHOLD 5
#define FAST9_HIGH_THRESHOLD 60

#ifndef OPTICFLOW_METHOD
#define OPTICFLOW_METHOD 0
#endif
PRINT_CONFIG_VAR(OPTICFLOW_METHOD)

#if OPTICFLOW_METHOD > 1
#error WARNING: Both Lukas Kanade and EdgeFlow are NOT selected
#endif

#ifndef OPTICFLOW_DEROTATION
#define OPTICFLOW_DEROTATION TRUE
#endif
PRINT_CONFIG_VAR(OPTICFLOW_DEROTATION)

#ifndef OPTICFLOW_MEDIAN_FILTER
#define OPTICFLOW_MEDIAN_FILTER TRUE
#endif
PRINT_CONFIG_VAR(OPTICFLOW_MEDIAN_FILTER)

//include state info for comparrison
#include "state.h"
#include "subsystems/datalink/telemetry.h"

//Include median filter
#include "filters/median_filter.h"
struct MedianFilterInt vel_x_filt, vel_y_filt;

//for kalman
#include "subsystems/imu.h"
#include "modules/state_autonomous_race/state_autonomous_race.h"

//for stereo optic flow
#include "modules/stereocam/stereocam2state/stereocam2state.h"


/* Functions only used here */
static uint32_t timeval_diff(struct timeval *starttime, struct timeval *finishtime);
static int cmp_flow(const void *a, const void *b);

//FOR FILE LOGGING
float vel_x = 0;
float vel_y = 0;
float body_vel_x = 0;
float body_vel_y = 0;
float opt_body_v_x = 0;
float opt_body_v_y = 0;
float optitrack_vel_x = 0;
float optitrack_vel_y = 0;   
float psi_ = 0;

//for derotation
struct FloatEulers prev_att = {0};

float vec_sum_x = 0;
float vec_sum_y = 0;
float vec_sum_z = 0;

struct timeval stop, start;
int16_t frame_counter = 0;


static void opticflow_debug_send(struct transport_tx *trans, struct link_device *dev)
    {
    pprz_msg_send_OPTIC_FLOW_DEBUG(trans, dev, AC_ID,&optitrack_vel_x,&optitrack_vel_y,&body_vel_x,&body_vel_y,&vel_x,&vel_y,&psi_);//
    }  


/**
 * Initialize the opticflow calculator
 * @param[out] *opticflow The new optical flow calculator
 * @param[in] *w The image width
 * @param[in] *h The image height
 */
void flow_speed_calc_init(struct opticflow_t *opticflow, uint16_t w, uint16_t h)
{

  register_periodic_telemetry(DefaultPeriodic, PPRZ_MSG_ID_OPTIC_FLOW_DEBUG, opticflow_debug_send);
  init_median_filter(&vel_x_filt);
  init_median_filter(&vel_y_filt);

  /* Create the image buffers */
  image_create(&opticflow->img_gray, w, h, IMAGE_GRAYSCALE);
  image_create(&opticflow->prev_img_gray, w, h, IMAGE_GRAYSCALE);

  /* Set the previous values */
  opticflow->got_first_img = false;
  FLOAT_RATES_ZERO(opticflow->prev_rates);

  /* Set the default values */
  opticflow->method = OPTICFLOW_METHOD; //0 = LK_fast9, 1 = Edgeflow
  opticflow->window_size = OPTICFLOW_WINDOW_SIZE;
  opticflow->search_distance = OPTICFLOW_SEARCH_DISTANCE;
  opticflow->derotation = OPTICFLOW_DEROTATION; //0 = OFF, 1 = ON

  opticflow->max_track_corners = OPTICFLOW_MAX_TRACK_CORNERS;
  opticflow->subpixel_factor = OPTICFLOW_SUBPIXEL_FACTOR;
  opticflow->max_iterations = OPTICFLOW_MAX_ITERATIONS;
  opticflow->threshold_vec = OPTICFLOW_THRESHOLD_VEC;
  opticflow->pyramid_level = OPTICFLOW_PYRAMID_LEVEL;
  opticflow->median_filter = OPTICFLOW_MEDIAN_FILTER;


  opticflow->fast9_adaptive = OPTICFLOW_FAST9_ADAPTIVE;
  opticflow->fast9_threshold = OPTICFLOW_FAST9_THRESHOLD;
  opticflow->fast9_min_distance = OPTICFLOW_FAST9_MIN_DISTANCE;
  opticflow->fast9_padding = OPTICFLOW_FAST9_PADDING;
  opticflow->fast9_rsize = 512;
  opticflow->fast9_ret_corners = malloc(sizeof(struct point_t) * opticflow->fast9_rsize);
  
  //temporary for timing
  gettimeofday(&start, 0);

}
/**
 * Run the optical flow with fast9 and lukaskanade on a new image frame
 * @param[in] *opticflow The opticalflow structure that keeps track of previous images
 * @param[in] *state The state of the drone
 * @param[in] *img The image frame to calculate the optical flow from
 * @param[out] *result The optical flow result
 */
void calc_fast9_lukas_kanade(struct opticflow_t *opticflow, struct opticflow_state_t *state, struct FloatEulers *att_eulers, struct image_t *img,
                             struct opticflow_result_t *result)
{
  if (opticflow->just_switched_method) {
    flow_speed_calc_init(opticflow, img->w, img->h);
  }

  // variables for size_divergence:
  float size_divergence; int n_samples;

  // variables for linear flow fit:
  float error_threshold;
  int n_iterations_RANSAC, n_samples_RANSAC, success_fit;
  struct linear_flow_fit_info fit_info;
  
  gettimeofday(&stop, 0);
  double curr_time = (double)(stop.tv_sec + stop.tv_usec / 1000000.0);
  double elapsed = curr_time - (double)(start.tv_sec + start.tv_usec / 1000000.0);
  gettimeofday(&start, 0);

  // Update FPS for information
  frame_counter += 1;
  img->ts.tv_sec + img->ts.tv_usec / 1000000.0;
  result->fps =  1 /(timeval_diff(&opticflow->prev_timestamp, &img->ts) / 1000.);
  opticflow->prev_timestamp = img->ts;

  // Convert image to grayscale
  image_to_grayscale(img, &opticflow->img_gray);

  // Copy to previous image if not set
  if (!opticflow->got_first_img) {
    image_copy(&opticflow->img_gray, &opticflow->prev_img_gray);
    opticflow->got_first_img = true;
  }

  // *************************************************************************************
  // Corner detection
  // *************************************************************************************

  // FAST corner detection
  // TODO: There is something wrong with fast9_detect destabilizing FPS. This problem is reduced with putting min_distance
  // to 0 (see defines), however a more permanent solution should be considered
  fast9_detect(img, opticflow->fast9_threshold, opticflow->fast9_min_distance,
               opticflow->fast9_padding, opticflow->fast9_padding, &result->corner_cnt,
               &opticflow->fast9_rsize,
               opticflow->fast9_ret_corners);

  // Adaptive threshold
  if (opticflow->fast9_adaptive) {
    // Decrease and increase the threshold based on previous values
    if (result->corner_cnt < 40
        && opticflow->fast9_threshold > FAST9_LOW_THRESHOLD) { // TODO: Replace 40 with OPTICFLOW_MAX_TRACK_CORNERS / 2
      opticflow->fast9_threshold--;
    } else if (result->corner_cnt > OPTICFLOW_MAX_TRACK_CORNERS * 2 && opticflow->fast9_threshold < FAST9_HIGH_THRESHOLD) {
      opticflow->fast9_threshold++;
    }
  }

#if OPTICFLOW_SHOW_CORNERS
  image_show_points(img, opticflow->fast9_ret_corners, result->corner_cnt);
#endif

  // Check if we found some corners to track
  if (result->corner_cnt < 1) {
    image_copy(&opticflow->img_gray, &opticflow->prev_img_gray);
    return;
  }

  // *************************************************************************************
  // Corner Tracking
  // *************************************************************************************

  // Execute a Lucas Kanade optical flow
  result->tracked_cnt = result->corner_cnt;
  struct flow_t *vectors = opticFlowLK(&opticflow->img_gray, &opticflow->prev_img_gray, opticflow->fast9_ret_corners,
                                       &result->tracked_cnt,
                                       opticflow->window_size / 2, opticflow->subpixel_factor, opticflow->max_iterations,
                                       opticflow->threshold_vec, opticflow->max_track_corners, opticflow->pyramid_level);

  
  struct FloatVect3 cross_vectors[result->tracked_cnt];
  struct FloatVect3 vel_vec_temp, vel_sum, v_est_vec;
  
  vec_sum_x = 0;
  vec_sum_y = 0;
  vec_sum_z = 0;
  
  derotate_flow_vectors(vectors, att_eulers,cross_vectors, result->tracked_cnt, opticflow->subpixel_factor);
  
  v_est_vec.x = vec_sum_x;
  v_est_vec.y = vec_sum_y;
  v_est_vec.z = vec_sum_z;
  
  vel_sum.x = 0;
  vel_sum.y = 0;
  vel_sum.z = 0;
  
  for(int i = 1; i < result->tracked_cnt;i++){
    VECT3_CROSS_PRODUCT(vel_vec_temp,cross_vectors[i],cross_vectors[i-1]);
    //sign = dot(foe,flow_sum)/(norm(flow_sum)*norm(foe));
    
    float sign = VECT3_DOT_PRODUCT(vel_vec_temp, v_est_vec);
      if(sign>0){
	VECT3_SMUL(vel_vec_temp, vel_vec_temp, -1) 
      }
    vel_sum.x += vel_vec_temp.x;
    vel_sum.y += vel_vec_temp.y;
    vel_sum.z += vel_vec_temp.z;
//     vec_sum_x += vel_vec_temp.x;
//     vec_sum_y += vel_vec_temp.y;
//     vec_sum_z += vel_vec_temp.z;
  }

  double norm = sqrt(VECT3_NORM2(vel_sum));
  VECT3_SDIV(vel_sum, vel_sum, norm);
  
  printf("vel_sum_x:%f\n",vel_sum.x);
  printf("vel_sum_y:%f\n",vel_sum.y);
  printf("vel_sum_z:%f\n",vel_sum.z);
  
    printf("vec_sum_x:%f\n",vec_sum_x);
    printf("vec_sum_y:%f\n",vec_sum_y);
     printf("vec_sum_z:%f\n",vec_sum_z);

  
#if OPTICFLOW_SHOW_FLOW
  //printf("show: n tracked = %d\n", result->tracked_cnt);
  image_show_flow(img, vectors, result->tracked_cnt, opticflow->subpixel_factor);
#endif

  // Estimate size divergence:
  if (SIZE_DIV) {
    n_samples = 100;
    size_divergence = get_size_divergence(vectors, result->tracked_cnt, n_samples);
    result->div_size = size_divergence;
  } else {
    result->div_size = 0.0f;
  }
  if (LINEAR_FIT) {
    // Linear flow fit (normally derotation should be performed before):
    error_threshold = 10.0f;
    n_iterations_RANSAC = 20;
    n_samples_RANSAC = 5;
    success_fit = analyze_linear_flow_field(vectors, result->tracked_cnt, error_threshold, n_iterations_RANSAC,
                                            n_samples_RANSAC, img->w, img->h, &fit_info);

    if (!success_fit) {
      fit_info.divergence = 0.0f;
      fit_info.surface_roughness = 0.0f;
    }

    result->divergence = fit_info.divergence;
    result->surface_roughness = fit_info.surface_roughness;
  } else {
    result->divergence = 0.0f;
    result->surface_roughness = 0.0f;
  }


  // Get the median flow
  qsort(vectors, result->tracked_cnt, sizeof(struct flow_t), cmp_flow);
  if (result->tracked_cnt == 0) {
    // We got no flow
    result->flow_x = 0;
    result->flow_y = 0;
  } else if (result->tracked_cnt > 3) {
    // Take the average of the 3 median points
    result->flow_x = vectors[result->tracked_cnt / 2 - 1].flow_x;
    result->flow_y = vectors[result->tracked_cnt / 2 - 1].flow_y;
    result->flow_x += vectors[result->tracked_cnt / 2].flow_x;
    result->flow_y += vectors[result->tracked_cnt / 2].flow_y;
    result->flow_x += vectors[result->tracked_cnt / 2 + 1].flow_x;
    result->flow_y += vectors[result->tracked_cnt / 2 + 1].flow_y;
    result->flow_x /= 3;
    result->flow_y /= 3;
  } else {
    // Take the median point
    result->flow_x = vectors[result->tracked_cnt / 2].flow_x;
    result->flow_y = vectors[result->tracked_cnt / 2].flow_y;
  }

  // Flow Derotation
  float diff_flow_x = 0;
  float diff_flow_y = 0;

  /*// Flow Derotation TODO:
  float diff_flow_x = (state->phi - opticflow->prev_phi) * img->w / OPTICFLOW_FOV_W;
  float diff_flow_y = (state->theta - opticflow->prev_theta) * img->h / OPTICFLOW_FOV_H;*/

  if (opticflow->derotation && result->tracked_cnt > 5) {
    diff_flow_x = (state->rates.p + opticflow->prev_rates.p) / 2.0f / result->fps * img->w /
                  OPTICFLOW_FOV_W;// * img->w / OPTICFLOW_FOV_W;
    diff_flow_y = (state->rates.q + opticflow->prev_rates.q) / 2.0f / result->fps * img->h /
                  OPTICFLOW_FOV_H;// * img->h / OPTICFLOW_FOV_H;
  }

  float correction_factor_x = 0.8;
  float correction_factor_y = 0.85;
  //diff_flow_y * opticflow->subpixel_factor * correction_factor_y;//
  //result->flow_y;
  result->flow_der_x = result->flow_x - diff_flow_x * opticflow->subpixel_factor * correction_factor_x;
  result->flow_der_y = result->flow_y - diff_flow_y * opticflow->subpixel_factor * correction_factor_y;
  opticflow->prev_rates = state->rates;

  // Velocity calculation
  // Right now this formula is under assumption that the flow only exist in the center axis of the camera.
  // TODO Calculate the velocity more sophisticated, taking into account the drone's angle and the slope of the ground plane.
   vel_x = result->flow_der_x * result->fps * state->agl / opticflow->subpixel_factor  / OPTICFLOW_FX;
   vel_y = result->flow_der_y * result->fps * state->agl / opticflow->subpixel_factor  / OPTICFLOW_FY;
   
   //Compare with optitrack speeds (used in file logger)
    optitrack_vel_x = stateGetSpeedNed_f()->x;
    optitrack_vel_y = stateGetSpeedNed_f()->y;
   
     psi_ = stateGetNedToBodyEulers_f()->psi;
    float s_psi = sinf(-psi_);
    float c_psi = cosf(-psi_);
    body_vel_y = (c_psi * optitrack_vel_x - s_psi * optitrack_vel_y);
    body_vel_x = -c_psi * optitrack_vel_y + s_psi * optitrack_vel_x; 

  //Apply a  median filter to the velocity if wanted
  if (opticflow->median_filter == true) {
    result->vel_x = (float)update_median_filter(&vel_x_filt, (int32_t)(vel_x * 1000)) / 1000;
    result->vel_y = (float)update_median_filter(&vel_y_filt, (int32_t)(vel_y * 1000)) / 1000;
  } else {
    result->vel_x = vel_x;
    result->vel_y = vel_y;
  }
  // Velocity calculation: uncomment if focal length of the camera is not known or incorrect.
  //  result->vel_x =  - result->flow_der_x * result->fps * state->agl / opticflow->subpixel_factor * OPTICFLOW_FOV_W / img->w
  //  result->vel_y =  result->flow_der_y * result->fps * state->agl / opticflow->subpixel_factor * OPTICFLOW_FOV_H / img->h


  // Determine quality of noise measurement for state filter
  //TODO develop a noise model based on groundtruth

  float noise_measurement_temp = (1 - ((float)result->tracked_cnt / ((float)opticflow->max_track_corners * 1.25)));
  result->noise_measurement = noise_measurement_temp;
  
  //DEBUFGGING
  //result->tracked_cnt = frame_counter;

  // *************************************************************************************
  // Next Loop Preparation
  // *************************************************************************************
  free(vectors);
  image_switch(&opticflow->img_gray, &opticflow->prev_img_gray);
}

void undistort_fisheye_point(float point_x, float point_y, float *undistorted_x, float *undistorted_y, int f, float k, float x_img_center, float y_img_center)
{
  /*
f = 168;
k = 1.085;
%k = 1.051;
%k = 1.118;

*/
  
  float x_mid = point_x - 157.0f;//-(float)(x_princip);
  float y_mid = point_y - 12.5;//32.0f;//-(float)(y_princip);
  
  //to polar coordinates
  float r = sqrtf((pow(x_mid,2))+(pow(y_mid,2)));
  float theta = atan2f(y_mid,x_mid);//atanf?
  
  //k = 1.085;
//k = 1.051;

//k = 1.080;//last k
   k = 1.118;
  //k = 1.218;
  
  //radial distortion correction
  float R = (float)f*tan(asin(sin( atan(r/(float)f))*k));
  
                                                  // +y
                                                  // ^
                                                  // |
                                                  // |
  *undistorted_x =  R * cos(theta);//+x_princip; in (0,0)--->+x
  *undistorted_y =  R * sin(theta);//+y_princip;
  
  
}

void vec_from_point_2(float point_x, float point_y, int f, struct FloatVect3 *vec)
{
  f = 168;
  float focus = (float)f;
  vec->x = 1.0;
  vec->y = point_x/focus;
  vec->z = point_y/focus;
}

void proj_vec(float *point_x, float *point_y, int f, struct FloatVect3 vec)
{
//   f = 168;
//   float focus = (float)f;
//   vec->x = 1.0;
//   vec->y = point_x/focus;
//   vec->z = point_y/focus;
  f = 168;
  float focus = (float)f;
  *point_x = vec.y*focus;
  *point_y = vec.z*focus;
}

void back_proj_points(struct FloatVect3 *gate_point, struct FloatVect3 *cam_pos, struct FloatMat33 *R_mat, float *x_res, float *y_res)
{
  
//   point_3d = R*(gate_point-cam_pos');
// 
//     hom_coord = [point_3d(2)/point_3d(1);
//                  point_3d(3)/point_3d(1);
//                  1];
// 
//     res = intr*hom_coord;
  
  struct FloatVect3 temp1, point_3d, hom_coord, res;
  struct FloatMat33 intr, R_t;
  
  VECT3_DIFF(temp1,*gate_point,*cam_pos);
  MAT33_TRANS(R_t,*R_mat);
  MAT33_VECT3_MUL(point_3d,R_t,temp1);
//   MAT33_VECT3_MUL(point_3d,*R_mat,temp1);
  hom_coord.x = point_3d.y/point_3d.x;
  hom_coord.y = -point_3d.z/point_3d.x;
  hom_coord.z = 1;
  
//   debug_1 = hom_coord.x;
//   debug_2 = hom_coord.y;
  
  MAT33_ELMT(intr, 0, 0) = 168;//(row,column)
  MAT33_ELMT(intr, 0, 1) = 0;
  MAT33_ELMT(intr, 0, 2) = 157.2;

  MAT33_ELMT(intr, 1, 0) = 0;//(row,column)
  MAT33_ELMT(intr, 1, 1) = 168;
  MAT33_ELMT(intr, 1, 2) = 12.5;

  MAT33_ELMT(intr, 2, 0) = 0;//(row,column)
  MAT33_ELMT(intr, 2, 1) = 0;
  MAT33_ELMT(intr, 2, 2) = 0;
  
  MAT33_VECT3_MUL(res,intr,hom_coord);
  
//   debug_1 = res.x;
//   debug_2 = res.y;
  
  *x_res = res.x;
  *y_res = res.y;
  
  
  //Rmat
  
  //image_error()//compare with undist_x and undist_y
}


void derotate_flow_vectors(struct flow_t *vectors, struct FloatEulers *att_eulers, struct FloatVect3 *cross_vecs, uint16_t points_cnt, uint8_t subpixel_factor)
{
    struct FloatVect3 temp_vec_1, temp_vec_2, temp_vec_3, rot_vec, rot_vec_1, rot_vec_2, der_flow,
    vec_curr, vec_prev, test_vec_1, test_vec_2;
    struct FloatRMat R, R_temp, R_curr, R_prev, R_20;
    struct FloatEulers attitude_delta, curr_att, cam_body;
    
  //get cur attitude and heading 
//     float curr_att_phi = stateGetNedToBodyEulers_f()->phi;
//     float curr_att_theta = stateGetNedToBodyEulers_f()->theta+20*(3.14/180);
//     float curr_att_psi = stateGetNedToBodyEulers_f()->psi;
    curr_att.phi = att_eulers->phi;
    curr_att.theta = att_eulers->theta;
    curr_att.psi = -att_eulers->psi;
    
    attitude_delta.phi = (prev_att.phi-curr_att.phi);
    attitude_delta.theta = prev_att.theta-curr_att.theta;
    attitude_delta.psi = prev_att.psi-curr_att.psi;
    
    printf("attitude_delta.phi:%f\n",attitude_delta.phi*(180/3.14));
    printf("attitude_delta.theta:%f\n",attitude_delta.theta*(180/3.14));
    printf("attitude_delta.psi:%f\n",attitude_delta.psi*(180/3.14));
    
  // Go through all the points
  for (uint16_t i = 0; i < points_cnt; i++) {
    // Draw a line from the original position with the flow vector
    struct point_t from = {
      vectors[i].pos.y / subpixel_factor,
      vectors[i].pos.x / subpixel_factor
    };
    struct point_t to = {
      (vectors[i].pos.y + vectors[i].flow_y) / subpixel_factor,
      (vectors[i].pos.x + vectors[i].flow_x) / subpixel_factor
    };
    
    float undist_point_x;
    float undist_point_y;
    float undist_flow_x;
    float undist_flow_y;
    
    float dist_point_x  = vectors[i].pos.y/subpixel_factor;
    float dist_point_y  = vectors[i].pos.x/subpixel_factor;
    float dist_flow_x = (vectors[i].pos.y + vectors[i].flow_y) / subpixel_factor;
    float dist_flow_y = (vectors[i].pos.x + vectors[i].flow_x) / subpixel_factor;
    
    float x_derotated;
    float y_derotated;
    
    float x_derotated_1;
    float y_derotated_1;
    float x_derotated_2;
    float y_derotated_2;
    
    //printf("dist_point_x:%f\n",dist_point_x);
    //undistort points and flow 
    undistort_fisheye_point(dist_point_x,dist_point_y, &undist_point_x, &undist_point_y, 168, 1.118, 157, 12.5);
    undistort_fisheye_point(dist_flow_x,dist_flow_y, &undist_flow_x, &undist_flow_y, 168, 1.118, 157, 12.5);
    
    //vector axis system x forward y left z up
    //vector from point
    vec_from_point_2(undist_point_x,undist_point_y, 168,&temp_vec_1);
//     double norm = sqrt(VECT3_NORM2(temp_vec_1));
//     VECT3_SDIV(temp_vec_1, temp_vec_1, norm);
    
    //vector from flow point 
    vec_from_point_2(undist_flow_x,undist_flow_y, 168,&temp_vec_2);
//     norm = sqrt(VECT3_NORM2(temp_vec_2));
//     VECT3_SDIV(temp_vec_2, temp_vec_2, norm);
    
    //camera to body rotation
    cam_body.phi = 0;
    cam_body.theta = -20*(3.14/180);
    cam_body.psi = 0;
    float_rmat_of_eulers_321(&R_20,&cam_body);
    MAT33_VECT3_MUL(vec_curr, R_20,temp_vec_1);
    MAT33_VECT3_MUL(vec_prev, R_20,temp_vec_2);
    
    //test rotations
    test_vec_1.x = 1;
    test_vec_1.y = 0.5;
    test_vec_1.z = 0;
    MAT33_VECT3_MUL(test_vec_2, R_20,test_vec_1);
    
    
    
//##################PREV and CUR ATT SWITCHED
    //rotation matrix of attitude and heading delta since last function call
    float_rmat_of_eulers_321(&R_curr,&prev_att);
    //MAT33_TRANS(R_curr,R_temp);
    
    float_rmat_of_eulers_321(&R_prev,&curr_att);
    //MAT33_TRANS(R_prev,R_temp);
    
    //rotate points to current attitude 
    MAT33_VECT3_MUL(rot_vec_1, R_curr,vec_curr);
//     norm = sqrt(VECT3_NORM2(rot_vec));
//     VECT3_SDIV(rot_vec, rot_vec, norm);
    MAT33_VECT3_MUL(rot_vec_2, R_prev,vec_prev);
    
    
    
    //VECT3_DIFF(der_flow,temp_vec_2,rot_vec);
    //VECT3_SUM(der_flow,temp_vec_2,rot_vec);
    VECT3_DIFF(der_flow,rot_vec_2,rot_vec_1);//mayby switch?
    
    VECT3_CROSS_PRODUCT(temp_vec_3,temp_vec_1,der_flow);
    
    cross_vecs[i] = temp_vec_3;
    
    //project rotated vector
    //back_proj_points(&rot_vec, struct FloatVect3 *cam_pos, struct FloatMat33 *R_mat, float *x_res, float *y_res)
    
    proj_vec(&x_derotated,&y_derotated, 168, rot_vec); 
    
    proj_vec(&x_derotated_1,&y_derotated_1, 168, rot_vec_1);
    proj_vec(&x_derotated_2,&y_derotated_2, 168, rot_vec_2);
    
    
    
//     vectors[i].pos.y = (undist_point_x+157.2)*subpixel_factor;
//     vectors[i].pos.x = (undist_point_y+12.5)*subpixel_factor;
//     //- undist_point_x /y
// //     vectors[i].flow_y = (undist_flow_x  - undist_point_x )*subpixel_factor;
// //     vectors[i].flow_x = (undist_flow_y  - undist_point_y )*subpixel_factor;
//     float flow_y = (undist_flow_x  - undist_point_x )*subpixel_factor;
//     float flow_x = (undist_flow_y  - undist_point_y )*subpixel_factor;
//     
//     
// //     vectors[i].flow_y = (x_derotated - undist_point_x)*subpixel_factor;
// //     vectors[i].flow_x = -(y_derotated - undist_point_y)*subpixel_factor;
//      float flow_y_rot = (x_derotated - undist_point_x)*subpixel_factor;
//      float flow_x_rot = (y_derotated - undist_point_y)*subpixel_factor;//was -
//      
// //      if(flow_y_rot < 0)flow_y_rot = 0;
// //      if(flow_x_rot < 0)flow_x_rot = 0;
//     
//     vectors[i].flow_y = flow_y + flow_y_rot;// (undist_flow_x  - x_derotated )*subpixel_factor;
//     vectors[i].flow_x = flow_x + flow_x_rot;//(undist_flow_y  + y_derotated )*subpixel_factor;
    //only plot vectors inside image
    if(x_derotated_1 < -157)x_derotated_1 = -156;
    if(y_derotated_1 < -12.5)y_derotated_1 = -12.5;
    vectors[i].pos.y = (x_derotated_1+157.2)*subpixel_factor;
    vectors[i].pos.x = (y_derotated_1+12.5)*subpixel_factor;
    vectors[i].flow_y = x_derotated_1 - x_derotated_2;
    vectors[i].flow_x = y_derotated_1 - y_derotated_2;
    
      vec_sum_x += der_flow.x;//vectors[i].flow_y;
      vec_sum_y += der_flow.y;//vectors[i].flow_x;
      vec_sum_z += der_flow.z;
     
    
    //printf("vectors[i].pos.y:%d\n",vectors[i].pos.y);
    
  }
  
    prev_att.phi = curr_att.phi;
    prev_att.theta = curr_att.theta;
    prev_att.psi = curr_att.psi;
    //printf("test_vec_2: x:%f y:%f z:%f\n",test_vec_2.x,test_vec_2.y,test_vec_2.z);
  
}

/**
 * Run the optical flow with EDGEFLOW on a new image frame
 * @param[in] *opticflow The opticalflow structure that keeps track of previous images
 * @param[in] *state The state of the drone
 * @param[in] *img The image frame to calculate the optical flow from
 * @param[out] *result The optical flow result
 */
void calc_edgeflow_tot(struct opticflow_t *opticflow, struct opticflow_state_t *state, struct image_t *img,
                       struct opticflow_result_t *result)
{
  // Define Static Variables
  static struct edge_hist_t edge_hist[MAX_HORIZON];
  static uint8_t current_frame_nr = 0;
  struct edge_flow_t edgeflow;
  static uint8_t previous_frame_offset[2] = {1, 1};

  // Define Normal variables
  struct edgeflow_displacement_t displacement;
  displacement.x = malloc(sizeof(int32_t) * img->w);
  displacement.y = malloc(sizeof(int32_t) * img->h);

  // If the methods just switched to this one, reintialize the
  // array of edge_hist structure.
  if (opticflow->just_switched_method == 1) {
    int i;
    for (i = 0; i < MAX_HORIZON; i++) {
      edge_hist[i].x = malloc(sizeof(int32_t) * img->w);
      edge_hist[i].y = malloc(sizeof(int32_t) * img->h);
      FLOAT_RATES_ZERO(edge_hist[i].rates);
    }
  }

  uint16_t disp_range;
  if (opticflow->search_distance < DISP_RANGE_MAX) {
    disp_range = opticflow->search_distance;
  } else {
    disp_range = DISP_RANGE_MAX;
  }

  uint16_t window_size;

  if (opticflow->window_size < MAX_WINDOW_SIZE) {
    window_size = opticflow->window_size;
  } else {
    window_size = MAX_WINDOW_SIZE;
  }

  uint16_t RES = opticflow->subpixel_factor;

  //......................Calculating EdgeFlow..................... //

  // Calculate current frame's edge histogram
  int32_t *edge_hist_x = edge_hist[current_frame_nr].x;
  int32_t *edge_hist_y = edge_hist[current_frame_nr].y;
  calculate_edge_histogram(img, edge_hist_x, 'x', 0);
  calculate_edge_histogram(img, edge_hist_y, 'y', 0);


  // Copy frame time and angles of image to calculated edge histogram
  edge_hist[current_frame_nr].frame_time = img->ts;
  edge_hist[current_frame_nr].rates = state->rates;

  // Calculate which previous edge_hist to compare with the current
  uint8_t previous_frame_nr[2];
  calc_previous_frame_nr(result, opticflow, current_frame_nr, previous_frame_offset, previous_frame_nr);

  //Select edge histogram from the previous frame nr
  int32_t *prev_edge_histogram_x = edge_hist[previous_frame_nr[0]].x;
  int32_t *prev_edge_histogram_y = edge_hist[previous_frame_nr[1]].y;
  
  
  float fps_x = 0;
  float fps_y = 0;
  float time_diff_x = (float)(timeval_diff(&edge_hist[previous_frame_nr[0]].frame_time, &img->ts)) / 1000.;
  float time_diff_y = (float)(timeval_diff(&edge_hist[previous_frame_nr[1]].frame_time, &img->ts)) / 1000.;
  fps_x = 1 / (time_diff_x);
  fps_y = 1 / (time_diff_y);

  result->fps = fps_x;

  //Calculate the corresponding derotation of the two frames
  int16_t der_shift_x = 0;
  int16_t der_shift_y = 0;

  if (opticflow->derotation) {
    der_shift_x = (int16_t)(((edge_hist[previous_frame_nr[0]].rates.p + edge_hist[current_frame_nr].rates.p) / 2.0f /
                            result->fps *
                            (float)img->w / (OPTICFLOW_FOV_W))*1);
    der_shift_y = (int16_t)(((edge_hist[previous_frame_nr[1]].rates.q + edge_hist[current_frame_nr].rates.q) / 2.0f /
                            result->fps *
                            (float)img->h / (OPTICFLOW_FOV_H))*1);
  }

  // Estimate pixel wise displacement of the edge histograms for x and y direction
  calculate_edge_displacement(edge_hist_x, prev_edge_histogram_x,
                              displacement.x, img->w,
                              window_size, disp_range,0);//  der_shift_x);
  calculate_edge_displacement(edge_hist_y, prev_edge_histogram_y,
                              displacement.y, img->h,
                              window_size, disp_range,0);// der_shift_y);

  // Fit a line on the pixel displacement to estimate
  // the global pixel flow and divergence (RES is resolution)
  line_fit(displacement.x, &edgeflow.div_x,
           &edgeflow.flow_x, img->w,
           window_size + disp_range, RES);
  line_fit(displacement.y, &edgeflow.div_y,
           &edgeflow.flow_y, img->h,
           window_size + disp_range, RES);

  /* Save Resulting flow in results
   * Warning: The flow detected here is different in sign
   * and size, therefore this will be multiplied with
   * the same subpixel factor and -1 to make it on par with
   * the LK algorithm of t opticalflow_calculator.c
   * */
  edgeflow.flow_x = -1 * edgeflow.flow_x;
  edgeflow.flow_y = -1 * edgeflow.flow_y;

  result->flow_x = ((int16_t)edgeflow.flow_x / previous_frame_offset[0])-(der_shift_x*opticflow->subpixel_factor);
  result->flow_y = ((int16_t)edgeflow.flow_y / previous_frame_offset[1])-(der_shift_y*opticflow->subpixel_factor);

  //Fill up the results optic flow to be on par with LK_fast9
  result->flow_der_x =  der_shift_x*opticflow->subpixel_factor;
  result->flow_der_y =  der_shift_y*opticflow->subpixel_factor;
  result->corner_cnt = getAmountPeaks(edge_hist_x, 500 , img->w);
  result->tracked_cnt = getAmountPeaks(edge_hist_x, 500 , img->w);
  result->divergence = (float)edgeflow.flow_x / RES;
  result->div_size = 0.0f;
  result->noise_measurement = 0.0f;
  result->surface_roughness = 0.0f;

  //......................Calculating VELOCITY ..................... //

  /*Estimate fps per direction
   * This is the fps with adaptive horizon for subpixel flow, which is not similar
   * to the loop speed of the algorithm. The faster the quadcopter flies
   * the higher it becomes
  */

  // Calculate velocity
  //float vel_x = edgeflow.flow_x * fps_x * state->agl * OPTICFLOW_FOV_W / (img->w * RES);
  //float vel_y = edgeflow.flow_y * fps_y * state->agl * OPTICFLOW_FOV_H / (img->h * RES);

  float vel_x = result->flow_x * fps_x * state->agl * OPTICFLOW_FOV_W / (img->w * RES);
  float vel_y = result->flow_y * fps_y * state->agl * OPTICFLOW_FOV_H / (img->h * RES);

  
  //Apply a  median filter to the velocity if wanted
  if (opticflow->median_filter == true) {
    result->vel_x = (float)update_median_filter(&vel_x_filt, (int32_t)(vel_x * 1000)) / 1000;
    result->vel_y = (float)update_median_filter(&vel_y_filt, (int32_t)(vel_y * 1000)) / 1000;
  } else {
    result->vel_x = vel_x;
    result->vel_y = vel_y;
  }

  result->noise_measurement = 0.2;



#if OPTICFLOW_SHOW_FLOW
  draw_edgeflow_img(img, edgeflow, prev_edge_histogram_x, edge_hist_x);
#endif
  // Increment and wrap current time frame
  current_frame_nr = (current_frame_nr + 1) % MAX_HORIZON;
}


/**
 * Run the optical flow on a new image frame
 * @param[in] *opticflow The opticalflow structure that keeps track of previous images
 * @param[in] *state The state of the drone
 * @param[in] *img The image frame to calculate the optical flow from
 * @param[out] *result The optical flow result
 */
void flow_speed_calc_frame(struct opticflow_t *opticflow, struct opticflow_state_t *state, struct FloatEulers *att_eulers,  struct image_t *img,
                          struct opticflow_result_t *result)
{

  // A switch counter that checks in the loop if the current method is similar,
  // to the previous (for reinitializing structs)
  static int8_t switch_counter = -1;
  if (switch_counter != opticflow->method) {
    opticflow->just_switched_method = true;
    switch_counter = opticflow->method;
  } else {
    opticflow->just_switched_method = false;
  }

  // Switch between methods (0 = fast9/lukas-kanade, 1 = EdgeFlow)
  if (opticflow->method == 0) {
    calc_fast9_lukas_kanade(opticflow, state, att_eulers, img, result);
  } else if (opticflow->method == 1) {
    calc_edgeflow_tot(opticflow, state, img, result);
  }
  
  //downward_facing optic flow
  result->vel_body_x = result->vel_y;
  result->vel_body_y = - result->vel_x;
  
  //stereo optic flow 
  //result->vel_body_x = -s_flow_vel_z;
  //result->vel_body_y = s_flow_vel_x;


  // KALMAN filter
  struct Int32Vect3 acc_meas_body;
  struct Int32RMat *body_to_imu_rmat = orientationGetRMat_i(&imu.body_to_imu);
  int32_rmat_transp_vmult(&acc_meas_body, body_to_imu_rmat, &imu.accel);

  static uint8_t wait_counter = 0;
  static float previous_state_x[2] = {0.0f, 0.0f};
  static float covariance_x[4] = {1.0f, 1.0f, 1.0f, 1.0f};
  float measurements_x[2];
  static float previous_state_y[2] = {0.0f, 0.0f};
  static float pitch[4] = {1.0f, 1.0f, 1.0f, 1.0f};
  static float covariance_y[4] = {1.0f, 1.0f, 1.0f, 1.0f};
  float measurements_y[2];
  float process_noise[2] = {0.01f, 0.01f};
  float measurement_noise[2] = {0.2,1.0};//result->noise_measurement, 1.0f};


  if (opticflow->just_switched_method == 1) {

    wait_counter = 0;
    previous_state_x[0] = 0.0f;
    previous_state_x[1] = 0.0f;
    covariance_x[0] = 1.0f;
    covariance_x[1] = 1.0f;
    covariance_x[2] = 1.0f;
    covariance_x[3] = 1.0f;

    previous_state_y[0] = 0.0f;
    previous_state_y[1] = 0.0f;
    covariance_y[0] = 1.0f;
    covariance_y[1] = 1.0f;
    covariance_y[2] = 1.0f;
    covariance_y[3] = 1.0f;

  }

  if (wait_counter > 100) {


    measurements_x[0] = result->vel_body_x;
    measurements_x[1] = ACCEL_FLOAT_OF_BFP(acc_meas_body.x) - 0.24;

    measurements_y[0] = result->vel_body_y;
    measurements_y[1] = ACCEL_FLOAT_OF_BFP(acc_meas_body.y);

    //printf("measurements %f\n",ACCEL_FLOAT_OF_BFP(acc_meas_body.y));
    
    //Adaptive measurement noise when turning flag is set from guidance module
    if(states_race.turning == 1)measurement_noise[1] = 5;//random large value so kalman filter wont use the distorted opticflow velocity

    kalman_filter(measurements_x, covariance_x,
                  previous_state_x, process_noise, measurement_noise, result->fps);
    kalman_filter(measurements_y, covariance_y,
                  previous_state_y, process_noise, measurement_noise, result->fps);


    result->vel_body_x = previous_state_x[0];
    result->vel_body_y =  previous_state_y[0];
   
    opt_body_v_x  = previous_state_x[0];
    opt_body_v_y  = previous_state_y[0];
    
    //TIMING HERE


  } else {
    wait_counter++;
  }


}

/**
 * Calculate the difference from start till finish
 * @param[in] *starttime The start time to calculate the difference from
 * @param[in] *finishtime The finish time to calculate the difference from
 * @return The difference in milliseconds
 */
static uint32_t timeval_diff(struct timeval *starttime, struct timeval *finishtime)
{
  uint32_t msec;
  msec = (finishtime->tv_sec - starttime->tv_sec) * 1000;
  msec += (finishtime->tv_usec - starttime->tv_usec) / 1000;
  return msec;
}

/**
 * Compare two flow vectors based on flow distance
 * Used for sorting.
 * @param[in] *a The first flow vector (should be vect flow_t)
 * @param[in] *b The second flow vector (should be vect flow_t)
 * @return Negative if b has more flow than a, 0 if the same and positive if a has more flow than b
 */
static int cmp_flow(const void *a, const void *b)
{
  const struct flow_t *a_p = (const struct flow_t *)a;
  const struct flow_t *b_p = (const struct flow_t *)b;
  return (a_p->flow_x * a_p->flow_x + a_p->flow_y * a_p->flow_y) - (b_p->flow_x * b_p->flow_x + b_p->flow_y *
         b_p->flow_y);
}

void kalman_filter(float *measurements, float *covariance, float *state
                   , float *process_noise, float *measurement_noise, float fps)
{
  // _______ Preparation kalman filter _______ \\

  // process model (linear)
  float _G[2][2];
  MAKE_MATRIX_PTR(G, _G, 2);
  G[0][0] = 1.0f;
  G[0][1] = 1.0f / fps;
  G[1][0] = 0.0f;
  G[1][1] = 1.0f;


  // transpose of G
  float _Gtrans[2][2];
  MAKE_MATRIX_PTR(Gtrans, _Gtrans, 2);
  float_mat_copy(Gtrans, G, 2, 2);
  float_mat_transpose(Gtrans, 2);

  // measurement model (linear)
  float _H[2][2];
  MAKE_MATRIX_PTR(H, _H, 2);
  H[0][0] = 1.0f;
  H[0][1] = 0.0f;
  H[1][0] = 0.0f;
  H[1][1] = 1.0f;

  //transpose of H
  float _Htrans[2][2];
  MAKE_MATRIX_PTR(Htrans, _Htrans, 2);
  float_mat_copy(Htrans, H, 2, 2);
  float_mat_transpose(Htrans, 2);

  //Previous state
  float _Xprev[1][2];
  MAKE_MATRIX_PTR(Xprev, _Xprev, 2);
  Xprev[0][0] = state[0];
  Xprev[1][0] = state[1]; //state[1];

  //Previous covariance
  float _Pprevious[2][2];
  MAKE_MATRIX_PTR(Pprevious, _Pprevious, 2);
  Pprevious[0][0] = covariance[0];
  Pprevious[0][1] = covariance[1];
  Pprevious[1][0] = covariance[2];
  Pprevious[1][1] = covariance[3];

  //measurements;
  float _Z[1][2];
  MAKE_MATRIX_PTR(Z, _Z, 2);
  Z[0][0] = measurements[0];
  Z[1][0] = measurements[1];

  //Process noise model
  float _Q[2][2];
  MAKE_MATRIX_PTR(Q, _Q, 2);
  Q[0][0] = process_noise[0];
  Q[0][1] = 0.0f;
  Q[1][0] = 0.0f;
  Q[1][1] =  process_noise[1];

  //measurement nosie model
  float _R[2][2];
  MAKE_MATRIX_PTR(R, _R, 2);
  R[0][0] = measurement_noise[0];
  R[0][1] = 0.0f;
  R[1][0] = 0.0f;
  R[1][1] = measurement_noise[1];

  //Variables during kalman computation:
  float _Xpredict[1][2];
  MAKE_MATRIX_PTR(Xpredict, _Xpredict, 2);
  float _Xnext[1][2];
  MAKE_MATRIX_PTR(Xnext, _Xnext, 2);

  float _Ppredict[2][2];
  MAKE_MATRIX_PTR(Ppredict, _Ppredict, 2);
  float _Pnext[2][2];
  MAKE_MATRIX_PTR(Pnext, _Pnext, 2);

  float _K[2][2];
  MAKE_MATRIX_PTR(K, _K, 2);

  float _eye[2][2];
  MAKE_MATRIX_PTR(eye, _eye, 2);
  eye[0][0] = 1.0f;
  eye[0][1] = 0.0f;
  eye[1][0] = 0.0f;
  eye[1][1] = 1.0f;

  float _temp_mat[2][2];
  MAKE_MATRIX_PTR(temp_mat, _temp_mat, 2);
  float _temp_mat2[2][2];
  MAKE_MATRIX_PTR(temp_mat2, _temp_mat2, 2)
  float _temp_mat3[2][2];
  MAKE_MATRIX_PTR(temp_mat3, _temp_mat3, 2)

  float _temp_vec[1][2];
  MAKE_MATRIX_PTR(temp_vec, _temp_vec, 2);
  float _temp_vec2[1][2];
  MAKE_MATRIX_PTR(temp_vec2, _temp_vec2, 2)


  //_______  KALMAN FILTER_______ \\
  //_______  calculate state predict  _______ \\

  //Xpredict = G* Xprev;
  float_mat_mul(Xpredict, G, Xprev, 2, 2, 1);

  // _______  calculate covariance predict _______ \\
  // Ppredict = G*Pprevious*Gtrans + Q
  //...Pprevious*Gtrans...
  float_mat_mul(temp_mat, Pprevious, Gtrans, 2, 2, 2);
  //G*Pprevious*Gtrans...
  float_mat_mul(temp_mat2, G, temp_mat, 2, 2, 2);
  //G*Pprevious*Gtrans+Q
  float_mat_sum(Ppredict, temp_mat2, Q, 2, 2);

  //_______ Calculate Kalman gain _______ \\
  // K = Ppredict * Htrans /( H * Ppredict * Htrans + R)
  // ... Ppredict * Htrans ...
  float_mat_mul(temp_mat, Ppredict, Htrans, 2, 2, 2);
  //... H * Predict * Htrans
  float_mat_mul(temp_mat2, H, temp_mat, 2, 2, 2);
  //..( H * Ppredict * Htrans + R)
  float_mat_sum(temp_mat3, temp_mat2, R, 2, 2);
  //...inv( H * Ppredict * Htrans + R)
  float det_temp2 = 1 / (temp_mat3[0][0] * temp_mat3[1][1] - temp_mat3[0][1] * temp_mat3[1][0]);
  temp_mat2[0][0] =  det_temp2 * (temp_mat3[1][1]);
  temp_mat2[0][1] =  det_temp2 * (-1 * temp_mat3[1][0]);
  temp_mat2[1][0] =  det_temp2 * (-1 * temp_mat3[0][1]);
  temp_mat2[1][1] =  det_temp2 * (temp_mat3[0][0]);
  // K = Ppredict * Htrans / *inv( H * Ppredict * Htrans + R)
  float_mat_mul(K, temp_mat, temp_mat2, 2, 2, 2);



  // _______  Calculate next state  _______ \\
  //Xnext = Xpredict + K *(Z - Htrans * Xpredict)
  // ... Htrans * Xpredict)
  float_mat_mul(temp_vec, Htrans, Xpredict, 2, 2, 1);

  //... (Z - Htrans * Xpredict)
  float_mat_diff(temp_vec2, Z, temp_vec, 2, 1);

  // ... K *(Z - Htrans * Xpredict)
  float_mat_mul(temp_vec, K, temp_vec2, 2, 2, 1);


  //Xnext = Xpredict + K *(Z - Htrans * Xpredict)
  float_mat_sum(Xnext, Xpredict, temp_vec, 2, 1);

  // _______ calculate next covariance matrix _______ \\
  // Pnext = (eye(2) - K*H)*P_predict
  // ...K*H...
  float_mat_mul(temp_mat, K, H, 2, 2, 2);
  //(eye(2) - K*H)
  float_mat_diff(temp_mat2, eye, temp_mat, 2, 2);
  // Pnext = (eye(2) - K*H)*P_predict
  float_mat_mul(Pnext, temp_mat2, Ppredict, 2, 2, 2);


  //save values for next state
  covariance[0] = Pnext[0][0];
  covariance[1] = Pnext[0][1];;
  covariance[2] = Pnext[1][0];;
  covariance[3] = Pnext[1][1];;


  state[0] = Xnext[0][0];
  state[1] = Xnext[1][0];
//
  /*
    printf("Xprev %f, %f\n",Xprev[0][0],Xprev[1][0]);
    printf("Xpredict %f, %f\n",Xpredict[0][0],Xpredict[1][0]);
    printf("Pprevious %f,%f,%f,%f\n", Pprevious[0][0],Pprevious[0][1],Pprevious[1][0],Pprevious[1][1]);
    printf("Ppredict %f,%f,%f,%f\n", Ppredict[0][0],Ppredict[0][1],Ppredict[1][0],Ppredict[1][1]);
    printf("Pnext %f,%f,%f,%f\n", Pnext[0][0],Pnext[0][1],Pnext[1][0],Pnext[1][1]);
    printf("Z %f, %f \n", Z[0][0],Z[1][0]);
    printf("Kalman %f,%f,%f,%f\n", K[0][0],K[0][1],K[1][0],K[1][1]);
    printf("Xnext %f, %f\n",Xnext[0][0],Xnext[1][0]);*/
}

