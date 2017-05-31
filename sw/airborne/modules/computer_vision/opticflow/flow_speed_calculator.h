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
 * @file modules/computer_vision/opticflow/opticflow_calculator.h
 * @brief Calculate velocity from optic flow.
 *
 * Using images from a vertical camera and IMU sensor data.
 */



#ifndef OPTICFLOW_CALCULATOR_H
#define OPTICFLOW_CALCULATOR_H

#include "std.h"
#include "inter_thread_data.h"
#include "lib/vision/image.h"
#include "lib/v4l/v4l2.h"
#include "state.h"

//FOR FILE LOGGING
float vel_x;
float vel_y;
extern float opt_body_v_x;
extern float opt_body_v_y;

struct opticflow_t {
  bool got_first_img;             ///< If we got a image to work with
  bool just_switched_method;
  struct FloatRates prev_rates;     ///< Gyro Rates from the previous image frame
  struct image_t img_gray;          ///< Current gray image frame
  struct image_t prev_img_gray;     ///< Previous gray image frame
  struct timeval prev_timestamp;    ///< Timestamp of the previous frame, used for FPS calculation

  uint8_t method;               ///< Method to use to calculate the optical flow
  uint16_t window_size;             ///< Window size for the blockmatching algorithm (general value for all methods)
  uint16_t search_distance;         ///< Search distance for blockmatching alg.
  bool derotation;             ///< Derotation switched on or off (depended on the quality of the gyroscope measurement)
  bool median_filter;          ///< Decides to use a median filter on the velocity

  uint16_t subpixel_factor;          ///< The amount of subpixels per pixel
  uint8_t max_iterations;           ///< The maximum amount of iterations the Lucas Kanade algorithm should do
  uint8_t threshold_vec;            ///< The threshold in x, y subpixels which the algorithm should stop
  uint8_t pyramid_level;        ///< Number of pyramid levels used in Lucas Kanade algorithm (0 == no pyramids used)

  uint8_t max_track_corners;        ///< Maximum amount of corners Lucas Kanade should track
  bool fast9_adaptive;            ///< Whether the FAST9 threshold should be adaptive
  uint8_t fast9_threshold;          ///< FAST9 corner detection threshold
  uint16_t fast9_min_distance;      ///< Minimum distance in pixels between corners
  uint16_t fast9_padding;           ///< Padding used in FAST9 detector

  uint16_t fast9_rsize;   ///< Amount of corners allocated
  struct point_t *fast9_ret_corners; ///< Corners
};


void flow_speed_calc_init(struct opticflow_t *opticflow, uint16_t w, uint16_t h);
void flow_speed_calc_frame(struct opticflow_t *opticflow, struct opticflow_state_t *state, struct FloatEulers *att_eulers,  struct image_t *img,
                          struct opticflow_result_t *result);

void calc_fast9_lukas_kanade(struct opticflow_t *opticflow, struct opticflow_state_t *state, struct FloatEulers *att_eulers, struct image_t *img,
                             struct opticflow_result_t *result);

void calc_edgeflow_tot(struct opticflow_t *opticflow, struct opticflow_state_t *state, struct image_t *img,
                       struct opticflow_result_t *result);

void kalman_filter(float *measurements, float *covariance, float *state
                   , float *process_noise, float *measurement_noise, float fps);

//Special vector vector operation
/* multiply _vin by _mat, store in _vout */
#define VECT3_VECT3_TRANS_MUL(_mat, _v_a,_v_b) {    \
    MAT33_ELMT((_mat),0,0) = (_v_a).x*(_v_b).x;  \
    MAT33_ELMT((_mat),0,1) = (_v_a).x*(_v_b).y;   \
    MAT33_ELMT((_mat),0,2) = (_v_a).x*(_v_b).z;   \
    MAT33_ELMT((_mat),1,0) = (_v_a).y*(_v_b).x;   \
    MAT33_ELMT((_mat),1,1) = (_v_a).y*(_v_b).y;   \
    MAT33_ELMT((_mat),1,2) = (_v_a).y*(_v_b).z;   \
    MAT33_ELMT((_mat),2,0) = (_v_a).z*(_v_b).x;   \
    MAT33_ELMT((_mat),2,1) = (_v_a).z*(_v_b).y;  \
    MAT33_ELMT((_mat),2,2) = (_v_a).z*(_v_b).z;   \
  }

  //custom linear algebra functions maybe later to paparazzi math
#define MAT33_TRANS(_mat1,_mat2) {     \
    MAT33_ELMT((_mat1),0,0) = MAT33_ELMT((_mat2),0,0);  \
    MAT33_ELMT((_mat1),0,1) = MAT33_ELMT((_mat2),1,0);  \
    MAT33_ELMT((_mat1),0,2) = MAT33_ELMT((_mat2),2,0);  \
    MAT33_ELMT((_mat1),1,0) = MAT33_ELMT((_mat2),0,1);  \
    MAT33_ELMT((_mat1),1,1) = MAT33_ELMT((_mat2),1,1);  \
    MAT33_ELMT((_mat1),1,2) = MAT33_ELMT((_mat2),2,1);  \
    MAT33_ELMT((_mat1),2,0) = MAT33_ELMT((_mat2),0,2);  \
    MAT33_ELMT((_mat1),2,1) = MAT33_ELMT((_mat2),1,2);  \
    MAT33_ELMT((_mat1),2,2) = MAT33_ELMT((_mat2),2,2);  \
  }

#endif /* OPTICFLOW_CALCULATOR_H */


