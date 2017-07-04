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
 * @file "modules/computer_vision//cv_active_random_filter.h"
 * @author Wilco Vlenterie (wv-tud)
 * Active random sampling colour filter
 */

#ifndef CV_BEBOP_CAMERA_STABILIZATION_H
#define CV_BEBOP_CAMERA_STABILIZATION_H

extern double     viewR;
extern double     default_k;
extern double     firstOrder_comp;
extern double     secondOrder_comp;
extern float      angleOfView;

extern uint16_t             cropCol;                                    ///< Column from which the ISP is cropped relative to MIN
extern int16_t              fillHeight;                                 ///< Column from which the ISP is cropped relative to MIN

void              cv_cam_stab_init(void);
/** Fisheye + Perspective correction **/
void              point2pixel        (double x_out, double y_out, double *x_in, double *y_in);
void              pixel2point        (double x_in, double y_in, double *x_out, double *y_out);
/** Conversion from angles to output frame point **/
void              angles2point        (double xAngle, double yAngle, double *x_out, double * y_out);
void              point2angles        (double x_out, double y_out, double *xAngle, double *yAngle);

#endif

