/*
 * Copyright (C) Wilco Vlenterie
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
 * @file "modules/computer_vision/autoswarm//autoswarm.h"
 * @author Wilco Vlenterie
 * Autonomous bebop swarming module based on vision
 */

#ifndef AS_H
#define AS_H

#include "std.h"
#include "stdbool.h"
#include <math.h>                                           ///< Used for sin/cos/tan/M_PI/abs
#include <time.h>                                           ///< Used to write time to results.txt
#include <stdio.h>                                          ///< Used for printing
#include <errno.h>                                          ///< Used for error handling

// Global options definitions
#define AS_POINT             0
#define AS_BUCKET            1
#define AS_CIRCLE_CW         2
#define AS_CIRCLE_CC         3
/// Camera options definitions
#define AS_CAM_FORWARD       1
#define AS_CAM_GLOBAL        2

extern double   settings_as_circle_radius;
extern double   settings_as_separation;
extern int      settings_as_attractor;
extern double   settings_as_vmax;
extern double   settings_as_global_strength;
extern int      settings_as_heading_mode;
extern double   settings_as_e;
extern double   settings_as_eps;

// Initialize global attractor
struct originPoint { double cx; double cy; double cz;};
struct originPoint globalOrigin;
static inline void setGlobalOrigin  (double x, double y, double z){ globalOrigin.cx = x; globalOrigin.cy = y; globalOrigin.cz = z;};
static inline bool setGlobalMode    (int mode){ settings_as_attractor = mode; return false; };
static inline bool setSwarmMode     (int mode){ settings_as_heading_mode = mode; return false; };

void autoswarm_init                 ( void );
void autoswarm_run                  ( void );
bool amIhome                        ( void );

#endif
