/*
 * Copyright (C) w.vlenterie
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
 * @file "modules/computer_vision/cv_image_pose.h"
 * @author w.vlenterie
 * Gets euler angles and rates at time of image capture
 */

#ifndef CV_IMAGE_POSE_H
#define CV_IMAGE_POSE_H

#include "pose_history/pose_history.h"

extern struct pose_t cv_image_pose;

void cv_image_pose_init( void );
struct image_t* cv_image_pose_func(struct image_t* img);

#endif

