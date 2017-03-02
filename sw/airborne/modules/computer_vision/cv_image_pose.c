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
 * @file "modules/computer_vision/cv_image_pose.c"
 * @author w.vlenterie
 * Gets euler angles and rates at time of image capture
 */

#include "modules/computer_vision/cv.h"
#include "modules/pose_history/pose_history.h"
#include "modules/computer_vision/cv_image_pose.h"

#include BOARD_CONFIG

#ifndef IMAGE_POSE_CAMERA
#define IMAGE_POSE_CAMERA front_camera
#endif
PRINT_CONFIG_VAR(IMAGE_POSE_CAMERA)

struct pose_t cv_image_pose;

void cv_image_pose_init( void ) {
    cv_add_to_device(&IMAGE_POSE_CAMERA, cv_image_pose_func);
}

struct image_t* cv_image_pose_func(struct image_t* img){
    cv_image_pose = get_rotation_at_timestamp(img->pprz_ts);
    return NULL;
}
