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
#include "modules/computer_vision/bebop_camera_stabilization.h"
#include "modules/computer_vision/cv_image_pose.h"

#include BOARD_CONFIG

#ifndef CAM_STAB_CAMERA
#define CAM_STAB_CAMERA front_camera
#endif
PRINT_CONFIG_VAR(CAM_STAB_CAMERA);

uint16_t frame_shift;

struct image_t* cv_cam_stab_func(struct image_t* img);
struct image_t* cv_cam_stab_func(struct image_t* img)
{
	bebop_camera_stabilization((char*) img->buf, (uint16_t) img->w, (uint16_t) img->h, &cv_image_pose.eulers,&frame_shift);
	return NULL;
}

void cv_cam_stab_init() {
    bebop_camera_stabilization_init();
    cv_add_to_device(&CAM_STAB_CAMERA, cv_cam_stab_func);
}
