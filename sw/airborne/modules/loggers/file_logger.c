/*
 * Copyright (C) 2014 Freek van Tienen <freek.v.tienen@gmail.com>
 *
 * This file is part of paparazzi.
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
 * along with paparazzi; see the file COPYING.  If not, write to
 * the Free Software Foundation, 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 *
 */

/** @file modules/loggers/file_logger.c
 *  @brief File logger for Linux based autopilots
 */

#include "file_logger.h"

#include <stdio.h>
#include "std.h"

#include "subsystems/imu.h"
#include "firmwares/rotorcraft/stabilization.h"
#include "firmwares/rotorcraft/stabilization/stabilization_indi.h"
#include "firmwares/rotorcraft/guidance/guidance_indi.h"
#include "state.h"
#include "subsystems/actuators.h"

/** Set the default File logger path to the USB drive */
#ifndef FILE_LOGGER_PATH
#define FILE_LOGGER_PATH /data/video/usb
#endif

/** The file pointer */
static FILE *file_logger = NULL;

/** Start the file logger and open a new file */
void file_logger_start(void)
{
  uint32_t counter = 0;
  char filename[512];

  // Check for available files
  sprintf(filename, "%s/%05d.csv", STRINGIFY(FILE_LOGGER_PATH), counter);
  while ((file_logger = fopen(filename, "r"))) {
    fclose(file_logger);

    counter++;
    sprintf(filename, "%s/%05d.csv", STRINGIFY(FILE_LOGGER_PATH), counter);
  }

  file_logger = fopen(filename, "w");

  if (file_logger != NULL) {
    fprintf(
      file_logger,
      "counter,pos_NED_x, pos_NED_y, pos_NED_z, filt_accel_ned_x, filt_accel_ned_y, filt_accel_ned_z, quat_i, quat_x, quat_y, quat_z,  sp_quat_i, sp_quat_x, sp_quat_y, sp_quat_z, sp_accel_x, sp_accel_y, sp_accel_z, accel_ned_x, accel_ned_y, accel_ned_z, speed_ned_x, speed_ned_y, speed_ned_z, imu_accel_unscaled_x, imu_accel_unscaled_y, imu_accel_unscaled_z, body_rates_p, body_rates_q, body_rates_r, actuaros_pprz_0, actuaros_pprz_1, actuaros_pprz_2, actuaros_pprz_3\n"
    );
  }
}

/** Stop the logger an nicely close the file */
void file_logger_stop(void)
{
  if (file_logger != NULL) {
    fclose(file_logger);
    file_logger = NULL;
  }
}

/** Log the values to a csv file */
void file_logger_periodic(void)
{
  if (file_logger == NULL) {
    return;
  }
  static uint32_t counter;
  counter =3 ;
  struct Int32Quat * quat       = stateGetNedToBodyQuat_i();
  struct NedCoor_f * pos        = stateGetPositionNed_f();
  struct NedCoor_f * accel_ned  = stateGetAccelNed_f();
  struct NedCoor_f * speed_ned  = stateGetSpeedNed_f();
  struct FloatRates* rates_body = stateGetBodyRates_f();

  fprintf(file_logger, "%d,%f,%f,%f,%f,%f,%f,%d,%d,%d,%d,%d,%d,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%d,%d,%d,%f,%f,%f,%d,%d,%d,%d\n",
          counter,
          pos->x,
          pos->y,
          pos->z,
          filt_accel_ned[0].o[0],
          filt_accel_ned[1].o[0],
          filt_accel_ned[2].o[0],
          quat->qi,
          quat->qx,
          quat->qy,
          quat->qz,
          stab_att_sp_quat.qi,
          stab_att_sp_quat.qx,
          stab_att_sp_quat.qy,
          stab_att_sp_quat.qz,
          sp_accel.x,
          sp_accel.y,
          sp_accel.z,
          accel_ned->x,
          accel_ned->y,
          accel_ned->z,
          speed_ned->x,
          speed_ned->y,
          speed_ned->z,
          imu.accel_unscaled.x,
          imu.accel_unscaled.y,
          imu.accel_unscaled.z,
          rates_body->p,
          rates_body->q,
          rates_body->r,
          actuators_pprz[0],
          actuators_pprz[1],
          actuators_pprz[2],
          actuators_pprz[3]
         );
  counter++;
}
