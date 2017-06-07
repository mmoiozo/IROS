/*
 * Copyright (C) 2014 Hann Woei Ho
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
 * @file modules/computer_vision/flow_speed.h
 * @brief optical-flow calculation for Parrot Drones
 *
 */

#ifndef OPTICFLOW_MODULE_H
#define OPTICFLOW_MODULE_H

// Include opticflow calculator
#include "opticflow/flow_speed_calculator.h"

// Needed for settings
extern struct opticflow_t opticflow;

// Module functions
extern void flow_speed_init(void);
extern void flow_speed_run(void);
extern void flow_speed_start(void);
extern void flow_speed_stop(void);

#endif /* OPTICFLOW_MODULE_H */
