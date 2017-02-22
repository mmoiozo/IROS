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
 * @file "modules/computer_vision/autoswarm/autoswarm.h"
 * @author Wilco Vlenterie
 * Autonomous bebop swarming module based on vision
 */

using namespace std;
#include "autoswarm_opencv.h"                       // Include own header file
#include <computer_vision/active_random_filter.h>   // Include active random filter header file
#include <cmath>                                    // Used for sin/cos/tan/M_PI/abs
#include <ctime>                                    // Used to write time to results.txt
#include <algorithm>                                // Used for min() and max()
#include <string>                                   // Used for strlen
#include <vector>                                   // Used for active random filter
#include <stdio.h>                                  // Used for printing

#define BOARD_CONFIG "boards/bebop.h"
#define RADIO_CONTROL_TYPE_H "radio_control/rc_datalink.h"

extern "C" {
    #include "state.h"                              // Used for accessing state variables
    #include "navigation.h"                         // Used for navigation functions
    #include <errno.h>                              // Used for error handling
    #include "subsystems/gps/gps_datalink.h"
    #include "generated/flight_plan.h"              // C header used for WP definitions (causes problems) TODO: fix this
}

#ifndef AUTOSWARM_GLOBAL_ATTRACTOR
#define AUTOSWARM_GLOBAL_ATTRACTOR AUTOSWARM_CIRCLE_CW
#endif

static void calcCamPosition(struct NedCoor_f *pos, double totV[3], double cPos[3], double gi[3]);
static void updateWaypoints(struct NedCoor_f *pos, double totV[3], double cPos[3]);
static void calcVelocityResponse(struct NedCoor_f *pos, double totV[3], double gi[3]);
static void calcGlobalVelocity(struct NedCoor_f *pos, double gi[3]);
static void calcLocalVelocity(struct NedCoor_f *pos, double li[3]);
static void limitYaw(struct NedCoor_f *pos, double cPos[3]);
static void limitNorm(double totV[3], double maxNorm);
static void calcDiffVelocity(double totV[3], double di[3]);
static void limitVelocityYaw(double totV[3]);
static void autoswarm_opencv_run_header(void);
static void autoswarm_opencv_run_trailer(void);
static void calcGlobalCoeff(struct NedCoor_f *pos, double gi[3], double* g_coeff);

// Debug options
#define AUTOSWARM_SHOW_WAYPOINT 1                   // Show the updated positions of the waypoints
#define AUTOSWARM_SHOW_MEM      0                   // Show the neighbours identified and their location
#define AUTOSWARM_WRITE_RESULTS 0                   // Write measurements to text file
#define AUTOSWARM_BENCHMARK     0                   // Print benchmark table

// Set up swarm parameters
int     AUTOSWARM_MODE          = 2;                // 0: follow (deprecated), 1: look in direction of flight, 2: look in direction of global component
double  AUTOSWARM_SEPERATION    = 1.0;              // m
double  AUTOSWARM_LATTICE_RATIO = 2.5;

double  AUTOSWARM_E             = 0.0299;           // Was 0.01x - 0.0005 at 12m/s OK (but close)
double  AUTOSWARM_EPS           = 0.03;           //
double  AUTOSWARM_LOGLO         = 2.0;

double  AUTOSWARM_GLOBAL        = 1.0;              // % of V_MAX
int     AUTOSWARM_FPS           = 18;               // Frames per second
double  AUTOSWARM_VMAX          = 2.0;              // m/s
double  AUTOSWARM_YAWRATEMAX    = 180;              // deg/s
double  AUTOSWARM_HOME          = 0.4;              // m

// Set up global attractor parameters
int     AUTOSWARM_ATTRACTOR     = AUTOSWARM_GLOBAL_ATTRACTOR;
double  AUTOSWARM_CIRCLE_R      = 2.0;
double  AUTOSWARM_DEADZONE      = 0.2;

// Initialize parameters to be assigned during runtime
extern uint8_t              nav_block;
extern memoryBlock          neighbourMem[AR_FILTER_MAX_OBJECTS];
extern uint8_t              neighbourMem_size;
static struct FloatEulers * eulerAngles;
static struct NedCoor_f *   pos;
static uint32_t             runCount        = 0;
static const char *         flight_blocks[] = FP_BLOCKS;
static double               prev_v_d[3]     = {0.0, 0.0, 0.0};

//Optional function declarations
#if AUTOSWARM_BENCHMARK
    static void initBenchmark(void);
    static void addBenchmark(const char * title);
    static void showBenchmark(void);
    vector<double>              benchmark_time;
    vector<double>              benchmark_avg_time;
    vector<const char*>         benchmark_title;
    clock_t                     bench_start, bench_end;
#endif

#if AUTOSWARM_WRITE_RESULTS
    static time_t               startTime;
    static time_t               currentTime;
    static int                  curT;
    FILE *                      pFile;
    char                        resultFile [50];
#endif

#define PRINT(string,...) fprintf(stderr, "[autoswarm->%s()] " string,__FUNCTION__ , ##__VA_ARGS__)
#define AUTOSWARM_VERBOSE FALSE

#if AUTOSWARM_VERBOSE
#define VERBOSE_PRINT PRINT
#else
#define VERBOSE_PRINT(...)
#endif

void autoswarm_opencv_run(){
    // Computer vision compatibility function used to call trackObjects and (optionally) parse modified data back to rtp stream as YUV
    if (strcmp("Swarm",flight_blocks[nav_block]) && strcmp("Swarm Home",flight_blocks[nav_block])){
        VERBOSE_PRINT("Not swarming yet.\n");
        return;
    }    // We aren't swarming yet, lets save some battery life and skip processing for now
    autoswarm_opencv_run_header();                  // Mainly code for printing and debugging
    runCount++;                                     // Update global run-counter
    eulerAngles = stateGetNedToBodyEulers_f();      // Get Euler angles
    pos         = stateGetPositionNed_f();          // Get your current position
    double totV[3]  = {0.0, 0.0, 0.0};
    double cPos[3]  = {0.0, 0.0, 0.0};
    double gi[3]    = {0.0, 0.0, 0.0};
    calcVelocityResponse(pos, totV, gi);            // Calculate the velocity response of the agent and store result in totV and gi
    calcCamPosition(pos, totV, cPos, gi);           // Calculate WP_CAM/heading based on pos,totV and gi and store in cPos
    updateWaypoints(pos, totV, cPos);               // Update waypoints based on velocity contribution
    autoswarm_opencv_run_trailer();                  // Mainly code for printing and debugging
    return;
}

void calcVelocityResponse(struct NedCoor_f *pos, double totV[3], double gi[3]){
    double li[3]    = {0.0, 0.0, 0.0};              // Initialize local contribution
    double di[3]    = {0.0, 0.0, 0.0};              // Initialize differential contribution
    calcLocalVelocity(pos, li);                     // Get the contribution due to the neighbours we see and have memorized
    calcGlobalVelocity(pos, gi);                    // Get the contribution due to the "attraction" towards global origin
    double g_coeff;
    calcGlobalCoeff(pos, gi, &g_coeff);
    totV[0]         = li[0] + g_coeff * gi[0];      // Average the X local contribution (#neighbours independent) and add the X global contribution
    totV[1]         = li[1] + g_coeff * gi[1];      // Do the same for Y
    limitNorm(totV, AUTOSWARM_VMAX);                // Check if ideal velocity exceeds AUTOSWARM_VMAX
    calcDiffVelocity(totV, di);                     // Calculate the differential component based on change in totV
    totV[0]         = totV[0] + di[0];              // Add differential x component
    totV[1]         = totV[1] + di[1];              // Add differential y component
    if (AUTOSWARM_MODE!=2){
        limitVelocityYaw(totV);   // Limit the velocity when relative angle gets larger
    }
    return;
}

void calcLocalVelocity(struct NedCoor_f *pos, double li[3]){
    double range, li_fac=0;                         // Initialize variables
    unsigned int r;                                 // Initialize r

    for (r=0; r < neighbourMem_size; r++){         // For all neighbours found
        range   = sqrt(pow((pos->x - neighbourMem[r].x_w), 2.0) + pow((pos->y - neighbourMem[r].y_w), 2.0));
        li_fac  = 12 * AUTOSWARM_E / range * (pow(AUTOSWARM_SEPERATION / range, 12.0) - pow(AUTOSWARM_SEPERATION / range, 6.0)); // Local contribution scaling factor
        li[0]   += li_fac * (pos->x - neighbourMem[r].x_w) / range; // Local contribution in x
        li[1]   += li_fac * (pos->y - neighbourMem[r].y_w) / range; // Local contribution in y
    }
    if (r > 0){                                     // Only if we found anyone
        li[0] = li[0] / r;                          // Average the X local contribution (#neighbours independent) and add the X global contribution
        li[1] = li[1] / r;                          // Do the same for Y
    }
    else{
                                                    // TODO: What if I don't see anyone?
    }
    return;
}

void calcGlobalCoeff(struct NedCoor_f *pos, double gi[3], double* g_coeff){
    if (neighbourMem_size > 0){
        double nDir[2], rel_dist, range, c_gi_coef, c_gi[2];
        double w_neighbours[2]  = {0.0, 0.0};
        double gi_n             = sqrt(pow(gi[0], 2.0) + pow(gi[1], 2.0));

        for (unsigned int r = 0; r < neighbourMem_size; r++){
            range            = sqrt(pow((pos->x - neighbourMem[r].x_w), 2.0) + pow((pos->y - neighbourMem[r].y_w), 2.0));
            if(range > AUTOSWARM_SEPERATION * AUTOSWARM_LATTICE_RATIO){
                continue;
            }
            nDir[0]          = (pos->x - neighbourMem[r].x_w) / range;
            nDir[1]          = (pos->y - neighbourMem[r].y_w) / range;
            rel_dist         = min(1.0, ( AUTOSWARM_LATTICE_RATIO * AUTOSWARM_SEPERATION - range) / (AUTOSWARM_SEPERATION * (AUTOSWARM_LATTICE_RATIO - 1)));
            w_neighbours[0] += pow(rel_dist, 2.0) * nDir[0];
            w_neighbours[1] += pow(rel_dist, 2.0) * nDir[1];
        }
        w_neighbours[0] *= gi_n;
        w_neighbours[1] *= gi_n;
        c_gi_coef   = (w_neighbours[0] * gi[0] + w_neighbours[1] * gi[1]) / pow(gi_n, 2.0); // dot(w_neighbours, gi) / dot(gi, gi);
        c_gi[0]     = c_gi_coef * gi[0];
        c_gi[1]     = c_gi_coef * gi[1];
        *g_coeff    = pow(min(1.0, max(0.0, 2.0 - sqrt(pow(gi[0] - c_gi[0], 2.0) + pow(gi[1] - c_gi[1], 2.0)) / gi_n )), AUTOSWARM_LOGLO);
    }else{
        *g_coeff = 1.0;
    }
}

void calcGlobalVelocity(struct NedCoor_f *pos, double gi[3]){
    switch (AUTOSWARM_ATTRACTOR){
    case AUTOSWARM_POINT :{
       double cr     = sqrt(pow(globalOrigin.cx - pos->x, 2.0) + pow(globalOrigin.cy - pos->y, 2.0));
        if (cr > AUTOSWARM_DEADZONE){
            gi[0]               = AUTOSWARM_GLOBAL * AUTOSWARM_VMAX * (globalOrigin.cx - pos->x) / cr;
            gi[1]               = AUTOSWARM_GLOBAL * AUTOSWARM_VMAX * (globalOrigin.cy - pos->y) / cr;
        }
        break;
    }
    case AUTOSWARM_BUCKET :{
        double cr     = sqrt(pow(globalOrigin.cx - pos->x, 2.0) + pow(globalOrigin.cy - pos->y, 2.0));
        if (cr > AUTOSWARM_DEADZONE){
            double gScalar      = 0.5 * AUTOSWARM_GLOBAL * (1 - 1 / (1 + exp(4 / AUTOSWARM_SEPERATION * (cr - AUTOSWARM_SEPERATION))));
            gi[0]               = gScalar * AUTOSWARM_VMAX * (globalOrigin.cx - pos->x) / cr;
            gi[1]               = gScalar * AUTOSWARM_VMAX * (globalOrigin.cy - pos->y) / cr;
        }
        break;
    }
    case AUTOSWARM_CIRCLE_CW :
    case AUTOSWARM_CIRCLE_CC :{
        double cr     = sqrt(pow(globalOrigin.cx - pos->x, 2.0) + pow(globalOrigin.cy - pos->y, 2.0));
        if (cr > AUTOSWARM_DEADZONE){
            double angle;
            if (cr >= AUTOSWARM_CIRCLE_R){
                angle   = 90 * AUTOSWARM_CIRCLE_R / cr;
            }
            else{
                angle   = 90 + 90 * (AUTOSWARM_CIRCLE_R - cr) / AUTOSWARM_CIRCLE_R;     // From 90 to 180 over the length of AUTOSWARM_CIRCLE_R
            }
            if (AUTOSWARM_ATTRACTOR != AUTOSWARM_CIRCLE_CC)     angle = -angle;         // Switch between cc (+) and cw (-)
            double xContrib     = AUTOSWARM_GLOBAL * AUTOSWARM_VMAX * (globalOrigin.cx - pos->x) / cr;
            double yContrib     = AUTOSWARM_GLOBAL * AUTOSWARM_VMAX * (globalOrigin.cy - pos->y) / cr;
            gi[0]               = cos(angle / 180 * M_PI) * xContrib - sin(angle / 180 * M_PI) * yContrib;
            gi[1]               = sin(angle / 180 * M_PI) * xContrib + cos(angle / 180 * M_PI) * yContrib;
        }
        break;
    }
    default : break;
    }
    return;
}

void calcDiffVelocity(double totV[3], double di[3]){
    di[0]       = - AUTOSWARM_EPS * (totV[0] - prev_v_d[0]);
    di[1]       = - AUTOSWARM_EPS * (totV[1] - prev_v_d[1]);
    prev_v_d[0] = totV[0];
    prev_v_d[1] = totV[1];
    prev_v_d[2] = totV[2];
    return;
}

void calcCamPosition(struct NedCoor_f *pos, double totV[3], double cPos[3], double gi[3]){
    if (AUTOSWARM_MODE == 1){                        // Watch where you're going mode. Point WP_CAM on unity circle towards totV
        double velR = sqrt(pow(totV[0], 2.0) + pow(totV[1], 2.0));
        if (velR > 0){
            cPos[0] = pos->x + totV[0] / velR;
            cPos[1] = pos->y + totV[1] / velR;
            cPos[2] = pos->z;
        }
        else{
            cPos[0] = pos->x + cos(eulerAngles->psi) * 1;
            cPos[1] = pos->y + sin(eulerAngles->psi) * 1;
            cPos[2] = pos->z;
        }
    }
    if (AUTOSWARM_MODE == 2){                        // Watch the global objective mode. Point WP_CAM on unity circle towards global component
        double velR = sqrt(pow(gi[0], 2.0) + pow(gi[1], 2.0));
        if (velR > 0){
            cPos[0]     = pos->x + gi[0] / velR;
            cPos[1]     = pos->y + gi[1] / velR;
            cPos[2]     = pos->z;
        }
        else{
            cPos[0]     = pos->x + cos(eulerAngles->psi) * 1;
            cPos[1]     = pos->y + sin(eulerAngles->psi) * 1;
            cPos[2]     = pos->z;
        }
    }
    limitYaw(pos, cPos);                            // Limit yaw
    return;
}

void limitYaw(struct NedCoor_f *pos, double cPos[3]){
    double cameraHeading    = atan2(cPos[1] - pos->y, cPos[0] - pos->x);
    double relHeading       = cameraHeading - eulerAngles->psi; // Z axis is defined downwards for psi so * -1 for the atan2

    if (relHeading > M_PI)       relHeading -= 2 * M_PI;
    if (relHeading < -M_PI)      relHeading += 2 * M_PI;

    if (abs(relHeading) > AUTOSWARM_YAWRATEMAX  / ((double) AUTOSWARM_FPS) / 180 * M_PI){
        double newHeading   = relHeading/abs(relHeading) * (AUTOSWARM_YAWRATEMAX / ((double) AUTOSWARM_FPS) / 180 * M_PI) + eulerAngles->psi;
        cPos[0]             = pos->x + cos(newHeading) * 1;
        cPos[1]             = pos->y + sin(newHeading) * 1;
        cPos[2]             = pos->z;
    }
    return;
}

void limitVelocityYaw(double totV[3]){
    double vTurn;
    double totVHeading  = atan2(totV[1], totV[0]);
    double relHeading   = totVHeading - eulerAngles->psi; // Z axis is defined downwards for psi so * -1 for the atan2

    if (abs(relHeading) > acos(0.01))    vTurn = 0.01            * sqrt(pow(totV[0], 2.0) + pow(totV[1], 2.0));
    else                                vTurn = cos(relHeading) * sqrt(pow(totV[0], 2.0) + pow(totV[1], 2.0));

    totV[0] = cos(totVHeading) * vTurn;
    totV[1] = sin(totVHeading) * vTurn;
    return;
}

void limitNorm(double totV[3], double maxNorm){
    double vR = sqrt(pow(totV[0],2.0) + pow(totV[1],2.0));      // Find the distance to our new desired position
    if (vR > maxNorm){                                          // Is the new desired position unachievable due to velocity constraint?
        totV[0] = maxNorm * totV[0] / vR;                       // Scale the totV X component so that ||totV|| = maxNorm
        totV[1] = maxNorm * totV[1] / vR;                       // Do the same for Y
        totV[2] = 0.0;
    }
    return;
}

void updateWaypoints(struct NedCoor_f *pos, double totV[3], double cPos[3]){
    /*
     *  BE CAREFUL! waypoint_set_xy_i is in ENU so N(x) and E(y) axis are changed compared to x_w and y_w
     */
    waypoint_set_xy_i((unsigned char) WP__GOAL, POS_BFP_OF_REAL(pos->y + totV[1]),  POS_BFP_OF_REAL(pos->x + totV[0]));     // Update WP_GOAL waypoint to add totV relative to our position
    waypoint_set_xy_i((unsigned char) WP__CAM,  POS_BFP_OF_REAL(cPos[1]),           POS_BFP_OF_REAL(cPos[0]));              // Update WP_GOAL waypoint to add totV relative to our position
    setGlobalOrigin(waypoint_get_y((unsigned char) WP_GLOBAL), waypoint_get_x((unsigned char) WP_GLOBAL), 1);

    if (!strcmp("Swarm",flight_blocks[nav_block]) || !strcmp("Swarm Home",flight_blocks[nav_block])){
        nav_set_heading_towards_waypoint(WP__CAM);      // Currently in block "Swarm" or "Swarm Home" so update heading
    }
    if (AUTOSWARM_SHOW_WAYPOINT) printf("WP_CAM (%0.2f m, %0.2f m) \tWP_GOAL (%0.2f m, %0.2f m) \tWP_GLOBAL (%0.2f m, %0.2f m)\n", cPos[0], cPos[1], pos->x + totV[0], pos->y + totV[1], globalOrigin.cx, globalOrigin.cy);    //, pos->z + totV[2]);    // Print to terminal
}

bool amIhome(void){
    struct EnuCoor_f *pos = stateGetPositionEnu_f();  // Get current position
    if (sqrt(pow(pos->x - waypoint_get_x((unsigned char) WP__TD), 2.0) + pow(pos->y - waypoint_get_y((unsigned char) WP__TD), 2.0)) < AUTOSWARM_HOME){
        return true;
    }
    else{
        return false;
    }
}

void autoswarm_opencv_init(int globalMode){
#if AUTOSWARM_WRITE_RESULTS
    startTime = time(0);
    tm * startTM    = localtime(&startTime);
    sprintf(resultFile, "/data/ftp/internal_000/AS_result-%d-%02d-%02d_%02d-%02d-%02d.txt", startTM->tm_year + 1900, startTM->tm_mon + 1, startTM->tm_mday, startTM->tm_hour, startTM->tm_min, startTM->tm_sec);
    pFile           = fopen(resultFile,"w");
    if (pFile == NULL){
        perror("[AS-ERROR] File error");
    }
    else{
        fprintf(pFile,"ID\tmem\ti\tt\tposX\tposY\tposZ\tobjX\tobjY\tobjZ\n");
        fclose(pFile);
        printf("[AS] Writing tracking results to: %s\n", resultFile);
    }
#endif
    AUTOSWARM_ATTRACTOR = globalMode;
    AR_FILTER_CAM_RANGE = AUTOSWARM_SEPERATION * AUTOSWARM_LATTICE_RATIO;
    printf("[AS] initialized, set ar_filter cam_range to %0.2f\n",AR_FILTER_CAM_RANGE);
    return;
}

void autoswarm_opencv_run_header(void){
#if AUTOSWARM_WRITE_RESULTS
    currentTime = time(0);                                                  // Get the current time
    curT        = difftime(currentTime,startTime);                          // Calculate time-difference between startTime and currentTime
#endif
}

void autoswarm_opencv_run_trailer(){
#if AUTOSWARM_WRITE_RESULTS                                                 // See if we want to write results
    pFile         = fopen("/data/ftp/internal_000/results.txt","a");        // Open file for appending
    if (pFile == NULL){
        perror("[AS-ERROR] File error");
    }
#endif
#if AUTOSWARM_WRITE_RESULTS || AUTOSWARM_SHOW_MEM
        for (unsigned int r=0; r < neighbourMem_size; r++){                // Print to file & terminal
#if AUTOSWARM_WRITE_RESULTS || AUTOSWARM_SHOW_MEM
            int memInt = 0;
            if (neighbourMem[r].lastSeen == runCount) memInt = 1;
#endif
#if AUTOSWARM_SHOW_MEM
            printf("%i - Object (see %i  id %i) at (%0.2f m, %0.2f m, %0.2f m) cur_pos (%0.2f, %0.2f, %0.2f)\n", runCount, memInt, neighbourMem[r].id, neighbourMem[r].x_w, neighbourMem[r].y_w, neighbourMem[r].z_w, pos->x, pos->y, pos->z);                                                         // Print to terminal
#endif
#if AUTOSWARM_WRITE_RESULTS
            if (pFile != NULL) fprintf(pFile, "%i\t%i\t%i\t%i\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\n", neighbourMem[r].id, memInt, runCount, curT, pos->x, pos->y, pos->z, neighbourMem[r].x_w, neighbourMem[r].y_w, neighbourMem[r].z_w); // if file writing is enabled, write to file
#endif
        }
#endif
#if AUTOSWARM_WRITE_RESULTS
        if (pFile != NULL) fclose(pFile);    // Close file
#endif
    if ((AUTOSWARM_SHOW_MEM==1 && neighbourMem_size > 0) || AUTOSWARM_SHOW_WAYPOINT || AUTOSWARM_BENCHMARK) printf("\n"); // Separate terminal output by newline
}

#if AUTOSWARM_BENCHMARK
void initBenchmark(){
    bench_start     = clock();
    benchmark_time.clear();
    benchmark_title.clear();
}

void addBenchmark(const char * title){
    bench_end       = clock();
    double time     = ((double) (bench_end - bench_start)) / 1000000;
    bench_start     = bench_end;
    benchmark_time.push_back(time);
    benchmark_title.push_back(title);
}

void showBenchmark(){
    double totalTime        = 0;
    unsigned int maxlength  = 0;
    for (unsigned int i=0; i < benchmark_time.size(); i++){
        totalTime += benchmark_time[i];                                     // Sum the total time elapsed
        if (runCount==1){
            benchmark_avg_time.push_back(0);                                // Prevent benchmark average time from being empty when t=0
        }
        if (strlen(benchmark_title[i]) > maxlength){
            maxlength = strlen(benchmark_title[i]);                         // This is the length of the longest title
        }
    }
    for (unsigned int i=0; i < benchmark_time.size(); i++){
        benchmark_avg_time[i]   = (benchmark_avg_time[i] * (runCount - 1) + benchmark_time[i]/totalTime*100) / runCount;    // Elapsed time in percentage
        if (strlen(benchmark_title[i]) < maxlength - 18){
            printf("%s\t\t\t\t%2.2f percent\n", benchmark_title[i], benchmark_avg_time[i]);     // 4 tabs
        }
        else{
            if (strlen(benchmark_title[i]) < maxlength - 12){
                printf("%s\t\t\t%2.2f percent\n", benchmark_title[i], benchmark_avg_time[i]);   // 3 tabs
            }
            else{
                if (strlen(benchmark_title[i]) < maxlength - 5){
                    printf("%s\t\t%2.2f percent\n", benchmark_title[i], benchmark_avg_time[i]); // 2 tabs
                }
                else{
                    printf("%s\t%2.2f percent\n", benchmark_title[i], benchmark_avg_time[i]);   // 1 tab
                }
            }
        }

    }
}
#endif
