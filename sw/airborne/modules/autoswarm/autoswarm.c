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

/// Header Files ///
#include "state.h"                                              ///< Used for accessing state variables
#include "navigation.h"                                         ///< Used for navigation functions
#include "generated/flight_plan.h"                              ///< Used for WP definitions
#include "autoswarm/autoswarm.h"                                ///< Include own header file
#include "computer_vision/active_random_filter.h"               ///< Include active random filter header file

#define PRINT(string,...) fprintf(stderr, "[autoswarm->%s()] " string,__FUNCTION__ , ##__VA_ARGS__)

#ifndef AS_VERBOSE
#define AS_VERBOSE FALSE
#endif
PRINT_CONFIG_VAR(AS_VERBOSE)

#if AS_VERBOSE
#define VERBOSE_PRINT PRINT
#else
#define VERBOSE_PRINT(...)
#endif

#ifndef AS_GLOBAL_ATTRACTOR
#define AS_GLOBAL_ATTRACTOR AS_CIRCLE_CW
#endif
PRINT_CONFIG_VAR(AS_GLOBAL_ATTRACTOR)

#ifndef AS_GLOBAL_STR
#define AS_GLOBAL_STR 1.0
#endif
PRINT_CONFIG_VAR(AS_GLOBAL_STR)

#ifndef AS_SEPARATION
#define AS_SEPARATION 1.3
#endif
PRINT_CONFIG_VAR(AS_SEPARATION)

#ifndef AS_CIRCLE_RADIUS
#define AS_CIRCLE_RADIUS 1.5
#endif
PRINT_CONFIG_VAR(AS_CIRCLE_RADIUS)

#ifndef AS_VMAX
#define AS_VMAX 2.5
#endif
PRINT_CONFIG_VAR(AS_VMAX)

#ifndef AS_HEADING_MODE
#define AS_HEADING_MODE AS_CAM_GLOBAL
#endif
PRINT_CONFIG_VAR(AS_CAM_GLOBAL)

#ifndef AS_E
#define AS_E 0.0025
#endif
PRINT_CONFIG_VAR(AS_E)

#ifndef AS_EPS
#define AS_EPS 0.03
#endif
PRINT_CONFIG_VAR(AS_EPS)

#ifndef AS_HOME
#define AS_HOME 0.4
#endif
PRINT_CONFIG_VAR(AS_HOME)

#ifndef AS_LOGLO
#define AS_LOGLO 1.75
#endif
PRINT_CONFIG_VAR(AS_LOGLO)

#ifndef AS_LATTICE_RATIO
#define AS_LATTICE_RATIO 2.0
#endif
PRINT_CONFIG_VAR(AS_LATTICE_RATIO)

#ifndef AS_DEADZONE
#define AS_DEADZONE 0.25
#endif
PRINT_CONFIG_VAR(AS_DEADZONE)

/// Debug options ///
#ifndef AS_PRINT_WAYPOINT
#define AS_PRINT_WAYPOINT 0                                     ///< Show the updated positions of the waypoints
#endif
PRINT_CONFIG_VAR(AS_PRINT_WAYPOINT)

#ifndef AS_WRITE_RESULTS
#define AS_WRITE_RESULTS 0                                      ///< Write measurements to text file
#endif
PRINT_CONFIG_VAR(AS_WRITE_RESULTS)

/// Static function declarations ///
static void calculateVelocityResponse(struct NedCoor_f *pos, struct FloatEulers* eulerAngles, double totV[3], double gi[3]);
static void calculateGlobalVelocity(struct NedCoor_f *pos, double gi[3]);
static void calculateLocalVelocity(struct NedCoor_f *pos, double li[3]);
static void calculateDiffVelocity(double totV[3], double di[3]);
static void calculateCamPosition(struct NedCoor_f *pos, struct FloatEulers* eulerAngles, double totV[3], double gi[3], double cPos[3]);
static void updateWaypoints(struct NedCoor_f *pos, double totV[3], double cPos[3]);
static void limitNorm(double totV[3], double maxNorm);
static void limitVelocityYaw(struct FloatEulers* eulerAngles, double totV[3]);
static void autoswarm_write_log(void);
static double calculateLocalGlobalCoeff(struct NedCoor_f *pos, double gi[3]);

/// Set up swarm parameters ///
int     settings_as_heading_mode    = AS_HEADING_MODE;          ///< Heading/Camera control modus
double  settings_as_separation      = AS_SEPARATION;            ///< Separation between neighbours
double  settings_as_lattice_ratio   = AS_LATTICE_RATIO;         ///< Local interactions are bounded within lattice_ratio * separation

double  settings_as_e               = AS_E;                     ///< Pinciroli E coefficient
double  settings_as_eps             = AS_EPS;                   ///< Differential component
double  settings_as_loglo           = AS_LOGLO;                 ///< Local-Global interaction exponent

double  settings_as_global_strength = AS_GLOBAL_STR;            ///< Global field strength (% of V_MAX)
double  settings_as_vmax            = AS_VMAX;                  ///< Maximum goal waypoint translation

/// Set up global attractor parameters ///
int     settings_as_attractor       = AS_GLOBAL_ATTRACTOR;      ///< The type of global attractor/field
double  settings_as_circle_radius   = AS_CIRCLE_RADIUS;         ///< Radius of the circle field
double  settings_as_deadzone        = AS_DEADZONE;              ///< Deadzone near the centre of the global field to prevent spinning around the centre

/// Initialize parameters to be assigned during runtime ///
static uint32_t             runCount        = 0;                ///< The current frame number
static const char *         flight_blocks[] = FP_BLOCKS;        ///< Array of flight blocks
static double               prev_v_d[3]     = {0.0, 0.0, 0.0};  ///< Previous totV used for differential component

/// Optional variables declarations ///
#if AS_WRITE_RESULTS
static time_t               startTime;                          ///< Start time of the program
static time_t               currentTime;                        ///< Current time
static int                  curT;                               ///< Current difference between start time and current time aka runtime
FILE*                       pFile;                              ///< File handle to write to
char                        resultFile [50];                    ///< File name to write to
#endif // AS_WRITE_RESULTS

/// Function definitions ///

/** initialize autoswarm
 *
 * @param[in] initial global mode (uint8_t)
 *
 * This function sets the default global mode, cam_range and log-file
 */
void autoswarm_init( void ){
#if AS_WRITE_RESULTS
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
    printf("[AS] initialized, set ar_filter cam_range to %0.2f\n",ARF_CAM_RANGE);
    return;
}

/** main function of autoswarm
 *
 * This function calculates the desired waypoints and sets the correct heading
 */
void autoswarm_run( void ){
    /*
    if (strcmp("Swarm",flight_blocks[nav_block]) && strcmp("Swarm Home",flight_blocks[nav_block])){
        ///< Don't run if we are not in flight block "Swarm" or "Swarm Home"
        return;
    }
    */
#if AS_WRITE_RESULTS
    currentTime     = time(0);                                              ///< Get the current time
    curT            = difftime(currentTime,startTime);                      ///< Calculate time-difference between startTime and currentTime
#endif
    struct FloatEulers* eulerAngles = stateGetNedToBodyEulers_f();          ///< Get vehicle body attitude euler angles (float).
    struct NedCoor_f* pos           = stateGetPositionNed_f();              ///< Get position in local NED coordinates (float).
    double totV[3]  = {0.0, 0.0, 0.0};                                      ///< Desired velocity
    double cPos[3]  = {0.0, 0.0, 0.0};                                      ///< Desired camera position in unit distance world coordinates
    double gi[3]    = {0.0, 0.0, 0.0};                                      ///< Global interaction vector
    calculateVelocityResponse(pos, eulerAngles, totV, gi);                  ///< Calculate the velocity response of the agent and store result in totV and gi
    calculateCamPosition(pos, eulerAngles, totV, gi, cPos);                 ///< Calculate WP_CAM/heading based on pos,totV and gi and store in cPos
    updateWaypoints(pos, totV, cPos);                                       ///< Update waypoints based on velocity contribution
    autoswarm_write_log();                                           ///< Write to log file
    runCount++;                                                             ///< Update global run-counter
    return;
}

/** Calulates the velocity response
 *
 * @param[in] Current position (struct NedCoor_f*)
 * @param[in] Array to store desired velocity (double[3])
 * @param[in] Array to store the global components (double[3])
 *
 * This function calculates velocity response based on the local, global, and differential components
 */
void calculateVelocityResponse(struct NedCoor_f *pos, struct FloatEulers* eulerAngles, double totV[3], double gi[3]){
    double li[3]    = {0.0, 0.0, 0.0};                      ///< Local contribution
    double di[3]    = {0.0, 0.0, 0.0};                      ///< Differential contribution
    calculateLocalVelocity(pos, li);                        ///< Get the local contribution based on our neighbours
    calculateGlobalVelocity(pos, gi);                       ///< Get the contribution due to the global attractor/field
    double g_coeff  = calculateLocalGlobalCoeff(pos, gi);   ///< Calculate the local-global interaction coefficient
    totV[0]         = li[0] + g_coeff * gi[0];              ///< Scale the global component with the local-global interaction coefficient and add them
    totV[1]         = li[1] + g_coeff * gi[1];
    limitNorm(totV, settings_as_vmax);                      ///< Check if the desired velocity exceeds settings_as_vmax
    calculateDiffVelocity(totV, di);                        ///< Calculate the differential component based on change in desired velocity
    totV[0]         = totV[0] + di[0];                      ///< Add differential component
    totV[1]         = totV[1] + di[1];
    if (settings_as_heading_mode!=AS_CAM_GLOBAL){
        limitVelocityYaw(eulerAngles, totV);                ///< Limit the velocity when heading change gets larger
    }
    return;
}

/** Calulates the local velocity response
 *
 * @param[in] Current position (struct NedCoor_f*)
 * @param[in] Array to store local interaction velocity (double[3])
 *
 * This function calculates local velocity response based on our neighbours
 */
void calculateLocalVelocity(struct NedCoor_f *pos, double li[3]){
#ifdef __linux__
  pthread_mutex_lock(&neighbourMem_mutex);
#endif
    for (uint8_t i=0; i < neighbourMem_size; i++){
        double r    = sqrt( pow( ( pos->x - neighbourMem[i].x_w ), 2.0 ) + pow( ( pos->y - neighbourMem[i].y_w ), 2.0 ) );                  ///< The range to this neighbour
        if(r < settings_as_lattice_ratio * settings_as_separation){
            double li_n = 12 * settings_as_e / r * ( pow( settings_as_separation / r, 12.0 ) - pow( settings_as_separation / r, 6.0 ) );    ///< Local contribution scaling factor
            li[0]      += li_n * (pos->x - neighbourMem[i].x_w) / r;                                                                        ///< Local contribution in x
            li[1]      += li_n * (pos->y - neighbourMem[i].y_w) / r;                                                                        ///< Local contribution in y
        }
    }
    if (neighbourMem_size > 0){
        li[0] = li[0] / neighbourMem_size;                                                                                                  ///< Average the local contribution (nr. neighbours independent)
        li[1] = li[1] / neighbourMem_size;
        li[2] = 0.0;
    }
    else{
        li[0] = 0.0;
        li[1] = 0.0;
        li[2] = 0.0;
    }
#ifdef __linux__
  pthread_mutex_unlock(&neighbourMem_mutex);
#endif
    return;
}

/** Calulates the local-global interaction coefficient
 *
 * @param[in] Current position (struct NedCoor_f*)
 * @param[in] Array containing global interaction velocity (double[3])
 * @param[out] Local-global interaction coefficient (double)
 *
 * This function calculates the local-global interaction coefficient
 */
double calculateLocalGlobalCoeff(struct NedCoor_f *pos, double gi[3]){
#ifdef __linux__
    pthread_mutex_lock(&neighbourMem_mutex);
#endif
    double result = 0.0;
    if (neighbourMem_size > 0){
        double w_neighbours[2]  = {0.0, 0.0};                                                                                               ///< weighed neighbours vector
        for (unsigned int r = 0; r < neighbourMem_size; r++){
            double range            = sqrt( pow( ( pos->x - neighbourMem[r].x_w ), 2.0 ) + pow( ( pos->y - neighbourMem[r].y_w ), 2.0 ) );  ///< Range of neighbour
            if(range < settings_as_separation * settings_as_lattice_ratio){
                double rel_dist         = fmin( 1.0, ( settings_as_lattice_ratio * settings_as_separation - range )
                        / ( settings_as_separation * ( settings_as_lattice_ratio - 1 ) ) );                             ///< Relative distance coefficient
                w_neighbours[0]        += pow( rel_dist, 2.0 ) * ( pos->x - neighbourMem[r].x_w ) / range;                                  ///< Calculate neighbour contribution to weighed neighbours
                w_neighbours[1]        += pow( rel_dist, 2.0 ) * ( pos->y - neighbourMem[r].y_w ) / range;
            }
        }
        double gi_n             = sqrt( pow( gi[0], 2.0 ) + pow( gi[1], 2.0 ) );                                                            ///< Norm of global interaction
        w_neighbours[0] *= gi_n;                                                                                                            ///< Scale weighed neighbours vector by global contribution
        w_neighbours[1] *= gi_n;
        double c_gi_coef   = ( w_neighbours[0] * gi[0] + w_neighbours[1] * gi[1] ) / pow( gi_n, 2.0 );                                      ///< dot product between weighed neighbours and global interaction, divided by dot produgt of global interaction and global interaction;
        result = pow( fmin( 1.0, fmax( 0.0,
                2.0 - sqrt( pow( gi[0] * ( 1 - c_gi_coef), 2.0 ) + pow( gi[1] * (1 - c_gi_coef), 2.0 ) ) / gi_n ) ), settings_as_loglo );   ///< Return local-global interaction coefficient
    }
#ifdef __linux__
    pthread_mutex_unlock(&neighbourMem_mutex);
#endif
    return result;
}

/** Calulates the global velocity component
 *
 * @param[in] Current position (struct NedCoor_f*)
 * @param[in] Array to store global interaction velocity (double[3])
 *
 * This function calculates the global velocity component
 */
void calculateGlobalVelocity(struct NedCoor_f *pos, double gi[3]){
    switch (settings_as_attractor){
    case AS_POINT :{
       double cr     = sqrt( pow( globalOrigin.cx - pos->x, 2.0 ) + pow( globalOrigin.cy - pos->y, 2.0 ) );
        if (cr > settings_as_deadzone){
            gi[0]               = settings_as_global_strength * settings_as_vmax * ( globalOrigin.cx - pos->x ) / cr;
            gi[1]               = settings_as_global_strength * settings_as_vmax * ( globalOrigin.cy - pos->y ) / cr;
        }
        break;
    }
    case AS_BUCKET :{
        double cr     = sqrt( pow( globalOrigin.cx - pos->x, 2.0 ) + pow( globalOrigin.cy - pos->y, 2.0 ) );
        if (cr > settings_as_deadzone){
            double gScalar      = 0.5 * settings_as_global_strength * ( 1 - 1 / ( 1 + exp( 4 / settings_as_separation * (cr - settings_as_separation ) ) ) );
            gi[0]               = gScalar * settings_as_vmax * ( globalOrigin.cx - pos->x ) / cr;
            gi[1]               = gScalar * settings_as_vmax * ( globalOrigin.cy - pos->y ) / cr;
        }
        break;
    }
    case AS_CIRCLE_CW :
    case AS_CIRCLE_CC :{
        double cr     = sqrt( pow( globalOrigin.cx - pos->x, 2.0 ) + pow( globalOrigin.cy - pos->y, 2.0 ) );
        if (cr > settings_as_deadzone){
            double angle;
            if (cr >= settings_as_circle_radius){
                angle   = 90 * settings_as_circle_radius / cr;
            }
            else{
                angle   = 90 + 90 * (settings_as_circle_radius - cr) / settings_as_circle_radius;     // From 90 to 180 over the length of settings_as_circle_radius
            }
            if (settings_as_attractor != AS_CIRCLE_CC)     angle = -angle;         // Switch between cc (+) and cw (-)
            double xContrib     = settings_as_global_strength * settings_as_vmax * ( globalOrigin.cx - pos->x ) / cr;
            double yContrib     = settings_as_global_strength * settings_as_vmax * ( globalOrigin.cy - pos->y ) / cr;
            gi[0]               = cos( angle / 180 * M_PI ) * xContrib - sin( angle / 180 * M_PI ) * yContrib;
            gi[1]               = sin( angle / 180 * M_PI ) * xContrib + cos( angle / 180 * M_PI ) * yContrib;
        }
        break;
    }
    default : break;
    }
    return;
}

/** Calulates the differential velocity component
 *
 * @param[in] Current position (struct NedCoor_f*)
 * @param[in] Array to store differential interaction velocity (double[3])
 *
 * This function calculates the differential velocity component
 */
void calculateDiffVelocity(double totV[3], double di[3]){
    di[0]       = - settings_as_eps * ( totV[0] - prev_v_d[0] );
    di[1]       = - settings_as_eps * ( totV[1] - prev_v_d[1] );
    prev_v_d[0] = totV[0];
    prev_v_d[1] = totV[1];
    prev_v_d[2] = totV[2];
    return;
}

/** Calulates the desired camera position / heading
 *
 * @param[in] Current position (struct NedCoor_f*)
 * @param[in] Array containing desired velocity (double[3])
 * @param[in] Array to store desired camera position as a unit vector from the current position (double[3])
 * @param[in] Array containing global velocity component (double[3])
 *
 * This function calculates the differential velocity component
 */
void calculateCamPosition(struct NedCoor_f *pos, struct FloatEulers* eulerAngles, double totV[3], double gi[3], double cPos[3]){
    if (settings_as_heading_mode == AS_CAM_FORWARD){                     ///< Point WP_CAM on unity circle towards totV
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
    if (settings_as_heading_mode == AS_CAM_GLOBAL){                        ///< Point WP_CAM on unity circle towards global component
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
    return;
}

/** Limits the norm of the desired velocity based on relative heading change
 *
 * @param[in] Array containing desired velocity (double[3])
 *
 * This function limits the velocity when increasing yaw changes are made
 */
void limitVelocityYaw(struct FloatEulers* eulerAngles, double totV[3]){
    double totVHeading  = atan2(totV[1], totV[0]);              ///< Heading of totV
    double relHeading   = totVHeading - eulerAngles->psi;       ///< Relative angle between totV and heading
    double vTurn;
    if (abs(relHeading) > acos(0.01)){
        vTurn = 0.01            * sqrt(pow(totV[0], 2.0) + pow(totV[1], 2.0));
    }
    else{
        vTurn = cos(relHeading) * sqrt(pow(totV[0], 2.0) + pow(totV[1], 2.0));
    }
    totV[0] = cos(totVHeading) * vTurn;
    totV[1] = sin(totVHeading) * vTurn;
    return;
}

/** Limits the norm the input vector based on the maximum given
 *
 * @param[in] Vector (double[3])
 * @param[in] maximum norm (double)
 *
 * This function limits the input vector to a maximum given
 */
void limitNorm(double vec[3], double maxNorm){
    double vR = sqrt(pow(vec[0],2.0) + pow(vec[1],2.0));      ///< Find the norm
    if (vR > maxNorm){
        vec[0] = maxNorm * vec[0] / vR;                       ///< vector so that ||vector|| = maxNorm
        vec[1] = maxNorm * vec[1] / vR;
        vec[2] = 0.0;
    }
    return;
}

/** Updates the waypoints and heading based on desired velocity and camera position
 *
 * @param[in] Current position (struct NedCoor_f*)
 * @param[in] Array containing desired velocity (double[3])
 * @param[in] Array containing desired camera position (double[3])
 *
 * This function updates waypoints and heading according to the given desired velocity and camera position
 */
void updateWaypoints(struct NedCoor_f *pos, double totV[3], double cPos[3]){
    /*
     *  BE CAREFUL! waypoint_set_xy_i is in ENU so N(x) and E(y) axis are changed compared to x_w and y_w
     */
    waypoint_set_xy_i((unsigned char) WP__GOAL, POS_BFP_OF_REAL(pos->y + totV[1]),  POS_BFP_OF_REAL(pos->x + totV[0]));     ///< Update WP_GOAL waypoint to add desired velocity to our position
    waypoint_set_xy_i((unsigned char) WP__CAM,  POS_BFP_OF_REAL(cPos[1]),           POS_BFP_OF_REAL(cPos[0]));              ///< Update WP_GOAL waypoint to add camera position to our position
    setGlobalOrigin(waypoint_get_y((unsigned char) WP_GLOBAL), waypoint_get_x((unsigned char) WP_GLOBAL), 1);               ///< Update global origin based on WP_GLOBAL location
    if (!strcmp("Swarm",flight_blocks[nav_block]) || !strcmp("Swarm Home",flight_blocks[nav_block])){
        nav_set_heading_towards_waypoint(WP__CAM);                                                                          ///< Currently in block "Swarm" or "Swarm Home" so update heading
    }
    if (AS_PRINT_WAYPOINT) PRINT("WP_CAM (%0.2f m, %0.2f m) \tWP_GOAL (%0.2f m, %0.2f m) \tWP_GLOBAL (%0.2f m, %0.2f m)\n", cPos[0], cPos[1], pos->x + totV[0], pos->y + totV[1], globalOrigin.cx, globalOrigin.cy);
}

/** Checks if position is close enough to waypoint _TD
 *
 * @param[in] Is our position close enough to waypoint _TD (bool)
 *
 * This function checks is our current position is close enough to waypoint _TD based on the AS_HOME distance
 */
bool amIhome(void){
    struct EnuCoor_f *pos = stateGetPositionEnu_f();  // Get current position
    if (sqrt(pow(pos->x - waypoint_get_x((unsigned char) WP__TD), 2.0) + pow(pos->y - waypoint_get_y((unsigned char) WP__TD), 2.0)) < AS_HOME){
        return true;
    }
    else{
        return false;
    }
}

/** Writes a log to the disk
 *
 * This function writes a log to the disk //TODO: fix
 */
void autoswarm_write_log(){
#if AS_WRITE_RESULTS                                                 // See if we want to write results
#ifdef __linux__
  pthread_mutex_lock(&neighbourMem_mutex);
#endif
    pFile         = fopen("/data/ftp/internal_000/results.txt","a");        // Open file for appending
    if (pFile == NULL){
        perror("[AS-ERROR] File error");
    }
    for (unsigned int r=0; r < neighbourMem_size; r++){                // Print to file & terminal
        int memInt = 0;
        if (neighbourMem[r].lastSeen == runCount) memInt = 1;
        if (pFile != NULL) fprintf(pFile, "%i\t%i\t%i\t%i\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\n", neighbourMem[r].id, memInt, runCount, curT, pos->x, pos->y, pos->z, neighbourMem[r].x_w, neighbourMem[r].y_w, neighbourMem[r].z_w); // if file writing is enabled, write to file
    }
    if (pFile != NULL) fclose(pFile);    // Close file
#ifdef __linux__
  pthread_mutex_unlock(&neighbourMem_mutex);
#endif
#endif
    if (AS_PRINT_WAYPOINT) printf("\n"); // Separate terminal output by newline
}
