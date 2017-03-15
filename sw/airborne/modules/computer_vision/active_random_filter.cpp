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
 * @file "modules/computer_vision//cv_active_random_filter.c"
 * @author Wilco Vlenterie (wv-tud)
 * Active random sampling colour filter
 */

#include "active_random_filter.h"
#include "bebop_camera_stabilization.h"
#include <vector>
#include <ctime>

#define BOARD_CONFIG "boards/bebop.h"               ///< Define which board

extern "C" {
    #include "boards/bebop.h"                       ///< C header used for bebop specific settings
    #include <state.h>                              ///< C header used for state functions and data
    #include <sys/time.h>                           ///< C header used for system time functions and data
    #include "mcu_periph/sys_time.h"                ///< C header used for PPRZ time functions and data
}

using namespace std;
#include <opencv2/core/core.hpp>                    ///< Load openCV
#include <opencv2/imgproc/imgproc.hpp>              ///< Load openCV image processing library
using namespace cv;


#define PRINT(string,...) fprintf(stderr, "[AR-FILTER->%s()] " string,__FUNCTION__ , ##__VA_ARGS__)

#define ARF_VERBOSE FALSE
#if ARF_VERBOSE
#define VERBOSE_PRINT PRINT
#else
#define VERBOSE_PRINT(...)
#endif

#define xSign(x) ( ( x ) >= ( 0 ) ? ( 1 ) : ( -1 ) )

#define ARF_MARK_CONTOURS   1                       ///< Mark all contour pixels green on sourceframe
#define ARF_BALL_CIRCLES    1                       ///< Draw circles around balls
#define ARF_GATE_CORNERS    1                       ///< Plot corner points of Gates
#define ARF_PLOT_COORDS     0                       ///< Plot the coordinates of balls on frame
#define ARF_DISTANCE_PLOT   0                       ///< Plot lines with distance on frame
#define ARF_CROSSHAIR       0                       ///< Plot horizon
#define ARF_SHOW_CAM_INFO   1                       ///< Show colour gains and exposure on frame
#define ARF_SHOW_STATS      0                       ///< Show statistics on the performance of the contour detection

#define ARF_MEASURE_FPS     1                       ///< Measure average FPS
#define ARF_TIMEOUT         50                      ///< Frames from start
#define ARF_WRITE_LOG       0                       ///< Write tracking results to logfile

#define ARF_WORLDPOS        0                       ///< Use world coordinates
#define ARF_NOYAW           1                       ///< Output in body horizontal XY
#define ARF_USE_ALTITUDE    0                       ///< Use own altitude for world pos

#define ARF_SHOW_REJECT     0                       ///< Print why shapes are rejected
#define ARF_MOD_VIDEO       1                       ///< Modify the frame to show relevant info
#define ARF_DRAW_BOXES 	    0                       ///< Draw boxes
#define ARF_SHOW_MEM        1                       ///< Print object locations to terminal

#define ARF_BALL            0
#define ARF_GATE            1
#define ARF_OBJECT          ARF_GATE//BALL

#define     ARF_MAX_OBJCONT_SIZE  10000                                 ///< Largest size of a contour allowed
#define     ARF_OBJ_X_OFFSET      0.0                                   ///< Offset x from object centre to object c.g. in world frame
#define     ARF_OBJ_Y_OFFSET      0.0                                   ///< Offset y from object centre to object c.g. in world frame
#define     ARF_OBJ_Z_OFFSET      0.1                                   ///< Offset z from object centre to object c.g. in world frame

extern void             plotHorizon         ( Mat& sourceFrameCrop, struct FloatEulers*  eulerAngles );

static void             active_random_filter_header ( void );
static void             active_random_filter_footer ( struct FloatEulers*  eulerAngles );
static Rect             setISPvars          ( uint16_t width, uint16_t height );
static void 			trackObjects	    ( Mat& sourceFrame, Mat& greyFrame );
#if ARF_OBJECT == ARF_BALL
static void             identifyObject      ( trackResults* trackRes );
#endif
#if ARF_OBJECT == ARF_GATE
void identifyObject(gateResults* gateRes);
#endif
static bool 			addContour			( vector<Point> contour, uint16_t offsetX, uint16_t offsetY, double minDist = 0.0, double maxDist = 0.0);
static void 			cam2body 			( trackResults* trackRes );
static void 			body2world 			( struct FloatEulers*  eulerAngles, trackResults* trackRes );
static void             estimatePosition    ( uint16_t xp, uint16_t yp, uint32_t area, double position[3]);
static bool             getNewPosition      ( uint8_t nextDir, uint16_t* newRow, uint16_t* newCol, int* maxRow, int* maxCol );
static void             eraseMemory         ( void );
static void             getYUVColours       ( Mat& sourceFrame, uint16_t row, uint16_t col, uint8_t* Y, uint8_t* U, uint8_t* V );
static void             createSearchGrid    ( uint16_t x_p, uint16_t y_p, Point searchGrid[], uint8_t searchLayer, uint16_t sGridSize, int* maxRow, int* maxCol);
/** Flood CW declarations **/
static bool             processImage_cw     ( Mat& sourceFrame, Mat& destFrame, uint16_t sampleSize );
static int              pixFindContour_cw   ( Mat& sourceFrame, Mat& destFrame, uint16_t row, uint16_t col, uint8_t prevDir, bool cascade );
static int              pixFollowContour_cw ( Mat& sourceFrame, Mat& destFrame, uint16_t row, uint16_t col, uint8_t prevDir );
static void             getNextDirection_cw ( uint8_t prevDir, uint8_t* nextDir, uint8_t* nextDirCnt );
static bool             objCont_add         ( double minDist = 0.0, double maxDist = 0.0);
static void             objCont_addPoint    ( uint16_t* row, uint16_t* col );
static Moments          objCont_moments     ( void );
/** Flood omni declarations **/
static bool             processImage_omni   ( Mat& sourceFrame, Mat& destFrame, uint16_t sampleSize );
static int              pixFindContour_omni ( Mat& sourceFrame, Mat& destFrame, uint16_t row, uint16_t col, uint8_t prevDir, bool cascade );
static void             getNextDirection_omni( uint8_t prevDir, uint8_t* nextDir, uint8_t* nextDirCnt );
static void             processCrops        ( Mat& frameGrey);
static void             addCrop             ( void );
static Rect             enlargeRectangle    ( Mat& sourceFrame, Rect rectangle, double scale );
static bool             inRectangle         ( Point pt, Rect rectangle );
/** Set up trackRes **/
static uint8_t          trackRes_size = 0;
static bool             trackRes_findMax    ( void );
static bool             trackRes_add        ( trackResults newRes, uint8_t overwriteId = trackRes_size);
static bool             trackRes_clear      ( void );
/** Set up neighbourMem **/
uint8_t                 neighbourMem_size = 0;
static bool             neighbourMem_findMax( void );
static bool             neighbourMem_add    ( memoryBlock newRes, uint8_t overwriteId = neighbourMem_size);
/** Optional functions **/
#if ARF_MOD_VIDEO
static void             mod_video           (Mat& sourceFrame, Mat& frameGrey);
#endif

/** Set up tracking parameters **/
uint16_t    default_calArea                     = 7650;                 ///< Area of a ball at 1m resolution on full sensor resolution
uint8_t     ARF_FLOOD_STYLE                     = ARF_FLOOD_CW;         ///< Flood style to search for contours
uint8_t     ARF_SAMPLE_STYLE                    = ARF_STYLE_RANDOM;     ///< Sample style to search for contours
double      ARF_CAM_RANGE                       = 10.0;                 ///< Maximum camera range of newly added objects
uint16_t    ARF_RND_PIX_SAMPLE                  = 2500;                 ///< Random pixel sample size
uint16_t    ARF_MAX_LAYERS                      = 5000;                 ///< Maximum recursive depth of CW flood
double 	    ARF_MAX_CIRCLE_DEF 	                = 0.15;                 ///< Maximum contour eccentricity
double      ARF_MIN_CIRCLE_PERC                 = 0.45;                 ///< Minimum percentage of circle in view
double      ARF_LARGE_SKIP_FACTOR               = 1.0 / 20.0;           ///< Percentage of length large contours are allowed to snap back to starting pos
/** Automatically calculated tracking parameters **/
uint16_t    ARF_MIN_POINTS;                                             ///< Mimimum contour length
double      ARF_MIN_CIRCLE_SIZE;                                        ///< Minimum contour area
uint16_t    ARF_MIN_LAYERS;                                             ///< Miminum recursive depth of CW flood
uint16_t    ARF_LARGE_LAYERS;                                           ///< Miminum recursive depth of CW flood before classified as a large contour

uint16_t    ARF_MIN_CROP_AREA                   = 100;                  ///< Minimal area of a crop rectangle

/** Set up colour filter **/
#if ARF_OBJECT == ARF_GATE
/* GATE CYBERZOO */
uint8_t     ARF_Y_MIN                           = 10;
uint8_t     ARF_Y_MAX                           = 250;
uint8_t     ARF_U_MIN                           = 100;
uint8_t     ARF_U_MAX                           = 170;
uint8_t     ARF_V_MIN                           = 150;
uint8_t     ARF_V_MAX                           = 220;
uint8_t     ARF_GREY_THRES                      = 0;
#endif
#if ARF_OBJECT == ARF_BALL
/* Cyberzoo */
uint8_t     ARF_Y_MIN                           = 50;                   ///< Minimum Y whilst searching and following contours
uint8_t     ARF_Y_MAX                           = 250;                  ///< Maximum Y whilst searching and following contours
uint8_t     ARF_U_MIN                           = 105;                  ///< Minimum U whilst searching and following contours
uint8_t     ARF_U_MAX                           = 170;                  ///< Maximum U whilst searching and following contours
uint8_t     ARF_V_MIN                           = 150;                  ///< Minimum V whilst searching and following contours
uint8_t     ARF_V_MAX                           = 210;                  ///< Maximum V whilst searching and following contours
uint8_t     ARF_GREY_THRES                      = 7;
/* FAKE LIGHT
uint8_t     ARF_Y_MIN             = 0;
uint8_t     ARF_Y_MAX             = 255;
uint8_t     ARF_U_MIN             = 0;
uint8_t     ARF_U_MAX             = 255;
uint8_t     ARF_V_MIN             = 158;
uint8_t     ARF_V_MAX             = 255;
*/

/* DAYLIGHT
uint8_t     ARF_Y_MIN           = 130;
uint8_t     ARF_Y_MAX           = 255;
uint8_t     ARF_U_MIN           = 95;
uint8_t     ARF_U_MAX           = 131;
uint8_t     ARF_V_MIN           = 145;
uint8_t     ARF_V_MAX           = 188;
*/
#endif

uint8_t     ARF_MAX_SEARCH_PIXEL_SKIP = 6;                              ///< Maximum nr of false pixels to skip whilst searching upwards for contours
/** Set up Remaining parameters **/
double 	    ARF_CROP_X 			                = 1.2;                  ///< Crop margin for blobs when using omni detection
uint8_t     ARF_MEMORY 			                = 40;                   ///< Frames to keep neighbours in memory
double      ARF_FPS                             = 17.0;                 ///< Estimated FPS to estimate lost neighbour decay
double      ARF_VMAX                            = 7.0;                  ///< Maximum estimated velocity of a neighbour (account for some noise)

/** Initialize parameters to be assigned during runtime **/
static uint16_t 	        pixCount            = 0;                    ///< Total pixels processed (resets to 0 before each frame)
static uint16_t 	        pixSucCount         = 0;                    ///< Total pixels passing pixTest whilst following contours (resets to 0 before each frame)
static uint16_t             pixDupCount         = 0;                    ///< Total pixels processed whilst being a duplicate contour (resets to 0 before each frame)
static uint16_t             pixSrcCount         = 0;                    ///< Total pixels processed during search for contours (resets to 0 before each frame)
static uint16_t             pixNofCount         = 0;                    ///< Total pixels failing pixTest (resets to 0 before each frame)
static uint16_t 	        layerDepth          = 0;                    ///< Recursive depth of pixel being processed (resets to 0 before each frame)
static uint16_t 	        sample              = 0;                    ///< Current sample being processed (resets to 0 before each frame)
static uint16_t 	        runCount            = 0;                    ///< Current frame number
static uint16_t 	        maxId 		        = 0;                    ///< Highest nr of unique neighbours trackes (used for assigning unique IDs)
static uint8_t              trackRes_maxId      = 0;                    ///< Index of current farthest tracking result
static double               trackRes_maxVal     = 0;                    ///< Range of current farthest tracking result
static uint8_t              trackRes_lastId     = 0;                    ///< Index of last processed tracking result
static uint8_t              neighbourMem_maxId  = 0;                    ///< Index of current farthest neighbour
static double               neighbourMem_maxVal = 0;                    ///< Range of current farthest neighbour
static uint8_t              neighbourMem_lastId = 0;                    ///< Index of last processed neighbour
static Point                objCont_store[ARF_MAX_OBJCONT_SIZE];        ///< Point in the current contour
static uint16_t             objCont_size        = 0;                    ///< Number of points in the current contour
static uint16_t             objCont_sRow        = 0;                    ///< Starting row of the current contour
static uint16_t             objCont_sCol        = 0;                    ///< Starting column of the current contour
extern uint16_t             ispWidth;                                   ///< Maximum width of ISP after applied scaling
extern uint16_t             ispHeight;                                  ///< Maximum height of ISP after applied scaling
extern uint16_t             initialWidth;                               ///< Initial width of ISP after applied scaling
extern uint16_t             initialHeight;                              ///< Initial height of ISP after applied scaling
extern uint16_t             cropCol;                                    ///< Column from which the ISP is cropped relative to MIN
extern int16_t              fillHeight;                                 ///< Column from which the ISP is cropped relative to MIN
extern double               ispScalar;                                  ///< Applied scalar by the ISP
static Rect 			    objCrop;                                    ///< Current recangle being processed when using omni search
static vector<Rect> 	    cropAreas;                                  ///< All the rectangles found
struct NedCoor_f*           pos;                                        ///< Current position of the UAV

#ifdef __linux__
pthread_mutex_t neighbourMem_mutex;
#endif

#if ARF_OBJECT == ARF_BALL
static trackResults         trackRes[ARF_MAX_OBJECTS];                  ///< Array to store the tracking results
memoryBlock                 neighbourMem[ARF_MAX_OBJECTS];              ///< Array to store the neighbours
#endif
#if ARF_OBJECT == ARF_GATE
static gateResults          trackRes[ARF_MAX_OBJECTS];                  ///< Array to store the tracking results
memoryGateBlock             neighbourMem[ARF_MAX_OBJECTS];              ///< Array to store the neighbours
gateResults hold_Gate[20]; //hold max 20 gates
int gate_count = 0;
#endif

#if ARF_MEASURE_FPS
    static struct timespec      time_now;                               ///< The current time
    static struct timespec      time_prev;                              ///< The time of the previous frame
    static struct timespec      time_init;                              ///< The time the processing began (after timeout)
    static uint32_t curT;                                               ///< The time in us between time_init and time_now
#endif

#if ARF_WRITE_LOG
    FILE *                      arf_File;                               ///< File handle of log file
    char                        arf_FileName[150];                      ///< File name of log file
#endif

void active_random_filter_init(void){
#if ARF_MEASURE_FPS || ARF_WRITE_LOG
    clock_gettime(CLOCK_MONOTONIC, &time_prev);
    time_init = time_prev;
#endif
#if ARF_WRITE_LOG
    time_t startTime    = time(0);
    tm * startTM        = localtime(&startTime);
    sprintf(arf_FileName, "/data/ftp/internal_000/ARF_result-%d-%02d-%02d_%02d-%02d-%02d.txt", startTM->tm_year + 1900, startTM->tm_mon + 1, startTM->tm_mday, startTM->tm_hour, startTM->tm_min, startTM->tm_sec);
    arf_File            = fopen(arf_FileName,"w");
    if (arf_File == NULL){
        perror("[AS-ERROR] File error");
    }
    else{
        fprintf(arf_File,"NR\tUS\tID\tMEM\tPOSX\tPOSY\tPOSZ\tPSI\tOBJX\tOBJY\tOBJZ\n");
        PRINT("Writing tracking results to: %s\n", arf_FileName);
    }
#endif
}

void active_random_filter(char* buff, uint16_t width, uint16_t height, struct FloatEulers* curEulerAngles){
    Mat sourceFrame (height, width, CV_8UC2, buff);                 // Initialize current frame in openCV (UYVY) 2 channel
    Mat frameGrey   (height, width, CV_8UC1, cvScalar(0.0));        // Initialize an empty 1 channel frame
#if ARF_SAVE_FRAME
    if(runCount == 25) saveBuffer(sourceFrame, "testBuffer.txt");   // (optional) save a raw UYVY frame
#endif // ARF_SAVE_FRAME
    active_random_filter_header();                                  // Mostly printing and storing
    Rect crop 	        = setISPvars( width, height); 	            // Calculate ISP related parameters
    if(runCount < ARF_TIMEOUT) {
        runCount++;
        PRINT("Timeout %d\n", ARF_TIMEOUT - runCount);
        frameGrey.release();
        sourceFrame.release();
        return;
    }
    Mat sourceFrameCrop = sourceFrame(crop); 				                // Crop the frame
	trackObjects(sourceFrameCrop, frameGrey);                       // Track objects in sourceFrame
	eraseMemory();
	uint8_t r;
	for(r=0; r < trackRes_size; r++){                             // Convert angles & Write/Print output
#if ARF_OBJECT == ARF_BALL
		cam2body(&trackRes[r]);						                // Convert from camera angles to body angles (correct for roll)
		body2world(curEulerAngles, &trackRes[r]); 		                            // Convert from body angles to world coordinates (correct yaw and pitch)
#endif
		identifyObject(&trackRes[r]);                               // Identify the spotted neighbours
	}
#if ARF_MOD_VIDEO
	mod_video(sourceFrameCrop, frameGrey);                              // Modify the sourceframesourceFrame.cols-1
	//PRINT("gate_count:%d\n",gate_count);
	gate_count = 0;
#endif // ARF_MOD_VIDEO
#if ARF_CROSSHAIR
	//circle(sourceFrame,Point(ispHeight/2 - cropCol + crop.x, ispWidth/2), CFG_MT9F002_FISHEYE_RADIUS * ispScalar, cvScalar(0,255), 1);
	plotHorizon(sourceFrameCrop, curEulerAngles);
#endif
	frameGrey.release(); 			                                // Release Mat
	sourceFrameCrop.release();
	sourceFrame.release();                                          // Release Mat
	active_random_filter_footer( curEulerAngles );
	return;
}

Rect setISPvars( uint16_t width, uint16_t height){
    // This function computes the cropping according to the desires FOV Y and the current euler angles
    ARF_MIN_CIRCLE_SIZE   = pow(sqrt((double) default_calArea  * pow(ispScalar,2.0)) / ARF_CAM_RANGE, 2.0);
    ARF_MIN_LAYERS        = (uint16_t) round(ARF_MIN_CIRCLE_PERC * 2 * M_PI * sqrt((double) default_calArea  * pow(ispScalar,2.0)/ M_PI) / ARF_CAM_RANGE);
    ARF_LARGE_LAYERS      = (uint16_t) round(ARF_MIN_CIRCLE_PERC * 2 * M_PI * sqrt((double) default_calArea  * pow(ispScalar,2.0)/ M_PI) / 1.0);
    ARF_MIN_POINTS        = (uint16_t) round(0.25 * ARF_MIN_LAYERS);

    Rect crop;
    crop                   = cvRect(fillHeight, 0, width - fillHeight, height);
    return crop;
}

void trackObjects(Mat& sourceFrame, Mat& frameGrey){
    // Main function for tracking multiple objects on the frame
    pixCount    = 0;
    pixSucCount = 0;
    pixSrcCount = 0;
    pixNofCount = 0;
    pixDupCount = 0;
    if(ARF_FLOOD_STYLE != ARF_FLOOD_CW){
        processImage_omni(sourceFrame, frameGrey, ARF_RND_PIX_SAMPLE);
        processCrops(frameGrey);
    }
    else{
        processImage_cw(sourceFrame, frameGrey, ARF_RND_PIX_SAMPLE);
    }
#if ARF_BENCHMARK
    addBenchmark("image Thresholded");
#endif // ARF_BENCHMARK
    return;
}

void processCrops(Mat& frameGrey){
    vector<vector<Point> > contours;
    for(unsigned int r=0; r < cropAreas.size(); r++)
    {
        if(cropAreas[r].x != 0 && cropAreas[r].width != 0)
        {
            contours.clear();
#if ARF_MOD_VIDEO && ARF_DRAW_BOXES
            findContours(frameGrey(cropAreas[r]).clone(), contours, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE);
#else
            findContours(frameGrey(cropAreas[r]), contours, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE);
#endif //ARF_MOD_VIDEO && ARF_DRAW_BOXES
            for(unsigned int tc=0; tc < contours.size(); tc++)
            {
                if(contours[tc].size() > ARF_MIN_POINTS){
                    addContour(contours[tc], (uint16_t) cropAreas[r].x, (uint16_t) cropAreas[r].y);
                }
            }
        }
    }
#if ARF_BENCHMARK
    addBenchmark("Contours found");
#endif //ARF_BENCHMARK
}

void eraseMemory(void){
#ifdef __linux__
  pthread_mutex_lock(&neighbourMem_mutex);
#endif
    uint8_t index   = 0;
    uint8_t i       = 0;
    for(i=0; i < neighbourMem_size; i++){
        if((runCount - neighbourMem[i].lastSeen) <= ARF_MEMORY){
            if(i != index){
                neighbourMem[index] = neighbourMem[i];
            }
            index++;
        }
    }
    neighbourMem_size -= i - index;
#ifdef __linux__
  pthread_mutex_unlock(&neighbourMem_mutex);
#endif
}

#if ARF_OBJECT == ARF_BALL
void identifyObject(trackResults* trackRes){
    bool identified = FALSE;
    uint8_t neighbourID;
    uint8_t closestNeighbour = ARF_CAM_RANGE;
    for(unsigned int i=0; i < neighbourMem_size; i++)
    {
        double radius	= (runCount - neighbourMem[i].lastSeen) * 1.0 / ARF_FPS * ARF_VMAX;
        double dx 		= trackRes->x_w - neighbourMem[i].x_w;
        double dy       = trackRes->y_w - neighbourMem[i].y_w;
        double nRadius  = sqrt(pow(dx, 2.0) + pow(dy, 2.0));
        if(nRadius <= radius)
        {
            if(nRadius < closestNeighbour){
                neighbourID         = i;
                closestNeighbour    = nRadius;
                identified          = TRUE;
            }
        }
    }
    if(identified){
#ifdef __linux__
  pthread_mutex_lock(&neighbourMem_mutex);
#endif
        VERBOSE_PRINT("Identified object %d at (%4d, %4d)p (%5.2f, %5.2f, %5.2f)w\n", neighbourMem[neighbourID].id, trackRes->x_p, trackRes->y_p, trackRes->x_w, trackRes->y_w, trackRes->z_w);
        neighbourMem[neighbourID].lastSeen    = runCount;
        neighbourMem[neighbourID].x_w         = trackRes->x_w;
        neighbourMem[neighbourID].y_w         = trackRes->y_w;
        neighbourMem[neighbourID].z_w         = trackRes->z_w;
        neighbourMem[neighbourID].x_p         = trackRes->x_p;
        neighbourMem[neighbourID].y_p         = trackRes->y_p;
        neighbourMem[neighbourID].area_p      = trackRes->area_p;
        neighbourMem[neighbourID].r_c         = trackRes->r_c;
#ifdef __linux__
  pthread_mutex_unlock(&neighbourMem_mutex);
#endif
    }
    else{
        // We haven't identified
        VERBOSE_PRINT("New object at (%4d, %4d)p (%6.2f, %6.2f)c (%5.2f, %5.2f)b (%5.2f, %5.2f, %5.2f)w\n", trackRes->x_p, trackRes->y_p, trackRes->x_c * 180 / M_PI, trackRes->y_c * 180 / M_PI, trackRes->x_b, trackRes->y_b, trackRes->x_w, trackRes->y_w,  trackRes->z_w);
        memoryBlock curN;
        curN.lastSeen 	= runCount;
        curN.id 		= maxId;
        curN.x_w 		= trackRes->x_w;
        curN.y_w 		= trackRes->y_w;
        curN.z_w 		= trackRes->z_w;
        curN.x_p 		= trackRes->x_p;
        curN.y_p 		= trackRes->y_p;
        curN.area_p     = trackRes->area_p;
        curN.r_c        = trackRes->r_c;
        neighbourMem_add(curN);
        maxId++;
    }
    return;
}
#endif

#if ARF_OBJECT == ARF_GATE
void identifyObject(gateResults* gateRes){
#ifdef __linux__
  pthread_mutex_lock(&neighbourMem_mutex);
#endif
    bool identified = FALSE;
    for(unsigned int i=0; i < neighbourMem_size; i++)
    {
        double radius   = (runCount - neighbourMem[i].lastSeen) * 1.0 / ARF_FPS * ARF_VMAX;
        double dx       = gateRes->x_p - neighbourMem[i].x_p;
        double dy       = gateRes->y_p - neighbourMem[i].y_p;
        if(dx <= radius && dy <= radius && sqrt(pow(dx, 2.0) + pow(dy, 2.0)) <= radius)
        {
            VERBOSE_PRINT("Identified object %d at (%4d, %4d)p (%5.2f, %5.2f, %5.2f)w\n", neighbourMem[i].id, gateRes->x_p, gateRes->y_p, gateRes->x_w, gateRes->y_w, gateRes->z_w);
            neighbourMem[i].lastSeen    = runCount;
            neighbourMem[i].x_p         = gateRes->x_p;
            neighbourMem[i].y_p         = gateRes->y_p;
            neighbourMem[i].area_p      = gateRes->area_p;
            neighbourMem[i].corners[0]  = gateRes->corners[0];
            neighbourMem[i].corners[1]  = gateRes->corners[1];
            neighbourMem[i].corners[2]  = gateRes->corners[2];
            neighbourMem[i].corners[3]  = gateRes->corners[3];
            identified                  = TRUE;
        }
    }
    if(!identified){
        // We haven't identified
        VERBOSE_PRINT("New object at (%4d, %4d)p (%6.2f, %6.2f)c (%5.2f, %5.2f)b (%5.2f, %5.2f, %5.2f)w\n", gateRes->x_p, gateRes->y_p, gateRes->x_c * 180 / M_PI, gateRes->y_c * 180 / M_PI, gateRes->x_b, gateRes->y_b, gateRes->x_w, gateRes->y_w,  gateRes->z_w);
        memoryGateBlock curN;
        curN.lastSeen   = runCount;
        curN.id         = maxId;
        curN.x_p        = gateRes->x_p;
        curN.y_p        = gateRes->y_p;
        curN.area_p     = gateRes->area_p;
        curN.corners[0] = gateRes->corners[0];
        curN.corners[1] = gateRes->corners[1];
        curN.corners[2] = gateRes->corners[2];
        curN.corners[3] = gateRes->corners[3];
        if(neighbourMem_size == ARF_MAX_OBJECTS){
            neighbourMem[ARF_MAX_OBJECTS - 1] = curN;
            neighbourMem_lastId                     = neighbourMem_maxId;
        }
        else{
            neighbourMem[neighbourMem_size]         = curN;
            neighbourMem_lastId                     = neighbourMem_size;
            neighbourMem_size++;
        }
        maxId++;
    }
#ifdef __linux__
  pthread_mutex_unlock(&neighbourMem_mutex);
#endif
    return;
}
#endif

void estimatePosition(uint16_t xp, uint16_t yp, uint32_t area, double position[3]){
    // This function estimates the 3D position (in camera  coordinate system) according to pixel position
    // Calculate corrected calibration parameters
    double x_in, y_in, xAngle, yAngle;
    double calArea          = round(default_calArea * pow(ispScalar,2.0));
    double corArea          = area;// * pow(areaCorFrac, 2.0); // TODO: Fix better
    double dist             = sqrt((double) calArea) / sqrt(corArea);
    pixel2point((double) yp, (double) xp, &x_in, &y_in);
    point2angles(x_in, y_in, &xAngle, &yAngle);
    VERBOSE_PRINT("pixel(%d, %d) point(%0.2f, %0.2f) angles(%0.2f, %0.2f)\n",xp, yp, x_in, y_in, xAngle / M_PI * 180, yAngle / M_PI * 180);
    position[0]             = yAngle;
    position[1]             = xAngle;
    position[2]             = dist;
    return;
}

void cam2body(trackResults* trackRes){
    // Neighbour position returned in 2 angles and a radius.
    // x_c is the angle wrt vertical camera axis.   Defined clockwise/right positive
    // y_c is angle wrt camera horizon axis.        Defined upwards positive
    // r_c is radial distance in m.
    trackRes->x_b = trackRes->r_c * cos( -trackRes->y_c) * cos( trackRes->x_c ) + MT9F002_X_OFFSET;
    trackRes->y_b = trackRes->r_c * cos( -trackRes->y_c) * sin( trackRes->x_c ) + MT9F002_Y_OFFSET;
    trackRes->z_b = trackRes->r_c * sin( -trackRes->y_c) + MT9F002_Z_OFFSET;
    VERBOSE_PRINT("camera (%0.2f deg, %0.2f deg, %0.2f m) -> body (%0.2f m, %0.2f m, %0.2f m)\n", trackRes->x_c * 180/M_PI, trackRes->y_c * 180/M_PI, trackRes->r_c, trackRes->x_b, trackRes->y_b, trackRes->z_b);
    return;
}

void body2world( struct FloatEulers*  eulerAngles, trackResults* trackRes){
#if ARF_WORLDPOS
    pos     = stateGetPositionNed_f();      // Get your current position
#else
    struct NedCoor_f fakePos;
    fakePos.x       = 0.0;
    fakePos.y       = 0.0;
    fakePos.z       = 0.0;
    pos             = &fakePos;
#endif
#if ARF_NOYAW
    double psi      = 0.0;
#else
    double psi      = eulerAngles->psi;
#endif
    Matx33f rotZ(    cos(psi),                 -sin(psi),               0,
                     sin(psi),                  cos(psi),               0,
                     0,                         0,                      1);
    Matx31f bPos(trackRes->x_b, trackRes->y_b, trackRes->z_b);
    Matx31f wPos    = rotZ * bPos;
    trackRes->x_w   = wPos(0,0) + pos->x + ARF_OBJ_X_OFFSET;
    trackRes->y_w   = wPos(1,0) + pos->y + ARF_OBJ_Y_OFFSET;
#if ARF_USE_ALTITUDE
    trackRes->z_w   = wPos(2,0) + pos->z + ARF_OBJ_Z_OFFSET;
#else
    trackRes->z_w   = wPos(2,0) + ARF_OBJ_Z_OFFSET;
#endif
    VERBOSE_PRINT("body (%0.2f m, %0.2f m, %0.2f m) + pos(%0.2f m, %0.2f m, %0.2f m) + euler (%0.2f deg, %0.2f deg, %0.2f deg) -> world (%0.2f m, %0.2f m, %0.2f m)\n", trackRes->x_b, trackRes->y_b, trackRes->z_b, pos->x, pos->y, pos->z, eulerAngles->phi, eulerAngles->theta, psi, trackRes->x_w, trackRes->y_w, trackRes->z_w);
    return;
}

#if ARF_OBJECT == ARF_BALL
bool addContour(vector<Point> contour, uint16_t offsetX, uint16_t offsetY, double minDist, double maxDist){
    Moments m;
    if(ARF_FLOOD_STYLE == ARF_FLOOD_CW){
        VERBOSE_PRINT("Analyzing contour of length %d\n", objCont_size);
        m = objCont_moments();
    }
    else{
        VERBOSE_PRINT("Analyzing contour of length %d\n", contour.size());
        m = moments(contour);
    }
    double semi_major   = sqrt( 2 * ( m.mu20 + m.mu02 + sqrt( pow(m.mu20 - m.mu02, 2.0) + 4 * pow( m.mu11, 2.0 ) ) ) / m.m00 );
    double semi_minor   = sqrt( 2 * ( m.mu20 + m.mu02 - sqrt( pow(m.mu20 - m.mu02, 2.0) + 4 * pow( m.mu11, 2.0 ) ) ) / m.m00 );
    //double e            = sqrt( 1 - pow( semi_minor, 2.0 ) / pow( semi_major, 2.0 ));
    double e          = (pow(m.mu20 - m.mu02,2.0) - 4 * pow(m.mu11,2.0)) / pow(m.mu20 + m.mu02,2.0);
    double corArea      = M_PI * pow( semi_major, 2.0 );
    VERBOSE_PRINT("m00: %7.2f  cA: %7.2f  x: %5.2f  y: %5.2f  ecc: %4.2f  fill: %4.2f\n",m.m00, corArea, m.m10 / m.m00, m.m01 / m.m00, e, m.m00 / corArea);
    if(fabs(e) < ARF_MAX_CIRCLE_DEF)
    {
        if ( m.m00 / corArea > ARF_MIN_CIRCLE_PERC && corArea >= ARF_MIN_CIRCLE_SIZE)
        {
            VERBOSE_PRINT("ecc: %f   smax: %f    smin: %f    corArea: %f\n", e, semi_major, semi_minor, corArea);
            trackResults curRes;
            double position[3];
            curRes.x_p 		    = m.m10 / m.m00 + cropCol + offsetX;
            curRes.y_p 		    = m.m01 / m.m00 + offsetY;
            //curRes.area_p       = (uint32_t) m.m00;
            curRes.area_p 	    = (uint32_t) M_PI * pow( semi_major, 2.0 );
            estimatePosition(curRes.x_p, curRes.y_p, curRes.area_p, position);  // Estimate position in camera reference frame based on pixel location and area
            curRes.x_c 		    = position[0];
            curRes.y_c 		    = position[1];
            curRes.r_c 		    = position[2];
            if((!minDist || minDist <= curRes.r_c) && (!maxDist || maxDist >= curRes.r_c)){
                if(curRes.r_c <= ARF_CAM_RANGE){
                    uint8_t overwriteId = trackRes_size;    // Invalid index, so won't overwrite
                    for(uint8_t tr = 0; tr < trackRes_size; tr++){
                        if( sqrt( pow( (double) (curRes.x_p - trackRes[tr].x_p), 2.0) + pow( (double) (curRes.y_p - trackRes[tr].y_p),2.0)) < sqrt( max( curRes.area_p, trackRes[tr].area_p ) / M_PI ) ){
                            if(curRes.area_p < trackRes[tr].area_p){
                                // We are inside another contour
                                return false;
                            }
                            else{
                                // There is another contour inside the current contour
                                overwriteId = tr;               // Mark this result for overwriting
                            }
                        }
                    }
                    trackRes_add(curRes, overwriteId);     // Save results and push into trackRes
                    return true;
                }
            }
        }else if(ARF_SHOW_REJECT){
            VERBOSE_PRINT("Rejected. fill: %4.2f (min %4.2f), area %7.2f (min %7.2f)\n", m.m00 / corArea, ARF_MIN_CIRCLE_PERC, corArea, ARF_MIN_CIRCLE_SIZE);
        }
    }else if(ARF_SHOW_REJECT) {
        VERBOSE_PRINT("Rejected. eccentricity: %4.2f (max %4.2f)\n", e, ARF_MAX_CIRCLE_DEF);
    }
	return false;
}
#else
#if ARF_OBJECT == ARF_GATE
bool addContour(vector<Point> contour, uint16_t offsetX, uint16_t offsetY, double minDist, double maxDist){
    if(objCont_size < 150)
        return false;

    char text[200];
    uint8_t skipSize = 35;
    uint16_t cornerSize = 0;

    //PRINT("Contour has %d points\n", objCont_size);
    Point gate[32]; // We probably won't detect more than 32 corners per contour
    double x1,x2,x3,y1,y2,y3, angle2, angle3, dAngle, tmpAngle = 0, tmpX = 0, tmpY = 0;
    uint8_t corner = 0;
    for( int i = skipSize; i < objCont_size - skipSize; i++)
    {
        if(corner >= 32){
            break;
        }
        x1      = objCont_store[i - skipSize].x;
        y1      = objCont_store[i - skipSize].y;
        x2      = objCont_store[i].x;
        y2      = objCont_store[i].y;
        x3      = objCont_store[i + skipSize].x;
        y3      = objCont_store[i + skipSize].y;
        angle2  = atan2(y1 - y2, x1 - x2);
        angle3  = atan2(y2 - y3, x2 - x3);
        dAngle  = min(fabs(angle2 - angle3), 2 * M_PI - fabs(angle2 - angle3));
        //printf("%6.2f\n", dAngle * 180 / M_PI);
        if(dAngle > 35.f / 180.f * M_PI){
            // Found point!
            // Remember this point, and see if the next angle is bigger
            cornerSize++;
            if(dAngle > tmpAngle){
                tmpX        = x2;
                tmpY        = y2;
                tmpAngle    = dAngle;
            }
        }
        else if(tmpAngle > 0 && cornerSize > 0.5 * skipSize){
            gate[corner].x = tmpX;
            gate[corner].y = tmpY;
         //PRINT("Corner point at %0.0f %0.0f (ang1: %0.2f  ang2: %0.2f) (x1: %0.0f y1: %0.0f, x2: %0.0f y2: %0.0f, x3: %0.0f y3: %0.0f)\n", gate[corner].x, gate[corner].y, angle2 / M_PI * 180, angle3 / M_PI * 180, x1, y1, x2, y2, x3, y3);
            //circle(frameForPlotting, cvPoint(gate[corner].x,gate[corner].y), 5, cvScalar(100,255), 1);

            double px1, py1, angleX, angleY;
            pixel2point((double) gate[corner].y, (double) gate[corner].x + cropCol, &px1, &py1);
            point2angles(px1, py1, &angleY, &angleX);
          //PRINT("pixel(%0.0f, %0.0f) point(%0.2f, %0.2f) angles(%0.2f, %0.2f)\n",gate[corner].x, gate[corner].y + cropCol, px1, py1, angleX / M_PI * 180, angleY / M_PI * 180);
            tmpAngle = 0.0;
            cornerSize = 0;
            corner++;
        }
        else{
            cornerSize = 0;
        }
    }
    if(corner == 4 ){//|| corner == 3){
        /* Storing */
        Moments m = objCont_moments();
//         gateResults curGate;
//         curGate.x_p          = m.m10 / m.m00 + cropCol;
//         curGate.y_p          = m.m01 / m.m00;
//         curGate.area_p       = m.m00;
//         curGate.corners[0].x = gate[0].x;
//         curGate.corners[0].y = gate[0].y;
//         curGate.corners[1].x = gate[1].x;
//         curGate.corners[1].y = gate[1].y;
//         curGate.corners[2].x = gate[2].x;
//         curGate.corners[2].y = gate[2].y;
//         curGate.corners[3].x = gate[3].x;
//         curGate.corners[3].y = gate[3].y;
//         hold_Gate.x_p          = m.m10 / m.m00 + cropCol;
//         hold_Gate.y_p          = m.m01 / m.m00;
//         hold_Gate.area_p       = m.m00;
        hold_Gate[gate_count].corners[0].x = gate[0].x;
        hold_Gate[gate_count].corners[0].y = gate[0].y;
        hold_Gate[gate_count].corners[1].x = gate[1].x;
        hold_Gate[gate_count].corners[1].y = gate[1].y;
        hold_Gate[gate_count].corners[2].x = gate[2].x;
        hold_Gate[gate_count].corners[2].y = gate[2].y;
        hold_Gate[gate_count].corners[3].x = gate[3].x;
        hold_Gate[gate_count].corners[3].y = gate[3].y;
	gate_count +=1;
	
	//print angle of first corner point, should be top left?
	 double p_x, p_y, angle_X, angle_Y;
         pixel2point_cam((double) gate[0].y, (double) gate[0].x + cropCol, &p_x, &p_y);
         point2angles(p_x, p_y, &angle_Y, &angle_X);
	PRINT("angle_X:%f angle_Y:%f\n",angle_X*57.324,angle_Y*57.324);
        return true;
    }
    return false;
}
#endif
#endif

/** This function is adapted from the openCV source code, please refer to their documentation **/
static Moments objCont_moments( void ){
    Moments m;
    if( objCont_size == 0 )
        return m;

    double a00 = 0, a10 = 0, a01 = 0, a20 = 0, a11 = 0, a02 = 0, a30 = 0, a21 = 0, a12 = 0, a03 = 0;
    double xi, yi, xi2, yi2, xi_1, yi_1, xi_12, yi_12, dxy, xii_1, yii_1;

    xi_1 = objCont_store[objCont_size-1].x;
    yi_1 = objCont_store[objCont_size-1].y;

    xi_12 = xi_1 * xi_1;
    yi_12 = yi_1 * yi_1;

    for( int i = 0; i < objCont_size; i++ )
    {
        xi = objCont_store[i].x;
        yi = objCont_store[i].y;

        xi2 = xi * xi;
        yi2 = yi * yi;
        dxy = xi_1 * yi - xi * yi_1;
        xii_1 = xi_1 + xi;
        yii_1 = yi_1 + yi;

        a00 += dxy;
        a10 += dxy * xii_1;
        a01 += dxy * yii_1;
        a20 += dxy * (xi_1 * xii_1 + xi2);
        a11 += dxy * (xi_1 * (yii_1 + yi_1) + xi * (yii_1 + yi));
        a02 += dxy * (yi_1 * yii_1 + yi2);
        a30 += dxy * xii_1 * (xi_12 + xi2);
        a03 += dxy * yii_1 * (yi_12 + yi2);
        a21 += dxy * (xi_12 * (3 * yi_1 + yi) + 2 * xi * xi_1 * yii_1 +
                   xi2 * (yi_1 + 3 * yi));
        a12 += dxy * (yi_12 * (3 * xi_1 + xi) + 2 * yi * yi_1 * xii_1 +
                   yi2 * (xi_1 + 3 * xi));
        xi_1 = xi;
        yi_1 = yi;
        xi_12 = xi2;
        yi_12 = yi2;
    }

    if( fabs(a00) > 1.19209289550781250000e-7F)
    {
        double db1_2, db1_6, db1_12, db1_24, db1_20, db1_60;

        if( a00 > 0 )
        {
            db1_2 = 0.5;
            db1_6 = 0.16666666666666666666666666666667;
            db1_12 = 0.083333333333333333333333333333333;
            db1_24 = 0.041666666666666666666666666666667;
            db1_20 = 0.05;
            db1_60 = 0.016666666666666666666666666666667;
        }
        else
        {
            db1_2 = -0.5;
            db1_6 = -0.16666666666666666666666666666667;
            db1_12 = -0.083333333333333333333333333333333;
            db1_24 = -0.041666666666666666666666666666667;
            db1_20 = -0.05;
            db1_60 = -0.016666666666666666666666666666667;
        }

        // spatial moments
        m.m00 = a00 * db1_2;
        m.m10 = a10 * db1_6;
        m.m01 = a01 * db1_6;
        m.m20 = a20 * db1_12;
        m.m11 = a11 * db1_24;
        m.m02 = a02 * db1_12;
        m.m30 = a30 * db1_20;
        m.m21 = a21 * db1_60;
        m.m12 = a12 * db1_60;
        m.m03 = a03 * db1_20;

        double cx = 0, cy = 0;
        double mu20, mu11, mu02;
        double inv_m00 = 0.0;

        if( fabs(m.m00) > double(2.22044604925031308085e-16L))
        {
            inv_m00 = 1. / m.m00;
            cx = m.m10 * inv_m00;
            cy = m.m01 * inv_m00;
        }

        mu20 = m.m20 - m.m10 * cx;
        mu11 = m.m11 - m.m10 * cy;
        mu02 = m.m02 - m.m01 * cy;

        m.mu20 = mu20;
        m.mu11 = mu11;
        m.mu02 = mu02;

        m.mu30 = m.m30 - cx * (3 * mu20 + cx * m.m10);
        mu11 += mu11;
        m.mu21 = m.m21 - cx * (mu11 + cx * m.m01) - cy * mu20;
        m.mu12 = m.m12 - cy * (mu11 + cy * m.m10) - cx * mu02;
        m.mu03 = m.m03 - cy * (3 * mu02 + cy * m.m01);

        double inv_sqrt_m00 = std::sqrt(std::abs(inv_m00));
        double s2 = inv_m00*inv_m00, s3 = s2*inv_sqrt_m00;

        m.nu20 = m.mu20*s2; m.nu11 = m.mu11*s2; m.nu02 = m.mu02*s2;
        m.nu30 = m.mu30*s3; m.nu21 = m.mu21*s3; m.nu12 = m.mu12*s3; m.nu03 = m.mu03*s3;
    }
    return m;
}

void createSearchGrid(uint16_t x_p, uint16_t y_p, Point searchGrid[], uint8_t searchLayer, uint16_t sGridSize, int* maxRow, int* maxCol){
    if(searchLayer > 0){
        int16_t curY   = y_p + 1 * searchLayer * sGridSize;
        int16_t curX   = x_p - 1 * searchLayer * sGridSize;
        int8_t  dX,dY, jitterX, jitterY;
        float jitterInd;
        for(unsigned int s = 0; s < 4; s++){
            switch(s){
            case 0 :    dX =  sGridSize;    dY =  0;            break; // Right
            case 1 :    dX =  0;            dY = -sGridSize;    break; // Down
            case 2 :    dX = -sGridSize;    dY =  0;            break; // Left
            case 3 :    dX =  0;            dY =  sGridSize;    break; // Up
            }
            for(unsigned int i = 0; i < searchLayer * 2; i++){
                curY += dY;
                curX += dX;
                jitterInd = rand() / ((float) RAND_MAX);
                if(jitterInd < 0.1){
                    jitterX = 1;
                    jitterY = 0;
                }
                else if(jitterInd < 0.2){
                    jitterX = 0;
                    jitterY = 1;
                }
                else if(jitterInd < 0.3){
                    jitterX = -1;
                    jitterY =  0;
                }
                else if(jitterInd < 0.4){
                    jitterX =  0;
                    jitterY = -1;
                }
                else{
                    jitterX = 0;
                    jitterY = 0;
                }
                if((curY+jitterY) >= 0 && (curX+jitterX) >= 0 && (curY+jitterY) < *maxRow && (curX+jitterX) < *maxCol){
                    searchGrid[ s * 2 * searchLayer + i] = Point(curY + jitterY , curX + jitterX);
                }
                else{
                    searchGrid[ s * 2 * searchLayer + i] = Point(y_p, x_p); // Add duplicate point
                }
            }
        }
    }
    else{
        searchGrid[0] = Point( y_p, x_p );
    }
}

bool processImage_cw(Mat& sourceFrame, Mat& destFrame, uint16_t sampleSize){
    bool obj_detected   = false;
    bool foundObj       = false;
    objCont_size        = 0;
    gate_count = 0;
    if (sourceFrame.cols > 0 && sourceFrame.rows > 0)
    {
        if (ARF_SAMPLE_STYLE > 0){
#if ARF_OBJECT == ARF_BALL
            const uint8_t maxLayer  = 30;
            Point searchGrid[8*maxLayer];
            for(unsigned int rnm=0; rnm < neighbourMem_size; rnm++)
            {
                if( ( ( (int16_t) neighbourMem[rnm].x_p ) - cropCol ) >= 0 && ( ( (int16_t) neighbourMem[rnm].x_p ) - cropCol ) < sourceFrame.cols && neighbourMem[rnm].y_p < sourceFrame.rows){
                    VERBOSE_PRINT("Looking for object %d\n",neighbourMem[rnm].id);
                    uint8_t searchLayer     = 0;
                    uint8_t searchPoints    = 1;
                    uint16_t sGridSize      = 0.1 * sqrt(((float) neighbourMem[rnm].area_p) / M_PI);
                    foundObj                = false;
                    while(foundObj == false && searchLayer < maxLayer){
                        VERBOSE_PRINT("Searching layer %d\n",searchLayer);
                        createSearchGrid(neighbourMem[rnm].x_p - cropCol, neighbourMem[rnm].y_p, searchGrid, searchLayer, sGridSize, &sourceFrame.rows, &sourceFrame.cols);
                        for(uint8_t rsg = 0; rsg < searchPoints; rsg++)
                        {
                            //VERBOSE_PRINT("Searching (%d,%d)\n",searchGrid[rsg].x, searchGrid[rsg].y);
                            layerDepth              = 0;
                            //sourceFrame.at<Vec2b>(searchGrid[rsg].x, searchGrid[rsg].y)[0] = 0;
                            //sourceFrame.at<Vec2b>(searchGrid[rsg].x, searchGrid[rsg].y)[0] = 255;
                            if(pixFindContour_cw(sourceFrame, destFrame, searchGrid[rsg].x, searchGrid[rsg].y, ARF_SEARCH, true) == ARF_FINISHED){
                                double r_margin = (runCount - neighbourMem[rnm].lastSeen) * ARF_VMAX * 1 / ARF_FPS;
                                if(objCont_add(neighbourMem[rnm].r_c - r_margin, neighbourMem[rnm].r_c + r_margin)){
                                    VERBOSE_PRINT("Found object %d from (%d, %d) at (%d, %d) after %d layers\n",neighbourMem[rnm].id, neighbourMem[rnm].x_p, neighbourMem[rnm].y_p, trackRes[trackRes_lastId].x_p, trackRes[trackRes_lastId].y_p, searchLayer);
                                    foundObj                = true;
                                    obj_detected            = true;
                                    break;
                                }
                                else{
                                    VERBOSE_PRINT("Object too small/big to be object %d\n", neighbourMem[rnm].id);
                                }
                            }
                        }
                        searchLayer++;
                        searchPoints = 8 * searchLayer;
                    }
                    if(!foundObj){
                        objCont_size = 0;
                        VERBOSE_PRINT("Could not find object %d\n",neighbourMem[rnm].id);
                    }
                }
                else{
                    VERBOSE_PRINT("Object %d is not within valid bounds (%d, %d) [0-%d][0-%d]\n",neighbourMem[rnm].id, ( (int16_t) neighbourMem[rnm].x_p ) - cropCol, neighbourMem[rnm].y_p, sourceFrame.cols, sourceFrame.rows);
                }
            }
#else
#if ARF_OBJECT == ARF_GATE
            const uint8_t maxLayer  = 30;
            Point searchGrid[8*maxLayer];
            for(unsigned int rnm=0; rnm < neighbourMem_size; rnm++)
            {
                if( ( ( (int16_t) neighbourMem[rnm].x_p ) - cropCol ) >= 0 && ( ( (int16_t) neighbourMem[rnm].x_p ) - cropCol ) < sourceFrame.cols && neighbourMem[rnm].y_p < sourceFrame.rows){
                    VERBOSE_PRINT("Looking for object %d\n",neighbourMem[rnm].id);
                    uint8_t searchLayer     = 0;
                    uint8_t searchPoints    = 1;
                    uint16_t sGridSize      = 0.1 * sqrt(((float) neighbourMem[rnm].area_p) / M_PI);
                    foundObj                = false;
                    while(foundObj == false && searchLayer < maxLayer){
                        VERBOSE_PRINT("Searching layer %d\n",searchLayer);
                        for(uint8_t corner = 0; corner < 4; corner++){
                            createSearchGrid(neighbourMem[rnm].corners[corner].x - cropCol, neighbourMem[rnm].corners[corner].y, searchGrid, searchLayer, sGridSize, &sourceFrame.rows, &sourceFrame.cols);
                            for(uint8_t rsg = 0; rsg < searchPoints; rsg++)
                            {
                                //VERBOSE_PRINT("Searching (%d,%d)\n",searchGrid[rsg].x, searchGrid[rsg].y);
                                layerDepth              = 0;
                                //sourceFrame.at<Vec2b>(searchGrid[rsg].x, searchGrid[rsg].y)[0] = 0;
                                //sourceFrame.at<Vec2b>(searchGrid[rsg].x, searchGrid[rsg].y)[0] = 255;
                                if(pixFindContour_cw(sourceFrame, destFrame, searchGrid[rsg].x, searchGrid[rsg].y, ARF_SEARCH, true) == ARF_FINISHED){
                                    if(objCont_add(0, 0)){
                                        VERBOSE_PRINT("Found object %d from (%d, %d) at (%d, %d) after %d layers\n",neighbourMem[rnm].id, neighbourMem[rnm].x_p, neighbourMem[rnm].y_p, trackRes[trackRes_lastId].x_p, trackRes[trackRes_lastId].y_p, searchLayer);
                                        foundObj                = true;
                                        obj_detected            = true;
                                        break;
                                    }
                                }
                            }
                            if(foundObj){
                                break;
                            }
                        }
                        searchLayer++;
                        searchPoints = 8 * searchLayer;
                    }
                    if(!foundObj){
                        objCont_size = 0;
                        VERBOSE_PRINT("Could not find object %d\n",neighbourMem[rnm].id);
                    }
                }
                else{
                    VERBOSE_PRINT("Object %d is not within valid bounds (%d, %d) [0-%d][0-%d]\n",neighbourMem[rnm].id, ( (int16_t) neighbourMem[rnm].x_p ) - cropCol, neighbourMem[rnm].y_p, sourceFrame.cols, sourceFrame.rows);
                }
            }
#endif
#endif
        }
        switch(ARF_SAMPLE_STYLE){
        case ARF_STYLE_FULL : {
            for(int r = 0; r < sourceFrame.rows; r++)
            {
                for(int c= 0; c < sourceFrame.cols; c++)
                {
                    layerDepth              = 0;
                    if(pixFindContour_cw(sourceFrame, destFrame, r, c, ARF_UP, false) == ARF_FINISHED)
                    {
                        obj_detected            = true;
                    }
                }
            }
            Rect fullCrop;
            fullCrop.x 		= 1;
            fullCrop.y 		= 0;
            fullCrop.width 	= sourceFrame.cols-1;
            fullCrop.height = sourceFrame.rows;
            cropAreas.push_back(fullCrop);
            break;
        }
        case ARF_STYLE_GRID : {
            int spacing     = (int) sqrt((sourceFrame.rows * sourceFrame.cols) / sampleSize);
            for(int r = spacing; r < sourceFrame.rows; r+=spacing)
            {
                for(int c=spacing; c < sourceFrame.cols; c+=spacing)
                {
                    sample++;
                    layerDepth      = 0;
                    //objCont.clear();
                    objCont_size    = 0;
                    if(pixFindContour_cw(sourceFrame, destFrame, r, c, ARF_SEARCH, true) == ARF_FINISHED)
                    {
                        if(objCont_add()){
                            obj_detected = true;
                        }
                    }
                }
            }
            break;
        }
        case ARF_STYLE_RANDOM : {
            int rndRow, rndCol;
            for(int i = 0; i<sampleSize; i++)
            {
                layerDepth  = 0;
                sample++;
                rndRow      = (int) round(((double) rand())/((double) RAND_MAX)*(sourceFrame.rows-1));
                rndCol      = (int) round(((double) rand())/((double) RAND_MAX)*(sourceFrame.cols-1));
                objCont_size    = 0;
                bool too_close  = false;
                for(unsigned int tr = 0; tr < trackRes_size; tr++){
                    if(sqrt(pow((double) (rndCol - (trackRes[tr].x_p - cropCol)),2.0)+pow((double) (rndRow - trackRes[tr].y_p),2.0)) <= 1.15 * sqrt(trackRes[tr].area_p / M_PI)){
                        too_close = true;
                        break;
                    }

                }
                if(too_close){
                    continue;
                }
                if(pixFindContour_cw(sourceFrame, destFrame, rndRow, rndCol, ARF_SEARCH, true) == ARF_FINISHED)
                {
                    if(objCont_add()){
                        obj_detected = true;
                    }
                }
            }
            break;
        }
        }
    }
    return obj_detected;
}

bool processImage_omni(Mat& sourceFrame, Mat& destFrame, uint16_t sampleSize){
    bool obj_detected   = false;
    bool foundObj       = false;
    cropAreas.clear();
    if (sourceFrame.cols > 0 && sourceFrame.rows > 0)
    {
        if(ARF_SAMPLE_STYLE > 0){
            const uint8_t maxLayer  = 5;
            Point searchGrid[8*maxLayer];
            for(unsigned int rnm=0; rnm < neighbourMem_size; rnm++)
            {
                objCrop.x               = neighbourMem[rnm].x_p - cropCol;
                objCrop.y               = neighbourMem[rnm].y_p;
                objCrop.width           = 0;
                objCrop.height          = 0;
                uint8_t searchLayer     = 0;
                uint8_t searchPoints    = 1;
                uint16_t sGridSize      =  0.5 * sqrt(((float) neighbourMem[rnm].area_p) / M_PI);
                foundObj                = false; // We're pessimistic that we can find the same object
                while(!foundObj && searchLayer < maxLayer){
                    createSearchGrid(neighbourMem[rnm].x_p - cropCol, neighbourMem[rnm].y_p, searchGrid, searchLayer, sGridSize, &sourceFrame.rows, &sourceFrame.cols);
                    for(uint8_t rsg=0; rsg < searchPoints; rsg++)
                    {
                        layerDepth              = 0;
                        if(pixFindContour_omni(sourceFrame, destFrame, searchGrid[rsg].x, searchGrid[rsg].y, ARF_SEARCH, true) == ARF_FINISHED){
                            if(layerDepth > ARF_MIN_LAYERS){
                                objCrop                 = enlargeRectangle(sourceFrame, objCrop, ARF_CROP_X);
                                addCrop();
                                foundObj                = true;
                                obj_detected            = true;
                                break;
                            }
                        }
                    }
                    searchLayer++;
                    searchPoints = 8 * searchLayer;
                }
                if(!foundObj){
                    VERBOSE_PRINT("Could not find object %d\n",neighbourMem[rnm].id);
                }
            }
        }
        switch(ARF_SAMPLE_STYLE){
        case ARF_STYLE_FULL : {
            for(int r = 0; r < sourceFrame.rows; r++)
            {
                for(int c= 0; c < sourceFrame.cols; c++)
                {
                    layerDepth              = 0;
                    if(pixFindContour_omni(sourceFrame, destFrame, r, c, ARF_UP, false) == ARF_FINISHED)
                    {
                        obj_detected            = true;
                    }
                }
            }
            Rect fullCrop;
            fullCrop.x      = 1;
            fullCrop.y      = 0;
            fullCrop.width  = sourceFrame.cols-1;
            fullCrop.height = sourceFrame.rows;
            cropAreas.push_back(fullCrop);
            break;
        }
        case ARF_STYLE_GRID : {
            int spacing     = (int) sqrt((sourceFrame.rows * sourceFrame.cols) / sampleSize);
            for(int r = spacing; r < sourceFrame.rows; r+=spacing)
            {
                for(int c=spacing; c < sourceFrame.cols; c+=spacing)
                {
                    sample++;
                    layerDepth      = 0;
                    objCrop.x       = c;
                    objCrop.y       = r;
                    objCrop.width   = 0;
                    objCrop.height  = 0;
                    if(pixFindContour_omni(sourceFrame, destFrame, r, c, ARF_SEARCH, true) == ARF_FINISHED)
                    {
                        if(layerDepth > ARF_MIN_LAYERS){
                            objCrop         = enlargeRectangle(sourceFrame, objCrop, ARF_CROP_X);
                            obj_detected    = true;
                            addCrop();
                        }
                    }
                }
            }
            break;
        }
        case ARF_STYLE_RANDOM : {
            int rndRow, rndCol;
            for(int i = 0; i<sampleSize; i++)
            {
                layerDepth  = 0;
                sample++;
                rndRow      = (int) round(((double) rand())/((double) RAND_MAX)*(sourceFrame.rows-1));
                rndCol      = (int) round(((double) rand())/((double) RAND_MAX)*(sourceFrame.cols-1));
                objCrop.x       = rndCol;
                objCrop.y       = rndRow;
                objCrop.width   = 0;
                objCrop.height  = 0;
                bool too_close  = false;
                for(unsigned int tr = 0; tr < trackRes_size; tr++){
                    if(sqrt(pow((double) (rndCol - (trackRes[tr].x_p - cropCol)),2.0)+pow((double) (rndRow - trackRes[tr].y_p),2.0)) <= 1.15 * sqrt(trackRes[tr].area_p / M_PI)){
                        too_close = true;
                        break;
                    }

                }
                if(too_close){
                    continue;
                }
                if(pixFindContour_omni(sourceFrame, destFrame, rndRow, rndCol, ARF_SEARCH, true) == ARF_FINISHED)
                {
                    if(layerDepth > ARF_MIN_LAYERS){
                        objCrop         = enlargeRectangle(sourceFrame, objCrop, ARF_CROP_X);
                        obj_detected    = true;
                        addCrop();
                    }
                }
            }
            break;
        }
        }
    }
    return obj_detected;
}

bool objCont_add( double minDist, double maxDist){
    if(layerDepth > ARF_MIN_LAYERS && objCont_size > ARF_MIN_POINTS){
        vector<Point> objCont;
        if(addContour(objCont, (uint16_t) 0, (uint16_t) 0, minDist, maxDist)){
            return true;
        }
    }
    objCont_size = 0;
    return false;
}

void objCont_addPoint(uint16_t* row, uint16_t* col){
    objCont_store[objCont_size] = Point(*col,*row);
    objCont_size++;
}

bool pixTest(uint8_t *Y, uint8_t *U, uint8_t *V, uint8_t *prevDir){
    if(*V > (*U + ARF_GREY_THRES) && *Y >= ARF_Y_MIN && *Y <= ARF_Y_MAX && *U >= ARF_U_MIN && *U <= ARF_U_MAX && *V >= ARF_V_MIN && *V <= ARF_V_MAX){
        return true;
    }
    else{
        //if(*prevDir != ARF_SEARCH){
            //if(abs(cmpY - *Y) <= ARF_CDIST_YTHRES && abs(cmpU - *U) <= ARF_CDIST_UTHRES && abs(cmpV - *V) <= ARF_CDIST_VTHRES){
                //PRINT("(cmpY: %d cmpU: %d cmpV: %d) (Y: %d U: %d V: %d) (dY: %d  dU: %d  dV: %d)\n", cmpY, cmpU, cmpV, *Y, *U, *V,(cmpY - *Y),(cmpU - *U),(cmpV - *V));
            //    return true;
            //}
            //else{
                return false;
            //}
        //}
        //else{
        //    return false;
        //}
    }//
}

int pixFindContour_cw(Mat& sourceFrame, Mat& destFrame, uint16_t row, uint16_t col, uint8_t prevDir, bool cascade){
    layerDepth++;
    pixCount++;
    if(destFrame.at<uint8_t>(row, col) >= 75){
            pixDupCount++;
            return ARF_DUPLICATE;
    }
    uint8_t U, Y, V;
    getYUVColours(sourceFrame, row, col, &Y, &U, &V);
    if(pixTest(&Y, &U, &V, &prevDir)){
        destFrame.at<uint8_t>(row, col) = 75;
        if(cascade){
            uint8_t nextDirCnt, nextDir[6];
            uint16_t newRow, newCol;
            bool success = false;
            uint8_t d = 0, edge = 0;
            getNextDirection_cw(prevDir, nextDir, &nextDirCnt);
            while(layerDepth < ARF_MAX_LAYERS && d < nextDirCnt && success == false){
                //cmpY                = Y;
                //cmpU                = U;
                //cmpV                = V;
                newRow              = row;
                newCol              = col;
                if(getNewPosition(nextDir[d], &newRow, &newCol, &sourceFrame.rows, &sourceFrame.cols)){
                    int result;
                    if(d == 0){
                        result = pixFindContour_cw(sourceFrame, destFrame, newRow, newCol, nextDir[d], true);
                        double depth = 0;
                        while( depth < ARF_MAX_SEARCH_PIXEL_SKIP && result == ARF_NO_FOUND && getNewPosition(nextDir[d], &newRow, &newCol, &sourceFrame.rows, &sourceFrame.cols) ){
                            result = pixFindContour_cw(sourceFrame, destFrame, newRow, newCol, nextDir[d], true);
                            ///< Check up another pixel to account for noise
                            depth++;
                        }
                    }
                    else{
                        destFrame.at<uint8_t>(row, col) = 76;
                        objCont_sRow                    = row;
                        objCont_sCol                    = col;
                        objCont_size                    = 0;
                        result                          = pixFollowContour_cw(sourceFrame, destFrame, newRow, newCol, nextDir[d]);
                    }
                    switch(result){ // Catch the proper response for the tested pixel
                    case ARF_FINISHED : {
                        pixSrcCount++;
                        return ARF_FINISHED;
                        break;
                    }
                    case ARF_NO_FOUND : {
                        edge++;
                        break;
                    }
                    case ARF_DUPLICATE : {
                        pixDupCount++;
                        //destFrame.at<uint8_t>(row, col) = 0;
                        return ARF_DUPLICATE;
                        break;
                    }
                    case ARF_ERROR : {
                        if(layerDepth > ARF_MAX_LAYERS){
                            //destFrame.at<uint8_t>(row, col) = 0;
                            return ARF_FINISHED;
                        }
                        else{
                            edge++;
                        }
                        break;
                    }
                    default : {
                        //destFrame.at<uint8_t>(row, col) = 0;
                        return ARF_ERROR;
                    }
                    }
                }
                else{
                    edge++;
                }
                d++;
            }
            pixNofCount++;
            //destFrame.at<uint8_t>(row, col) = 0;
            return ARF_NO_FOUND; // Dead end
        }
        else{
            pixSucCount++;
            return ARF_FINISHED;
        }
    }
    else{
        pixNofCount++;
        return ARF_NO_FOUND;
    }
}

int pixFollowContour_cw(Mat& sourceFrame, Mat& destFrame, uint16_t row, uint16_t col, uint8_t prevDir){
    layerDepth++;
    pixCount++;
    if(destFrame.at<uint8_t>(row, col) == 76  // Arrived neatly back at starting pos
            || (destFrame.at<uint8_t>(row, col) == 250 && layerDepth > 5 * ARF_MIN_LAYERS)){
            //|| ( layerDepth > ARF_LARGE_LAYERS && abs(objCont_sRow - row) < round( layerDepth * ARF_LARGE_SKIP_FACTOR ) && abs(objCont_sCol - col) < round( layerDepth * ARF_LARGE_SKIP_FACTOR ) )){  // Close enough for large contours
        // This is my starting position, finished!
        if(layerDepth > ARF_MIN_LAYERS){
            destFrame.at<uint8_t>(row, col) = 255;
            objCont_addPoint(&row,&col);
            VERBOSE_PRINT("ARF_FINISHED back at (%d, %d) after %d pixels\n",row, col, layerDepth);
            pixSucCount++;
            return ARF_FINISHED;
        }
        else{
            VERBOSE_PRINT("ARF_NO_FOUND back at (%d, %d) after only %d pixels\n",row, col, layerDepth);
            pixNofCount++;
            objCont_size = 0;
            return ARF_NO_FOUND; // TODO: Should this be ARF_NO_FOUND?
        }

    }
    if(destFrame.at<uint8_t>(row, col) >= 75){
        // I've already searched this pixel before
        pixDupCount++;
        return ARF_DUPLICATE;
    }
    uint8_t U, Y, V;
    getYUVColours(sourceFrame, row, col, &Y, &U, &V);
    if(pixTest(&Y, &U, &V, &prevDir)){
        if(layerDepth < 4 * ARF_MIN_LAYERS){
            destFrame.at<uint8_t>(row, col) = 250;
        }
        else{
            destFrame.at<uint8_t>(row, col) = 255;
        }
        uint8_t nextDirCnt, nextDir[6];
        uint16_t newRow, newCol;
        bool success = false;
        uint8_t d = 0, edge = 0;
        getNextDirection_cw(prevDir, nextDir, &nextDirCnt);
        while(layerDepth < ARF_MAX_LAYERS && d < nextDirCnt && success == false){
            newRow              = row;
            newCol              = col;
            if(getNewPosition(nextDir[d], &newRow, &newCol, &sourceFrame.rows, &sourceFrame.cols)){
                switch(pixFollowContour_cw(sourceFrame, destFrame, newRow, newCol, nextDir[d])){ // Catch the proper response for the tested pixel
                case ARF_FINISHED : {
#if ARF_OBJECT == ARF_BALL
                    if(prevDir != nextDir[d]){
                        objCont_addPoint(&row,&col);
                    }
#else
                    objCont_addPoint(&row,&col);
#endif
#if ARF_MARK_CONTOURS
                    sourceFrame.at<Vec2b>(row, col)[0] = 0;
                    sourceFrame.at<Vec2b>(row, col)[1] = 127;
#endif
                    pixSucCount++;
                    return ARF_FINISHED;
                    break;
                }
                case ARF_NO_FOUND : {
                    edge++;
                    break;
                }
                case ARF_DUPLICATE : {
                    pixDupCount++;
                    destFrame.at<uint8_t>(row, col) = 0;
                    return ARF_DUPLICATE;
                    break;
                }
                case ARF_ERROR : {
                    if(layerDepth > ARF_MAX_LAYERS){
                        destFrame.at<uint8_t>(row, col) = 0;
                        return ARF_FINISHED;
                    }
                    else{
                        edge++;
                    }
                    break;
                }
                default : {
                    destFrame.at<uint8_t>(row, col) = 0;
                    return ARF_ERROR;
                }
                }
            }
            else{
                edge++;
            }
            d++;
        }
        pixNofCount++;
        destFrame.at<uint8_t>(row, col) = 0;
        return ARF_NO_FOUND; // Dead end
    }
    else{
        pixNofCount++;
        return ARF_NO_FOUND;
    }
}

int pixFindContour_omni(Mat& sourceFrame, Mat& destFrame, uint16_t row, uint16_t col, uint8_t prevDir, bool cascade){
    layerDepth++;
    pixCount++;
    if(prevDir == ARF_SEARCH){
        prevDir = ARF_UP;
    }
    if(destFrame.at<uint8_t>(row, col) == 255){
        pixDupCount++;
        return ARF_DUPLICATE;
    }
    uint8_t U, Y, V;
    getYUVColours(sourceFrame, row, col, &Y, &U, &V);
    if(pixTest(&Y, &U, &V, &prevDir))
    {
        if(prevDir != ARF_SEARCH)
        {
            destFrame.at<uint8_t>(row, col) = 255;
        }
        if(cascade)
        {
            pixSucCount++;
            uint8_t nextDir[3];
            uint8_t nextDirCnt  = 3;
            getNextDirection_omni(prevDir, nextDir, &nextDirCnt);
            bool success        = false;
            uint16_t d          = 0;
            uint16_t newRow, newCol;
            int res[4]          = {ARF_NO_FOUND,ARF_NO_FOUND,ARF_NO_FOUND,ARF_NO_FOUND};
            while(d < nextDirCnt && success == false)
            {
                newRow              = row;
                newCol              = col;
                if(getNewPosition(nextDir[d], &newRow, &newCol, &sourceFrame.rows, &sourceFrame.cols)){
                    res[d]              = pixFindContour_omni(sourceFrame, destFrame, newRow, newCol, nextDir[d], true);
                }
                else{
                    res[d]              = ARF_NO_FOUND;
                }
                d++;
            }
            if(res[0] <= ARF_NO_FOUND && res[1] <= ARF_NO_FOUND && res[2] <= ARF_NO_FOUND && res[3] <= ARF_NO_FOUND)
            {
                return ARF_FINISHED;
            }
            objCrop.width       = max<uint16_t>(objCrop.x + objCrop.width,  col) - min<uint16_t>(objCrop.x, col);
            objCrop.height      = max<uint16_t>(objCrop.y + objCrop.height, row) - min<uint16_t>(objCrop.y, row);
            objCrop.x           = min<uint16_t>(objCrop.x, col);
            objCrop.y           = min<uint16_t>(objCrop.y, row);
            return ARF_FINISHED;
        }else{
            pixSucCount++;
            return ARF_FINISHED;
        }
    }else{
        pixNofCount++;
        return ARF_NO_FOUND;
    }
}

void getYUVColours(Mat& sourceFrame, uint16_t row, uint16_t col, uint8_t* Y, uint8_t* U, uint8_t* V){
    if (((col & 1) == 0 && (cropCol & 1) == 0) || ((col & 1) == 1 && (cropCol & 1) == 1))
    {
        // Even col number
        *U = sourceFrame.at<Vec2b>(row, col)[0]; // U1
        *Y = sourceFrame.at<Vec2b>(row, col)[1]; // Y1
        if(col + 1 < sourceFrame.cols)
        {
            *V = sourceFrame.at<Vec2b>(row, col + 1)[0]; // V2
        }else{
            *V = sourceFrame.at<Vec2b>(row, col - 1)[0]; // V2
        }
    }else{
        // Uneven col number
        *V = sourceFrame.at<Vec2b>(row, col)[0]; // V2
        *Y  = sourceFrame.at<Vec2b>(row, col)[1]; // Y2
        if(col > 0)
        {
            *U = sourceFrame.at<Vec2b>(row, col - 1)[0]; // U1
        }else{
            *U = sourceFrame.at<Vec2b>(row, col + 1)[0]; // U1
        }
    }
}

void getNextDirection_cw(uint8_t prevDir, uint8_t* nextDir, uint8_t* nextDirCnt){
	*nextDirCnt = 5;
	switch(prevDir) // Find out which directions to try next
	{
	default:
	case ARF_SEARCH :
	    nextDir[0] = ARF_SEARCH;
	    nextDir[1] = ARF_UP_RIGHT;
	    nextDir[2] = ARF_RIGHT;
	    nextDir[3] = ARF_RIGHT_DOWN;
	    *nextDirCnt = 4;
	    break;
	case ARF_UP :
	    //nextDir[0] = ARF_LEFT;
	    nextDir[0] = ARF_LEFT_UP;
	    nextDir[1] = ARF_UP;
	    nextDir[2] = ARF_UP_RIGHT;
	    nextDir[3] = ARF_RIGHT;
	    nextDir[4] = ARF_RIGHT_DOWN;
	    break;
	case ARF_UP_RIGHT :
	    nextDir[0] = ARF_LEFT_UP;
	    nextDir[1] = ARF_UP;
	    nextDir[2] = ARF_UP_RIGHT;
	    nextDir[3] = ARF_RIGHT;
	    nextDir[4] = ARF_RIGHT_DOWN;
	    nextDir[5] = ARF_DOWN;
	    *nextDirCnt = 6;
	    break;
	case ARF_RIGHT :
	    //nextDir[0] = ARF_UP;
	    nextDir[0] = ARF_UP_RIGHT;
	    nextDir[1] = ARF_RIGHT;
	    nextDir[2] = ARF_RIGHT_DOWN;
	    nextDir[3] = ARF_DOWN;
	    nextDir[4] = ARF_DOWN_LEFT;
	    break;
	case ARF_RIGHT_DOWN :
	    nextDir[0] = ARF_UP_RIGHT;
	    nextDir[1] = ARF_RIGHT;
	    nextDir[2] = ARF_RIGHT_DOWN;
	    nextDir[3] = ARF_DOWN;
	    nextDir[4] = ARF_DOWN_LEFT;
	    nextDir[5] = ARF_LEFT;
	    *nextDirCnt = 6;
	    break;
	case ARF_DOWN :
	    //nextDir[0] = ARF_RIGHT;
	    nextDir[0] = ARF_RIGHT_DOWN;
	    nextDir[1] = ARF_DOWN;
	    nextDir[2] = ARF_DOWN_LEFT;
	    nextDir[3] = ARF_LEFT;
	    nextDir[4] = ARF_LEFT_UP;
	    break;
	case ARF_DOWN_LEFT :
	    nextDir[0] = ARF_RIGHT_DOWN;
	    nextDir[1] = ARF_DOWN;
	    nextDir[2] = ARF_DOWN_LEFT;
	    nextDir[3] = ARF_LEFT;
	    nextDir[4] = ARF_LEFT_UP;
	    nextDir[5] = ARF_UP;
	    *nextDirCnt = 6;
	    break;
	case ARF_LEFT :
	    //nextDir[0] = ARF_DOWN;
	    nextDir[0] = ARF_DOWN_LEFT;
	    nextDir[1] = ARF_LEFT;
	    nextDir[2] = ARF_LEFT_UP;
	    nextDir[3] = ARF_UP;
	    nextDir[4] = ARF_UP_RIGHT;
	    break;
	case ARF_LEFT_UP :
	    nextDir[0] = ARF_DOWN_LEFT;
	    nextDir[1] = ARF_LEFT;
	    nextDir[2] = ARF_LEFT_UP;
	    nextDir[3] = ARF_UP;
	    nextDir[4] = ARF_UP_RIGHT;
	    nextDir[5] = ARF_RIGHT;
	    *nextDirCnt = 6;
	    break;
	}
	return;
}

void getNextDirection_omni(uint8_t prevDir, uint8_t* nextDir, uint8_t* nextDirCnt){
   switch(prevDir)
    {
    case ARF_UP :
        nextDir[0] = ARF_LEFT;
        nextDir[1] = ARF_UP;
        nextDir[2] = ARF_RIGHT;
        *nextDirCnt = 3;
        break;
    case ARF_RIGHT :
        nextDir[0] = ARF_UP;
        nextDir[1] = ARF_RIGHT;
        nextDir[2] = ARF_DOWN;
        *nextDirCnt = 3;
        break;
    case ARF_DOWN :
        nextDir[0] = ARF_RIGHT;
        nextDir[1] = ARF_DOWN;
        nextDir[2] = ARF_LEFT;
        *nextDirCnt = 3;
        break;
    case ARF_LEFT :
        nextDir[0] = ARF_DOWN;
        nextDir[1] = ARF_LEFT;
        nextDir[2] = ARF_UP;
        *nextDirCnt = 3;
        break;
    }
    return;
}

bool getNewPosition(uint8_t nextDir, uint16_t* newRow, uint16_t* newCol, int* maxRow, int* maxCol){
	switch(nextDir) // Set the location of the next pixel to test
	{
	case ARF_SEARCH :
	    if(*newRow > 0){
	        *newRow += -1;
	    }
	    else{
	        return false;
	    }
	    break;
	case ARF_UP :
	    if(*newRow > 0){
	        *newRow += -1;
	    }
	    else{
            return false;
        }
	    break;
	case ARF_UP_RIGHT :
	    if(*newRow > 0){
	        *newRow += -1;
	    }
	    else{
	        return false;
	    }
	    if(*newCol < *maxCol - 1)
	    {
	        *newCol += 1;
	    }
	    else{
	        return false;
	    }
	    break;
	case ARF_RIGHT :
	    if(*newCol < *maxCol - 1){
	        *newCol += 1;
	    }
	    else{
	        return false;
	    }
	    break;
	case ARF_RIGHT_DOWN :
	    if(*newRow < *maxRow - 1){
	        *newRow += 1;
	    }
	    else{
	        return false;
	    }
	    if(*newCol < *maxCol - 1){
	        *newCol += 1;
	    }
	    else{
	        return false;
	    }
	    break;
	case ARF_DOWN :
	    if(*newRow < *maxRow - 1){
	        *newRow += 1;
	    }
	    else{
	        return false;
	    }
		break;
	case ARF_DOWN_LEFT :
	    if(*newRow < *maxRow - 1){
	        *newRow += 1;
	    }
	    else{
	        return false;
	    }
		if(*newCol > 0){
		    *newCol += -1;
		}
		else{
		    return false;
		}
		break;
	case ARF_LEFT :
	    if(*newCol > 0){
	        *newCol += -1;
	    }
	    else{
	        return false;
	    }
	    break;
	case ARF_LEFT_UP :
	    if(*newRow > 0){
	        *newRow += -1;
	    }
	    else{
	        return false;
	    }
	    if(*newCol > 0){
	        *newCol += -1;
	    }
	    else{
	        return false;
	    }
	    break;
	default:
	    VERBOSE_PRINT("[AR_FILTER-ERR] Invalid next-dir: %i\n", nextDir);
		return false;
		break;
	}
	return true;
}



//line check
int line_check(Mat& checkFrame,Point P1,Point P2,float quality_threshhold, int min_sides)
{
uint8_t y_c = 0;
uint8_t u_c = 0;
uint8_t v_c = 0;

//getYUVColours(checkFrame, checkFrame.rows / 2.0, checkFrame.cols / 2.0, &y_c, &u_c, &v_c);
//PRINT("y_c:%d u_c:%d v_c:%d\n",y_c,u_c,v_c);

  float n_points = 0;
  float n_colored_points = 0;

  float t_step = 0.05;
  float x, y;
  float t;
//   PRINT("P1.x:%d P1.y:%d\n",P1.x,P1.y);
//   PRINT("P2.x:%d P2.y:%d\n",P2.x,P2.y);
  // go from Q1 to Q2 in 1/t_step steps:
  for (t = 0.0f; t < 1.0f; t += t_step) {
    // determine integer coordinate on the line:
    float x1_temp = (float)P1.x;
    float x2_temp = (float)P2.x;
    float y1_temp = (float)P1.y;
    float y2_temp = (float)P2.y;
    
    x = (t * x1_temp + (1.0f - t) * x2_temp);
    y = (t * y1_temp + (1.0f - t) * y2_temp);
//PRINT("x:%f y:%f t:%f x1_temp:%f x2_temp:%f\n",x,y,t,x1_temp,x2_temp);
    // if (x >= 0 && x < im->w && y >= 0 && y < im->h) {
    if (x >= 0 && x < checkFrame.cols && y >= 0 && y < checkFrame.rows) {
      // augment number of checked points:
      n_points += 1;

      //Point(x,y) -> Point(col,row) x = col y = row 
      //getYUVColours(checkFrame, row, col, &Y, &U, &V); hence:
      getYUVColours(checkFrame,y,x, &y_c, &u_c, &v_c);
      if (y_c > ARF_Y_MIN && y_c < ARF_Y_MAX && u_c > ARF_U_MIN && u_c < ARF_Y_MAX && v_c > ARF_V_MIN && v_c < ARF_V_MAX) {
        // the point is of the right color:
        n_colored_points += 1;
      }
    }
  }
  PRINT("y_c:%d u_c:%d v_c:%d\n",y_c,u_c,v_c);
  return n_colored_points;
}

void sort_corners()
{
  
}

#if ARF_MOD_VIDEO
void mod_video(Mat& sourceFrame, Mat& frameGrey){
	char text[200];
#if ARF_MEASURE_FPS
    sprintf(text,"%5.2f %5.d %8.2fs", ARF_FPS,(runCount), curT / 1000000.0);
    putText(sourceFrame, text, Point(10,sourceFrame.rows-40), FONT_HERSHEY_PLAIN, 1, Scalar(0,255,255), 1);
#endif // CAM_STAB_MEASURE_FPS
	if(ARF_FLOOD_STYLE != ARF_FLOOD_CW){
#if ARF_DRAW_BOXES
		for(unsigned int r=0; r < cropAreas.size(); r++)
		{
			if(cropAreas[r].x != 0 && cropAreas[r].width != 0)
			{
				vector<Mat> channels;
				Mat thr_frame(cropAreas[r].height, cropAreas[r].width, CV_8UC2, cvScalar(0.0,0.0));
				Mat emptyCH(cropAreas[r].height, cropAreas[r].width, CV_8UC1, cvScalar(127.0));
				channels.push_back(emptyCH);
				channels.push_back(frameGrey(cropAreas[r]));
				merge(channels, thr_frame);
				thr_frame.copyTo(sourceFrame(cropAreas[r])); 			               // Copy threshold result to black frame
				emptyCH.release();
				thr_frame.release();
				rectangle(sourceFrame, cropAreas[r], Scalar(0,255), 2);
			}
		}
#endif //ARF_DRAW_BOXES
	}
#if ARF_BALL_CIRCLES && ARF_OBJECT == ARF_BALL
	for(unsigned int r=0; r < trackRes_size; r++)         // Convert angles & Write/Print output
	{
		circle(sourceFrame,cvPoint(trackRes[r].x_p - cropCol, trackRes[r].y_p), sqrt(trackRes[r].area_p / M_PI), cvScalar(100,255), 1);
	}
#endif //ARF_BALL_CIRCLES
#if ARF_DISTANCE_PLOT
	for(unsigned int r=0; r < trackRes_size; r++)         // Convert angles & Write/Print output
	{
	    line(sourceFrame, Point(0,sourceFrame.rows / 2.0), Point(trackRes[r].x_p - cropCol, trackRes[r].y_p), Scalar(0,255), 1);
	    sprintf(text,"%4.2f", trackRes[r].r_c);
	    putText(sourceFrame, text, Point((trackRes[r].x_p - cropCol + 0) / 2.0, (trackRes[r].y_p + sourceFrame.rows / 2.0) / 2.0), FONT_HERSHEY_PLAIN, 1, Scalar(0,255), 1);
	}
#endif
#if ARF_GATE_CORNERS && ARF_OBJECT == ARF_GATE
//printf("trackRes_size:%d---------------------\n",trackRes_size);
    Point p1, p2;
    
    uint8_t y_ = 0;
    uint8_t u_ = 0;
    uint8_t v_ = 0;
    
    //getYUVColours(sourceFrame, sourceFrame.rows / 2.0, sourceFrame.cols / 2.0, &y_, &u_, &v_);
    //PRINT("y_:%d\nu_:%d\nv_:%d\n",y_,u_,v_);
    
    int num_color_lines;
    
    //Point coordinates x positive up
    //y positive left
    
    int top_shift[gate_count];
    
      for(unsigned int r=0; r < gate_count; r++)         // Convert angles & Write/Print output
    {
      num_color_lines = 0;
      int max_av = 0;
        p1.x = hold_Gate[r].corners[3].x;
        p1.y = hold_Gate[r].corners[3].y;
        for(uint8_t i = 0; i < 4; i++){
            p2.x = hold_Gate[r].corners[i].x;
            p2.y = hold_Gate[r].corners[i].y;
	  
	    int av_x = (p1.x+p2.x)/2;
	    if(av_x > max_av){
	      top_shift[r] = i;
	      max_av = av_x;
	    }
	    
            p1 = p2;
        }
        
        PRINT("top_shift[%d]:%d\n",r,top_shift[r]);
    }
    
    for(unsigned int r=0; r < gate_count; r++)         // Convert angles & Write/Print output
    {
	num_color_lines = 0;
        p1.x = hold_Gate[r].corners[3].x;
        p1.y = hold_Gate[r].corners[3].y;
       // for(uint8_t i = 0; i < 4; i++){
            p2.x = hold_Gate[r].corners[0].x;
            p2.y = hold_Gate[r].corners[0].y;
            if(line_check(sourceFrame,p1,p2,0.6, 4)>10){
	      num_color_lines +=1;
	      line(sourceFrame, p1, p2, Scalar(250,128), 5);
	    }
            p1 = p2;
        //}
	  p2.x = hold_Gate[r].corners[1].x;
          p2.y = hold_Gate[r].corners[1].y;
            if(line_check(sourceFrame,p1,p2,0.6, 4)>10){
	      num_color_lines +=1;
	      line(sourceFrame, p1, p2, Scalar(100,100), 5);
	     }
          p1 = p2;
	  
	   p2.x = hold_Gate[r].corners[2].x;
          p2.y = hold_Gate[r].corners[2].y;
            if(line_check(sourceFrame,p1,p2,0.6, 4)>10){
	      num_color_lines +=1;
	      line(sourceFrame, p1, p2, Scalar(100,100), 5);
	     }
          p1 = p2;
	  
	   p2.x = hold_Gate[r].corners[3].x;
          p2.y = hold_Gate[r].corners[3].y;
            if(line_check(sourceFrame,p1,p2,0.6, 4)>10){
	      num_color_lines +=1;
	      line(sourceFrame, p1, p2, Scalar(100,100), 5);
	     }
          p1 = p2;
        
        PRINT("num_color_lines:%d\n",num_color_lines);
    }
    
    
    
    for(unsigned int r=0; r < gate_count; r++)         // Convert angles & Write/Print output
    {
      num_color_lines = 0;
        p1.x = hold_Gate[r].corners[3].x;
        p1.y = hold_Gate[r].corners[3].y;
        for(uint8_t i = 0; i < 4; i++){
            p2.x = hold_Gate[r].corners[i].x;
            p2.y = hold_Gate[r].corners[i].y;
	    if(i == 0) {
	      circle(sourceFrame, p2, 5, cvScalar(100,255), 5);
	    }
	    else{
	      circle(sourceFrame, p2, 5, cvScalar(110,140), 5);
	    }
//             if(line_check(sourceFrame,p1,p2,0.6, 4)>10){
// 	      num_color_lines +=1;
// 	      line(sourceFrame, p1, p2, Scalar(100,100), 5);
// 	    }
            p1 = p2;
        }
        
        //PRINT("num_color_lines:%d\n",num_color_lines);
    }
    
//     PRINT("p1.x:%d p1.y:%d\n",p1.x,p1.y);
//     int num_color_points = line_check(sourceFrame,p1,p2,0.6, 4);
//     PRINT("num_color_points:%d\n",num_color_points);
    
//     PRINT("half cols:%f\n",(sourceFrame.cols / 2.0));
//     PRINT("half rows:%f\n",(sourceFrame.rows / 2.0));
    //cv Point(x,y)
    
    
    line(sourceFrame,Point((sourceFrame.cols / 2.0),sourceFrame.rows / 2.0),Point((sourceFrame.cols / 2.0)+15,sourceFrame.rows / 2.0),Scalar(0,255), 1);
    line(sourceFrame,Point(sourceFrame.cols / 2.0,sourceFrame.rows / 2.0),Point(sourceFrame.cols / 2.0,(sourceFrame.rows / 2.0)+25),Scalar(0,255), 1);
    int half_crop = (sourceFrame.cols / 2.0);
    
    
    
    
    //PRINT("cropCol:%d half_crop:%d\n",cropCol,half_crop);
#endif //ARF_BALL_CIRCLES
#if ARF_DISTANCE_PLOT
    for(unsigned int r=0; r < trackRes_size; r++)         // Convert angles & Write/Print output
    {
        line(sourceFrame, Point(0,sourceFrame.rows / 2.0), Point(trackRes[r].x_p - cropCol, trackRes[r].y_p), Scalar(0,255), 1);
        sprintf(text,"%4.2f", trackRes[r].r_c);
        putText(sourceFrame, text, Point((trackRes[r].x_p - cropCol + 0) / 2.0, (trackRes[r].y_p + sourceFrame.rows / 2.0) / 2.0), FONT_HERSHEY_PLAIN, 1, Scalar(0,255), 1);
    }
#endif
#if ARF_PLOT_COORDS
	for(unsigned int r=0; r < neighbourMem_size; r++)         // Convert angles & Write/Print output
	{
	    if(neighbourMem[r].lastSeen < runCount){
	        uint8_t tColor = (uint8_t) round(255 * (ARF_MEMORY - (runCount - neighbourMem[r].lastSeen)) / ((float) ARF_MEMORY));
	        sprintf(text,"x%5.2f", neighbourMem[r].x_w);
	        putText(sourceFrame, text, Point(neighbourMem[r].x_p - cropCol + sqrt(neighbourMem[r].area_p / M_PI) + 10, neighbourMem[r].y_p - 15), FONT_HERSHEY_PLAIN, 1, Scalar(0,tColor), 1);
	        sprintf(text,"y%5.2f", neighbourMem[r].y_w);
	        putText(sourceFrame, text, Point(neighbourMem[r].x_p - cropCol + sqrt(neighbourMem[r].area_p / M_PI) + 10, neighbourMem[r].y_p), FONT_HERSHEY_PLAIN, 1, Scalar(0,tColor), 1);
	        sprintf(text,"z%5.2f", neighbourMem[r].z_w);
	        putText(sourceFrame, text, Point(neighbourMem[r].x_p - cropCol + sqrt(neighbourMem[r].area_p / M_PI) + 10, neighbourMem[r].y_p + 15), FONT_HERSHEY_PLAIN, 1, Scalar(0,tColor), 1);

	    }
	}
	for(unsigned int r=0; r < trackRes_size; r++)         // Convert angles & Write/Print output
	{
	    sprintf(text,"x%5.2f", trackRes[r].x_w);
	    putText(sourceFrame, text, Point(trackRes[r].x_p - cropCol + sqrt(trackRes[r].area_p / M_PI) + 10, trackRes[r].y_p - 15), FONT_HERSHEY_PLAIN, 1, Scalar(0,255,255), 1);
	    sprintf(text,"y%5.2f", trackRes[r].y_w);
	    putText(sourceFrame, text, Point(trackRes[r].x_p - cropCol + sqrt(trackRes[r].area_p / M_PI) + 10, trackRes[r].y_p), FONT_HERSHEY_PLAIN, 1, Scalar(0,255,255), 1);
	    sprintf(text,"z%5.2f", trackRes[r].z_w);
	    putText(sourceFrame, text, Point(trackRes[r].x_p - cropCol + sqrt(trackRes[r].area_p / M_PI) + 10, trackRes[r].y_p + 15), FONT_HERSHEY_PLAIN, 1, Scalar(0,255,255), 1);
	}
#endif
#if ARF_SHOW_STATS
	sprintf(text,"t:%4.1f%% o:%4.1f%%", pixCount/((float) ispHeight * ispWidth) * 100, pixSucCount/((float) pixCount) * 100);
	putText(sourceFrame, text, Point(10,sourceFrame.rows-80), FONT_HERSHEY_PLAIN, 1, Scalar(0,255,255), 1);
	sprintf(text,"d:%4.1f%% n:%4.1f%% s:%4.1f%%", pixDupCount/((float) pixCount) * 100, pixNofCount/((float) pixCount) * 100, pixSrcCount/((float) pixCount) * 100);
	putText(sourceFrame, text, Point(10,sourceFrame.rows-60), FONT_HERSHEY_PLAIN, 1, Scalar(0,255,255), 1);
#endif
#if ARF_SHOW_CAM_INFO
	sprintf(text,"R:%4.1f B:%4.1f G1:%4.1f G2:%4.1f Exp: %4.1f / %4.1f", mt9f002.gain_red, mt9f002.gain_blue, mt9f002.gain_green1, mt9f002.gain_green2, mt9f002.real_exposure, mt9f002.target_exposure);
	putText(sourceFrame, text, Point(10 , 40), FONT_HERSHEY_PLAIN, 1, Scalar(0,255,255), 1);
#endif
	return;
}
#endif // ARF_MOD_VIDEO

void active_random_filter_header( void ){
#if ARF_MEASURE_FPS
    clock_gettime(CLOCK_MONOTONIC, &time_now);
    curT            = sys_time_elapsed_us(&time_init, &time_now);
    uint32_t dt_us  = sys_time_elapsed_us(&time_prev, &time_now);
    ARF_FPS   = 0.975 * ARF_FPS + 0.025 * 1000000.f / dt_us;
    time_prev       = time_now;
    VERBOSE_PRINT("Measured FPS: %0.2f\n", ARF_FPS);
#endif
    trackRes_clear();
}

void active_random_filter_footer( struct FloatEulers*  eulerAngles ){
#if ARF_SHOW_MEM
    for(unsigned int r=0; r < neighbourMem_size; r++)        // Print to file & terminal
    {
        PRINT("%i - Object %d at (%0.2f m, %0.2f m, %0.2f m)\n", runCount, neighbourMem[r].id, neighbourMem[r].x_w, neighbourMem[r].y_w, neighbourMem[r].z_w);                                                        // Print to terminal
#if ARF_WRITE_LOG
        clock_gettime(CLOCK_MONOTONIC, &time_now);
        curT            = sys_time_elapsed_us(&time_init, &time_now);
        fprintf(arf_File,"%d\t%0.6f\t%d\t%d\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\n", runCount, curT / 1000000.f, neighbourMem[r].id, neighbourMem[r].lastSeen != runCount, pos->x, pos->y, pos->z, eulerAngles->psi, neighbourMem[r].x_w, neighbourMem[r].y_w, neighbourMem[r].z_w);
#endif
    }
    printf("\n");
#endif // ARF_SHOW_MEM
#if ARF_CALIBRATE_CAM
    if(runCount >= (ARF_TIMEOUT + 100)) calibrateEstimation();
#endif // ARF_CALIBRATE_CAM
    VERBOSE_PRINT("pixCount: %d  (%.2f%%), pixSucCount: %d  (%.2f%%), pixDupCount: %d  (%.2f%%), pixNofCount: %d  (%.2f%%), pixSrcCount: %d  (%.2f%%)\n", pixCount, pixCount/((float) ispHeight * ispWidth) * 100, pixSucCount, pixSucCount/((float) pixCount) * 100, pixDupCount, pixDupCount/((float) pixCount) * 100, pixNofCount, pixNofCount/((float) pixCount) * 100, pixSrcCount, pixSrcCount/((float) pixCount) * 100);
    runCount++;                                                // Increase counter
}

void addCrop(void){
    for(unsigned int r=0; r < cropAreas.size(); r++)
    {
        bool overlap = false;
        if(!overlap && (inRectangle(Point(objCrop.x, objCrop.y), cropAreas[r]) || inRectangle(Point(objCrop.x + objCrop.width, objCrop.y), cropAreas[r]) || inRectangle(Point(objCrop.x + objCrop.width, objCrop.y + objCrop.height), cropAreas[r]) || inRectangle(Point(objCrop.x, objCrop.y + objCrop.height), cropAreas[r])))
        {
            overlap = true; // One of the corner points is inside the cropAreas rectangle
        }
        if(!overlap && objCrop.x >= cropAreas[r].x && (objCrop.x + objCrop.width) <= (cropAreas[r].x + cropAreas[r].width) && objCrop.y <= cropAreas[r].y && (objCrop.y + objCrop.height) >= (cropAreas[r].y + cropAreas[r].height))
        {
            overlap = true; // less wide, yet fully overlapping in height
        }
        if(!overlap && objCrop.y >= cropAreas[r].y && (objCrop.y + objCrop.height) <= (cropAreas[r].y + cropAreas[r].height) && objCrop.x <= cropAreas[r].x && (objCrop.x + objCrop.width) >= (cropAreas[r].x + cropAreas[r].width))
        {
            overlap = true; // less high, yet fully overlapping in width
        }
        if(!overlap && (inRectangle(Point(cropAreas[r].x, cropAreas[r].y), objCrop) || inRectangle(Point(cropAreas[r].x + cropAreas[r].width, cropAreas[r].y), objCrop) || inRectangle(Point(cropAreas[r].x + cropAreas[r].width, cropAreas[r].y + cropAreas[r].height), objCrop) || inRectangle(Point(cropAreas[r].x, cropAreas[r].y + cropAreas[r].height), objCrop)))
        {
            overlap = true; // One of the corner points is inside the objCrop rectangle
        }
        if(overlap == true)
        {
            objCrop.width       = max(objCrop.x + objCrop.width, cropAreas[r].x + cropAreas[r].width) - min(objCrop.x, cropAreas[r].x);
            objCrop.height      = max(objCrop.y + objCrop.height, cropAreas[r].y + cropAreas[r].height) - min(objCrop.y, cropAreas[r].y);
            objCrop.x           = min(objCrop.x, cropAreas[r].x);
            objCrop.y           = min(objCrop.y, cropAreas[r].y);
            cropAreas[r].x      = 0;
            cropAreas[r].y      = 0;
            cropAreas[r].width  = 0;
            cropAreas[r].height = 0;
        }
    }
    if(objCrop.width * objCrop.height >= ARF_MIN_CROP_AREA * ispScalar * ispScalar)
    {
        cropAreas.push_back(objCrop);
    }
    return;
}

bool inRectangle(Point pt, Rect rectangle){
    if(pt.x >= rectangle.x && pt.x <= (rectangle.x + rectangle.width) && pt.y >= rectangle.y && pt.y <= (rectangle.y + rectangle.height))
    {
        return true;
    }else{
        return false;
    }
}

Rect enlargeRectangle(Mat& sourceFrame, Rect rectangle, double scale){
    int Hincrease       = round(scale / 2 * rectangle.width);
    int Vincrease       = round(scale / 2 * rectangle.height);
    rectangle.width     = min(sourceFrame.cols - 1, rectangle.x + rectangle.width + Hincrease) - max(0, rectangle.x - Hincrease);
    rectangle.height    = min(sourceFrame.rows - 1, rectangle.y + rectangle.height + Vincrease) - max(0, rectangle.y - Vincrease);
    rectangle.x         = max(0, rectangle.x - Hincrease);
    rectangle.y         = max(0, rectangle.y - Vincrease);
    return rectangle;
}

#if ARF_OBJECT == ARF_BALL 
bool trackRes_findMax( void ){
    trackRes_maxVal = 0.0;
    trackRes_maxId  = 0;
    for(uint8_t i=0; i < trackRes_size; i++){
        if(trackRes[i].r_c >= trackRes_maxVal){
            trackRes_maxVal = trackRes[i].r_c;
            trackRes_maxId  = i;
        }
    }
    return true;
}

bool trackRes_add( trackResults newRes, uint8_t overwriteId){
    if(overwriteId < trackRes_size){
        trackRes[overwriteId]           = newRes;
        trackRes_lastId                 = overwriteId;
        if(overwriteId == trackRes_maxId){
            if(trackRes_size == ARF_MAX_OBJECTS){
                trackRes_findMax();
            }
        }
    }
    else if(trackRes_size == ARF_MAX_OBJECTS){
        if(newRes.r_c < trackRes_maxVal){
            trackRes[trackRes_maxId]    = newRes;
            trackRes_lastId             = trackRes_maxId;
            trackRes_findMax();
        }
        else{
            return false;
        }
    }
    else{
        trackRes[trackRes_size]         = newRes;
        trackRes_lastId                 = trackRes_size;
        trackRes_size++;
        if(trackRes_size == ARF_MAX_OBJECTS){
            trackRes_findMax();
        }
    }
    return true;
}
#endif

bool trackRes_clear( void ){
    trackRes_size = 0;
    return true;
}

#if ARF_OBJECT == ARF_BALL
bool neighbourMem_findMax( void ){
#ifdef __linux__
  pthread_mutex_lock(&neighbourMem_mutex);
#endif
    neighbourMem_maxVal = 0.0;
    neighbourMem_maxId  = 0;
    for(uint8_t i=0; i<neighbourMem_size; i++){
        if(neighbourMem[i].r_c >= neighbourMem_maxVal){
            neighbourMem_maxVal = neighbourMem[i].r_c;
            neighbourMem_maxId  = i;
        }
    }
#ifdef __linux__
  pthread_mutex_unlock(&neighbourMem_mutex);
#endif
    return true;
}

bool neighbourMem_add( memoryBlock newRes, uint8_t overwriteId){
#ifdef __linux__
  pthread_mutex_lock(&neighbourMem_mutex);
#endif
    bool result = TRUE;
    if(overwriteId < neighbourMem_size){
        neighbourMem[overwriteId]           = newRes;
        neighbourMem_lastId                 = overwriteId;
    }
    else if(neighbourMem_size == ARF_MAX_OBJECTS){
        if(newRes.r_c < neighbourMem_maxVal){
            neighbourMem[neighbourMem_maxId]    = newRes;
            neighbourMem_lastId                 = neighbourMem_maxId;
        }
        else{
            result = FALSE;
        }
    }
    else{
        neighbourMem[neighbourMem_size]         = newRes;
        neighbourMem_lastId                     = neighbourMem_size;
        neighbourMem_size++;
    }
#ifdef __linux__
  pthread_mutex_unlock(&neighbourMem_mutex);
#endif
    if(result && neighbourMem_size == ARF_MAX_OBJECTS){
        neighbourMem_findMax();
    }
    return result;
}
#endif
