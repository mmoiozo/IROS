/*
 * Copyright (C) Freek van Tienen
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
 * @file "modules/computer_vision/cv_ae_awb.c"
 * @author Freek van Tienen
 * Auto exposure and Auto white balancing for the Bebop 1 and 2
 */

#include "modules/computer_vision/cv_ae_awb.h"
#include "lib/isp/libisp.h"
#include <stdio.h>

#ifndef CV_AE_AWB_CAMERA
#define CV_AE_AWB_CAMERA front_camera
#endif

#ifndef CV_AE_AWB_FTOLERANCE
#define CV_AE_AWB_FTOLERANCE 0.05 // 0.3
#endif

#ifndef CV_AE_AWB_TARGET_AWB
#define CV_AE_AWB_TARGET_AWB 0.0
#endif

#ifndef CV_AE_AWB_MU
#define CV_AE_AWB_MU 0.1 // 0.0312
#endif

#ifndef CV_AE_AWB_VERBOSE
#define CV_AE_AWB_VERBOSE 0
#endif

#ifndef CV_AE_AWB_MIN_GAINS
#define CV_AE_AWB_MIN_GAINS 1
#endif

#ifndef CV_AE_AWB_MAX_GAINS
#define CV_AE_AWB_MAX_GAINS 64
#endif

#define PRINT(string,...) fprintf(stderr, "[cv_ae_awb->%s()] " string,__FUNCTION__ , ##__VA_ARGS__)

#if CV_AE_AWB_VERBOSE
#define VERBOSE_PRINT PRINT
#else
#define VERBOSE_PRINT(...)
#endif

#define MAX_HIST_Y 256-20
#define MIN_HIST_Y 16

#define awb_muAbs(x) ( ( x ) >= ( 0 ) ? ( x ) : ( -x ))
#define awb_muSign(x) ( ( x ) > ( 0 ) ? ( 1 ) : ( ( x ) < ( 0 ) ? ( -1 ) : ( 0 ) ) )
#define awb_muK(x) ( awb_muAbs(x) >= ( 2.5 ) ? ( 0.25 * awb_muSign( x ) ) : (awb_muAbs(x) >= ( 1.0 ) ? ( 0.125 * awb_muSign( x ) ) : ( awb_muAbs(x) >= ( 0.5 ) ? ( 0.075 * awb_muSign( x ) ) : ( awb_muAbs( x ) >= ( 0.2 ) ? 0.015 * awb_muSign( x ) : ( 0 )) )) )

bool gains_maxed            = false;
bool gains_minned           = false;
int8_t adjustDir            = 0;

float awb_fTolerance        = CV_AE_AWB_FTOLERANCE;
float awb_targetAWB         = CV_AE_AWB_TARGET_AWB;
float awb_mu                = CV_AE_AWB_MU;

uint8_t  bright_bin             = 211;
uint8_t  sat_bin                = 25;
uint8_t  middle_index           = 40;//85;

#include "boards/bebop/mt9f002.h"
struct image_t* cv_ae_awb_periodic(struct image_t* img);
struct image_t* cv_ae_awb_periodic(struct image_t* img) {
    struct isp_yuv_stats_t yuv_stats;
    if(isp_get_statistics_yuv(&yuv_stats) == 0 && yuv_stats.nb_valid_Y != 0) {
        // Calculate the CDF based on the histogram
        uint32_t cdf[MAX_HIST_Y];
        cdf[MIN_HIST_Y-1] = 0;
        VERBOSE_PRINT("nb_valid_Y: %d\n\n",yuv_stats.nb_valid_Y);
        for(int i = MIN_HIST_Y; i < MAX_HIST_Y; i++) {
            cdf[i] = cdf[i-1] + yuv_stats.ae_histogram_Y[i];
        }
        float adjustment                = 1.0f;
        int l = MIN_HIST_Y + 6;
        uint32_t total_pixels = 0;
        while (total_pixels < ( (cdf[MAX_HIST_Y - 1 - 6] - cdf[MIN_HIST_Y + 6]) / 2.0 ) && l < MAX_HIST_Y - 1 - 6) {
            total_pixels += yuv_stats.ae_histogram_Y[l];
            l++;
        }
        // that level is supposed to be in the middle of the bright_pixels bins
        adjustment = 1 + 0.5 * (( (float) (middle_index) ) / ( (float) l ) - 1);
        VERBOSE_PRINT("des_middle_index: %d,  actual_middle: %d,  adjustment: %f\n", middle_index, l, adjustment);
        // Calculate exposure
        Bound(adjustment, 1/4.0f, 2.0);
        float desiredExposure       = mt9f002.real_exposure * adjustment;
        mt9f002.target_exposure     = desiredExposure;
        mt9f002_set_exposure(&mt9f002);
        VERBOSE_PRINT("Desired exposure: %5.2f ms (real: %5.2f ms)\r\n", desiredExposure, mt9f002.real_exposure);
        // Calculate AWB (Robust Automatic White Balance Algorithm using Gray Color Points in Images - Huo et al.)
        if(yuv_stats.awb_nb_grey_pixels > 0){
            float avgU          = (((float) yuv_stats.awb_sum_U) / ((float) yuv_stats.awb_nb_grey_pixels) - 128);
            float avgV          = (((float) yuv_stats.awb_sum_V) / ((float) yuv_stats.awb_nb_grey_pixels) - 128);

            VERBOSE_PRINT("avgU = %d / %d - 127 = %f   avgU = %d / %d - 127 = %f\n", yuv_stats.awb_sum_U, yuv_stats.awb_nb_grey_pixels, avgU, yuv_stats.awb_sum_V, yuv_stats.awb_nb_grey_pixels, avgV);
            if(fabs(avgU) > (fabs(avgV) + awb_fTolerance) || ( fabs( fabs( avgU ) - fabs( avgV ) ) < awb_fTolerance && fabs( avgU ) > awb_fTolerance ) ){
                // filter output: avgU
                float error = awb_targetAWB - avgU;
                VERBOSE_PRINT("Adjust blue gain (error: %f)  blue_gain = %f + %f * %d\n", error, mt9f002.gain_blue, awb_mu, awb_muK(error));
                mt9f002.gain_blue   *= 1 + awb_mu * awb_muK(error);
                if(mt9f002.gain_blue > CV_AE_AWB_MAX_GAINS || mt9f002.gain_blue < CV_AE_AWB_MIN_GAINS){
                    mt9f002.gain_red    *= 1 - 0.5 * awb_mu * awb_muK(error);
                    mt9f002.gain_green1 *= 1 - 0.5 * awb_mu * awb_muK(error);
                    mt9f002.gain_green2 *= 1 - 0.5 * awb_mu * awb_muK(error);
                }
            }else if((fabs(avgU) + awb_fTolerance) < fabs(avgV)){
                // filter output: avgV
                float error = awb_targetAWB - avgV;
                VERBOSE_PRINT("Adjust red gain (error: %f)  red_gain = %f + %f * %d\n", error, mt9f002.gain_red, awb_mu, awb_muK(error));
                mt9f002.gain_red    *= 1 + awb_mu * awb_muK(error);
                if(mt9f002.gain_red > CV_AE_AWB_MAX_GAINS || mt9f002.gain_red < CV_AE_AWB_MIN_GAINS){
                    mt9f002.gain_blue   *= 1 - 0.5 * awb_mu * awb_muK(error);
                    mt9f002.gain_green1 *= 1 - 0.5 * awb_mu * awb_muK(error);
                    mt9f002.gain_green2 *= 1 - 0.5 * awb_mu * awb_muK(error);
                }
            }else{
                VERBOSE_PRINT("White balance achieved\n");
            }
            if((mt9f002.target_exposure > 37.5) && !gains_maxed)
            {
                // Let's try to decrease the exposure to be able to account for darker conditions better
                gains_minned           = false;
                mt9f002.gain_blue     *= 1 + 0.5 * awb_mu;
                mt9f002.gain_red      *= 1 + 0.5 * awb_mu;
                mt9f002.gain_green1   *= 1 + 0.5 * awb_mu;
                mt9f002.gain_green2   *= 1 + 0.5 * awb_mu;
            }
            else if((mt9f002.target_exposure < 1.5) && !gains_minned)
            {
                // Let's try to increase the exposure to be able to account for lighter conditions better
                gains_maxed            = false;
                mt9f002.gain_blue     *= 1 - 0.5 * awb_mu;
                mt9f002.gain_red      *= 1 - 0.5 * awb_mu;
                mt9f002.gain_green1   *= 1 - 0.5 * awb_mu;
                mt9f002.gain_green2   *= 1 - 0.5 * awb_mu;
            }
            else{
                // Let's try to center the exposure, that way we'll be able to account for fast changes in lightness
                if(((mt9f002.target_exposure > 15.0 || (mt9f002.target_exposure > 12.5 && adjustDir > 0))
                        && mt9f002.gain_blue   < (CV_AE_AWB_MAX_GAINS * 0.2)
                        && mt9f002.gain_red    < (CV_AE_AWB_MAX_GAINS * 0.2)
                        && mt9f002.gain_green1 < (CV_AE_AWB_MAX_GAINS * 0.2)
                        && mt9f002.gain_green2 < (CV_AE_AWB_MAX_GAINS * 0.2))
                        || ((mt9f002.gain_blue < (CV_AE_AWB_MIN_GAINS * 3.0)
                        || mt9f002.gain_red    < (CV_AE_AWB_MIN_GAINS * 3.0)
                        || mt9f002.gain_green1 < (CV_AE_AWB_MIN_GAINS * 3.0)
                        || mt9f002.gain_green2 < (CV_AE_AWB_MIN_GAINS * 3.0)) && mt9f002.target_exposure > 2.0)){
                    mt9f002.gain_blue     *= 1 + 0.01 * awb_mu;
                    mt9f002.gain_red      *= 1 + 0.01 * awb_mu;
                    mt9f002.gain_green1   *= 1 + 0.01 * awb_mu;
                    mt9f002.gain_green2   *= 1 + 0.01 * awb_mu;
                    gains_minned = false;
                    adjustDir = +1;
                }
                else if (((mt9f002.target_exposure < 10.0  || (mt9f002.target_exposure < 12.5 && adjustDir < 0))
                        && mt9f002.gain_blue   > (CV_AE_AWB_MIN_GAINS * 3.0)
                        && mt9f002.gain_red    > (CV_AE_AWB_MIN_GAINS * 3.0)
                        && mt9f002.gain_green1 > (CV_AE_AWB_MIN_GAINS * 3.0)
                        && mt9f002.gain_green2 > (CV_AE_AWB_MIN_GAINS * 3.0))
                        || ((mt9f002.gain_blue > (CV_AE_AWB_MAX_GAINS * 0.2)
                        || mt9f002.gain_red    > (CV_AE_AWB_MAX_GAINS * 0.2)
                        || mt9f002.gain_green1 > (CV_AE_AWB_MAX_GAINS * 0.2)
                        || mt9f002.gain_green2 > (CV_AE_AWB_MAX_GAINS * 0.2)) && mt9f002.target_exposure < 35.0 )){
                    mt9f002.gain_blue     *= 1 - 0.01 * awb_mu;
                    mt9f002.gain_red      *= 1 - 0.01 * awb_mu;
                    mt9f002.gain_green1   *= 1 - 0.01 * awb_mu;
                    mt9f002.gain_green2   *= 1 - 0.01 * awb_mu;
                    gains_maxed = false;
                    adjustDir = -1;
                }
                else{
                    adjustDir = 0;
                }
            }
            if(mt9f002.gain_blue >= (CV_AE_AWB_MAX_GAINS * 0.8) || mt9f002.gain_red >= (CV_AE_AWB_MAX_GAINS * 0.8)  || mt9f002.gain_green1 >= (CV_AE_AWB_MAX_GAINS * 0.8)  || mt9f002.gain_green2 >= (CV_AE_AWB_MAX_GAINS * 0.8)){
                gains_maxed = true;
            }
            if(mt9f002.gain_blue <= (CV_AE_AWB_MIN_GAINS * 1.2) || mt9f002.gain_red <= (CV_AE_AWB_MIN_GAINS * 1.2)  || mt9f002.gain_green1 <= (CV_AE_AWB_MIN_GAINS * 1.2)  || mt9f002.gain_green2 <= (CV_AE_AWB_MIN_GAINS * 1.2)){
                gains_minned = true;
            }
            Bound(mt9f002.gain_blue,    CV_AE_AWB_MIN_GAINS, CV_AE_AWB_MAX_GAINS);
            Bound(mt9f002.gain_red,     CV_AE_AWB_MIN_GAINS, CV_AE_AWB_MAX_GAINS);
            Bound(mt9f002.gain_green1,  CV_AE_AWB_MIN_GAINS, CV_AE_AWB_MAX_GAINS);
            Bound(mt9f002.gain_green2,  CV_AE_AWB_MIN_GAINS, CV_AE_AWB_MAX_GAINS);
            mt9f002_set_gains(&mt9f002);
        }
        else{
            VERBOSE_PRINT("Error: nb_grey_pixels = %d\n", yuv_stats.awb_nb_grey_pixels);
        }
    }
    else{
        VERBOSE_PRINT("Error: 0 valid Y pixels\n");
    }
    return img;
}

void cv_ae_awb_init(void) {
    cv_add_to_device(&CV_AE_AWB_CAMERA, cv_ae_awb_periodic);
}
