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

#ifndef CV_AE_AWB_AV
#define CV_AE_AWB_AV 1
#endif

#ifndef CV_AE_AWB_TARGET_BRIGHT
#define CV_AE_AWB_TARGET_BRIGHT 0.05
#endif

#ifndef CV_AE_AWB_MAX_SATURATED
#define CV_AE_AWB_MAX_SATURATED 0.01
#endif

#ifndef CV_AE_AWB_FTOLERANCE
#define CV_AE_AWB_FTOLERANCE 0.05 // 0.3
#endif

#ifndef CV_AE_AWB_TARGET_AWB
#define CV_AE_AWB_TARGET_AWB 0.0
#endif

#ifndef CV_AE_AWB_MU
#define CV_AE_AWB_MU 0.0312 // 0.0312
#endif

#ifndef CV_AE_AWB_VERBOSE
#define CV_AE_AWB_VERBOSE 0
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
#define awb_muK(x) (awb_muAbs(x) >= ( 1.6 ) ? ( 4.0 * awb_muSign( x ) ) : ( awb_muAbs(x) >= ( 0.8 ) ? ( 2.0 * awb_muSign( x ) ) : ( awb_muAbs( x ) >= ( 0.05 ) ? awb_muSign( x ) : ( 0 )) ))

bool gains_maxed            = false;
float setting_target_bright = CV_AE_AWB_TARGET_BRIGHT;
float setting_max_saturated = CV_AE_AWB_MAX_SATURATED;
float awb_fTolerance        = CV_AE_AWB_FTOLERANCE;
float awb_targetAWB         = CV_AE_AWB_TARGET_AWB;
float awb_mu                = CV_AE_AWB_MU;

uint8_t  bright_bin             = 55;
uint8_t  sat_bin                = 25;

#include "boards/bebop/mt9f002.h"
struct image_t* cv_ae_awb_periodic(struct image_t* img);
struct image_t* cv_ae_awb_periodic(struct image_t* img) {
    struct isp_yuv_stats_t yuv_stats;
    if(isp_get_statistics_yuv(&yuv_stats) == 0 && yuv_stats.nb_valid_Y > 0.5 * img->h * img->w) {
        // Calculate the CDF based on the histogram
        uint32_t cdf[MAX_HIST_Y];
        cdf[MIN_HIST_Y-1] = 0;
        printf("nb_valid_Y: %d\n",yuv_stats.nb_valid_Y);
        for(int i = MIN_HIST_Y; i < MAX_HIST_Y; i++) {
            cdf[i] = cdf[i-1] + yuv_stats.ae_histogram_Y[i];
        }
        // Calculate bright and saturated pixels
        uint32_t bright_pixels          = cdf[MAX_HIST_Y - sat_bin] - cdf[MAX_HIST_Y - bright_bin];           // Top 30-15 bins
        uint32_t saturated_pixels       = cdf[MAX_HIST_Y - 1] - cdf[MAX_HIST_Y - sat_bin];            // top 15 bins
        uint32_t target_bright_pixels   = yuv_stats.nb_valid_Y * setting_target_bright;
        uint32_t max_saturated_pixels   = yuv_stats.nb_valid_Y * setting_max_saturated;
        float adjustment                = 1.0f;
        // Fix saturated pixels
        if(saturated_pixels > max_saturated_pixels) {
            adjustment = 1.0f - ( (float) (saturated_pixels - max_saturated_pixels) ) / yuv_stats.nb_valid_Y;
        }
        // Fix bright pixels
        else if (bright_pixels < target_bright_pixels) {
            // increase brightness to try and hit the desired number of well exposed pixels
            int l = MAX_HIST_Y - bright_bin;
            while (bright_pixels < target_bright_pixels && l > 0) {
                bright_pixels += yuv_stats.ae_histogram_Y[l];
                l--;
            }
            // that level is supposed to be at MAX_HIST_Y-31;
            adjustment = ( (float) (MAX_HIST_Y - bright_bin + 1) ) / (l+1);
        }
        else if (bright_pixels > target_bright_pixels) {
            // decrese brightness to try and hit the desired number of well exposed pixels
            int l = MAX_HIST_Y - bright_bin + 1;
            while (bright_pixels > target_bright_pixels && l < MAX_HIST_Y) {
                bright_pixels -= yuv_stats.ae_histogram_Y[l];
                l++;
            }
            // that level is supposed to be at MAX_HIST_Y-31;
            adjustment = ( (float) (MAX_HIST_Y - bright_bin + 1) ) / (l+1);
        }
        // Calculate exposure
        Bound(adjustment, 1/16.0f, 4.0);
        float desiredExposure       = mt9f002.real_exposure * adjustment * CV_AE_AWB_AV;
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
                mt9f002.gain_blue   += awb_mu * awb_muK(error);
            }else if((fabs(avgU) + awb_fTolerance) < fabs(avgV)){
                // filter output: avgV
                float error = awb_targetAWB - avgV;
                VERBOSE_PRINT("Adjust red gain (error: %f)  red_gain = %f + %f * %d\n", error, mt9f002.gain_red, awb_mu, awb_muK(error));
                mt9f002.gain_red    += awb_mu * awb_muK(error);
            }else{
                VERBOSE_PRINT("White balance achieved\n");
            }
            // filter output: 0
            //error = awb_targetAWB - 0.0;
            if((mt9f002.target_exposure > 25.0) && !gains_maxed)
            {
                mt9f002.gain_blue     += awb_mu;
                mt9f002.gain_red      += awb_mu;
                mt9f002.gain_green1   += awb_mu;
                mt9f002.gain_green2   += awb_mu;
            }
            else if((mt9f002.target_exposure < 0.4))
            {
                mt9f002.gain_blue     -= awb_mu;
                mt9f002.gain_red      -= awb_mu;
                mt9f002.gain_green1   -= awb_mu;
                mt9f002.gain_green2   -= awb_mu;
                /*if(mt9f002.gain_blue >= CV_AE_AWB_MAX_GAINS || mt9f002.gain_red >= CV_AE_AWB_MAX_GAINS  || mt9f002.gain_green1 >= CV_AE_AWB_MAX_GAINS  || mt9f002.gain_green2 >= CV_AE_AWB_MAX_GAINS){
                    gains_maxed = true;
                }
                VERBOSE_PRINT("Exposure saturated: new gains(B %f, R %f, G %f, G %f)\r\n", mt9f002.gain_blue, mt9f002.gain_red, mt9f002.gain_green1, mt9f002.gain_green2);
                */
            }
            Bound(mt9f002.gain_blue,    1, 75);
            Bound(mt9f002.gain_red,     1, 75);
            Bound(mt9f002.gain_green1,  1, 75);
            Bound(mt9f002.gain_green2,  1, 75);
            mt9f002_set_gains(&mt9f002);
        }
        else{
            PRINT("Error: nb_grey_pixels = %d\n", yuv_stats.awb_nb_grey_pixels);
        }
    }
    return img;
}

void cv_ae_awb_init(void) {
    cv_add_to_device(&CV_AE_AWB_CAMERA, cv_ae_awb_periodic);
}
