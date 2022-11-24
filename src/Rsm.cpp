#include "Rsm.h"



void RSM::PROX(u_int16_t PROX_MAX)
    {
                        
    }


void RSM::JMLE(u_int16_t JMLE_MAX)
    {
        u_int16_t iter=0;

        
        for(iter; iter<JMLE_MAX;iter++)
            {
                
                estimate_counts();
                estimate_thresholds();
                estimate_full_probability();
                estimate_model_moments();
                calculate_residuals();

                
                estimate_difficulty();
                estimate_ability();

            }
    }

