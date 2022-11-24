#include "Rsm.h"



void RSM::PROX(int PROX_MAX)
    {
                        
    }


void RSM::JMLE(int JMLE_MAX)
    {
        int iter=0;

        
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

