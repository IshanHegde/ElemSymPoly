#ifndef PROB_H

#define PROB_H
#include "Eigen/Dense"


double PCM_pr(const double & difficulty, const Eigen::VectorXd & RA_threshold, const int &val)
    {
        double neum =0.0;
        double dnom =1.0;
        double temp_dnom;
        int k;
        int j;
        int m =RA_threshold.size();

        for (k=0;k<val;k++)
            {
                neum+=difficulty-RA_threshold(k);
                temp_dnom=0.0;
                for (j=0;j<k;j++)
                    {
                        temp_dnom+=difficulty-RA_threshold(j);
                    }
                dnom+=exp(temp_dnom);

            }
        
        for (k=val;k<m;k++)
            {
                temp_dnom=0.0;
                for (j=0;j<k;j++)
                    {
                        temp_dnom+=difficulty-RA_threshold(j);
                    }
                dnom+=exp(temp_dnom);    
            }
        
        neum= exp(neum);

        return neum/dnom;

    }

Eigen::MatrixXd PCM_E()


#endif