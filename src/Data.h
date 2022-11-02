#ifndef DATA_H

#define DATA_H
#include "Eigen/Dense"
#include <vector>
#include <unordered_map>

struct Rasch
    {
        const Eigen::MatrixXd data;
        const int N;
        const int I;

        const int MAX_PERSON_RAW_SCORE;
        const Eigen::VectorXd MAX_ITEM_RAW_SCORE;

        std::vector<std::vector<double>> RA_Thresholds;
        std::vector<std::vector<int>> observed_counts;

        Eigen::VectorXd difficulty;
        Eigen::VectorXd ability;

        std::vector<std::vector<std::vector<double>>> data_probability;

        Rasch(const Eigen::MatrixXd & t_data);

        const int find_max_person_score();

        const Eigen::VectorXd find_max_item_scores();

        std::vector<std::vector<double>> estimate_thresholds();

        void PROX(int PROX_MAX);

        double estimate_expected_value(const Eigen::VectorXd & probabilities, const int max_item_score);

        double estimate_probability(const double & ability_n, const Eigen::VectorXd & RA_threshold_i, const int &val);

        std::vector<std::vector<std::vector<double>>> estimate_full_probability();

        double estimate_expected_frequency();        

    };

const int Rasch::find_max_person_score()
    {
        return int(data.rowwise().sum().maxCoeff());
    }

const Eigen::VectorXd Rasch::find_max_item_scores()
    {
        return data.colwise().maxCoeff()*N;
    }

Rasch::Rasch(const Eigen::MatrixXd & t_data):data(t_data),N(data.rows()),I(data.cols()),MAX_PERSON_RAW_SCORE(find_max_person_score()),MAX_ITEM_RAW_SCORE(find_max_item_scores())
    {
        //MAX_PERSON_RAW_SCORE= find_max_person_score();

        //MAX_ITEM_RAW_SCORE = find_max_item_scores();

        difficulty = Eigen::VectorXd::Constant(I,0.0);
        ability = Eigen::VectorXd::Constant(N,0.0);

        RA_Thresholds.reserve(I);
        observed_counts.reserve(I);
        
        int j;
        int i=0;
        int m_i;

        
        for (auto col : data.colwise())
            {
                m_i=MAX_ITEM_RAW_SCORE.coeff(i);
                std::vector<int> temp;
                temp.reserve(m_i+1);
                std::vector<int> temp_map(m_i+1,0);
                for (j=0;j<N;j++)
                    {
                        int val= temp_map.at(int(col.coeff(j)));
                        val++;
                        temp_map[int(col.coeff(j))]=val;
                    }
                observed_counts.emplace_back(temp_map);
                RA_Thresholds.emplace_back(std::vector<double>(m_i+1,0.0));
                i++;
            }


    }

double Rasch::estimate_probability(const double & ability_n, const Eigen::VectorXd & RA_threshold_i, const int &val)
    {
        double neum =0.0;
        double dnom =1.0;
        double temp_dnom;
        int k;
        int j;
        int m =RA_threshold_i.size();

        for (k=0;k<val;k++)
            {
                neum+=ability_n-RA_threshold_i(k);
                temp_dnom=0.0;
                for (j=0;j<k;j++)
                    {
                        temp_dnom+=ability_n-RA_threshold_i(j);
                    }
                dnom+=exp(temp_dnom);

            }
        
        for (k=val;k<m;k++)
            {
                temp_dnom=0.0;
                for (j=0;j<k;j++)
                    {
                        temp_dnom+=ability_n-RA_threshold_i(j);
                    }
                dnom+=exp(temp_dnom);    
            }
        
        neum= exp(neum);

        return neum/dnom;
    }

double Rasch::estimate_expected_value(const Eigen::VectorXd & probabilities, const int max_item_score)
    {
        double sum=0.0;
        int i;
        for (i=0;i<max_item_score+1;i++)
            {
                sum+=probabilities.coeff(i)*i;
            }
        
        return sum;
    }

std::vector<std::vector<std::vector<double>>> Rasch::estimate_full_probability()
    {
        std::vector<std::vector<std::vector<double>>> out;
        out.reserve(N);
        int n;
        int i;
        int j;
        int k;
        int l;

        double ability_n;
        std::vector<double> RA_threshold_i;
        int m_i;

        double dnom;
        double neom;
        double temp_dnom;
        double prob;

        for (n=0;n<N;n++)
            {
                std::vector<std::vector<double>> person_n_probability;
                person_n_probability.reserve(I);
                ability_n= ability.coeff(n);

                for(i=0;i<I;i++)
                    {
                        RA_threshold_i = RA_Thresholds[i];
                        m_i=MAX_ITEM_RAW_SCORE.coeff(i);

                        std::vector<double> person_n_item_i_probability;
                        person_n_item_i_probability.reserve(m_i+1);

                        dnom=1.0;
                        for(j=0;j<m_i+1;j++)
                            {
                                temp_dnom=0.0;
                                for (k=0;k<j+1;k++)
                                    {
                                        temp_dnom+=ability_n-RA_threshold_i[k];
                                    }
                                dnom+=exp(temp_dnom);
                                
                            }
                        
                        for(j=0;j<m_i+1;j++)
                            {
                                neom=0.0;
                                for(l=0;l<j;l++)
                                    {
                                        neom+=ability_n-RA_threshold_i[l];
                                    }
                                prob= exp(neom)/dnom;
                                person_n_item_i_probability.emplace_back(prob);
                            }

                        person_n_probability.emplace_back(person_n_item_i_probability);                        

                    }
                out.emplace_back(person_n_probability);
            }
        
        return out;
    }

#endif