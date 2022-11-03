#ifndef DATA_H

#define DATA_H
#include "Eigen/Dense"
#include <vector>
#include <iostream>

struct Rasch
    {
        const Eigen::MatrixXd data;
        const int N;
        const int I;

        const int MAX_PERSON_RAW_SCORE;
        const Eigen::VectorXd MAX_ITEM_RAW_SCORE;
        const Eigen::VectorXi MAX_ITEM_SCORES;

        std::vector<std::vector<double>> RA_Thresholds;
        std::vector<std::vector<double>> observed_counts;
        std::vector<std::vector<double>> estimated_counts;

        Eigen::VectorXd difficulty;
        Eigen::VectorXd ability;

        Eigen::MatrixXd expected_value;
        Eigen::MatrixXd variance;

        std::vector<std::vector<std::vector<double>>> data_probability;

        Rasch(const Eigen::MatrixXd & t_data);

        const int find_max_person_score();

        const Eigen::VectorXd find_max_item_raw_scores();

        const Eigen::VectorXi find_max_item_scores();

        std::vector<std::vector<double>> estimate_thresholds();

        std::vector<std::vector<double>> estimate_counts();

        void PROX(int PROX_MAX);

        void JMLE(int JMLE_MAX);

        void estimate_model_moments();
        
        Eigen::MatrixXd calculate_residuals();

        std::vector<std::vector<std::vector<double>>> estimate_full_probability();      

    };

const int Rasch::find_max_person_score()
    {
        return int(data.rowwise().sum().maxCoeff());
    }

const Eigen::VectorXd Rasch::find_max_item_raw_scores()
    {
        return data.colwise().maxCoeff()*N;
    }

const Eigen::VectorXi Rasch::find_max_item_scores()
    {
        return data.cast<int>().colwise().maxCoeff();
    }

Rasch::Rasch(const Eigen::MatrixXd & t_data):data(t_data),N(data.rows()),I(data.cols()),MAX_PERSON_RAW_SCORE(find_max_person_score()),MAX_ITEM_RAW_SCORE(find_max_item_raw_scores()),MAX_ITEM_SCORES(find_max_item_scores())
    {

        difficulty = Eigen::VectorXd::Constant(I,0.0);
        ability = Eigen::VectorXd::Constant(N,0.0);

        RA_Thresholds.reserve(I);
        observed_counts.reserve(I);
        
        int j;
        int i=0;
        int m_i;

        
        for (auto col : data.colwise())
            {
                m_i=int(MAX_ITEM_SCORES.coeff(i));
                std::cout<<m_i<<std::endl;
                std::vector<double> temp;
                temp.reserve(m_i+1);
                std::vector<double> temp_map(m_i+1,0.0);
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

        estimated_counts= observed_counts;

        data_probability=estimate_full_probability();

    }

Eigen::MatrixXd Rasch::estimate_expected_values()
    {
        Eigen::MatrixXd out(N,I);
        int i;
        int j;
        int n;
        int m_i;
        double sum;
        

        for (n=0;n<N;n++)
            {
                for (i=0;i<I;i++)
                    {
                        m_i=MAX_ITEM_SCORES.coeff(i);
                        sum=0.0;
                        for (j=1;j<m_i+1;j++)
                            {
                                sum+=data_probability.at(n).at(i).at(j)*j;
                            }
                        out(n,i)=sum;
                    }
            }
            
        return out;
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
                        m_i=MAX_ITEM_SCORES.coeff(i);

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

std::vector<std::vector<double>> Rasch::estimate_counts()
    {
        std::vector<std::vector<double>> out;
        out.reserve(I);
        int i;
        int j;
        int n;
        double sum;

        for (i=0;i<I;i++)
            {
                std::vector<double> item_prob;
                item_prob.reserve(MAX_ITEM_SCORES.coeff(i)+1);
                for (j=0;j<MAX_ITEM_SCORES.coeff(i)+1;j++)
                    {
                        sum=0;
                        for (n=0;n<N;n++)
                            {
                                sum+=data_probability.at(n).at(i).at(j);
                            }
                        item_prob.emplace_back(sum);
                    }
                out.emplace_back(item_prob);
            }
        
        return out;
        
    }

std::vector<std::vector<double>> Rasch::estimate_thresholds()
    {
        std::vector<std::vector<double>> out;
        out.reserve(I);
        int i;
        int j;
        int k;
        int m_i;
        double observed_lower;
        double observed_higher;
        double estimated_lower;
        double estimated_higher;

        for (i=0;i<I;i++)
            {
                m_i=MAX_ITEM_SCORES.coeff(i);
                std::vector<double> temp;
                temp.reserve(m_i+1);

                temp.emplace_back(0.0);
                for (j=1;j<m_i+1;j++)
                    {
                        observed_lower=observed_counts.at(i).at(j-1);
                        observed_higher=observed_counts.at(i).at(j);
                        estimated_lower=estimated_counts.at(i).at(j-1);
                        estimated_higher=estimated_counts.at(i).at(j);
                            
                        temp.emplace_back(estimated_counts.at(i).at(j)+log(observed_lower/observed_higher)-log(estimated_lower/estimated_higher));


                    }
                
                out.emplace_back(temp);

            }
        return out;
    }

Eigen::MatrixXd Rasch::calculate_residuals()
    {
        Eigen::MatrixXd out(N,I);

        return estimate_expected_values()-data;
    }
#endif