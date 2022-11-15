#ifndef DATA_H

#define DATA_H
#include "Eigen/Dense"
#include <vector>
#include <iostream>
#include "RunningStats.h"

double LIMIT =2;
double CORRECTION =0.3;
double Minimim_Variance;

struct Rasch
    {
        Eigen::MatrixXd data;
        const int N;
        const int I;

        const int MAX_PERSON_RAW_SCORE;
        const Eigen::VectorXd MAX_ITEM_RAW_SCORE;
        const Eigen::VectorXi MAX_ITEM_SCORES;

        std::vector<std::vector<double>> RA_Thresholds;
        std::vector<std::vector<double>> observed_counts;
        std::vector<std::vector<double>> estimated_counts;

        std::vector<double> threshold_average;

        Eigen::VectorXd difficulty;
        Eigen::VectorXd ability;

        Eigen::MatrixXd expected_value;
        Eigen::MatrixXd variance;

        Eigen::MatrixXd residuals;

        std::vector<std::vector<std::vector<double>>> data_probability;

        Rasch(const Eigen::MatrixXd & t_data);

        const int find_max_person_score();

        const Eigen::VectorXd find_max_item_raw_scores();

        const Eigen::VectorXi find_max_item_scores();

        void estimate_thresholds();

        void estimate_counts();

        void PROX(int PROX_MAX);

        void JMLE(int JMLE_MAX);

        void estimate_model_moments();

        void estimate_ability();

        void estimate_difficulty();
        
        void calculate_residuals();

        void estimate_full_probability();      

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

        PROX(30);

        observed_counts.reserve(I);
        RA_Thresholds.reserve(I);
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

        //estimated_counts= observed_counts;
        estimate_full_probability();
        estimate_counts();

        threshold_average =std::vector<double>(I,0.0);

        //estimate_thresholds();

        estimate_model_moments();

        calculate_residuals();

    }

void Rasch::estimate_model_moments()
    {
        Eigen::MatrixXd out_expected_value(N,I);
        Eigen::MatrixXd out_variance(N,I);
        int i;
        int n;
        int m_i;
        int j;
        double mean;
        double var;

        for (n=0;n<N;n++)
            {   
                
                for(i=0;i<I;i++)
                    {
                        m_i= MAX_ITEM_SCORES.coeff(i);
                        Eigen::VectorXd scores;
                        scores.setLinSpaced(m_i,1,m_i);
                        
                        mean=0.0;
                        for (j=0;j<m_i-1;j++)
                            {
                                mean+=data_probability.at(n)[i].at(j+1)*scores.coeff(j);
                            }
                        var=0.0;
                        for (j=0;j<m_i;j++)
                            {
                                var+=(data.coeff(n,i)-mean)*(data.coeff(n,i)-mean);
                            }
                        var/=m_i;

                        //std::cout<<"mean :"<<mean<<std::endl;
                        //RunningStats stats = get_stats_obj(std::vector<double>(data_probability.at(n)[i].begin() + 1, data_probability.at(n)[i].end()),scores);
                        out_expected_value(n,i)=mean;
                        out_variance(n,i)=var;
                    }
            }
        
        expected_value= out_expected_value;
        variance = out_variance;

    }

void Rasch::estimate_full_probability()
    {
        std::vector<std::vector<std::vector<double>>> out;
        out.reserve(N);
        int n;
        int i;
        int j;
        int k;
        int l;

        double ability_n;
        double difficulty_i;
        
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
                        m_i=MAX_ITEM_SCORES.coeff(i);
                        difficulty_i=difficulty.coeff(i);

                        std::vector<double> person_n_item_i_probability;
                        person_n_item_i_probability.reserve(m_i+1);

                        dnom=0.0;
                        
                        for (j =0;j<m_i+1;j++)
                            {
                                double aux_sum=0.0;
                                for (l=0;l<j;l++)
                                    {
                                        aux_sum+=RA_Thresholds.at(i).at(l);
                                    }
                                dnom+=exp(j*(ability_n-difficulty_i)-aux_sum);
                                //std::cout<<"gere"<<std::endl;
                                
                            }
                        //dnom=1+exp((m_i)*(ability_n-difficulty_i)- aux_sum);

                        for(j=0;j<m_i+1;j++)
                            {
                                double aux_neom;
                                aux_neom=0.0;
                                for(l=0;l<j;l++)
                                    {
                                        aux_neom+=RA_Thresholds.at(i).at(l);
                                        //neom+=l*(ability_n-RA_Thresholds.at(i).at(l));
                                    }
                                neom=j*(ability_n-difficulty_i)-aux_neom;
                                prob= exp(neom)/dnom;
                                person_n_item_i_probability.emplace_back(prob);
                                //std::cout<<"prob "<<n<<" "<<i<<" "<<j<<" "<<"ability: "<<ability_n<<" "<<prob<<std::endl;
                            }
                        
                        person_n_probability.emplace_back(person_n_item_i_probability); 
                        person_n_item_i_probability.clear();                      

                    }
                
                out.emplace_back(person_n_probability);
                person_n_probability.clear();
            }
        
        data_probability = out;
        out.clear();
        
    }

void Rasch::estimate_counts()
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
        

        
        
        estimated_counts =out;
        
    }

void Rasch::estimate_thresholds()
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
        double temp_val;
        double temp_average;
        double temp_sum;

        for (i=0;i<I;i++)
            {
                m_i=MAX_ITEM_SCORES.coeff(i);
                std::vector<double> temp;
                temp.reserve(m_i+1);
                temp.emplace_back(0.0);
                temp_sum=0.0;

                for (j=1;j<m_i+1;j++)
                    {
                        observed_lower=observed_counts.at(i).at(j-1);
                        observed_higher=observed_counts.at(i).at(j);
                        estimated_lower=estimated_counts.at(i).at(j-1);
                        estimated_higher=estimated_counts.at(i).at(j);

                        if (observed_lower==0.0)
                            {
                                observed_lower=1;
                            }
                        
                        if (observed_higher==0.0)
                            {
                                observed_higher=1;
                            }

                        if (estimated_lower==0.0)
                            {
                                estimated_lower=1;
                            }
                        if (estimated_higher==0.0)
                            {
                                estimated_higher=1;
                            }

                        temp_val=RA_Thresholds.at(i).at(j)+log(observed_lower/observed_higher)-log(estimated_lower/estimated_higher);
                            
                        temp.emplace_back(temp_val);
                        temp_sum+=temp_val;

                    }
                
                temp_average=temp_sum/m_i;

                for (j=1;j<m_i+1;j++)
                    {
                        temp[j]-=temp_average;
                    }
                
                out.emplace_back(temp);

            }
        
        

        RA_Thresholds =out;
    }

void Rasch::calculate_residuals()
    {
        residuals= expected_value-data;
    }

void Rasch::estimate_difficulty()
    {
        Eigen::VectorXd temp;
        double temp_mean;
        temp =difficulty - (residuals.colwise().sum().array()/variance.colwise().sum().array()).matrix().transpose().unaryExpr([](double x) {return std::min(std::max(x,-LIMIT),LIMIT);});
        temp_mean = temp.mean();
        difficulty= temp.unaryExpr([&temp_mean](double x){return x-temp_mean;});

    }

void Rasch::estimate_ability()
    {
        Eigen::VectorXd temp;
        double temp_mean;
        temp = ability + (residuals.rowwise().sum().array()/variance.rowwise().sum().array()).matrix().unaryExpr([](double x) {return std::min(std::max(x,-LIMIT),LIMIT);});
        temp_mean= temp.mean();
        ability=temp.unaryExpr([&temp_mean](double x){return x-temp_mean;});

    }

#endif