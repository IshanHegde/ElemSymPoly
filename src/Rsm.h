#ifndef RSM_H

#define RSM_H
#include "Eigen/Dense"
#include <vector>
#include <map>
#include <iostream>

double LIMIT =2;
double CORRECTION =0.3;
double MINIMUM_VARIANCE=0;

// rating scale model
struct RSM
    {
        Eigen::MatrixXd data;
        const unsigned long long N;
        const unsigned long long I;

        const unsigned int MAX_PERSON_RAW_SCORE;
        const Eigen::VectorXd MAX_ITEM_RAW_SCORE;
        const Eigen::VectorXi MAX_ITEM_SCORES;

        std::map<double,double> RA_Thresholds;
        std::map<double,double> observed_counts;
        std::map<double,double> estimated_counts;

        unsigned int max_score;

        
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

const int RSM::find_max_person_score()
    {
        return int(data.rowwise().sum().maxCoeff());
    }

const Eigen::VectorXd RSM::find_max_item_raw_scores()
    {
        return data.colwise().maxCoeff()*N;
    }

const Eigen::VectorXi RSM::find_max_item_scores()
    {
        return data.cast<int>().colwise().maxCoeff();
    }

RSM::RSM(const Eigen::MatrixXd & t_data):data(t_data),N(data.rows()),I(data.cols()),MAX_PERSON_RAW_SCORE(find_max_person_score()),MAX_ITEM_RAW_SCORE(find_max_item_raw_scores()),MAX_ITEM_SCORES(find_max_item_scores())
    {

        difficulty = Eigen::VectorXd::Constant(I,0.0);
        ability = Eigen::VectorXd::Constant(N,0.0);

        size_t i;
        size_t n;

        max_score= MAX_ITEM_SCORES.maxCoeff();
        std::map<double,double> counts;
        std::map<double,double> estimated_counts;

        for ( i=0;i<max_score+1;i++)
            {
                counts[i]=0;
                RA_Thresholds[i]=0;
                estimated_counts[i]=0;
            }

        for (i=0;i<data.cols();i++)
            {
                for (n=0; n<data.rows();n++)
                    {
                        counts[data.coeff(n,i)]+=1;
                    }
            }

        observed_counts=counts;

        estimate_full_probability();
        estimate_counts();


        estimate_thresholds();

        estimate_model_moments();

        calculate_residuals();

    }

void RSM::estimate_model_moments()
    {
        Eigen::MatrixXd out_expected_value(N,I);
        Eigen::MatrixXd out_variance(N,I);
        size_t i;
        size_t n;
        size_t j;
        double mean;
        double var;

        for (n=0;n<N;n++)
            {
                for (i=0;i<I;i++)
                    {

                        mean=0.0;
                        for (j=0;j<max_score+1;j++)
                            {
                                mean+=data_probability.at(n).at(i).at(j)*j;
                            }
                        var=0.0;
                        for (j=0;j<max_score+1;j++)
                            {
                                var+=data_probability.at(n).at(i).at(j)*j*j;
                            }

                        var-=mean*mean;

                        out_expected_value(n,i)=mean;
                        out_variance(n,i)=var;
                    }
            }

        
        
        expected_value= out_expected_value;
        variance = out_variance;

    }

void RSM::estimate_full_probability()
    {
        std::vector<std::vector<std::vector<double>>> out;
        out.reserve(N);
        size_t n;
        size_t i;
        size_t j;
        size_t l;

        double ability_n;
        double difficulty_i;
        

        double dnom;
        double neom;
        double prob;

        for (n=0;n<N;n++)
            {
                std::vector<std::vector<double>> person_n_probability;
                person_n_probability.reserve(I);
                ability_n= ability.coeff(n);

                for(i=0;i<I;i++)
                    {
                        
                        difficulty_i=difficulty.coeff(i);

                        std::vector<double> person_n_item_i_probability;
                        person_n_item_i_probability.reserve(int(max_score)+1);

                        dnom=0.0;
                        
                        for (j =0;j<max_score+1;j++)
                            {
                                double aux_sum=0.0;
                                for (l=0;l<j+1;l++)
                                    {
                                        aux_sum+=RA_Thresholds.at(l);
                                    }
                                
                                dnom+=exp(j*(ability_n-difficulty_i)-aux_sum);
                                
                                
                            }
                        

                        for(j=0;j<max_score+1;j++)
                            {
                                double aux_neom;
                                aux_neom=0.0;
                                for(l=0;l<j+1;l++)
                                    {
                                        aux_neom+=RA_Thresholds.at(l);
                                        
                                    }
                                
                                neom=j*(ability_n-difficulty_i)-aux_neom;
                                prob= exp(neom)/dnom;
                                person_n_item_i_probability.emplace_back(prob);
                                
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

void RSM::estimate_counts()
    {
        std::map<double,double> out;
        
        size_t i;
        size_t j;
        size_t n;
        double sum;

        for (i=0;i<max_score+1;i++)
            {
                out[i]=0;
            }
        

                
        for (j=0;j<max_score+1;j++)
            {
                sum=0;
                for (i=0;i<I;i++)
                    {
                        for (n=0;n<N;n++)
                        {
                            sum+=data_probability.at(n).at(i).at(j);
                        }
                    }

                
                out[j]=sum;
            }
                
            
        

        estimated_counts =out;
        
    }

void RSM::estimate_thresholds()
    {
        std::vector<double> out;
        out.reserve(max_score+1);
        out.emplace_back(0.0);
        size_t i;
        size_t j;
        double observed_lower;
        double observed_higher;
        double estimated_lower;
        double estimated_higher;
        double temp_val;
        double temp_average;
        double temp_sum;

        temp_sum=0.0;
        for (j=1;j<max_score+1;j++)
            {
                observed_lower=observed_counts.at(j-1);
                observed_higher=observed_counts.at(j);
                estimated_lower=estimated_counts.at(j-1);
                estimated_higher=estimated_counts.at(j);
                
                if (observed_lower==0.0)
                    {
                        observed_lower=0.3;
                    }
                
                if (observed_higher==0.0)
                    {
                        observed_higher=1;
                    }

                if (estimated_lower==0.0)
                    {
                        estimated_lower=0.3;
                    }
                if (estimated_higher==0.0)
                    {
                        estimated_higher=1;
                    }
                
                temp_val=RA_Thresholds.at(j)+log(observed_lower/observed_higher)-log(estimated_lower/estimated_higher);
                    
                out.emplace_back(temp_val);
                temp_sum+=temp_val;
            }
        temp_average=temp_sum/(out.size()-1);
        for (i=1;i<max_score+1;i++)
            {
                
                RA_Thresholds[i]=out.at(i)-temp_average;
                
            }

        
        
    }

void RSM::calculate_residuals()
    {
        residuals= data-expected_value;
    }

void RSM::estimate_difficulty()
    {

        Eigen::VectorXd temp;
        double temp_mean;
        Eigen::VectorXd temp_variance = variance.colwise().sum().unaryExpr([MINIMUM_VARIANCE](double x){return std::max(x,MINIMUM_VARIANCE);});


        temp =difficulty - (residuals.colwise().sum().array()/temp_variance.transpose().array()).matrix().transpose().unaryExpr([](double x) {return std::min(std::max(x,-LIMIT),LIMIT);});
        temp_mean = temp.mean();
        difficulty= temp.unaryExpr([&temp_mean](double x){return x-temp_mean;});

    }

void RSM::estimate_ability()
    {
        Eigen::VectorXd temp;
        double temp_mean;
        Eigen::VectorXd temp_variance = variance.rowwise().sum().unaryExpr([MINIMUM_VARIANCE](double x){return std::max(x,MINIMUM_VARIANCE);});

        temp = ability + (residuals.rowwise().sum().array()/temp_variance.array()).matrix().unaryExpr([](double x) {return std::min(std::max(x,-LIMIT),LIMIT);});
        temp_mean= temp.mean();
        ability=temp.unaryExpr([&temp_mean](double x){return x-temp_mean;});

    }

#endif