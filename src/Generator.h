#ifndef GENERATOR_H


#define GENERATOR_H
#include "Eigen/Dense"
#include <vector>
#include <random>
#include <algorithm>
#include <functional>
#include <numeric>

inline float normal_pdf(float x, float m, float s)
{
    static const float inv_sqrt_2pi = 0.3989422804014327;
    float a = (x - m) / s;

    return inv_sqrt_2pi / s * std::exp(-0.5f * a * a);
}

struct Generator
    {
        const int NUM_PERSON;
        const int NUM_ITEMS;
        
        Eigen::VectorXd difficulty;
        Eigen::VectorXd ability;

        std::vector<std::vector<double>> pcc_means;
        Eigen::VectorXd max_scores;

        Eigen::MatrixXd data;

        Generator(const int  t_num_person,const int t_num_items);

    };

Generator::Generator(const int  t_num_person,const int t_num_items):NUM_PERSON(t_num_person),NUM_ITEMS(t_num_items)
    {
        std::default_random_engine generator;
        std::normal_distribution<double> norm_distribution(0.0,3.0);
        std::uniform_int_distribution<int> uniform_distribution(0,5);

        std::random_device rd;
        std::mt19937 gen(rd());

        difficulty= Eigen::VectorXd::LinSpaced(NUM_ITEMS,-3.0,3.0);

        Eigen::VectorXd temp_ability(NUM_PERSON);
        Eigen::MatrixXd temp_data(NUM_PERSON,NUM_ITEMS);
        
        Eigen::VectorXd max_scores(NUM_ITEMS);
        
        std::vector<std::vector<double>> pcc_means;
        pcc_means.reserve(NUM_ITEMS);
        
        unsigned int i;
        unsigned int j;
        unsigned int n;
        unsigned int max_size;
        double max_score;
        double ability_n;
        double mean_ij;


        for (i=0;i<NUM_PERSON;i++)
            {
                temp_ability(i)=norm_distribution(generator);
            }
        
        ability =temp_ability;
        
        for (i=0;i<NUM_ITEMS;i++)
            {
                max_score=uniform_distribution(generator);
                max_scores(i)=max_score;

                

                std::vector<double> mean_points;
                mean_points.reserve(max_score+1);
                
                for (j=0;j<max_score+1;j++)
                    {
                        mean_points.emplace_back(norm_distribution(generator));
                    }
                
                std::sort(mean_points.begin(),mean_points.end());
                
                pcc_means.emplace_back(mean_points);
                
            }
        
        for (n=0;n<NUM_PERSON;n++)
            {
                ability_n= ability.coeff(n);
                for (i=0;i<NUM_ITEMS;i++)
                    {
                        max_size=pcc_means.at(i).size();
                        std::vector<double> raw_probs;
                        raw_probs.emplace_back(max_size);

                        

                        for (j=0;j<max_size;j++)
                            {
                                mean_ij=pcc_means.at(i).at(j);
                                raw_probs.emplace_back(normal_pdf(ability_n,mean_ij,1.0)*100);
                            }
                        
                        
                        std::discrete_distribution<> d{raw_probs.begin(),raw_probs.end()};
                        temp_data(n,i)=d(gen);
                    }
            }

        data=temp_data;
        

    }


#endif