#include "RunningStats.h"
#include <pybind11/pybind11.h>
#include <iostream>
namespace py = pybind11;




double find_max_person_score(const Eigen::MatrixXd & data)
    {
        return data.rowwise().sum().maxCoeff();
    }

Eigen::VectorXd find_max_item_scores(const Eigen::MatrixXd & data)
    {
        return data.colwise().maxCoeff()*data.rows();
    }


struct parameters
    {
        Eigen::VectorXd ability;
        Eigen::VectorXd difficulty;
    };


parameters PROX_iter(const Eigen::MatrixXd & data,int PROX_MAX)
    {
        
        bool dichotomous=false;
        
        Eigen::VectorXd ability = Eigen::VectorXd::Constant(data.rows(), 0.0);
        
        Eigen::VectorXd difficulty = Eigen::VectorXd::Constant(data.cols(), 0.0);
        
        double MAX_PERSON_SCORE = find_max_person_score(data);
        
        Eigen::VectorXd MAX_ITEM_SCORE = find_max_item_scores(data);
        
        printf("%f",MAX_PERSON_SCORE);

        
        
        int iter = 0;

        
        for (iter;iter<PROX_MAX;iter++)
            {
                int n =0;
                int i=0;

                double difficulty_average;
                double difficulty_std;
                if (iter==0)
                {
                    difficulty_average=0;
                    difficulty_std = 0;
                }
                else
                {

                    RunningStats stats  = get_stats_obj(difficulty);
                    difficulty_average = stats.Mean();
                    difficulty_std = stats.StandardDeviation();
                }

                //printf("iter: %d, res : %f",iter,difficulty_average);

                double ability_const= sqrt(1+(difficulty_std*difficulty_std)/2.9);

                for (auto row : data.rowwise())
                {
                    double person_n_raw_score = row.sum();
                    //printf(" p score: %f",person_n_raw_score);

                    if (person_n_raw_score==0)
                    {
                        person_n_raw_score+=1;
                    } 
                    else if (person_n_raw_score==MAX_PERSON_SCORE)
                    {
                        person_n_raw_score-=1;
                    }
                    
                    ability(n)= difficulty_average+ability_const*log(person_n_raw_score/(MAX_PERSON_SCORE-person_n_raw_score));
                    //printf(" p ablity: %f ",ability(n));
                    n++;
                }

                RunningStats ability_stats  = get_stats_obj(ability);
                double temp_ability_average = ability_stats.Mean();
                ability = ability.array()-temp_ability_average;


                double ability_average;
                double ability_std;
                if (iter==0)
                {
                    ability_average=0;
                    ability_std = 0;
                }
                else
                {
                    RunningStats stats  = get_stats_obj(ability);
                    ability_average = stats.Mean();
                    ability_std = stats.StandardDeviation();
                }

                double difficulty_const= sqrt(1+(ability_std*ability_std)/2.9);
                
                for (auto column : data.colwise())
                {
                    double item_i_raw_score = column.sum();
                    //printf(" i score: %f",item_i_raw_score);

                    if (item_i_raw_score==0)
                    {
                        item_i_raw_score++;
                    } 
                    else if (item_i_raw_score==MAX_ITEM_SCORE(i))
                    {
                        item_i_raw_score--;
                    }
                    
                    difficulty(i)= ability_average - difficulty_const*log(item_i_raw_score/(MAX_ITEM_SCORE(i)-item_i_raw_score));
                    //printf(" i diff: %f",difficulty(i));
                    i++;
                }

                RunningStats diff_stats = get_stats_obj(difficulty);
                double temp_difficulty_average = diff_stats.Mean();
                difficulty = difficulty.array() -temp_difficulty_average;
                
            }
        
        parameters result = {ability,difficulty};

        return result;
        
    }



PYBIND11_MODULE(PROX,m)
    {
        m.doc() = "pybind11 example plugin";

        m.def("PROX_iter",&PROX_iter,py::arg("data"),py::arg("PROX_MAX"));

        py::class_<parameters>(m,"parameters")
        .def_readwrite("ability",&parameters::ability,py::return_value_policy::reference_internal)
        .def_readwrite("difficulty",&parameters::difficulty,py::return_value_policy::reference_internal);
        

    }
