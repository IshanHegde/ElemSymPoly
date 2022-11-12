#include "RunningStats.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <iostream>
#include "Data.h"
#include <pybind11/eigen.h>
namespace py = pybind11;



void Rasch::PROX(int PROX_MAX)
    {
                
        int iter = 0;

        Eigen::VectorXd x;
        

        
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
                    x.setLinSpaced(I,1,1);
                    RunningStats stats  = get_stats_obj(difficulty,x);
                    difficulty_average = stats.Mean();
                    difficulty_std = stats.StandardDeviation();
                }

  
                double ability_const= sqrt(1+(difficulty_std*difficulty_std)/2.9);

                for (auto row : data.rowwise())
                {
                    double person_n_raw_score = row.sum();
                    

                    if (person_n_raw_score==0)
                    {
                        person_n_raw_score+=0.3;
                    } 
                    else if (person_n_raw_score==MAX_PERSON_RAW_SCORE)
                    {
                        person_n_raw_score-=0.3;
                    }
                    
                    ability(n)= difficulty_average+ability_const*log(person_n_raw_score/(MAX_PERSON_RAW_SCORE-person_n_raw_score));
                    
                    n++;
                }
                
                x.setLinSpaced(N,1,1);
                RunningStats ability_stats  = get_stats_obj(ability,x);
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
                    x.setLinSpaced(N,1,1);
                    RunningStats stats  = get_stats_obj(ability,x);
                    ability_average = stats.Mean();
                    ability_std = stats.StandardDeviation();
                }

                double difficulty_const= sqrt(1+(ability_std*ability_std)/2.9);
                
                for (auto column : data.colwise())
                {
                    double item_i_raw_score = column.sum();

                    if (item_i_raw_score==0)
                    {
                        item_i_raw_score++;
                    } 
                    else if (item_i_raw_score==MAX_ITEM_RAW_SCORE(i))
                    {
                        item_i_raw_score--;
                    }
                    
                    difficulty(i)= ability_average - difficulty_const*log(item_i_raw_score/(MAX_ITEM_RAW_SCORE(i)-item_i_raw_score));

                    i++;
                }


                x.setLinSpaced(I,1,1);
                RunningStats diff_stats = get_stats_obj(difficulty,x);
                double temp_difficulty_average = diff_stats.Mean();
                difficulty = difficulty.array() -temp_difficulty_average;
                
            }
        
    }


void Rasch::JMLE(int JMLE_MAX)
    {
        int iter=0;

        for(iter; iter<JMLE_MAX;iter++)
            {
                estimate_difficulty();
                estimate_ability();
                estimate_thresholds();
                estimate_counts();
                estimate_model_moments();
                estimate_full_probability();
                calculate_residuals();
            }
    }

PYBIND11_MODULE(pyrasch,m)
    {
        m.doc() = "pybind11 example plugin";
        
        py::class_<Rasch>(m,"Rasch")
        .def(py::init<const Eigen::MatrixXd &>())
        .def("PROX",&Rasch::PROX)
        .def("JMLE",&Rasch::JMLE)
        .def("estimate_thresholds",&Rasch::estimate_thresholds)
        .def("estimate_counts",&Rasch::estimate_counts)
        .def_readwrite("expected_value",&Rasch::expected_value,py::return_value_policy::reference_internal)
        .def_readwrite("variance",&Rasch::variance,py::return_value_policy::reference_internal)
        .def_readwrite("data_probability",&Rasch::data_probability,py::return_value_policy::reference_internal)
        .def_readwrite("RA_Thresholds",&Rasch::RA_Thresholds,py::return_value_policy::reference_internal)
        .def_readwrite("observed_counts",&Rasch::observed_counts,py::return_value_policy::reference_internal)
        .def_readwrite("ability",&Rasch::ability,py::return_value_policy::reference_internal)
        .def_readwrite("difficulty",&Rasch::difficulty,py::return_value_policy::reference_internal);

    }
