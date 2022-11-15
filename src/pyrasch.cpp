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
                std::cout<<exp(1.001)<<std::endl;
                std::cout<<"_________________________"<<std::endl;
                for (int n=0;n<N;n++)
                    {
                        for (int i=0;i<I;i++)
                            {
                                int m_i= MAX_ITEM_SCORES.coeff(i);
                                for (int m=0;m<m_i+1;m++)
                                    {
                                        std::cout<<data_probability.at(n).at(i).at(m)<<" ";
                                    }
                                
                                std::cout<<'\t';
                            }
                        
                        std::cout<<'\n';
                    }                
                std::cout<<"-------------------------"<<iter<<"---------------------"<<std::endl;
                for (int i=0;i<estimated_counts.size();i++)
                    {
                        for (int j=0;j<estimated_counts[i].size();j++)
                            {
                                std::cout<<estimated_counts[i][j]<<" ";
                            }
                        std::cout<<'\n';
                    }
                std::cout<<"--------"<<std::endl;
                for (int i=0;i<observed_counts.size();i++)
                    {
                        for (int j=0;j<observed_counts[i].size();j++)
                            {
                                std::cout<<observed_counts[i][j]<<" ";
                            }
                        std::cout<<'\n';
                    }
                std::cout<<"--------"<<std::endl;
                for (int i=0;i<RA_Thresholds.size();i++){
                    for (int j=0;j<RA_Thresholds[i].size();j++)
                        {
                            std::cout<<RA_Thresholds[i][j]<<" ";
                        }
                    
                    std::cout<<'\n';
                }

                std::cout<<"metric "<<residuals.rowwise().sum().colwise().sum()<<std::endl;
                std::cout<<ability<<std::endl;
                std::cout<<difficulty<<std::endl;
                //std::cout<<"threshold "<<RA_Thresholds<<std::endl;


                //std::cout<<"estimated_counts "<<estimated_counts<<std::endl;

                
                std::cout<<"EXPEc "<<expected_value<<std::endl;
                std::cout<<"variance "<<variance<<std::endl;
                //std::cout<<"probs "<<data_probability<<std::endl;
                std::cout<<"resuduals "<<residuals<<std::endl;

                estimate_counts();
                estimate_thresholds();
                estimate_full_probability();
                estimate_model_moments();
                calculate_residuals();

                
                estimate_difficulty();
                estimate_ability();

                
                
                
                
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
