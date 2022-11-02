#ifndef DATA_H

#define DATA_H
#include "Eigen/Dense"

struct RaschData
    {
        const Eigen::MatriXd data;
        const int N;
        const int I;

        const int MAX_PERSON_RAW_SCORE;
        const Eigen::VectorXd MAX_ITEM_RAW_SCORE;

        Eigen::VectorXd<Eigen::VectorXd> RA_Threshold;

        Eigen::VectorXd difficulty;
        Eigen::VectorXd ability;

        RaschData(const Eigen::MatriXd & t_data);

        const int find_max_person_score(const Eigen::MatriXd & t_data);

        const Eigen::VectorXd find_max_item_scores(const Eigen::MatriXd & t_data);


    };

double RaschData::find_max_person_score(const Eigen::MatriXd & t_data)
    {
        return t_data.rowwise().sum().maxCoeff();
    }

Eigen::VectorXd RaschData::find_max_item_scores(const Eigen::MatriXd & t_data)
    {
        return t_data.colwise().maxCoeff()*t_data.rows();
    }

RaschData::RaschData(const Eigen::MatriXd & t_data):data(t_data),N(data.rows()),I(data.cols())
    {
        MAX_PERSON_RAW_SCORE= find_max_person_score(data);

        MAX_ITEM_RAW_SCORE = find_max_item_scores(data);

        difficulty = Eigen::VectorXd::Constant(I,0.0);
        ability = Eigen::VectorXd::Constant(N,0.0);

    }


#endif