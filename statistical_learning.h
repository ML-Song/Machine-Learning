//
//  statistical_learning.h
//  ML
//
//  Created by Bayes on 23/08/2017.
//  Copyright © 2017 Bayes. All rights reserved.
//

#ifndef statistical_learning_h
#define statistical_learning_h

#include "matrix.h"

using std::vector;
using std::cout;
using std::endl;

class Perceptron
{
public:
    Perceptron(Matrix<double> train, Matrix<double>test);
    vector<double> GetWeights();
    Status Train(long iteration=1000, double rate=0.1);
    vector<double> Evaluate();
    bool Predict(Matrix<double> input);
    //感知机 二分类器
protected:
    Matrix<double> train_;
    Matrix<double> test_;
    Matrix<double> weigths_;
};


#endif /* statistical_learning_h */
