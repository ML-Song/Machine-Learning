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

class Perceptron
{
public:
    Perceptron(Matrix<double> train, Matrix<double>test);
    Matrix<double> GetWeights();
    Status PrintWeights();
    Status Train(long iteration=1000, double rate=0.01, double stop=1e-3);
    double Evaluate();
    bool Predict(Matrix<double> input);
    //感知机 二分类器
protected:
    Matrix<double> train_x_;
    Matrix<double> train_y_;
    Matrix<double> test_x_;
    Matrix<double> test_y_;
    Matrix<double> weigths_;
};

Perceptron::Perceptron(Matrix<double> train, Matrix<double>test)
{
    Matrix<double> ones(train.GetColumn(), 1, 1.0);
    //ones.Show();
    train_x_ = train.SubMatrix(1, train.GetColumn(), 1, train.GetRow() - 1).Cat(ones);
    train_y_ = train.SubMatrix(1, train.GetColumn(), train.GetRow(), train.GetRow());
    test_x_ = test.SubMatrix(1, test.GetColumn(), 1, test.GetRow() - 1).Cat(ones);
    test_y_ = test.SubMatrix(1, test.GetColumn(), test.GetRow(), test.GetRow());
    //train_x_.Show();
    //train_y_.Show();
    weigths_ = Matrix<double>(train.GetRow(), 1, 1.0);//列向量
    //(train_x_.Dot(weigths_) - train_y_).Show();
    //weigths_.Show();
}

Status Perceptron::PrintWeights()
{
    weigths_.Show();
    return OK;
}

Matrix<double> Perceptron::GetWeights()
{
    return weigths_;
}

bool Perceptron::Predict(Matrix<double> input)
{
    //input.Dot(weigths_).Show();
    input = input.Cat(Matrix<double>(input.GetColumn(), 1, 1.0));
    //input.Show();
    return input.Dot(weigths_).Get(0, 0) > 0;
}

double Perceptron::Evaluate()
{
    double result = 0;
    //train_x_.Show();
    //train_y_.Show();
    vector<double> delta = (train_x_.Dot(weigths_) * train_y_).Transpose().GetMat()[0];
    for (auto i: delta)
        result -= i;
    return result;
}

Status Perceptron::Train(long iteration, double rate, double stop)
{
    long count = 0;
    //train_x_.Show();
    //train_y_.Show();
    while (count < iteration) {
        Matrix<double> delta(train_x_.Transpose().Dot((train_y_ - train_x_.Dot(weigths_))));
        if (delta.NormVector_2() < stop)
            break;
        weigths_ = weigths_ + rate * delta;
        ++count;
    }
    return OK;
}
#endif /* statistical_learning_h */
