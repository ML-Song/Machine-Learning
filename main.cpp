//
//  main.cpp
//  ML
//
//  Created by Bayes on 17/08/2017.
//  Copyright © 2017 Bayes. All rights reserved.
//

#include <iostream>
#include <time.h>

#include "matrix.h"
#include "statistical_learning.h"

using namespace std;

int main(int argc, const char * argv[]) {
    clock_t start = clock();
    
//    vector<double> test1{0, 0, 0};
//    vector<double> test2{0, 1, 0};
//    vector<double> test3{1, 0, 0};
//    vector<double> test4{1, 1, 1};
//    vector<vector<double>> test5{test1, test2, test3, test4};
//    Matrix<double> test(test5);
//    Matrix<double> tmp(test);
//    Matrix<double> tmp_v(test1);
    
//    Matrix<double>(3).Show();
//    Matrix<double>(3, 4, 4).Show();
//    Matrix<double>(test1).Show();
//    cout << tmp.Get(1, 1) << endl;
//    tmp.Set(0, 0, 3);
//    tmp.Show();
//    cout << tmp.GetColumn() << tmp.GetRow() << endl;
//    tmp.SubMatrix(2, 3, 2, 3).Show();
//    tmp.Cofactor(2, 2).Show();
//    tmp_v = tmp;
//    tmp_v.Show();
//    for (auto &i: tmp[1])
//        cout << i << " ";
//    tmp.SwitchRow(1, 2);
//    tmp.SwitchColumn(3, 2);
//    tmp.MultiplyColumn(1, 0.2);
//    tmp.MultiplyRow(2, 0.8);
//    tmp.AddColumn(1, 2, 2);
//    tmp.AddRow(2, 0.5, 3);
//    tmp.RowEchelon().Show();
//    tmp.Cat(tmp).Show();
//    cout << tmp.IsInversible() << " " << tmp.Det() << endl;
//    tmp.Transpose().Show();
//    (-tmp).Show();
//    tmp.Dot(tmp).Show();
//    tmp.Inverse().Dot(tmp).Show();
//    tmp.RowEchelon().Show();
//    cout << tmp.Tr() << endl;
//    cout << tmp.Rank() << endl;
//    cout << Matrix<double>(test1).NormVector_2() << endl;
//    tmp.GramSchmidt().Show();
//    tmp.QR().Transpose().Show();
//    tmp.QR().Inverse().Show();
//    for (auto &i: tmp.Eigen())
//        cout << i << " ";
//    tmp.Transpose().Shuffle().Show();
//    cout << (tmp == tmp_v) << endl;
//    (tmp / tmp).Show();
//    tmp.Show();
    vector<vector<double>> tmp;
    for (auto i: {0.0, 1.0}) {
        for (auto j: {0.0, 1.0})
        {
            for (auto k: {0.0, 1.0})
            {
                tmp.push_back(vector<double>{i, j, k, double (((int) i) && ((int) j) && ((int) k))});
            }
        }
    }
    Matrix<double> test(tmp);
    //test.Show();
    Perceptron perceptron(test, test);
    //cout << perceptron.Predict(Matrix<double>(test4));
    //cout << perceptron.Evaluate();
    //cout << Matrix<double>(test1).NormVector_1();
    perceptron.Train(2000, 0.01, 1e-5);
    perceptron.PrintWeights();
    cout << perceptron.Predict(Matrix<double>(vector<double>{0.5, 0.5, 0.5}));
    
    
    
    clock_t end = clock();
    cout << endl << "Used Time: " << double (end - start) / CLOCKS_PER_SEC << "s" << endl;
    return OK;
}
