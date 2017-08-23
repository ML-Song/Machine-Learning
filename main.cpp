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
    
    vector<double> test1{5, -3, 2};
    vector<double> test2{6, -4, 4};
    vector<double> test3{4, -4, 7};
    //vector<double> test4{1, -1, 2};
    vector<vector<double>> test5{test1, test2, test3};
    Matrix<double> test(test5);
    Matrix<double> tmp(test);
    
    
//    tmp.SwapColumn(2, 3);
//    tmp.SwapRow(1, 2);
//    tmp.MultiplyColumn(1, -1, 2);
//    tmp.MultiplyRow(1, -1, 2);
//    tmp.RowEchelon();
//    tmp.SubMatrix(2, 3, 2, 3).Show();
//    tmp.Cofactor(2, 1).Show();
//    cout << tmp.Det() << endl;
//    tmp.Transpose().Show();
//    tmp = test * test;
//    Dot(tmp, test).Show();
//    cout << tmp.IsInversible() << endl;
//    tmp.Inverse().Show();
//    Dot(tmp.Inverse(), tmp).Show();
//    tmp.Show();
//    Matrix<double> tmp1(tmp.Cofactor(1, 2));
//    tmp1.Show();
//    cout << tmp1.Det() << endl;
//    tmp1.RowEchelon();
//    tmp1.Show();
//    (test - test).Show();
//    tmp.Show();
//    cout << tmp.Det() << endl;
//    tmp.Transpose().Schmidt().Transpose().Show();
//    tmp.MultiplyRow(2, 2);
//    tmp.GramSchmidt().Show();
//    tmp.QR().Transpose().Show();
//    tmp.QR().Inverse().Show();
//    Dot(tmp.QR().Inverse(), tmp).Show();
//    for (auto i: tmp.Eigen())
//        cout << i << endl;
//    Cat(tmp.Transpose(), tmp.Transpose()).Transpose().Show();
    
    
    clock_t end = clock();
    cout << endl << "Used Time: " << double (end - start) / CLOCKS_PER_SEC << "s" << endl;
    return OK;
}
