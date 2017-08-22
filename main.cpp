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

using namespace std;

int main(int argc, const char * argv[]) {
//    int input[2][2] = {1, 2, 3, 4};
//    Matrix<int> test(len(input), len(input[0]), &input[0][0]);
    clock_t start = clock();
    vector<double> test1{1, 2, 3};
    vector<double> test2{4, 5, 6};
    vector<double> test3{7, 8, 10};
    vector<vector<double>> test4{test1, test2, test3};
    Matrix<double> test(test4);
    Matrix<double> tmp(test);
    
    //tmp.SwapColumn(2, 3);
    //tmp.SwapRow(1, 2);
    //tmp.MultiplyColumn(1, -1, 2);
    //tmp.MultiplyRow(1, -1, 2);
    //tmp.RowEchelon();
    //tmp.SubMatrix(2, 3, 2, 3).Show();
    //tmp.Cofactor(2, 1).Show();
    //cout << tmp.Det() << endl;
    //tmp.Transpose();
    //tmp = test * test;
    //test.Transpose();
    //Dot(tmp, test).Show();
    //cout << tmp.IsInversible() << endl;
    //tmp.Inverse().Show();
    Dot(tmp.Inverse(), tmp).Show();
    //tmp.Show();
    //Matrix<double> tmp1(tmp.Cofactor(1, 2));
    //tmp1.Show();
    //cout << tmp1.Det() << endl;
    //tmp1.RowEchelon();
    //tmp1.Show();
    
    //cout << tmp.Det() << endl;
    
    clock_t end = clock();
    cout << double (end - start) / CLOCKS_PER_SEC << endl;
    return OK;
}
