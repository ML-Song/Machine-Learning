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
    
//    vector<double> test1{5, -3, 2};
//    vector<double> test2{6, -4, 4};
//    vector<double> test3{4, -4, 7};
//    //vector<double> test4{1, -1, 2};
//    vector<vector<double>> test5{test1, test2, test3};
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
//    cout << Matrix<double>(test1).NormVector() << endl;
//    tmp.GramSchmidt().Show();
//    tmp.QR().Transpose().Show();
//    tmp.QR().Inverse().Show();
//    for (auto &i: tmp.Eigen())
//        cout << i << " ";
//    tmp.Transpose().Shuffle().Show();
//    cout << (tmp == tmp_v) << endl;
//    (tmp / tmp).Show();
//    tmp.Show();
    
    
    
    
    
    clock_t end = clock();
    cout << endl << "Used Time: " << double (end - start) / CLOCKS_PER_SEC << "s" << endl;
    return OK;
}
