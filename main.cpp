//
//  main.cpp
//  ML
//
//  Created by Bayes on 17/08/2017.
//  Copyright © 2017 Bayes. All rights reserved.
//

#include <iostream>
#include <vector>

using namespace std;

typedef int Status;

#define OK 1;
#define ERR 0;
#define len(arr) (end(arr) - begin(arr))

template <typename T>
class Matrix {
public:
    Matrix();
    Matrix(long col, long row, T *input);
    Matrix(vector<vector<T>> input);
    //~Matrix();
    Status Show();
    Matrix<T> Dot(Matrix<T> m);
    Matrix<T> Inverse();
    Matrix<T> Transpose();
    Matrix<T> Tr();
    Matrix<T> Det();
    Matrix<T> operator-();
    Matrix<T> operator=(Matrix<T> m);
    friend Matrix<T> operator+(Matrix<T> m1, Matrix<T> m2);
    friend Matrix<T> operator-(Matrix<T> m1, Matrix<T> m2);
    friend Matrix<T> operator*(Matrix<T> m1, Matrix<T> m2);
protected:
    long column_ = 0;
    long row_ = 0;
    vector<vector<T>> content_;
};


int main(int argc, const char * argv[]) {
//    int input[2][2] = {1, 2, 3, 4};
//    Matrix<int> test(len(input), len(input[0]), &input[0][0]);
    vector<int> test1{1, 2};
    vector<int> test2{3, 4};
    vector<vector<int>> test3{test1, test2};
    Matrix<int> test(test3);
    test.Show();
    return OK;
}

template <typename T>
Matrix<T>::Matrix()
{
    column_ = 1;
    row_ = 1;
    content_.push_back(vector<T>(1));
}

template <typename T>
Matrix<T>::Matrix(long col, long row, T *input)
{
    column_ = col;
    row_ = row;
    for (int i = 0; i < column_; i++) {
        vector<T> tmp;
        for (int j = 0; j < row_; j++) {
            tmp.push_back(input[i * row_ + j]);
        }
        content_.push_back(tmp);
    }
}

template <typename T>
Matrix<T>::Matrix(vector<vector<T>> input)
{
    column_ = input.end() - input.begin();
    row_ = input[0].end() - input[0].begin();
    for (auto i: input)
    {
        content_.push_back(i);
    }
}

template <typename T>
Status Matrix<T>::Show()
{
    for (auto i: content_) {
        for (auto j: i) {
            cout << j << ' ';
        }
        cout << endl;
    }
    return OK;
}
