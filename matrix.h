//
//  matrix.h
//  ML
//
//  Created by Bayes on 22/08/2017.
//  Copyright © 2017 Bayes. All rights reserved.
//

#ifndef matrix_h
#define matrix_h

#include <vector>
#include <math.h>

typedef int Status;

#define OK 1;
#define ERR 0;
#define EPSILON 1e-12

#define len(arr) (end(arr) - begin(arr))

using std::vector;
using std::cout;
using std::endl;

template <typename T> class Matrix;
template <typename T> Matrix<T> Dot(Matrix<T> const &m1, Matrix<T> const &m2);
//约束友元函数

template <typename T>
class Matrix {
public:
    Matrix();
    Matrix(long col, long row, T ele = T());
    Matrix(long col, long row, T *input);
    Matrix(vector<T> const &input);
    Matrix(vector<vector<T>> const &input);
    Matrix(Matrix<T> const &m);
    //~Matrix();
    Status Show();
    T Get(long i, long j);
    Status Set(long i, long j, T ele);
    Matrix<T> SubMatrix(long start_col, long end_col, long start_row, long end_row);
    Matrix<T> Cofactor(long x, long y);
    Matrix<T> operator=(Matrix<T> const &m);
    vector<T> operator[](long index);
    Status SwitchColumn(long i1, long i2);
    Status SwitchRow(long j1, long j2);
    Status MultiplyColumn(long i, double multiplier);
    Status MultiplyRow(long j, double multiplier);
    Status AddColumn(long i1, double multiplier, long i2);
    Status AddRow(long j1, double multiplier, long j2);
    Status RowEchelon(long *swap = nullptr);
    bool IsInversible();
    double Det();
    Matrix<T> Transpose();
    Matrix<T> operator-();
    friend Matrix<T> Dot<T>(Matrix<T> const &m1, Matrix<T> const &m2);
    Matrix<T> Inverse();
    T Tr();
    long Rank();
    double NormVector();
    Matrix<T> Schmidt();
    Matrix<T> GramSchmidt();
    Matrix<T> QR();// return Q
    vector<T> Eigen();
    
    friend Matrix<T> operator+(Matrix<T> const &m1, Matrix<T> const &m2)
    {
        if (m1.column_ != m2.column_ || m1.row_ != m2.row_)
        {
            throw std::invalid_argument("\ntwo matrix have different dimensions\n");
        }
        vector<vector<T>> tmp(m1.mat_);
        for (long i = 0; i < m1.column_; ++i) {
            for (long j = 0; j < m1.row_; ++j)
            {
                tmp[i][j] += m2.mat_[i][j];
            }
        }
        return Matrix<T>(tmp);
    }
    
    friend Matrix<T> operator*(T const &m1, Matrix<T> const &m2)
    {
        vector<vector<T>> tmp(m2.mat_);
        for (long i = 0; i < m2.column_; ++i) {
            for (long j = 0; j < m2.row_; ++j)
            {
                tmp[i][j] = m1 * m2.mat_[i][j];
            }
        }
        return Matrix<T>(tmp);
    }
    
    friend Matrix<T> operator-(Matrix<T> const &m1, Matrix<T> const &m2)
    {
        if (m1.column_ != m2.column_ || m1.row_ != m2.row_)
        {
            throw std::invalid_argument("\ntwo matrix have different dimensions\n");
        }
        vector<vector<T>> tmp(m1.mat_);
        for (long i = 0; i < m1.column_; ++i) {
            for (long j = 0; j < m1.row_; ++j)
            {
                tmp[i][j] -= m2.mat_[i][j];
            }
        }
        return Matrix<T>(tmp);
    }
    
    friend Matrix<T> operator*(Matrix<T> const &m1, Matrix<T> const &m2)
    {
        if (m1.column_ != m2.column_ || m1.row_ != m2.row_)
        {
            throw std::invalid_argument("\ntwo matrix have different dimensions\n");
        }
        vector<vector<T>> tmp(m1.mat_);
        for (long i = 0; i < m1.column_; ++i) {
            for (long j = 0; j < m1.row_; ++j)
            {
                tmp[i][j] *= m2.mat_[i][j];
            }
        }
        return Matrix<T>(tmp);
    }
    
    friend Matrix<T> operator/(Matrix<T> const &m1, Matrix<T> const &m2)
    {
        if (m1.column_ != m2.column_ || m1.row_ != m2.row_)
        {
            throw std::invalid_argument("\ntwo matrix have different dimensions\n");
        }
        vector<vector<T>> tmp(m1.mat_);
        for (long i = 0; i < m1.column_; ++i) {
            for (long j = 0; j < m1.row_; ++j)
            {
                tmp[i][j] /= m2.mat_[i][j];
            }
        }
        return Matrix<T>(tmp);
    }
protected:
    long column_ = 0;
    long row_ = 0;
    vector<vector<T>> mat_;
};


template <typename T>
Matrix<T>::Matrix()
{
    column_ = 1;
    row_ = 1;
    mat_.push_back(vector<T>(1));
}

template <typename T>
Matrix<T>::Matrix(long col, long row, T ele)
{
    column_ = col;
    row_ = row;
    mat_ = vector<vector<T>>(col, vector<T>(row, ele));
}

template <typename T>
Matrix<T>::Matrix(long col, long row, T *input)
{
    column_ = col;
    row_ = row;
    for (long i = 0; i < column_; ++i) {
        vector<T> tmp;
        for (long j = 0; j < row_; ++j) {
            tmp.push_back(input[i * row_ + j]);
        }
        mat_.push_back(tmp);
    }
}

template <typename T>
Matrix<T>::Matrix(vector<T> const &input)
{
    column_ = 1;
    row_ = input.end() - input.begin();
    mat_ = vector<vector<T>>{input};
}

template <typename T>
Matrix<T>::Matrix(vector<vector<T>> const &input)
{
    column_ = input.end() - input.begin();
    row_ = input[0].end() - input[0].begin();
    for (auto i: input)
    {
        mat_.push_back(i);
    }
}

template <typename T>
Matrix<T>::Matrix(Matrix<T> const &m)
{
    column_ = m.column_;
    row_ = m.row_;
    mat_ = m.mat_;
}

template <typename T>
Status Matrix<T>::Show()
{
    for (auto i: mat_) {
        for (auto j: i) {
            if (fabs(j) > EPSILON)
                printf("%-8.5lf ", j);
                //cout << setw(6) << j << ' ';
            else
                printf("0        ");
        }
        cout << endl;
    }
    return OK;
}

template <typename T>
T Matrix<T>::Get(long i, long j)
{
    return mat_[i][j];
}

template <typename T>
Status Matrix<T>::Set(long i, long j, T ele)
{
    mat_[i][j] = ele;
    return OK;
}

template <typename T>
Matrix<T> Matrix<T>::SubMatrix(long start_col, long end_col, long start_row, long end_row)
{
    start_row -= 1;
    start_col -= 1;
    //end_row -= 1;
    //end_col -= 1;
    vector<vector<T>> tmp(end_col - start_col, vector<T>(end_row - start_row));
    for (long i = start_col; i < end_col; ++i) {
        for (long j = start_row; j < end_row; ++j) {
            tmp[i - start_col][j - start_row] = mat_[i][j];
        }
    }
    return Matrix<T>(tmp);
}

template <typename T>
Matrix<T> Matrix<T>::operator=(Matrix<T> const &m)
{
    column_ = m.column_;
    row_ = m.row_;
    mat_ = m.mat_;
    return m;
}

template <typename T>
vector<T> Matrix<T>::operator[](long index)
{
    return mat_[index];
}

template <typename T>
Matrix<T> Matrix<T>::Cofactor(long x, long y)
{
    --x;
    --y;
    vector<vector<T>> tmp(column_ - 1, vector<T>(row_ - 1));
    for (long i = 0; i < column_; ++i) {
        long tmp_i;
        if (i < x)
        {
            tmp_i = i;
        }
        else if (i == x)
        {
            continue;
        }
        else
        {
            tmp_i = i - 1;
        }
        for (long j = 0; j < row_; ++j) {
            long tmp_j;
            if (j < y)
            {
                tmp_j = j;
            }
            else if (j == y)
            {
                continue;
            }
            else
            {
                tmp_j = j - 1;
            }
            tmp[tmp_i][tmp_j] = mat_[i][j];
        }
    }
    return Matrix<T>(tmp);
}

template <typename T>
Status Matrix<T>::SwitchColumn(long i1, long i2)
{
    --i1;
    --i2;
    vector<T> tmp(mat_[i1]);
    mat_[i1] = mat_[i2];
    mat_[i2] = tmp;
    return OK;
}

template <typename T>
Status Matrix<T>::SwitchRow(long j1, long j2)
{
    --j1;
    --j2;
    T tmp;
    for (long i = 0; i < column_; ++i)
    {
        tmp = mat_[i][j1];
        mat_[i][j1] = mat_[i][j2];
        mat_[i][j2] = tmp;
    }
    return OK;
}

template <typename T>
Status Matrix<T>::MultiplyColumn(long i, double multiplier)
{
    --i;
    for (long j = 0; j < row_; ++j)
    {
        mat_[i][j] *= multiplier;
    }
    return OK;
}

template <typename T>
Status Matrix<T>::MultiplyRow(long j, double multiplier)
{
    --j;
    for (long i = 0; i < row_; ++i)
    {
        mat_[i][j] *= multiplier;
    }
    return OK;
}

template <typename T>
Status Matrix<T>::AddColumn(long i1, double multiplier, long i2)
{
    --i1;
    --i2;
    for (long i = 0; i < row_; ++i)
    {
        mat_[i2][i] += multiplier * mat_[i1][i];
        if (fabs(mat_[i2][i]) < EPSILON)
            mat_[i2][i] = 0;
    }
    return OK;
}

template <typename T>
Status Matrix<T>::AddRow(long j1, double multiplier, long j2)
{
    --j1;
    --j2;
    for (long i = 0; i < column_; ++i)
    {
        mat_[i][j2] += multiplier * mat_[i][j1];
        if (fabs(mat_[i][j2]) < EPSILON)
            mat_[i][j2] = 0;
    }
    return OK;
}

template <typename T>
Status Matrix<T>::RowEchelon(long *swap)
{
//    vector<vector<T>> tmp, return_tmp;
//    for (auto i: mat_)
//    {
//        tmp.push_back(i);
//    }
    for (long i = 0; i < column_ - 1 && i < row_; ++i)
    {
        long max = i;
        for (long j = i; j < column_; ++j)
        {
            if (mat_[j][i] > mat_[max][i]) {
                max = j;
            }
        }
        if (i != max)
        {
            if (swap != nullptr)
                ++*swap;
            SwitchColumn(i + 1, max + 1);
        }
        if (fabs(mat_[i][i]) > EPSILON )
        {
            for (long j = i  + 1; j < column_; ++j) {
                AddColumn(i + 1, - mat_[j][i] / mat_[i][i], j + 1);
            }
        }
        else
        {
            ;//continue;
        }
        //Show();
    }
    //return_tmp = mat_;
    //mat_ = tmp;
    return OK;//Matrix<T>(return_tmp);
}

template <typename T>
bool Matrix<T>::IsInversible()
{
    if (column_ != row_)
        return false;
    vector<vector<T>> tmp(mat_);
    RowEchelon();
    for (long i = 0; i < column_; ++i)
    {
        if (fabs(mat_[i][i]) < EPSILON)
        {
            mat_ = tmp;
            return false;
        }
    }
    mat_ = tmp;
    return true;
}

template <typename T>
double Matrix<T>::Det()
{
    if (column_ != row_) {
        throw std::invalid_argument("\nNot Square Matrix\n");
    }
    vector<vector<T>> tmp(mat_);
    double det = 1;
    long swap = 0;
    RowEchelon(&swap);
    for (long i = column_ - 1; i >= 0; --i)
    {
        //cout << det << endl;
        if (fabs(mat_[i][i]) < EPSILON)
        {
            det = 0;
            break;
        }
        det *= mat_[i][i];
    }
    mat_ = tmp;
    return det * (swap % 2 ? -1 : 1);
}

template <typename T>
Matrix<T> Matrix<T>::Transpose()
{
    vector<vector<T>> tmp(row_, vector<T>(column_));
    for (long i = 0; i < column_; ++i) {
        for (long j = 0; j < row_; ++j) {
            tmp[j][i] = mat_[i][j];
        }
    }
    return Matrix<T>(tmp);
}

template <typename T>
Matrix<T> Matrix<T>::operator-()
{
    vector<vector<T>> tmp(mat_);
    for (auto &i: tmp) {
        for (auto &j: i)
        {
            j = -j;
        }
    }
    return Matrix<T>(tmp);
}

template <typename T>
Matrix<T> Dot(Matrix<T> const &m1, Matrix<T> const &m2)
{
    if (m1.row_ != m2.column_)
        throw std::invalid_argument("\nm1->row != m2->column\n");
    vector<vector<T>> tmp(m1.column_, vector<T>(m2.row_, 0));
    for (long x = 0; x < m1.column_; ++x)
        for (long y = 0; y < m2.row_; ++y)
            for (long i = 0; i < m1.row_; ++i)
            {
                tmp[x][y] += m1.mat_[x][i] * m2.mat_[i][y];
            }
    return Matrix<T>(tmp);
}

template <typename T>
long Matrix<T>::Rank()
{
    vector<vector<T>> tmp(mat_);
    RowEchelon();
    long rank = column_;
    for (long i = column_ - 1; i > 0; --i)
    {
        if (fabs(mat_[i][row_ - 1]) < EPSILON)
            --rank;
        else
            break;
    }
    mat_ = tmp;
    return rank;
}

template <typename T>
T Matrix<T>::Tr()
{
    if (column_ != row_)
        throw std::invalid_argument("\nNot Square Matrix\n");
    T tr = 0;
    for (long i = 0; i < column_; ++i)
    {
        tr += mat_[i][i];
    }
    return tr;
}

template <typename T>
Matrix<T> Matrix<T>::Inverse()
{
    if (!IsInversible())
        throw std::invalid_argument("\nNot Inversible\n");
    vector<vector<T>> tmp(mat_);
    for (long i = 0; i < column_; ++i)
        for (long j = 0; j < row_; ++j)
        {
            tmp[i][j] = Cofactor(i + 1, j + 1).Det() * ((i + j) % 2 ? -1 : 1);
        }
    Matrix<T> result(tmp);
    //result.Show();
    result = result.Transpose();
    return (1.0 / Det()) * result;
}

template <typename T>
double Matrix<T>::NormVector()
{
    if (column_ != 1 && row_ != 1)
    {
        throw std::invalid_argument("\nNot a Vector\n");
    }
    if (column_ == 1)
    {
        return sqrt(Dot(*this, this->Transpose()).Get(0, 0));
    }
    else
    {
        return sqrt(Dot(this->Transpose(), *this).Get(0, 0));
    }
}

template <typename T>
Matrix<T> Matrix<T>::Schmidt()
{
    vector<vector<T>> tmp(column_, vector<T>(row_));
    for (long i = 0; i < column_; ++i) {
        Matrix<T> eta(vector<T>{mat_[i]});
        Matrix<T> sum(eta);
        for (long j = 0; j < i; ++j) {
            Matrix<T> tmp_m(vector<T>{tmp[j]});
            sum = sum - ((Dot(eta, tmp_m.Transpose())) / Dot(tmp_m, tmp_m.Transpose())).Get(0, 0) * tmp_m;
            //tmp_m.Show();
            //eta.Show();
            //sum.Show();
        }
        tmp[i] = sum[0];
    }
    return Matrix<T>(tmp);
}

template <typename T>
Matrix<T> Matrix<T>::GramSchmidt()
{
    Matrix<T> tmp(*this);
    tmp = tmp.Transpose().Schmidt();
    for (long i = 1; i <= row_; ++i)
    {
        double norm = Matrix<T>(tmp[i - 1]).NormVector();
        tmp.MultiplyColumn(i, 1.0 / norm);
    }
    return tmp.Transpose();
}

template <typename T>
Matrix<T> Matrix<T>::QR()
{
    Matrix<T> tmp(column_, row_);
    for (long i = 0; i < column_; ++i) {
        for (long j = 0; j < row_; ++j) {
            tmp.Set(i, j, mat_[i][j]);
        }
    }
    //tmp.MultiplyRow(2, 2);
    return tmp.GramSchmidt();
}

template <typename T>
vector<T> Matrix<T>::Eigen()
{
    if (column_ != row_)
    {
        throw std::invalid_argument("\nNot Square Mat\n");
    }
    Matrix<T> Q, A = *this;
    vector<T> result;
    for (int i = 0; i < (column_ + row_) * 5; ++i) {
        Q = A.QR();
        A = Dot(Dot(Q.Transpose(), A), Q);
    }
    //A.Show();
    for (long i = 0; i < column_; ++i)
    {
        result.push_back(A.Get(i, i));
    }
    return result;
}
#endif /* matrix_h */
