#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <sstream>
#include <cmath>

using namespace std;

double Ro_inf(vector<double> a, vector<double> b) {
    double max = 0;
    for (int i = 0; i < a.size(); i++)
        if (max < fabs(a[i] - b[i])) //Функция fabs вычисляет абсолютное значение (модуль) и возвращает его |х|
            max = fabs(a[i] - b[i]);

    return max;
}

double Ro_1(vector<double> a, vector<double> b) {
    double sum = 0;
    for (int i = 0; i < a.size(); i++)
        sum += fabs(a[i] - b[i]);
    return sum;
}

double Ro_2(vector<double> a, vector<double> b) {
    double sum = 0;
    for (int i = 0; i < a.size(); i++)
        sum += (a[i] - b[i]) * (a[i] - b[i]);
    return sqrt(sum);
}

double Alpha(vector<vector<double>>& matrix) {
    double max_alpha = 0;
    for (int i = 0; i < matrix.size(); i++) {
        double cur_alpha = 0;
        for (int j = 0; j < matrix.size(); j++)
            cur_alpha += fabs(matrix[i][j]);

        if (cur_alpha > max_alpha)
            max_alpha = cur_alpha;
    }

    return max_alpha;
}

double Betta(vector<vector<double>>& matrix) {
    double max_betta = 0;
    for (int j = 0; j < matrix.size(); j++) {
        double cur_Betta = 0;
        for (int i = 0; i < matrix.size(); i++)
            cur_Betta += fabs(matrix[i][j]);

        if (cur_Betta > max_betta)
            max_betta = cur_Betta;
    }

    return max_betta;
}

double Gamma(vector<vector<double>>& matrix) {
    double gamma = 0;
    for (auto i : matrix)
        for (int j = 0; j < i.size() - 1; j++)
            gamma += (i[j] * i[j]);

    return gamma;
}

void outputMatrix(ostream& os, vector<vector<double>>& matrix) {
    stringstream tmp;
    int outWidth = 20;
    for (int i = 0; i < matrix.size(); i++) {
        for (int j = 0; j < matrix[i].size() - 1; j++) {
            tmp << (matrix[i][j] < 0 ? "-" : "+")
                << setprecision(5) << fabs(matrix[i][j])
                << " X_" << j + 1 << ' ';
            os << setw(outWidth) << tmp.str();
            tmp.str("");
        }

        tmp << (matrix[i][matrix[i].size() - 1] < 0 ? "-" : "+")
            << setprecision(5) << fabs(matrix[i][matrix[i].size() - 1]);
        os << setw(outWidth - 4) << tmp.str() << endl;
        tmp.str("");

    }
}

void Mat_A(vector<vector<double>>& matrix) {

    for (int i = 0; i < matrix.size(); i++) {
        double division = matrix[i][i];
        matrix[i][i] = 0;
        for (int j = 0; j < matrix[i].size(); j++)
            matrix[i][j] /= division;
        for (int j = 0; j < matrix[i].size() - 1; j++)
            matrix[i][j] *= (-1);
    }
}

void Metod_Iter(const char* fname) {
    ifstream fin(fname);
    int N; // kolvo peremennih
    fin >> N;

    double eps; // epsilon
    fin >> eps;

    vector<vector<double>> matrix(N, vector<double>(N + 1));
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N + 1; j++)
            fin >> matrix[i][j];

    fin.close();
    ofstream fout("output.txt");

    fout << "Ishodnay Matrica" << endl;

    outputMatrix(fout, matrix);

    Mat_A(matrix);

    fout << endl << "Matrica A" << endl;

    outputMatrix(fout, matrix);

    vector<double> cur_X(N);
    for (int i = 0; i < N; i++)
        cur_X[i] = matrix[i][N];

    double ABG[3];
    ABG[0] = Alpha(matrix); 
    ABG[1] = Betta(matrix);  
    ABG[2] = Gamma(matrix);   

    int outWidth = 20;

    fout << endl << setw(outWidth) << "alpha" << setw(outWidth) << "betta" << setw(outWidth) << "gamma" << endl;
    for (int i = 0; i < 3; i++)
        fout << setprecision(5) << setw(outWidth) << ABG[i];
    fout << endl;

    double mods[3];

    vector<double> zero_vec(N, 0);

    mods[0] = Ro_inf(zero_vec, cur_X);
    mods[1] = Ro_1(zero_vec, cur_X);
    mods[2] = Ro_2(zero_vec, cur_X);

    bool enable[3] = { 0 };
    for (int i = 0; i < 3; i++)
        enable[i] = (ABG[i] < 1);

    int iterations = 0;

    fout << endl << "Kolichestvo iteraci dla shogdenia metoda:" << endl << setw(outWidth) << "alpha" << setw(outWidth) << "betta" << setw(outWidth) << "gamma" << endl;
    for (int i = 0; i < 3; i++)
        if (enable[i]) {
            double n = log(eps * ((1 - ABG[i]) / mods[i])) / log(ABG[i]);
            n = ((int)n == n) ? n : (int)n + 1;

            if (iterations < n) iterations = n;

            fout << setw(outWidth) << (int)n;
        }
        else
            fout << setw(outWidth) << "-";

    fout << endl;
    fout << "Iteracii:" << endl;

    cur_X = vector<double>(N, 0);

    fout << "X_0" << ":\t";
    for (int i = 0; i < N; i++)
        fout << setprecision(5) << setw(outWidth) << 0;

    fout << endl;


    for (int i = 0; i < iterations; i++) {
        fout << "X_" << i + 1 << ":\t";
        vector<double> new_X(N, 0);
        for (int k = 0; k < N; k++) {
            for (int j = 0; j < N; j++)
                new_X[k] += matrix[k][j] * cur_X[j];
            new_X[k] += matrix[k][N];
        }

        for (int k = 0; k < N; k++)
            fout << setprecision(5) << setw(outWidth) << new_X[k];

        stringstream tmp;

        fout << setw(10) << "| Ro_inf:" << setprecision(5) << setw(outWidth) << Ro_inf(cur_X, new_X) << '\t';
        fout << setw(10) << "Ro_1:" << setprecision(5) << setw(outWidth) << Ro_1(cur_X, new_X) << '\t';
        fout << setw(10) << "Ro_2:" << setprecision(5) << setw(outWidth) << Ro_2(cur_X, new_X);

        fout << endl;

        cur_X = new_X;


    }

    fout.close();

}

int main() {
    Metod_Iter("input.txt");
    return 0;
}
