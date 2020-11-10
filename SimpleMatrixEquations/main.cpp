#include <iostream>
#include <cassert>
int swap_count = 0;

class Row{
private:
    int size;
    double* numbers;
public:
    explicit Row(const double* initial_numbers, int elements_number){
        numbers = new double[elements_number];
        for(int i=0; i < elements_number; i++)
            numbers[i] = initial_numbers[i];
        size = elements_number;
    }
    explicit Row(double value, int elements_number){
        numbers = new double[elements_number];
        for(int i=0; i < elements_number; i++)
            numbers[i] = value;
        size = elements_number;
    }
    explicit Row(int elements_number){
        numbers = new double[elements_number];
        for(int i=0; i < elements_number; i++)
            numbers[i] = 0;
        size = 0;
    }

    void print(){
        for(int i=0; i < size; i++) {
            printf("%.3lf", numbers[i]);
            for(int j=0; j < 3; j++)
                printf(" ");
        }
    }

    double get(int j){
        return numbers[j];
    }

    friend Row operator+(const Row &row, double value);
    friend Row operator+(const Row &row1, const Row &row2);
    friend Row operator-(const Row &row1, double value);
    friend Row operator-(const Row &row1, const Row &row2);
    friend Row operator*(const Row &row1, double value);
    friend Row operator/(const Row &row1, double value);
    double& operator[](int index);

    void operator+=(double value){
        for(int i=0; i < size; i++)
            numbers[i] += value;
    }
    void operator+=(const Row &row){
        for(int i=0; i < size; i++)
            numbers[i] += row.numbers[i];
    }
    void operator-=(double value){
        for(int i=0; i < size; i++)
            numbers[i] -= value;
    }
    void operator-=(const Row &row){
        for(int i=0; i < size; i++)
            numbers[i] -= row.numbers[i];
    }
    void operator*=(double value){
        for(int i=0; i < size; i++)
            numbers[i] *= value;
    }
    void operator/=(double value){
        for(int i=0; i < size; i++)
            numbers[i] /= value;
    }


};
double& Row::operator[](int index){
    return numbers[index];
}

Row operator+(const Row &row, double value){
    auto *new_numbers = new double[row.size];
    for(int i=0; i < row.size; i++)
        new_numbers[i] = row.numbers[i] + value;
    Row new_row = Row(new_numbers, row.size);
    delete[] new_numbers;
    return new_row;
}
Row operator+(const Row &row1, const Row &row2){
    auto *new_numbers = new double [row1.size];
    for(int i=0; i < row1.size; i++)
        new_numbers[i] = row1.numbers[i] + row2.numbers[i];
    Row new_row = Row(new_numbers, row1.size);
    delete[] new_numbers;
    return new_row;
}
Row operator-(const Row &row, double value){
    auto *new_numbers = new double[row.size];
    for(int i=0; i < row.size; i++)
        new_numbers[i] = row.numbers[i] - value;
    Row new_row = Row(new_numbers, row.size);
    delete[] new_numbers;
    return new_row;
}
Row operator-(const Row &row1, const Row &row2){
    auto *new_numbers = new double[row1.size];
    for(int i=0; i < row1.size; i++)
        new_numbers[i] = row1.numbers[i] - row2.numbers[i];
    Row new_row = Row(new_numbers, row1.size);
    delete[] new_numbers;
    return new_row;
}
Row operator*(const Row &row, double value){
    auto *new_numbers = new double[row.size];
    for(int i=0; i < row.size; i++)
        new_numbers[i] = row.numbers[i] * value;
    Row new_row = Row(new_numbers, row.size);
    delete[] new_numbers;
    return new_row;
}
Row operator/(const Row &row, double value){
    auto *new_numbers = new double[row.size];
    for(int i=0; i < row.size; i++)
        new_numbers[i] = row.numbers[i] / value;
    Row new_row = Row(new_numbers, row.size);
    delete[] new_numbers;
    return new_row;
}


class Matrix{
private:
    Row** rows;
    int number_of_lines;
    int line_size;
    static int find_next_nonzero(Matrix matrix, int current_row, int current_element, int size){
    int n = current_row + 1;
    while (n < size && matrix[n][current_element] == 0)
        n++;
    if (n == size && matrix[n][current_element] == 0)
        return 0; //stop gauss, answer is 0
    else {
        matrix.swap_rows(current_row, n);
        swap_count++;
        return 1; //success, continue gauss
    }
}
    int gauss_for_det(Matrix matrix, int current_row, int current_element, int m_size){
    if (current_row == m_size && matrix[current_row] [current_element] != 0)
        return 1;
    if (current_row == m_size)
        return 0;
    if (matrix[current_row][current_element] == 0 && find_next_nonzero(matrix, current_row, current_element, m_size) == 0)
        return 0;

    for(int i = current_row + 1; i <= m_size; i++){
        if (matrix[i][current_element] != 0) {
            double ratio = matrix[i][current_element] / matrix[current_row][current_element];
            matrix[i] -= matrix[current_row] * ratio;
        }
    }
    return gauss_for_det(matrix, current_row + 1, current_element + 1, m_size);
}
public:
    explicit Matrix(Row** mat, int m, int n){
        rows = new Row*[m];
        for(int i=0; i < m; i++)
            rows[i] = mat[i];
        line_size = n;
        number_of_lines = m;
    }
    explicit Matrix(Row** mat, int n){
        rows = new Row*[n];
        for(int i=0; i < n; i++)
            rows[i] = mat[i];
        line_size = n;
        number_of_lines = n;
    }
    explicit Matrix(int m, int n, double value){
        rows = new Row*[m];
        for(int i=0; i < n; i++)
            rows[i] = new Row(value, n);
        line_size = n;
        number_of_lines = m;
    }
    explicit Matrix (int m, int n){
        rows = new Row*[m];
        for(int i=0; i < n; i++)
            rows[i] = new Row(n);
        line_size = n;
        number_of_lines = m;
    }
    explicit Matrix (int n){
        rows = new Row*[n];
        for(int i=0; i < n; i++)
            rows[i] = new Row(n);
        line_size = n;
        number_of_lines = n;
    }

    void print(){
        for(int i=0; i < number_of_lines; i++) {
            rows[i]->print();
            std::cout << "\n";
        }
    }

    void swap_rows(int i, int j){
        std::swap(rows[i], rows[j]);
    }

    double determinator(){
        assert(number_of_lines == line_size);
        Row** r = new Row*[line_size];
        for(int i = 0; i < line_size; i++){
            auto *row = new double [line_size];
            for(int j=0; j < line_size; j++)
                row[j] = rows[i]->get(j);
            r[i] = new Row(row, line_size);
        }
        auto mat = Matrix(r, line_size);
        delete[] r;
        swap_count = 0;
        if (gauss_for_det(mat, 0, 0, line_size - 1) == 1){
            double determinator = 1;
            for(int i=0; i < line_size; i++)
                determinator *= mat[i][i];
            return swap_count % 2 == 0 ? determinator : -determinator;
        }
        return 0;
    }

    void operator+=(double value){
        for(int i=0; i < number_of_lines; i++)
            *rows[i] += value;
    }
    void operator+=(const Matrix &matrix){
        for(int i=0; i < number_of_lines; i++)
            *rows[i] += *matrix.rows[i];
    }
    void operator-=(double value){
        for(int i=0; i < number_of_lines; i++)
            *rows[i] -= value;
    }
    void operator-=(const Matrix &matrix){
        for(int i=0; i < number_of_lines; i++)
            *rows[i] -= *matrix.rows[i];
    }
    void operator*=(double value){
        for(int i=0; i < number_of_lines; i++)
            *rows[i] *= value;
    }
    void operator*=(const Matrix &matrix){
        assert(line_size == matrix.number_of_lines);
        for(int i = 0; i < number_of_lines; i++){
            auto *row = new double [matrix.line_size];
            for(int j=0; j < matrix.line_size; j++)
                for(int k=0; k < line_size; k++)
                    row[j] += rows[i]->get(k) * matrix.rows[k]->get(j);
            rows[i] = new Row(row, matrix.line_size);
        }
        line_size = matrix.line_size;
    }
    void operator/=(double value){
        for(int i=0; i < number_of_lines; i++)
            *rows[i] /= value;
    }

    Row& operator[](int index);
    friend Matrix operator+(const Matrix &matrix, double value);
    friend Matrix operator+(const Matrix &matrix1, const Matrix &matrix2);
    friend Matrix operator-(const Matrix &matrix, double value);
    friend Matrix operator-(const Matrix &matrix1, const Matrix &matrix2);
    friend Matrix operator*(const Matrix &matrix, double value);
    friend Matrix operator*(const Matrix &matrix1, const Matrix &matrix2);
    friend Matrix operator/(const Matrix &matrix, double value);
    friend Matrix operator!(const Matrix &matrix); //reverse
    friend Matrix operator~(const Matrix &matrix); //transpose
};
Row& Matrix::operator[](int index) {
    return *rows[index];
}

int straight_find_first_nonzero(Matrix init_matrix, Matrix rev_matrix, int curr_row, int curr_elem, int size){
    int n = curr_row+1;
    while (n < size && init_matrix[n][curr_elem] == 0)
        n++;
    if (n == size && init_matrix[n][curr_elem] == 0)
        return 0;
    init_matrix.swap_rows(n, curr_row);
    rev_matrix.swap_rows(n, curr_row);
    return 1;
}
int reverse_find_first_nonzero(Matrix init_matrix, Matrix rev_matrix, int curr_row, int curr_elem){
    int n = curr_row - 1;
    while(n > -1 && init_matrix[n][curr_elem] == 0)
        n--;
    if (n == -1 && init_matrix[n][curr_elem] == 0)
        return 0;
    init_matrix.swap_rows(curr_row, n);
    rev_matrix.swap_rows(curr_row, n);
    return 1;
}
int straight_gauss(Matrix init_matrix, Matrix rev_matrix, int curr_row, int curr_elem, int size){
    if (curr_row == size + 1 || curr_elem == size + 1)
        return 0;
    if (init_matrix[curr_row][curr_elem] == 0 && straight_find_first_nonzero(init_matrix, rev_matrix, curr_row, curr_elem, size) == 0)
        return straight_gauss(init_matrix, rev_matrix, curr_row, curr_elem + 1, size);
    for(int i = curr_row + 1; i <= size; i++){
        if (init_matrix[i][curr_elem] != 0) {
            double ratio = init_matrix[i][curr_elem] / init_matrix[curr_row][curr_elem];
            //init_matrix.row(i) -= init_matrix.row(curr_row) * ratio;
            //rev_matrix.row(i) -= rev_matrix.row(curr_row) * ratio;
            init_matrix[i] -= init_matrix[curr_row] * ratio;
            rev_matrix[i] -= rev_matrix[curr_row] * ratio;
        }
    }
    return straight_gauss(init_matrix, rev_matrix, curr_row + 1, curr_elem + 1, size);
}
int reverse_gauss(Matrix init_matrix, Matrix rev_matrix, int curr_row, int curr_elem, int size){
    if (curr_row == -1 || curr_elem == -1)
        return 0;
    if (init_matrix[curr_row][curr_elem] == 0 && reverse_find_first_nonzero(init_matrix, rev_matrix, curr_row, curr_elem) == 0)
        return reverse_gauss(init_matrix, rev_matrix, curr_row, curr_elem - 1, size);
    for(int i = curr_row - 1; i >= 0; i--){
        if (init_matrix[i][curr_elem] != 0){
            double ratio = init_matrix[i][curr_elem] / init_matrix[curr_row][curr_elem];
            init_matrix[i] -= init_matrix[curr_row] * ratio;
            rev_matrix[i] -= rev_matrix[curr_row] * ratio;
        }
    }
    return  reverse_gauss(init_matrix, rev_matrix, curr_row - 1, curr_elem - 1, size);
}

Matrix operator+(const Matrix &matrix1, double value){
    int n = matrix1.line_size;
    int m = matrix1.number_of_lines;
    auto **new_rows = new Row*[m];
    auto *row = new double[n];
    for(int i = 0; i < m; i++) {
        for(int j=0; j < n; j++)
            row[j] = matrix1.rows[i]->get(j) + value;
        new_rows[i] = new Row(row, n);
    }
    Matrix new_matrix = Matrix(new_rows, m, n);
    delete[] new_rows;
    delete[] row;
    return new_matrix;
}
Matrix operator+(const Matrix &matrix1, const Matrix &matrix2){
    assert(matrix1.number_of_lines == matrix2.number_of_lines && matrix1.line_size == matrix2.line_size);
    int n = matrix1.line_size;
    int m = matrix1.number_of_lines;
    auto **new_rows = new Row*[m];
    auto *row = new double[n];
    for(int i = 0; i < m; i++) {
        for(int j=0; j < n; j++)
            row[j] = matrix1.rows[i]->get(j) + matrix2.rows[i]->get(j) ;
         new_rows[i] = new Row(row, n);
    }
    Matrix new_matrix = Matrix(new_rows, m, n);
    delete[] new_rows;
    delete[] row;
    return new_matrix;
}
Matrix operator-(const Matrix &matrix1, double value){
    int n = matrix1.line_size;
    int m = matrix1.number_of_lines;
    auto **new_rows = new Row*[m];
    auto *row = new double[n];
    for(int i = 0; i < m; i++) {
        for(int j=0; j < n; j++)
            row[j] = matrix1.rows[i]->get(j) - value;
        new_rows[i] = new Row(row, n);
    }
    Matrix new_matrix = Matrix(new_rows, m, n);
    delete[] new_rows;
    delete[] row;
    return new_matrix;
}
Matrix operator-(const Matrix &matrix1, const Matrix &matrix2){
    assert(matrix1.number_of_lines == matrix2.number_of_lines && matrix1.line_size == matrix2.line_size);
    int n = matrix1.line_size;
    int m = matrix1.number_of_lines;
    auto **new_rows = new Row*[m];
    auto *row = new double[n];
    for(int i = 0; i < m; i++) {
        for(int j=0; j < n; j++)
            row[j] = matrix1.rows[i]->get(j) - matrix2.rows[i]->get(j) ;
        new_rows[i] = new Row(row, n);
    }
    Matrix new_matrix = Matrix(new_rows, m, n);
    delete[] new_rows;
    delete[] row;
    return new_matrix;
}
Matrix operator*(const Matrix &matrix1, double value){
    int n = matrix1.line_size;
    int m = matrix1.number_of_lines;
    auto **new_rows = new Row*[m];
    auto *row = new double[n];
    for(int i = 0; i < m; i++) {
        for(int j=0; j < n; j++)
            row[j] = matrix1.rows[i]->get(j) * value;
        new_rows[i] = new Row(row, n);
    }
    Matrix new_matrix = Matrix(new_rows, m, n);
    delete[] new_rows;
    delete[] row;
    return new_matrix;
}
Matrix operator*(const Matrix &matrix1, const Matrix &matrix2){
    assert(matrix1.line_size == matrix2.number_of_lines);

    //matrix1.number_of_lines x matrix2.line_size
    int m1 = matrix1.number_of_lines;
    int n1 = matrix1.line_size;
    int m2 = matrix2.number_of_lines;
    int n2 = matrix2.line_size;

    auto **new_rows = new Row*[m1];
    for(int i = 0; i < m1; i++){
        auto *row = new double[n2];
        for(int j = 0; j < n2; j++)
            for(int k=0; k < n1; k++)
                row[j] += matrix1.rows[i]->get(k) * matrix2.rows[k]->get(j);
        new_rows[i] = new Row(row, n2);
        //leak mb
    }
    Matrix new_matrix = Matrix(new_rows, m1, n2);
    delete[] new_rows;
    return new_matrix;
}
Matrix operator/(const Matrix &matrix1, double value){
    int n = matrix1.line_size;
    int m = matrix1.number_of_lines;
    auto **new_rows = new Row*[m];
    auto *row = new double[n];
    for(int i = 0; i < m; i++) {
        for(int j=0; j < n; j++)
            row[j] = matrix1.rows[i]->get(j) / value;
        new_rows[i] = new Row(row, n);
    }
    Matrix new_matrix = Matrix(new_rows, m, n);
    delete[] new_rows;
    delete[] row;
    return new_matrix;
}
Matrix operator~(const Matrix &matrix){
    int m = matrix.number_of_lines;
    int n = matrix.line_size;
    Row **new_rows = new Row*[m];
    for(int i = 0; i < m; i++){
        auto *row = new double[n];
        for(int j=0; j < n; j++)
            row[j] = matrix.rows[j]->get(i);
        new_rows[i] = new Row(row, n);
        delete[] row;
    }
    Matrix new_matrix = Matrix(new_rows, m, n);
    delete[] new_rows;
    return new_matrix;
}
Matrix operator!(const Matrix &matrix){
    assert(matrix.line_size == matrix.number_of_lines);
    int n = matrix.number_of_lines;
    Row **new_rows = new Row*[n];
    for(int i = 0; i < n; i++){
        auto *row = new double[n];
        for(int j=0; j < n; j++)
            row[j] = 0;
        row[i] = 1;
        new_rows[i] = new Row(row, n);
    }
    Matrix rev_mat = Matrix(new_rows, n);
    Matrix init_mat = Matrix(matrix.rows, n);
    straight_gauss(init_mat, rev_mat, 0, 0, n - 1);
    reverse_gauss(init_mat, rev_mat, n - 1, n - 1, n - 1);
    for(int i=0; i < n; i++)
        rev_mat[i] /= init_mat[i][i];
    return rev_mat;
} //reverse


Matrix scan_matrix(int n){
    auto *row = new double[n];
    auto **rows = new Row*[n];
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++)
            std::cin >> row[j];
        rows[i] = new Row(row, n);
    }
    Matrix mat = Matrix(rows, n);
    delete[] rows;
    delete[] row;
    return mat;
}

Matrix scan_matrix(int m, int n){
    auto **rows = new Row*[m];
    auto *row = new double[n];
    for(int i = 0; i < m; i++){
        for(int j = 0; j < n; j++)
            std::cin >> row[j];
        rows[i] = new Row(row, n);
    }
    Matrix mat = Matrix(rows, m, n);
    delete[] rows;
    delete[] row;
    return mat;
}

int main() {
    std::cout.precision(3);
    int m, n;
    std::cout << "Enter size for A: ";
    std::cin >> n;
    std::cout << "\nEnter matrix A:\n";
    Matrix A = scan_matrix(n);
    std::cout << "\nEnter number of rows for B: ";
    std::cin >> m;
    std::cout << "\nEnter number of cols for B: ";
    std::cin >> n;
    std::cout << "\nEnter matrix B:\n";
    Matrix B = scan_matrix(m, n);
    std::cout << "\nAnswer:\n";
    (!A * B).print();
}
