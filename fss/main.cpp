#include <iostream>
#include <cassert>
#include <cmath>

class Row{
private:
    double* numbers;
public:
    int size;
    explicit Row(const double* initial_numbers, int elements_number){
        assert(elements_number >= 0);
        size = elements_number;
        numbers = new double[size];
        for(int i=0; i < size; i++)
            numbers[i] = initial_numbers[i];
    }
    explicit Row(double value, int elements_number){
        assert(elements_number >= 0);
        size = elements_number;
        numbers = new double[size];
        for(int i=0; i < size; i++)
            numbers[i] = value;
    }
    explicit Row(int elements_number){
        assert(elements_number >= 0);
        size = elements_number;
        numbers = new double[size];
        for(int i=0; i < size; i++)
            numbers[i] = 0;

    }
    Row(const Row &row){
        size = row.size;
        numbers = new double[size];
        for (int i = 0; i < size; i++)
            numbers[i] = row.numbers[i];
    }
    void print(){
        for(int i=0; i < size; i++)
            printf("%.3lf   ", numbers[i]);
    }
    double& get(int j){
        assert(j >= 0 && j < size);
        return numbers[j];
    }
    void insert(double elem, int idx){
        assert(idx >= 0 && idx <= size);
        auto *new_row = new double[size + 1];
        for (int i = 0; i < idx; i++)
            new_row[i] = numbers[i];
        new_row[idx] = elem;
        for(int i = idx + 1; i < size + 1; i++)
            new_row[i] = numbers[i - 1];
        delete[] numbers;
        numbers = new_row;
        size++;
    }
    void push_back(double elem){
        insert(elem, size);
    }
    int isNull(){
        for(int i=0; i < size; i++)
            if (fabs(numbers[i]) > 1E-6)
                return 0;
        return 1;
    }

    Row& operator=(const Row &row);
    Row operator+(const Row &row) const;
    Row operator+(double value) const;
    Row operator-(const Row &row) const;
    Row operator-(double value) const;
    Row operator*(double value) const;
    Row operator/(double value) const;
    Row& operator+=(double value);
    Row& operator+=(const Row &addition_row);
    Row& operator-=(double value);
    Row& operator-=(const Row &subtraction_row);
    Row& operator*=(double value);
    Row& operator/=(double value);
    double& operator[](int index);
};
double& Row::operator[](int index){
    assert(index >= 0 && index < size);
    return numbers[index];
}

Row& Row::operator+=(double value){
    for(int i = 0; i < size; i++)
        numbers[i] += value;
    return *this;
}
Row& Row::operator+=(const Row &addition_row){
    assert(size == addition_row.size);
    for(int i = 0; i < size; i++)
        numbers[i] += addition_row.numbers[i];
    return *this;
}
Row& Row::operator-=(double value){
    for(int i = 0; i < size; i++)
        numbers[i] -= value;
    return *this;
}
Row& Row::operator-=(const Row &subtraction_row){
    assert(size == subtraction_row.size);
    for(int i = 0; i < size; i++)
        numbers[i] -= subtraction_row.numbers[i];
    return *this;
}
Row& Row::operator*=(double value){
    for(int i = 0; i < size; i++)
        numbers[i] *= value;
    return *this;
}
Row& Row::operator/=(double value){
    assert(value != 0);
    for(int i = 0; i < size; i++)
        numbers[i] /= value;
    return *this;
}
Row Row::operator+(double value) const {
    return  Row(*this) += value;
}
Row Row::operator+(const Row &row) const {
    return Row(*this) += row;
}
Row Row::operator-(double value) const {
    return Row(*this) -= value;
}
Row Row::operator-(const Row &row) const {
    return Row(*this) -= row;
}
Row Row::operator*(double value) const {
    return Row(*this) *= value;
}
Row Row::operator/(double value) const {
    return Row(*this) /= value;
}
Row &Row::operator=(const Row &row) {
    numbers = new double[row.size];
    for(int i=0; i < row.size; i++)
        numbers[i] = row.numbers[i];
    return *this;
}

class Matrix{
private:
    Row** rows;
    int swap_count = 0;
    static int find_next_nonzero(Matrix matrix, int current_row, int current_element, int size){
    int n = current_row + 1;
    while (n < size && matrix[n][current_element] == 0)
        n++;
    if (n >= size || matrix[n][current_element] == 0)
        return 0; //stop gauss, answer is 0
    else {
        matrix.swap_rows(current_row, n);
        matrix.swap_count++;
        return 1; //success, continue gauss
    }
}
    int gauss_for_det(Matrix &matrix, int current_row, int current_element, int m_size){
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
    int gauss_for_rank(Matrix &matrix, int current_row, int current_element, int row_size, int line_size) {
        if (current_element == line_size + 1 || current_row == row_size + 1)
            return current_row;
        if (matrix[current_row][current_element] == 0 &&
            find_next_nonzero(matrix, current_row, current_element, row_size) == 0)
            return gauss_for_rank(matrix, current_row, current_element + 1, row_size, line_size);
        for (int i = current_row + 1; i <= row_size; i++) {
            if (matrix[i][current_element] != 0) {
                double ratio = matrix[i][current_element] / matrix[current_row][current_element];
                matrix[i] -= matrix[current_row] * ratio;
            }
        }
        return gauss_for_rank(matrix, current_row + 1, current_element + 1, row_size, line_size);
    }
public:
    int number_of_lines;
    int line_size;
    explicit Matrix(Row** mat, int m, int n){
        assert(m >= 0);
        assert(n >= 0);
        rows = new Row*[m];
        for(int i=0; i < m; i++)
            rows[i] = mat[i];
        line_size = n;
        number_of_lines = m;
    }
    explicit Matrix(Row** mat, int n){
        assert(n >= 0);
        rows = new Row*[n];
        for(int i=0; i < n; i++)
            rows[i] = mat[i];
        line_size = n;
        number_of_lines = n;
    }
    explicit Matrix(int m, int n, double value){
        assert(m >= 0);
        assert(n >= 0);
        rows = new Row*[m];
        for(int i=0; i < n; i++)
            rows[i] = new Row(value, n);
        line_size = n;
        number_of_lines = m;
    }
    explicit Matrix (int m, int n){
        assert(m >= 0);
        assert(n >= 0);
        rows = new Row*[m];
        for(int i=0; i < n; i++)
            rows[i] = new Row(n);
        line_size = n;
        number_of_lines = m;
    }
    explicit Matrix (int n){
        assert(n >= 0);
        rows = new Row*[n];
        for(int i=0; i < n; i++)
            rows[i] = new Row(n);
        line_size = n;
        number_of_lines = n;
    }
    Matrix(const Matrix &matrix){
        number_of_lines = matrix.number_of_lines;
        line_size = matrix.line_size;
        rows = new Row*[number_of_lines];
        for (int i = 0; i < number_of_lines; i++)
            rows[i] = new Row(*matrix.rows[i]);
    }
    int isNull(){
        for(int i=0; i < number_of_lines; i++)
            if (rows[i]->isNull() == 0)
                return 0;
        return 1;
    }
    void print(){
        for(int i=0; i < number_of_lines; i++) {
            rows[i]->print();
            std::cout << "\n";
        }
    }
    void swap_rows(int first_row, int second_row){
        assert(first_row >= 0 && first_row < number_of_lines);
        assert(second_row >= 0 && second_row < number_of_lines);
        std::swap(rows[first_row], rows[second_row]);
    }
    void add_row(Row *row, int idx){
        assert(idx >= 0 && idx <= number_of_lines);
        assert(row->size == line_size);
        Row** new_rows = new Row*[number_of_lines+1];
        for(int i=0; i < idx; i++)
            new_rows[i] = rows[i];
        new_rows[idx] = row;
        for(int i=idx + 1; i < number_of_lines + 1; i++)
            new_rows[i] = rows[i - 1];
        delete[] rows;
        rows = new_rows;
        number_of_lines++;
    }
    void add_row(Row *row){
        add_row(row, number_of_lines);
    }
    void delete_row(int idx){
        assert(idx >= 0 && idx < number_of_lines);
        for (int i = idx + 1; i < number_of_lines; i++)
            rows[i - 1] = rows[i];
        number_of_lines--;
    }
    void delete_column(int idx){
        assert(idx >= 0 && idx < line_size);
        for (int i = 0; i < number_of_lines; i++) {
            for (int j = idx + 1; j < line_size; j++)
                rows[i]->get(j - 1) = rows[i]->get(j);
            rows[i]->size--;
        }
        line_size--;
    }
    void add_column(Matrix column, int idx){
        assert(idx >= 0 && idx <= line_size);
        assert(column.number_of_lines == number_of_lines);
        for(int i=0; i < number_of_lines; i++)
            rows[i]->insert(column[i].get(0), idx);
        line_size++;
    }
    void add_column(Matrix column){
        for(int i=0; i < number_of_lines; i++)
            rows[i]->insert(column[i].get(0), line_size);
        line_size++;
    }
    void change_row(Row *row, int idx){
        assert(idx >= 0 && idx < number_of_lines);
        rows[idx] = row;
    }
    void change_column(Matrix column, int idx){
        assert(idx >= 0 && idx < line_size);
        for(int i=0; i < number_of_lines; i++)
            rows[i]->get(idx) = column[i][0];
    }
    void extend_matrix(Matrix matrix){
        assert(matrix.number_of_lines == number_of_lines);
        for(int i=0; i < matrix.line_size; i++)
            add_column(matrix.get_column(i));
    }
    void delete_null_rows(){
        int chk = 0;
        while(chk == 0) {
            chk = 1;
            for (int i = 0; i < number_of_lines; i++) {
                if (rows[i]->isNull()) {
                    delete_row(i);
                    chk = 0;
                }
            }
        }
    }
    void add_column(){
        for(int i=0; i < number_of_lines; i++)
            rows[i]->insert(0, line_size);
        line_size++;
    }
    void delete_null_columns(){
        for(int i = line_size - 1; i >= 0; i--) {
            int chk = 0;
            for (int j = 0; j < number_of_lines; j++) {
                if (rows[j]->get(i) != 0) {
                    chk = 1;
                    break;
                }
            }
            if (chk == 0) {
                delete_column(i);
                delete_column(i);
            }
        }
    }
    void add_to_column(double value, int idx){
        for(int i=0; i < number_of_lines; i++)
            rows[idx]->get(i) += value;
    }
    void add_to_column(const Matrix &matrix, int idx){
        for(int i=0; i < number_of_lines; i++)
            rows[i]->get(idx) += matrix.rows[i]->get(0);
    }
    void delete_kr(){
        for(int i=0; i < number_of_lines - 1; i++){
            for(int j = i + 1; j < number_of_lines; j++){
                if (rows[i]->isNull() || rows[j]->isNull())
                    continue;
                double t = INT32_MAX;
                for(int k=0; k < line_size; k++){
                    if ((fabs(rows[j]->get(k)) < 1E-6 && fabs(rows[i]->get(k)) > 1E-6) ||
                        (fabs(rows[i]->get(k)) < 1E-6 && fabs(rows[j]->get(k)) > 1E-6))
                        break;
                    if (fabs(rows[i]->get(k)) < 1E-6 && fabs(rows[i]->get(k)) < 1E-6)
                        continue;
                    if (t != INT32_MAX){
                        if (rows[j]->get(k) / rows[i]->get(k) == t)
                            continue;
                        else{
                            break;
                        }
                    }else{
                        t = rows[j]->get(k) / rows[i]->get(k);
                    }
                    *rows[j] -= *rows[i] * t;
                }
            }
        }
    }
    Matrix get_column(int idx){
        assert(idx >= 0 && idx < line_size);
        auto **new_rows = new Row*[number_of_lines];
        auto *row = new double[1];
        for(int i=0; i < number_of_lines; i++){
            row[0] = rows[i]->get(idx);
            new_rows[i] = new Row(row, 1);
        }
        Matrix column = Matrix(new_rows, number_of_lines, 1);
        delete[] new_rows;
        return column;
    }
    double determinator(){
        assert(number_of_lines == line_size);
        auto mat = Matrix(*this);
        swap_count = 0;
        if (gauss_for_det(mat, 0, 0, line_size - 1) == 1){
            double determinator = 1;
            for(int i=0; i < line_size; i++)
                determinator *= mat[i][i];
            return swap_count % 2 == 0 ? determinator : -determinator;
        }
        return 0;
    }
    int rank(){
        auto mat = Matrix(*this);
        return gauss_for_rank(mat, 0, 0, number_of_lines-1, line_size-1);
    }
    Row& operator[](int index);
    Matrix& operator=(const Matrix &matrix);
    friend Matrix operator!(const Matrix &matrix);
    friend Matrix operator~(const Matrix &matrix);
    Matrix operator+(Matrix &addition_matrix) const;
    Matrix operator+(double value) const;
    Matrix operator-(Matrix &subtraction_matrix) const;
    Matrix operator-(double value) const;
    Matrix operator*(double value) const;
    Matrix operator*(Matrix &mult_matrix) const;
    Matrix operator/(double value) const;
    Matrix& operator+=(double value);
    Matrix& operator+=(const Matrix &addition_matrix);
    Matrix& operator-=(double value);
    Matrix& operator-=(const Matrix &subtraction_matrix);
    Matrix& operator*=(double value);
    Matrix& operator*=(const Matrix &mult_matrix);
    Matrix& operator/=(double value);

};
Row& Matrix::operator[](int index) {
    assert(index >= 0 && index < number_of_lines);
    return *rows[index];
}

int straight_find_first_nonzero(Matrix &first_matrix, Matrix &second_matrix, int curr_row, int curr_elem, int size){
    int n = curr_row+1;
    while (n < size && first_matrix[n][curr_elem] == 0)
        n++;
    if (n >= size || first_matrix[n][curr_elem] == 0)
        return -1;
    first_matrix.swap_rows(n, curr_row);
    second_matrix.swap_rows(n, curr_row);
    return n;
}
int reverse_find_first_nonzero(Matrix &first_matrix, Matrix &second_matrix, int curr_row, int curr_elem){
    int n = curr_row - 1;
    while(n > -1 && first_matrix[n][curr_elem] == 0)
        n--;
    if (n == -1 || first_matrix[n][curr_elem] == 0)
        return -1;
    first_matrix.swap_rows(curr_row, n);
    second_matrix.swap_rows(curr_row, n);
    return n;
}
int straight_gauss(Matrix &first_matrix, Matrix &second_matrix, int curr_row, int curr_elem, int row_size, int line_size){
    if (curr_row == row_size + 1 || curr_elem == line_size + 1)
        return 0;
    if (first_matrix[curr_row][curr_elem] == 0 && straight_find_first_nonzero(first_matrix, second_matrix, curr_row, curr_elem, row_size) == -1)
        return straight_gauss(first_matrix, second_matrix, curr_row, curr_elem + 1, row_size, line_size);
    for(int i = curr_row + 1; i <= row_size; i++){
        if (first_matrix[i][curr_elem] != 0) {
            double ratio = first_matrix[i][curr_elem] / first_matrix[curr_row][curr_elem];
            first_matrix[i] -= first_matrix[curr_row] * ratio;
            second_matrix[i] -= second_matrix[curr_row] * ratio;
        }
    }
    return straight_gauss(first_matrix, second_matrix, curr_row + 1, curr_elem + 1, row_size, line_size);
}

int reverse_gauss(Matrix &first_matrix, Matrix &second_matrix, int curr_row, int curr_elem, int row_size, int line_size){
    if (curr_row == -1 || curr_elem == -1)
        return 0;
    if (first_matrix[curr_row][curr_elem] == 0 && reverse_find_first_nonzero(first_matrix, second_matrix, curr_row, curr_elem) == -1)
        return reverse_gauss(first_matrix, second_matrix, curr_row, curr_elem - 1, row_size, line_size);
    for(int i = curr_row - 1; i >= 0; i--){
        if (first_matrix[i][curr_elem] != 0){
            double ratio = first_matrix[i][curr_elem] / first_matrix[curr_row][curr_elem];
            first_matrix[i] -= first_matrix[curr_row] * ratio; //   here
            second_matrix[i] -= second_matrix[curr_row] * ratio;
        }
    }
    return reverse_gauss(first_matrix, second_matrix, curr_row - 1, curr_elem - 1, row_size, line_size);
}
int straight_find_first_nonzero(Matrix &matrix, int curr_row, int curr_elem, int size){
    int n = curr_row+1;
    while (n < size && matrix[n][curr_elem] == 0)
        n++;
    if (n >= size || matrix[n][curr_elem] == 0)
        return -1;
    matrix.swap_rows(n, curr_row);
    return n;
}
int reverse_find_first_nonzero(Matrix &matrix, int curr_row, int curr_elem){
    int n = curr_row - 1;
    while(n > -1 && matrix[n][curr_elem] == 0)
        n--;
    if (n == -1 || matrix[n][curr_elem] == 0)
        return -1;
    matrix.swap_rows(curr_row, n);
    return n;
}
int straight_gauss(Matrix &matrix, int curr_row, int curr_elem, int row_size, int line_size){
    //если будут баги с гауссом, они в этой строке
    //matrix.print();
    //std::cout << "\n\n";
    if (curr_row == row_size || curr_elem == line_size)
        return 0;
    if (matrix[curr_row][curr_elem] == 0 && straight_find_first_nonzero(matrix, curr_row, curr_elem, row_size) == -1)
        return straight_gauss(matrix, curr_row, curr_elem + 1, row_size, line_size);
    for(int i = curr_row + 1; i <= row_size; i++){
        if (fabs(matrix[i][curr_elem]) > 1E-6 && fabs(matrix[curr_row][curr_elem]) > 1E-6) {
            double ratio = matrix[i][curr_elem] / matrix[curr_row][curr_elem];
            matrix[i] -= matrix[curr_row] * ratio;
        }
    }
    return straight_gauss(matrix, curr_row + 1, curr_elem + 1, row_size, line_size);
}
int reverse_gauss(Matrix &matrix, int curr_row, int curr_elem, int row_size, int line_size){
    //matrix.print();
    //std::cout << "\n\n";
    if (curr_row == 0 || curr_elem == 0)
        return 0;
    if (matrix[curr_row][curr_elem] == 0 && reverse_find_first_nonzero(matrix, curr_row, curr_elem) == -1)
        return reverse_gauss(matrix, curr_row, curr_elem - 1, row_size, line_size);
    for(int i = curr_row - 1; i >= 0; i--){
        if (fabs(matrix[i][curr_elem]) > 1E-6 && fabs(matrix[curr_row][curr_elem]) > 1E-6){
            double ratio = matrix[i][curr_elem] / matrix[curr_row][curr_elem];
            matrix[i] -= matrix[curr_row] * ratio;
        }
    }
    return  reverse_gauss(matrix, curr_row - 1, curr_elem - 1, row_size, line_size);
}

Matrix &Matrix::operator+=(double value) {
    for(int i = 0;i < number_of_lines; i++)
        *rows[i] += value;
    return *this;
}
Matrix &Matrix::operator+=(const Matrix &addition_matrix){
    assert(number_of_lines == addition_matrix.number_of_lines && line_size == addition_matrix.line_size);
    for(int i=0; i < number_of_lines; i++)
        *rows[i] += *addition_matrix.rows[i];
    return *this;
}
Matrix &Matrix::operator-=(double value) {
    for(int i = 0;i < number_of_lines; i++)
        *rows[i] -= value;
    return *this;
}
Matrix &Matrix::operator-=(const Matrix &subtraction_matrix){
    assert(number_of_lines == subtraction_matrix.number_of_lines && line_size == subtraction_matrix.line_size);
    for(int i=0; i < number_of_lines; i++)
        *rows[i] -= *subtraction_matrix.rows[i];
    return *this;
}
Matrix &Matrix::operator*=(double value) {
    for(int i = 0;i < number_of_lines; i++)
        *rows[i] *= value;
    return *this;
}
Matrix &Matrix::operator*=(const Matrix &mult_matrix){
    assert(line_size == mult_matrix.number_of_lines);
    for(int i = 0; i < number_of_lines; i++){
        auto *row = new double [mult_matrix.line_size];
        for(int j=0; j < mult_matrix.line_size; j++)
            for(int k=0; k < line_size; k++)
                row[j] += rows[i]->get(k) * mult_matrix.rows[k]->get(j);
        rows[i] = new Row(row, mult_matrix.line_size);
    }
    line_size = mult_matrix.line_size;
    return *this;
}
Matrix &Matrix::operator/=(double value) {
    for(int i = 0;i < number_of_lines; i++)
        *rows[i] /= value;
    return *this;
}
Matrix Matrix::operator+(double value) const {
    return Matrix(*this) += value;
}
Matrix Matrix::operator+(Matrix &addition_matrix) const {
    return Matrix(*this) += addition_matrix;
}
Matrix Matrix::operator-(double value) const {
    return Matrix(*this) -= value;
}
Matrix Matrix::operator-(Matrix &subtraction_matrix) const {
    return Matrix(*this) -= subtraction_matrix;
}
Matrix Matrix::operator*(double value) const {
    return Matrix(*this) *= value;
}
Matrix Matrix::operator*(Matrix &mult_matrix) const {
    return Matrix(*this) *= mult_matrix;
}
Matrix Matrix::operator/(double value) const {
    return Matrix(*this) /= value;
}
Matrix &Matrix::operator=(const Matrix &matrix) {
    delete[] rows;
    line_size = matrix.line_size;
    number_of_lines = matrix.number_of_lines;
    rows = new Row*[line_size];
    for (int i = 0; i < number_of_lines; i++)
        rows[i] = new Row(*matrix.rows[i]);
    return *this;
}
Matrix operator~(const Matrix &matrix){
    int m = matrix.number_of_lines;
    int n = matrix.line_size;
    Row **new_rows = new Row*[n];
    for(int i = 0; i < n; i++){
        auto *row = new double[m];
        for(int j=0; j < m; j++)
            row[j] = matrix.rows[j]->get(i);
        new_rows[i] = new Row(row, m);
        delete[] row;
    }
    Matrix new_matrix = Matrix(new_rows, n, m);
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
    auto init_mat = matrix;
    straight_gauss(init_mat, rev_mat, 0, 0, n - 1, n - 1);
    reverse_gauss(init_mat, rev_mat, n - 1, n - 1, n - 1, n - 1);
    for(int i=0; i < n; i++)
        rev_mat[i] /= init_mat[i][i];
    return rev_mat;
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
Matrix scan_matrix(int n){
    return scan_matrix(n, n);
}
Matrix cramer_rule(Matrix A){
    Matrix ans = A.get_column(A.line_size - 1);
    Matrix system = A;
    system.delete_column(system.line_size - 1);
    auto *determinators = new double[system.line_size+1];
    determinators[0] = system.determinator();
    Matrix changed_system = system;
    for(int i=1; i <= system.line_size; i++){
        changed_system.change_column(ans, i - 1);
        determinators[i] = changed_system.determinator();
        changed_system.change_column(system.get_column(i - 1), i - 1);
    }
    for(int i=1 ; i <= ans.number_of_lines; i++)
        ans[i-1][0] = determinators[i] / determinators[0];
    ans.print();
    return ans;
}

void test_mult_matrix(Matrix A, Matrix B){
    std::cout << "\nMatrix mult:\n";
    (A * B).print();
}
void test_mult_value(Matrix A, double value){
    std::cout << "\nvalue mult:\n";
    (A * value).print();
}
void test_add_matrix(Matrix A, Matrix B){
    std::cout << "\nMatrix add:\n";
    (A + B).print();
}
void test_add_value(Matrix A, double value){
    std::cout << "\nValue add:\n";
    (A + value).print();
}
void test_sub_matrix(Matrix A, Matrix B){
    std::cout << "\nMatrix sub:\n";
    (A - B).print();
}
void test_sub_value(Matrix A, double value){
    std::cout << "\nValue sub:\n";
    (A - value).print();
}
void test_divide_value(Matrix A, double value){
    std::cout << "\nValue div:\n";
    (A / value).print();
}
void test_rank(Matrix A){
    std::cout << "\nMatrix rank:\n";
    std::cout << A.rank();
}
void test_determinator(Matrix A){
    std::cout << "\nMatrix det:\n";
    std::cout << A.determinator();
}
void test_reverse(Matrix A){
    assert(A.determinator() != 0);
    std::cout << "\nMatrix reverse:\n";
    (!A).print();
}
void test_trans(Matrix A){
    std::cout << "\nMatrix trans:\n";
    (~A).print();
}
void test_extend(Matrix A, Matrix B){
    std::cout << "\nMatrix extend:\n";
    A.extend_matrix(B);
    A.print();
}
void test_operations(){
    int m1, n1, m2, n2, value;
    std::cin >> m1 >> n1;
    Matrix A = scan_matrix(m1, n1);
    std::cin >> m2 >> n2;
    Matrix B = scan_matrix(m2, n2);
    std::cin >> value;
    std::cout << "\n";
    test_add_value(A, value);
    test_add_matrix(A, B);
    test_sub_value(A, value);
    test_sub_matrix(A, B);
    test_divide_value(A, value);
    test_mult_value(A, value);
    test_mult_matrix(A, B);
    test_rank(A);
    test_determinator(A);
    test_reverse(A);
    test_trans(A);
    test_extend(A, B);
}


Matrix fsr(Matrix A){
    straight_gauss(A, 0, 0, A.number_of_lines - 1, A.line_size - 1);
    reverse_gauss(A, A.number_of_lines - 1, A.line_size - 1, A.number_of_lines-1, A.line_size-1);
    A.delete_kr();
    straight_gauss(A, 0, 0, A.number_of_lines - 1, A.line_size - 1);
    Matrix ans = Matrix(A.line_size, A.line_size);
    for(int i=0; i < A.number_of_lines; i++){
        for(int j = 0; j < A.line_size; j++){
            if (fabs(A[i][j])>1E-6){
                ans[j] = A[i] * -1 / A[i][j];
                for(int k = 0; k < j; k++)
                    ans[k] += ans[j] * ans[k][j];
                ans[j][j] = 0;
                if (ans[j].isNull())
                    ans[j][0] = INT32_MAX;
                break;
            }
        }
    }
    //ans.print();
    return ans;
}

int main() {
    setlocale(LC_ALL, "Russian");
    freopen("input.txt", "r", stdin);
    //test_operations();
    int a,b;
    std::cin >> a >> b;
    Matrix A = scan_matrix(a, b);
    Matrix ans = fsr(A);
    bool chk = false;
    for(int i=0; i < ans.number_of_lines; i++){
        if (ans[i][0] != INT32_MAX)
            chk = true;
    }
    if (chk == false){
        std::cout << "NO";
        exit(0);
    }
    bool *nulls = new bool[ans.line_size];
    std::cout << "Общее решение:\n";
    for(int i=0; i < ans.number_of_lines; i++){
        if (ans[i].isNull()) {
            std::cout << "x" << i + 1 << " is free\n";
            nulls[i] = true;
        }
        else if (ans[i][0] == INT32_MAX){
            std::cout << "x" << i + 1 << " = 0\n";
        }
        else{
            nulls[i] = false;
            std::cout << "x" << i + 1 << " = ";
            int j = i + 1;
            for(; j < ans.line_size; j++){
                if (fabs(ans[i][j]) > 1E-6) {
                    std::cout << ans[i][j] << "x" << j+1;
                    break;
                }
            }
            j++;
            for(;j < ans.line_size; j++){
                if (fabs(ans[i][j]) > 1E-6)
                    std::cout << (ans[i][j] > 0 ? " + " : " - ") << fabs(ans[i][j]) << "x" << j + 1;
            }
            std::cout << "\n";
        }
    }
    std::cout << "\nФСР:\n";
    for(int i=0; i < ans.line_size; i++){
        if (nulls[i] && ans[i][0] != INT32_MAX){
            for(int j = 0; j < ans.line_size; j++)
                std::cout << (j == i ? 1 : fabs(ans[j][i]) > 1E-6 ? ans[j][i] : 0) << " ";
            std::cout << "\n";
        }
    }

}
