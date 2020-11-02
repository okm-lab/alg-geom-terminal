#include <iostream>
#include <cstdlib>

int swap_count = 0;

class Row{
private:
    double* numbers;
    int size;
public:
    Row(const double* initial_numbers, int elements_number){
        numbers = (double*)malloc(sizeof(double) * elements_number);
        for(int i=0; i < elements_number; i++)
            numbers[i] = initial_numbers[i];
        size = elements_number;
    }

    void print(){
        for(int i=0; i < size; i++)
            std::cout << numbers[i] << " ";
    }

    double get_number(int j){
        return numbers[j];
    }

    friend Row operator*(const Row &row1, double value);
    friend Row operator-(const Row &row1, const Row &row2);
    void operator-=(const Row &row){
        for(int i=0; i < size; i++)
            numbers[i] -= row.numbers[i];
    }
    void operator*=(double value){
        for(int i=0; i < size; i++)
            numbers[i] *= value;
    }
};


Row operator-(const Row &row1, const Row &row2){
    auto *new_numbers = (double*)malloc(sizeof(double)*row1.size);
    for(int i=0; i < row1.size; i++)
        new_numbers[i] = row1.numbers[i] - row2.numbers[i];
    return Row(new_numbers, row1.size);
}

Row operator*(const Row &row, double value){
    auto *new_numbers = (double*)malloc(sizeof(double)*row.size);
    for(int i=0; i < row.size; i++)
        new_numbers[i] = row.numbers[i] * value;
    return Row(new_numbers, row.size);
}

class Matrix{
private:
    Row** rows;
    int size;
public:
    Matrix(Row** m, int n){
        size = n;
        rows = m;
    }

    void print(){
        for(int i=0; i < size; i++){
            rows[i]->print();
            std::cout << "\n";
        }
    }

    double get_element(int i, int j){
        return rows[i]->get_number(j);
    }

    Row row(int i){
        return *rows[i];
    }

    void swap_rows(int i, int j){
        std::swap(rows[i], rows[j]);
    }

};


//The function finds the next row with nonzero coordinate (current_element)
int find_next_nonzero(Matrix matrix, int current_row, int current_element, int size){
    int n = current_row + 1;
    while (n < size && matrix.get_element(n, current_element) == 0)
        n++;
    if (n == size && matrix.get_element(n, current_element) == 0)
        return 0; //stop gauss, answer is 0
    else {
        matrix.swap_rows(current_row, n);
        swap_count++;
        return 1; //success, continue gauss
    }
}

int gauss_for_det(Matrix matrix, int current_row, int current_element, int size){
    //return 1 - success, we can determine determinator
    //return 0 - answer is 0, there is a zero element on the diagonal
    if (current_row == size && matrix.get_element(current_row, current_element) != 0)
        return 1;
    if (current_row == size)
        return 0;

    /*
        if the current number of the matrix is 0:
        check - if this number is 0 for all rows:
            return 0
        otherwise:
            replace the rows of the matrix and continue the algorithm
     */
    if (matrix.get_element(current_row, current_element) == 0)
        if (find_next_nonzero(matrix, current_row, current_element, size) == 0)
            return 0;

    for(int i= current_row + 1; i <= size; i++){
        if (matrix.get_element(i, current_element) != 0) {
            double ratio = matrix.get_element(i, current_element) / matrix.get_element(current_row, current_element);
            matrix.row(i) -= matrix.row(current_row) * ratio;
        }
    }
    return gauss_for_det(matrix, current_row + 1, current_element + 1, size);
}

double determine_determinator(Matrix matrix, int size){
    if (gauss_for_det(matrix, 0, 0, size - 1) == 1){
        double determinator = 1;
        for(int i=0; i < size; i++)
            determinator *= matrix.get_element(i, i);
        return swap_count % 2 == 0 ? determinator : -determinator;
    }
    return 0;
}


int main() {
    int n;
    std::cout << "Enter matrix size: ";
    std::cin >> n;
    std::cout << "Enter matrix:\n";
    auto *row = new double[n];
    auto **rows = new Row*[n];
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++)
            std::cin >> row[j];
        rows[i] = new Row(row, n);
    }
    Matrix matrix = Matrix(rows, n);
    std::cout << "Determinator is: " << determine_determinator(matrix, n);
    free(row);
    free(rows);
    return 0;
}
