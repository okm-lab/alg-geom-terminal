#include <iostream>
#include <vector>

class Vector{
private:
    int dims; // number of coordinates
public:
    std::vector<double> coords;
    Vector(const double* coordinates, int dimensions){
        for(int i=0; i < dimensions; i++)
            coords.push_back(coordinates[i]);
        dims = dimensions;
    }
    //This function outputs the coordinates of the vector
    void print_coords(){
        for(double n: coords)
            std::cout << n << " ";
    }

    //This function multiplies a vector by a scalar
    void mult_by_scalar(double scalar){
        for(int i=0; i < dims; i++)
            coords[i] *= scalar;
    }

    //Function for finding the difference of vectors
    void subtract_coords(std::vector<double> minus_vect){
        for(int i=0; i < dims; i++)
            coords[i] -= minus_vect[i];
    }
};

//The function finds the next vector with nonzero coordinate (current_coord)
int find_next_nonzero(Vector** vectors, int current_vector, int current_coord, int lastvector){
    int n = current_vector + 1;
    while (n < lastvector && vectors[n]->coords[current_coord] == 0)
        n++;
    if (n == lastvector)
        return 0;
    else{
        std::swap(vectors[n], vectors[current_vector]);
        return 1;
    }
}

int gauss(Vector** vectors, int current_vector, int current_coord, int last_vector, int last_coord){
    //if the last coordinate is 0, then the rank is determined
    if (current_coord == last_coord)
        return current_vector;

    if (current_vector == last_vector)
        return current_vector + 1;

    /*
        if the current coordinate of the vector is 0:
        check - if this coordinate for all vectors is 0:
            skip current coordinate
        otherwise:
            replace the rows of the matrix and continue the algorithm
     */
    if (vectors[current_vector]->coords[current_coord] == 0){
        if (find_next_nonzero(vectors, current_vector, current_coord, last_vector) == 0) {
            return gauss(vectors, current_vector, current_coord + 1, last_vector, last_coord);
        }
    }

    for(int i = current_vector + 1; i <= last_vector; i++){
        if (vectors[i]->coords[current_coord] == 0)
            continue;
        double ratio = vectors[current_vector]->coords[current_coord] / vectors[i]->coords[current_coord];
        vectors[i]->mult_by_scalar(ratio);
        vectors[i]->subtract_coords(vectors[current_vector]->coords);
    }
    return gauss(vectors, current_vector + 1, current_coord + 1, last_vector, last_coord);
}


int main(){
    int dimension, num_of_vectors;
    std::cin >> dimension >> num_of_vectors;
    auto **vectors = new Vector*[num_of_vectors];
    auto *coords = new double[dimension];
    for(int i=0; i < num_of_vectors; i++){
        std::cout << "Enter " << i+1 << " vector: ";
        for(int j = 0; j < dimension; j++)
            std::cin >> coords[j];
        vectors[i] = new Vector(coords, dimension);
    }
    int rank = gauss(vectors, 0, 0, num_of_vectors - 1, dimension);
    std::cout << "Basis:\n";
    for(int i=0; i < rank; i++) {
        vectors[i]->print_coords();
        std::cout << "\n";
    }
    std::cout << "Rank: " << rank << "\n";
    if (rank == num_of_vectors)
        std::cout << "This vector system is not linearly dependent";
    else
        std::cout << "This vector system is linearly dependent";
    free(coords);
}