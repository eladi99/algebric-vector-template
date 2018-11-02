#include <cmath>
#include <string>
#include <sstream>

template <unsigned int dim>
class Vector {
    double *_arr;
    static constexpr double epsilon = 0.001;

public:
    explicit Vector (double arr[] = nullptr);
    Vector (const Vector& vec);
    Vector (const Vector&& vec);
    ~Vector ();
    
    Vector (double) = delete;
    operator double () = delete;

    static Vector zero_vec ();
    static Vector unit_vec (unsigned int ind = 0);
    
    unsigned int get_dim () const;
    
    Vector add (const Vector& vec) const;
    Vector subtract (const Vector& vec) const;
    Vector multiply_scalar (const double alpha) const;
    double dot_product (const Vector& vec) const;
    bool is_equal (const Vector& vec) const;
    
    double norm () const;
    Vector normalize () const;
    
    double& operator[] (unsigned int ind);
    const double& operator[] (unsigned int ind) const;
    Vector& operator= (const Vector& vec);
    Vector operator-() const;

    static double angle (const Vector& vec1, const Vector& vec2);
    Vector operator+ (const Vector& vec) const;
    Vector operator- (const Vector& vec) const;
    Vector operator* (double alpha) const;
    double operator* (const Vector& vec) const;
    Vector& operator+= (const Vector& vec);
    Vector& operator-= (const Vector& vec);
    Vector& operator*= (double alpha);
    bool operator== (const Vector& vec) const;
    bool operator!= (const Vector& vec) const;
    
    std::string to_string () const;
};

template <unsigned int dim>
Vector<dim>::Vector (double arr[]): _arr (new double[dim]) {
    if (arr == nullptr) {
        for (unsigned int i = 0; i < dim; i++) {
            _arr[i] = 0.0;
            return;
        }
    }
    
    for (unsigned int i = 0; i < dim; i++) {
        _arr[i] = arr[i];
    }
}

template <unsigned int dim>
Vector<dim>::Vector (const Vector& vec): _arr (new double[dim]) {
    for (unsigned int i = 0; i < dim; i++) {
        _arr[i] = vec._arr[i];
    }
}

template <unsigned int dim>
Vector<dim>::Vector (const Vector&& vec): _arr (new double[dim]) {
    for (unsigned int i = 0; i < dim; i++) {
        _arr[i] = vec._arr[i];
    }
}

template <unsigned int dim>
Vector<dim>::~Vector () {
    delete[] _arr;
}

template <unsigned int dim>
Vector<dim> Vector<dim>::zero_vec () {
    return Vector ();
}

template <unsigned int dim>
Vector<dim> Vector<dim>::unit_vec (unsigned int ind) {
    Vector vec0 = Vector ();
    vec0._arr[ind] = 1.0;
    return vec0;
}

template <unsigned int dim>
unsigned int Vector<dim>::get_dim () const {
    return dim;
}

template <unsigned int dim>
Vector<dim> Vector<dim>::add (const Vector& vec) const {
    Vector sum_vec = Vector ();
    for (unsigned int i = 0; i < dim; i++) {
        sum_vec._arr[i] = _arr[i] + vec._arr[i];
    }

    return sum_vec;
}

template <unsigned int dim>
Vector<dim> Vector<dim>::subtract (const Vector& vec) const {
    Vector sum_vec = Vector ();
    for (unsigned int i = 0; i < dim; i++) {
        sum_vec._arr[i] = _arr[i] - vec._arr[i];
    }

    return sum_vec;
}

template <unsigned int dim>
Vector<dim> Vector<dim>::multiply_scalar (const double alpha) const {
    Vector mult_vec (*this);
    for (unsigned int i = 0; i < dim; i++) {
        mult_vec._arr[i] *= alpha;
    }

    return mult_vec;
}

template <unsigned int dim>
double Vector<dim>::dot_product (const Vector& vec) const {
    double mult = 0;
    for (unsigned int i = 0; i < dim; i++) {
        mult += _arr[i] * vec._arr[i];
    }

    return mult;
}

template <unsigned int dim>
double Vector<dim>::norm () const {
    return sqrt (dot_product (*this));
}

template <unsigned int dim>
Vector<dim> Vector<dim>::normalize () const {
    return Vector (*this).multiply_scalar (1 / norm ());
}

template <unsigned int dim>
bool Vector<dim>::is_equal (const Vector& vec) const {
    for (unsigned int i = 0; i < dim; i++) {
        if (abs (_arr[i] - vec._arr[i]) > epsilon) {
            return false;
        }
    }

    return true;
}

template <unsigned int dim>
Vector<dim>& Vector<dim>::operator= (const Vector& vec) {
    if (this == &vec) {
        return *this;
    }
    
    delete[] _arr;
    _arr = new double[dim];
    for (unsigned int i = 0; i < dim; i++) {
        _arr[i] = vec._arr[i];
    }
}

template <unsigned int dim>
Vector<dim> Vector<dim>::operator- () const {
    return multiply_scalar (-1.0);
}

template <unsigned int dim>
double& Vector<dim>::operator[] (unsigned int ind) {
    return _arr[ind];
}

template <unsigned int dim>
const double& Vector<dim>::operator[] (unsigned int ind) const {
    return _arr[ind];
}

template <unsigned int dim>
double Vector<dim>::angle (const Vector& vec1, const Vector& vec2) {
    return acos (vec1.normalize () * vec2.normalize ());
}

template <unsigned int dim>
Vector<dim> Vector<dim>::operator+ (const Vector& vec) const {
    return add (vec);
}

template <unsigned int dim>
Vector<dim> Vector<dim>::operator- (const Vector& vec) const {
    return subtract (vec);
}

template <unsigned int dim>
Vector<dim> Vector<dim>::operator* (double alpha) const {
    return multiply_scalar (alpha);
}

template <unsigned int dim>
double Vector<dim>::operator* (const Vector& vec) const {
    return dot_product (vec);
}

template <unsigned int dim>
Vector<dim>& Vector<dim>::operator+= (const Vector& vec) {
    for (unsigned int i = 0; i < dim; i++) {
        _arr[i] += vec._arr[i];
    }

    return *this;
}

template <unsigned int dim>
Vector<dim>& Vector<dim>::operator-= (const Vector& vec) {
    for (unsigned int i = 0; i < dim; i++) {
        _arr[i] -= vec._arr[i];
    }

    return *this;
}

template <unsigned int dim>
Vector<dim>& Vector<dim>::operator*= (double alpha) {
    for (unsigned int i = 0; i < dim; i++) {
        _arr[i] *= alpha;
    }

    return *this;
}

template <unsigned int dim>
bool Vector<dim>::operator== (const Vector& vec) const {
    return this == &vec || is_equal (vec);
}

template <unsigned int dim>
bool Vector<dim>::operator!= (const Vector& vec) const {
    return !is_equal (vec);

}

template <unsigned int dim>
std::string Vector<dim>::to_string () const {
    std::ostringstream out_str;
    out_str << "(" << _arr[0];
    for (unsigned int i = 1; i < dim; i++) {
        out_str << ", " << _arr[i];
    }

    out_str << ")";
    return out_str.str();
}