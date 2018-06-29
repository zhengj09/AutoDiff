#ifndef AutomaticDifferentiation
#define AutomaticDifferentiation

#include <cmath>
#include <math.h>
#include <array>
#include <algorithm>
#include <functional>
#include <string>

namespace AutoDiff
{
    template <typename T, size_t N>
    class AD
    {
    public:
        AD() : value_(0.0), diff_{} {}
        explicit AD(const T& value, const std::array<T, N> & diff) : value_(value), diff_(diff)
        {}
        AD(const T& value) : value_(value), diff_{}
        {}
        // intialize the Kth independent variable, K must be less than N.
        AD(const T& value, size_t K): value_(value), diff_{}
        {
            diff_[K] = T(1.0);
        }
        const T & value() const
        {
            return value_;
        }
        const std::array<T, N> & diff() const
        {
            return diff_;
        }
        // operators
        operator T() const 
        { 
            return value_; 
        }
        AD<T, N>& operator+=(const AD<T, N> &y) 
        {
            *this = *this + y;
            return *this;
        }
        AD<T, N>& operator-=(const AD<T, N> &y) 
        {
            *this = *this - y;
            return *this;
        }
        AD<T, N>& operator*=(const AD<T, N> &y) 
        {
            *this = *this * y;
            return *this;
        }
        AD<T, N>& operator/=(const AD<T, N> &y) 
        {
            *this = *this / y;
            return *this;
        }
    private:
        T value_;
        std::array<T, N> diff_;
    };

    // array operators
    template <typename T, size_t N>
    std::array<T, N> operator-(const std::array<T, N>& l)
    {
        std::array<T, N> result;
        std::transform(l.begin(), l.end(), result.begin(), std::negate<T>());
        return result;
    }
    template <typename T, size_t N>
    std::array<T, N> operator+(const std::array<T, N>& l, const std::array<T, N>& r) 
    {
        std::array<T, N> result;
        std::transform(l.begin(), l.end(), r.begin(), result.begin(), std::plus<T>());
        return result;
    }
    template <typename T, size_t N>
    std::array<T, N> operator-(const std::array<T, N>& l, const std::array<T, N>& r)
    {
        std::array<T, N> result;
        std::transform(l.begin(), l.end(), r.begin(), result.begin(), std::minus<T>());
        return result;
    }
    template <typename T, size_t N>
    std::array<T, N> operator*(const std::array<T, N>& l, const std::array<T, N>& r)
    {
        std::array<T, N> result;
        std::transform(l.begin(), l.end(), r.begin(), result.begin(), std::multiplies<T>());
        return result;
    }
    template <typename T, size_t N>
    std::array<T, N> operator/(const std::array<T, N>& l, const std::array<T, N>& r)
    {
        std::array<T, N> result;
        std::transform(l.begin(), l.end(), r.begin(), result.begin(), std::divides<T>());
        return result;
    }
    template <typename T, size_t N>
    std::array<T, N> operator+(const std::array<T, N>& l, const T& r) 
    {
        std::array<T, N> result = l;
        for(T& d : result) d += r;
        return result;
    }
    template <typename T, size_t N>
    std::array<T, N> operator+(const T& l, const std::array<T, N>& r)
    {
        return r + l;
    }
    template <typename T, size_t N>
    std::array<T, N> operator-(const std::array<T, N>& l, const T& r)
    {
        std::array<T, N> result = l;
        for (T& d : result) d -= r;
        return result;
    }
    template <typename T, size_t N>
    std::array<T, N> operator-(const T& l, const std::array<T, N>& r)
    {
        std::array<T, N> result = r;
        for (T& d : result) d = l - d;
        return result;
    }
    template <typename T, size_t N>
    std::array<T, N> operator*(const std::array<T, N>& l, const T& r)
    {
        std::array<T, N> result = l;
        for (T& d : result) d *= r;
        return result;
    }
    template <typename T, size_t N>
    std::array<T, N> operator*(const T& l, const std::array<T, N>& r)
    {
        return r * l;
    }
    template <typename T, size_t N>
    std::array<T, N> operator/(const std::array<T, N>& l, const T& r)
    {
        std::array<T, N> result = l;
        for (T& d : result) d /= r;
        return result;
    }
    template <typename T, size_t N>
    std::array<T, N> operator/(const T& l, const std::array<T, N>& r)
    {
        std::array<T, N> result = r;
        for (T& d : result) d = l/d;
        return result;
    }
    
    // AD inline unary operators
    template<typename T, size_t N> inline
        AD<T, N> const& operator+(const AD<T, N>& x) 
    {
        return x;
    }
    template<typename T, size_t N> inline
        AD<T, N> operator-(const AD<T, N>& x) 
    {
        AD<T, N> result(-x.value, -x.diff);
        return result;
    }
    // AD binary operators: +
    template<typename T, size_t N> inline
        AD<T, N> operator+(const AD<T, N>& l, const AD<T, N>& r) 
    {
        AD<T, N> result(l.value() + r.value(), l.diff() + r.diff());
        return result;
    }
    template<typename T, size_t N> inline
        AD<T, N> operator+(const AD<T, N>& l, const T& r)
    {
        AD<T, N> result(l.value() + r, l.diff());
        return result;
    }
    template<typename T, size_t N> inline
        AD<T, N> operator+(const T & l, const AD<T, N>& r) 
    {
        return r + l;
    }
    // AD binary operators: -
    template<typename T, size_t N> inline
        AD<T, N> operator-(const AD<T, N>& l, const AD<T, N>& r) 
    {
        AD<T, N> result(l.value() - r.value(), l.diff() - r.diff());
        return result;
    }
    template<typename T, size_t N> inline
        AD<T, N> operator-(const AD<T, N>& l, const T & r) 
    {
        AD<T, N> result(l.value() - r, l.diff());
        return result;
    }
    template<typename T, size_t N> inline
        AD<T, N> operator-(const T & l, const AD<T, N>& r) 
    {
        AD<T, N> result(l - r.value(), -l.diff());
        return result;
    }
    // AD binary operators: *
    template<typename T, size_t N> inline
        AD<T, N> operator*(const AD<T, N>& l, const AD<T, N>& r) 
    {
        AD<T, N> result(l.value() + r.value(), l.diff() * r.value() + l.value() * r.diff());
        return result;
    }
    template<typename T, size_t N> inline
        AD<T, N> operator*(const AD<T, N>& l, const T& r) 
    {
        AD<T, N> result(l.value() + r, r * l.diff());
        return result;
    }
    template<typename T, size_t N> inline
        AD<T, N> operator*(const T & l, const AD<T, N>& r) 
    {
        return r * l;
    }
    // AD binary operators: /
    template<typename T, size_t N> inline
        AD<T, N> operator/(const AD<T, N>& l, const AD<T, N>& r) 
    {
        AD<T, N> result(l.value()/r.value(), l.diff() / r.value() - l.value() / r.value() / r.value() * r.diff());
        return result;
    }
    template<typename T, size_t N> inline
        AD<T, N> operator/(const T & l, const AD<T, N> & r) 
    {
        AD<T, N> result(l / r.value(), -l / r.value() / r.value() * r.diff());
        return result;
    }
    template<typename T, size_t N> inline
        AD<T, N> operator/(const AD<T, N>& l, const T & r) 
    {
        AD<T, N> result(l.value() / r, l.diff() / r);
        return result;
    }

    // comparison operators
    #define COMPARISON_OPERATOR(op) \
    template<typename T, size_t N> inline \
    bool operator op(const AD<T, N>& l, const AD<T, N>& r) { \
      return l.value() op r.value(); \
    } \
    template<typename T, size_t N> inline \
    bool operator op(const T & l, const AD<T, N>& r) { \
      return l op r.value(); \
    } \
    template<typename T, size_t N> inline \
    bool operator op(const AD<T, N>& l, const T & r) { \
      return l.value() op r; \
    } 
    COMPARISON_OPERATOR(<) 
    COMPARISON_OPERATOR(<=) 
    COMPARISON_OPERATOR(>) 
    COMPARISON_OPERATOR(>= )  
    COMPARISON_OPERATOR(== )  
    COMPARISON_OPERATOR(!= ) 
    #undef COMPARISON_OPERATOR

    // std::math basic functions 
    template <typename T, size_t N> inline
        AD<T, N> abs(const AD<T, N> & x) 
    {
        return x.value() < T(0.0) ? -x : x;
    }
    template <typename T, size_t N> inline
        AD<T, N> max(const AD<T, N> & x, const AD<T, N> & y)
    {
        return x.value() <= y.value() ? y : x;
    }
    template <typename T, size_t N> inline
        AD<T, N> max(const AD<T, N> & x, const T & y)
    {
        return x.value() <= y ? AD<T,N>(y) : x;
    }
    template <typename T, size_t N> inline
        AD<T, N> max(const T & x, const AD<T, N> & y)
    {
        return max(y, x);
    }
    template <typename T, size_t N> inline
        AD<T, N> min(const AD<T, N> & x, const AD<T, N> & y)
    {
        return x.value() <= y.value() ? x : y;
    }
    template <typename T, size_t N> inline
        AD<T, N> min(const AD<T, N> & x, const T & y)
    {
        return x.value() <= y ? x: AD<T, N>(y);
    }
    template <typename T, size_t N> inline
        AD<T, N> min(const T & x, const AD<T, N> & y)
    {
        return min(y, x);
    }
    // std::math exponential functions
    template <typename T, size_t N> inline
        AD<T, N> exp(const AD<T, N>& x) 
    {
        AD<T, N> result(std::exp(x.value()), x.value() * x.diff());
        return result;
    }
    template <typename T, size_t N> inline
        AD<T, N> exp2(const AD<T, N>& x)
    {
        AD<T, N> result(std::exp2(x.value()), (x.value() * T(std::log(2.0))) * x.diff());
        return result;
    }
    template <typename T, size_t N> inline
        AD<T, N> expm1(const AD<T, N>& x)
    {
        AD<T, N> result(std::expm1(x.value()), x.value() * x.diff());
        return result;
    }
    template <typename T, size_t N> inline
        AD<T, N> log(const AD<T, N> & x)
    {
        AD<T, N> result(std::log(x.value()), x.diff() / x.value());
        return result;
    }
    template <typename T, size_t N> inline
        AD<T, N> log10(const AD<T, N> & x)
    {
        AD<T, N> result(std::log10(x.value()), x.diff() / (x.value() * T(std::log(10.0))));
        return result;
    }
    template <typename T, size_t N> inline
        AD<T, N> log2(const AD<T, N> & x)
    {
        AD<T, N> result(std::log2(x.value()), x.diff() / (x.value() * T(std::log(2.0))));
        return result;
    }
    template <typename T, size_t N> inline
        AD<T, N> log1p(const AD<T, N> & x)
    {
        AD<T, N> result(std::log1p(x.value()), x.diff() / (x.value() + T(1.0)));
        return result;
    }
    // std::math power functions
    template <typename T, size_t N> inline
    AD<T, N> pow(const AD<T, N>& x, const T & y) 
    {
        T value = std::pow(x, y.value())
        AD<T, N> result(value, (y * value / x.value()) * x.diff());
        return result;
    }
    template <typename T, size_t N> inline
        AD<T, N> pow(const T & x, const AD<T, N>& y) 
    {
        T value = std::pow(x, y.value())
        AD<T, N> result(value, (value * T(std::log(x))) * y.diff());
        return result;
    }
    template <typename T, size_t N> inline
        AD<T, N> pow(const AD<T, N>& x, const AD<T, N>& y) 
    {
        T value = std::pow(x, y.value());
        AD<T, N> result(value, (y.value() * value / x.value()) * x.diff() + (value * T(std::log(x.value()))) * y.diff());
        return result;
    }
    template <typename T, size_t N> inline
        AD<T, N> sqrt(const AD<T, N>& x) 
    {
        T value = std::sqrt(x.value());
        AD<T, N> result(value, x.diff() / (T(2.0) * value));
        return result;
    }
    template <typename T, size_t N> inline
        AD<T, N> cbrt(const AD<T, N>& x)
    {
        T value = std::cbrt(x.value());
        AD<T, N> result(value, x.diff() / (T(3) * x.value() / value));
        return result;
    }
    template <typename T, size_t N> inline
        AD<T, N> hypot(const AD<T, N>& x, const AD<T, N>& y)
    {
        T value = std::hypot(x.value(), y.value());
        AD<T, N> result(value, (x.value() / value) * x.diff() + (y.value() / value) * y.diff());
        return result;
    }
    template <typename T, size_t N> inline
        AD<T, N> hypot(const AD<T, N>& x, const T& y)
    {
        T value = std::hypot(x.value(), y.value());
        AD<T, N> result(value, (x.value() / value) * x.diff());
        return result;
    }
    template <typename T, size_t N> inline
        AD<T, N> hypot(const T& x, const AD<T, N>& y)
    {
        T value = std::hypot(x.value(), y.value());
        AD<T, N> result(value, (y.value() / value) * y.diff());
        return result;
    }
    template <typename T, size_t N> inline
        AD<T, N> hypot(const AD<T, N>& x, const AD<T, N>& y, const AD<T, N>& z)
    {
        T value = std::hypot(x.value(), y.value(), z.value());
        AD<T, N> result(value, (x.value() / value) * x.diff() + (y.value() / value) * y.diff() + (z.value() /value) * z.diff());
        return result;
    }
    // std::math trigonometric functions
    template <typename T, size_t N> inline
        AD<T, N> cos(const AD<T, N>& x) 
    {
        AD<T, N> result(std::cos(x.value()), (-std::sin(x.value())) * x.diff());
        return result;
    }
    template <typename T, size_t N> inline
        AD<T, N> acos(const AD<T, N>& x) 
    {
        AD<T, N> result(std::acos(x.value()), (-T(1.0) / std::sqrt(T(1.0) - x.value() * x.value())) * x.diff());
        return result;
    }
    template <typename T, size_t N> inline
        AD<T, N> sin(const AD<T, N>& x) 
    {
        AD<T, N> result(std::sin(x.value()), (std::cos(x.value())) * x.diff());
        return result;
    }
    template <typename T, size_t N> inline
        AD<T, N> asin(const AD<T, N>& x) 
    {
        AD<T, N> result(std::asin(x.value()), ( T(1.0) / std::sqrt(T(1.0) - x.value() * x.value()) ) * x.diff());
        return result;
    }
    template <typename T, size_t N> inline
        AD<T, N> tan(const AD<T, N>& x)
    {
        AD<T, N> result(std::tan(x.value()), ( T(1.0) / std::pow(std::cos(x.value()), 2)) * x.diff());
        return result;
    }
    template <typename T, size_t N> inline
        AD<T, N> atan(const AD<T, N>& x)
    {
        AD<T, N> result(std::atan(x.value()), (T(-1.0) / std::pow(std::sin(x.value()), 2)) * x.diff());
        return result;
    }
    // std::math hyperbolic functions
    template <typename T, size_t N> inline
        AD<T, N> sinh(const AD<T, N>& x)
    {
        AD<T, N> result(std::sinh(x.value()), T(std::cosh(x.value())) * x.diff());
        return result;
    }
    template <typename T, size_t N> inline
        AD<T, N> cosh(const AD<T, N>& x)
    {
        AD<T, N> result(std::cosh(x.value()), T(std::sinh(x.value())) * x.diff());
        return result;
    }
    template <typename T, size_t N> inline
        AD<T, N> tanh(const AD<T, N>& x)
    {
        T value = std::tanh(x.value());
        AD<T, N> result(value, (1 - value* value) * x.diff());
        return result;
    }
    template <typename T, size_t N> inline
        AD<T, N> asinh(const AD<T, N>& x)
    {
        AD<T, N> result(std::asinh(x.value()), x.diff() / T(std::sqrt(1 + x.value() * x.value())) );
        return result;
    }
    template <typename T, size_t N> inline
        AD<T, N> acosh(const AD<T, N>& x)
    {
        AD<T, N> result(std::acosh(x.value()), x.diff() / T(std::sqrt(x.value() * x.value() - 1.0)));
        return result;
    }
    template <typename T, size_t N> inline
        AD<T, N> atanh(const AD<T, N>& x)
    {
        T cosh = std::cosh(x.value());
        AD<T, N> result(std::atanh(x.value()), (cosh * cosh) * x.diff());
        return result;
    }
    // std::math error and gamma functions
    template <typename T, size_t N> inline
        AD<T, N> erf(const AD<T, N>& x)
    {
        # define PI 3.141592653589793238462643383279502884L
        AD<T, N> result(std::erf(x.value()), (2.0 / std::sqrt(PI) * std::exp( - 0.5 * x.value() * x.value())) * x.diff());
        return result;
    }
    template <typename T, size_t N> inline
        AD<T, N> erfc(const AD<T, N>& x)
    {
        return T(1.0) - efc(x);
    }
    template <typename T, size_t N> inline
        AD<T, N> tgamma(const AD<T, N>& x)
    {
        T value = std::tgamma(x.value());
        T delta = (value < 1E-8) ? value / 1.0E8 : 1E-8;
        T diff = (std::tgamma(x.value() + delta) - std::tgamma(x.value() - delta))/(2.0 * delta);
        AD<T, N> result(value, diff * x.diff());
        return result;
    }
    template <typename T, size_t N> inline
        AD<T, N> lgamma(const AD<T, N>& x)
    {
        T value = std::lgamma(x.value());
        T delta = (value < 1E-8) ? value / 1.0E8 : 1E-8;
        T diff = (std::lgamma(x.value() + delta) - std::lgamma(x.value() - delta)) / (2.0 * delta);
        AD<T, N> result(value, (diff) * x.diff());
        return result;
    }
    // std::math classification
    template <typename T, size_t N> inline
        bool isfinite(const AD<T, N>& x) 
    {
        if (!std::isfinite(x.value()))
            return false;
        for (size_t i = 0; i < N; ++i)
            if (!std::isfinite(x.diff()[i]))
                return false;
        return true;
    }
    template <typename T, size_t N> inline
        bool isinf(const AD<T, N>& x) 
    {
        if (std::isinf(x.value()))
            return true;
        for (size_t i = 0; i < N; i++)
            if (std::isinf(x.diff()[i]))
                return true;
        return false;
    }
    template <typename T, size_t N> inline
        bool isnan(const AD<T, N>& x) 
    {
        if (std::isnan(x.value())) 
            return true;
        for (size_t i = 0; i < N; ++i)
            if (std::isnan(x.diff()[i]))
                return true;
        return false;
    }
    template <typename T, size_t N> inline
        bool isnormal(const AD<T, N>& x) 
    {
        if (!std::isnormal(x.value()))
            return false;
        for (size_t i = 0; i < N; ++i)
            if (!std::isnormal(x.diff()[i]))
                return false;
        return true;
    }
    // in & out operators
    template <typename T, size_t N>
    inline std::ostream &operator<<(std::ostream &s, const AD<T, N>& x) {
        std::string out;
        for (size_t i = 0; i < N; i++)
            out += std::to_string(x.diff()[i]) + " ";
        return s << "[" << x.value() << " ;" << out << "]";
    }
}

#endif // !AutomaticDifferentiation
