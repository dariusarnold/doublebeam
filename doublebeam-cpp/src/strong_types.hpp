#ifndef DOUBLEBEAM_CPP_STRONG_TYPES_HPP
#define DOUBLEBEAM_CPP_STRONG_TYPES_HPP

#include <ostream>
#include <typeinfo>
#include <utility>

#include <boost/operators.hpp>


/**
 * Base class for strong named types.
 * Use to declare new strong types holding a value of type T:
 * using Meter = NamedType<double, MeterParameter>;
 * @tparam T Base type of the thing which NamedType is holding.
 * @tparam Parameter Type that makes instantiations of NamedType unique even when T is the same.
 */
template <typename T, typename Parameter>
class NamedType : boost::totally_ordered<NamedType<T, Parameter>>,
                  boost::additive<NamedType<T, Parameter>>,
                  boost::multiplicative<NamedType<T, Parameter>, T> {
public:
    NamedType() : value_m(0.) {}
    explicit NamedType(const T& value) : value_m(value) {}
    explicit NamedType(T&& value) : value_m(std::move(value)) {}

    T& get() {
        return value_m;
    }
    const T& get() const {
        return value_m;
    }

    friend std::ostream& operator<<(std::ostream& os, const NamedType& quantity) {
        return os << quantity.get();
    }

    friend std::istream& operator>>(std::istream& is, NamedType& quantity) {
        is >> quantity.get();
        return is;
    }

    // comparison of strong types
    bool operator==(const NamedType& other) const {
        return value_m == other.get();
    }

    bool operator<(const NamedType& other) const {
        return value_m < other.get();
    }

    // addition/subtraction of strong types
    NamedType& operator+=(const NamedType& rhs) {
        value_m += rhs.get();
        return *this;
    }

    NamedType& operator-=(const NamedType& rhs) {
        value_m -= rhs.get();
        return *this;
    }

    // division/multiplication by another scalar
    NamedType& operator*=(T factor) {
        value_m *= factor;
        return *this;
    }

    NamedType& operator/=(T factor) {
        value_m /= factor;
        return *this;
    }

    // unary negative/positive operator
    NamedType operator-() const {
        return NamedType(-value_m);
    }

    NamedType& operator+() const {
        return *this;
    }

    friend NamedType abs(const NamedType& quantity) {
        return NamedType{std::abs(quantity.value)};
    }

    friend bool isnan(const NamedType& quantity) {
        return std::isnan(quantity.get());
    }

private:
    T value_m;
};

/**
 * Helper macro to define a new strong type.
 * @param NameOfType Name of the new strong type.
 * @param BaseType Underlying type of the thing which the strong type is holding.
 */
#define DEFINE_STRONG_TYPE(NameOfType, BaseType)                                                   \
    using NameOfType = NamedType<BaseType, struct PhantomType##NameOfType>;

/**
 * Helper macro to create user defined literals for a custom type.
 * @param NameOfType Name of an existing strong type for which to define literals.
 * @param Suffix Suffix to use for user defined literal.
 */
#define DEFINE_TYPE_LITERAL(NameOfType, Suffix)                                                    \
    inline NameOfType operator""##Suffix(long double param) {                                      \
        return NameOfType(param);                                                                  \
    }                                                                                              \
    inline NameOfType operator""##Suffix(unsigned long long param) {                               \
        return NameOfType(param);                                                                  \
    }

/**
 * Helper macro to create user defined literals for a custom type that multiply input parameter by a
 * factor.
 * @param NameOfType Name of an existing strong type for which to define literals.
 * @param Suffix Suffix to use for user defined literal.
 * @param Factor The value of the literal will be multiplied by this constant.
 */
#define DEFINE_TYPE_LITERAL_WITH_FACTOR(NameOfType, Suffix, Factor)                                \
    inline NameOfType operator""##Suffix(long double param) {                                      \
        return NameOfType(param * Factor);                                                         \
    }                                                                                              \
    inline NameOfType operator""##Suffix(unsigned long long param) {                               \
        return NameOfType(param * Factor);                                                         \
    }


/** Examples of the macro usages:
 * Define a strong type called MyType, underlying type double
DEFINE_STRONG_TYPE(MyType, double);

 * Define type literal so that 1_mytype will be converted to MyType(1)
 * and 2.3_mytype will be converted to MyType(2.3)
DEFINE_TYPE_LITERAL(MyType, _mytype);

 * Define type literal with factor of 1000, resulting in 1.4_kmytpye == MyType(1400.)
DEFINE_TYPE_LITERAL_WITH_FACTOR(MyType, _kmytype, 1000);

*/

#endif // DOUBLEBEAM_CPP_STRONG_TYPES_HPP
