#ifndef DOUBLEBEAM_CPP_STRONG_TYPES_HPP
#define DOUBLEBEAM_CPP_STRONG_TYPES_HPP

#include <typeinfo>
#include <utility>
#include <ostream>


/**
 * Base class for strong named types.
 * Use to declare new strong types holding a value of type T:
 * using Meter = NamedType<double, MeterParameter>;
 * @tparam T Base type of the thing which NamedType is holding.
 * @tparam Parameter Type that makes instantiations of NamedType unique even when T is the same.
 */
template <typename T, typename Parameter>
class NamedType {
public:
    explicit NamedType(const T& value) : value(value) {}
    explicit NamedType(T&& value) : value(std::move(value)) {}

    T& get() {
        return value;
    }
    const T& get() const {
        return value;
    }

    friend std::ostream& operator<<(std::ostream& os, const NamedType& quantity) {
        return os << quantity.value;
    }

    bool operator==(const NamedType& quantity) const{
        return value == quantity.value;
    }

private:
    T value;
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