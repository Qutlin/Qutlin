package Qutlin

import org.hipparchus.complex.Complex
import java.lang.RuntimeException


typealias ComplexArray = Array<Complex>

fun ComplexArray.str(): String {
    var res = ""
    this.forEach { res += "$it, " }
    return res
}


inline fun ComplexArray.norm() = (this.bra() * this).real


fun complexArrayOf(vararg numbers: Any): ComplexArray {
    val res = ComplexArray(numbers.size) { Complex.ZERO }
    for (i in numbers.indices) {
        when (val x = numbers[i]) {
            is Complex -> res[i] = x
            is Double -> res[i] = Complex(x.toDouble())
            is Int -> res[i] = Complex(x.toDouble())
            else -> throw RuntimeException("No valid number chosen! ${x}")
        }
    }
    return res
}

fun List<Complex>.toComplexArray(): ComplexArray = this.toTypedArray()


operator fun ComplexArray.plus(other: ComplexArray): ComplexArray {
    val res = copyOf()
    for (i in res.indices) res[i] += other[i]
    return res
}

operator fun ComplexArray.times(other: Complex): ComplexArray {
    val res = copyOf()
    for (i in res.indices) res[i] *= other
    return res
}

inline operator fun ComplexArray.times(other: Double) =
    this.times(Complex(other))

inline operator fun ComplexArray.minus(other: ComplexArray) =
    this.plus(other.times(Complex(-1.0)))

inline operator fun ComplexArray.div(other: Complex) =
    this.times(Complex(1.0) / other)


operator fun ComplexArray.times(other: ComplexMatrix): ComplexArray {
    assert(size == other.dimensions.second)

    val dim = other.dimensions.first
    val res = ComplexArray(dim) { Complex.ZERO }

    for (d in 0 until dim)
        for (i in 0 until other.dimensions.second) {
            res[d] += this[i] * other[i, d]
        }
    return res
}

operator fun ComplexMatrix.times(other: ComplexArray): ComplexArray {
    assert(dimensions.first == other.size)

    val res = ComplexArray(dimensions.second) { Complex.ZERO }

    for (d in 0..dimensions.second)
        for (i in other.indices) {
            res[d] += this[d, i] * other[i]
        }

    return res
}

operator fun ComplexArray.times(other: ComplexArray): Complex {
    assert(size == other.size)
    var res = Complex.ZERO
    for (i in 0 until size) res += this[i] * other[i]
    return res
}

operator fun ComplexArray.div(other: ComplexArray): ComplexMatrix {
    val res = ComplexMatrix(Pair(size, other.size))
    for (i in indices)
        for (j in other.indices) res[i, j] = this[i] * other[j]
    return res
}