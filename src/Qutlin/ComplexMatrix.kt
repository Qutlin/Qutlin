package Qutlin

import org.hipparchus.complex.Complex
import org.hipparchus.linear.Array2DRowRealMatrix
import org.hipparchus.linear.EigenDecomposition
import org.hipparchus.util.FastMath.*


data class ComplexMatrix(
    val dimensions: Pair<Int, Int>,
    val values: ComplexArray = ComplexArray(dimensions.first * dimensions.second) { Complex(0.0, 0.0) }
) {
    companion object {
        fun eye(dimensions: Pair<Int, Int>): ComplexMatrix {
            val res = ComplexMatrix(dimensions)
            for (i in 0 until min(dimensions.first, dimensions.second)) {
                res[i, i] = 1.0.toComplex()
            }
            return res;
        }

        fun ones(dimensions: Pair<Int, Int>): ComplexMatrix = ComplexMatrix(
            dimensions,
            Array(dimensions.first * dimensions.second) { Complex(1.0, 0.0) }
        )

        fun fromDoubleArray(dimensions: Pair<Int, Int>, data: DoubleArray): ComplexMatrix {
            assert(data.size == dimensions.first * dimensions.second * 2)
            return ComplexMatrix(
                dimensions,
                ComplexArray(dimensions.first * dimensions.second) { Complex(data[2 * it], data[2 * it + 1]) }
            )
        }

        fun pauliX() = ComplexMatrix(Pair(2,2), complexArrayOf(
            0.0, R,
            R, 0.0
        ))
        fun pauliY() = ComplexMatrix(Pair(2,2), complexArrayOf(
            0.0, -I,
            I, 0.0
        ))
        fun pauliZ() = ComplexMatrix(Pair(2,2), complexArrayOf(
            R, 0.0,
            0.0, -R
        ))

        /**
         * calculates `Rx(θ) = e^(-i θ/2 pauliX)`
         */
        fun rotPauliX(θ: Double) = ComplexMatrix(Pair(2,2), complexArrayOf(
            cos(θ/2.0), -I*sin(θ/2.0),
            -I*sin(θ/2.0), cos(θ/2.0)
        ))
        /**
         * calculates `Ry(θ) = e^(-i θ/2 pauliY)`
         */
        fun rotPauliY(θ: Double) = ComplexMatrix(Pair(2,2), complexArrayOf(
            cos(θ/2.0), -sin(θ/2.0),
            sin(θ/2.0), cos(θ/2.0)
        ))
        /**
         * calculates `Rz(θ) = e^(-i θ/2 pauliZ)`
         */
        fun rotPauliZ(θ: Double) = ComplexMatrix(Pair(2,2), complexArrayOf(
            exp(-I*θ/2.0), 0.0,
            0.0, exp(I*θ/2.0)
        ))
    }


    override fun equals(other: Any?): Boolean {
        if (this === other) return true
        if (javaClass != other?.javaClass) return false

        other as ComplexMatrix

        if (dimensions != other.dimensions) return false
        if (!values.contentEquals(other.values)) return false

        return true
    }

    override fun hashCode(): Int {
        var result = dimensions.hashCode()
        result = 31 * result + values.contentHashCode()
        return result
    }

    fun toDoubleArray() = DoubleArray(values.size * 2) {
        if (it % 2 == 0) values[it / 2].real else values[(it - 1) / 2].imaginary
    }

    fun toRealMatrix(): Array2DRowRealMatrix {
        val res = Array2DRowRealMatrix(dimensions.first * 2, dimensions.second * 2)
        for (i in 0 until dimensions.first)
            for (j in 0 until dimensions.second) {
                val c = this[i, j]
                res.setEntry(2 * i, 2 * j, c.real)
                res.setEntry(2 * i, 2 * j + 1, -c.imaginary)
                res.setEntry(2 * i + 1, 2 * j, c.imaginary)
                res.setEntry(2 * i + 1, 2 * j + 1, c.real)
            }
        return res
    }

    /**
     * Calculates the eigenvalues and eigenvectors of this ComplexMatrix.
     * @return A List of Pairs of eigenvalue and corresponding eigenvector. The eigenvalues are sorted from lowest to highest.
     */
    fun eigenSystem(): List<Pair<Complex, ComplexArray>> {
//        val solver = ComplexEigenDecomposition(toRealMatrix())
//        val eigenvalues: Array<Complex> = solver.eigenvalues;

//        println("Calculating EigenDecomposition of ")
//        println(this)
        val solver = EigenDecomposition(toRealMatrix())

        val eigenvalues: Array<Complex> = solver.realEigenvalues.zip(solver.imagEigenvalues) {
                r,i -> Complex(r,i)
        }.toComplexArray()

        val eigenvectors = List(eigenvalues.size) { ComplexArray(dimensions.first) { Complex.ZERO } }

        for (i in eigenvalues.indices) {
            val e0 = solver.getEigenvector(i).toArray()
            for (j in eigenvectors[i].indices)
                eigenvectors[i][j] = e0[2 * j] + I * e0[2 * j + 1]
        }

        // ? sort results by real parts of eigenvalues
        return eigenvalues.mapIndexed { i, v ->
            Pair(v, eigenvectors[i])
        }.sortedBy { p -> p.first.real }.chunked(2) { l -> l[0] }
    }

    fun trace(): Complex {
        var res = Complex.ZERO
        for (i in 0 until min(dimensions.first, dimensions.second))
            res += this[i, i]
        return res
    }


    fun conjugateTransposed() {
        if (dimensions.first != dimensions.second) throw Error("matrix not square!");
        for (i in 0 until dimensions.first) {
            for (j in i until dimensions.first) {
                if (i == j) {
                    this[i, i] = this[i, i].conjugate();
                    continue;
                }

                val tmp = this[i, j];
                this[i, j] = this[j, i].conjugate();
                this[j, i] = tmp.conjugate();
            }
        }
    }

    fun conjugateTranspose(): ComplexMatrix {
        val res = ComplexMatrix(dimensions = dimensions, values = values.clone());
        res.conjugateTransposed();
        return res;
    }

}

inline operator fun ComplexMatrix.get(i: Int, j: Int): Complex = values[dimensions.second * i + j]
inline operator fun ComplexMatrix.set(i: Int, j: Int, value: Complex) {
    values[dimensions.second * i + j] = value
}

operator fun ComplexMatrix.plus(other: ComplexMatrix): ComplexMatrix {
    val res = this.copy()
    for (i in 0 until dimensions.first)
        for (j in 0 until dimensions.second)
            res[i, j] = res[i, j] + other[i, j]
    return res
}

operator fun ComplexMatrix.minus(other: ComplexMatrix): ComplexMatrix {
    val res = this.copy()
    for (i in 0 until dimensions.first)
        for (j in 0 until dimensions.second)
            res[i, j] = res[i, j] - other[i, j]
    return res
}

operator fun ComplexMatrix.times(other: ComplexMatrix): ComplexMatrix {
    val res = ComplexMatrix(dimensions = Pair(dimensions.first, other.dimensions.second))
    for (j in 0 until other.dimensions.second)
        for (k in 0 until dimensions.second) {
            val tmp = other[k, j];
            for (i in 0 until dimensions.first) {
                res[i, j] = res[i, j] + this[i, k] * tmp
            }
        }
    return res
}


operator fun ComplexMatrix.times(other: Complex): ComplexMatrix {
    val res = this.copy()
    for (i in values.indices) res.values[i] *= other
    return res
}

inline operator fun ComplexMatrix.times(other: Double) = this * other.toComplex()
inline operator fun Double.times(other: ComplexMatrix) = other * this
inline operator fun Complex.times(other: ComplexMatrix) = other * this



operator fun ComplexMatrix.div(other: Complex): ComplexMatrix {
    val res = this.copy()
    for (i in values.indices) res.values[i] /= other
    return res
}
inline operator fun ComplexMatrix.div(other: Double) = this / other.toComplex()
inline operator fun Double.div(other: ComplexMatrix) = other / this
inline operator fun Complex.div(other: ComplexMatrix) = other / this


inline fun commutator(a: ComplexMatrix, b: ComplexMatrix) = a * b - b * a
inline fun antiCommutator(a: ComplexMatrix, b: ComplexMatrix) = a * b + b * a


fun main() {
    var m = ComplexMatrix.eye(Pair(2, 2))
    m[0, 1] = 3.0 + 4.0 * I
    println(m)
    m = 0.5.toComplex() * (m + m.conjugateTranspose())
    println(m)

    val eig = m.eigenSystem()
    eig.forEach {
        println("eigenvalue: ${it.first}")
        it.second.forEach { v -> print("$v, ") }
        println()
    }
}