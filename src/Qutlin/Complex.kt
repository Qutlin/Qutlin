package Qutlin

import org.hipparchus.complex.Complex
import org.hipparchus.special.Erf

val I: Complex = Complex.I
val R: Complex = Complex(1.0)

inline fun Double.toComplex() = Complex(this)

inline operator fun Complex.unaryMinus(): Complex = this.multiply(-1.0)

inline operator fun Complex.plus(c: Complex): Complex = this.add(c)
inline operator fun Complex.plus(c: Double): Complex = this.add(c)
inline operator fun Double.plus(c: Complex): Complex = c.add(this)

inline operator fun Complex.minus(c: Complex): Complex = this.subtract(c)
inline operator fun Complex.minus(c: Double): Complex = this.subtract(c)
inline operator fun Double.minus(c: Complex): Complex = c.subtract(this)

inline operator fun Complex.times(c: Complex): Complex = this.multiply(c)
inline operator fun Complex.times(c: Double): Complex = this.multiply(c)
inline operator fun Double.times(c: Complex): Complex = c.multiply(this)

inline operator fun Complex.div(c: Complex): Complex = this.divide(c)
inline operator fun Double.div(c: Complex): Complex = c.divide(this)
inline operator fun Complex.div(c: Double): Complex = this.divide(c)

inline fun exp(c: Complex) = c.exp()!!


inline fun Complex.absSquare() = this.real * this.real + this.imaginary * this.imaginary
