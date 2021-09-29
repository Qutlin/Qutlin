package Qutlin

import org.hipparchus.ode.ODEIntegrator
import org.hipparchus.ode.ODEState
import org.hipparchus.ode.OrdinaryDifferentialEquation
import org.hipparchus.ode.nonstiff.GraggBulirschStoerIntegrator
import org.hipparchus.util.FastMath.sqrt


typealias State = ComplexArray
typealias Operator = ComplexMatrix

fun State.normalize() {
    val norm = sqrt((this.bra() * this).real)
    for(i in this.indices) set(i, get(i) / norm)
}

fun State.normalized(): State {
    val res = copyOf()
    res.normalize()
    return res
}

fun State.bra(): State {
    val res = copyOf()
    for(i in res.indices) res[i] = res[i].conjugate()
    return res
}

fun Operator.dagger(): Operator {
    return conjugateTranspose();
}

fun State.ketBra() = this / this.bra()


/**
 * Returns the result of a collapse operator `O` acting on a density matrix `ρ`:
 *
 * `collapse(ρ) = O ρ O^† - 1/2 [O^† O, ρ]`
 */
fun collapse(O: Operator, ρ: Operator): Operator {
    val opDagger = O.dagger();
    return O * ρ * opDagger - 0.5 * antiCommutator(opDagger * O, ρ);
}



fun solveEvolution(
    t0: Double = 0.0,
    t: Double,
    initial: ComplexMatrix,
    derivative: (Double, ComplexMatrix) -> ComplexMatrix,
    integrator: ODEIntegrator = GraggBulirschStoerIntegrator(Double.MIN_VALUE, 1.0, 1e-8, 1e-8)
) : Pair<Double, ComplexMatrix> {

    val ode = object: OrdinaryDifferentialEquation {
        override fun computeDerivatives(t: Double, rho: DoubleArray) =
            derivative(t, ComplexMatrix.fromDoubleArray(initial.dimensions, rho)).toDoubleArray()
        override fun getDimension(): Int = initial.dimensions.first * initial.dimensions.second * 2
    }

    val initialState = ODEState(t0, initial.toDoubleArray())

    val res = integrator.integrate(ode, initialState, t )
    return Pair(res.time, ComplexMatrix.fromDoubleArray(initial.dimensions, res.primaryState))
}



fun main() {
    val rho = ComplexMatrix.ones(Pair(2, 2))
    println("rho = $rho")
    val v0 = ComplexArray(2) { 0.0 * R }
    v0[0] = 1.0.toComplex()
    println("v0 = ${v0.str()}")

    val v1 = ComplexArray(2) { 0.0 * R }
    v1[1] = 1.0.toComplex()
    println("v1 = ${v1.str()}")


    val c = v0.bra() / v1
    println("c = $c")
    val coll = collapse(c, rho)
    println("collapse = $coll")
}