
import Qutlin.*
import kotlinx.coroutines.async
import kotlinx.coroutines.awaitAll
import kotlinx.coroutines.coroutineScope
import org.hipparchus.analysis.integration.SimpsonIntegrator
import org.hipparchus.analysis.integration.UnivariateIntegrator
import org.hipparchus.complex.Complex
import org.hipparchus.migration.ode.FirstOrderDifferentialEquations
import org.hipparchus.ode.nonstiff.GraggBulirschStoerIntegrator
import org.hipparchus.util.FastMath
import org.hipparchus.util.FastMath.sqrt


/*
8888888888                                    d8b                                 888
888                                           Y8P                                 888
888                                                                               888
8888888    888  888 88888b.   .d88b.  888d888 888 88888b.d88b.   .d88b.  88888b.  888888
888        `Y8bd8P' 888 "88b d8P  Y8b 888P"   888 888 "888 "88b d8P  Y8b 888 "88b 888
888          X88K   888  888 88888888 888     888 888  888  888 88888888 888  888 888
888        .d8""8b. 888 d88P Y8b.     888     888 888  888  888 Y8b.     888  888 Y88b.
8888888888 888  888 88888P"   "Y8888  888     888 888  888  888  "Y8888  888  888  "Y888
                    888
                    888
                    888
*/



abstract class Experiment (
    open val f_epsilon : (Double, Double) -> Double,
    open val f_d_epsilon : (Double, Double) -> Double,
    val name: String = "Experiment"
) {
    abstract fun dtheta(t: Double, tau: Double): Double
    abstract fun gamma(ep: Double): Double
    abstract fun bx(ep: Double): Double
    abstract fun bz(ep: Double): Double

    fun phi(t: Double, tau: Double): Double {
        if(t <= 0.0) return 0.0
        val integrator: UnivariateIntegrator = SimpsonIntegrator(
            1e-8, 1e-8,2,64
        )

        val b = fun(x: Double): Double {
            val ep = f_epsilon(x, tau)
            val bx = bx(ep)
            val bz = bz(ep)
            return sqrt(bx*bx + bz*bz)
        }

        return integrator.integrate(10000, b, 0.0, t)
    }

    private fun computeDerivatives(t: Double, S: DoubleArray, dS: DoubleArray, tau: Double) {
        val ep = f_epsilon(t, tau)
        val dtheta = dtheta(t,tau)
        val phi = phi(t, tau)
//        val B = sqrt(bx(ep).e2() + bz(ep).e2())
        val i = Complex.I

        val Sp = S[1] + i*S[2]

        dS[0] = (-dtheta * Sp * exp(i * phi)).real
        val c = exp(-i * phi) * S[0] * dtheta - Sp * gamma(ep)
        dS[1] = c.real
        dS[2] = c.imaginary
    }


    fun initRelB(alongB: DoubleArray) : DoubleArray {
        val ep0 = f_epsilon(0.0, 1.0)
        val bx = bx(ep0)
        val bz = bz(ep0)
        val b = sqrt(bx.e2() + bz.e2())
        println("initial bz = $bz, bx = $bx, b = $b")
        return doubleArrayOf(
            bz/b * alongB[0], bx/b * alongB[1], 0.0
        )
    }


    fun integrate(initial: DoubleArray, tau: Double) : DoubleArray {

        val integrator = GraggBulirschStoerIntegrator(
            1e-16,
            1000.0,
            1e-8,
            1e-8
        )

        val eqs = object : FirstOrderDifferentialEquations {
            override fun computeDerivatives(t: Double, y: DoubleArray, yDot: DoubleArray) =
                computeDerivatives(t, y, yDot, tau)
            override fun getDimension(): Int = 3
        }

        val res = DoubleArray(3)
        integrator.integrate(eqs, 0.0, initial, tau, res)

        return res
    }
}


/*
8888888b.  8888888b.         .d8888b.   .d8888b.
888  "Y88b 888  "Y88b       d88P  Y88b d88P  Y88b
888    888 888    888       Y88b.      Y88b.
888    888 888    888        "Y888b.    "Y888b.
888    888 888    888           "Y88b.     "Y88b.
888    888 888    888 888888      "888       "888
888  .d88P 888  .d88P       Y88b  d88P Y88b  d88P
8888888P"  8888888P"         "Y8888P"   "Y8888P"
*/



class DonorDotSS(
    var Om: Double, var A: Double, var g: Double,
    override val f_epsilon: (Double, Double) -> Double,
    override val f_d_epsilon: (Double, Double) -> Double
    ) : Experiment(f_epsilon, f_d_epsilon, "DonorDotSS") {

    override fun dtheta(t: Double, tau: Double): Double {
        val ep = f_epsilon(t, tau)
        val dep = f_d_epsilon(t, tau)
        return Om * dep / (ep.e2() + Om.e2())
    }
    override fun gamma(ep: Double): Double =
        g * ep.e2()  / ( ep.e2() + Om.e2())
    override fun bx(ep: Double): Double = Om
    override fun bz(ep: Double): Double = ep
}


/*
8888888b.  8888888b.         .d8888b. 88888888888
888  "Y88b 888  "Y88b       d88P  Y88b    888
888    888 888    888       Y88b.         888
888    888 888    888        "Y888b.      888
888    888 888    888           "Y88b.    888
888    888 888    888 888888      "888    888
888  .d88P 888  .d88P       Y88b  d88P    888
8888888P"  8888888P"         "Y8888P"     888

*/



class DonorDotST(
    var Om: Double, var A: Double, var g: Double,
    override val f_epsilon: (Double, Double) -> Double,
    override val f_d_epsilon: (Double, Double) -> Double
) : Experiment(f_epsilon, f_d_epsilon, "DonorDotST") {

    override fun dtheta(t: Double, tau: Double): Double {
        val ep = f_epsilon(t, tau)
        val dep = f_d_epsilon(t, tau)
        val beta = sqrt(Om * Om + ep * ep)

        return -A * Om * Om * (ep - 3 * beta) * dep / (
                (Om * Om + ep * ep) * sqrt(2 + 2 * ep/beta) *
                        ((A * A - 4 * Om * Om) * ep - 4 * ep * ep * ep
                                + (A * A + 2 * Om * Om) * beta
                                + 4 * ep * ep * beta)
                )
    }
    override fun gamma(ep: Double): Double =
        g / 4.0 * FastMath.pow(
            (sqrt(ep * ep + Om * Om) - ep) / sqrt(ep * ep + Om * Om)
        , 2)
    override fun bx(ep: Double): Double {
        val beta = sqrt(Om * Om + ep * ep)
        return A / 2.0 * sqrt((ep + beta)/(2.0 * beta))
    }
    override fun bz(ep: Double): Double =
        0.5 * (sqrt(ep * ep + Om * Om) - ep)
}


/*
 .d8888b.                                  888 d8b
d88P  Y88b                                 888 Y8P
Y88b.                                      888
 "Y888b.    8888b.  88888b.d88b.  88888b.  888 888 88888b.   .d88b.
    "Y88b.     "88b 888 "888 "88b 888 "88b 888 888 888 "88b d88P"88b
      "888 .d888888 888  888  888 888  888 888 888 888  888 888  888
Y88b  d88P 888  888 888  888  888 888 d88P 888 888 888  888 Y88b 888
 "Y8888P"  "Y888888 888  888  888 88888P"  888 888 888  888  "Y88888
                                  888                            888
                                  888                       Y8b d88P
                                  888                        "Y88P"
*/



suspend fun sampling(
    tau: List<Double>,
    initial: DoubleArray,
    f_integrated: (DoubleArray, Double) -> DoubleArray
): List<DoubleArray> {
    return tau.pmap { f_integrated(initial, it) }
}