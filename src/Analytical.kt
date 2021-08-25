//import Qutlin.*
//import org.hipparchus.analysis.UnivariateFunction
//import org.hipparchus.analysis.integration.BaseAbstractUnivariateIntegrator.DEFAULT_MAX_ITERATIONS_COUNT
//import org.hipparchus.analysis.integration.RombergIntegrator
//import org.hipparchus.analysis.solvers.BracketingNthOrderBrentSolver
//import org.hipparchus.complex.Complex
//import org.hipparchus.util.FastMath
//import org.hipparchus.util.FastMath.cos
//import org.hipparchus.util.FastMath.sin
//import java.lang.Error
//import kotlin.math.exp
//import kotlin.math.pow
//import kotlin.math.sqrt
//
//fun main() {
//    calculateAnalyticalError(16.6 * _ns)
//}
//
//
//fun noiseComponent(
//    tf: Double = 16.6 * _ns,
//    smoothStd: Double = 5.0 * _ns,
//    spectrum: (Double) -> Double,
//    sigma: Double,
//    cx: (Double) -> Double,
//    cy: (Double) -> Double,
//    B: (Double) -> Double,
//    Bd: (Double) -> Double,
//    Bdd: (Double) -> Double,
//): Double {
//    val integrator = RombergIntegrator(1e-6, 1e-10, 3, 32)
//    val solver = BracketingNthOrderBrentSolver(1e-10, 1e-8, 5)
//
//
//    fun cPerp(t: Double) = cy(t) - I * cx(t)
//
//    fun innerSum(omega: Double): Double {
//        // * find zeros (T tildes)
//        // ! how to make sure to find all zeros?
//        fun findInbetween(zeroFun: (Double) -> Double, left: Double, right: Double): List<Double> {
//
//            val T0 = try {
//                solver.solve(100, zeroFun, left, right)
//            } catch (error: Error) {
//                // * solution not found
//                return emptyList()
//            }
//
//            val res = mutableListOf(T0)
//
//            res.addAll(findInbetween(left, T0))
//            res.addAll(findInbetween(T0, right))
//
//            return res
//        }
//
////        fun xi(Double: T) = erf
//
//        // * left term
//        fun zeroLeft(T: Double) = Bd(T) * T + B(T) + omega
//        val leftZeros = findInbetween(::zeroLeft, -1000 * tf, 1001 * tf)
//        val resLeft = leftZeros.map { zero ->
//            zero
//        }.sum()
//
//
//    }
//
//    fun omegaIntegralRe(omega: Double) = 1.0 / TAU * spectrum(omega) * innerSum(omega)
////    fun omegaIntegralIm(omega: Double) = 1.0 / TAU * spectrum(omega) * innerIntegralIm(omega)
//
//    val omegaIntegralRe = integrator.integrate(
//        DEFAULT_MAX_ITERATIONS_COUNT, ::omegaIntegralRe, -5 * sigma, 5 * sigma
//    )
//
//    return -2 * omegaIntegralRe
//}
//}
//
//
//
//
//
//
//
//fun calculateAnalyticalError(
//    tf: Double = 16.6 * _ns,
//    smoothStd: Double = 5.0 * _ns,
//    noiseTau: Double = 1.0 * _ns,
//) {
//    val a = 2.0 * FastMath.PI * 100.0 * 1e6 / _s
//    val deltaNm = 1.0
//    val dbz = -a / 2.0 * 0.5 * deltaNm
//    val om = 20.0 * _μeV / _ħ
//    val epsBounds = Pair(-200 * _μeV / _ħ, 1700 * _μeV / _ħ)
//
////    val smoothStd = 5.0 * _ns
//    val overhang = 3 * smoothStd
//    val start = -overhang
//    val end = tf + overhang
//
//    fun deps(ep: Double): Double {
//        val beta = FastMath.sqrt(ep * ep + om * om)
//        val x = beta.e(3) * FastMath.sqrt(1.0 + ep / beta)
//        val y = a * a + 2 * om * om + 4 * ep * ep + ep / beta * (a * a - 4 * om * om - 4 * ep * ep)
//        val z = 2 * a * om * om * (ep - 3 * beta)
//        return -x * y.e(3.0 / 2.0) / z
//    }
//
//    val pulse = calculatePulseShape(epsBounds, ::deps)
//    val pulse0 = pulse.first()
//    val pulse1 = pulse.last()
//
//    fun epsShaped(t: Double, tau: Double) = clamp(
//        epsBounds.first,
//        interpolate(pulse, pulse1.first * t / tau, epsBounds),
//        epsBounds.second
//    )
//
//    val sm = Smooth(smoothStd)
//    val dt = FastMath.min(smoothStd / 10.0, tf / 1000.0)
//    val smoothEps = List(((end - start) / dt + 1).toInt()) {
//        Pair(
//            it * dt + start,
//            sm(it * dt + start) { t -> epsShaped(clamp(0.0, t, tf), tf) }
//        )
//    }
//
//    fun smoothEps(t: Double) = interpolate(smoothEps, t, epsBounds)
//
//
//    fun err(sig: Double): Double {
//
//        val integrator = RombergIntegrator(1e-6, 1e-10, 3, 32)
//
//
//        fun r0(t: Double): Double {
//            val ep = smoothEps(t)
//            return sqrt(ep * ep + om * om)
//        }
//
//        val r0UF = UnivariateFunction { t -> r0(t) }
//
//        fun rho(t1: Double): Double {
//            if (t1 >= 0.0) {
//                if (t1 == 0.0) return 0.0
//                return integrator.integrate(DEFAULT_MAX_ITERATIONS_COUNT, r0UF, 0.0, t1)
//            } else {
//                return -integrator.integrate(DEFAULT_MAX_ITERATIONS_COUNT, r0UF, t1, 0.0)
//            }
//        }
//
//        var data = List(1000) {
//            listOf(it.toDouble(), rho((it / 500.0 - 1.0) * 100.0))
//        }
//
//        val p = Plot(Plot.DataSet(data, symbol = Plot.Symbols.None))
//
//        fun inner(omega: Double): Complex {
//            fun fRe(x: Double, om: Double): Double {
//                val r0 = r0(x)
//                return 1.0 / (r0 * r0) * cos(om * x - rho(x))
//            }
//
//            fun fIm(x: Double, om: Double): Double {
//                val r0 = r0(x)
//                return 1.0 / (r0 * r0) * sin(om * x - rho(x))
//            }
//
//            val nn = 10000
//            val dRe = List(nn) {
//                val x = (it / (nn.toDouble() / 2.0) - 1.0) * 20
//                listOf(x, fRe(x, 1.0))
//            }
//            val pp = Plot(Plot.DataSet(dRe, symbol = Plot.Symbols.None))
//
//
//            val fReUF = UnivariateFunction { x -> fRe(x, omega) }
//            val fImUF = UnivariateFunction { x -> fIm(x, omega) }
//
//            println("integrating Re omega = $omega")
//            val resRe = integrator.integrate(DEFAULT_MAX_ITERATIONS_COUNT, fReUF, -1e3, 1e3)
//            println("integrating Im omega = $omega")
//            val resIm = integrator.integrate(DEFAULT_MAX_ITERATIONS_COUNT, fImUF, -1e3, 1e3)
//
//            return resRe + I * resIm
//        }
//
//        val N = 4
//        data = List(N) {
//            listOf(it.toDouble(), inner((it / N.toDouble() * 2.0 - 1.0) * 100.0).abs())
//        }
//        val p2 = Plot(Plot.DataSet(data, symbol = Plot.Symbols.None))
//
//        fun spectrum(omega: Double) = exp(-(omega / sig).pow(2))
//
//        fun f(omega: Double) = 1.0 / TAU * spectrum(omega) * inner(omega).abs().pow(2)
//        val fUF = UnivariateFunction { om -> f(om) }
//
//        return integrator.integrate(DEFAULT_MAX_ITERATIONS_COUNT, fUF, -1e3, 1e3)
//    }
//
//    for (i in 0..100) {
//        println("$i")
//        println("$i : ${err(10.0.pow(i / 100.0 - 0.5))}")
//    }
//}