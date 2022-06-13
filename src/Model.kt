import Qutlin.*
import org.hipparchus.util.FastMath.*

/**
 * The abstract Model defines properties that every simulation run of a concrete setup should have, such as
 * the final time `tf`, the dimensions of the Hilbert space `dimensions`, the overhang of the time variable before and
 * after the targeted time `t ∈ [0-overhang,tf+overhang]`, as well as preparation and measurement unitaries
 * `U_p` & `U_m` when using the _"generalized"_ strategy.
 *
 * The noise-free Hamiltonian `H_0` and the noisy Hamiltonian `H_η`, and the `collapse` operators are to be initialized
 * by the concrete implementation in `build()`, which is to be overridden.
 *
 * The `collapse` operator is given as a pair of the scalar prefactor `Γ` and a function `D(t)` returning the collapse
 * operator at time `t`.
 */
abstract class Model(
    val initial: Int, // ? the initially prepare eigenstate
    val tf: Double,
    val dimensions: Int,
    var overhang: Double = 0.0,
    val U_p: Operator? = null,
    val U_m: Operator? = null,
) {
    abstract fun build()
    abstract val maxIntegrationStep: Double

    var H_η: ((Double) -> Operator)? = null
    var H_0: ((Double) -> Operator)? = null
    var collapse: Pair<Double, (Double) -> Operator>? = null

    operator fun invoke(t: Double) = H_η!!.invoke(t)

    fun free() {
        H_0 = null
        H_η = null
    }
}


// ########   #######  ########
// ##     ## ##     ## ##     ##
// ##     ## ##     ## ##     ##
// ##     ## ##     ## ##     ##
// ##     ## ##  ## ## ##     ##
// ##     ## ##    ##  ##     ##
// ########   ##### ## ########

class DoubleQuantumDotModel(
    initial: Int,
    tf: Double,
    δbz: Double = 1.0 * _μeV / _ħ, // ?  ~1.5/ns
    Ω: Double = 20.0 * _μeV / _ħ,
    useShapedPulse: Boolean = false,
    useSmoothPulse: Boolean = true,
    τ: Double = 5.0 * _ns,
    ε_max: Double = 1700 * _μeV / _ħ,
    ε_min: Double = -200 * _μeV / _ħ,
//    σ: Double = 1.0 * _μeV / _ħ,
//    τ_c: Double = 1.0 * _ns,
    Γ: Double = 0.0 / _ns,
    noiseType: NoiseType,
) : DonorDotModel(
    initial,
    tf,
    4 * δbz, // ?  ~6.0/ns
    Ω,
    useShapedPulse,
    useSmoothPulse,
    τ,
    ε_max,
    ε_min,
//    σ,
//    τ_c,
    Γ,
    noiseType,
) {}

// ########  ########
// ##     ## ##     ##
// ##     ## ##     ##
// ##     ## ##     ##
// ##     ## ##     ##
// ##     ## ##     ##
// ########  ########

open class DonorDotModel(
        initial: Int,
        tf: Double,
        private val a: Double = TAU * 100.0 * 1e6 / _s, // ?  ~0.6/ns
        private val Ω: Double = 20.0 * _μeV / _ħ,
        private val useShapedPulse: Boolean = false,
        private val useSmoothPulse: Boolean = true,
        private val τ: Double = 5.0 * _ns,
        private val ε_max: Double = 1700 * _μeV / _ħ,
        private val ε_min: Double = -200 * _μeV / _ħ,
        private val Γ: Double = 0.0 / _ns,
        private val noiseType: NoiseType,
) : Model(initial, tf, 3) {

    var pulse: List<Pair<Double, Double>>? = null

    override val maxIntegrationStep : Double = π2 / max(noiseType.ω_sampling, abs(ε_max), abs(ε_min), Ω).toDouble() * 0.1

    override fun build() {
        val deltaNm = 1.0 // ? depending on the nuclear spin transition
        val δbz = -a / 2.0 * 0.5 * deltaNm

        var ε = if (useShapedPulse) {
            if (pulse == null) {
                // ? derivative of the detuning depending on the detuning `dε/dt[ε(t)]`
                fun dε(ε_: Double): Double {
                    val β = sqrt(ε_ * ε_ + Ω * Ω)
                    val x = β.e(3) * sqrt(1.0 + ε_ / β)
                    val y = a * a + 2 * Ω * Ω + 4 * ε_ * ε_ + ε_ / β * (a * a - 4 * Ω * Ω - 4 * ε_ * ε_)
                    val z = 2 * a * Ω * Ω * (ε_ - 3 * β)
                    return -x * y.e(3.0 / 2.0) / z
                }

                pulse = calculatePulseShape(Pair(ε_min, ε_max), ::dε)
            }

            fun(t: Double) = clamp(
                ε_min,
                interpolate(pulse!!, pulse!!.last().first * clamp(0.0, t, tf) / tf, Pair(ε_min, ε_max)),
                ε_max,
            )
        } else { // ? linear pulse
            fun(t: Double) = clamp(
                ε_min,
                ε_max + clamp(0.0, t, tf) / tf * (ε_min - ε_max),
                ε_max,
            )
        }

        if (useSmoothPulse) {
            overhang = 3.0 * τ

            val sm = Smooth(τ)
            val dt = min(τ / 10.0, tf / 1000.0)

            val start = -overhang
            val end = tf + overhang

            val smoothEps = List(((end - start) / dt + 1).toInt()) {
                Pair(
                    it * dt + start,
                    sm(it * dt + start) { t -> ε(t) }
                )
            }

            ε = fun(t: Double) = interpolate(smoothEps, t, Pair(ε_min, ε_max))
        }


        // * the highest frequency is given by max(|ε_max|, |ε_min|) !
        val ω_max = max(abs(ε_max), abs(ε_min), Ω, noiseType.ω_sampling)
        val t_total = if(useShapedPulse) tf+6*τ else tf

        val η = Noise(
            t_total,
            ω_max * 10.0 * 10.0,
            π2/t_total * 0.1,
            ω_max * 10.0,
        )
        η.generate(noiseType)


        // * noisy Hamiltonian
        H_η = if(useShapedPulse) fun(t: Double) = Operator(
            Pair(3, 3), complexArrayOf(
                0.0, δbz, Ω / 2.0,
                δbz, 0.0, 0.0,
                Ω / 2.0, 0.0, ε(t) + η(t+3*τ),
            )
        ) else fun(t: Double) = Operator(
            Pair(3, 3), complexArrayOf(
                0.0, δbz, Ω / 2.0,
                δbz, 0.0, 0.0,
                Ω / 2.0, 0.0, ε(t) + η(t),
            )
        )

        // * noise-free Hamiltonian
        H_0 = fun(t: Double) = Operator(
            Pair(3, 3), complexArrayOf(
                0.0, δbz, Ω / 2.0,
                δbz, 0.0, 0.0,
                Ω / 2.0, 0.0, ε(t),
            )
        )


        if (Γ != 0.0) {
            collapse = Pair(Γ,
                fun(t: Double): Operator {
                    // ? changed from H_η to H_0 since that's what we write in the text
                    val sys = H_0!!(t).eigenSystem()
                    val singletNeg = sys[0].second.normalized()
                    val singletPos = sys[2].second.normalized()
                    return singletNeg / singletPos.bra()
                })
        }

    }
}

//  ######   ######
// ##    ## ##    ##
// ##       ##
// ##       ##   ####
// ##       ##    ##
// ##    ## ##    ##
//  ######   ######

class ConstantGapModel(
    initial: Int,
    tf: Double,
    private val gap: Double = 10.0,
    val noiseType: NoiseType,
    initialTransformation: Operator? = null,
    finalTransformation: Operator? = null,
) : Model(
    initial,
    tf,
    2,
    U_p = initialTransformation,
    U_m = finalTransformation,
) {
    override val maxIntegrationStep : Double = π2 / max(noiseType.ω_sampling, gap) * 0.1
    override fun build() {
        val ε = fun(t: Double) = gap * cos(π * clamp(0.0, t, tf) / tf)
        val Ω = fun(t: Double) = gap * sin(π * clamp(0.0, t, tf) / tf)


        // The purpose of the cutoff frequency is to generate a "smooth" signal η.
        // If there is no cutoff frequency, the high frequency components lead to
        // "jaggy" lines of η(t) and the integrator will be slowed down unnecessarily.
        // By choosing a cutoff frequency that is much higher than the relevant
        // frequency (in this case, the energy gap between the states), we cover
        // the relevant physical aspects, and we can neglect higher frequencies.
        // Since E=ħω, the cutoff frequency is given by
        //     ω_max = gap / ħ
        // Using h=1 and ħ=h/2π, we get
        //     ω_max = 2π gap
        //
        // On the other hand, we want to sample frequencies much higher than the
        // width of the Lorenzian governing the Ornstein-Uhlenbeck noise
        //     S(ω) = 2σ²γ/(γ²+ω²)

        val ω_max = max(gap, noiseType.ω_sampling)

        val η = Noise(
            tf,
            ω_max * 10.0 * 10.0,
            π2/tf * 0.1,
            ω_max * 10,
        )
        η.generate(noiseType)

        H_η = fun(t: Double) = Operator(
            Pair(2, 2),
            complexArrayOf(
                ε(t) + η(t), Ω(t),
                Ω(t), -(ε(t) + η(t))
            ) * 0.5
        )

        H_0 = fun(t: Double) = Operator(
            Pair(2, 2),
            complexArrayOf(
                ε(t), Ω(t),
                Ω(t), -ε(t)
            ) * 0.5
        )
    }
}

// ##       ########
// ##            ##
// ##           ##
// ##          ##
// ##         ##
// ##        ##
// ######## ########

class LandauZenerModel(
    initial: Int,
    tf: Double,
    private val Ω: Double = 10.0,
    private val ε0: Double = -10.0 * Ω,
    private val ε1: Double =  10.0 * Ω,
    private val useShapedPulse: Boolean = false,
    private val noiseType: NoiseType,
    initialTransformation: Operator? = null,
    finalTransformation: Operator? = null,
) : Model(
    initial,
    tf,
    2,
    U_p = initialTransformation,
    U_m = finalTransformation,
) {
    override val maxIntegrationStep : Double = π2 / max(noiseType.ω_sampling, abs(ε0), abs(ε1), Ω).toDouble() * 0.1
    override fun build() {

        val ε = if (!useShapedPulse) {
            fun(t: Double) = ε0 + (ε1-ε0) * t/tf
        } else {
            val δ = -1.0/(Ω*tf) * ( ε1/sqrt(Ω*Ω + ε1*ε1) - ε0/sqrt(Ω*Ω + ε0*ε0) )
            val t0 = -1.0/(Ω*δ) * ε0/sqrt(Ω*Ω + ε0*ε0)

            fun(t: Double): Double {
                val tt = clamp(0.0, t, tf)
                val x = (tt+t0)*Ω*δ
                return -x * Ω / sqrt(1 - x*x)
            }
        }

        val ω_max = max(π2/tf ,abs(ε0), abs(ε1), Ω, noiseType.ω_sampling)

        val η = Noise(
            tf,
            ω_max * 10.0, // ? safety factor 10.0 to have high enough time resolution - Fehse, 2022-06-13
            π2/tf * 0.1,       // ? safety factor 10.0 to include all relevant frequencies during the evolution - Fehse, 2022-06-13
            ω_max * 10.0,      // ? safety factor 10.0 to include all relevant frequencies during the evolution - Fehse, 2022-06-13
        )
        η.generate(noiseType)

        H_η = fun(t: Double) = Operator(
            Pair(2, 2),
            complexArrayOf(
                ε(t) + η(t), Ω,
                Ω, -(ε(t) + η(t))
            ) * 0.5
        )

        H_0 = fun(t: Double) = Operator(
            Pair(2, 2),
            complexArrayOf(
                ε(t), Ω,
                Ω, -ε(t)
            ) * 0.5
        )
    }
}

fun main() {
    val cg = ConstantGapModel(
        1,
        1.0,
        1.0/10000.0,
        OUNoise(1.0,1.0),
    )
    cg.build()

    println(cg.H_0!!(0.0))
    val eigenvalues = cg.H_0!!(0.0).eigenSystem()
    eigenvalues.forEach{ println("${it.first} -> ${it.second.str()}") }
}
