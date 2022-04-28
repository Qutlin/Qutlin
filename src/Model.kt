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
    val maxIntegrationStep: Double,
    val dimensions: Int,
    var overhang: Double = 0.0,
    val U_p: Operator? = null,
    val U_m: Operator? = null,
) {
    abstract fun build()
    lateinit var H_η: (Double) -> Operator
    lateinit var H_0: (Double) -> Operator
    var collapse: Pair<Double, (Double) -> Operator>? = null

    operator fun invoke(t: Double) = H_η(t)
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
    maxIntegrationStep: Double,
    δbz: Double = 1.0 * _μeV / _ħ, // ?  ~1.5/ns
    Ω: Double = 20.0 * _μeV / _ħ,
    useShapedPulse: Boolean = false,
    useSmoothPulse: Boolean = true,
    τ: Double = 5.0 * _ns,
    ε_max: Double = 1700 * _μeV / _ħ,
    ε_min: Double = -200 * _μeV / _ħ,
    σ: Double = 1.0 * _μeV / _ħ,
    τ_c: Double = 1.0 * _ns,
    Γ: Double = 0.0 / _ns,
) : DonorDotModel(
    initial,
    tf,
    maxIntegrationStep,
    4 * δbz, // ?  ~6.0/ns
    Ω,
    useShapedPulse,
    useSmoothPulse,
    τ,
    ε_max,
    ε_min,
    σ,
    τ_c,
    Γ,
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
        maxIntegrationStep: Double,
        private val a: Double = TAU * 100.0 * 1e6 / _s, // ?  ~0.6/ns
        private val Ω: Double = 20.0 * _μeV / _ħ,
        private val useShapedPulse: Boolean = false,
        private val useSmoothPulse: Boolean = true,
        private val τ: Double = 5.0 * _ns,
        private val ε_max: Double = 1700 * _μeV / _ħ,
        private val ε_min: Double = -200 * _μeV / _ħ,
        private val σ: Double = 1.0 * _μeV / _ħ, // TODO correct values!?
        private val τ_c: Double = 1.0 * _ns,
        private val Γ: Double = 0.0 / _ns,
) : Model(
    initial, tf, maxIntegrationStep, 3
) {
//    companion object {
//        // ? shaped pulse
//        var pulse: List<Pair<Double, Double>>? = null
//    }
    var pulse: List<Pair<Double, Double>>? = null

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


        val γ = 1.0 / τ_c
        // * the highest frequency is given by max(|ε_max|, |ε_min|) !
        val ΔEmax = max(abs(ε_max), abs(ε_min))
        val cutoff = max(2.0 * π * ΔEmax, γ) * 10.0
        val noiseType = OUNoise(
            σ,
            γ,
            cutoff,
            initialSpacing = 2 * π / cutoff,
        )

        val η = if(useShapedPulse) Noise(
            max(tf+6*τ, 10.0*τ_c),
            min(noiseType.initialSpacing, (tf+6*τ) / 100.0),
            noiseType.wnDeltaRate,
            tf + 6*τ
        ) else Noise(
            max(tf, 10.0*τ_c),
            min(noiseType.initialSpacing, tf / 100.0),
            noiseType.wnDeltaRate,
            tf
        )
        η.generate(noiseType::envelope, rescaleWN = true)


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
                    val sys = H_0(t).eigenSystem()
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
    maxIntegrationStep: Double,
    private val gap: Double = 10.0,
    private val σ: Double = 1.0,
    private val γ: Double = 1.0,
    initialTransformation: Operator? = null,
    finalTransformation: Operator? = null,
) : Model(
    initial,
    tf,
    maxIntegrationStep,
    2,
    U_p = initialTransformation,
    U_m = finalTransformation,
) {
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
        val cutoffFrequency = max(2.0 * π * gap, γ) * 10

        val noiseType = OUNoise(
            σ,
            γ,
            cutoffFrequency,
            initialSpacing = 2 * π / cutoffFrequency
        )
        val τ_c = 1.0/γ;
        println("tf = $tf, τ_c = $τ_c, γ = $γ")
        val η = Noise(
            max(tf, 10 * τ_c), // ? total time; make sure that the smallest, non-zero frequency ω_1 = 2π/tf is small compared to the width of the spectrum given by γ = 1/τ_c; The smalles frequency is given by ω_1 = 
            min(noiseType.initialSpacing, tf/100.0), // ? see above; make sure that the time-resolution of the noise is high enough to result in a smooth curve
            noiseType.wnDeltaRate, // ? wnVariance
            tf, // ? maximum time to consider
        )
        η.generate(noiseType::envelope, rescaleWN = true)

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
    maxIntegrationStep: Double,
    private val gap: Double = 10.0,
    private val σ: Double = 1.0,
    private val γ: Double = 1.0,
    private val ε_max: Double = 10.0 * gap,

    private val useShapedPulse: Boolean = false,

    initialTransformation: Operator? = null,
    finalTransformation: Operator? = null,
) : Model(
    initial, 
    tf, 
    maxIntegrationStep, 
    2,
    U_p = initialTransformation,
    U_m = finalTransformation,
) {
    override fun build() {

        val Ω = gap

        val ε = if (!useShapedPulse) {
            fun(t: Double) = (1.0 - 2.0 * t / tf) * ε_max
        } else {
            val α = sqrt(Ω * Ω + ε_max * ε_max) * tf * Ω / (2 * ε_max)

            fun(t: Double): Double {
                val tt = clamp(0.0, t, tf) - tf / 2.0
                return -tt * Ω * Ω / sqrt(α * α - tt * tt * Ω * Ω)
            }
        }

        val cutoff = max(2.0 * π * gap, γ) * 10.0
        val τ_c = 1.0/γ

        val noiseType = OUNoise(
            σ,
            γ,
            cutoff,
            initialSpacing = 2 * π / cutoff,
        )
        val η = Noise(
            max(tf, 10.0 * τ_c),
            min(noiseType.initialSpacing, tf / 100.0),
            noiseType.wnDeltaRate,
            tf,
        )
        η.generate(noiseType::envelope, rescaleWN = true)


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
    )
    cg.build()

    println(cg.H_0(0.0));
    val eigenvalues = cg.H_0(0.0).eigenSystem()
    eigenvalues.forEach{ println("${it.first} -> ${it.second.str()}") }
}
