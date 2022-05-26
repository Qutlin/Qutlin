package Qutlin

import NoiseType
import Plot
import kotlinx.atomicfu.atomic
import org.hipparchus.random.GaussianRandomGenerator
import org.hipparchus.random.RandomDataGenerator
import org.hipparchus.random.UniformRandomGenerator
import org.hipparchus.transform.DftNormalization
import org.hipparchus.transform.FastFourierTransformer
import org.hipparchus.transform.TransformType
import org.hipparchus.util.FastMath.*
import java.time.LocalDateTime
import kotlin.math.log2


/**
 * Noise Class to generate Gaussian noise with a given envelope of the spectrum.
 * Use objects as functions, e.g. `val noise = Noise(...)` used as `noise(t)`.
 */
class Noise(
    val time: Double,
    val ω_sampling: Double,
    val ω_min: Double,
    val ω_max: Double,
    val ω_0: Double? = null,
) {
    companion object {
        // * make sure the seed is new in every run
        var seed = atomic(LocalDateTime.now().nano.toLong()/100)
    }

    val N: Int = (time * ω_sampling/π2).toInt()
    val dt: Double = time/N.toDouble()
    lateinit var values: DoubleArray

    val tf: Double = π2/ω_min
    val Nt: Int = (tf * ω_max/π2).toInt()

    /**
     * from the Hipparchus docs:
     *
     * forward transform: `yn =     ∑_k=0^N-1 x_k  exp(-2πi n k / N)`,
     *
     * inverse transform: `xk = 1/N ∑_n=0^N-1 y_n exp( 2πi n k / N)`
     *
     * in terms of `ω_k` and `t_n`: `a(ω_k) =     ∑_n=1^N η(t_n)   exp(-i ω_k t_n)`,
     *
     * inverse transform: `η(t_n)   = 1/N ∑_k=1^N a(ω_k) exp( i ω_k t_n)`
     *
     * @param envelope function in radian frequencies (`ω`)
     * */
    fun generate(noise_type: NoiseType, σ_w: Double = sqrt(1.0/dt)) {
        val seed = seed.addAndGet(LocalDateTime.now().nano.toLong()/100)
        println("seed = $seed")

        // * generator for uniform normalized random numbers (0 mean, 1 std)
        // * Hipparchus: Since it is a normalized random generator, it generates values
        // * from a uniform distribution with mean equal to 0 and standard deviation equal
        // * to 1. Generated values fall in the range [-√3, +√3].
        val u_generator = UniformRandomGenerator( RandomDataGenerator(seed))
        val isqrt3 = 1.0/sqrt(3.0)
        val g_generator = GaussianRandomGenerator(RandomDataGenerator(seed))

        values = DoubleArray(N)

        val Nω = (ω_max/ω_min).toInt()
        println("$Nω Fourier frequencies to sum over.")
        val σ_r = σ_w * sqrt(Nω.toDouble())

        for (k in 0..Nω) {
            if (k % (Nω/10) == 0)
                println("$k/$Nω")

            val α = π2 * 0.5*(u_generator.nextNormalizedDouble() * isqrt3 + 1) // see comment above
            val p =      0.5*(u_generator.nextNormalizedDouble() * isqrt3 + 1) // see comment above
            val r = sqrt(-2.0*log(1.0-p)) * σ_r

            val ωk = ω_min * k

            val s = sqrt(noise_type.envelope(ωk))
//            val z = r * exp(I*α)
            for (i in 0 until N) {
                // ? Re[z * exp(ωk*t)] * sqrt(S(ωk))
//                values[i] += (z.real*cos(ωk*(i*dt)) - z.imaginary*sin(ωk*(i*dt))) * s
                values[i] += r * cos(ωk*(i*dt) + α) * s
            }
        }

        val offset = if (ω_0 != null) g_generator.nextNormalizedDouble() * noise_type.variance(ω_0, ω_min) else 0.0
        val sqrt2 = sqrt(2.0) // ? since we're only summing from 0 to Nω, not from -Nω to Nω; see notes from 2022-05-26 - Fehse 2022-05-26
        for (i in 0 until Nt) {
            values[i] = values[i]/Nω.toDouble() * sqrt2 + offset
        }

        println("noise generated with mean ${values.average()} and std ${values.std()}")
    }


    /**
     * Generate interpolated values from generated noise data.
     */
    private fun interpolate(t: Double): Double {
        val tt = t/dt
        val index0 = floor(tt).toInt()
        val index1 = floor(tt + 1.0).toInt()
        val tr = tt - index0.toDouble()

        if(index0 >= Nt) return values[values.lastIndex]
        if(index1 >= Nt) return values[index0]
        if(index0 < 0)   return values[0]

        return values[index0] * (1.0 - tr) + values[index1] * tr
    }

    /**
     * Makes `Noise` callable, e.g. `val noise = Noise(...)` can be used as `noise(t)`
     */
    operator fun invoke(t: Double) = this.interpolate(t)//.real



    override fun toString(): String {
        val avg      = values.average()
        val variance = values.variance()

        return """
            Noise
                avg         : $avg
                std         : ${sqrt(variance)}
                var         : $variance
                var * delta : ${variance*dt}
            """.trimIndent()
    }
}






fun plotNoise(noise: Noise, take: Int = 0) {

    val v = if(take > 0) noise.values.take(take).toDoubleArray() else noise.values;
    val nRe = v.mapIndexed { i, n -> listOf(i.toDouble()*noise.dt, n) }

    val p1 = Plot(
        Plot.DataSet(
            nRe,
            color = Plot.colorPalette[0],
            symbol = Plot.Symbols.None
        )
    )
}




//
//
//fun main() {
//    fftTest()
//}
//
//
//
//
//fun fftTest() {
//    val seed = LocalDateTime.now().nano.toLong()
//    val generator = GaussianRandomGenerator(RandomDataGenerator(seed))
//    val fftTransformer = FastFourierTransformer(DftNormalization.STANDARD)
//
//    val N = pow(2.0,10).toInt()
//    val T = 1.0
//    val t = List(N) { it.toDouble()/N.toDouble() * T}
//
//    val freq = 8
//
////    val eta = List(N) { generator.nextNormalizedDouble() }
//    val eta = t.map { cos(it * freq * TAU - TAU/8.0) }
//    println("eta: var = ${eta.variance()}, avg = ${eta.average()}")
//
//    val p0 = Plot(Plot.DataSet(
//        t.zip(eta){ it, ieta -> listOf(it,ieta)},
//        color = Plot.colorPalette[0],
//        symbol = Plot.Symbols.None
//    ), xScale = Mapper.Companion.Scales.Linear)
//
//
//    val frequencies = List(N) {
//        if(it <= N/2) TAU/T * it
//        else -TAU/T * (N - it)
//    }
//    val amplitudes = fftTransformer.transform( eta.toDoubleArray(), TransformType.FORWARD )
////    println(amplitudes.toList())
//
//    val p1 = Plot(Plot.DataSet(
//        frequencies.zip(amplitudes) {f,a -> listOf(f, a.real)},
//        color = Plot.colorPalette[0],
//        symbol = Plot.Symbols.None
//    ))
//    p1.addDataSet(Plot.DataSet(
//        frequencies.zip(amplitudes) {f,a -> listOf(f, a.imaginary)},
//        color = Plot.colorPalette[1],
//        symbol = Plot.Symbols.None
//    ))
//
//    val etaBack = fftTransformer.transform( amplitudes, TransformType.INVERSE ).map { it.real }
//    println("etaBack: var = ${etaBack.variance()}, avg = ${etaBack.average()}")
//}
