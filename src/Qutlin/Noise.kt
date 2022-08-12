package Qutlin

import NoiseType
import Plot
import kotlinx.atomicfu.atomic
import org.hipparchus.random.GaussianRandomGenerator
import org.hipparchus.random.RandomDataGenerator
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
//    val ω_sampling: Double,
//    val ω_min: Double? = null,
//    val ω_max: Double? = null,
) {
    companion object {
        // * make sure the seed is new in every run
        var seed = atomic(LocalDateTime.now().nano.toLong()/100)
    }

    // ? The FFT algorithms needs arrays with a size of a power of 2
    var Nt: Int = 0// = pow(2.0, ceil(log2(time * ω_sampling/π2)).toInt()).toInt()
    var dt: Double = 0.0// = time/Nt.toDouble()
    lateinit var values: DoubleArray

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
    fun generate(noise_type: NoiseType) {
        Nt = pow(2.0, ceil(log2(time * noise_type.ω_sampling/π2)).toInt()).toInt()
        dt = time/Nt.toDouble()
        val σ_wn = sqrt(1.0/dt)

        val seed = seed.addAndGet(LocalDateTime.now().nano.toLong()/100)
        println("seed = $seed")

        val g_rng = GaussianRandomGenerator(RandomDataGenerator(seed))

        println("time = $time, ω_sampling = $noise_type.ω_sampling")
        println("time*ω_sampling/π2 = ${time * noise_type.ω_sampling/π2}")
        println("log2(time*ω_sampling/π2) = ${log2(time * noise_type.ω_sampling/π2)}")
        println("Nt = $Nt, log2(Nt) = ${log2(Nt.toDouble())}")
        val whiteNoise = DoubleArray(Nt) { g_rng.nextNormalizedDouble() * σ_wn }
        val fftTransformer = FastFourierTransformer(DftNormalization.STANDARD)
        val amplitudes = fftTransformer.transform(whiteNoise, TransformType.FORWARD)

        for (i in amplitudes.indices) {
            val ω = if(i <= Nt/2) TAU/time*i else -TAU/time*(Nt-i)
//            if      (ω_max != null && abs(ω) > ω_max) amplitudes[i] = 0.0.toComplex()
//            else if (ω_min != null && abs(ω) < ω_min) amplitudes[i] = 0.0.toComplex()
//            else amplitudes[i] = amplitudes[i] *  sqrt(noise_type.envelope(ω))
            amplitudes[i] = amplitudes[i] * sqrt(noise_type.envelope(ω))
        }

        val cvalues = fftTransformer.transform(amplitudes, TransformType.INVERSE)
        values = DoubleArray(Nt) { cvalues[it].real }

        println("noise generated with Nt = $Nt, mean ${values.average()} and std ${values.std()}")
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
