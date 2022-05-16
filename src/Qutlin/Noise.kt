package Qutlin

import Plot
import kotlinx.atomicfu.atomic
import org.hipparchus.complex.Complex
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
    initialSpacing: Double, 
    val wnVariance: Double = 1.0, 
    val maxTime: Double = time
) {
    companion object {
        private val unit = fun(_: Double) = 1.0

        // * make sure the seed is new in every run
        var seed = atomic(LocalDateTime.now().nano.toLong()/100)
    }

    private val fftTransformer = FastFourierTransformer(DftNormalization.STANDARD)
    val N: Int
    val realSpacing: Double
    val omegaFreq: List<Double>
    lateinit var amplitudes: ComplexArray
    lateinit var values: DoubleArray

    init {
        // ? Generate a List of `omegaFreq` with 2^N elements
        val n = time/initialSpacing
        N = round(pow(2.0, ceil(log2(n)).toInt())).toInt()
        println("Noise of time $time initialized with $N~$n frequencies")
        realSpacing = time/N.toDouble()
        // * omegaFreq in units of TAU
        omegaFreq = List(N) {
            if(it <= N/2) TAU/time * it
            else -TAU/time * (N-it)
        }
    }

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
    fun generate(envelope: (Double) -> Double = unit, rescaleWN: Boolean = true, targetVariance: Double = -1.0) {

        val localdatetime = LocalDateTime.now().nano.toLong()/100
        seed.plusAssign(localdatetime)
        println("seed = $seed")
        // * generator for gaussian normalized random numbers (0 mean, 1 std)

        val generator = GaussianRandomGenerator(RandomDataGenerator(seed.value))

        val factor = if (rescaleWN) sqrt(wnVariance/realSpacing) else 1.0
        val whiteNoise = List(N) { generator.nextNormalizedDouble() * factor }

        amplitudes = fftTransformer.transform(whiteNoise.toDoubleArray(), TransformType.FORWARD )

        amplitudes = amplitudes.zip(omegaFreq) {a, o -> envelope(o) * a}.toComplexArray()

        val cvalues = fftTransformer.transform(amplitudes, TransformType.INVERSE )
        values = cvalues.map {it.real}.toDoubleArray()

        if(targetVariance > 0.0) {
//            val variance = values.map { it.real }.variance()
            val variance = values.toList().variance()
            val f = sqrt(targetVariance / variance)
            values = values.map { it * f }.toDoubleArray() //.toComplexArray()
        }
        println("noise generated with mean ${values.toList().average()} and std ${values.toList().std()}")

        val maxN = ceil(maxTime / realSpacing).toInt() + 1
        values = values.take(maxN).toDoubleArray()//.toComplexArray()

    }


    /**
     * Generate interpolated values from generated noise data.
     */
    private fun interpolate(t: Double): Double {
        val tt = t/realSpacing
        val index0 = floor(tt).toInt()
        val index1 = floor(tt + 1.0).toInt()
        val tr = tt - index0.toDouble()

        val N = values.size

        if(index0 >= N) return values[values.lastIndex]
        if(index1 >= N) return values[index0]
        if(index0 < 0) return values[0]

        return values[index0] * (1.0 - tr) + values[index1] * tr
    }

    /**
     * Makes `Noise` callable, e.g. `val noise = Noise(...)` can be used as `noise(t)`
     */
    operator fun invoke(t: Double) = this.interpolate(t)//.real



    override fun toString(): String {
        val realValues = values//.map { it.real }

        val avg = realValues.average()
        val std = realValues.toList().std()
        val variance = realValues.toList().variance()

        return """
            Noise
                avg         : $avg
                std         : $std
                var         : $variance
                realSpacing : $realSpacing
                var * delta : ${std*std*realSpacing}
            """.trimIndent()
    }
}






fun plotNoise(noise: Noise, take: Int = -1, plotFrequencies: Boolean = true, plotTimeSeries: Boolean = true) {

    if(plotFrequencies) {
        val p0 = Plot(
            Plot.DataSet(
                noise.omegaFreq.drop(1).zip(noise.amplitudes.drop(1)) { f, a -> listOf(f, a.real) },
                color = Plot.colorPalette[0],
                symbol = Plot.Symbols.None
            ), xScale = Mapper.Companion.Scales.Linear
        )
        p0.addDataSet(
            Plot.DataSet(
                noise.omegaFreq.drop(1).zip(noise.amplitudes.drop(1)) { f, a -> listOf(f, a.imaginary) },
                color = Plot.colorPalette[1],
                symbol = Plot.Symbols.None
            )
        )
    }

    if(plotTimeSeries) {
        var nRe = noise.values.mapIndexed { i, n -> listOf(i.toDouble(), n) }
//        var nRe = noise.values.mapIndexed { i, n -> listOf(i.toDouble(), n.real) }
//        var nIm = noise.values.mapIndexed { i, n -> listOf(i.toDouble(), n.imaginary) }

        if (take > 0) {
            nRe = nRe.take(take)
//            nIm = nIm.take(take)
        }


        val p1 = Plot(
            Plot.DataSet(
                nRe,
                color = Plot.colorPalette[0],
                symbol = Plot.Symbols.None
            )
        )
//        p1.addDataSet(
//            Plot.DataSet(
//                nIm,
//                color = Plot.colorPalette[1],
//                symbol = Plot.Symbols.None
//            )
//        )
    }
}






fun main() {
    fftTest()
}




fun fftTest() {
    val seed = LocalDateTime.now().nano.toLong()
    val generator = GaussianRandomGenerator(RandomDataGenerator(seed))
    val fftTransformer = FastFourierTransformer(DftNormalization.STANDARD)

    val N = pow(2.0,10).toInt()
    val T = 1.0
    val t = List(N) { it.toDouble()/N.toDouble() * T}

    val freq = 8

//    val eta = List(N) { generator.nextNormalizedDouble() }
    val eta = t.map { cos(it * freq * TAU - TAU/8.0) }
    println("eta: var = ${eta.variance()}, avg = ${eta.average()}")

    val p0 = Plot(Plot.DataSet(
        t.zip(eta){ it, ieta -> listOf(it,ieta)},
        color = Plot.colorPalette[0],
        symbol = Plot.Symbols.None
    ), xScale = Mapper.Companion.Scales.Linear)


    val frequencies = List(N) {
        if(it <= N/2) TAU/T * it
        else -TAU/T * (N - it)
    }
    val amplitudes = fftTransformer.transform( eta.toDoubleArray(), TransformType.FORWARD )
//    println(amplitudes.toList())

    val p1 = Plot(Plot.DataSet(
        frequencies.zip(amplitudes) {f,a -> listOf(f, a.real)},
        color = Plot.colorPalette[0],
        symbol = Plot.Symbols.None
    ))
    p1.addDataSet(Plot.DataSet(
        frequencies.zip(amplitudes) {f,a -> listOf(f, a.imaginary)},
        color = Plot.colorPalette[1],
        symbol = Plot.Symbols.None
    ))

    val etaBack = fftTransformer.transform( amplitudes, TransformType.INVERSE ).map { it.real }
    println("etaBack: var = ${etaBack.variance()}, avg = ${etaBack.average()}")
}
