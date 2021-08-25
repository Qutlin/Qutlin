import Qutlin.*
import org.hipparchus.util.FastMath.*


abstract class NoiseType {
    abstract val name: String
    abstract val initialSpacing: Double
    abstract val wnDeltaRate: Double
    abstract val initialVariance: Double
    abstract fun envelope(ω: Double): Double
}


class GaussianSpectrumNoise(
    var correlationTime: Double = 1.0,
    override var initialVariance: Double = 1.0
) : NoiseType() {
    val omegaStd: Double
        get() = 2.0 / correlationTime // ? CHECK

    override val name = "gauss"
    override val wnDeltaRate: Double
        get() = sqrt(PI) * correlationTime * initialVariance // ? second sqrt already in noise generation // ? CHECK
    override val initialSpacing
        get() = 0.1 / omegaStd


    /**
     * omega: in TAU/s
     */
    override fun envelope(omega: Double): Double { // ? CHECK
        val x = omega/omegaStd
        return exp(-0.5 * x * x)
    }


    override fun toString() = "${name}[corr$correlationTime var%.3f]".format(initialVariance)
}

class FNoise(
    var highCutoff: Double = 1000.0 * _μeV,
    var lowCutoff: Double = 0.1 * _μeV,
    var amplitude: Double = 1.0,
    override var initialVariance: Double = 0.0 // TODO fix
) : NoiseType() {
    override val name: String = "1f"

    override val wnDeltaRate: Double = 1.0

    override val initialSpacing: Double
        get() = 0.1 / highCutoff

    override fun envelope(omega: Double): Double {
        if (omega <= lowCutoff || omega >= highCutoff) return 0.0
        return amplitude / sqrt(omega)
    }


    override fun toString() =
        "${name}[a$amplitude lowC$lowCutoff highC$highCutoff]"
}


class WhiteNoise: NoiseType() {
    override var initialSpacing = 0.01
    override val name: String = "white"
    override val wnDeltaRate: Double = 1.0
    override val initialVariance: Double = 1.0 // TODO fix
    override fun envelope(omega: Double): Double = 1.0
}











/**
 * Ornstein-Uhlenbeck Noise
 * 
 * stores the properties of the noise:
 * 
 * σ     : std deviation of the noise
 * 
 * γ     : inverse correlation time γ = 1/τ_c
 * 
 * cutoff: 
 */
class OUNoise(
    var σ: Double,
    var γ: Double,
    var cutoff: Double = 20.0 * γ,
    override var initialSpacing: Double = 2*PI/(20.0 * γ) * 0.1 // ? make sure time-resolution is high enough
): NoiseType() {
    override var initialVariance: Double = σ * σ
    override val name: String = "OU"
    override val wnDeltaRate = 1.0

    override fun envelope(ω: Double) =
        if (abs(ω) < cutoff)
            σ * Math.sqrt( 2.0 * γ / (γ*γ + ω*ω) )
        else 0.0

    override fun toString() =
        "${name}[g$γ s%.1e spacing%.1e co%.1e]".format(σ, initialSpacing, cutoff)

}













fun main() {
    noiseTypesTest()
}




fun noiseTypesTest() {
    val time = 2000000.0
//    val spacing = 0.0025

    println("Ornstein-Uhlenbeck noise -------------------------------------")
    val γ = 0.00001
    val variance = 1e-8
    val σ = sqrt(variance)
    println("σ = $σ")

    val ouSetup = OUNoise(σ = σ, γ = γ, cutoff = Double.MAX_VALUE)
    println("$ouSetup")

    val varianceCorrection = 2.0 * atan(ouSetup.cutoff/γ)/PI
    println("variance correction factor: $varianceCorrection")
    println("target var * varCorr = ${variance * varianceCorrection}")

    val noise = Noise(time, ouSetup.initialSpacing)
    noise.generate(ouSetup::envelope, rescaleWN = true)//, targetVariance = variance)

    plotNoise(noise, take = 200)
    println(noise)
    noise.values.take(100).forEach { println("noise values: ${it.real}") }

//    // * testing interpolation
//    val ts = List(10000) { 0.005 * time/10000.0 * it.toDouble() }
//    val nvals = ts.map { ouNoiseData(it).real }
//
//    val nsPlot = Plot(
//        Plot.DataSet(ts.zip(nvals){t,v -> listOf(t,v)}, symbol = Plot.Symbols.None)
//    )

    // * calculate correlation function
    val noiseValues = noise.values.map{it.real}
    val autoCorrelation = noiseValues.pAutoCorrelation(normalized = false)
    val timeValues = List(autoCorrelation.size) {it.toDouble() * noise.realSpacing}

    // * plot correlation functions
    val dataPairs = timeValues.zip(autoCorrelation) { t, c -> listOf(t,abs(c))}
    dataPairs.take(100).forEach { println(it) }

    val plot = Plot(
        Plot.DataSet(
            dataPairs,
            symbol = Plot.Symbols.None,
        ),
//        xScale = Mapper.Companion.Scales.Log10,
        yScale = Mapper.Companion.Scales.Log10,
        minX = 0.0,
        maxX = 4.0/γ, // ? this corresponds to 4 * tau_c (correlation time)
        maxY = variance,
        minY = 1e-2 * variance,
    )

    // * theoretical correlation curve
    val theo = timeValues.map { listOf(it, variance * exp(-γ * it)) }
    plot.addDataSet(Plot.DataSet(theo, symbol = Plot.Symbols.None, color = Plot.colorPalette[1]))


//
//
//
//    println("1/f noise ----------------------------------------------------")
//    val fNoise = FNoise()
////    noise = Noise(time, spacing = fNoise.spacing, wnVariance = fNoise.kappa)
//
//
//    println("gaussian -----------------------------------------------------")
//    val gNoise = GaussianSpectrumNoise()
//
//    gNoise.correlationTime = 0.05 * _ns
//    val gamma = 1.0 / _ns
//    gNoise.variance = gamma / (TAU * gNoise.correlationTime) // ? CHECK
//    println("target variance: ${gNoise.variance}")
//
//    val noise2 = Noise(time = time, spacing = gNoise.spacing, wnVariance = gNoise.kappa)
//    noise2.generate(gNoise::envelope)
//    plotNoise(noise2)
//    println(noise2)
//
//    val ns = noise2.noise.map{it.real}
//
//    val variance = ns.variance()
//    val avg = ns.average()
//    println("kappa = ${gNoise.kappa}")
//    println("kappa/delta = ${gNoise.kappa / noise2.delta}")
//    println("gamma = ${TAU * variance * gNoise.correlationTime}")
//    println("gamma * delta = ${TAU * variance * gNoise.correlationTime * noise2.delta}")
//
//
//    // * correlation time
//
//    val ct = ns.pAutoCorrelation()
//    val ts = List(ct.size) {(it+1).toDouble() * gNoise.spacing}
//    val tmp = ts.zip(ct) {t, c -> listOf(t,c)}
//    tmp.take(100).forEach { println(it) }
//    val p = Plot(
//        Plot.DataSet(
//            tmp,
//            symbol = Plot.Symbols.None
//        ),
//        xScale = Mapper.Companion.Scales.Log10
//    )
//
//
//
//    println("white -----------------------------------------------------")
//
//    val whiteNoise = WhiteNoise()
//    val noise3 = Noise(time = time, spacing = whiteNoise.spacing, wnVariance = whiteNoise.kappa)
//    noise3.generate(whiteNoise::envelope)
////    plotNoise(noise3)
////    println(noise3)
//
//
//    // ! check, if the noise is generated correctly
////    save("results/noise.csv", noise.noise.mapIndexed{i,n -> listOf(i.toDouble() * noise.dt, n.real)})
//
////    val generator = GaussianRandomGenerator(RandomDataGenerator(42))
////    save("results/rng.csv", List(1000){
////        listOf(it.toDouble(), generator.nextNormalizedDouble())
////    })
}