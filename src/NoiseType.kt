import Qutlin.*
import org.hipparchus.util.FastMath.*


abstract class NoiseType {
    abstract val name: String
    abstract val initialSpacing: Double
    abstract val wnDeltaRate: Double
    abstract val initialVariance: Double
    abstract fun envelope(ω: Double): Double
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
            σ * sqrt( 2.0 * γ / (γ*γ + ω*ω) )
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
}