import Qutlin.*
import org.hipparchus.util.FastMath.*


abstract class NoiseType {
    abstract val name: String
    abstract val initialSpacing: Double
    abstract val wnDeltaRate: Double
    abstract val min_tf: Double
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
    override val name: String = "OU"
    override val wnDeltaRate = 1.0
    override val min_tf = 10.0/γ

    override fun envelope(ω: Double) =
        if (abs(ω) < cutoff)
            σ * sqrt( 2.0 * γ / (γ*γ + ω*ω) )
        else 0.0

    override fun toString() =
        "${name}[g$γ s%.1e spacing%.1e co%.1e]".format(σ, initialSpacing, cutoff)

}




/**
 * 1/f Noise
 * 
 * stores the properties of the noise:
 * 
 * S0     : factor
 * 
 * ω0     : low-frequency cutoff
 * 
 * cutoff: 
 */
class f_inv_Noise(
    var S0: Double,
    var ω0: Double,
    var cutoff: Double,
    override val min_tf: Double = 2.0*π/ω0,
    override var initialSpacing: Double,
): NoiseType() {
    override val name: String = "1_f"
    override val wnDeltaRate = 1.0
//    override val minTime = 10.0/cutoff

    override fun envelope(ω: Double) =
        if (abs(ω) < ω0 || abs(ω) > cutoff) 0.0
        else sqrt(S0/abs(ω))

    override fun toString() =
        "${name}[S0$S0 ω0$ω0 spacing%.1e co%.1e]".format(initialSpacing, cutoff)

}










//
//fun main() {
////    noiseTypesTest()
//    fnoiseTest()
//}
//
//
//
//
//
//fun fnoiseTest() {
//    val time = 1e3
////    val spacing = 0.0025
//
//    println("1/f noise -------------------------------------")
//    val S0 = 2.67*1e6*_Hz*_Hz * 4.0*π*π * _ħ
//    val ω0 = 2.0*π/(60*60*_s) // period given by max tf
//    println("S0 = $S0, ω0 = $ω0")
//
//    val noiseType = f_inv_Noise(S0, ω0, cutoff = Double.MAX_VALUE, 4.0/ω0, 1.0e-3)
//    println("$noiseType")
//
//    val noise = Noise(time, noiseType.initialSpacing, noiseType.wnDeltaRate, time)
//    noise.generate(noiseType::envelope, rescaleWN = true)//, targetVariance = variance)
//
//    plotNoise(noise, take = 2000)
//    println(noise)
////    noise.values.take(100).forEach { println("noise values: ${it.real}") }
//    noise.values.take(100).forEach { println("noise values: $it") }
//
////    // * testing interpolation
////    val ts = List(10000) { 0.005 * time/10000.0 * it.toDouble() }
////    val nvals = ts.map { ouNoiseData(it).real }
////
////    val nsPlot = Plot(
////        Plot.DataSet(ts.zip(nvals){t,v -> listOf(t,v)}, symbol = Plot.Symbols.None)
////    )
//
//    // * calculate correlation function
//    val noiseValues = noise.values//.map{it.real}
//    val autoCorrelation = noiseValues.take(100000).pAutoCorrelation(normalized = false)
//    val timeValues = List(autoCorrelation.size) {it.toDouble() * noise.realSpacing}
//
//    // * plot correlation functions
//    val dataPairs = timeValues.zip(autoCorrelation) { t, c -> listOf(t,abs(c))}
//    dataPairs.take(100).forEach { println(it) }
//
//    val plot = Plot(
//        Plot.DataSet(
//            dataPairs,
//            symbol = Plot.Symbols.None,
//        ),
////        xScale = Mapper.Companion.Scales.Log10,
//        yScale = Mapper.Companion.Scales.Log10,
//        minX = 0.0,
//        maxX = 100 * _ns,
//        maxY = 10*S0,
//        minY = 1e-2 * S0,
//    )
//
//    // * theoretical correlation curve
////    val theo = timeValues.map { listOf(it, variance * exp(-γ * it)) }
////    plot.addDataSet(Plot.DataSet(theo, symbol = Plot.Symbols.None, color = Plot.colorPalette[1]))
//}
//fun noiseTypesTest() {
//    val time = 2000000.0
////    val spacing = 0.0025
//
//    println("Ornstein-Uhlenbeck noise -------------------------------------")
//    val γ = 0.00001
//    val variance = 1e-8
//    val σ = sqrt(variance)
//    println("σ = $σ")
//
//    val ouSetup = OUNoise(σ = σ, γ = γ, cutoff = Double.MAX_VALUE)
//    println("$ouSetup")
//
//    val varianceCorrection = 2.0 * atan(ouSetup.cutoff/γ)/PI
//    println("variance correction factor: $varianceCorrection")
//    println("target var * varCorr = ${variance * varianceCorrection}")
//
//    val noise = Noise(time, ouSetup.initialSpacing)
//    noise.generate(ouSetup::envelope, rescaleWN = true)//, targetVariance = variance)
//
//    plotNoise(noise, take = 200)
//    println(noise)
////    noise.values.take(100).forEach { println("noise values: ${it.real}") }
//    noise.values.take(100).forEach { println("noise values: $it") }
//
////    // * testing interpolation
////    val ts = List(10000) { 0.005 * time/10000.0 * it.toDouble() }
////    val nvals = ts.map { ouNoiseData(it).real }
////
////    val nsPlot = Plot(
////        Plot.DataSet(ts.zip(nvals){t,v -> listOf(t,v)}, symbol = Plot.Symbols.None)
////    )
//
//    // * calculate correlation function
//    val noiseValues = noise.values//.map{it.real}
//    val autoCorrelation = noiseValues.toList().pAutoCorrelation(normalized = false)
//    val timeValues = List(autoCorrelation.size) {it.toDouble() * noise.realSpacing}
//
//    // * plot correlation functions
//    val dataPairs = timeValues.zip(autoCorrelation) { t, c -> listOf(t,abs(c))}
//    dataPairs.take(100).forEach { println(it) }
//
//    val plot = Plot(
//        Plot.DataSet(
//            dataPairs,
//            symbol = Plot.Symbols.None,
//        ),
////        xScale = Mapper.Companion.Scales.Log10,
//        yScale = Mapper.Companion.Scales.Log10,
//        minX = 0.0,
//        maxX = 4.0/γ, // ? this corresponds to 4 * tau_c (correlation time)
//        maxY = variance,
//        minY = 1e-2 * variance,
//    )
//
//    // * theoretical correlation curve
//    val theo = timeValues.map { listOf(it, variance * exp(-γ * it)) }
//    plot.addDataSet(Plot.DataSet(theo, symbol = Plot.Symbols.None, color = Plot.colorPalette[1]))
//}