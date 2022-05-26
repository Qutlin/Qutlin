import Qutlin.*
import org.hipparchus.util.FastMath.*


abstract class NoiseType {
    abstract val name: String
    abstract val ω_max: Double
    abstract fun envelope(ω: Double): Double
    abstract fun variance(ω_min: Double, ω_max: Double): Double
}



data class WhiteNoise(
    val σ: Double,
    override val ω_max: Double
) : NoiseType() {
    override val name: String = "WN"

    override fun envelope(ω: Double) = σ*σ

    override fun variance(ω_min: Double, ω_max: Double): Double = 0.0 // TODO correct variance

    override fun toString() = name
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
data class OUNoise(
    val σ: Double,
    val γ: Double,
): NoiseType() {
    override val name: String = "OU"
    override val ω_max: Double = π2*γ

    override fun envelope(ω: Double) =
        σ*σ * 2.0 * γ / (γ*γ + ω*ω)

    override fun variance(ω_min: Double, ω_max: Double): Double {
        return 0.0
    }


    override fun toString() = "${name}[γ%.1e σ%.1e]".format(γ, σ)

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
data class f_inv_Noise(
    val S0: Double,
    val ω_min: Double,
    override val ω_max: Double,
): NoiseType() {
    override val name: String = "1_f"
    override fun envelope(ω: Double) =
        if (abs(ω) < ω_min || abs(ω) > ω_max) 0.0
        else S0/abs(ω)

    override fun variance(ω_min: Double, ω_max: Double): Double {
        return 0.0
    }

    override fun toString() = "${name}[S0%.1e ω_min%.1e ω_max%.1e]".format(S0, ω_min, ω_max)
}











fun main() {
    noiseTypesTest()
//    fnoiseTest()
}





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
fun noiseTypesTest() {
    val time = 1.0
    val N = 20000
    val dt = time/N.toDouble()

    println("Ornstein-Uhlenbeck noise -------------------------------------")
    val γ = 40.0
    val variance = 1.0
    val σ = sqrt(variance)
    println("σ = $σ")

//    val ouSetup = WhiteNoise(σ, π2/dt)
    val ouSetup = OUNoise(σ = σ, γ = γ)
    println("$ouSetup")

    val noise = Noise(time, π2/dt, π2/time*1.0, π2/dt*1.0)
    noise.generate(ouSetup, 1.0/sqrt(dt))

    plotNoise(noise)
    println(noise)
//    noise.values.take(100).forEach { println("noise values: $it") }
    val σ2_noise = noise.values.variance()

//    // * testing interpolation
//    val ts = List(10000) { 0.005 * time/10000.0 * it.toDouble() }
//    val nvals = ts.map { ouNoiseData(it).real }
//    val nsPlot = Plot(
//        Plot.DataSet(ts.zip(nvals){t,v -> listOf(t,v)}, symbol = Plot.Symbols.None)
//    )

    // * calculate correlation function
    val autoCorrelation = noise.values.pAutoCorrelation(normalized = false)
    val timeValues = List(autoCorrelation.size) {it.toDouble() * dt}

    // * plot correlation functions
    val dataPairs = timeValues.zip(autoCorrelation) { t, c -> listOf(t,abs(c))}
//    dataPairs.take(100).forEach { println(it) }

    val plot = Plot(
        Plot.DataSet(
            dataPairs,
            symbol = Plot.Symbols.None,
        ),
//        xScale = Mapper.Companion.Scales.Log10,
        yScale = Mapper.Companion.Scales.Log10,
        minX = 0.0,
        maxX = 1.0/γ * 4.0, // ? this corresponds to 4 * tau_c (correlation time)
        maxY = variance,
        minY = 1e-2 * variance,
    )

    // * theoretical correlation curve
    val theo = timeValues.map { listOf(it, variance * exp(-γ * it)) }
    plot.addDataSet(Plot.DataSet(theo, symbol = Plot.Symbols.None, color = Plot.colorPalette[1]))
}