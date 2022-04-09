import Qutlin.*
import Qutlin.ComplexMatrix.Companion.rotPauliX
import Qutlin.ComplexMatrix.Companion.rotPauliY
import kotlinx.coroutines.*
import org.hipparchus.util.FastMath.*
import java.util.concurrent.Executors
import kotlin.math.pow


fun main() {
    // constant_gap()
    // landau_zener()
    // donor_dot()
    double_quantum_dot()
}


/**
 * Calculate all the data for the constant-gap model.
 */
fun constant_gap() {

    // transfer error depending on variance `σ`
    completeSet_ConstantGap(
        x = linspace(-4.0, 1.0, 100).map { 10.0.pow(it) },
        samples = 20,
        gap = 10.0,
        γ = 1.0,
        σ = 0.1,
        tf = 100.0,
        variable = "σ",
        saveName = "2021 05 04",
        superAdiabatic = false,
    )

    // transfer error depending on `tf`
    completeSet_ConstantGap(
        x = linspace(-2.0, 3.0, 200).map { 10.0.pow(it) },
        samples = 20,
        gap = 10.0,
        γ = 1.0,
        σ = 0.1,
        tf = 100.0,
        variable = "tf",
        saveName = "2021 05 31",
        superAdiabatic = false,
    )

    // transfer error depending on rate `γ`
    completeSet_ConstantGap(
        x = linspace(-3.0, 4.0, 50).map { 10.0.pow(it) },
        samples = 20,
        variable = "γ",
        saveName = "2021 02 09"
    )
}

/**
 * Calculate all the data for the Landau-Zener model.
 */
fun landau_zener() {

    // transfer error depending on `tf`
    completeSet_LandauZener(
        x = concatenate(
            linspace(0.0, 1.0, 50).map { 10.0.pow(it) },
            linspace(1.0, 3.0, 50, skipFirst = true).map { 10.0.pow(it) },
        ),
        samples = 20,
        σ = 0.1,
        γ = 1.0,
        gap = 1.0,
        variable = "tf",
        useShapedPulse = false,
        saveName = "2021 05 04 LZ",
    )

    // transfer error depending on noise variance `σ`
    completeSet_LandauZener(
        x = linspace(-3.0, 1.0, 100).map { 10.0.pow(it) },
        samples = 20,
        variable = "σ",
        saveName = "2021 05 04 LZ",
    )

    // transfer error depending on rate `γ`
    completeSet_LandauZener(
        x = linspace(-4.0, 4.0, 100).map { 10.0.pow(it) },
        samples = 20,
        variable = "γ",
        saveName = "2021 05 04 LZ",
    )
}

/**
 * Calculate all the data for the donor-dot model.
 */
fun donor_dot() {

    // transfer error depending on `tf` for the linear pulse
    completeSet_DonorDot(
        x = linspace(0.0, 2.0, 100).map { 10.0.pow(it) / _ns },
        variable = "tf",
        samples = 20,

        a = TAU * 100.0 * 1e6 / _s,
        Ω = 20.0 * _μeV / _ħ,

        σ = 1.0 * _μeV / _ħ,
        τ_c = 1.0 * _ns,

        ε_max = 2000 * _μeV / _ħ,
        ε_min = -200 * _μeV / _ħ,
        τ = 4 * _ns,

        useShapedPulse = false,
        useSmoothPulse = true,

        saveData = true,
        saveName = "2021 09 26 DonorDot",
    )

    // transfer error depending on `tf` for the fast-QUAD pulse
    completeSet_DonorDot(
        x = linspace(0.0, 2.0, 100).map { 10.0.pow(it) / _ns },
        variable = "tf",
        samples = 20,

        a = TAU * 100.0 * 1e6 / _s,
        Ω = 20.0 * _μeV / _ħ,

        σ = 1.0 * _μeV / _ħ,
        τ_c = 1.0 * _ns,

        ε_max = 2000 * _μeV / _ħ,
        ε_min = -200 * _μeV / _ħ,
        τ = 4 * _ns,

        useShapedPulse = true,
        useSmoothPulse = true,

        saveData = true,
        saveName = "2021 09 26 DonorDot",
    )

    // transfer error depending on relaxation rate `Γ` for the linear pulse
    completeSet_DonorDot(
        x = linspace(-3.0, 4.0, 25).map { 10.0.pow(it) / _ns },
        variable = "Γ",
        samples = 20,

        tf = 18.0 * _ns,
        a = TAU * 100.0 * 1e6 / _s,
        Ω = 20.0 * _μeV / _ħ,

        σ = 1.0 * _μeV / _ħ,
        τ_c = 1.0 * _ns,

        τ = 4.0 * _ns,
        ε_max = 2000 * _μeV / _ħ,
        ε_min = -200 * _μeV / _ħ,

        useShapedPulse = true,
        useSmoothPulse = true,

        saveData = true,
        saveName = "2021 09 28 DonorDot",
    )

    // transfer error depending on relaxation rate `Γ` for the fast-QUAD pulse
    completeSet_DonorDot(
        x = linspace(-3.0, 4.0, 25).map { 10.0.pow(it) / _ns },
        variable = "Γ",
        samples = 20,

        tf = 18.0 * _ns,
        a = TAU * 100.0 * 1e6 / _s,
        Ω = 20.0 * _μeV / _ħ,

        σ = 1.0 * _μeV / _ħ,
        τ_c = 1.0 * _ns,

        τ = 4.0 * _ns,
        ε_max = 2000 * _μeV / _ħ,
        ε_min = -200 * _μeV / _ħ,

        useShapedPulse = false,
        useSmoothPulse = true,

        saveData = true,
        saveName = "2021 09 28 DonorDot",
    )

    // transfer error depending on noise correlation time `τ_c = 1/γ` for the fast-QUAD pulse
    completeSet_DonorDot(
        x = linspace(-3.0, 1.0, 25).map { 10.0.pow(it) * _ns },
        variable = "τ_c",
        samples = 20,

        tf = 18.0 * _ns,
        a = TAU * 100.0 * 1e6 / _s,
        Ω = 20.0 * _μeV / _ħ,

        σ = 1.0 * _μeV / _ħ,

        τ = 4.0 * _ns,
        ε_max = 2000 * _μeV / _ħ,
        ε_min = -200 * _μeV / _ħ,

        useShapedPulse = true,
        useSmoothPulse = true,

        saveData = true,
        saveName = "2021 09 28 DonorDot",
    )

    // transfer error depending on noise correlation time `τ_c = 1/γ` for the linear pulse
    completeSet_DonorDot(
        x = linspace(-3.0, 1.0, 25).map { 10.0.pow(it) * _ns },
        variable = "τ_c",
        samples = 20,

        tf = 18.0,
        a = TAU * 100.0 * 1e6 / _s,
        Ω = 20.0 * _μeV / _ħ,

        σ = 1.0 * _μeV / _ħ,

        τ = 4.0 * _ns,
        ε_max = 2000 * _μeV / _ħ,
        ε_min = -200 * _μeV / _ħ,

        useShapedPulse = false,
        useSmoothPulse = true,

        saveData = true,
        saveName = "2021 09 28 DonorDot",
    )


    // transfer error depending on `tf` for the following list of values
    // for the time constant `τ` of the pulse-smoothing Gaussian.
    val time_constant_τ = concatenate(
        listOf(0.1, 0.25, 0.5, 0.75),
        linsteps(1.0, 1.0, 10.0).map { it * _ns }
    )

    val dateName = "2021 09 24v2 DonorDot"

    time_constant_τ.forEach { smooth ->
        println("smooth = $smooth")

        completeSet_DonorDot(
            x = concatenate(
                linsteps(5.0, 1.0, 40.0).map { it * _ns }
            ),
            variable = "tf",
            samples = 5,

            τ = smooth,

            ε_max = 2000 * _μeV / _ħ,
            ε_min = -200 * _μeV / _ħ,

            a = TAU * 100.0 * 1e6 / _s,
            Ω = 20.0 * _μeV / _ħ,

            σ = 1.0 * _μeV / _ħ,
            τ_c = 1.0 * _ns,

            Γ = 0.0 / _ns,

            useShapedPulse = true, // fast-QUAD pulse
            useSmoothPulse = true,

            saveData = true,
            saveName = dateName,
        )

        completeSet_DonorDot(
            x = concatenate(
                linsteps(5.0, 1.0, 40.0).map { it * _ns }
            ),
            variable = "tf",
            samples = 5,

            τ = smooth,

            ε_max = 2000 * _μeV / _ħ,
            ε_min = -200 * _μeV / _ħ,
            a = TAU * 100.0 * 1e6 / _s,
            Ω = 20.0 * _μeV / _ħ,

            σ = 1.0 * _μeV / _ħ,
            τ_c = 1.0 * _ns,

            Γ = 0.0 / _ns,

            useShapedPulse = false, // linear pulse
            useSmoothPulse = true,

            saveData = true,
            saveName = dateName,
        )
    }
}


/**
 * Calculate all the data for the donor-dot model.
 */
fun double_quantum_dot() {

    val setup = DqdSetup(
        δbz = 0.1 * _μeV / _ħ // very close to the donor_dot value, rounded to next full μeV
    )
    val saveName = "2022 03 10 DQD"


//    // transfer error depending on `tf` for the following list of values
//    // for the time constant `τ` of the pulse-smoothing Gaussian.
//    // ! This should be the FIRST run of calculation you do, before settling on a time constant τ
//    val time_constant_τ = concatenate(
//        listOf(0.1, 0.25, 0.5, 0.75),
//        linsteps(1.0, 1.0, 10.0).map { it * _ns }
////         linsteps(3.0, 1.0, 10.0).map { it * _ns }
//    )


//    val x = linsteps(0.7, 0.025, 1.7).map { 10.0.pow(it) }
//    val x = concatenate(
////                linsteps(0.25, 0.25, 5.0).map { it * _ns }
////        listOf(0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5),
//        linsteps(6.0, 1.0, 20.0).map { it * _ns },
//        linsteps(20.0, 5.0, 60.0).map { it * _ns },
//    )

//    time_constant_τ.forEach { smooth ->
//        println("smooth = $smooth")
//
//        completeSet_DQD(
//            x = x,
//            variable = "tf",
//            setup = setup.copy(τ = smooth),
//
//            useShapedPulse = true, // fast-QUAD pulse
//            useSmoothPulse = true,
//
//            saveData = true,
//            saveName = saveName,
//        )
//
//        completeSet_DQD(
//            x = x,
//            variable = "tf",
//            setup = setup.copy(τ = smooth),
//
//            useShapedPulse = false, // linear pulse
//            useSmoothPulse = true,
//
//            saveData = true,
//            saveName = saveName,
//        )
//    }

//
//    // ! This run is to analyze the effect of the effective magnetic field gradient δbz on the error
//    // ? units _μeV/_ħ will be multiplied down below
////    val δbz_range = linsteps(-2.0,0.25,1.0).map{10.0.pow(it)} // * original datapoints
//    val δbz_range_calculated = linsteps(-2.0,0.25,1.0)
//    val δbz_range = linsteps(-2.0,0.125,-0.5).filterNot{ δbz_range_calculated.asList().contains(it) }.map{10.0.pow(it)} // * additional datapoints
//    // ? relevant values are the first 6 for δbz, so 10^{-2.0,-1.75,-1.5,-1.25,-1,-0.75,-0.5}
//    // ? and tf is most relevant in [10, 50] ns
//
//    δbz_range.forEach { δbz ->
//        println("δbz = $δbz")
//        val saveName = "${saveName} δbz%.2e".format(δbz)
////        val tf = linsteps(0.75, 0.025, 1.75).map { 10.0.pow(it) * _ns }
//        val tf = linsteps(1.0, 0.025, 1.75).map { 10.0.pow(it) * _ns }
//        val setup = setup.copy(δbz=δbz*_μeV/_ħ, samples=5) // ! don't forget the units
//
//        completeSet_DQD(
//            x = tf,
//            variable = "tf",
//            setup = setup,
//
//            useShapedPulse = true, // fast-QUAD pulse
//            useSmoothPulse = true,
//
//            saveData = true,
//            saveName = saveName,
//        )
//
//        completeSet_DQD(
//            x = tf,
//            variable = "tf",
//            setup = setup,
//
//            useShapedPulse = false, // linear pulse
//            useSmoothPulse = true,
//
//            saveData = true,
//            saveName = saveName,
//        )
//    }


    // ! the following block should be evaluated AFTER the optimal time constant τ has been determined
    // transfer error depending on `tf` for the linear pulse
//    completeSet_DQD(
//        x = linspace(0.0, 2.0, 100).map { 10.0.pow(it) / _ns },
//        variable = "tf",
//        setup = setup,
//
//        useShapedPulse = false,
//        useSmoothPulse = true,
//
//        saveData = true,
//        saveName = saveName,
//    )
//
//    // transfer error depending on `tf` for the fast-QUAD pulse
//    completeSet_DQD(
//        x = linspace(0.0, 2.0, 100).map { 10.0.pow(it) / _ns },
//        variable = "tf",
//        setup = setup,
//
//        useShapedPulse = true,
//        useSmoothPulse = true,
//
//        saveData = true,
//        saveName = saveName,
//    )

    // transfer error depending on relaxation rate `Γ` for the linear pulse
    completeSet_DQD(
        x = linspace(-3.0, 4.0, 25).map { 10.0.pow(it) / _ns },
        variable = "Γ",
        setup = setup,

        useShapedPulse = true,
        useSmoothPulse = true,

        saveData = true,
        saveName = saveName,
    )

    // transfer error depending on relaxation rate `Γ` for the fast-QUAD pulse
    completeSet_DQD(
        x = linspace(-3.0, 4.0, 25).map { 10.0.pow(it) / _ns },
        variable = "Γ",
        setup = setup,

        useShapedPulse = false,
        useSmoothPulse = true,

        saveData = true,
        saveName = saveName,
    )

    // transfer error depending on noise correlation time `τ_c = 1/γ` for the fast-QUAD pulse
    completeSet_DQD(
        x = linspace(-3.0, 1.0, 25).map { 10.0.pow(it) * _ns },
        variable = "τ_c",
        setup = setup,

        useShapedPulse = true,
        useSmoothPulse = true,

        saveData = true,
        saveName = saveName,
    )

    // transfer error depending on noise correlation time `τ_c = 1/γ` for the linear pulse
    completeSet_DQD(
        x = linspace(-3.0, 1.0, 25).map { 10.0.pow(it) * _ns },
        variable = "τ_c",
        setup = setup,

        useShapedPulse = false,
        useSmoothPulse = true,

        saveData = true,
        saveName = saveName,
    )
}


/**
 * This class simply stores the default values of the double-quantum dot system parameters.
 */
data class DqdSetup(
    val samples: Int = 20,

    // system parameters
    val tf: Double = 18.0 * _ns,
    val Ω: Double = 20.0 * _μeV / _ħ,
    val δbz: Double = 2.0 * _μeV / _ħ,
    val ε_max: Double = 2000.0 * _μeV / _ħ,
    val ε_min: Double = -200.0 * _μeV / _ħ,
    val τ: Double = 4.0 * _ns,

    //  detuning noise
    val σ: Double = 1.0 * _μeV / _ħ,
    val τ_c: Double = 1.0 * _ns,

    // relaxation
    val Γ: Double = 0.0 / _ns,
)


/**
 * This function generates a list of parameters for which the population transfer is being calculated.
 * `x` is the variable over the list is generated over. By setting `variable`, you can choose what parameter to
 * iterate over. Available are `"Γ"`, `"τ_c"`, `"noiseVariance"`, `"tf"`, `"smooth"`. The default is `"tf"`.
 * It will solve the master equation for a number of `samples` independent (noise) realizations.
 */
fun completeSet_DQD(
    x: List<Double> = linspace(-1.0, 3.0, 100).map { 10.0.pow(it) },
    variable: String = "tf",

    setup: DqdSetup,

    useSmoothPulse: Boolean = true,
    useShapedPulse: Boolean = false,
    saveData: Boolean = true,
    saveName: String = "2021 12 03 DQD",
    plotData: Boolean = true,
) {

    val (samples, tf, Ω, δbz, ε_max, ε_min, τ, σ, τ_c, Γ) = setup // deconstruct the variables

    val γ = 1.0 / τ_c
    runBlocking(Executors.newSingleThreadExecutor().asCoroutineDispatcher()) { // ? single threaded
//    runBlocking(Executors.newFixedThreadPool(8).asCoroutineDispatcher()) { // ? 8 threads
//    runBlocking(Dispatchers.Default) { // ? system default (as many as possible usually)
        // ? only interested in initializing the lowest two eigenstates, not the excited hybridized singlet state
        List(2) { initial ->
            x.map {
                val cutoff = max(2.0 * π * ε_max, γ) * 10.0
                val initialSpacing = 2 * π / (cutoff * 2)

                when (variable) {
                    "Γ" -> DoubleQuantumDotModel(
                        initial,
                        tf,
                        initialSpacing,
                        δbz,
                        Ω,
                        useShapedPulse,
                        useSmoothPulse,
                        τ,
                        ε_max,
                        ε_min,
                        σ,
                        τ_c,
                        it,
                    )
                    "τ_c" -> DoubleQuantumDotModel(
                        initial,
                        tf,
                        initialSpacing,
                        δbz,
                        Ω,
                        useShapedPulse,
                        useSmoothPulse,
                        τ,
                        ε_max,
                        ε_min,
                        σ,
                        it,
                        Γ,
                    )
                    "noiseVariance" -> DoubleQuantumDotModel(
                        initial,
                        tf,
                        initialSpacing,
                        δbz,
                        Ω,
                        useShapedPulse,
                        useSmoothPulse,
                        τ,
                        ε_max,
                        ε_min,
                        it,
                        τ_c,
                        Γ,
                    )
                    "tf" -> DoubleQuantumDotModel(
                        initial,
                        it,
                        initialSpacing,
                        δbz,
                        Ω,
                        useShapedPulse,
                        useSmoothPulse,
                        τ,
                        ε_max,
                        ε_min,
                        σ,
                        τ_c,
                        Γ,
                    )
                    "smooth" -> DoubleQuantumDotModel(
                        initial,
                        tf,
                        initialSpacing,
                        δbz,
                        Ω,
                        useShapedPulse,
                        useSmoothPulse,
                        it,
                        ε_max,
                        ε_min,
                        σ,
                        τ_c,
                        Γ,
                    )
                    "δbz" -> DoubleQuantumDotModel(
                        initial,
                        tf,
                        initialSpacing,
                        it,
                        Ω,
                        useShapedPulse,
                        useSmoothPulse,
                        τ,
                        ε_max,
                        ε_min,
                        σ,
                        τ_c,
                        Γ,
                    )
                    else -> DoubleQuantumDotModel(
                        initial,
                        it,
                        initialSpacing,
                        δbz,
                        Ω,
                        useShapedPulse,
                        useSmoothPulse,
                        τ,
                        ε_max,
                        ε_min,
                        σ,
                        τ_c,
                        Γ,
                    )
                }
            }
        }.mapIndexed { initial, models ->
            async {
                val res = sampleSweeps(
                    models,
                    samples,
                )
                if (saveData) {
                    val smooth = if (useSmoothPulse) "smooth(%.2e) ".format(τ) else ""
                    val shaped = if (useShapedPulse) "shaped " else ""
                    val filename =
                        "_results_/$saveName $variable ${smooth}${shaped}i$initial Ω%.2e s%.2e tc%.2e coll%.2e tf%.2e n$samples.csv"
                            .format(Ω, σ, τ_c, Γ, tf)
                    saveSweep(filename, x, res)
                }
                if (plotData)
                    try {
                        plotSweep(x, res)
                    } catch (e: Exception) {
                        println("could not plot")
                    }
            }
        }.awaitAll()
    }
}


/**
 * This function generates a list of parameters for which the population transfer is being calculated.
 * `x` is the variable over the list is generated over. By setting `variable`, you can choose what parameter to
 * iterate over. Available are `"Γ"`, `"τ_c"`, `"noiseVariance"`, `"tf"`, `"smooth"`. The default is `"tf"`.
 * It will solve the master equation for a number of `samples` independent (noise) realizations.
 */
fun completeSet_DonorDot(
    samples: Int = 20,
    x: List<Double> = linspace(-1.0, 3.0, 100).map { 10.0.pow(it) },
    variable: String = "tf",

    tf: Double = 100.0,
    a: Double = TAU * 100.0 * 1e6 / _s,
    Ω: Double = 20.0 * _μeV / _ħ,

    τ: Double = 5.0 * _ns,
    ε_max: Double = 1700 * _μeV / _ħ,
    ε_min: Double = -200 * _μeV / _ħ,

    σ: Double = 1.0 * _μeV,
    τ_c: Double = 1.0 * _ns,
    Γ: Double = 0.0,

    useSmoothPulse: Boolean = true,
    useShapedPulse: Boolean = false,
    saveData: Boolean = true,
    saveName: String = "2021 02 15 DonorDot",
    plotData: Boolean = true,
) {
    val γ = 1.0 / τ_c
    runBlocking(Executors.newSingleThreadExecutor().asCoroutineDispatcher()) { // ? single threaded
//    runBlocking(Executors.newFixedThreadPool(8).asCoroutineDispatcher()) { // ? 8 threads
//    runBlocking(Dispatchers.Default) { // ? system default (as many as possible usually)
        // ? only interested in initializing the lowest two eigenstates, not the excited hybridized singlet state
        List(2) { initial ->
            x.map {
                val cutoff = max(2.0 * π * ε_max, γ) * 10.0
                val initialSpacing = 2 * π / (cutoff * 2)

                when (variable) {
                    "Γ" -> DonorDotModel(
                        initial,
                        tf,
                        initialSpacing,
                        a,
                        Ω,
                        useShapedPulse,
                        useSmoothPulse,
                        τ,
                        ε_max,
                        ε_min,
                        σ,
                        τ_c,
                        it,
                    )
                    "τ_c" -> DonorDotModel(
                        initial,
                        tf,
                        initialSpacing,
                        a,
                        Ω,
                        useShapedPulse,
                        useSmoothPulse,
                        τ,
                        ε_max,
                        ε_min,
                        σ,
                        it,
                        Γ,
                    )
                    "noiseVariance" -> DonorDotModel(
                        initial,
                        tf,
                        initialSpacing,
                        a,
                        Ω,
                        useShapedPulse,
                        useSmoothPulse,
                        τ,
                        ε_max,
                        ε_min,
                        it,
                        τ_c,
                        Γ,
                    )
                    "tf" -> DonorDotModel(
                        initial,
                        it,
                        initialSpacing,
                        a,
                        Ω,
                        useShapedPulse,
                        useSmoothPulse,
                        τ,
                        ε_max,
                        ε_min,
                        σ,
                        τ_c,
                        Γ,
                    )
                    "smooth" -> DonorDotModel(
                        initial,
                        tf,
                        initialSpacing,
                        a,
                        Ω,
                        useShapedPulse,
                        useSmoothPulse,
                        it,
                        ε_max,
                        ε_min,
                        σ,
                        τ_c,
                        Γ,
                    )
                    else -> DonorDotModel(
                        initial,
                        it,
                        initialSpacing,
                        a,
                        Ω,
                        useShapedPulse,
                        useSmoothPulse,
                        τ,
                        ε_max,
                        ε_min,
                        σ,
                        τ_c,
                        Γ,
                    )
                }
            }
        }.mapIndexed { initial, models ->
            async {
                val res = sampleSweeps(
                    models,
                    samples,
                )
                if (saveData) {
                    val smooth = if (useSmoothPulse) "smooth(%.2e) ".format(τ) else ""
                    val shaped = if (useShapedPulse) "shaped " else ""
                    val filename =
                        "_results_/$saveName $variable ${smooth}${shaped}i$initial Ω%.2e s%.2e tc%.2e coll%.2e tf%.2e n$samples.csv"
                            .format(Ω, σ, τ_c, Γ, tf)
                    saveSweep(filename, x, res)
                }
                if (plotData)
                    try {
                        plotSweep(x, res)
                    } catch (e: Exception) {
                        println("could not plot")
                    }
            }
        }.awaitAll()
    }
}


/**
 * This function generates a list of parameters for which the population transfer is being calculated.
 * `x` is the variable over the list is generated over. By setting `variable`, you can choose what parameter to
 * iterate over. Available are `"σ"`, `"γ"`, (implicitly) `"tf"`. The default is `"tf"`.
 * It will solve the master equation for a number of `samples` independent (noise) realizations.
 */
fun completeSet_ConstantGap(
    samples: Int = 20,
    x: List<Double> = linspace(-1.0, 3.0, 100).map { 10.0.pow(it) },
    variable: String = "tf",

    gap: Double = 10.0,
    γ: Double = 1.0,
    σ: Double = 0.1,
    tf: Double = 100.0,

    superAdiabatic: Boolean = false,

    saveData: Boolean = true,
    saveName: String = "2021 02 01",
    plotData: Boolean = true,
) {

    // The preparation and measurement unitaries for the generalized strategy
    fun transformations(tf: Double): Pair<Operator?, Operator?> {
        if (!superAdiabatic) return Pair(null, null);

        val ϕ = atan(π / (tf * gap))
        println("ϕ(tf = $tf) = $ϕ")

        val initTrans = rotPauliX(ϕ);
        val finalTrans = (rotPauliY(π) * rotPauliX(ϕ).dagger() * rotPauliY(π).dagger());

        return Pair(initTrans, finalTrans)
    }




    runBlocking(Dispatchers.Default) { // ? multi-threaded
//    runBlocking(Executors.newSingleThreadExecutor().asCoroutineDispatcher()) { // ? single threaded

        List(2) { initial ->
            x.map {
                val cutoff = 2.0 * π * gap * 10 //max(2.0 * π * gap * 10, 10 * γ)
                val maxIntegrationStep = 2 * π / (cutoff * 10)

                when (variable) {
                    "σ" -> {
                        val trans = transformations(tf)
                        ConstantGapModel(
                            initial,
                            tf,
                            maxIntegrationStep,
                            gap,
                            it,
                            γ,
                            trans.first,
                            trans.second,
                        )
                    }
                    "γ" -> {
                        val trans = transformations(tf)
                        ConstantGapModel(
                            initial,
                            tf,
                            maxIntegrationStep,
                            gap,
                            σ,
                            it,
                            trans.first,
                            trans.second,
                        )
                    }
                    else -> {
                        val trans = transformations(it);
                        ConstantGapModel(
                            initial,
                            it,
                            maxIntegrationStep,
                            gap,
                            σ,
                            γ,
                            trans.first,
                            trans.second,
                        )
                    }
                }
            }
        }.mapIndexed { initial, models ->
            async {
                val res = sampleSweeps(
                    models,
                    samples,
                )
                if (saveData) {
                    val filename =
                        "_results_/$saveName CG $variable i$initial B%.2e s%.2e g%.2e tf%.2e n$samples.csv".format(
                            gap,
                            σ,
                            γ,
                            tf
                        )
                    saveSweep(filename, x, res)
                }
                if (plotData)
                    try {
                        plotSweep(x, res)
                    } catch (e: Exception) {
                        println("could not plot")
                    }
            }
        }.awaitAll()
    }
}


/**
 * This function generates a list of parameters for which the population transfer is being calculated.
 * `x` is the variable over the list is generated over. By setting `variable`, you can choose what parameter to
 * iterate over. Available are `"σ"`, `"γ"`, (implicitly) `"tf"`. The default is `"tf"`.
 * It will solve the master equation for a number of `samples` independent (noise) realizations.
 */
fun completeSet_LandauZener(
    samples: Int = 20,
    x: List<Double> = linspace(-4.0, 3.0, 200).map { 10.0.pow(it) },
    variable: String = "tf",

    gap: Double = 1.0,
    γ: Double = 1.0,
    σ: Double = 0.1,
    ε_max: Double = 10.0 * gap,
    tf: Double = 100.0,
    saveData: Boolean = true,
    saveName: String = "2021 02 01 LZ",
    useShapedPulse: Boolean = false,

    plotData: Boolean = true,
) {
    val name = if (!useShapedPulse) saveName else "$saveName shaped"

//    runBlocking(Dispatchers.Default) {
    runBlocking(Executors.newSingleThreadExecutor().asCoroutineDispatcher()) {

        List(2) { initial ->
            x.map {
                val cutoff = max(2.0 * π * ε_max * 10, 10 * γ)
                val initialSpacing = 2 * π / (cutoff * 10)
//                if(useShapedPulse) initialSpacing *= 0.1

                when (variable) {
                    "σ" -> LandauZenerModel(
                        initial,
                        tf,
                        initialSpacing,
                        gap,
                        it,
                        γ,
                        ε_max,
                        useShapedPulse,
                    )
                    "γ" -> LandauZenerModel(
                        initial,
                        tf,
                        initialSpacing,
                        gap,
                        σ,
                        it,
                        ε_max,
                        useShapedPulse,
                    )
                    else -> LandauZenerModel(
                        initial,
                        it,
                        initialSpacing,
                        gap,
                        σ,
                        γ,
                        ε_max,
                        useShapedPulse,
                    )
                }
            }
        }.mapIndexed { initial, models ->
            async {
                val res = sampleSweeps(
                    models,
                    samples,
                )
                if (saveData) {
                    val filename = "_results_/$name $variable i$initial B%.2e s%.2e g%.2e tf%.2e n$samples.csv".format(
                        gap,
                        σ,
                        γ,
                        tf
                    )
                    saveSweep(filename, x, res)
                }
                if (plotData)
                    try {
                        plotSweep(x, res)
                    } catch (e: Exception) {
                        println("could not plot")
                    }
            }
        }.awaitAll()
    }
}