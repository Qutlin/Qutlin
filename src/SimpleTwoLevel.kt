import Qutlin.*
import Qutlin.ComplexMatrix.Companion.rotPauliX
import Qutlin.ComplexMatrix.Companion.rotPauliY
import kotlinx.coroutines.*
import org.hipparchus.util.FastMath.*
import java.util.concurrent.Executors
import kotlin.math.pow


fun main() {
//    completeSet_ConstantGap(
//        x = linspace(-4.0, 1.0, 100).map { 10.0.pow(it) },
////        x = linspace(0.01, 10.0, 50).map{ it },
//        samples = 20,
//        gap = 10.0,
//        γ = 1.0,
//        σ = 0.1,
//        tf = 100.0,
//        variable = "σ",
//        saveName = "2021 05 04",
//        superAdiabatic = false,
//    )

    // completeSet_ConstantGap(
    //     x = linspace(-2.0, 3.0, 200).map { 10.0.pow(it) },
    //     samples = 20,
    //     gap = 10.0,
    //     γ = 1.0,
    //     σ = 0.1,
    //     tf = 100.0,
    //     variable = "tf",
    //     saveName = "2021 05 31",
    //     superAdiabatic = false,
    // )

//    completeSet_ConstantGap(
//        x = linspace(-3.0, 4.0, 50).map{10.0.pow(it)},
//        samples = 20,
//        variable = "γ",
//        saveName = "2021 02 09"
//    )

//    completeSet_LandauZener(
//        x = concatenate(
//            linspace(0.0, 1.0, 50).map { 10.0.pow(it) },
//            linspace(1.0, 3.0, 50, skipFirst = true).map { 10.0.pow(it) },
//        ),
//        samples=20,
//        σ = 0.1,
//        γ = 1.0,
//        gap = 1.0,
//        variable = "tf",
//        useShapedPulse = false,
//        saveName = "2021 05 04 LZ",
//    )

//    completeSet_LandauZener(
//        x = linspace(-3.0, 1.0, 100).map { 10.0.pow(it) },
//        samples = 20,
//        variable = "σ",
//        saveName = "2021 05 04 LZ",
//    )
//
//    completeSet_LandauZener(
//        x = linspace(-4.0, 4.0, 100).map { 10.0.pow(it) },
//        samples = 20,
//        variable = "γ",
//        saveName = "2021 05 04 LZ",
//    )


//    val dateName = "2021 05 06 DonorDot"

    completeSet_DonorDot(
        x = linspace(0.0, 2.0, 100).map { 10.0.pow(it)/_ns },
        variable = "tf",
        samples = 20,

        tf = 16.6, // * does not matter here, will be overwritten
        a = TAU * 100.0 *1e6 / _s,
        Ω = 20.0 * _μeV / _ħ,

        σ = 1.0 * _μeV / _ħ, // TODO correct values!?
        τ_c = 1.0 * _ns,

        useShapedPulse = false,

        useSmoothPulse = true,
        saveData = true,
        saveName = "2021 08 24 DonorDot",
    )

    completeSet_DonorDot(
        x = linspace(0.0, 2.0, 100).map { 10.0.pow(it)/_ns },
        variable = "tf",
        samples = 20,

        tf = 16.6, // * does not matter here, will be overwritten
        a = TAU * 100.0 *1e6 / _s,
        Ω = 20.0 * _μeV / _ħ,

        σ = 1.0 * _μeV / _ħ, // TODO correct values!?
        τ_c = 1.0 * _ns,

        useShapedPulse = true,

        useSmoothPulse = true,
        saveData = true,
        saveName = "2021 08 24 DonorDot",
    )


//    completeSet_DonorDot(
//        x = linspace(-3.0, 4.0, 25).map { 10.0.pow(it)/_ns },
//        variable = "Γ",
//        samples = 20,
//
//        tf = 16.6,
//        a = TAU * 100.0 * 1e6 / _s,
//        Ω = 20.0 * _μeV / _ħ,
//
//        σ = 1.0 * _μeV / _ħ, // TODO correct values!?
//        τ_c = 1.0 * _ns,
//
//        useShapedPulse = true,
//
//        useSmoothPulse = true,
//        saveData = true,
//        saveName = "2021 06 09 DonorDot",
//    )
//
//    completeSet_DonorDot(
//        x = linspace(-3.0, 4.0, 25).map { 10.0.pow(it)/_ns },
//        variable = "Γ",
//        samples = 20,
//
//        tf = 16.6,
//        a = TAU * 100.0 * 1e6 / _s,
//        Ω = 20.0 * _μeV / _ħ,
//
//        σ = 1.0 * _μeV / _ħ, // TODO correct values!?
//        τ_c = 1.0 * _ns,
//
//        useShapedPulse = false,
//
//        useSmoothPulse = true,
//        saveData = true,
//        saveName = "2021 06 09 DonorDot",
//    )



//   completeSet_DonorDot(
//       x = linspace(-3.0, 1.0, 100).map { 10.0.pow(it) * _ns },
//       variable = "τ_c",
//       samples = 20,
//
//       tf = 16.6,
//       a = TAU * 100.0 *1e6 / _s,
//       Ω = 20.0 * _μeV / _ħ,
//
//       σ = 1.0 * _μeV / _ħ, // TODO correct values!?
////        τ_c = 1.0 * _ns,
//
//       useShapedPulse = true,
//
//       useSmoothPulse = true,
//       saveData = true,
//       saveName = "2021 06 03 DonorDot",
//   )
//    completeSet_DonorDot(
//        x = linspace(-3.0, 1.0, 100).map { 10.0.pow(it) * _ns },
//        variable = "τ_c",
//        samples = 20,
//
//        tf = 16.6,
//        a = TAU * 100.0 *1e6 / _s,
//        Ω = 20.0 * _μeV / _ħ,
//
//        σ = 1.0 * _μeV / _ħ, // TODO correct values!?
// //        τ_c = 1.0 * _ns,
//
//        useShapedPulse = false,
//
//        useSmoothPulse = true,
//        saveData = true,
//        saveName = "2021 06 03 DonorDot",
//    )



//
//    val smooth = concatenate(
//        listOf(0.1, 0.25, 0.5, 0.75),
//        linsteps(1.0, 1.0, 20.0).map { it * _ns }
//    )
//
//    smooth.forEach { smooth ->
//        println("smooth = $smooth")
//
//        completeSet_DonorDot(
//            x = concatenate(
//                linsteps(5.0, 1.0, 40.0).map { it * _ns }
//            ),
//            variable = "tf",
//            samples = 5,
//
//            τ = smooth,
//
//            tf = 16.6,
//            a = TAU * 100.0 * 1e6 / _s,
//            Ω = 20.0 * _μeV / _ħ,
//
//            noiseVariance = 1.0 * _μeV / _ħ, // TODO correct values!?
//            τ_c = 1.0 * _ns,
//
//            Γ = 0.0 / _ns,
//
//            useShapedPulse = true,
//
//            useSmoothPulse = true,
//            saveData = true,
//            saveName = dateName,
//        )
//
//        completeSet_DonorDot(
//            x = concatenate(
//                linsteps(5.0, 1.0, 40.0).map { it * _ns }
//            ),
//            variable = "tf",
//            samples = 5,
//
//            τ = smooth,
//
//            tf = 16.6,
//            a = TAU * 100.0 * 1e6 / _s,
//            Ω = 20.0 * _μeV / _ħ,
//
//            noiseVariance = 1.0 * _μeV / _ħ, // TODO correct values!?
//            τ_c = 1.0 * _ns,
//
//            Γ = 0.0 / _ns,
//
//            useShapedPulse = false,
//
//            useSmoothPulse = true,
//            saveData = true,
//            saveName = dateName,
//        )
//    }

}







fun completeSet_DonorDot(
        samples: Int = 20,
        x: List<Double> = linspace(-1.0, 3.0, 100).map { 10.0.pow(it) },
        variable: String = "tf",

        tf: Double = 100.0,
        a: Double = TAU * 100.0 *1e6 / _s,
        Ω: Double = 20.0 * _μeV / _ħ,
        
        useShapedPulse: Boolean = false,
        useSmoothPulse: Boolean = true,
        τ: Double = 5.0 * _ns,

        ε_max: Double = 1700 * _μeV / _ħ,
        ε_min: Double = -200 * _μeV / _ħ,
        σ: Double = 1.0 * _μeV, // TODO correct values!?
        τ_c: Double = 1.0 * _ns,
        Γ: Double = 0.0,

        saveData: Boolean = true,
        saveName: String = "2021 02 15 DonorDot",
        plotData: Boolean = true,
) {
    val γ = 1.0/τ_c
    runBlocking(Executors.newSingleThreadExecutor().asCoroutineDispatcher()) {
//    runBlocking(Executors.newFixedThreadPool(8).asCoroutineDispatcher()) {
//    runBlocking(Dispatchers.Default) {
        // ? only interested in initializing the lowest two eigenstates, not the excited hybridized singlet state
        List(2) { initial ->
            x.map {
                val cutoff = max(2.0 * π * ε_max, γ) * 10.0
//                val initialSpacing = 2*π/(cutoff * 10)
                val initialSpacing = 2*π/(cutoff * 2)

                when(variable) {
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
                if(saveData) {
                    val smooth = if(useSmoothPulse) "smooth(%.2e) ".format(τ) else ""
                    val shaped = if(useShapedPulse) "shaped " else ""
                    val filename = "results/$saveName $variable ${smooth}${shaped}i$initial Ω%.2e s%.2e tc%.2e coll%.2e tf%.2e n$samples.csv"
                        .format(Ω, σ, τ_c, Γ, tf)
                    saveSweep(filename, x, res)
                }
                if(plotData)
                    try { plotSweep(x, res) }
                    catch (e: Exception) { println("could not plot") }
            }
        }.awaitAll()
    }
}











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


    fun transformations(tf: Double) : Pair<Operator?, Operator?> {
        if(!superAdiabatic) return Pair(null, null);

        val ϕ = atan(π / (tf * gap))
        println("ϕ(tf = $tf) = $ϕ")

//        val initTrans = ComplexMatrix(Pair(2, 2), complexArrayOf(
//            cos(ϕ/2.0), -I*sin(ϕ/2.0),
//            -I*sin(ϕ/2.0), cos(ϕ/2.0)
//        ))
        val initTrans = rotPauliX(ϕ);
        val finalTrans = (rotPauliY(π) * rotPauliX(ϕ).dagger() * rotPauliY(π).dagger());

        // the final transformation acts on the measurement states, not the time-evolved states!
        // meaning we measure Tr[ (Um^+ |1(tf)><1(tf)| Um) rho(t_f) ]
        // and we have Up = Um^+, so initTrans == finalTrans in this specific case (CG model, theta_f = π)
//        val finalTrans = initTrans//.dagger()

        return Pair(initTrans, finalTrans)
    }




    runBlocking(Dispatchers.Default) {
//    runBlocking(Executors.newSingleThreadExecutor().asCoroutineDispatcher()) {

        List(2) { initial ->
            x.map {
                val cutoff = 2.0*π*gap*10 //max(2.0 * π * gap * 10, 10 * γ)
                val maxIntegrationStep = 2*π/(cutoff * 10)

                when(variable) {
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
                if(saveData) {
                    val filename = "results/$saveName CG $variable i$initial B%.2e s%.2e g%.2e tf%.2e n$samples.csv".format(gap, σ, γ, tf)
                    saveSweep(filename, x, res)
                }
                if(plotData)
                    try { plotSweep(x, res) }
                    catch (e: Exception) { println("could not plot") }
            }
        }.awaitAll()
    }
}



fun completeSet_LandauZener(
    samples: Int = 20,
    x: List<Double> = linspace(-4.0, 3.0, 200).map { 10.0.pow(it) },
    variable: String = "tf",

    gap: Double = 1.0,
    γ: Double = 1.0,
    σ: Double = 0.1,
    ε_max: Double = 10.0 * gap,
    tf: Double = 100.0,
    useShapedPulse: Boolean = false,

    saveData: Boolean = true,
    saveName: String = "2021 02 01 LZ",
    plotData: Boolean = true,
) {
    val name = if(!useShapedPulse) saveName else "$saveName shaped"

//    runBlocking(Dispatchers.Default) {
    runBlocking(Executors.newSingleThreadExecutor().asCoroutineDispatcher()) {

        List(2) { initial ->
            x.map {
                val cutoff = max(2.0 * π * ε_max * 10, 10 * γ)
                val initialSpacing = 2*π/(cutoff * 10)
//                if(useShapedPulse) initialSpacing *= 0.1

                when(variable) {
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
                if(saveData) {
                    val filename = "results/$name $variable i$initial B%.2e s%.2e g%.2e tf%.2e n$samples.csv".format(gap, σ, γ, tf)
                    saveSweep(filename, x, res)
                }
                if(plotData)
                    try { plotSweep(x, res) }
                    catch (e: Exception) { println("could not plot") }
            }
        }.awaitAll()
    }
}