import Qutlin.*
import kotlinx.coroutines.*
import org.hipparchus.ode.nonstiff.AdamsMoultonIntegrator
import java.util.concurrent.Executors
import kotlin.math.min


/**
 * Solve the Lindblad master equation for `t ∈ [-overhang, tf + overhang]`
 *
 * `d/dt ρ(t) = -i [H_η(t), ρ(t)] + Γ_0 D[Γ[1]] ρ(t)`,
 *
 * given a noisy Hamiltonian `H_η(t)` and a collapse operator `Γ.second` with factor `Γ.first`.
 *
 * The initial density matrix `ρ(0) = ρ_init = |initialState><initialState|`.
 *
 * `U_p` and `U_m` are the preparation and measurement unitaries used for the "generalized" strategy.
 *
 *
 */
fun integrate(
    initialState: Int,
    tf: Double,
    overhang: Double = 0.0,
    H_η: (Double) -> Operator,
    H_0: (Double) -> Operator,
    Γ: Pair<Double, (Double) -> Operator>? = null,
    maxIntegrationStep: Double = 1.0,
    U_p: Operator? = null,
    U_m: Operator? = null,
): List<Double> {

    val t0 = -overhang
    val t1 = tf + overhang

    // * generate initial density matrix depending on the chosen initial eigenstate
    val initSys = H_0(t0).eigenSystem()
//    val initSys = H_η(t0).eigenSystem()

    var ρ_init = initSys[initialState].second.normalized().ketBra()
    if (U_p != null) {
        ρ_init = U_p * ρ_init * U_p.dagger()
    }

    // * generate the density matrices of the final eigenstates
    val finalSys = H_0(t1).eigenSystem()
    val projection_ops = List(initSys.size) {
        finalSys[it].second.normalized().ketBra()
    }
//    if(U_m != null)
//        ρ = ρ.map { U_m * it * U_m.dagger() }

//    val ρ = initSys.map { rotPauliY(PI) * it.second.normalized().ketBra() * rotPauliY(PI).dagger() }

    // * define the ODE for the von Neumann equation
    val ode = if (Γ == null)
        fun(t: Double, ρ: ComplexMatrix) = -I * commutator(H_η(t), ρ)
    else
        fun(t: Double, ρ: ComplexMatrix) =
            -I * commutator(H_η(t), ρ) + Γ.first * collapse(Γ.second(t), ρ)


    val (_, ρ_res) = solveEvolution(
        t0, t1, ρ_init, ode,
        integrator = AdamsMoultonIntegrator(
            12, Double.MIN_VALUE, min(maxIntegrationStep, tf / 10000.0),
            1e-10, 1e-8
        )
    )

     return if (U_m != null) projection_ops.map {Π ->
         (Π * U_m * ρ_res * U_m.dagger()).trace().real
     } else projection_ops.map {Π ->
         (Π * ρ_res).trace().real
     }
//    return ρ.map { (it * ρ_res).trace().real } // I already applied U_m on ρ above - Fehse, 2022-04-29
}

/**
 * Helper function to calculate multiple quantum systems given by the `models` parameters with `samples` number of independent noise-samples each.
 */
fun sampleSweeps(
    models: MutableList<Model>,
    samples: Int = 20,
    parallel_over_samples: Boolean = false,
): List<List<List<Double>>> {

    // * prepare results storage
    val nd = models[0].dimensions
    val p = List(samples) { MutableList(nd) { MutableList(models.size) { 0.0 } } }

    val N_THREADS = 8

    println("start solver...")
    if (parallel_over_samples) {
//    runBlocking(Executors.newSingleThreadExecutor().asCoroutineDispatcher()) { // ? single threaded execution
        runBlocking(Executors.newFixedThreadPool(N_THREADS).asCoroutineDispatcher()) {
            val futures = List(samples) { sample ->
                println("sample $sample")

                async {
                    models.mapIndexed { i_model: Int, model: Model ->
                        println("%3d / $samples : %3d / ${models.size}".format(sample+1, i_model+1))
                        model.build()

                        val res = integrate(
                            initialState = model.initial,
                            tf = model.tf,
                            overhang = model.overhang,
                            H_η = model.H_η!!,
                            H_0 = model.H_0!!,
                            Γ = model.collapse,
                            maxIntegrationStep = model.maxIntegrationStep,
                            U_p = model.U_p,
                            U_m = model.U_m,
                        )

                        res.forEachIndexed { i, r -> p[sample][i][i_model] = r }
                    }
                }
            }
            futures.awaitAll()//.forEach { it.awaitAll() }
        }
    } else {
//        runBlocking(Executors.newSingleThreadExecutor().asCoroutineDispatcher()) {
        runBlocking(Executors.newFixedThreadPool(N_THREADS).asCoroutineDispatcher()) {
            val futures = models.mapIndexed { i_model: Int, model: Model ->
                async {
                    (0 until samples).forEach { sample ->
                        println("%3d / $samples : %3d / ${models.size}".format(sample+1, i_model+1))
                        model.build()

                        val res = integrate(
                            initialState = model.initial,
                            tf = model.tf,
                            overhang = model.overhang,
                            H_η = model.H_η!!,
                            H_0 = model.H_0!!,
                            Γ = model.collapse,
                            maxIntegrationStep = model.maxIntegrationStep,
                            U_p = model.U_p,
                            U_m = model.U_m,
                        )

                        res.forEachIndexed { i, r -> p[sample][i][i_model] = r }
                    }
                }
            }
            futures.awaitAll()
        }
    }
    println("solver finished.")
    return p
}


fun plotSweep(
    x: List<Double>,
    p: List<List<List<Double>>>, // * samples, dimensions, p(x)
) {
    val nd = p[0].size
    val samples = p.size

    val datasets = List(samples * nd) {
        val dim = it % nd
        val sample = (it - dim) / nd

        Plot.DataSet(
            x.zip(p[sample][dim]) { x, y -> listOf(x, y) },
            symbol = Plot.Symbols.None,
            color = Plot.colorPalette[dim]
        )
    }

    val plot = Plot(
        datasets[0],
        xScale = Mapper.Companion.Scales.Log10, yScale = Mapper.Companion.Scales.Log10, minY = 1e-6, maxY = 1.0
    )

    datasets.takeLast(datasets.size - 1).forEach { plot.addDataSet(it) }
}

fun saveSweep(
    filename: String,
    x: List<Double>,
    p: List<List<List<Double>>>, // * samples, dimensions, p(x)
) {
    val res = mutableListOf(x)
    p.forEach { sample -> sample.forEach { dim -> res.add(dim) } }
    save(filename, res)
    println("saved data as: $filename")
}