import Qutlin.*
import Qutlin.ComplexMatrix.Companion.rotPauliY
import kotlinx.coroutines.*
import org.hipparchus.ode.nonstiff.AdamsMoultonIntegrator
import org.hipparchus.ode.nonstiff.GraggBulirschStoerIntegrator
import java.util.concurrent.Executors
import kotlin.math.PI
import kotlin.math.min


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
//    val t1 = t0 * 1.01;

    // * generate initial density matrix depending on the chosen initial eigenstate
    val initSys = H_0(t0).eigenSystem()
//    val initSys = H_η(t0).eigenSystem()
    var ρ_init = initSys[initialState].second.normalized().ketBra()
    if (U_p != null) {
//        println("initial transformation = $U_p")
        ρ_init = U_p * ρ_init * U_p.dagger()
    }

    // * generate the density matrices of the final eigenstates
    val finalSys = H_0(t1).eigenSystem()
    var ρ = List(initSys.size) {
        finalSys[it].second.normalized().ketBra()
    }
    if(U_m != null)
        ρ = ρ.map { U_m * it * U_m.dagger() }
//    val ρ = initSys.map { rotPauliY(PI) * it.second.normalized().ketBra() * rotPauliY(PI).dagger() }

    // * define the ODE for the von Neumann equation
    val ode = if (Γ == null)
//        fun(t: Double, ρ: ComplexMatrix) = -I * commutator(ρ, H_η(t))
        fun(t: Double, ρ: ComplexMatrix) = -I * commutator(H_η(t), ρ)
    else
//        fun(t: Double, ρ: ComplexMatrix) = -I * commutator(ρ, H_η(t)) + Γ.first * collapse(Γ.second(t), ρ)
        fun(t: Double, ρ: ComplexMatrix) =
            -I * commutator(H_η(t), ρ) + Γ.first * collapse(Γ.second(t), ρ)


    val (_, ρ_res) = solveEvolution(
        t0, t1, ρ_init, ode,
        integrator = AdamsMoultonIntegrator(
            12, Double.MIN_VALUE, min(maxIntegrationStep, tf / 10000.0),
            1e-10, 1e-8
        )
//        integrator = AdamsBashforthIntegrator(
//            (tf/maxIntegrationStep / 100).toInt(), Double.MIN_VALUE, maxIntegrationStep,
//            1e-6, 1e-4
//        )
//        integrator = GraggBulirschStoerIntegrator(
//            Double.MIN_VALUE,min(maxIntegrationStep, tf/10000.0),
//            1e-10, 1e-8
//        )
    )

//    val ρ_res = ρ_init;


    return if (U_m != null) ρ.map {
        (it * U_m * ρ_res * U_m.dagger()).trace().real
    } else ρ.map {
        (it * ρ_res).trace().real
    }
}


fun sampleSweeps(
    models: List<Model>,
    samples: Int = 20,
): List<List<List<Double>>> {

    // * prepare results storage
    val nd = models[0].dimensions
    val p = List(samples) { MutableList(nd) { MutableList(models.size) { 0.0 } } }

    println("start solver...")
//    runBlocking(Dispatchers.Default) {
//    runBlocking(Executors.newSingleThreadExecutor().asCoroutineDispatcher()) {
    runBlocking(Executors.newFixedThreadPool(14).asCoroutineDispatcher()) {
        val futures = List(samples) { sample ->
            println("sample $sample")
            models.mapIndexed { i_model: Int, model: Model ->
                async {
                    println("%3d / $samples : %3d / ${models.size}".format(sample, i_model))

                    model.build()

                    val res = integrate(
                        initialState = model.initial,
                        tf = model.tf,
                        overhang = model.overhang,
                        H_η = model.H_η,
                        H_0 = model.H_0,
                        Γ = model.collapse,
                        maxIntegrationStep = model.maxIntegrationStep,
                        U_p = model.U_p,
                        U_m = model.U_m,
                    )
                    res.forEachIndexed { i, r -> p[sample][i][i_model] = r }
                }//.await()
            }
        }
        futures.forEach { it.awaitAll() }
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