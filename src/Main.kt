import Mapper.Companion.Scales.Linear
import Mapper.Companion.Scales.Log10
import Plot.DataSet
import Plot.Symbols
import Qutlin.*
import kotlinx.coroutines.Dispatchers
import kotlinx.coroutines.async
import kotlinx.coroutines.awaitAll
import kotlinx.coroutines.runBlocking
import org.hipparchus.complex.Complex
import org.hipparchus.util.FastMath.*

/*
888b     d888          d8b
8888b   d8888          Y8P
88888b.d88888
888Y88888P888  8888b.  888 88888b.
888 Y888P 888     "88b 888 888 "88b
888  Y8P  888 .d888888 888 888  888
888   "   888 888  888 888 888  888
888       888 "Y888888 888 888  888



*/

enum class Shape { LINEAR, SHAPED }


fun main() {
//    sweep_tau()
//    pulseShape()
//    singleSweep()

//    tauSweeps(
//        0,
//        Shape.LINEAR,
//        true,
//        1,
//        10,
//        100.0,
//        false
//    )

//    completeSet(plotData = false)

//    completeSet(   0.01, plotData = false)
//    completeSet(   0.05, plotData = false)
//    completeSet(   0.50, plotData = false)
//    completeSet(   5.00, plotData = false)
//    completeSet(  50.00, plotData = false)
//    completeSet( 500.00, plotData = false)
//    completeSet(5000.00, plotData = false)

//    completeSet(0.001, plotData = false)
//    completeSet(0.005, plotData = false)
//    completeSet(10000.00, plotData = false)

//    completeSet(0.0)


//    completeSet(collapseGamma = 0.0, smoothStd = 0.0, plotData = false)
//    completeSet(collapseGamma = 0.0, smoothStd = 0.1, plotData = false)
//    completeSet(collapseGamma = 0.0, smoothStd = 0.5, plotData = false)
//    completeSet(collapseGamma = 0.0, smoothStd = 1.0, plotData = false)
//    completeSet(collapseGamma = 0.0, smoothStd = 2.0, plotData = false)
//    completeSet(collapseGamma = 0.0, smoothStd = 3.0, plotData = false)
//    completeSet(collapseGamma = 0.0, smoothStd = 4.0, plotData = false)
//    completeSet(collapseGamma = 0.0, smoothStd = 6.0, plotData = false)
//    completeSet(collapseGamma = 0.0, smoothStd = 7.0, plotData = false)
//    completeSet(collapseGamma = 0.0, smoothStd = 8.0, plotData = false)
//    completeSet(collapseGamma = 0.0, smoothStd = 9.0, plotData = false)
//    completeSet(collapseGamma = 0.0, smoothStd = 10.0, plotData = false)

//    completeSet(collapseGamma = 0.0, correlationTime =    0.01, plotData = false)
//    completeSet(collapseGamma = 0.0, correlationTime =    0.02, plotData = false)
//    completeSet(collapseGamma = 0.0, correlationTime =    0.03, plotData = false)
//    completeSet(collapseGamma = 0.0, correlationTime =    0.04, plotData = false)
//    completeSet(collapseGamma = 0.0, correlationTime =    0.05, plotData = false)
//    completeSet(collapseGamma = 0.0, correlationTime =    0.06, plotData = false)
//    completeSet(collapseGamma = 0.0, correlationTime =    0.07, plotData = false)
//    completeSet(collapseGamma = 0.0, correlationTime =    0.08, plotData = false)
//    completeSet(collapseGamma = 0.0, correlationTime =    0.09, plotData = false)
//    completeSet(collapseGamma = 0.0, correlationTime =    0.1 , plotData = false)
    completeSet_ConstantGap_tf(collapseGamma = 0.0, correlationTime =    0.11 , plotData = false)
    completeSet_ConstantGap_tf(collapseGamma = 0.0, correlationTime =    0.12 , plotData = false)
    completeSet_ConstantGap_tf(collapseGamma = 0.0, correlationTime =    0.13 , plotData = false)
    completeSet_ConstantGap_tf(collapseGamma = 0.0, correlationTime =    0.14 , plotData = false)
    completeSet_ConstantGap_tf(collapseGamma = 0.0, correlationTime =    0.15 , plotData = false)
    completeSet_ConstantGap_tf(collapseGamma = 0.0, correlationTime =    0.16 , plotData = false)
    completeSet_ConstantGap_tf(collapseGamma = 0.0, correlationTime =    0.17 , plotData = false)
    completeSet_ConstantGap_tf(collapseGamma = 0.0, correlationTime =    0.18 , plotData = false)
    completeSet_ConstantGap_tf(collapseGamma = 0.0, correlationTime =    0.19 , plotData = false)
//    completeSet(collapseGamma = 0.0, correlationTime =    0.2 , plotData = false)
//    completeSet(collapseGamma = 0.0, correlationTime =    0.3 , plotData = false)
//    completeSet(collapseGamma = 0.0, correlationTime =    0.4 , plotData = false)
//    completeSet(collapseGamma = 0.0, correlationTime =    0.5 , plotData = false)
//    completeSet(collapseGamma = 0.0, correlationTime =    0.6 , plotData = false)
//    completeSet(collapseGamma = 0.0, correlationTime =    0.7 , plotData = false)
//    completeSet(collapseGamma = 0.0, correlationTime =    0.8 , plotData = false)
//    completeSet(collapseGamma = 0.0, correlationTime =    0.9 , plotData = false)
//    completeSet(collapseGamma = 0.0, correlationTime =    1.0 , plotData = false)
//    completeSet(collapseGamma = 0.0, correlationTime =    5.0 , plotData = false)
//    completeSet(collapseGamma = 0.0, correlationTime =   10.0 , plotData = false)
//    completeSet(collapseGamma = 0.0, correlationTime =   50.0 , plotData = false)
//    completeSet(collapseGamma = 0.0, correlationTime =  100.0 , plotData = false)
//    completeSet(collapseGamma = 0.0, correlationTime =  500.0 , plotData = false)
//    completeSet(collapseGamma = 0.0, correlationTime = 1000.0 , plotData = false)
//    completeSet(collapseGamma = 0.0, correlationTime = 5000.0 , plotData = false)
}


fun completeSet_ConstantGap_tf(
    collapseGamma: Double = 0.1 / _ns,
    smoothStd: Double = 5.0 * _ns,
    correlationTime: Double = 1.0 * _ns,
    plotData: Boolean = true
) {
    runBlocking(Dispatchers.Default) {

        val useSmooth = smoothStd > 0.0
        val samples = 10
        val tfSubdivision = 100
//        val collapseGamma = 0.1 / _ns
        val saveData = true
        val saveName = "2020 09 01"


        val gamma = 1.0 / _ns
//        val correlationTime = 1.0 * _ns

//        val noise = GaussianSpectrumNoise(
//            correlationTime,
//            variance = gamma / (TAU * correlationTime)
//        )

        val noise = FNoise(
            20.0 * _μeV / _ħ,
            TAU / (40.0 * _ns / 10), // ? given by 1/10th of slowest frequency given by tf

        )


        val futures = listOf(
            async {
                tauSweeps(
                    0,
                    Shape.SHAPED,
                    useSmooth,
                    smoothStd,
                    noise,
                    samples,
                    tfSubdivision,
                    collapseGamma,
                    saveData,
                    saveName,
                    plotData
                )
            },
            async {
                tauSweeps(
                    1,
                    Shape.SHAPED,
                    useSmooth,
                    smoothStd,
                    noise,
                    samples,
                    tfSubdivision,
                    collapseGamma,
                    saveData,
                    saveName,
                    plotData
                )
            },
            async {
                tauSweeps(
                    0,
                    Shape.LINEAR,
                    useSmooth,
                    smoothStd,
                    noise,
                    samples,
                    tfSubdivision,
                    collapseGamma,
                    saveData,
                    saveName,
                    plotData
                )
            },
            async {
                tauSweeps(
                    1,
                    Shape.LINEAR,
                    useSmooth,
                    smoothStd,
                    noise,
                    samples,
                    tfSubdivision,
                    collapseGamma,
                    saveData,
                    saveName,
                    plotData
                )
            }
        )
        futures.awaitAll()
    }
}


fun tauSweeps(
    initial: Int = 0,
    shape: Shape = Shape.LINEAR,
    useSmooth: Boolean = true,
    smoothStd: Double = 5.0 * _ns,
    noiseType: NoiseType,
    samples: Int = 20,
    tfSubdivision: Int = 500,
    collapseGamma: Double = 0.0,
    saveData: Boolean = true,
    saveName: String = "",
    plotData: Boolean = true
) {
    val a = 2.0 * PI * 100.0 * 1e6 / _s
    val deltaNm = 1.0
    val dbz = -a / 2.0 * 0.5 * deltaNm
    println("dbz = ${dbz * _ħ}")
    val om = 20.0 * _μeV / _ħ
    val epsBounds = Pair(-200 * _μeV / _ħ, 1700 * _μeV / _ħ)


    // * SETUP /////////////////////////////////////////////////////////////////////////////////////////////////////////
//    val tfMin = 0.0
//    val tfMax = 1.6
//    val tf = linspace(tfMin, tfMax, tfSubdivision).map { pow(10.0, it) }
    val tfMin = 1.0
    val tfMax = 40.0
    val tf = linspace(tfMin, tfMax, tfSubdivision).toList()


    // * pulse smoothing
//    val smoothStd = 5.0 * _ns
    val overhang = 3 * smoothStd


    val filename = run {
        var shp = when (shape) {
            Shape.LINEAR -> "lin"
            Shape.SHAPED -> "shd"
        }
        if (useSmooth) shp += " smth"

        "results/${saveName} $noiseType i${initial}_n${samples} v$smoothStd tf${tfMin}_${tfMax}_n$tfSubdivision g$collapseGamma $shp.csv"
    }


    // * linear epsilon function
    fun epsLinear(t: Double, tau: Double) = clamp(
        epsBounds.first,
        epsBounds.second + t / tau * (epsBounds.first - epsBounds.second),
        epsBounds.second
    )

    // * calculate shaped epsilon function
    @Suppress("DuplicatedCode")
    fun deps(ep: Double): Double {
        val beta = sqrt(ep * ep + om * om)
        val x = beta.e(3) * sqrt(1.0 + ep / beta)
        val y = a * a + 2 * om * om + 4 * ep * ep + ep / beta * (a * a - 4 * om * om - 4 * ep * ep)
        val z = 2 * a * om * om * (ep - 3 * beta)
        return -x * y.e(3.0 / 2.0) / z
    }

    val pulse = calculatePulseShape(epsBounds, ::deps)
//    println(pulse)
    val pulse0 = pulse.first()
    val pulse1 = pulse.last()
//    println("t0 = $pulse0, t1 = $pulse1")
//    Plot(DataSet(pulse.map { listOf(it.first, it.second) }))

    fun epsShaped(t: Double, tau: Double) = clamp(
        epsBounds.first,
        interpolate(pulse, pulse1.first * t / tau, epsBounds),
        epsBounds.second
    )


    // * define Hamiltonian of the system (3-level)
//    val em = -min(epsBounds.first, epsBounds.second) * 2.0

    //    println("$em ? $epsBounds")
    fun hamiltonian(ep: Double) = Operator(
        Pair(3, 3), complexArrayOf(
            0.0, dbz, om / 2.0,
            dbz, 0.0, 0.0,
            om / 2.0, 0.0, ep
        )
    )


    // * choose epsilon function
    fun rawEps(t: Double, tau: Double) = when (shape) {
        Shape.LINEAR -> epsLinear(clamp(0.0, t, tau), tau)// smooth(t){ x -> epsLinear(x,tau) }
        Shape.SHAPED -> epsShaped(clamp(0.0, t, tau), tau)// smooth(t){ x -> epsLinear(x,tau) }
    }


    val p0 = List(samples) { MutableList(tf.size) { 0.0 } }
    val p1 = List(samples) { MutableList(tf.size) { 0.0 } }
    val p2 = List(samples) { MutableList(tf.size) { 0.0 } }

    println("start solver...")
//    runBlocking {
    runBlocking(Dispatchers.Default) {
        println("run blocking...")
        val futures = List(samples) { sample ->
            println("sample $sample")
            tf.mapIndexed { itf, tf ->
                async {
                    println("$itf/ $tfSubdivision : $tf")
                    val start = -overhang
                    val end = tf + overhang


                    // * generate data of smoothed out pulse
                    val smoothEps = if (useSmooth) {
                        val sm = Smooth(smoothStd)
                        val dt = min(smoothStd / 10.0, tf / 1000.0)

                        List(((end - start) / dt + 1).toInt()) {
                            Pair(
                                it * dt + start,
                                sm(it * dt + start) { t -> rawEps(t, tf) }
                            )
                        }
                    } else null

                    // * function to call smoothed pulse shape using (linear) interpolation of smoothEps data
                    fun smoothEps(t: Double, tau: Double): Double {
                        return if (!useSmooth) rawEps(t, tau)
                        else interpolate(smoothEps!!, t, epsBounds)
                    }




                    // * generate the noise for the simulation
                    // * the spacing is chosen such that the discretization frequency is 10 x the noise high-cutoff
                    // * frequency
//                    val noise = Noise(time = end - start, spacing = 1.0 / (noiseHighCutoff * 10))
                    val noise = Noise(time = end - start, initialSpacing = noiseType.initialSpacing, wnVariance = noiseType.wnDeltaRate)
//                    noise.generate(envelope)
                    noise.generate(noiseType::envelope)



                    // * generate initial density matrix depending on the chose initial eigenstate
                    val initSys = hamiltonian(smoothEps(start, tf) + noise(0.0)).eigenSystem()
                    val rhoInitial = initSys[initial].second.normalized().ketBra()

                    // * generate the density matrices of the final eigenstates
                    val finalSys = hamiltonian(smoothEps(end, tf) + noise(end-start)).eigenSystem()
                    val rho0 = finalSys[0].second.normalized().ketBra()
                    val rho1 = finalSys[1].second.normalized().ketBra()
                    val rho2 = finalSys[2].second.normalized().ketBra()




                    // * define the ODE for the von Neumann equation
                    val ode =
                        if (collapseGamma == 0.0) {
                            fun(t: Double, rho: ComplexMatrix) =
                                -I * commutator(rho, hamiltonian(smoothEps(t, tf) + noise(t - start)))
                        } else {
                            // ? using approximate eigenstates
//                            fun(t: Double, rho: ComplexMatrix): ComplexMatrix {
////                                val sys = hamiltonian(eps(t, tf)).eigenSystem()
////                                val sMin = sys[0].second.normalized()
////                                val sPls = sys[2].second.normalized()
//
//                                val ep = eps0(t,tf) + noise(t - start).real
//                                val beta = sqrt(ep*ep + om*om)
//                                val sNeg = complexArrayOf(-(ep+beta)/om, 0.0, 1.0)
//                                val sPos = complexArrayOf(-(ep-beta)/om, 0.0, 1.0)
//
//                                val op = sNeg.ket() / sPos
//
//                                return -I * commutator(rho, hamiltonian(eps0(t, tf) + noise(t - start).real)) + collapse(op, rho) * collapseGamma.toComplex()
//                            }

                            // ? using exact eigenstates
                            fun(t: Double, rho: ComplexMatrix): ComplexMatrix {
                                val H = hamiltonian(smoothEps(t, tf) + noise(t - start))

                                val sys = H.eigenSystem()
                                val sNeg = sys[0].second.normalized()
                                val sPos = sys[2].second.normalized()

                                val op = sNeg.bra() / sPos

                                return -I * commutator(rho, H) + collapseGamma * collapse(op, rho)
                            }
                        }

                    val (_, rRho) = solveEvolution(start, end, rhoInitial, ode)

                    p0[sample][itf] = (rho0 * rRho).trace().real
                    p1[sample][itf] = (rho1 * rRho).trace().real
                    p2[sample][itf] = (rho2 * rRho).trace().real
                }
            }
        }
        futures.forEach { it.awaitAll() }
    }
    println("solver finished.")

    val res = mutableListOf(tf)

    for (i in 0 until samples) {
        res.add(p0[i])
        res.add(p1[i])
        res.add(p2[i])
        if (i == 0) continue // ! what does this?
    }


    if (saveData) {
        save(filename, res)
        println("saved as: $filename")
    }


    if (plotData) {
        val plot2 = Plot(
            DataSet(tf.zip(p0[0]) { x, y -> listOf(x, y) }, symbol = Symbols.None, color = Plot.colorPalette[0]),
            xScale = Linear, yScale = Log10, minY = 1e-6, maxY = 1.0
        )
        plot2.addDataSet(
            DataSet(
                tf.zip(p1[0]) { x, y -> listOf(x, y) },
                symbol = Symbols.None,
                color = Plot.colorPalette[1]
            )
        )
        plot2.addDataSet(
            DataSet(
                tf.zip(p2[0]) { x, y -> listOf(x, y) },
                symbol = Symbols.None,
                color = Plot.colorPalette[2]
            )
        )

        for (i in 0 until samples) {
            plot2.addDataSet(
                DataSet(
                    tf.zip(p0[i]) { x, y -> listOf(x, y) },
                    symbol = Symbols.None,
                    color = Plot.colorPalette[0]
                )
            )
            plot2.addDataSet(
                DataSet(
                    tf.zip(p1[i]) { x, y -> listOf(x, y) },
                    symbol = Symbols.None,
                    color = Plot.colorPalette[1]
                )
            )
            plot2.addDataSet(
                DataSet(
                    tf.zip(p2[i]) { x, y -> listOf(x, y) },
                    symbol = Symbols.None,
                    color = Plot.colorPalette[2]
                )
            )
        }
    }

}


fun singleSweep() {
    val initial = zeroState(3)
    initial[0] = Complex(1.0, 0.5)
    initial[1] = Complex(0.2, 0.6)
    var rhoInit = initial.normalized().ketBra()

    println("basis states")
    val rho0 = state(1.0, 0.0, 0.0).ketBra()
    val rho1 = state(0.0, 1.0, 0.0).ketBra()
    val rho2 = state(0.0, 0.0, 1.0).ketBra()

    println("Hamiltonian")
    val dbz = 0.1
    val om = 1.0
    fun eps(t: Double, tau: Double) = 100 - t / tau * 200

    val hamiltonian = Operator(
        Pair(3, 3), complexArrayOf(
            0.0, dbz, 0.0,
            dbz, 0.0, om,
            0.0, om, 0.0
        )
    )

    val ode = fun(t: Double, rho: ComplexMatrix): ComplexMatrix {
        hamiltonian[2, 2] = eps(t, 10.0).toComplex()
        return -I * commutator(rho, hamiltonian)
    }

    val tSubdivision = 10000
    val t = List(tSubdivision) {
        it.toDouble() / tSubdivision.toDouble() * 10.0
    }

    val p0 = MutableList(tSubdivision) { List(2) { 0.0 } }
    val p1 = MutableList(tSubdivision) { List(2) { 0.0 } }
    val p2 = MutableList(tSubdivision) { List(2) { 0.0 } }


    println("starting solver...")
    for (i in 1 until t.size) {
//        val prev = t[i-1]
//        val t = (i.toDouble() + 1.0)/N.toDouble() * 10

        val (tf, new) = solveEvolution(t[i - 1], t[i], rhoInit, ode)

        p0[i] = listOf(tf, (rho0 * new).trace().real)
        p1[i] = listOf(tf, (rho1 * new).trace().real)
        p2[i] = listOf(tf, (rho2 * new).trace().real)

        rhoInit = new
    }

    println("solver finished.")
//    res.forEach { println(it) }

    val plot = Plot(DataSet(p0, symbol = Symbols.None))
    plot.addDataSet(DataSet(p1, symbol = Symbols.None))
    plot.addDataSet(DataSet(p2, symbol = Symbols.None))
    plot.addDataSet(DataSet(listOf(listOf(0.0, 1.0), listOf(0.0, 0.0))))
}


fun pulseShape() {

    val om = 1.0

    val a = 0.1 * om
    val epsBounds = Pair(-10.0 * om, 10.0 * om)

    val deps = fun(ep: Double): Double {
        val beta = sqrt(ep * ep + om * om)
        val x = beta.e(3) * sqrt(1.0 + ep / beta)
        val y = a * a + 2 * om * om + 4 * ep * ep + ep / beta * (a * a - 4 * om * om - 4 * ep * ep)
        val z = 2 * a * om * om * (ep - 3 * beta)
        return -x * y.e(3.0 / 2.0) / z
    }

    val pulse = calculatePulseShape(epsBounds, deps)
    println(pulse)

    val plot = Plot(DataSet(pulse.map { listOf(it.first, it.second) }))

    val tInterpolated = linspace(0.0, pulse.last().first, 10000)
    val interpolated = tInterpolated.map { listOf(it, interpolate(pulse, it, epsBounds)) }
    plot.addDataSet(DataSet(interpolated, color = Plot.colorPalette[1], symbol = Symbols.None))
}


fun sweepTau() {


    val om = 1.0

    val a = 0.1 * om
    val g = 0.0 * om
    val epsMax = 10.0 * om
    val epsMin = -10.0 * om

    fun eps(t: Double, tau: Double): Double {
        val tt = when {
            t < 0.0 -> 0.0
            t > tau -> tau
            else -> t
        }

        return (epsMax - epsMin) * tt / tau + epsMin
    }

    fun deps(t: Double, tau: Double): Double {
        if (t < 0 || t > tau) return 0.0
        return (epsMax - epsMin) / tau
    }

    val ss = DonorDotSS(om, a, g, ::eps, ::deps)

    val tau = 1 * om
    val t = linspace(0.0, tau, 20)
    for (tt in t) {
        println("$tt : ${ss.phi(tt, tau)}")
    }

    val phiData = t.map { listOf(it, ss.phi(it, tau)) }
    Plot(DataSet(phiData))

    val initialSS = ss.initRelB(doubleArrayOf(1.0, 0.0, 0.0))
    println("initial state = ${initialSS.toList()}")
    val res = ss.integrate(initialSS, tau)
    println("res = ${res.toList()}")

    val aTau = linspace(-2.0, 3.0, 500).map { 10.0 e it }

    val st = DonorDotST(om, a, g, ::eps, ::deps)
    val initialST = st.initRelB(doubleArrayOf(1.0, 0.0, 0.0))
//    val initialST = st.initRelB(doubleArrayOf(1.0,0.0,0.0))


    val (resSS, resST) = runBlocking(Dispatchers.Default) {
        val rSS = aTau.pmap { ss.integrate(initialSS, it) }
//        val rSS = sampling(
//            aTau,
//            initial = initialSS,
//            f_integrated = ss::integrate
//        )
        val rST = aTau.pmap { st.integrate(initialST, it) }
//        val rST = sampling(
//            aTau,
//            initial = initialST,
//            f_integrated = st::integrate
//        )
        Pair(rSS, rST)
    }
    val zErrSS = resSS.map { 0.5 * (1 - it[0]) }
    println(zErrSS)
    val zDataSS = aTau.mapIndexed { index, d -> listOf(d, zErrSS[index]) }

    val zErrST = resST.map { 0.5 * (1 - it[0]) }
    val zDataST = aTau.mapIndexed { index, d -> listOf(d, zErrST[index]) }

    val plot0 = Plot(DataSet(zDataSS, symbol = Symbols.None), xScale = Log10, yScale = Log10)
    plot0.addDataSet(DataSet(zDataST, color = Plot.colorPalette[1], symbol = Symbols.None))


//    val zss = resSS.map{it[0]}
//    val zdatass = aTau.mapIndexed{index, d -> listOf(d, zss[index])}
//    val plot0 = Plot(Plot.DataSet(zdatass, symbol = Plot.Symbols.None))
}

