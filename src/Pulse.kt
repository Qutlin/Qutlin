import org.hipparchus.migration.ode.FirstOrderDifferentialEquations
import org.hipparchus.ode.nonstiff.GraggBulirschStoerIntegrator
import org.hipparchus.util.FastMath.floor

fun calculatePulseShape(
    eps: Pair<Double,Double>,
    depsilon: (Double) -> Double) : List<Pair<Double, Double>> {

    val integrator = GraggBulirschStoerIntegrator(
        1e-16,
        1e-3,
        1e-8,
        1e-8
    )

    val eqs = object : FirstOrderDifferentialEquations {
        override fun computeDerivatives(t: Double, y: DoubleArray, yDot: DoubleArray) {
            yDot[0] = depsilon(y[0])
        }
        override fun getDimension(): Int = 1
    }

    var tStep = 0.001
    var tNow = 0.0

    val res0 = doubleArrayOf(eps.second)
    val res = DoubleArray(1)

    val resData = mutableListOf(Pair(tNow, eps.second))

    while(res0[0] > eps.first) {
        try {
            val tNext = integrator.integrate(eqs, tNow, res0, tNow-tStep, res)
            resData.add(Pair(tNext, res[0]))
            tNow = tNext
            res0[0] = res[0]
        } catch(e: Exception) {
            // ? here, an error in the integrator occured. This usually happens when the step size is too large.
            // ? By halving the time interval of the integration, the integrator is usually better able to handle the integration.
            // ? This is recursive, meaning `calculatePulseShape` will continuously divide the intervals by 2 until it's small enough for the integrator to succeed.
//            println("error occured. Splitting time steps. tStep = $tStep")
            tStep /= 2
        }
    }

    return resData
}


/**
 * Interpolates the data given at the chosen point, abiding by the outer bounds for the interpolated value
 * @param data monotone (in first parameter) increasing List of data-points
 * @param p point to calculate interpolation at
 * @param bounds Pair of lower and upper bound for the interpolated value
 */
fun interpolate(data: List<Pair<Double, Double>>, p: Double, bounds: Pair<Double, Double>): Double {
    fun bound(d: Double) = clamp(bounds.first, d, bounds.second)

    var point = p

    val d0t = data.first().first
    val d1t = data.last().first

    if(point == d0t) return bound(data.first().second)
    if(point == d1t) return bound(data.last().second)

    val (mi, ma) = if(d0t < d1t) Pair(d0t, d1t) else Pair(d1t, d0t)

    point = clamp(mi, point, ma)

    val rp = (point - d0t)/(d1t-d0t) * data.size
    val i0 = floor(rp).toInt()
    val i1 = floor(rp + 1.0).toInt()
    val r = rp - i0.toDouble()

    if(i1 >= data.size) return data.last().second
    return data[i0].second * r + data[i1].second * (1.0 - r)
}