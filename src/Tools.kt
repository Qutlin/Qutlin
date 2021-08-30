import kotlinx.coroutines.async
import kotlinx.coroutines.awaitAll
import kotlinx.coroutines.coroutineScope
import org.hipparchus.special.Erf.erfInv
import org.hipparchus.util.FastMath.*
import java.io.File




suspend fun <A, B> Iterable<A>.pmap(f: suspend (A) -> B): List<B> = coroutineScope {
    map { async { f(it) } }.awaitAll()
}

suspend fun <A, B> Iterable<A>.pmapIndexed(f: suspend (Int, A) -> B): List<B> = coroutineScope {
    mapIndexed { i,v -> async { f(i, v) } }.awaitAll()
}



fun <T> concatenate(vararg lists: List<T>) = listOf(*lists).flatten()



fun save(path: String, data: DoubleArray, newline: Int) {
    try {
        File(path).bufferedWriter().use { out ->
            data.forEachIndexed() {i, v ->
                if(i%newline == 0) out.newLine()
                out.write(v.toString())
                if( (i+1)%newline != 0) out.write(", ")
            }
        }
    } catch (e: Exception) {
        e.printStackTrace()
    }
}

fun save(path: String, data: List<List<Double>>) {
    try {
        File(path).bufferedWriter().use { out ->
            data.forEachIndexed { i,d ->
                if(i > 0) out.newLine()
                d.forEachIndexed {ii, dd ->
                    if(ii < d.lastIndex) out.write("$dd, ")
                    else out.write(dd.toString())
                }
            }
        }
    } catch (e: Exception) {
        e.printStackTrace()
    }
}






fun linsteps(min: Double, step: Double, max: Double): DoubleArray {
    val res = mutableSetOf(min)
    var x = min
    while(x < max) {
        x += step
        res.add(x)
    }
    if (x != max) {res.add(max)}
    return res.toDoubleArray()
}

fun linspace(min: Double, max: Double, subdivisions: Int, skipFirst: Boolean = false, skipLast: Boolean = false): DoubleArray {
    val res = if(skipFirst) mutableListOf<Double>() else mutableSetOf(min)
    var x = min
    for(i in 0 until subdivisions-1) {
        x += (max - min)/subdivisions.toDouble()
        res.add(x)
    }
    if(!skipLast) res.add(max)
    return res.toDoubleArray()
}

fun clamp(min: Double, x: Double, max: Double): Double {
    if(x < min) return min
    if(x > max) return max
    return x
}

/**
 * Normalized Gaussian function
 *
 * `x0` : average
 *
 * `σ` : standard deviation
 */
fun gaussian(x: Double, x0: Double = 0.0, σ: Double = 1.0)
    = exp(-0.5 * pow((x-x0)/σ, 2.0)) / (σ * sqrt(2.0*PI))


/**
 * Helper class to smoothen out a function ___f(x)=y___ via a Gaussian convolution.
 * 
 * Analytically, the result is given by ___f'(x) = ∫f(x₀+x)g(x)dx___ (given ___g(x)___ is symmetric), where ___dx∈(-∞,∞)___.
 * 
 * The numerical version is ___f'(x) ≈ ∑ᵢf(x₀+xᵢ)gᵢΔᵢ___ Here, ___xᵢ___ is the ___i___-th discrete point on the ___x___-axis, ___Δᵢ = xᵢ₊₁-xᵢ___, and ___gᵢ=g(xᵢ)/∑ₙg(xₙ)___, see below.
 * 
 * The points ___xᵢ___ are distributed with higher density at the center of the gaussian than at the tails which allows for smaller numerical errors (using the inverse error function `errInv`):
 * 
 * ___αᵢ = erfInv(2 i/(N+1) - 1) √2 δx___,
 * 
 * ___xᵢ = g(αᵢ)___,
 * 
 * where ___δx___ is the standard deviation, ___N___ the number of subdivisions, and ___i=1..N___.
 * 
 * The values ___g(xᵢ)___ are finally normalized ___gᵢ=g(xᵢ)/∑ₙg(xₙ)∆ₙ___.
 *
 * `erf` is defined as ___erf(x) = 2/√π ∫₀ˣ exp(-t²) dt___
 */
class Smooth(
    private val std: Double,
    private val subdivisions: Int = 1001,
) {
    private var gx: List<Double>
    private var gy: List<Double>
    init{
        gx = List(subdivisions){ i ->
            val r = (i+1).toDouble() / (subdivisions.toDouble() + 1.0)
            erfInv(2 * r - 1) * sqrt(2.0) * std
        }
        gy = gx.take(gx.size - 1).mapIndexed { i,v -> gaussian(v, 0.0, std) * (gx[i+1]-v) }
        val sum = gy.sum()
        gy = gy.map { it/sum }
    }

    operator fun invoke(x: Double, f: (Double) -> Double)
        = gx.zip(gy) {ix,iy -> f(x+ix) * iy}.sum()
}




fun main() {
    val x = linspace(-10.0, 10.0, 200)
    fun f(x: Double) = signum(x)
    val y0 = x.map { f(it) }

    val s1 = Smooth(2.0)
    val y1 = x.map { s1(it, ::f) }

    val p1 = Plot(
        Plot.DataSet(
            x.zip(y0) {x,y -> listOf(x,y)},
            color = Plot.colorPalette[0],
            symbol = Plot.Symbols.None
        ),
        minY = -1.5,
        maxY = 1.5,
    )
    p1.addDataSet(Plot.DataSet(
        x.zip(y1) {x,y -> listOf(x,y)},
        color = Plot.colorPalette[1],
        symbol = Plot.Symbols.None
    ))

    // ? analytical check: A step function of height 2 (-1 to 1)
    // ? leads to a gradient of 2*gaussian(0.0, 0.0, std)
    val yCheck = x.map { it * 2.0 * gaussian(0.0, 0.0, 2.0) }
    p1.addDataSet(Plot.DataSet(
        x.zip(yCheck) { x, y -> listOf(x,y)},
        color = Plot.colorPalette[2],
        symbol = Plot.Symbols.None
    ))
}