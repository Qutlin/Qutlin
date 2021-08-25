package Qutlin

import kotlinx.coroutines.Dispatchers
import kotlinx.coroutines.runBlocking
import org.hipparchus.analysis.UnivariateFunction
import org.hipparchus.complex.Complex
import org.hipparchus.util.FastMath
import org.hipparchus.util.FastMath.PI
import org.hipparchus.util.FastMath.sqrt
import pmap
import java.util.function.Function

const val _ns : Double = 1.0
const val _s = 1.0e9 * _ns
const val _Hz = 1.0/_s
const val _kHz = 1e3 * _Hz
const val _MHz = 1e3 * _kHz
const val _μeV : Double = 1.0
const val _eV = 1e6 * _μeV
const val _h = 4.135667 * 1e-15 * _eV * _s
const val _ħ = 6.582119 * 1e-16 * _eV * _s
const val _T : Double = 1.0
const val _μB = 5.788382 * 1e-5 * _eV / _T

const val TAU = 2.0 * PI
const val π = PI


infix fun Double.e(b: Double) = FastMath.pow(this, b)
fun Double.e2() = FastMath.pow(this, 2)
fun Double.e3() = FastMath.pow(this, 3)
infix fun Double.e(b: Int) = FastMath.pow(this, b)




fun List<Complex>.average(): Complex =
    this.fold(Complex.ZERO){acc, c -> acc + c} / size.toDouble()




fun List<Double>.variance(): Double {
    val avg = average()
    return fold(0.0){acc, v -> acc + (v-avg)*(v-avg)} / size //(size - 1.0)
}

fun List<Double>.std(): Double = sqrt(variance())

fun List<Double>.pAutoCorrelation(normalized: Boolean = true) = runBlocking(Dispatchers.Default) {
    val tmp = List(size) {it}
    val avg = average()
    val σ2 = if(normalized) variance() else 1.0
    tmp.pmap {
        mapIndexed { index: Int, d: Double ->
            if(index+it<size)
                    (d-avg) * (get(index+it)-avg)
            else
                0.0
        }.sum() / (σ2 * size)//(size-1))
    }
}



fun ((Double) -> Double).uvFun() = UnivariateFunction{ x -> this.invoke(x) }