import org.hipparchus.util.FastMath
import org.hipparchus.util.FastMath.tanh
import java.awt.*
import java.awt.event.WindowEvent
import java.awt.event.WindowListener
import java.awt.geom.AffineTransform
import javax.swing.JFrame


/*
888b     d888
8888b   d8888
88888b.d88888
888Y88888P888  8888b.  88888b.  88888b.   .d88b.  888d888
888 Y888P 888     "88b 888 "88b 888 "88b d8P  Y8b 888P"
888  Y8P  888 .d888888 888  888 888  888 88888888 888
888   "   888 888  888 888 d88P 888 d88P Y8b.     888
888       888 "Y888888 88888P"  88888P"   "Y8888  888
                       888      888
                       888      888
                       888      888
*/



class Mapper(
    private val data: List<Plot.DataSet>,
    var width: Int, var height: Int
) {
    companion object {
        enum class Scales {
            Linear, Log10, Tanh//, SymLog10
        }
    }

    var xScale = Scales.Linear
    var yScale = Scales.Linear
    var minX: Double? = null
    var maxX: Double? = null
    var minY: Double? = null
    var maxY: Double? = null

    var flipy = true

    private val xBound = mutableListOf(0.0, 0.0)
    private val yBound = mutableListOf(0.0, 0.0)

    init {
        updateBounds()
    }

    fun updateBounds() {
        xBound[0] = data[0].points[0][0]
        xBound[1] = data[0].points[0][0]
        yBound[0] = data[0].points[0][1]
        yBound[1] = data[0].points[0][1]

        data.forEach { outer ->
            outer.points.forEach {
                if (it[0] < xBound[0]) xBound[0] = it[0]
                if (it[0] > xBound[1]) xBound[1] = it[0]
                if (it[1] < yBound[0]) yBound[0] = it[1]
                if (it[1] > yBound[1]) yBound[1] = it[1]
            }
        }
//        println("Mapper: bounds updated to: $xBound, $yBound")
        if(minX != null) xBound[0] = minX!!
        if(maxX != null) xBound[1] = maxX!!
        if(minY != null) yBound[0] = minY!!
        if(maxY != null) yBound[1] = maxY!!
    }

    fun map(x: Double, y: Double): Pair<Int, Int> {
        val xres = when (xScale) {
            Scales.Linear -> (width * (x - xBound[0]) / (xBound[1] - xBound[0])).toInt()
            Scales.Log10 -> {
                val mx = FastMath.log10(xBound[1])
                val mn = FastMath.log10(xBound[0])
                val lx = FastMath.log10(x)
                (width * (lx - mn) / (mx - mn)).toInt()
            }
            Scales.Tanh -> {
                val mn = tanh(xBound[0])
                val mx = tanh(xBound[1])
                (width * (tanh(x) - mn) / (mx - mn)).toInt()
            }
        }
        var yres = when (yScale) {
            Scales.Linear -> (height * (y - yBound[0]) / (yBound[1] - yBound[0])).toInt()
            Scales.Log10 -> {
                val mx = FastMath.log10(yBound[1])
                val mn = FastMath.log10(yBound[0])
                val ly = FastMath.log10(y)
//                println("${yBound[1]}, ${yBound[0]} -> mx = $mx, my = $mn")
//                println("mapped $y -> $ly")
                (height * (ly - mn) / (mx - mn)).toInt()
            }
            Scales.Tanh -> {
                val mn = tanh(xBound[0])
                val mx = tanh(xBound[1])
                (width * (tanh(y) - mn) / (mx - mn)).toInt()
            }
        }
        if(flipy) yres = height - yres

//        println("mapping ($x, $y) -> ($xres, $yres)")

        return Pair(xres, yres)
    }
}







/*
8888888b.  888          888
888   Y88b 888          888
888    888 888          888
888   d88P 888  .d88b.  888888
8888888P"  888 d88""88b 888
888        888 888  888 888
888        888 Y88..88P Y88b.
888        888  "Y88P"   "Y888



*/



class Plot(
    points: DataSet,
    var xScale: Mapper.Companion.Scales = Mapper.Companion.Scales.Linear,
    var yScale: Mapper.Companion.Scales = Mapper.Companion.Scales.Linear,
    var minX: Double? = null,
    var maxX: Double? = null,
    var minY: Double? = null,
    var maxY: Double? = null
    ) : JFrame("Plot") {

    companion object {
        var counter = 0
        var colorPalette = mutableListOf( //* matplotlib colors
            Color(31, 119, 180),
            Color(255, 127, 14),
            Color(44, 160, 44),
            Color(214, 39, 40),
            Color(148, 103, 189)
        )
    }

    enum class Symbols {
        None, CircleFilled, Circle
    }

    enum class Lines {
        None, Solid
    }

    enum class ColorData {
        None, Data
    }

    class DataSet(
        val points: List<List<Double>>,
        val colordata: ColorData = ColorData.None,
        val color: Color = colorPalette[0],
        val symbol: Symbols = Symbols.Circle,
        val symbolSize: Int = 3,
        val line: Lines = Lines.Solid,
        val lineWidth: Float = 1f) {
        lateinit var displayPoints: List<Pair<Int, Int>>
    }


    private val data = mutableListOf<DataSet>()

    private var initialized = false
    private var mapper: Mapper

    private var psize = Pair(0,0)
    private val padding = Pair(50,50)

    private val superSampling = 1
    private lateinit var offscreen: Image
    private lateinit var bg: Graphics2D

    init {
        data.add(points)

        addWindowListener(object : WindowListener {
            override fun windowDeiconified(e: WindowEvent?) {}
            override fun windowClosed(e: WindowEvent?) {}
            override fun windowActivated(e: WindowEvent?) {}
            override fun windowDeactivated(e: WindowEvent?) {}
            override fun windowOpened(e: WindowEvent?) {}
            override fun windowIconified(e: WindowEvent?) {}
            override fun windowClosing(e: WindowEvent?) {
                println("closing: $title")
                dispose()
            }
        })

        title = "Plot $counter"
        counter++

        setSize(800, 600)
        isVisible = true

        val d = contentPane.size
        psize = Pair(d.width, d.height)

        mapper = Mapper(data, psize.first - 2 * padding.first, psize.second - 2 * padding.second)
        mapper.xScale = xScale
        mapper.yScale = yScale
        mapper.minX = minX
        mapper.maxX = maxX
        mapper.minY = minY
        mapper.maxY = maxY

//        offscreen = createImage(psize.first, psize.second)
//        bg = offscreen.graphics

        initialized = true
        repaint()
    }



    fun addDataSet(dset: DataSet) {
        data.add(dset)
        repaint()
    }



    private fun update() {
        val d = contentPane.size
        psize = Pair(d.width, d.height)

        mapper.width = psize.first - 2 * padding.first
        mapper.height = psize.second - 2 * padding.second
        mapper.updateBounds()
    }

    private fun setupGraphics(): Pair<Graphics2D, Graphics2D> {
        offscreen = createImage(psize.first * superSampling, psize.second * superSampling)
        bg = offscreen.graphics as Graphics2D

        if(superSampling != 1)
            bg.transform = AffineTransform.getScaleInstance(superSampling.toDouble(), superSampling.toDouble())

        val cpg = contentPane.graphics as Graphics2D

//        val g = canvas.graphics
        bg.setRenderingHint(
            RenderingHints.KEY_ANTIALIASING,
            RenderingHints.VALUE_ANTIALIAS_ON
        )
//        (bg as Graphics2D).setRenderingHint(
//            RenderingHints.KEY_ALPHA_INTERPOLATION,
//            RenderingHints.VALUE_ALPHA_INTERPOLATION_QUALITY
//        )
        bg.setRenderingHint(
            RenderingHints.KEY_RENDERING,
            RenderingHints.VALUE_RENDER_QUALITY
        )
//
        cpg.setRenderingHint(
            RenderingHints.KEY_RENDERING,
            RenderingHints.VALUE_RENDER_QUALITY
        )
        cpg.setRenderingHint(
            RenderingHints.KEY_ANTIALIASING,
            RenderingHints.VALUE_ANTIALIAS_ON
        )
        cpg.setRenderingHint(
            RenderingHints.KEY_INTERPOLATION,
            RenderingHints.VALUE_INTERPOLATION_BICUBIC
        )
//        (cpg as Graphics2D).setRenderingHint(
//            RenderingHints.KEY_ALPHA_INTERPOLATION,
//            RenderingHints.VALUE_ALPHA_INTERPOLATION_QUALITY
//        )
        return Pair(bg, cpg)
    }

    private fun applyImage(cpg: Graphics2D) =
        cpg.drawImage(offscreen, AffineTransform.getScaleInstance(1.0/superSampling, 1.0/superSampling), this)

    private fun drawData(bg: Graphics2D) {
//        println("drawData")
        data.forEach { outer ->
            outer.displayPoints = outer.points.map {
                val p = mapper.map(it[0], it[1])
                Pair(p.first + padding.first, p.second + padding.second)
            }
        }

        data.forEach { dset ->
            bg.color = dset.color
            dset.displayPoints.forEachIndexed { index, p ->
                if(index > 0) {
                    val pp = dset.displayPoints[index-1]
                    if(dset.line != Lines.None) {
                        bg.stroke = BasicStroke(dset.lineWidth)
                    }
                    when(dset.line) {
                        Lines.None -> {}
                        Lines.Solid -> bg.drawLine(pp.first, pp.second, p.first, p.second)
                    }
                }

                bg.stroke = BasicStroke()
                val symbolSize = dset.symbolSize
                when(dset.symbol) {
                    Symbols.None -> {}
                    Symbols.CircleFilled -> bg.fillOval(p.first - symbolSize, p.second - symbolSize, 2*symbolSize, 2*symbolSize)
                    Symbols.Circle -> bg.drawOval(p.first - symbolSize, p.second - symbolSize, 2*symbolSize, 2*symbolSize)
                }

            }
        }
    }

    private fun drawFrame(bg: Graphics2D) {
        bg.color = Color.black
        bg.drawRect(padding.first, padding.second,
            psize.first - 2 * padding.first, psize.second - 2 * padding.second)
    }

    override fun paint(g: Graphics) {
        if(!initialized) return

        update()
        val (bg, cpg) = setupGraphics()

        drawFrame(bg)
        drawData(bg)

        applyImage(cpg)
    }
}