using PythonCall
ENV["MPLBACKEND"] = "Qt5Agg"
using PythonPlot
pyplot.ion()

function test_plot()
    # Create some data
    x = 0:0.1:2π
    y = sin.(x)
    fig, axes = pyplot.subplots(2,2)
    figurename = "Test Plot"
    fig.suptitle(figurename, fontsize=16)
    for ax in axes.flat
    # Create a plot using Matplotlib functions
        ax.plot(x, y, label="sin(x)")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_title("Sine Wave Plot")
        ax.legend()
        ax.grid(true)
    end
    pyplot.tight_layout()
    pyplot.savefig("sine_wave_plot.png")
    # Display the plot (if not in an interactive environment like IJulia)
    pyplot.show()
    println("finis")
end

test_plot()