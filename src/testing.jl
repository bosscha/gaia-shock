## testing function (for notebooks e.g.)

function __plot_check(dfcart,plotdir)
    cart= DataFrame(X=dfcart.data[1,:], Y=dfcart.data[2,:], Z=dfcart.data[3,:])

    println("## check plot subtraction ...")

    PyPlot.plt.figure(figsize=(9.0,8.0))
    PyPlot.plt.subplot(1, 1, 1 , ylim=[100,450], xlim=[-50,50])
    PyPlot.plt.scatter(cart.Y, cart.X, s = 0.1 )
    PyPlot.plt.xlabel("Y (pc)")
    PyPlot.plt.ylabel("X (pc)")
    PyPlot.plt.grid(true)

    PyPlot.plt.savefig(plotdir*"/"*plotfile)
    PyPlot.plt.show()

end

function __plot_nstars(nstarh,plotfile="test-stats-votable.png", plotdir= ".")
    println("## plotting distribution...")
    PyPlot.plt.figure(figsize=(9.0,8.0))
    PyPlot.plt.subplot(1, 1, 1 )
    nbins = 50
    PyPlot.plt.hist(nstarh,nbins, range = [0,3e5],  color = "g", alpha=0.8 , label = "Votable stars",density=false)
    PyPlot.plt.xlabel("Stars")
    PyPlot.plt.ylabel("N")
    PyPlot.plt.grid(true)
    PyPlot.plt.savefig(plotdir*"/"*plotfile)
    PyPlot.plt.show()
end
