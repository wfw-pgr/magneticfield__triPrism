import sys
import numpy                      as np
import nkUtilities.load__config   as lcf
import nkUtilities.plot1D         as pl1
import nkUtilities.configSettings as cfs


# ========================================================= #
# ===  display                                          === #
# ========================================================= #
def display():

    xs_,ys_ = 2, 3
    xt_,yt_ = 1, 2

    # ------------------------------------------------- #
    # --- [1] Arguments                             --- #
    # ------------------------------------------------- #
    config   = lcf.load__config()
    datFile1 = "dat/subdivE.dat"
    datFile2 = "dat/triangles.dat"
    pngFile  = "png/subdivE.png"

    # ------------------------------------------------- #
    # --- [2] Fetch Data                            --- #
    # ------------------------------------------------- #
    import nkUtilities.load__pointFile as lpf
    Data1    = lpf.load__pointFile( inpFile=datFile1, returnType="point" )
    Data2    = lpf.load__pointFile( inpFile=datFile2, returnType="point" )
    
    # ------------------------------------------------- #
    # --- [3] config Settings                       --- #
    # ------------------------------------------------- #
    cfs.configSettings( configType="plot.def"   , config=config )
    cfs.configSettings( configType="plot.marker", config=config )
    config["xTitle"]         = "X (m)"
    config["yTitle"]         = "Y (m)"
    config["plt_xAutoRange"] = False
    config["plt_yAutoRange"] = False
    config["plt_xRange"]     = [-0.6,+0.6]
    config["plt_yRange"]     = [-0.6,+0.6]
    config["plt_linewidth"]  = 1.0
    config["xMajor_Nticks"]  = 7
    config["yMajor_Nticks"]  = 7

    # ------------------------------------------------- #
    # --- [4] plot Figure                           --- #
    # ------------------------------------------------- #
    fig = pl1.plot1D( config=config, pngFile=pngFile )
    fig.add__plot( xAxis=Data1[:,xs_], yAxis=Data1[:,ys_], marker=".", linestyle="none", color="RoyalBlue" )
    fig.add__plot( xAxis=Data2[:,xt_], yAxis=Data2[:,yt_], marker=".", linestyle="none", color="Magenta"   )
    fig.set__axis()
    fig.save__figure()


# ======================================== #
# ===  実行部                          === #
# ======================================== #
if ( __name__=="__main__" ):
    display()

