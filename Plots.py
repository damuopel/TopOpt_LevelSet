import gif
# Classes
class Plots():   
    @gif.frame
    def MaterialDistribution(self,Top,):
        # Plot material distribution
        plt.pcolor(top.x.reshape(ny,nx),cmap='cool',edgecolors='k',linewidth=0.1)
        plt.colorbar()
        ax = plt.gca()
        ax.axis('equal')
        ax.axis('off')