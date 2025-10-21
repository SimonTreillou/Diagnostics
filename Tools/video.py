import matplotlib.animation as animation
from IPython.display import Video
import numpy as np
import matplotlib.pyplot as plt

def write_map_to_video(var,time,path,title=r'$\omega$ ($s^{-1})$',fps=5,vmin=-0.1,vmax=0.1,cmap="RdBu_r",label="Vorticity"):
    FFMpegWriter = animation.writers['ffmpeg']
    metadata = dict(title="Video", artist='Simon Treillou', comment='Allez le TFC') 
    writer = FFMpegWriter(fps=fps, metadata=metadata, codec='mpeg4')

    # --- Figure setup ---
    fig, ax = plt.subplots(figsize=(8,8), dpi=300)
    pcm = ax.pcolormesh(var[0], cmap=cmap, shading="auto", vmin=vmin, vmax=vmax)
    cbar = fig.colorbar(pcm, ax=ax, label=label)
    cbar.set_clim(vmin, vmax)  # fix limits for the colorbar
    ax.set_xlabel(r'$x$ (m)')
    ax.set_ylabel(r'$y$ (m)')
    tit=ax.set_title(title + f' at t = {time[0]:.2f} s',loc='left')

    def update(frame):
        pcm.set_array(var[frame].ravel())
        tit.set_text(title + f' at t = {time[frame]:.2f} s')
        return pcm, tit

    ani = animation.FuncAnimation(fig, update, frames=var.shape[0], interval=100, blit=True)

    # Save animation to file
    ani.save('./Videos/'+path+'.mp4', writer=writer)