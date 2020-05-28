import cut_plane_plot as cp

#cp.plot([2003, 11, 20, 7, 0, 0, 176., 57.50], 'p', png=False)

cp.plot([2003, 11, 20, 7, 0, 0, 176., 57.50], 'p',
        nx=100, ny=100, xlim=[-30, 30], ylim=[-20, 20], png=False)

# plot(time, [r, mlat, mlong], None, 'p')
# plot(time, [GSMx,GSMy,GSMz], [v1, v2], 'p')

# plot(time, [0,0,0], [[1, 0, 0], [0, 0, 1]], 'p') # GSM y-z plane
