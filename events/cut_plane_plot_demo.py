import cut_plane_plot as cp

#cp.plot([2003, 11, 20, 7, 0, 0, 176., 57.50], 'p', png=False)

#cp.plot([2003, 11, 20, 7, 0, 0, 176., 57.50], 'p',
#        nx=100, ny=100, xlim=[-30, 30], ylim=[-20, 20], png=False)

# plot(time, [r, mlat, mlong], None, 'p')
# plot(time, [GSMx,GSMy,GSMz], [v1, v2], 'p')

# plot(time, [0,0,0], [[1, 0, 0], [0, 0, 1]], 'p') # GSM x-z plane

time = [2003, 11, 20, 7, 0, 0]
r = 1.01
mlat = 57.50
mlong = 176.
cp.plot(time, [r, mlat, mlong], None, 'p',
        nx=100, ny=100, xlim=[-30, 30], ylim=[-20, 20], png=False)

#cp.plot(time, [0 ,0 ,0 ], [[0, 1, 0], [0, 0, 1]], 'p',
#        nx=100, ny=100, xlim=[-30, 30], ylim=[-20, 20], png=False)

# plot(time, 'xy', None, 'p')
# plot(time, 'xz', None, 'p')
# plot(time, 'yz', None, 'p')


