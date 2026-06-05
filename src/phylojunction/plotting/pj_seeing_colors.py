import matplotlib
import numpy as np
import matplotlib.colors as colors
from matplotlib import pyplot as plt


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap


# matplotlib color palette name, n colors
cmap = plt.get_cmap('seismic', 5)
for i in range(cmap.N):
    rgb = cmap(i)[:3]  # will return rgba, we take only first 3 so we get rgb
    print(matplotlib.colors.rgb2hex(rgb))

qual_cmap = plt.get_cmap('tab20')
qual_color_list = \
    [matplotlib.colors.rgb2hex(qual_cmap(i)[:3]) for i in range(qual_cmap.N)]

# n_colors = 10
n_colors = 20
# n_colors = 120
mv = 120 * 2.08 / n_colors  # I found this 2.08 empirically... by trying!

# this palette does up to 250 colors, not more
cmap = plt.get_cmap('terrain', n_colors)

# for up to 120
# new_cmap = truncate_colormap(cmap, minval=0.0, maxval=2.08, n=n_colors)

# for (120, 250)
new_cmap = truncate_colormap(cmap, minval=0.0, maxval=mv, n=n_colors)

# arr = np.linspace(0, 50, 100).reshape((10, 10))
# fig, ax = plt.subplots(ncols=2)
# ax[0].imshow(arr, interpolation='nearest', cmap=cmap)
# ax[1].imshow(arr, interpolation='nearest', cmap=new_cmap)
# plt.show()

# use this for n_colors = 10
# color_list = \
#     [matplotlib.colors.rgb2hex(new_cmap(i)[:3])
#      for i in [0,9,19,29,39,49,59,69,79,89,99]]

# use this for n_colors = 20
# color_list = \
#     [matplotlib.colors.rgb2hex(new_cmap(i)[:3])
#      for i in [0,4,9,14,19,24,29,34,39,44,49,54,59,64,69,74,79,84,89,94,99]]

# use this for n_colors >= 120
# color_list = \
#     [matplotlib.colors.rgb2hex(new_cmap(i)[:3]) for i in range(cmap.N)]

# print(color_list) # to see the hex codes

if __name__ == "__main__":
    fig, ax = matplotlib.pyplot.subplots(figsize=(8, 5))
    # for x in range(n_colors):
    #     # ax.axvline(x, color=color_list[x], linewidth=3)
    #     ax.axvline(x, color=qual_color_list[x], linewidth=3)
    # plt.show()

    # colors = plt.cm.Vega20c((4. / 3 * np.arange(20 * 3 / 4)).astype(int))
    plt.scatter(np.arange(20), np.ones(20), c=qual_color_list, s=180)
    print(qual_color_list)

    plt.show()
