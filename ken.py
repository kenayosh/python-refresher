import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

fig = plt.figure(figsize=(10, 10))
plt.plot()

plt.gca().add_patch(
    Rectangle((25, 50), 15, 15, fill=True, color="g", alpha=0.5, zorder=100, figure=fig)
)

plt.gca().add_patch(
    Rectangle((50, 100), 40, 80, angle=30, edgecolor="red", facecolor="none", lw=4)
)

plt.savefig("plot.png")
