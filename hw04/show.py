import sys
import matplotlib.pyplot as plt

v = [l.strip() for l in sys.stdin]
plt.title(v[0])
plt.hist(
    [float(x) for i, x in enumerate(v) if i > 0], bins=100, density=True, color="b"
)
plt.xlabel("Random x")
plt.ylabel("Probability density")
plt.show()
