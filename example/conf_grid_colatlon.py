import numpy as np

# Generate 600 points between 0 and 17.2512 (inclusive)
lon_points = np.linspace(0, 17.2512, 600)

# Generate 200 points between 84.2688 and 90 (inclusive)
colat_points = np.linspace(84.2688, 90, 200)

# Save the points to files
np.savetxt('./grid_lon.dat', lon_points, fmt='%f')
np.savetxt('./grid_colat.dat', colat_points, fmt='%f')

# Verify the first few points for confirmation
lon_points[:5], colat_points[:5]
