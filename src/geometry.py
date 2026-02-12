# src/geometry.py
import numpy as np

class BluntBody:
    """Starship-like blunt body nose geometry."""
    
    def __init__(self, R_nose=4.5, theta_max=np.pi/2,
                 n_theta=50, n_phi=40):
        self.R_nose = R_nose
        self.theta_max = theta_max
        self.n_theta = n_theta
        self.n_phi = n_phi
        self.points, self.normals = (
            self._generate_mesh()
        )
    
    def _generate_mesh(self):
        """Generate surface mesh points and 
        outward normals on the spherical cap."""
        theta = np.linspace(
            0, self.theta_max, self.n_theta
        )
        phi = np.linspace(0, 2*np.pi, self.n_phi,
                          endpoint=False)
        TH, PH = np.meshgrid(theta, phi,
                              indexing='ij')
        
        x = self.R_nose * np.sin(TH) * np.cos(PH)
        y = self.R_nose * np.sin(TH) * np.sin(PH)
        z = self.R_nose * np.cos(TH)
        
        pts = np.stack([x, y, z], axis=-1)
        # Normals point radially outward 
        # for a sphere
        norms = pts / self.R_nose
        return pts.reshape(-1, 3), \
               norms.reshape(-1, 3)
    
    def surface_area(self):
        """Analytic area of spherical cap."""
        return (2 * np.pi * self.R_nose**2 
                * (1 - np.cos(self.theta_max)))
    
    def standoff_at(self, coil_radius, 
                    coil_offset_z):
        """Distance from each surface point
        to nearest coil location."""
        coil_pos = np.array([0, 0, coil_offset_z])
        dists = np.linalg.norm(
            self.points - coil_pos, axis=1
        )
        return dists - coil_radius
    