import numpy as np


class Quaternion:
    
    def __init__(self, w=1.0, x=0.0, y=0.0, z=0.0):
        """
        Initialize quaternion as w + xi + yj + zk
        Default is identity (no rotation)
        """
        self.q = np.array([w, x, y, z], dtype=np.float64)
        self._normalize()
    
    def _normalize(self):
        """Normalize to unit quaternion"""
        norm = np.linalg.norm(self.q)
        if norm > 1e-10:
            self.q /= norm
        else:
            self.q = np.array([1.0, 0.0, 0.0, 0.0])
    
    @classmethod
    def from_axis_angle(cls, axis, angle):
        """
        Create quaternion from axis-angle representation
        
        Parameters
        ----------
        axis : array-like (3,)
            Rotation axis vector (will be normalized)
        angle : float
            Rotation angle in radians
        
        Returns
        -------
        Quaternion
        """
        axis = np.asarray(axis, dtype=np.float64)
        axis_norm = np.linalg.norm(axis)
        
        if axis_norm < 1e-10:
            return cls()  # Identity quaternion
        
        axis = axis / axis_norm
        half_angle = angle / 2.0
        w = np.cos(half_angle)
        xyz = axis * np.sin(half_angle)
        
        return cls(w, xyz[0], xyz[1], xyz[2])
    
    @classmethod
    def from_euler(cls, roll, pitch, yaw, sequence='xyz'):
        """
        Create quaternion from Euler angles (for user convenience)
        
        Parameters
        ----------
        roll, pitch, yaw : float
            Euler angles in radians
        sequence : str
            Rotation sequence (default 'xyz')
        
        Returns
        -------
        Quaternion
        
        Note: While Euler angles have gimbal lock issues, the resulting
        quaternion does not. This is just a convenient input method.
        """
        # For 'xyz' sequence (most common)
        cr, cp, cy = np.cos([roll/2, pitch/2, yaw/2])
        sr, sp, sy = np.sin([roll/2, pitch/2, yaw/2])
        
        w = cr*cp*cy + sr*sp*sy
        x = sr*cp*cy - cr*sp*sy
        y = cr*sp*cy + sr*cp*sy
        z = cr*cp*sy - sr*sp*cy
        
        return cls(w, x, y, z)
    
    @classmethod
    def from_rotation_matrix(cls, R):
        """
        Create quaternion from 3x3 rotation matrix
        Uses Shepperd's method for numerical stability
        
        Parameters
        ----------
        R : array-like (3, 3)
            Rotation matrix
        
        Returns
        -------
        Quaternion
        """
        R = np.asarray(R, dtype=np.float64)
        
        trace = np.trace(R)
        
        if trace > 0:
            s = 0.5 / np.sqrt(trace + 1.0)
            w = 0.25 / s
            x = (R[2, 1] - R[1, 2]) * s
            y = (R[0, 2] - R[2, 0]) * s
            z = (R[1, 0] - R[0, 1]) * s
        elif R[0, 0] > R[1, 1] and R[0, 0] > R[2, 2]:
            s = 2.0 * np.sqrt(1.0 + R[0, 0] - R[1, 1] - R[2, 2])
            w = (R[2, 1] - R[1, 2]) / s
            x = 0.25 * s
            y = (R[0, 1] + R[1, 0]) / s
            z = (R[0, 2] + R[2, 0]) / s
        elif R[1, 1] > R[2, 2]:
            s = 2.0 * np.sqrt(1.0 + R[1, 1] - R[0, 0] - R[2, 2])
            w = (R[0, 2] - R[2, 0]) / s
            x = (R[0, 1] + R[1, 0]) / s
            y = 0.25 * s
            z = (R[1, 2] + R[2, 1]) / s
        else:
            s = 2.0 * np.sqrt(1.0 + R[2, 2] - R[0, 0] - R[1, 1])
            w = (R[1, 0] - R[0, 1]) / s
            x = (R[0, 2] + R[2, 0]) / s
            y = (R[1, 2] + R[2, 1]) / s
            z = 0.25 * s
        
        return cls(w, x, y, z)
    
    @classmethod
    def align_vectors(cls, vec_from, vec_to):
        """
        Create quaternion that rotates vec_from to vec_to
        
        Parameters
        ----------
        vec_from : array-like (3,)
            Initial vector
        vec_to : array-like (3,)
            Target vector
        
        Returns
        -------
        Quaternion
        """
        v1 = np.asarray(vec_from, dtype=np.float64)
        v2 = np.asarray(vec_to, dtype=np.float64)
        
        v1 = v1 / np.linalg.norm(v1)
        v2 = v2 / np.linalg.norm(v2)
        
        dot = np.dot(v1, v2)
        
        # Vectors are already aligned
        if dot > 0.99999:
            return cls()
        
        # Vectors are opposite
        if dot < -0.99999:
            # Find perpendicular axis
            axis = np.cross([1, 0, 0], v1)
            if np.linalg.norm(axis) < 1e-6:
                axis = np.cross([0, 1, 0], v1)
            return cls.from_axis_angle(axis, np.pi)
        
        # General case
        axis = np.cross(v1, v2)
        w = 1.0 + dot
        
        q = cls(w, axis[0], axis[1], axis[2])
        q._normalize()
        return q
    
    def to_rotation_matrix(self):
        """
        Convert quaternion to 3x3 rotation matrix
        
        Returns
        -------
        ndarray (3, 3)
            Rotation matrix
        """
        w, x, y, z = self.q
        
        R = np.array([
            [1-2*(y*y+z*z),   2*(x*y-w*z),   2*(x*z+w*y)],
            [  2*(x*y+w*z), 1-2*(x*x+z*z),   2*(y*z-w*x)],
            [  2*(x*z-w*y),   2*(y*z+w*x), 1-2*(x*x+y*y)]
        ], dtype=np.float64)
        
        return R
    
    def rotate_points(self, points):
        """
        Rotate point(s) using this quaternion
        
        Parameters
        ----------
        points : ndarray (..., 3)
            Points to rotate in Cartesian coordinates
        
        Returns
        -------
        ndarray (..., 3)
            Rotated points
        """
        points = np.asarray(points, dtype=np.float64)
        original_shape = points.shape
        
        # Reshape to (N, 3) for processing
        points_2d = points.reshape(-1, 3)
        
        # Use rotation matrix for efficiency with multiple points
        R = self.to_rotation_matrix()
        rotated = points_2d @ R.T
        
        return rotated.reshape(original_shape)
    
    def __mul__(self, other):
        """Quaternion multiplication (composition of rotations)"""
        if not isinstance(other, Quaternion):
            raise TypeError("Can only multiply Quaternion with Quaternion")
        
        w1, x1, y1, z1 = self.q
        w2, x2, y2, z2 = other.q
        
        w = w1*w2 - x1*x2 - y1*y2 - z1*z2
        x = w1*x2 + x1*w2 + y1*z2 - z1*y2
        y = w1*y2 - x1*z2 + y1*w2 + z1*x2
        z = w1*z2 + x1*y2 - y1*x2 + z1*w2
        
        return Quaternion(w, x, y, z)
    
    def conjugate(self):
        """Return conjugate (inverse for unit quaternions)"""
        w, x, y, z = self.q
        return Quaternion(w, -x, -y, -z)
    
    def __repr__(self):
        w, x, y, z = self.q
        return f"Quaternion(w={w:.4f}, x={x:.4f}, y={y:.4f}, z={z:.4f})"


class Transform:
    """
    Transformation class combining translation, rotation, and scaling
    Applied in order: scale -> rotate -> translate
    """
    
    def __init__(self, translation=None, rotation=None, scale=None):
        """
        Parameters
        ----------
        translation : array-like (3,), optional
            Translation vector (default: no translation)
        rotation : Quaternion, optional
            Rotation quaternion (default: no rotation)
        scale : float or array-like (3,), optional
            Uniform or per-axis scaling (default: no scaling)
        """
        self.translation = np.zeros(3, dtype=np.float64) if translation is None \
                          else np.asarray(translation, dtype=np.float64)
        
        self.rotation = Quaternion() if rotation is None else rotation
        
        if scale is None:
            self.scale = np.ones(3, dtype=np.float64)
        elif np.isscalar(scale):
            self.scale = np.full(3, scale, dtype=np.float64)
        else:
            self.scale = np.asarray(scale, dtype=np.float64)
    
    @classmethod
    def from_offset(cls, offset):
        """
        Create Transform from old-style Point offset (backward compatibility)
        
        Parameters
        ----------
        offset : Point or array-like (3,)
            Translation offset
        
        Returns
        -------
        Transform
        """
        if hasattr(offset, 'get_cart_crds'):
            translation = offset.get_cart_crds()
        else:
            translation = np.asarray(offset)
        
        return cls(translation=translation)
    
    def apply(self, points):
        """
        Apply transformation to point(s): scale -> rotate -> translate
        
        Parameters
        ----------
        points : ndarray (..., 3)
            Points in Cartesian coordinates
        
        Returns
        -------
        ndarray (..., 3)
            Transformed points
        """
        points = np.asarray(points, dtype=np.float64)
        
        # Apply scaling
        transformed = points * self.scale
        
        # Apply rotation
        transformed = self.rotation.rotate_points(transformed)
        
        # Apply translation
        transformed = transformed + self.translation
        
        return transformed
    
    def apply_inverse(self, points):
        """
        Apply inverse transformation: untranslate -> unrotate -> unscale
        
        Parameters
        ----------
        points : ndarray (..., 3)
            Transformed points
        
        Returns
        -------
        ndarray (..., 3)
            Original points
        """
        points = np.asarray(points, dtype=np.float64)
        
        # Reverse translation
        result = points - self.translation
        
        # Reverse rotation
        result = self.rotation.conjugate().rotate_points(result)
        
        # Reverse scaling
        result = result / self.scale
        
        return result
    
    def __repr__(self):
        return (f"Transform(\n"
                f"  translation={self.translation},\n"
                f"  rotation={self.rotation},\n"
                f"  scale={self.scale}\n"
                f")")


# Convenience functions for common transformations

def rotation_x(angle):
    """Rotation about X-axis"""
    return Quaternion.from_axis_angle([1, 0, 0], angle)


def rotation_y(angle):
    """Rotation about Y-axis"""
    return Quaternion.from_axis_angle([0, 1, 0], angle)


def rotation_z(angle):
    """Rotation about Z-axis"""
    return Quaternion.from_axis_angle([0, 0, 1], angle)


def rotation_about_axis(axis, angle):
    """Rotation about arbitrary axis"""
    return Quaternion.from_axis_angle(axis, angle)


def rotation_align(vec_from, vec_to):
    """Rotation that aligns one vector to another"""
    return Quaternion.align_vectors(vec_from, vec_to)
