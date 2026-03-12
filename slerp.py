import numpy as np
import matplotlib.pyplot as plt
from stl import mesh
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import KDTree
# =========================
# 3. Visualization
# =========================

# fig = plt.figure()
# ax = fig.add_subplot(projection='3d')

# # Draw mesh
# for triangle in triangles:
#     x = [triangle[0][0], triangle[1][0], triangle[2][0], triangle[0][0]]
#     y = [triangle[0][1], triangle[1][1], triangle[2][1], triangle[0][1]]
#     z = [triangle[0][2], triangle[1][2], triangle[2][2], triangle[0][2]]
#     ax.plot(x, y, z)

# # Draw normals (only some to avoid visual clutter)
# for i in range(0, len(centers), max(1, len(centers)//100)):
#     ax.quiver(
#         centers[i][0],
#         centers[i][1],
#         centers[i][2],
#         normals[i][0],
#         normals[i][1],
#         normals[i][2],
#         length=0.1
#     )

# ax.set_xlabel("X")
# ax.set_ylabel("Y")
# ax.set_zlabel("Z")

# plt.show()


# Convert each normal n into a quaternion q that rotates u = (0,0,1) towards n
def quat_from_two_vectors(u, v):
    u = u / np.linalg.norm(u)
    v = v / np.linalg.norm(v)
    dot = np.dot(u, v)           # dot = 1: parallel vectors (same direction)
                                # dot = 0: perpendicular vectors
                                # dot = -1: opposite vectors (opposite directions)
    if np.isclose(dot, -1.0):
        # Opposite vectors, rotate 180 degrees around any axis perpendicular to u
        axis = np.cross(u, np.array([1.0, 0.0, 0.0]))
        if np.isclose(np.linalg.norm(axis), 0.0):  # If u is parallel to (1,0,0), use another axis
            axis = np.cross(u, np.array([0, 1.0, 0]))
        axis = axis / np.linalg.norm(axis)
        return np.array([0.0, axis[0], axis[1], axis[2]])  # 180-degree quaternion
        
    w = 1.0 + dot
    xyz = np.cross(u, v)
    q = np.concatenate(([w], xyz))
    return q / np.linalg.norm(q)

def slerp_quat(q0, q1, t):
    """
    q0, q1: np.array(4,) quaternions [w,x,y,z] (unit)
    t: scalar in [0,1]
    """
    # normalize just in case
    q0 = q0 / np.linalg.norm(q0)
    q1 = q1 / np.linalg.norm(q1)

    # dot product
    dot = np.dot(q0, q1)

    # If dot < 0, flip q1 to take the short path in S^3
    # (q and -q represent the same rotation)
    if dot < 0.0:
        q1 = -q1
        dot = -dot

    # Avoid numerical issues
    dot = np.clip(dot, -1.0, 1.0)

    # If very close, approximate with LERP (and renormalize)
    if np.isclose(dot, 1.0):
        result = (1.0 - t) * q0 + t * q1
        return result / np.linalg.norm(result)

    # Angle between quats
    theta_0 = np.arccos(dot)      # original angle
    theta = theta_0 * t           # scaled angle
    sin_theta_0 = np.sin(theta_0)
    sin_theta = np.sin(theta)

    s0 = np.sin(theta_0 - theta) / sin_theta_0
    s1 = sin_theta / sin_theta_0

    return (s0 * q0) + (s1 * q1)



def nearest_neighbor(centers, quats):
    """
    Sorts centers and their respective quaternions following
    the nearest neighbor strategy.
    """
    remaining_points = list(range(len(centers)))
    order = []
    
    # Start from the first point in the array
    current_index = 0
    order.append(current_index)
    remaining_points.remove(current_index)
    
    while remaining_points:
        current_point = centers[current_index]
        
        # Compute Euclidean distances to all remaining points
        # (Using numpy vectorization for speed)
        distances = np.linalg.norm(centers[remaining_points] - current_point, axis=1)
        
        # Find the index of the minimum distance
        next = np.argmin(distances)
        
        # Get the real index from the original array
        current_index = remaining_points[next]
        
        order.append(current_index)
        remaining_points.remove(current_index)
        
    return centers[order], quats[order]


def KDTree_sort(centers, quats):
    n_points = len(centers)
    visited = np.zeros(n_points, dtype=bool)
    order = []
    tree = KDTree(centers)

    actual_index = 0
    order.append(actual_index)
    visited[actual_index] = True
    
    for i in range(n_points - 1):
        actual_point = centers[actual_index]
        # Search for the nearest unvisited neighbor
        distances, index = tree.query(actual_point, k=min(20, n_points)) 
        founded=False
        for idx in index:
            if not visited[idx]:
                actual_index = idx
                order.append(actual_index)
                visited[actual_index] = True
                founded=True
                break
        # If no unvisited neighbor found (rare), just take the next unvisited one
        if not founded:
            _,all_index = tree.query(actual_point, k=n_points)
            for id in all_index:
                if not visited[id]:
                    actual_index = id
                    founded=True
                    break
        order.append(actual_index)
        visited[actual_index] = True

    return centers[order], quats[order]

def print_comparison(centers_orig, centers_sorted, n):
    print(f"{'--- BEFORE SORTING (First {n}) ---':^50}")
    print(f"{'Index':<10} | {'Coordinates (X, Y, Z)':<30} | {'Dist. to next'}")
    print("-" * 65)
    for i in range(n):
        dist = np.linalg.norm(centers_orig[i+1] - centers_orig[i]) if i < n-1 else 0
        print(f"{i:<10} | {str(np.round(centers_orig[i], 2)):<30} | {dist:.4f}")

    print("\n" + "="*65 + "\n")

    print(f"{'--- AFTER SORTING (First {n}) ---':^50}")
    print(f"{'Index':<10} | {'Coordinates (X, Y, Z)':<30} | {'Dist. to next'}")
    print("-" * 65)
    for i in range(n):
        dist = np.linalg.norm(centers_sorted[i+1] - centers_sorted[i]) if i < n-1 else 0
        print(f"{i:<10} | {str(np.round(centers_sorted[i], 2)):<30} | {dist:.4f}")



def interpolate_full_trajectory(centers_ord, quats_ord, factor_densidad=10):
    """
    Generates intermediate points between each triangle center.
    factor_densidad: how many extra points to create between each original pair.
    """
    puntos_interp = []
    quats_interp = []
    
    for i in range(len(centers_ord) - 1):
        p0, p1 = centers_ord[i], centers_ord[i+1]
        q0, q1 = quats_ord[i], quats_ord[i+1]
        
        # Create t from 0 to 1 (excluding 1 to avoid duplicating with the next pair)
        for t in np.linspace(0, 1, factor_densidad, endpoint=False):
            # 1. LERP for position (straight line between centers)
            pos = (1 - t) * p0 + t * p1
            
            # 2. SLERP for orientation (smooth spherical rotation)
            ori = slerp_quat(q0, q1, t)
            
            puntos_interp.append(pos)
            quats_interp.append(ori)
            
    return np.array(puntos_interp), np.array(quats_interp)


# This function allows us to visualize the interpolated quaternion in matplotlib
def rotate_vector_by_quat(v, q):
    """ Rotation of a vector v by a quaternion q [w, x, y, z] """
    w, x, y, z = q
    q_vec = np.array([x, y, z])
    # Simplified rotation formula
    return v + 2.0 * np.cross(q_vec, np.cross(q_vec, v) + w * v)



def visualize_final_result(pos_interp, quat_interp):
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # 1. Plot the continuous trajectory line (Blue)
    ax.plot(pos_interp[:, 0], pos_interp[:, 1], pos_interp[:, 2], 
            label='Interpolated Trajectory (LERP+SLERP)', color='royalblue', lw=2, alpha=0.8)
    
    # 2. Plot the interpolated normals (Red)
    # Show 1 every 5 to avoid cluttering the view
    paso_flechas = 5
    for i in range(0, len(pos_interp), paso_flechas):
        p = pos_interp[i]
        q = quat_interp[i]
        
        # The vector (0,0,1) is our tool originally pointing "up"
        # We rotate it according to the computed quaternion
        n_dir = rotate_vector_by_quat(np.array([0, 0, 0.8]), q) # 0.8 is the visual length
        
        ax.quiver(p[0], p[1], p[2], 
                  n_dir[0], n_dir[1], n_dir[2], 
                  color='red', alpha=0.5, arrow_length_ratio=0.2)

    ax.set_title("Trajectory Simulation: Orientation via SLERP", fontsize=14)
    ax.set_xlabel("X (mm)")
    ax.set_ylabel("Y (mm)")
    ax.set_zlabel("Z (mm)")
    ax.legend()
    
    # Adjust view so the ring looks good
    ax.view_init(elev=30, azim=45)
    plt.show()


#---MAIN---

# Load STL
model = mesh.Mesh.from_file('ring.stl')
# Extract triangles
triangles = model.vectors
print(f"Number of triangles: {len(triangles)}")

# Compute normals manually
normals = []
centers = []  # We need the centers to better visualize the normals


for tri in triangles:
    v1, v2, v3 = tri
    # Triangle edge vectors
    edge1 = v2 - v1
    edge2 = v3 - v1
    # Compute normal with cross product
    normal = np.cross(edge1, edge2)
    # Normalize to get only direction, not magnitude
    normal = normal / np.linalg.norm(normal)
    # Triangle center
    center = (v1 + v2 + v3) / 3

    normals.append(normal)
    centers.append(center)

normals = np.array(normals)
centers = np.array(centers)

print("First 5 normals:")
print(normals[:5])

# Convert each normal into a quaternion that rotates (0,0,1) towards that normal
u_ref = np.array([0.0, 0.0, 1.0])
quaternions = np.array([quat_from_two_vectors(u_ref, n) for n in normals])
print(f"Generated quaternions: {quaternions.shape}")  # (N, 4)

# Sort centers and quaternions using KDTree
print("Sorting centers and quaternions with KDTree...")
centers_sorted, quats_sorted = KDTree_sort(centers, quaternions)
print("Centers sorted")
# Generate intermediate points with LERP for positions and SLERP for orientations
print("Generating interpolated trajectory...")
interpolated_points, interpolated_quats = interpolate_full_trajectory(centers_sorted, quats_sorted, factor_densidad=20)
# Visualize the final result
visualize_final_result(interpolated_points, interpolated_quats)

# Compare the first 10 points
#print_comparison(centers, centers_sorted, 20)