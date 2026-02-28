import numpy as np
import matplotlib.pyplot as plt

def slerp(v0, v1, t):
    # Producto punto
    dot = np.dot(v0, v1)
    
    # Evitar errores numéricos
    dot = np.clip(dot, -1.0, 1.0)

    # Ángulo entre vectores
    theta = np.arccos(dot)

    # Si el ángulo es muy pequeño, usar LERP
    if np.isclose(theta, 0):
        return (1 - t) * v0 + t * v1
    
    # Fórmula SLERP
    sin_theta = np.sin(theta)
    
    factor1 = np.sin((1 - t) * theta) / sin_theta
    factor2 = np.sin(t * theta) / sin_theta
    
    return factor1 * v0 + factor2 * v1


# Vectores inicial y final

v0 = np.array([5, 4, 1])
v1 = np.array([2, 1, 3])


# Normalizar vectores
v0 = v0 / np.linalg.norm(v0)
v1 = v1 / np.linalg.norm(v1)

print (f"v0 normalizado: {v0}")
print (f"v1 normalizado: {v1}")

# Generar trayectoria completa
ts = np.linspace(0, 1, 100)
trajectory = np.array([slerp(v0, v1, t) for t in ts])

# imprimir algunos puntos de la trayectoria
for t in np.linspace(0, 1, 5):
    point = slerp(v0, v1, t)
    print(f"t={t:.2f} -> {point}")

fig = plt.figure()
ax = fig.add_subplot(projection='3d')

# Trayectoria SLERP
ax.plot(trajectory[:,0], trajectory[:,1], trajectory[:,2])

# Vector inicial
ax.quiver(0, 0, 0, v0[0], v0[1], v0[2], color = 'red')

# Vector final
ax.quiver(0, 0, 0, v1[0], v1[1], v1[2], color = 'blue')

ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")

ax.set_xlim([-1, 1])
ax.set_ylim([-1, 1])
ax.set_zlim([-1, 1])

plt.show()


# =========================
# Imprimir algunos puntos
# =========================

for t in np.linspace(0, 1, 5):
    point = slerp(v0, v1, t)
    print(f"t={t:.2f} -> {point}")