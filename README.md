# Simulación del estómago

Simulación de la digestión gástrica en Python, construida **paso a paso**: primero cada proceso por separado y al final un modelo integrado. Cubre la química ácida, el flujo con paredes peristálticas, la rotura y erosión de las partículas de comida y el vaciado por el píloro. Cada fase escribe sus ecuaciones gobernantes, verifica la numérica y lleva la cuenta del error.

Todos los cuadernos están pensados para **Google Colab** y funcionan en CPU con las librerías que Colab trae instaladas.

![Ocho momentos de la digestión: ingesta, mezcla, trituración y vaciado](media/fase7_momentos_clave.png)

Película completa, 4 h de digestión: [`media/fase7_digestion_2D.mp4`](media/fase7_digestion_2D.mp4) (también en [GIF](media/fase7_digestion_2D.gif)).

## Fases

| Fase | Cuaderno | Contenido |
|---|---|---|
| 1 | [`fase1_quimica_estomago`](notebooks/fase1_quimica_estomago.ipynb) | Química 0D: secreción de ácido, activación autocatalítica de la pepsina, Michaelis–Menten dependiente del pH. Euler, Heun, RK4 y Radau, con extrapolación de Richardson y soluciones manufacturadas. |
| 1b | [`fase1b_validacion_carne_pan`](notebooks/fase1b_validacion_carne_pan.ipynb) | Validación con alimentos reales (carne frente a pan): tampón de varios pKa en equilibrio y vaciado por energía. Abre el registro de validación. |
| 2 | [`fase2_transporte_1D_2D`](notebooks/fase2_transporte_1D_2D.ipynb) | Advección–difusión–reacción en 1D y 2D con volúmenes finitos (MUSCL / van Leer, SSP-RK2). |
| 3 | [`fase3_fluido_LBM`](notebooks/fase3_fluido_LBM.ipynb) | Navier–Stokes 2D con Lattice Boltzmann D2Q9 escrito a mano en NumPy y paredes móviles peristálticas, validado con el bombeo de Shapiro. Versión opcional para GPU (JAX o CuPy). |
| 4 | [`fase4_particulas`](notebooks/fase4_particulas.ipynb) | Partículas: arrastre (Maxey–Riley, Schiller–Naumann), rotura, erosión y criba pilórica, con un modelo lagrangiano (Monte Carlo) y otro de balance de poblaciones. |
| 5 | [`fase5_visualizacion`](notebooks/fase5_visualizacion.ipynb) | Visualización: campos, `quiver`, `streamplot`, criterio Q y animaciones con Matplotlib (Plotly opcional). |
| 6 | [`fase6_modelo_integrado`](notebooks/fase6_modelo_integrado.ipynb) | Un único bucle de tiempo con flujo (promediado en fase desde el LBM), química local, partículas y píloro abierto. Incluye una tabla que compara cada versión del modelo con métricas clínicas, separando el error numérico, el paramétrico y el de modelo. |
| 7 | [`fase7_representacion_2D`](notebooks/fase7_representacion_2D.ipynb) | Representación anatómica 2D de la ingesta, los movimientos peristálticos y la dilución progresiva de la comida. |

Los cuadernos se pueden ejecutar en orden o por separado. Los que necesitan el flujo de la fase 3 (6 y 7) lo regeneran si no encuentran `campo_antro.npz`.

## Cómo ejecutarlo

**En Colab**: abre `https://colab.research.google.com/github/<usuario>/<repositorio>/blob/main/notebooks/<cuaderno>.ipynb` y ejecuta todas las celdas.

**En local**:
```bash
pip install -r requirements.txt
jupyter lab notebooks/
```
Para guardar las animaciones en MP4 hace falta `ffmpeg`. Sin él, se guardan en GIF.

## Limitaciones

Es un proyecto **docente y exploratorio**, no una herramienta clínica.
- Muchos parámetros son ilustrativos o vienen de rangos de la literatura. Cada cuaderno indica cuáles.
- La geometría es una losa 2D. La forma anatómica de la fase 7 es un dibujo sobre el que se proyectan los resultados del modelo.
- El código de GPU (JAX o CuPy) y de Plotly está protegido para que el cuaderno funcione sin esas librerías. Se ha probado menos que el de NumPy.

Las fuentes de las métricas de validación están citadas en el cuaderno 1b.
