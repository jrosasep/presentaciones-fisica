# Procedencia y reproducción de las figuras

Este documento distingue los cálculos numéricos propios de los recursos
gráficos provenientes de publicaciones o sitios externos. Las atribuciones y
referencias científicas completas se conservan en la presentación.

## Figuras numéricas reproducibles

| Archivo | Sistema o procedimiento | Programa |
|---|---|---|
| `images/figures/doble_pendulo_region_regular.*` | Doble péndulo, masas y longitudes unitarias; control `(theta1, theta2)=(0.25, 0.35)` | `scripts/generar_doble_pendulo.py` |
| `images/figures/doble_pendulo_region_caotica.*` | Mismo sistema; control `(theta1, theta2)=(2.20, 0.40)` | `scripts/generar_doble_pendulo.py` |
| `images/figures/poincare_intro_real.png` | Hamiltoniano cuártico reducido, comparación de energías | `scripts/generar_figuras_caos.py` |
| `images/figures/salasnich_poincare_6E_contraste.png` | Secciones de Poincaré del Hamiltoniano cuártico para seis energías | `scripts/generar_figuras_caos.py` |
| `images/figures/electroweak_poincare_contraste.png` | Reducción homogénea del sector SU(2)-Higgs con `g=v_EW=1` | `scripts/generar_figuras_caos.py` |
| `images/figures/canfora_poincare_contraste.png` | Sistema efectivo de Canfora–Grandi–Oyarzo–Oliva | `scripts/generar_figuras_caos.py` |

El programa general también genera `lyapunov_exponente.png` y
`reducciones_adicionales_poincare.png`, aunque esas imágenes no forman parte de
la versión final de la presentación.

## Figura numérica adaptada de una fuente publicada

`dyons_poincare_2panel_contraste.png` es un recorte y ajuste de contraste de
`dyons_poincare_original.png`. El programa `generar_figuras_caos.py` realiza
solamente ese posprocesamiento: no inventa ni vuelve a calcular los puntos. La
fuente científica se identifica en la diapositiva correspondiente.

## Recursos externos

Las fotografías, logotipos, recortes de artículos y figuras experimentales que
no aparecen en la tabla anterior no son resultados numéricos propios. Se
mantienen como recursos de apoyo de la presentación y deben conservar su
atribución original.

## Ejecución

Desde esta carpeta:

```bash
python -m venv .venv
python -m pip install -r requirements.txt
python scripts/generar_doble_pendulo.py
python scripts/generar_figuras_caos.py
```

Los programas escriben las figuras en `images/figures/`. El generador del doble
péndulo también conserva las series numéricas en `data/` como archivos CSV.
