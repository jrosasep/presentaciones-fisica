# ¿Cómo encontrar caos partiendo desde el Modelo Estándar?

Repositorio de la presentación académica **¿Cómo encontrar caos partiendo desde
el Modelo Estándar?**, preparada en LaTeX Beamer por José Ignacio Rosas
Sepúlveda para el curso Sistemas Dinámicos e Introducción al Caos, en agosto de
2026.

La presentación recorre el paso desde una teoría clásica de campos hasta un
sistema hamiltoniano no lineal de pocos grados de libertad. El análisis utiliza
secciones de Poincaré para estudiar la aparición de dinámica regular, mixta y
caótica en reducciones relacionadas con teorías gauge y el sector Higgs.

## Contenido del repositorio

- `Como_encontrar_caos_partiendo_desde_el_Modelo_Estandar.tex`: fuente editable
  de la presentación.
- `Como_encontrar_caos_partiendo_desde_el_Modelo_Estandar.pdf`: versión
  compilada.
- `images/`: fotografías, logotipo y recursos visuales generales.
- `images/figures/`: figuras de dinámica no lineal y secciones de Poincaré.
- `images/papers/`: recortes de los artículos citados en la presentación.
- `scripts/`: programas utilizados para generar o reproducir las figuras
  numéricas.
- `data/`: series numéricas del doble péndulo exportadas como CSV.
- `FIGURE_PROVENANCE.md`: procedencia de cada figura y parámetros de
  reproducción.
- `requirements.txt`: dependencias de Python utilizadas por los programas.

El repositorio contiene únicamente la versión final y los recursos que ésta
necesita para compilar. No incluye versiones antiguas ni archivos auxiliares de
LaTeX.

## Compilación

El documento puede compilarse con una distribución moderna de LaTeX que incluya
Beamer, TikZ, SVG y los paquetes matemáticos habituales:

```bash
pdflatex Como_encontrar_caos_partiendo_desde_el_Modelo_Estandar.tex
pdflatex Como_encontrar_caos_partiendo_desde_el_Modelo_Estandar.tex
```

## Autor

José Ignacio Rosas Sepúlveda  
Universidad de Concepción

Las referencias científicas y atribuciones de las imágenes se indican en las
diapositivas correspondientes y en la sección final de referencias.

## Reproducción numérica

Las figuras numéricas no se publican como resultados aislados: se incluyen sus
programas, condiciones iniciales y dependencias. Consulta
`FIGURE_PROVENANCE.md` antes de reutilizarlas o interpretar sus resultados.
