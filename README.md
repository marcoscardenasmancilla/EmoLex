# EmoLex (ES)

**Author / Autor**: Dr. Marcos H. Cárdenas Mancilla  
**E-mail**: marcoscardenasmancilla@gmail.com  
**Creation date / Fecha de creación**: 2025-07-25  
**License / Licencia**: AGPL V3  
**Copyright (c) 2025** Marcos H. Cárdenas Mancilla

---

## Description

This Python script implements a **Random Forest classifier** to automatically assign Spanish words to affective-semantic subgroups based on psycholinguistic and emotional variables. It integrates unsupervised clustering results (sub-clusters) with supervised classification to improve scalability and accuracy in lexical profiling.

### Key Features:

1. **Data Input**: Loads preprocessed data from `long_format_sub-clustering.csv`, containing affective ratings and sub-cluster labels.
2. **Data Cleaning**: Removes rows with missing values in predictor or target variables.
3. **Label Encoding**: Encodes string labels (if necessary) for classification.
4. **Training/Testing Split**: 80% training, 20% testing.
5. **Model Training**: Trains a `RandomForestClassifier` with 100 estimators.
6. **Evaluation**: Prints precision, recall, f1-score, and confusion matrix.
7. **Model Export**: Saves the trained model with a timestamp using `joblib`.
8. **Visualization**: Plots feature importances using `matplotlib` and `seaborn`.

### Predictors:

- `Valence_Mean`
- `Arousal_Mean`
- `Concreteness_Mean`
- `Emotionality`
- `Zipf_EsPal`
- `Balanced_Integration_Score`

### Objective:

To automate and enhance the classification of emotional words in Spanish by leveraging machine learning techniques that combine statistical and semantic cues.

---

## Descripción

Este script en Python implementa un **clasificador Random Forest** para asignar automáticamente palabras en español a subgrupos afectivo-semánticos, basándose en variables psicolingüísticas y emocionales. Integra resultados de clasificación no supervisada (sub-clústeres) con aprendizaje supervisado para mejorar la escalabilidad y precisión del perfilamiento léxico.

### Características principales:

1. **Entrada de datos**: Carga el archivo `long_format_sub-clustering.csv` con etiquetas de subagrupamiento y puntuaciones afectivas.
2. **Limpieza**: Elimina filas con valores faltantes en predictores o variable objetivo.
3. **Codificación de etiquetas**: Convierte etiquetas no numéricas en enteros si es necesario.
4. **División del conjunto**: 80% entrenamiento, 20% prueba.
5. **Entrenamiento del modelo**: Utiliza `RandomForestClassifier` con 100 árboles.
6. **Evaluación**: Imprime métricas de precisión, recall, f1-score y matriz de confusión.
7. **Exportación del modelo**: Guarda el modelo entrenado con `joblib` y timestamp.
8. **Visualización**: Grafica la importancia de los atributos predictivos con `matplotlib` y `seaborn`.

### Predictores utilizados:

- `Valence_Mean`
- `Arousal_Mean`
- `Concreteness_Mean`
- `Emotionality`
- `Zipf_EsPal`
- `Balanced_Integration_Score`

### Objetivo:

Automatizar y mejorar la clasificación de palabras emocionales en español utilizando técnicas de aprendizaje automático que combinan información estadística y semántica.

---

## References / Referencias

- Liesefeld, H. R., & Janczyk, M. (2019). Combining speed and accuracy to control for speed–accuracy trade-offs. *Behavior Research Methods, 51*(1), 40–60. https://doi.org/10.3758/s13428-018-1076-x  
- Pedregosa, F., Varoquaux, G., Gramfort, A., Michel, V., Thirion, B., Grisel, O., ... & Duchesnay, E. (2011). Scikit-learn: Machine learning in Python. *Journal of Machine Learning Research, 12*, 2825–2830.  
- Pérez-Sánchez, M. Á., Stadthagen-Gonzalez, H., Guasch, M., Hinojosa, J. A., Fraga, I., Marín, J., & Ferré, P. (2021). EmoPro: Emotional prototypicality for 1,286 Spanish words: Relationships with affective and psycholinguistic variables. *Behavior Research Methods, 53*(5), 1857–1875. https://doi.org/10.3758/s13428-020-01519-9
- Warriner, A. B., Kuperman, V., & Brysbaert, M. (2013). Norms of valence, arousal, and dominance for 13,915 English lemmas. *Behavior Research Methods, 45*, 1191–1207. https://doi.org/10.3758/s13428-012-0314-x  

---

## Cross-validation output log / Registro de salida de validación cruzada

<img width="755" height="373" alt="imagen" src="https://github.com/user-attachments/assets/d1206ac5-ab86-4daa-afdd-14e1b68b4eae" />
