## Presentation Notes:
### Slide 1
Title Slide.
### Slide 4:
Reduced order modeling is an unsupervised learning problem, where the high correlations between parater and/or responses of interest data streams, we exploit to find a reduced complexity representation, that involves fewer number of active variables. 
The premise is that the number of active variables (Intrinsic dimensionality) is significantly less than the nominal dimensionality and hence executing the reduced model is less taxing.
This is usually achieved by projection methods such as POD which is great for steady state problems, but less significant in time-dependent once since POD averages out the temporal behaviour through the snapshots.
In this work, Dynamic mode decomposition, is used as a soley data-driven approach to extract the temporal effect and decouple it from the spatial one, and thus suits the time dependednt applications.
### Slide 5-7:
Introduce LRA Problem.
### Slide 8:
DMD has the ability to provide the frequencies associated to each dynamic mode on the same plate. It can be looked at as a hybridization of Principle componant analysis (The discrete version of POD) responsible for extracting the spatial behaviour, and the Discrete Fourier Transform responsible for extracting the decay, growth or oscillatory frequencies.
### Slide 9: DMD Methodology
- Sequential dataset containing $m$ $n$-dimentional snapsots, is aggregated in two data matrices, $X_1$ the past data, and $X_2$ the future data, since each column in $X_2$ represents the corresponding state vector in $X_1$ after $\Delta t$.
- The algorithm baisically searches for the push forward linear operator ${\texfbf{A}}$ that maps data from $X_1$ tp $X_2$. 
### Slide 12
- This is achieved via minimizing the frobenous norm $\|textbf{X_2} - \textbf{AX}_1\|$
- It can be shown that this is similar to compiting the eigen values of the covariance of the data $\frac{1}{m-1}\textbf{XX^T}$ or directly computing the Singular value decomposition of the past data $\textbf{X_1 = U \Sigma V^H}$.
- Computing the low-rank approximation of $\textbf{A}$ is a lot less expensive and facilitates the computation of the original eigenpairs since $\textbf{\tilde{A}} = \textbf{U_{r}^HAU_r}$ is a similarity transformation, (then the eigenvalues are similar and the original eigenvectors are the projection of eigen vectors)
### Slide 13
- There are several methods to compute the amplitudes (mode contributions) such as computing them based on the first snapshot (representing the I.C.), or an optimized version introduced by Jovanovi'c, or a time dependent version based on the inner product with right (adjoint) DMD modes \cite{Wei}. In this work we only used the first two alternatives.
### Slide 14
- Unlike POD and since DMD modes are non-orthogonal/optimal, increasing the rank does not necessarily enhance accuracy.
- In transient problems some dynamics might be faster then others, some behavious might stabilize earlier, then others or maybe disappear before others, and hence the required rank is dependent on the time window. This is usually considered an undesirable feature in DMD, yet it was the real motivation for the idea of Partitioned DMD!.    
