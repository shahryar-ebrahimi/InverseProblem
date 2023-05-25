# InverseProblem
a matlab code to solve inverse problem for EEG and MEG Electrodes

A) Let's consider a spherical model with three concentric spheres of radii r2=8.5 cm, r1=8 cm, and r3=9 cm, respectively. The conductivities in accordance with the spherical model introduced in class are s1=0, s2=1, and s3=2/80. In this section, we deal with the process of filling the matrices, data simulation, and solving the inverse problem.

B) Suppose the electrode positions for MEG (Magnetoencephalography) and EEG (Electroencephalography) on the scalp are arranged as follows: 33 electrodes are placed on the head, where the first electrode is on the z-axis and the others are distributed as follows: on each of the 4 arcs, φi=i*(45°) for i=0,1,2,3,...,7, and there are four electrodes at angles θj=j*(22.5°) for j=1,2,3,4.

C) Let's assume the possible locations of the bipolar current dipoles on a sphere with radius r0=7 cm. There are 65 locations on the surface of this sphere, with one on the z-axis and the rest distributed as follows: on each of the 16 arcs, φi=i*(22.5°) for i=0,1,2,3,...,15, and there are four points at angles θj=j*(22.5°) for j=1,2,3,4.

-1) Calculate the direct lead field matrix for MEG and plot its column 33 in your report. Each column of this matrix represents something. What does each row represent?

-2) Assuming at a specific point r0=7 cm, q0=45°, and j0=180°, there exists a bipolar current dipole with a magnitude of 1. Calculate the magnetic field magnitude radially at all the electrodes. Represent this magnitude as a color map in your report, or you can use a two-dimensional image with theta and phi axes.

-3) To measure the obtained magnitudes in the previous section, we want to solve the inverse problem. Consider the electrode positions as the 65 locations described in part C. Solve the inverse problem using the minimum norm method. Calculate the reconstructed magnitudes of the bipolar dipoles at all the points and represent them as an image or mapping. Are the results close to the actual values?

-4) Model the signal disturbance as Gaussian noise. Generate 33 random numbers with Gaussian distribution, with mean 0 and standard deviation sigma. Add these values to the measured amplitudes on all the electrodes. Choose sigma such that the signal-to-noise ratio (SNR) is Log 20 50 A_s. Solve the inverse problem again using the minimum norm method and plot the reconstructed bipolar dipoles.

-5) Apply the regularized minimum norm method to the given data in the previous section. Report the results for different regularization parameter values. What role does this parameter play in the solution?

-6) Consider the formula for calculating the potential in the spherical model for EEG. Derive the equation in matrix/vector form for a single source and a bipolar dipole. Separate the matrix components and write them in your report.

-7) For the same simplified case as before, if the bipolar location is specified as (r0, 0, 0) and (q, j), modify the equation by performing the necessary rotations in the coordinate axes. Note that since the bipolar location can be in any octant of the space, the rotations may differ.

Using the equation mentioned in the previous section, if a bipolar electrode with the specified magnitude and location in section 2 exists in the head volume, the calculation of the potential on all electrodes on the scalp is desired. Report the mapping (image) of the resulting potential.

In section 7, calculate the direct interaction matrix for EEG for 65 bipolar locations and 33 electrodes using the results from the previous section. Plot the 23rd row of this matrix.

In section 9, using the calculated data from section 8 and the matrix from section 7, solve the inverse problem of EEG using the minimum norm method. Plot and report the amplitudes of the bipolar sources in the mapping.

In section 10, to consider the time element in the simulation, modulate the assumed bipolar amplitudes from section 2 in a triangular waveform with an amplitude of one and a phase of one second. Assume a sampling time of 2 milliseconds. Another bipolar source is placed at r0 = 7cm, q0 = 45 degrees, and j0 = 90 degrees in the first half of the second, and its amplitude is [0, 0, 1] in the next half in the opposite direction (180 degrees). Add noise to the signals at all time points, assuming the noise variance is constant. Finally, solve the inverse problem using the MUSIC method, report the obtained amplitudes of the bipolar sources at ./25 seconds in the mapping.
