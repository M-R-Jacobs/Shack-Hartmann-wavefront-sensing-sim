# Shack-Hartmann Wavefront Sensor Simulator

Shack-Hartmann wavefront sensor (SHWFS)
Closed adaptive optics (AO) loop: known aberrated wavefront goes in, a sensor measures local slopes, reconstructor estimates the wavefront, and a deformable mirror corrects it. Every stage is validated against true wavefront data.

---

## 1. Wavefronts

A perfect point source emitting isotropically arrives at distant telescopes as a plane wave. Every point across the aperture has the same phase. Complications like atmospheric turbulence, misaligned optics, thermal distortion, and figure errors on mirrors all delay some parts of the beam more than others. The wavefront error, here called $W(x,y)$ is the optical path difference (OPD) from a theoretical reference surface (units of length). The corresponding phase error is:

$$\phi(x,y) = \frac{2\pi}{\lambda} W(x,y)$$

The phase error spreads light out of the diffraction-limited core of the point spread function (PSF). For small aberrations, the peak intensity relative to a perfect system (the Strehl ratio) follows the Maréchal approximation:

$$S \approx e^{-\sigma_\phi^2}$$

where $\sigma_\phi$ is the RMS phase error over the pupil in radians. A system is conventionally "diffraction limited" at $S \geq 0.8$, which corresponds to an RMS wavefront error of about $\lambda/14$. Adaptive optics exists to drive $\sigma_\phi$ down in real time. This is the same quantity an FTIR measures in a different form - a Michelson interferometer records intensity as a function of OPD between two arms. A wavefront sensor measures OPD as a function of position across one beam.

## 2. Zernike polynomials

These describe a wavefront. Over a circular pupil of radius $R$, any wavefront can be expanded in Zernike polynomials, which are orthonormal on the unit disk. In normalized polar coordinates $\rho = r/R$ and $\theta$:

$$W(\rho,\theta) = \sum_j a_j Z_j(\rho,\theta)$$

Each polynomial has a radial order $n$ and azimuthal frequency $m$ (with $n - |m|$ even and $|m| \leq n$):

$$Z_n^m(\rho,\theta) = N_n^m R_n^{|m|}(\rho) \cos(m\theta) \quad (m \geq 0)$$

$$Z_n^{m}(\rho,\theta) = N_n^m R_n^{|m|}(\rho) \sin(|m|\theta) \quad (m < 0)$$

with radial polynomial:

$$R_n^{m}(\rho) = \sum_{k=0}^{(n-m)/2} \frac{(-1)^k (n-k)!}{k! \, ((n+m)/2 - k)! \, ((n-m)/2 - k)!} \rho^{n-2k}$$

and normalization $N_n^m = \sqrt{2(n+1)}$ for $m \neq 0$, $N_n^0 = \sqrt{n+1}$. This sim uses Noll's single-index ordering $j$, under which the orthonormality condition is

$$\frac{1}{\pi} \int_0^{2\pi} \int_0^1 Z_i Z_j \, \rho \, d\rho \, d\theta = \delta_{ij}$$

The first few modes are the classical aberrations: $Z_1$ piston, $Z_2$ and $Z_3$ tip/tilt, $Z_4$ defocus, $Z_5$ and $Z_6$ astigmatism, $Z_7$ and $Z_8$ coma, $Z_{11}$ spherical. Because of the normalization, the RMS wavefront error is simply the root-sum-square of the coefficients, excluding piston (which has no effect on image quality):

$$\sigma_W^2 = \sum_{j \geq 2} a_j^2$$

This is considered the main reason to work in a Zernike basis: the error budget comes directly from the coefficient vector.

## 3. Shack-Hartmann sensors

A Shack-Hartmann sensor places a lenslet array at a plane conjugate to the pupil, with a detector at the lenslets' focal plane. The array divides the pupil into subapertures of pitch $d$, and each lenslet of focal length $f$ forms its own image of the source.

If the wavefront across a subaperture is flat, the spot lands on the lenslet's optical axis. If the wavefront is tilted, the spot moves. Geometrically, the displacement is the focal length times the wavefront slope averaged over the subaperture $A_i$:

$$\Delta x_i = \frac{f}{A_i} \iint_{A_i} \frac{\partial W}{\partial x} \, dA \qquad \Delta y_i = \frac{f}{A_i} \iint_{A_i} \frac{\partial W}{\partial y} \, dA$$

Remember that the sensor never measures the wavefront directly. It measures a grid of local gradients, and the wavefront has to be reconstructed from them. More on this below.

Each spot is the far-field diffraction pattern of its subaperture, so the intensity on the detector is the squared magnitude of the Fourier transform of the subaperture field:

$$I_i(u,v) \propto \left| \mathcal{F}\{ P_i(x,y) \, e^{i\phi(x,y)} \} \right|^2$$

where $P_i$ is the subaperture's aperture function. The simulator computes spots this way with an FFT rather than just shifting a Gaussian, which captures spot blurring from higher-order aberrations within a subaperture. A diffraction-limited spot has a width of order $\lambda f / d$.

The spot position is estimated with a center-of-mass centroid over the detector pixels $p$ in each subaperture window:

$$x_c = \frac{\sum_p x_p I_p}{\sum_p I_p} \qquad y_c = \frac{\sum_p y_p I_p}{\sum_p I_p}$$

The measured slope is the centroid minus a reference centroid recorded with a flat wavefront, divided by $f$.

There are some design tradeoffs here:
- If a spot moves out of its subaperture window, it is lost or assigned to the wrong lenslet. The maximum measurable tilt is roughly $\theta_{max} \approx d/(2f)$, so a short focal length gives more range.
- A longer focal length moves the spot more per unit slope, so the same centroid noise corresponds to a smaller slope error. Range and sensitivity trade directly against each other through $f$.
- More subapertures resolve higher-order aberrations but split the available light more ways. Aberration structure finer than the subaperture pitch cannot be sensed and aliases into lower-order modes.
- Noise. With photon-limited light, centroid error scales roughly as the spot width over the square root of the photon count, $\sigma_c \sim (\lambda f/d)/\sqrt{N_{ph}}$. Detector read noise and background add to this, and they hurt most in faint subapertures at the pupil edge.

## 4. The measurement model

Stack all $x$ and $y$ slopes from $M$ valid subapertures into a vector $\mathbf{s}$ of length $2M$, and the first $J$ Zernike coefficients into $\mathbf{a}$. Because slope averaging is linear in $W$, the sensor is a linear operator plus noise:

$$\mathbf{s} = G \mathbf{a} + \mathbf{n}$$

Each column of the interaction matrix $G$ is the slope pattern produced by one unit Zernike mode. Each element is the subaperture-averaged derivative of that mode:

$$G_{(i,x),j} = \frac{1}{A_i} \iint_{A_i} \frac{\partial Z_j}{\partial x} \, dA$$

Piston is excluded since it produces zero slope everywhere, so no slope sensor can see it.

## 5. Wavefront reconstruction

Reconstruction inverts the measurement model. With more slope measurements than modes ($2M > J$), the least-squares estimate is

$$\hat{\mathbf{a}} = (G^T G)^{-1} G^T \mathbf{s} = G^{+} \mathbf{s}$$

where $G^{+}$ is the Moore-Penrose pseudoinverse, computed by singular value decomposition (SVD). This is modal reconstruction, solving for a set of global modes. The alternative is zonal reconstruction, which solves for the wavefront height at grid points directly from neighbouring slope differences. Modal reconstruction is used here because the Zernike coefficients feed straight into the error budget.

The SVD also shows which modes the geometry can and cannot see. Small singular values belong to modes that produce almost no slope signal. Inverting them amplifies noise, so the reconstructor truncates or regularizes them.

Noise propagation needs to be accounted for. If slope noise is white with variance $\sigma_s^2$, the covariance of the reconstructed coefficients is:

$$C_{\hat{a}} = \sigma_s^2 \, (G^T G)^{-1}$$

Which is a predicted error, derived from the sensor geometry alone, and is the basis of the consistency test in Section 7.

## 6. Closing the loop: adaptive optics

A deformable mirror (DM) applies a correcting surface. Each actuator $k$ deforms the mirror with an influence function, modeled here as a Gaussian with inter-actuator coupling $\omega$ at actuator pitch $p$:

$$I_k(x,y) = \exp(\ln(\omega) \, r_k^2 / p^2)$$

where $r_k$ is the distance from actuator $k$. With actuator commands $\mathbf{c}$, the mirror surface is $W_{DM} = \sum_k c_k I_k$, and the light reflects off it, so the wavefront correction is twice the surface.

The DM interaction matrix $D$ is measured the way a real system is calibrated: poke each actuator by a known amount, record the slope vector, and store it as a column. The command matrix is the truncated pseudoinverse $D^{+}$. Truncating small singular values removes unsensed patterns, such as the checkerboard "waffle" mode that a Shack-Hartmann geometry is blind to.

The loop runs as a discrete integrator. At each frame $k$, the sensor measures the residual wavefront $W_{res} = W_{in} - W_{corr}$, and the commands update by

$$\mathbf{c}_{k+1} = \mathbf{c}_k - g \, D^{+} \mathbf{s}_k$$

with loop gain $g$. With no latency the loop is stable for $0 < g < 2$. With a one-frame delay between measurement and correction, which every real system has, the characteristic equation becomes $z^2 - z + g = 0$ and stability requires $g < 1$. Low gain rejects noise but responds slowly; high gain tracks quickly but amplifies noise and approaches instability. Residual RMS against gain is the central performance plot of the loop.

## 7. Validation

The simulator keeps the truth model strictly separate from the sensor and reconstructor. The truth wavefront is generated, the sensor only sees noisy pixels, and the reconstructor only sees slopes. Results are scored against the truth. Inject single Zernike modes of known amplitude and confirm $\hat{a}_j \approx a_j$ with negligible cross-talk into other modes, then measure where cross-talk grows (aliasing from modes beyond $J$). Over many Monte Carlo trials, the normalized estimation error:

$$\epsilon = (\hat{\mathbf{a}} - \mathbf{a})^T C_{\hat{a}}^{-1} (\hat{\mathbf{a}} - \mathbf{a})$$

should follow a $\chi^2$ distribution with $J$ degrees of freedom. If it does, the predicted covariance (from Section 5) matches reality. This is the same consistency logic as the NEES test for Kalman filters!

For a Kolmogorov phase screen with Fried parameter $r_0$ over aperture diameter $D$, the piston-removed phase variance is computed analytically (Noll, 1976):

$$\sigma_\phi^2 \approx 1.0299 \, (D/r_0)^{5/3} \ \text{rad}^2$$

along with the variance remaining after removing the first $J$ modes. Simulated phase screens are checked against these values before being used as loop inputs.

Residual RMS and Strehl ratio are measured against loop gain, photon count, and number of corrected modes, and compared with the Maréchal estimate.

## References for theory/data

- R. J. Noll, "Zernike polynomials and atmospheric turbulence," *J. Opt. Soc. Am.* 66, 207–211 (1976).
- B. C. Platt and R. Shack, "History and principles of Shack-Hartmann wavefront sensing," *J. Refract. Surg.* 17, S573–S577 (2001).
- J. W. Hardy, *Adaptive Optics for Astronomical Telescopes*, Oxford University Press (1998).
- R. K. Tyson, *Principles of Adaptive Optics*, CRC Press.
- J. W. Goodman, *Introduction to Fourier Optics*, for the subaperture diffraction model.
