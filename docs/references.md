---
title: References
nav_order: 12
description: "Bibliography for methods, materials, and software used in COSMO."
---

# References
{: .no_toc }

---

## Governing equations and compressible flow

**Euler, L.** (1757). "Principes généraux du mouvement des fluides."
*Mémoires de l'Académie Royale des Sciences et des Belles-Lettres de Berlin*,
**11**, 274–315.

**Anderson, J. D.** (2003).
*Modern Compressible Flow: With Historical Perspective*, 3rd ed.
McGraw-Hill. ISBN 978-0-07-242443-0.

**LeVeque, R. J.** (2002).
*Finite Volume Methods for Hyperbolic Problems*.
Cambridge University Press. ISBN 978-0-521-00924-9.

**Toro, E. F.** (2009).
*Riemann Solvers and Numerical Methods for Fluid Dynamics: A Practical
Introduction*, 3rd ed.
Springer. ISBN 978-3-540-25202-3.
[doi:10.1007/b79761](https://doi.org/10.1007/b79761)

---

## HLLC Riemann solver

**Toro, E. F., Spruce, M., and Speares, W.** (1994).
"Restoration of the contact surface in the HLL Riemann solver."
*Shock Waves*, **4**(1), 25–34.
[doi:10.1007/BF01414629](https://doi.org/10.1007/BF01414629)

**Davis, S. F.** (1988).
"Simplified second-order Godunov-type methods."
*SIAM Journal on Scientific and Statistical Computing*, **9**(3), 445–473.
[doi:10.1137/0909030](https://doi.org/10.1137/0909030)

**Harten, A., Lax, P. D., and van Leer, B.** (1983).
"On upstream differencing and Godunov-type schemes for hyperbolic
conservation laws."
*SIAM Review*, **25**(1), 35–61.
[doi:10.1137/1025002](https://doi.org/10.1137/1025002)

---

## MUSCL reconstruction and slope limiters

**van Leer, B.** (1979).
"Towards the ultimate conservative difference scheme. V. A second-order
sequel to Godunov's method."
*Journal of Computational Physics*, **32**(1), 101–136.
[doi:10.1016/0021-9991(79)90145-1](https://doi.org/10.1016/0021-9991(79)90145-1)

**Sweby, P. K.** (1984).
"High resolution schemes using flux limiters for hyperbolic conservation
laws."
*SIAM Journal on Numerical Analysis*, **21**(5), 995–1011.
[doi:10.1137/0721062](https://doi.org/10.1137/0721062)

**Harten, A.** (1983).
"High resolution schemes for hyperbolic conservation laws."
*Journal of Computational Physics*, **49**(3), 357–393.
[doi:10.1016/0021-9991(83)90136-5](https://doi.org/10.1016/0021-9991(83)90136-5)

---

## Equation of state — JWL

**Lee, E. L., Finger, M., and Collins, W.** (1973).
*JWL Equation of State Coefficients for High Explosives*
(UCID-16189). Lawrence Livermore National Laboratory.

**Dobratz, B. M., and Crawford, P. C.** (1985).
*LLNL Explosives Handbook: Properties of Chemical Explosives and Explosive
Simulants* (UCRL-52997, Rev. 2).
Lawrence Livermore National Laboratory.

**Cooper, P. W.** (1996).
*Explosives Engineering*.
Wiley-VCH. ISBN 978-0-471-18636-6.

---

## Detonation and programmed burn

**Fickett, W., and Davis, W. C.** (1979).
*Detonation*.
University of California Press. ISBN 978-0-520-03587-0.

**Mader, C. L.** (2007).
*Numerical Modeling of Explosives and Propellants*, 3rd ed.
CRC Press. ISBN 978-0-8493-7439-4.

**Tarver, C. M., and Hallquist, J. O.** (1981).
"Modeling short pulse duration shock initiation of solid explosives."
*Proceedings of the 7th Symposium on Detonation*, Annapolis, MD.

**Brode, H. L.** (1955).
"Numerical solutions of spherical blast waves."
*Journal of Applied Physics*, **26**(6), 766–775.
[doi:10.1063/1.1722085](https://doi.org/10.1063/1.1722085)

---

## Strong-stability-preserving time integration

**Shu, C.-W., and Osher, S.** (1988).
"Efficient implementation of essentially non-oscillatory shock-capturing
schemes."
*Journal of Computational Physics*, **77**(2), 439–471.
[doi:10.1016/0021-9991(88)90177-5](https://doi.org/10.1016/0021-9991(88)90177-5)

**Gottlieb, S., Shu, C.-W., and Tadmor, E.** (2001).
"Strong stability-preserving high-order time discretization methods."
*SIAM Review*, **43**(1), 89–112.
[doi:10.1137/S003614450036757X](https://doi.org/10.1137/S003614450036757X)

**Heun, K.** (1900).
"Neue Methoden zur approximativen Integration der Differentialgleichungen
einer unabhängigen Veränderlichen."
*Zeitschrift für Mathematik und Physik*, **45**, 23–38.

---

## Software

**Bezanson, J., Edelman, A., Karpinski, S., and Shah, V. B.** (2017).
"Julia: A fresh approach to numerical computing."
*SIAM Review*, **59**(1), 65–98.
[doi:10.1137/141000671](https://doi.org/10.1137/141000671)

**Byrne, S., Wilcox, L. C., and Churavy, V.** (2021).
"MPI.jl: Julia bindings for the Message Passing Interface."
*Proceedings of the JuliaCon Conferences*, **1**(1), 68.
[doi:10.21105/jcon.00068](https://doi.org/10.21105/jcon.00068)

**Polanco, J. I.** (2024).
"WriteVTK.jl — Julia package for writing VTK XML files."
[github.com/JuliaIO/WriteVTK.jl](https://github.com/JuliaIO/WriteVTK.jl)

---

## Blast and weapons effects

**Kingery, C. N., and Bulmash, G.** (1984).
*Airblast Parameters from TNT Spherical Air Burst and Hemispherical
Surface Burst* (ARBRL-TR-02555).
U.S. Army Ballistic Research Laboratory.

**UFC 3-340-02** (2008).
*Structures to Resist the Effects of Accidental Explosions*.
Unified Facilities Criteria, U.S. Department of Defense.

**Baker, W. E., Cox, P. A., Westine, P. S., Kulesz, J. J., and Strehlow, R. A.**
(1983).
*Explosion Hazards and Evaluation*.
Elsevier. ISBN 978-0-444-42094-3.
