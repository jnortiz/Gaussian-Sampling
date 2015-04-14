# PKE-PEKS
<br> University of Campinas | Institute of Computing </br>
<br> Laboratory of Security and Applied Cryptography (LASCA) </br>

<br> Author: Jheyne N. Ortiz </br>
<br> Advisor: Ricardo Dahab </br>
<br> Co-advisor: Diego F. Aranha </br>

This implementation is intended to have three C++ classes: HIBE, OT and PKE-PEKS. Different from Chen and Lin, we replace the HIBE system of Agrawal et al. by the HIBE construction of (Mochetti and Dahab, 2014) using ideal lattices in order to reduce the public key size.

So far, we are considering only the HIBE class, which has four main algorithms: SetUp, KeyDerive, Encrypt and Decrypt. The SetUp algorithm is already able to generate the system master keys, e.g. two basis for an orthogonal lattice and it uses the IdealTrapGen algorithm (Stehlé et al., 2009).

In focus, we have the KeyDerive algorithm that requires samples from a discrete Gaussian distribution. The discrete Ziggurat algorithm (Buchmann et al., 2014) was adopted in this work due to its simplicity and efficiency (about one million samples per second). It is a sampler by rejection that outputs integers statistically close to a discrete Gaussian distribution centered in zero.
