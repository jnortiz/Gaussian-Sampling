# PKE-PEKS
<br> University of Campinas </br>
<br> Institute of Computing </br>
<br> Laboratory of Security and Applied Cryptography (LASCA) </br>

<br><b> Project name: PKE-PEKS</b> </br>
<br> Authors: Jheyne N. Ortiz </br>
<br> Advisor: Ricardo Dahab </br>
<br> Co-advisor: Diego F. Aranha </br>
<br> Brief description: Software implementation of the PKE-PEKS scheme proposed by (Chen and Lin, 2014) in C++ language using NTL. </br>

This implementation includes three C++ classes: HIBE, OT and PKE-PEKS. So far, our focus is the HIBE class and, more specifically, the Setup algorithm. The Setup algorithm creates a pair of master keys (public and private) and it used the IdealTrapGen algorithm proposed by Stehlé et al. in 2009 in order to build two bases for the same orthogonal lattice. Different from Chen and Lin, we replace the HIBE system of Agrawal et al. by the HIBE construction of Mochetti and Dahab using ideal lattices in order to reduce the public key size.
