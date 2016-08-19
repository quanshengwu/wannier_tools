# Wannier_tools

[**More examples in Wiki for Wannier_tools**](https://github.com/quanshengwu/wannier_tools/wiki)

![](https://lh3.googleusercontent.com/-NGkPcF7iUDY/Vy-34BbICBI/AAAAAAAAASY/e2YiWSnQJD4jpHh-kDWceThf2jKKSGAxwCCo/s526/wannier_tools-logo-purple.jpg)

**Authorship**

Written by QuanSheng Wu in Fortran 90 (wuquansheng@gmail.com  wuq@phys.ethz.ch)

Copyright (c) 2016 QuanSheng Wu and ShengNan Zhang. All rights reserved.

**Pull down the package**

git clone https://github.com/quanshengwu/wannier_tools.git

**Brief introductions**

Use tight binding model to get the surface states of slab systems
or edge states of nanowire systems or just bulk bands. Especially
usefull for topological novel systems, including topological insulator, Dirac semimetal, Weyl semimetal, nodal line systems, nodal chain systems, triple point systems and unknown topological systems.

With given tight binding model well written in the format of wannier90_hr.dat in [Wannier90](http://wannier.org), 

0. There are several examples in the examples folder, including Bi2Se3, WTe2, IrF4. The user guide is in the doc folder. Please read the user guide first. 

1. Identify topological class by calculating the Wilson loop (Wannier charge center)

2. Get surface state spectrum.

3. Identify Weyl points and nodal-line structure. 

**License and agreement**

If you use our code, please cite this website  “ The surface spectrums are calculated by the software package Wannier_tools~\cite{wannier_tools}, which is based on the iterative Green’s function~\cite{Sancho1985}.” 


Reference 

wannier_tools:  Q.S.Wu, S.N.Zhang, https://github.com/quanshengwu/wannier_tools

Sancho1985: [Highly convergent schemes for the calculation of bulk and surface Green functions, M P Lopez Sancho, J M Lopez Sancho, J M L Sancho and J Rubio, J.Phys.F.Met.Phys.15(1985)851-858](http://iopscience.iop.org/article/10.1088/0305-4608/15/4/009/meta;jsessionid=A349A81FE38B2B55DB42032F6792B275.c1)
 
If you have good ideas to improve this code, do not hesitate to contact me. Your contribution will be recorded.

**Publications**
* Heavy Weyl fermion state in CeRu4Sn6, Yuanfeng Xu, Changming Yue, Hongming Weng, Xi Dai, [arXiv:1608.04602](http://arxiv.org/abs/1608.04602) 
* Triple Point Topological Metals Ziming Zhu, Georg W. Winkler, **QuanSheng Wu**, Ju Li, Alexey A. Soluyanov  [arXiv:1605.04653](http://arxiv.org/abs/1605.04653) Phys. Rev. X 6, 031003 – Published 7 July 2016.
* Fermi arcs and their topological character in the candidate type-II Weyl semimetal MoTe2, A. Tamai, **Q. S. Wu**, I. Cucchi, F. Y. Bruno, S. Ricco, T.K. Kim, M. Hoesch, C. Barreteau, E. Giannini, C. Bernard, A. A. Soluyanov, F. Baumberger  [arXiv:1604.08228](http://arxiv.org/abs/1604.08228) 
* Nodal chain metals, Tomáš Bzdušek, **QuanSheng Wu**, Andreas Rüegg, Manfred Sigrist, Alexey A. Soluyanov [arXiv:1604.03112, 2016 ](https://arxiv.org/abs/1604.03112) Shown in Nature soon.
* Surface states and bulk electronic structure in the candidate type-II Weyl semimetal WTe2, F. Y. Bruno, A. Tamai, **Q. S. Wu**, I. Cucchi, C. Barreteau, A. de la Torre, S. McKeown Walker, S. Riccò, Z. Wang, T. K. Kim, M. Hoesch, M. Shi, N. C. Plumb, E. Giannini, A. A. Soluyanov, F. Baumberger [arXiv:1604.02411, 2016](https://arxiv.org/abs/1604.02411)
* Topological Phases in InAs1−xSbx: From Novel Topological Semimetal to Majorana Wire, Georg W. Winkler, **QuanSheng Wu**, Matthias Troyer, Peter Krogstrup, Alexey A. Soluyanov [arxiv:1602.07001, 2016](https://arxiv.org/abs/1602.07001)
* Type-II Weyl semimetals, Alexey A. Soluyanov,	Dominik Gresch,	Zhijun Wang,	**QuanSheng Wu**,	Matthias Troyer,	Xi Dai	& B. Andrei Bernevig, [Nature 527, 495–498 (26 November 2015)](http://www.nature.com/nature/journal/v527/n7579/full/nature15768.html) 


