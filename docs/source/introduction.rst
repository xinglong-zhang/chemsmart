##############
 Introduction
##############

.. image:: _static/chemsmart_logo.png
   :width: 400
   :align: center

Chemsmart is a python-based toolkit for the automatic creation of input and submission script files, the submission and
the analysis of quantum chemistry simulation jobs.

It uses the same submission command regardless of the queueing systems (SLURM, Torque or SLF) used by any High
Performance Computing (HPC) cluster.

Users can customize their own HPC server settings and project settings to run different jobs, without modifying the
codes in this package.

##########
 Citation
##########

If you use **CHEMSMART** in your work, please follow good scholarly practice and kindly cite our work

📄 **Paper URL**: https://arxiv.org/abs/2508.20042

-  ACS Style

      Zhang, X.; Tan, H.; Liu, J.; Li, Z.; Wang, L.; Chen, B. W. J. CHEMSMART: Chemistry Simulation and Modeling
      Automation Toolkit for High-Efficiency Computational Chemistry Workflows. *arXiv* **2025**, arXiv:2508.20042.
      https://doi.org/10.48550/arXiv.2508.20042.

-  BibTeX

      .. code:: console

         @misc{zhang2025chemsmartchemistrysimulationmodeling,
           title        = {CHEMSMART: Chemistry Simulation and Modeling Automation Toolkit for High-Efficiency Computational Chemistry Workflows},
           author       = {Xinglong Zhang and Huiwen Tan and Jingyi Liu and Zihan Li and Lewen Wang and Benjamin W. J. Chen},
           year         = {2025},
           eprint       = {2508.20042},
           archivePrefix= {arXiv},
           primaryClass = {physics.chem-ph},
           url          = {https://arxiv.org/abs/2508.20042},
           doi          = {10.48550/arXiv.2508.20042}
           }

In addition, if you use **ASE** Atoms object in **CHEMSMART**, please cite:

-  ACS Style

      Ask Hjorth Larsen et al The atomic simulation environment—a Python library for working with atoms. J. Phys.:
      Condens. Matter, 2017, 29, 273002.

-  BibTeX

      .. code:: console

         @article{Hjorth Larsen_2017,
            doi = {10.1088/1361-648X/aa680e},
            url = {https://dx.doi.org/10.1088/1361-648X/aa680e},
            year = {2017},
            month = {jun},
            publisher = {IOP Publishing},
            volume = {29},
            number = {27},
            pages = {273002},
            author = {Hjorth Larsen, Ask and Jørgen Mortensen, Jens and Blomqvist, Jakob and Castelli, Ivano E and Christensen, Rune and Dułak, Marcin and Friis, Jesper and Groves, Michael N and Hammer, Bjørk and Hargus, Cory and Hermes, Eric D and Jennings, Paul C and Bjerre Jensen, Peter and Kermode, James and Kitchin, John R and Leonhard Kolsbjerg, Esben and Kubal, Joseph and Kaasbjerg, Kristen and Lysgaard, Steen and Bergmann Maronsson, Jón and Maxson, Tristan and Olsen, Thomas and Pastewka, Lars and Peterson, Andrew and Rostgaard, Carsten and Schiøtz, Jakob and Schütt, Ole and Strange, Mikkel and Thygesen, Kristian S and Vegge, Tejs and Vilhelmsen, Lasse and Walter, Michael and Zeng, Zhenhua and Jacobsen, Karsten W},
            title = {The atomic simulation environment—a Python library for working with atoms},
            journal = {Journal of Physics: Condensed Matter},
            abstract = {The atomic simulation environment (ASE) is a software package written in the Python programming language with the aim of setting up, steering, and analyzing atomistic simulations. In ASE, tasks are fully scripted in Python. The powerful syntax of Python combined with the NumPy array library make it possible to perform very complex simulation tasks. For example, a sequence of calculations may be performed with the use of a simple ‘for-loop’ construction. Calculations of energy, forces, stresses and other quantities are performed through interfaces to many external electronic structure codes or force fields using a uniform interface. On top of this calculator interface, ASE provides modules for performing many standard simulation tasks such as structure optimization, molecular dynamics, handling of constraints and performing nudged elastic band calculations.}
            }

If you use RDKit functionalities in **CHEMSMART**, please cite:

-  ACS Style

      ARDKit: Open-source cheminformatics. https://www.rdkit.org

-  BibTeX

      .. code:: console

         @article{Landrum2016RDKit2016_09_4,
            added-at = {2017-04-11T06:11:47.000+0200},
            author = {Landrum, Greg},
            biburl = {https://www.bibsonomy.org/bibtex/28d01fceeccd6bf2486e47d7c4207b108/salotz},
            description = {Release 2016_09_4 (Q3 2016) Release · rdkit/rdkit},
            interhash = {ee9a4ddeff3121aa622cf35709fa6e21},
            intrahash = {8d01fceeccd6bf2486e47d7c4207b108},
            keywords = {chemoinformatics drug-design pharmacophores software},
            timestamp = {2017-04-11T06:11:47.000+0200},
            title = {RDKit: Open-Source Cheminformatics Software},
            url = {https://github.com/rdkit/rdkit/releases/tag/Release_2016_09_4},
            year = 2016
            }

Our package has minimal dependencies on **pymatgen**, but if you convert **CHEMSMART** molecule into pymatgen
**AseAtomsAdaptor**, please cite:

-  ACS Style

      A. Jain, S.P. Ong, G. Hautier, W. Chen, W.D. Richards, S. Dacek, S. Cholia, D. Gunter, D. Skinner, G. Ceder, K.A.
      Persson The Materials Project: A materials genome approach to accelerating materials innovation. *APL Materials*,
      2013, 1(1), 011002.

-  BibTeX
      .. code:: console

         @article{Jain2013,
            author = {Jain, Anubhav and Ong, Shyue Ping and Hautier, Geoffroy and Chen, Wei and Richards, William Davidson and Dacek, Stephen and Cholia, Shreyas and Gunter, Dan and Skinner, David and Ceder, Gerbrand and Persson, Kristin a.},
            doi = {10.1063/1.4812323},
            issn = {2166532X},
            journal = {APL Materials},
            number = {1},
            pages = {011002},
            title = {{The Materials Project: A materials genome approach to accelerating materials innovation}},
            url = {http://link.aip.org/link/AMPADS/v1/i1/p011002/s1\&Agg=doi},
            volume = {1},
            year = {2013}
            }

If you use **scikit-learn**, please cite

-  ACS Style

Pedregosa et al., Scikit-learn: Machine Learning in Python, *J. Mach. Learn. Res* 2011, 12, 2825-2830.

-  BibTeX
      .. code:: console

         @article{scikit-learn,
            title={Scikit-learn: Machine Learning in {P}ython},
            author={Pedregosa, F. and Varoquaux, G. and Gramfort, A. and Michel, V.
                  and Thirion, B. and Grisel, O. and Blondel, M. and Prettenhofer, P.
                  and Weiss, R. and Dubourg, V. and Vanderplas, J. and Passos, A. and
                  Cournapeau, D. and Brucher, M. and Perrot, M. and Duchesnay, E.},
            journal={Journal of Machine Learning Research},
            volume={12},
            pages={2825--2830},
            year={2011}
            }

**Please also cite other relavant software (e.g., Gaussian, ORCA, NCIPLOT, PyMOL) and DFT functionals and basis sets you
use in your research accordingly.**
