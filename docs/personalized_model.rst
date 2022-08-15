Personalized signaling models
=============================

Breast cancer
-------------

.. raw:: html

    <a href="https://www.cell.com/iscience/fulltext/S2589-0042(22)00214-0">
    <img src="https://els-jbs-prod-cdn.jbs.elsevierhealth.com/cms/attachment/f4a8bd2a-0cae-4e62-ae32-b2d6a5d7575c/fx1_lrg.jpg" width="400px" hspace="30px" align="left">
    </a>

Patient heterogeneity precludes cancer treatment and drug development; hence, development of methods for finding prognostic markers for individual treatment is urgently required.

Using Pasmopy, we built a mechanistic model of ErbB receptor signaling network, trained with protein quantification data obtained from cultured cell lines, and performed *in silico* simulation of the pathway activities on breast cancer patients using `The Cancer Genome Atlas (TCGA) <https://www.cancer.gov/about-nci/organization/ccg/research/structural-genomics/tcga>`_ transcriptome datasets.
The temporal dynamics of Akt, extracellular signal-regulated kinase (ERK), and c-Myc in each patient were able to accurately predict the difference in prognosis and sensitivity to kinase inhibitors in triple-negative breast cancer (TNBC).

Protocols
---------

* **Paper:** https://doi.org/10.1016/j.xpro.2022.101619
* **Code:** https://github.com/pasmopy/breast_cancer