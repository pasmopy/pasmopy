Personalized signaling models
=============================

Breast cancer
-------------

.. raw:: html

    <a href="https://www.cell.com/iscience/fulltext/S2589-0042(22)00214-0">
    <img src="https://els-jbs-prod-cdn.jbs.elsevierhealth.com/cms/attachment/f4a8bd2a-0cae-4e62-ae32-b2d6a5d7575c/fx1_lrg.jpg" width="400px" hspace="30px" align="left">
    </a>

The temporal activation dynamics of signaling pathways play important roles for cell fate decisions. Therefore, we hypothesized that signaling dynamics can be further utilized as prognostic biomarkers for human diseases. However, the majority of available data obtained from patients represent static snapshots taken at a single point in time, and not time-resolved dynamics. To overcome this problem, we developed a Patient-Specific Modeling in Python (Pasmopy), an open-source package for the development of dynamic pathway models that are individualized to patient-specific data.

Using Pasmopy, we built a mechanistic model of ErbB receptor signaling network, trained with protein quantification data obtained from cultured cell lines, and performed *in silico* simulation of the pathway activities on breast cancer patients using `The Cancer Genome Atlas (TCGA) <https://www.cancer.gov/about-nci/organization/ccg/research/structural-genomics/tcga>`_ transcriptome datasets.
The temporal activation dynamics of Akt, extracellular signal-regulated kinase (ERK), and c-Myc in each patient were able to accurately predict the difference in prognosis and sensitivity to kinase inhibitors in triple-negative breast cancer (TNBC).

Protocol
--------

This protocol describes in detail the step-by-step method for model construction of the ErbB signaling network, parameterization of the models, integration of transcriptomic data, and stratification of breast cancer patients based on signaling dynamics.

* **Paper:** https://doi.org/10.1016/j.xpro.2022.101619

* **Code:** https://github.com/pasmopy/breast_cancer