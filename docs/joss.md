---
title: 'Staring down direct infusion metabolomics datasets with MeDUSA'
tags:
  - R
  - Single Cell Metabalomics
  - Mass Spectrometry
  - Metabolomics
  - Data Analytics

authors:
  - name: Laura Ann Hetzel
    orcid: 0000-0002-6922-9423
    equal-contrib: true
    affiliation: 1
  - name: Eric Hetzel
    equal-contrib: true
    affiliation: 2
  - name: Ahmed Ali
    orcid: 0000-0003-2157-4399
    affiliation: 1
    corresponding: true
affiliations:
 - name: Leiden University, Leiden, Netherlands
   index: 1
 - name: Independent Researcher, Netherlands
   index: 2
date: 4 June 2024
bibliography: joss.bib
---

# Summary
Over ten trillion cells are hard at work in the human body[@B:2013] and there can be significant heterogeneity amongst them affecting biological development, disease progression, and treatment response.[@ZV:2018]

One technique to categorize this heterogeneity is single cell analysis using mass spectrometry. However this technique introduces new challenges. One of which being the small sample volume limits separation possibilities, this makes the data non-compatible with traditional analysis pipelines. We introduce the R package `MeDUSA` for Metabolomic Direct-infusion Untargeted Single-cell Analysis.

`MeDUSA` is a start-to-finish analysis package allowing metabolomics researchers to focus on analytical content rather than R proficiency. `MeDUSA` handles the suggested workflow by Southam [@SWEJV:2017] from data extraction to filtering without a chromatogram, carrying the data through statistical analysis and identification of features for biological interpretation.

# Statement of Need
Due to the small volume of a single cell, direct infusion(DI), nano-electrospray ionization (nESI) is a highly suitable technique. However liquid chromatography mass spectrometry (LC-MS) is significantly more common than DI mass spectrometry, and therefore drives software development. This leads to most software being dependent on a separation chromatogram. For instance, XCMS offers filtering and statistical analysis with visualization, however, all of the filtering is reliant on the presence of a chromatogram. [@SWOS:2006] Similarly, MetaboAnalyst has an impressive interface that allows for metabolomic pathway analysis, but the full potential is only unlocked with a chromatogram.[@P:2021] Even paid subscription based software such as Thermo Fisher Scientific’s Compound Discoverer is designed to align and filter peaks based on a chromatogram.[@CY:2024] Furthermore, the mentioned software does not allow the user to define preprocessing methods such as centroiding or alignment. In contrast, the modularity of `MeDUSA` is built specifically for direct infusion data, and it lays the framework for method specification. The modularity also offers the user the ability to bypass functions, introduce external functions to filter, and process the data as they see fit, an option that other software does not offer. Currently, if researchers want this level of modularity and customization designed for single cell data, they must write their own scripts which requires proficiency in a programming language and a significant time allotment. Therefore the metabolomic community needs a software option that will enable complete, modular, and customizable processing of mass spectra without a chromatogram.

# Description
The goal of `MeDUSA` is to provide a toolset that is modular, customizable, and user friendly. There are five major sections along this standard flow as shown in \autoref{fig:summary}.

Modularity is achieved by using standardized interoperable data-objects. This allows the user to choose any collection and order of functions as they see fit.  See the README for a list and description of the standard objects.

Customization is achieved by being greatly parameterized and leveraging modularity. This allows the user to dial in their variables, such as thresholds, aggregation methods, tolerances, etc. The user may also interrupt the suggested flow to perform any custom logic to their needs, and reintroduce their updated data into the `MeDUSA` flow, so long as the data structure is maintained. The three primary data-objects have the same structure and are differentiated by name for human readability along the suggested flow. However, they are technically interchangeable, thus increasing customization.

User-friendliness is achieved with “magic” functions, readability, and containerization. The magic functions leverage suggested parameter values and simplify many functions within each of the five major sections. These methods can allow a user to go from mzML files to a list of compounds in few commands. \autoref{fig:detailed} illustrates a detailed list of functions, and how the magic functions simplify them. To ease readability, files and methods are prepended with the expected input data-object type. This naming convention helps users identify what methods are available to them at different stages of the suggested flow. Containerization manages dependencies and provides HMDB and lipid data for convenient m/z to compound mapping.

![Map of the five sections of MeDUSA and the capabilities of each section. The bold arrows indicate a suggested workflow. The dashed arrows indicate references. The circled text indicated the object data type. The plot symbol indicates the function may output a plot.\label{fig:summary}](medusa-medusa.png)

![Detailed suggested flow and function map. Note the “magic” functions which aggregate similar functions for user ease. The bold arrows indicate a suggested workflow. The circled text indicated the object data type. The plot symbol indicates the function may output a plot. \label{fig:detailed}](medusa-medusa_detailed.png)

# Research projects using the software

Current research projects relating to single cell metabolic profiling of the cell cycle of FUCCI cells, stem cell differentiation, hypoxic organoids, and metastatic organoids would benefit from this package. The projects currently use two different commercially available software as well as R scripts that lack robustness for the data analysis; shifting the analysis to MeDUSA will enable timely and reliable results. Upon the release of the package, it will be implemented lab wide for live single cell metabolomics.

# Acknowledgements

This research was funded in part by the European Union’s Horizon 2020 research and innovation program under the Marie Skladowska-Curie grant agreement number 861196.

Rosa (species = Canis Familiaris, breed = Catahoula Leopard) was essential to the mental well-being and creativity of the authors, and therefore essential to the project.

We would like to thank Alida Kindt and Thomas Hankemeier for their support and direction on this project.

# References
