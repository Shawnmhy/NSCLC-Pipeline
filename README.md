<!--
*** Thanks for checking out the Best-README-Template. If you have a suggestion
*** that would make this better, please fork the repo and create a pull request
*** or simply open an issue with the tag "enhancement".
*** Thanks again! Now go create something AMAZING! :D
***
***
***
*** To avoid retyping too much info. Do a search and replace for the following:
*** github_username, repo_name, twitter_handle, email, project_title, project_description
-->



<!-- PROJECT SHIELDS -->
<!--
*** I'm using markdown "reference style" links for readability.
*** Reference links are enclosed in brackets [ ] instead of parentheses ( ).
*** See the bottom of this document for the declaration of the reference variables
*** for contributors-url, forks-url, etc. This is an optional, concise syntax you may use.
*** https://www.markdownguide.org/basic-syntax/#reference-style-links
-->

<!-- PROJECT LOGO -->
<br />
<p align="center">
  <a href="https://github.com/Shawnmhy/NSCLC-Pipeline/">
    <img src="./images/logo.png" alt="Logo" width="250" height="180">
  </a>

  <h3 align="center">NSCLC-mIF processing pipeline</h3>
</p>



<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary><h2 style="display: inline-block">Table of Contents</h2></summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
    </li>
    <li>
      <a href="#system-requirement">System requirement</a>
      <ul>
        <li><a href="#operating-systems">Operating systems</a></li>
        <li><a href="#software-dependencies">Software dependencies</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgements">Acknowledgements</a></li>
  </ol>
</details>


<!-- ABOUT THE PROJECT -->
## About The Project

**This project proposed a computational pipeline to analyze the tumor microenvironment in patients with non-small cell lung cancer (NSCLC) using archival pathological samples and multiplexed immunofluorescence. 



<!-- SYSTEM REQUIREMENT -->
## System requirement

<!-- OPERATING SYSTEMS -->
### Operating systems

The computational pipeline was developed on the following operating system:
* Windows 10 Pro.
* Processor: 12th Gen Intel(R) Core(TM) i7-12700K   3.61 GHz
* 32GB RAM
* 64-bit operating system, x64-based processor.



### Software dependencies

The computational pipeline was developed using the following softwares:
* R version 4.2.0.
* Rstudio Desktop version 1.4.
* Pycharm Python IDE version 2022.12.0.
* R packages listed at the beginning of each R script.
* BioRender


<!-- INSTALLATION GUIDE -->
## Installation guide

### Instructions
To install the pipeline, simply download the codes and run from local R compiler.

## Installation time
* Installation time all dependencies should take no longer than 30 minutes.
* Codes can be readily used upon downloading and do not require extra installations.

<!-- USAGE EXAMPLES -->
### Usage

The pipeline consists of five components:
* DataPreparation.R (Read single-cell data for each core and combine as one, compute cell densities by core, read cell boundaries data)
* NSCLC_Pipeline1.R (First-order characterizations, tSNE plot, population fractions comparison and clustering, survival analysis, patient stratification)
* NSCLC_Pipeline2.R (Cell density comparison betweeen survival group, Shannon's entropy, ratios, )
* NSCLC_Pipeline3.R (Network analysis, minimum spanning tree, pairwise Gcross function, voronoi tesselation)
* NSCLC_Pipeline4.R (RiskScale)
* NSCLC_Pipeline5.R (PD-1/PD-L1 statistics)
* NSCLC_Pipeline6.R (Analysis on validation cohort)
* Function.R (Custom R functions defined for computations)


<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE` for more information.



<!-- CONTACT -->
## Contact
Haoyang Mi - hmi1@jhmi.edu







<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/github_username/repo.svg?style=for-the-badge
[contributors-url]: https://github.com/github_username/repo_name/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/github_username/repo.svg?style=for-the-badge
[forks-url]: https://github.com/github_username/repo_name/network/members
[stars-shield]: https://img.shields.io/github/stars/github_username/repo.svg?style=for-the-badge
[stars-url]: https://github.com/github_username/repo_name/stargazers
[issues-shield]: https://img.shields.io/github/issues/github_username/repo.svg?style=for-the-badge
[issues-url]: https://github.com/github_username/repo_name/issues
[license-shield]: https://img.shields.io/github/license/github_username/repo.svg?style=for-the-badge
[license-url]: https://github.com/github_username/repo_name/blob/master/LICENSE.txt
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://linkedin.com/in/github_username
