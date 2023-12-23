<p align="center">
  <a href="" rel="noopener">
 <img width=200px height=200px src="https://i.imgur.com/6wj0hh6.jpg" alt="Project logo"></a>
</p>

<h3 align="center">Project Title</h3>

<div align="center">

[![Status](https://img.shields.io/badge/status-active-success.svg)]()
[![GitHub Issues](https://img.shields.io/github/issues/kylelobo/The-Documentation-Compendium.svg)](https://github.com/kylelobo/The-Documentation-Compendium/issues)
[![GitHub Pull Requests](https://img.shields.io/github/issues-pr/kylelobo/The-Documentation-Compendium.svg)](https://github.com/kylelobo/The-Documentation-Compendium/pulls)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](/LICENSE)

</div>

---

<p align="center"> Few lines describing your project.
    <br> 
</p>

## 📝 Table of Contents

- [About](#about)
- [Getting Started](#getting_started)
- [Deployment](#deployment)
- [Usage](#usage)
- [Built Using](#built_using)
- [TODO](../TODO.md)
- [Contributing](../CONTRIBUTING.md)
- [Authors](#authors)
- [Acknowledgments](#acknowledgement)

## 🧐 About <a name = "about"></a>

Welcome to tR3D! This tool enables the efficient exploration of gene neighborhoods surrounding enzymes of interest, for the goal of Biosynthetic Gene Cluster discovery. This enzyme family-based approach uses data directly from the Enzyme Function Institute (https://efi.igb.illinois.edu/) to identify new biosynthetic pathways and new enzyme functions. 

This tool was developed by Doug Millar at the University of California, Berkeley in the lab of Michelle Chang, with the help of Shannon Millar (B.S. Industrial Engineering, UC-Berkeley). Yes, we're married :D

## 🏁 Getting Started <a name = "getting_started"></a>

Start by generating a SSN of your enzyme family of choice, with a sequence alignment cutoff score that separates known enzymes from unknown (to the best of your knowledge). Submit this SSN to the EFI-GNT tool with 10 neighboring genes and 0% co-occurrence cutoff. Upload the Cluster SSN hub node file that is generated to Cytoscape, and download the table. This CSV is directly loaded to tR3D for visualization. That’s pretty much it. Look for high co-occurrences > 75% and low median distances (<5).

A note: If your SSN is very large, the co-occurrences will be very low and lack meaningful information. To address this, you can generate subcluster files through uploading all of the Uniprot IDs from a particular cluster for their own SSN, increase the alignment score, and repeat the process above. This "cluster of cluster" approach allows a deeper dive into clades of enzymes. 


### Prerequisites

What things you need to install the software and how to install them.

```
Give examples
```

### Installing

A step by step series of examples that tell you how to get a development env running.

Say what the step will be

```
Give the example
```

And repeat

```
until finished
```

End with an example of getting some data out of the system or using it for a little demo.

## 🔧 Running the tests <a name = "tests"></a>

Explain how to run the automated tests for this system.

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## 🎈 Usage <a name="usage"></a>

Add notes about how to use the system.

## 🚀 Deployment <a name = "deployment"></a>

Add additional notes about how to deploy this on a live system.

## ⛏️ Built Using <a name = "built_using"></a>

- [MongoDB](https://www.mongodb.com/) - Database
- [Express](https://expressjs.com/) - Server Framework
- [VueJs](https://vuejs.org/) - Web Framework
- [NodeJs](https://nodejs.org/en/) - Server Environment

## ✍️ Authors <a name = "authors"></a>

- [@kylelobo](https://github.com/kylelobo) - Idea & Initial work

See also the list of [contributors](https://github.com/kylelobo/The-Documentation-Compendium/contributors) who participated in this project.

## 🎉 Acknowledgements <a name = "acknowledgement"></a>

- Hat tip to anyone whose code was used
- Inspiration
- References
