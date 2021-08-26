# Tumor-Cell-Extravasation-Simulation

The adhesion of tumor cells to the endothelium is a crucial step in the extravasation process of cancer migration and metastasis. Tumor cell adhesion involves several complex mechanisms and factors, including receptor-ligand kinetics, rupture forces, loading rates, and statistical mechanics, to name just a few. Computer simulations are an invaluable method of modeling these systems and understanding the relationships between the various factors affecting the process of adhesion to the endothelium. This paper discusses important and current research into the adhesion process within tumor cell extravasation, details a computer-simulated parametric case study of the factors and statistical mechanics involved in adhesion, and highlights potential future directions for research on the adhesion process involving the usage of computer modeling and simulations.


Interactions between receptors and ligands are a key component of tumor cell adhesion to the endothelial cell layer. In particular, the zero-force unbinding rate, also called the initial unbinding rate, plays a major role in the process. The zero-force unbinding rate can be viewed as the natural rate at which the receptor-ligand bond will be broken due to random chance as opposed to being broken by external forces applied to the bond. The value of the zero-force unbinding rate directly affects several other important factors and statistics, such as the most likely rupture force and the bond survival probability (Evans and Calderwood). However, the value of the zero-force unbinding rate is itself widely varied across different types of tumor cells and bonds, as its value is affected by a multitude of factors. These factors include the thermal energy, temperature, pressure, pH, flow characteristics, and hydrodynamic shear, as well as others that may impact the structural integrity of the bond between the tumor cell and the endothelial cell layer (Barsegov and Thirumalai). 


Unfortunately, this means that experimentally determined values for the zero-force unbinding rate will vary depending on the specific conditions of the local microenvironment under which the testing was performed, as well as the methods used to measure the value of this rate. For example, using highly similar experimental testing setups, different groups of researchers used varying methods to measure the zero-force unbinding rate of the LFA-1/ICAM-1 bond. Three different values were found by the three different groups of researchers, likely due to slight differences in the local microenvironments as well as various inaccuracies or effects of the measurement methods used. These values for the zero-force unbinding rate of this bond was measured to be 0.1 s-1 using a surface plasmon resonance assay technique, 0.17 s-1 using atomic force microscopy, and 0.3 s-1 using a parallel plate flow chamber (Hoskins and Dong). These differing results demonstrate the complexity of the effects of receptor-ligand kinetics on tumor cell adhesion.


One of the most heavily researched areas around the study of tumor cell adhesion is the determination of bond rupture forces. The rupture force, also called the detachment force or deadhesion force, is the force necessary to stretch bond between the tumor cell and the endothelial cell until the bond detaches or ruptures. Stronger bonds will therefore have, on average, higher rupture forces (Laurent, et al.). 
There are several methods used to experimentally measure the rupture force of specific types of bonds. One of the most common methods is known as atomic force microscopy. In this method, a tumor cell is attached to the tip of a cantilever and the cantilever is moved downward at a constant velocity toward a constructed layer of endothelial cells. The tumor cell is allowed to come into contact with the endothelial cell layer, and once this occurs, the cantilever stops moving downward. The tumor cell is held in place for a predetermined amount of time, usually on the order of a few seconds, in order to give adequate time for complete or at least strong adhesion of the tumor cell to the endothelial cell layer through formation of bonds between receptors and ligands. After this predetermined period of time, the cantilever is moved upward at a constant velocity, and the deflection of the cantilever as it moves upward is recorded. The cantilever is moved upward until the bonds between the tumor cell and the endothelial cell layer have ruptured, causing deadhesion. These force data are used to determined the force curve for the experiment, which can then be used to determine the force at which the bonds ruptured and the tumor cell was completely detached from the endothelial cell layer (Khalili and Ahmad). A simple diagram illustrating how an atomic force microscopy experiment would be conducted to measure the rupture force of bonds between a tumor cell and an endothelial cell layer is shown in Figure 1. 

![Figure 1](https://i.imgur.com/ygaan7n.png)

Figure 1: The three phases of an atomic force microscopy experiment (Benoit, et al.).

The value of the average rupture force varies with the type of bond, which is not always as simple as merely the specific receptor-ligand pair involved in the bond between the tumor cell and endothelial cell. In fact, there are multiple types of adhesive bonds, including persistent connectors, transient connectors, and catch bonds, which are a special subtype of transient connector. Persistent connectors, such as ICAM-1/αLβ2 bonds, have relatively high rupture force values even at lower force loading rates. This is not the case for transient connectors, such as VCAM-1/VLA-4 bonds,which require significantly higher force loading rates in order to exhibit larger rupture forces, suggesting that these types of bonds are relatively weak at low force loading rates. Catch bonds, such as PSGL-1/P-selectin bonds, are unique in that these bonds have essentially no strength under conditions of slow force loading or low force values, leading to a quick release of the bond (Evans and Calderwood). A graphical representation of these types of adhesive bonds using plots of rupture force vs. loading rate is shown in Figure 2.

![Figure 2](https://i.imgur.com/s3qLaoI.png).

Figure 2: Plots of bond strength vs. force loading rate illustrating three different types of adhesive bonds (Evans and Calderwood).

The rupture force values for bonds comprised of receptor-ligand interactions can range from approximately 20 pN to 200 pN. For comparison, the forces generated by molecular and cytoskeletal motor proteins range from about 2 pN to 10 pN, while the rupture force of a covalent bond has been found to be approximately 4 nN, or 4000 pN (Müller, et al.). 

A computer simulation of the adhesion process using MATLAB software was performed to illustrate the concepts of statistical mechanics that underlie the process. The simulation was a parametric case study based on the results of an atomic force microscopy experiment attempting to determine the average rupture forces of single bonds comprised of interactions between the VLA-4 receptor of B16 melanoma cells and the VCAM-1 counter-receptor of bEnd.3 endothelial cell monolayers performed by R.H. Eibl and M. Benoit, as well as research performed by Xiaohui Zhang and colleagues. 

The experiment performed by R.H. Eibl and M. Benoit determined that under physiological conditions of a temperature of 37 °C and a pH of 7.4, single VCAM-1/VLA-4 bonds have an average rupture force of 33 pN. The standard deviation of this experimental measurement was 12 pN. Eibl and Benoit measured the rupture forces of 960 VCAM-1/VLA-4 bonds. They included a histogram plot of the rupture forces, which exhibited a standard Gaussian distribution (Eibl and Benoit).
 
Xiaohui Zhang and colleagues performed a much more in-depth study of the effects of various mutations of the wild type VCAM-1/VLA-4 bond on rupture forces, zero-force unbinding rates, and characteristic barrier lengths. It was understood that the energy landscape of the VCAM-1/VLA-4 complex features two energy barriers. However, these researchers concluded that when dealing with forces around the order of 50 pN or less, the first energy barrier has little effect on the overall energy landscape and system characteristics. The researchers made comparisons of the zero-force unbinding rate and characteristic barrier length, which exist in two forms related both to the first and second energy barrier, but the second energy barrier’s zero-force unbinding rate and characteristic barrier length dominate at lower forces (on the order of roughly 50 pN), so these factors are taken into consideration while the characteristic barrier length and zero-force unbinding rate of the first energy barrier in the VCAM-1/VLA-4 complex are largely ignored, as the overall effects on the system are negligible. The researchers performed experiments using different mutations or conditions, including the presence of 5 mM ethylenediaminetetraacetic acid (EDTA), a mutation of the aspartic acid-40 residue of VCAM-1 to a neutral, uncharged alanine residue (D40A), and a mutation of the aspartic acid-40 residue of VCAM-1 to a negatively charged glutamate residue (D40E). 

Due to the complex system that is the human physiology, numerous assumptions were made in order to create a useful simulation. Firstly, the simulation was only based on the adhesion step in the extravasation process. For cancer to metastasize, it must intravasate, diffuse through the blood, extravasation out of the endothelium, and proliferate in another location. The simulation was chosen to only simulate the adhesion process because it is the critical step in cancer metastasis. If tumor cells cannot adhere to the endothelium, it cannot proliferate and become lethal. 

Another limitation on the simulation is the assumption made on the average rupture forces. This value will change dramatically between physiologies. Things that will vary this value are the type of tumor cell, temperature, pH, tumor cell concentrations, endothelium wall strength, presence of certain proteins, and much more. The value of 33 pN was chosen because it can be considered a typical value. Research shows the average value can be from 20 pN to 200 pN [Muller et al.].

Time is a variable that vary on several orders of magnitude. For this simulation, it was assumed to take place over 10 seconds. Cancer metastasis is a process that occurs over a period of months or years, but the adhesion process is a quicker process. 10 seconds was assumed to preserve computational resources. Due to the simplicity of the simulation, the length of time it runs doesn’t fundamentally change the results. The rupture forces still have a normal distribution about the set average, and the survivability of bonds still approximately linear. 

The number of tumor cells used in the simulation is a value that severely limited the capabilities of the simulation. Again, 1000 tumor cells was assumed in order to preserve computational resources. The concentration of tumor cells on the endothelium alters the entire system drastically. Due to the properties of tumor cell endothelium interactions, an increase in concentration acts as a positive feedback that allows tumor cells to extravasate more quickly. These dynamics were not taken into account for this simulation, due to the not fully understood mechanics behind it. It is reasonable to assume, however, that an increase in the number of tumor cells would decrease the percentage of adhered tumor cells and decrease the average rupture force, as can be shown by receptor ligand kinetics.

The future of cancer research is uncertain. There are thousands of efforts worldwide to reduce the lethality of cancer. Every type of cancer is a unique disease, with its own complications and dynamics, so the area of research is vast. An exciting branch of research is the design of nanotechnology specifically for cancer. There are several types of nanotechnologies already available that are in a limited use: liposomes, synthetic polymers, chitosan, buckyballs, nanotubes, mesoporous silica and more. Their designs generally are for drug delivery, due to the fact that nanotechnology allows for structures impossible with traditional medicine. Other practical uses of nanotechnology on cancer treatment are for diagnostics, biomarkers, therapeutic delivery, and immunotherapeutics.


### References

Barsegov, V., and D. Thirumalai. “Dynamics of Unbinding of Cell Adhesion Molecules: Transition from Catch to Slip Bonds.” Proceedings of the National Academy of Sciences, vol. 102, no. 6, 2005, pp. 1835–1839., doi:10.1073/pnas.0406938102. 

Benoit M, Gaub H, E: Measuring Cell Adhesion Forces with the Atomic Force Microscope at the Molecular Level. Cells Tissues Organs 2002;172:174-189. doi: 10.1159/000066964

Eibl, R H, and M Benoit. “Molecular resolution of cell adhesion forces.” IEE proceedings. Nanobiotechnology vol. 151,3 (2004): 128-32. doi:10.1049/ip-nbt:20040707

Evans, Evan A, and David A Calderwood. “Forces and bond dynamics in cell adhesion.” Science (New York, N.Y.) vol. 316,5828 (2007): 1148-53. doi:10.1126/science.1137592

Hoskins, Meghan H, and Cheng Dong. “Kinetics analysis of binding between melanoma cells and neutrophils.” Molecular & cellular biomechanics : MCB vol. 3,2 (2006): 79-87.

Khalili AA, Ahmad MR. A Review of Cell Adhesion Studies for Biomedical and Biological Applications. International Journal of Molecular Sciences. 2015; 16(8):18149-18184. https://doi.org/10.3390/ijms160818149

Laurent, Valérie M et al. “Atomic force microscopy reveals a role for endothelial cell ICAM-1 expression in bladder cancer cell adherence.” PloS one vol. 9,5 e98034. 23 May. 2014, doi:10.1371/journal.pone.0098034

Müller, Daniel J et al. “Force probing surfaces of living cells to molecular resolution.” Nature chemical biology vol. 5,6 (2009): 383-90. doi:10.1038/nchembio.181

Zhang, Xiaohui et al. “Molecular basis for the dynamic strength of the integrin alpha4beta1/VCAM-1 interaction.” Biophysical journal vol. 87,5 (2004): 3470-8. doi:10.1529/biophysj.104.045690
