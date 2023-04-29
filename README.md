# ScPectral
A tool for automated discovery of developmental pathways.<br>
<pre>
Environment setup:
    conda create -n scPectral python=3.7.16
    conda activate scPectral
    pip install hypernetx==1.2.5
    conda install -c anaconda ipykernel==6.15.2
    conda install -c anaconda scikit-learn==1.0.2
</pre>

<pre>
Analysis of mouse embryonic stem cells (mESCs) <sup>1</sup> in mESCsResult notebook.
Analysis of murine hematapoietic stem and progenitor cells (HSPCs)<sup>2</sup> in HSPCsResult notebook. 

HyperNetX package used for clustering and visualisation of hypergraph<sup>3</sup>.
LocaTE <sup>4</sup> was used to obtain input networks for ScPectral. 

References:
[1] Hayashi, T., Ozaki, H., Sasagawa, Y. et al. Single-cell full-length total RNA sequencing uncovers dynamics of recursive splicing and enhancer RNAs.
    Nat Commun 9, 619 (2018). https://doi.org/10.1038/s41467-018-02866-0.

[2] Marot-Lassauzaie V, Bouman BJ, Donaghy FD, Demerdash Y, Essers MAG, Haghverdi L (2022) Towards reliable quantification of cell state velocities. 
    PLoS Comput Biol 18(9): e1010031. https://doi.org/10.1371/journal.pcbi.1010031.

[3] https://github.com/pnnl/HyperNetX

[4] Stephen Y Zhang, Michael P H Stumpf. Learning cell-specific networks from dynamical single cell data. 
    bioRxiv 2023.01.08.523176; doi: https://doi.org/10.1101/2023.01.08.523176.
</pre>
