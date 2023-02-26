# Twistronics

Julia codes for computing electronic bandstructure in twisted 2D systems

# Installation
Download the repository, which is compatible with Julia 1.7

# Disclaimer
This package is still under heavy development

# Examples
Codes can be found in the examples folder
## Twisted bilayer graphene (TBG)
### Supercell and bandstructure
![Alt text](figs/TBG.png?raw=true "TBG" )
(a) Supercell of TBG at twist angle $1.89^\circ$. (b)-(d) Bandstructure of TBG at twist angle (b) $\theta=21.79^\circ$, (c) $\theta=9.43^\circ$ and (d) $\theta=1.89^\circ$.

### Interlayer bias and layer operator
![Alt text](figs/tbg_bias.png?raw=true "TBG_bias" )
Bandstructure of TBG with interlayer bias $V=0.1t$ at twist angle (a) $\theta=9.43^\circ$ and (b) $\theta=1.89^\circ$. The color indicates either the bands belong to the first/second layer ($\langle L\rangle=\pm 1$) or mixed ($\langle L\rangle=0$).

### Interlayer bias and valley operator
![Alt text](figs/tbg_valley.png?raw=true "TBG_valley" )
(a) Bandstructure of TBG with a small valley potential $0.05\mathcal{V}_z$. (b) Bandstructure of TBG with an interlayer bias $V=0.1t$. The color indicates the valley polarization.

## Twisted Dirac quantum spin liquid
![Alt text](figs/TBQSL.png?raw=true "TBQSL")
(a) Part of the supercell of the twisted $\pi$-flux model at twist angle $3.48^\circ$ under a fixed gauge. (b)-(d) Bandstructure of the twisted $\pi$-flux model at twist angle (b) $\theta=21.79^\circ$, (c) $\theta=9.43^\circ$ and (d) $\theta=3.48^\circ$.
