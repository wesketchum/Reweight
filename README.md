# GENIE Reweight

The GENIE Reweight product is a selection of tools to propagate model uncertainties and to support generator-related 
analysis tasks. Users should note the inherent limitations of the reweighting procedure and be aware that this product
does not include weight calculators for several important systematics, and that important modelling aspects are not
reweightable in principle. The GENIE Reweight product does *not* provide the full systematic error for any GENIE 
comprehensive model or tune and, indeed, the **GENIE tuning procedure makes no use of the ReWeight product** 
but it relies on weighting functions from brute-force parameter scans made with the aid of the Professor tool [https://professor.hepforge.org/].

The GENIE Collaboration has medium-term plans to release Professor/YODA Generator response functions for all important 
modelling uncertainties, as well as to release all covariance matrices from the GENIE global fits to scattering data.
In the mean time, please be aware of the numerous caveats with the use of the Reweight product.

For more information, visit http://www.genie-mc.org

<pre>
                                                   .oooooo.    oooooooooooo ooooo      ooo ooooo oooooooooooo  
                                                  d8P'  `Y8b   `888'     `8 `888b.     `8' `888' `888'     `8  
                                                 888            888          8 `88b.    8   888   888           
                        Ndyooym          dN      888            888oooo8     8   `88b.  8   888   888oooo8     
                     Nds//+sdmoy       d+m       888     ooooo  888    "     8     `88b.8   888   888    "      
                   Nh+//ohN  m+s      N//syyN    `88.    .88'   888       o  8       `888   888   888       o  
                 Ny+//od   Nh+oN       o///+      `Y8bood8P'   o888ooooood8 o8o        `8  o888o o888ooooood8  
               Nh+//om   Nh+/yN       o///s                                                                    
              d+//+d   my+/smmyhN    m///h                                                           REWEIGHT    
            Ns///yN NdyoshNNs///d   h////yN                                                                    
           mo//om        ms///+m   d///////oyhmN                                                               
          N+//yN       ms////+N    h////////////oym                                                         
          s//h       ho/+///sN     N///////////////od                                                          
         N+/h     my++yh+//y       s/////////////////oN                                                        
         Ns/m Nmy++ymNh//+d        s//////////////////+m                                                       
          NhssoshN  Ny//sN         m///////////////////+                                                       
                  Nmo/+ohdmN        mo//////////////////h            NmmN                                      
              ms/+s//o/-----:+sdN     mhso+ooys/////////o      mhs+/------/ohN                                 
    Nhd     N+-/o+/o+------------/shN     Ndy+//////////+ mhs+:---------------om                               
  mo/sN    m:/oo/oo:-----------------:://////////////////-----------------------y                              
  y//d    Noo++o+:----------:------------:://///////////:-------:/+osyso/:-------h                             
  Nyo+syysooo+/--------------::--------------:://///////:----::////+ossyso/------:                             
     NNNNNs--------------------:/::-------------:://///:-:://////////oosyo+:------                             
          y---------------------:////:::-----------:////////////////+o++++/:-----/                             
          m-----------------------://////////////////////////////+so/--:---------y                             
           +------------------------://////////////////////////+yy:---+oh/--:y+-/N                             
           N/-------------------------://////////////////////ohy/-------yy--os-/m                              
            No--------------------------://///////////////+shs/---------/d++/-oN                               
              mo--------------------------:////////////+shyo:------------sy/sm                                 
                Ny+:------------------------::://///oyhyo:--------------/sd                                    
                    mhso+/:------------------:/+oyhyo/----------:/+syhm                                        
                           NmddhyyyssssyyyhdmmNNNmhhhyyyyhhddmN                                                
</pre>

## Current authors:

- Luis Alvarez-Ruso (*IFIC*)
- Costas Andreopoulos (+) (*Liverpool*)
- Adi Ashkenazi (*Tel Aviv*)
- Joshua Barrow (*Tel Aviv; MIT*)
- Steve Dytman (*Pittsburgh*)
- Hugh Gallagher (*Tufts*)
- Alfonso Andres Garcia Soto (*Harvard and IFIC*)
- Steven Gardiner (*Fermilab*)
- Matan Goldenberg (*Tel Aviv*)
- Robert Hatcher (*Fermilab*)
- Or Hen (*MIT*)
- Igor Kakorin (*JINR*)
- Konstantin Kuzmin (*ITEP and JINR*)
- Weijun Li (*Oxford*)
- Xianguo Lu (*Warwick*)
- Anselmo Meregaglia (*Bordeaux, CNRS/IN2P3*)
- Vadim Naumov (*JINR*)
- Afroditi Papadopoulou (*Argonne*)
- Gabriel Perdue (*Fermilab*)
- Komninos-John Plows (*Oxford*)
- Marco Roda (*Liverpool*)
- Beth Slater (*Liverpool*)
- Alon Sportes (*Tel Aviv*)
- Noah Steinberg (*Fermilab*)
- Vladyslav Syrotenko (*Tufts*)
- Júlia Tena Vidal (*Tel Aviv*)
- Jeremy Wolcott (*Tufts*)
- Qiyu Yan (*UCAS and Warwick*)

---
(+) Corresponding Author:

**Prof. Costas Andreopoulos < c.andreopoulos \at cern.ch >**

University of Liverpool, Department of Physics, Oliver Lodge Lab 316,  Liverpool L69 7ZE, UK  

 
## Past authors and other key contributors

Past authors: 
- Christopher Barry (*Liverpool*)
- Steve Dennis (*Liverpool*)
- Walter Giele (*Fermilab*)
- Timothy Hobbs (*Fermilab*)
- Libo Jiang (*Pittsburgh*)
- Rhiannon Jones (*Liverpool*)
- Donna Naples (*Pittsburgh*)
- Julia Yarba (*Fermilab*) 


## Copyright

Copyright (c) 2003-2023, The GENIE Collaboration. For information, visit http://copyright.genie-mc.org


## Physics & User manual

For installation and usage information, as well as information on the GENIE framework, event generator modules and tuning, 
see the GENIE Physics & User Manual in the public section of the GENIE Document Database:
https://genie-docdb.pp.rl.ac.uk/cgi-bin/ShowDocument?docid=2


## Public releases and physics tunes

For a list of public releases and a summary information, see http://releases.genie-mc.org

For a list of model configurations and tunes supported in each release, see http://tunes.genie-mc.org

The naming conventions for releases, model configurations and tunes are outlined in the GENIE web page and in the Physics & User manual.


## Contribution guidelines

GENIE welcomes community contributions through its Incubator. An Incubator Project is the unique route for any physics or software development into any of the GENIE suite products (Generator, Comparisons, Tuning). Details on the Incubator Project Phases (Launch, Research and Development, Graduation, Integration and Deployment) can be found in the GENIE Bylaws in the public section of the GENIE Document Database: https://genie-docdb.pp.rl.ac.uk/cgi-bin/ShowDocument?docid=1


## Citing GENIE

If you use GENIE, please **always** cite the following reference: 

<pre>
@article{Andreopoulos:2009rq,
      author         = "Andreopoulos, C. and others",
      title          = "{The GENIE Neutrino Monte Carlo Generator}",
      journal        = "Nucl. Instrum. Meth.",
      volume         = "A614",
      year           = "2010",
      pages          = "87-104",
      doi            = "10.1016/j.nima.2009.12.009",
      eprint         = "0905.2517",
      archivePrefix  = "arXiv",
      primaryClass   = "hep-ph",
      reportNumber   = "FERMILAB-PUB-09-418-CD",
      SLACcitation   = "%%CITATION = ARXIV:0905.2517;%%"
}
</pre>

If you used any of the standard GENIE applications, built-in flux and geometry drivers, or if you used any of its event reweightng and error propagation tools, please **add the following reference**:
<pre>
@article{Andreopoulos:2015wxa,
      author         = "Andreopoulos, Costas and Barry, Christopher and Dytman,
                        Steve and Gallagher, Hugh and Golan, Tomasz and Hatcher,
                        Robert and Perdue, Gabriel and Yarba, Julia",
      title          = "{The GENIE Neutrino Monte Carlo Generator: Physics and
                        User Manual}",
      year           = "2015",
      eprint         = "1510.05494",
      archivePrefix  = "arXiv",
      primaryClass   = "hep-ph",
      reportNumber   = "FERMILAB-FN-1004-CD",
      SLACcitation   = "%%CITATION = ARXIV:1510.05494;%%"
}
</pre>

Finally, if you used any of the new model configurations and tunes provided in the GENIE v3* series, please consider adding any of the following references is relevant:

<pre>
@article{GENIE:2021npt,
    author = "Alvarez-Ruso, Luis and others",
    collaboration = "GENIE",
    title = "{Recent highlights from GENIE v3}",
    eprint = "2106.09381",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    reportNumber = "FERMILAB-PUB-21-266-SCD-T",
    doi = "10.1140/epjs/s11734-021-00295-7",
    journal = "Eur. Phys. J. ST",
    volume = "230",
    number = "24",
    pages = "4449--4467",
    year = "2021"
}
</pre>

<pre>
@article{GENIE:2021zuu,
    author = "Tena-Vidal, J\'ulia and others",
    collaboration = "GENIE",
    title = "{Neutrino-nucleon cross-section model tuning in GENIE v3}",
    eprint = "2104.09179",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    reportNumber = "FERMILAB-PUB-20-531-SCD-T",
    doi = "10.1103/PhysRevD.104.072009",
    journal = "Phys. Rev. D",
    volume = "104",
    number = "7",
    pages = "072009",
    year = "2021"
}
</pre>

<pre>
@article{GENIE:2021wox,
    author = "Tena-Vidal, J\'ulia and others",
    collaboration = "GENIE",
    title = "{Hadronization model tuning in genie v3}",
    eprint = "2106.05884",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    reportNumber = "FERMILAB-PUB-21-024-QIS-SCD-T",
    doi = "10.1103/PhysRevD.105.012009",
    journal = "Phys. Rev. D",
    volume = "105",
    number = "1",
    pages = "012009",
    year = "2022"
}
</pre>

<pre>
@article{GENIE:2022qrc,
    author = "Tena-Vidal, Julia and others",
    collaboration = "GENIE",
    title = "{Neutrino-nucleus CC0$\pi$ cross-section tuning in GENIE v3}",
    eprint = "2206.11050",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    reportNumber = "FERMILAB-PUB-22-296-ND-QIS-SCD",
    doi = "10.1103/PhysRevD.106.112001",
    journal = "Phys. Rev. D",
    volume = "106",
    number = "11",
    pages = "112001",
    year = "2022"
}
</pre>

Please notice that the GENIE authors endorse the **MCNET guidelines for fair academic use** which can be found in http://www.montecarlonet.org/GUIDELINES. We invite users to consider which GENIE components are important for a particular analysis and cite them, in addition to the main references. A list of such references in maintained in the official GENIE web page.
